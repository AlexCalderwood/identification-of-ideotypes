library(here) #  here_1.0.1  
library(data.table) # data.table_1.13.2
library(ggplot2) # ggplot2_3.3.2 
library(bayesplot) # bayesplot_1.7.2 
library(cowplot) # cowplot_1.0.0 
library(rstan) # rstan_2.21.2 
stan_version() # 2.21.0
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source('trait_modification_functions.R')

#library(GGally) # ggpairs plot
#library(bnlearn)
#library('Rgraphviz') # allows plot() graphNEL objects
#library('rbmn') # for the conditional distributions stuff / dependencies
#source('stan_functions.R') # functions for making and plotting Stan models and data


# set RNG 
set.seed(123)


# setup ----
in.dir <- here('id_trait_structure_output') # output of "id_trait_structure.R"
out.dir <- here('trait_modification_output')
PERT.TRAIT <- 'height' # trait(s) to simulate modifying
num.perturbations <- 30 # number of values to test
num.iterations.per.chain <- 30 # number of sampled points for each perturbed value 
TESTING <- TRUE # controls whether calculates perturbation predictions, or just loads
                # provided example precomputed values for spring ecotype, "height" trait.
                # (very slow to calculate new values!)


model.dir <- paste0(in.dir, '/data/') # where to find bnlearn models, and data used to generate them
kernel_function_file <- './GP_ARD_functions.stan' # function to define covariance kernel used in GP models

graph.dir <- paste0(out.dir, '/graphs/')
pert.graph.dir <- graph.dir
pert.model.dir <- paste0(out.dir, '/model_files/')

dir.create(pert.model.dir, recursive=T)
dir.create(graph.dir)


# load model structure and data ----
D.data2.cultivars <- readRDS(paste0(model.dir, 'bnlearn_data.rds'))
D.data2.cultivars[, rep:=seq(1,nrow(.SD)), by=.(cultivar)] # add column to keep track of individual plants

results <- readRDS(paste0(model.dir, 'bnlearn_results.rds'))



# Scale the data for model fitting ----
data.columns <- names(D.data2.cultivars)[!(names(D.data2.cultivars) %in% c('cultivar', 'rep'))]
data.means <- apply(subset(D.data2.cultivars, select=data.columns), 2, mean) # store means and sd for unscaling after modelling
data.sds <- apply(subset(D.data2.cultivars, select=data.columns), 2, sd)
D.data.unscaled <- subset(D.data2.cultivars, select=data.columns)
D.data.sc <- data.table(data.frame(scale(data.frame(subset(D.data2.cultivars, select=data.columns)))))

D.data.sc.cultivar <- D.data.sc
D.data.sc.cultivar$cultivar <- D.data2.cultivars$cultivar
D.data.sc.cultivar$rep <- D.data2.cultivars$rep



##### NULL MODEL -----
##' Will need NULL predictions for any trait which is not predicted from the other traits in a particular 
##' k-fold model (ie has no parent nodes), but which are predicted (ie has parents) in at least one other. 
##' This is so that can average predictions over the k-fold models. 
##' 
##' In the models in which they don't have parents, hese NULL prediction traits are not predicted to vary 
##' with any perturbation of the other traits,  so will just be sampled from a normal distribution with 
##' mean and sd estimated from the observed values for that trait.

# parse k-fold models to get all traits which are not root nodes in at least one model, and the subset of 
# these traits which are descendents of PERT.TRAIT
O <- get_all_predicted_and_perturbed_traits(results, PERT.TRAIT)
all.pred.traits <- O[[1]]
all.pert.pred.traits <- O[[2]]

all.data.df <- subset(D.data.sc.cultivar, select=all.pred.traits)

stan_data <- make_stan_data_NULL(all.data.df) # convert all.data.df to a List so can pass to stan

# write stan source code for null model to file - each trait is sampled from a normal distribution based
# on observed values
write.path.NULL <- paste0(pert.model.dir, PERT.TRAIT, '_NULL.stan')
make_stan_model_NULL(stan_data, write.path.NULL)
gc()

# sample using the just written model -----
NULL.fit.file <- paste0(pert.model.dir, PERT.TRAIT, '_NULL_fitted.rds') # nb shouldn't just be "_NULL.rds", as an rds file with the same name as the .stan file confuses stan!
# if don't remove previously compiled files, rstan uses old version
if (file.exists(paste0(pert.model.dir, PERT.TRAIT, '_NULL.rds'))) {
  file.remove(paste0(pert.model.dir, PERT.TRAIT, '_NULL.rds'))
}
if (file.exists(NULL.fit.file)) {
  file.remove(NULL.fit.file)
}
numChains <- 4
NULL.fit <- rstan::stan(file=write.path.NULL, data=stan_data, 
            warmup =1000, iter = 2000, chains = numChains, refresh=50)
saveRDS(NULL.fit, file=NULL.fit.file)
extracted.NULL.fitObj <- extract(NULL.fit)

# downsample NULL estimates so same dimension as real predictions will have ----
extracted.NULL.fitObj <- downsample(extracted.NULL.fitObj,
                                    numChains, num.iterations.per.chain)

# convert the "_rep" null predictions generated, to be in the same format that the real "_pred"
# predictions will be from the non NULL model. By multiplying each element by the numper of perturbations 
# to be applied (so will be in a,a,a,a,a,b,b,b,b,b,c,c,c,c,c,c format where is repeated the number of
# perturbations times.
for (trait in names(stan_data)) {
  if (trait != 'N') {
    m1<- extracted.NULL.fitObj[[paste0(trait, '_rep')]]
    tmp <- matrix(data = apply(m1, 2, function(x) rep(x, num.perturbations)), ncol = ncol(m1)*num.perturbations)
    extracted.NULL.fitObj[[paste0(trait, '_pred')]] <- tmp
  }
}
gc()

#### End of NULL modelling ---

#### Real predictions ----

# for each k-fold DAG model found
i <- 1
first.loop <- TRUE # flag to track if first k-fold
for (i in seq(from=1, to=5)) {
  print(paste0('k: ', i))
  # create k-fold directories
  curr.fold.graph.dir <- paste0(graph.dir, 'k', i, '-fold/')
  dir.create(curr.fold.graph.dir)
  curr.fold.model.dir <- paste0(pert.model.dir, 'k', i, '-fold/')
  dir.create(curr.fold.model.dir)
  
  k_model = results$model.strings[i] # get the DAG model string
  
  # iterate over perturbed traits 
  pert.trait <- 'height'
  for (pert.trait in c(PERT.TRAIT)) {
    # define where results will be written.
    curr.pert.plt.dir <- paste0(curr.fold.graph.dir, pert.trait, '/') # plots of effects of perturbations  
    dir.create(curr.pert.plt.dir)
    curr.pert.model.dir <- paste0(curr.fold.model.dir, pert.trait, '/')
    dir.create(curr.pert.model.dir)
    
    # Make input data with altered values for pert.trait ----
    # get the traits want to use real (measured) values for the perturbed data predictions.
    # should be the perturbed trait itself, and any other nodes, which are not direct or indirect descendents
    # of the perturbed trait. Produces _4pred datas in tot.stan.data, which are then used by 
    # bnmodel2stanmodel_GP_pred() for making the data, and generative sections of the stan model
    perturbed.independent.traits <- get_perturbed_independent_traits(k_model, pert.trait) 
    data.traits <- c(perturbed.independent.traits, pert.trait) 
    # make the perturbed data points. Each observed plant is repeated num.perturbations times, 
    # and each repeat is given a new value for the perturned trait. 
    pred.data.df <- make_stan_prediction_data_df_GP(D.data.sc, pert.trait, num.perturbations)
    # add plant and cultivar information to pred.data.df
    pred.data.df.cultivars <- combine_pred_data_with_plant_id(pred.data.df, D.data.sc.cultivar)
    # just get the trait data for the modelling
    pred.data.df.fixed <- subset(pred.data.df.cultivars, select=data.traits)
    tot.stan.data <- make_stan_data_GP(D.data.sc, pred.data.df.fixed, k_model, seperate.predictions = T)
    
    
    # Optimise GP prior hyperparameters ------
    # make the optimising stan model source code file
    write.path.opt <- paste0(curr.pert.model.dir, pert.trait, '_opt.stan')
    bnmodel2stanmodel_GP_opt(k_model, tot.stan.data, kernel_function_file, write.path.opt)
    # where will save fitted model object
    # get the optimal hyperparameters, and add to the stan data
    gp_opt <- stan_model(file=write.path.opt) # compile the stan source code
    opt_fit <- optimizing(gp_opt, data=tot.stan.data, seed=95848338, hessian=FALSE)
    write_opt_parameters_to_file(opt_fit, paste0(curr.pert.plt.dir, 'optimal_hyperparameters.csv'))
    tot.stan.data <- add_opt_parameters_to_stan_data(opt_fit, tot.stan.data)
    gc()
    
  
    # Predict the values of descendent traits of the perturbed parameter using optimised hyperparams  -------
    # Make the prediction stan model source code file. 
    # Direct children of perturbed trait are predicted, these predictions are used in the prediction 
    # of indirect descendents of perturbed trait etc. etc. Predictions are made accounting for order 
    # of dependencies of traits.
    write.path.pred <- paste0(curr.pert.model.dir, pert.trait, '_pred.stan')
    fit.file <- paste0(curr.pert.model.dir, pert.trait, '_fitted.rds')
    bnmodel2stanmodel_GP_pred_v2(k_model, pert.trait, tot.stan.data, kernel_function_file, write.path.pred)
    
    # make the predictions using the perturbed traits, and optimized hyperparams ------
    if (!(TESTING)) { # if not testing, calculate predictions
      fit <- stan(file=write.path.pred, data=tot.stan.data, 
                  warmup = 0, iter = num.iterations.per.chain, chains = 4, refresh=1, algorithm='Fixed_param')
    } else { # otherwise load example precomputed values for "height" trait perturbation
      print('loading example pred model fit...')
      eg.file <- here('data', 'example_trait_modification_data', 
                       paste0('k', i, '/height_fitted.rds'))
      fit <- readRDS(file=eg.file)
    }
    saveRDS(fit, file=fit.file)
    write_convergence_check(fit, paste0(curr.pert.plt.dir, 'fit_summary.csv')) # write the nEff, and Rhat to file
    extractedfitObj <- extract(fit)
    
    # augment the current k-fold predictions with any NULL ones required because are predicted in another model----
    aug.extractedFitObj <- extractedfitObj
    for (rep.trait in all.pred.traits) {
      if (!(paste0(rep.trait, '_rep') %in% names(aug.extractedFitObj))) {
        aug.extractedFitObj[[paste0(rep.trait, '_rep')]] <- extracted.NULL.fitObj[[paste0(rep.trait, '_rep')]]
      }
    }
    for (pred.trait in all.pert.pred.traits) {
      if (!(paste0(pred.trait, '_pred') %in% names(aug.extractedFitObj))) {
        aug.extractedFitObj[[paste0(pred.trait, '_pred')]] <- extracted.NULL.fitObj[[paste0(pred.trait, '_pred')]]
      }
    }
    
    # store the results across k-folds so can average predictions later ----
    if (first.loop==TRUE) {
      all.extractedFitObj <- aug.extractedFitObj
      all.tot.stan.data <- tot.stan.data
      first.loop <- FALSE
    } else {
      for (n in names(aug.extractedFitObj)) {
        all.extractedFitObj[[n]] <- rbind(all.extractedFitObj[[n]], aug.extractedFitObj[[n]])
      }
      for (n in names(tot.stan.data)) {
        if (!(n %in% all.tot.stan.data)) {
          all.tot.stan.data[[n]] <- tot.stan.data[[n]]
        }
      }
    }
    gc() 
    
  } # close pert.trait loop
} # close current k-model loop


## Combine the predictions across the k models -----------------------------------------------------

# undo the mean=0, sd=1 scaling ----
unscaled.fit.extracted <- unscale.extracted.fit(all.extractedFitObj, k_model, data.means, data.sds)
unscaled.tot.stan.data <- unscale.stan.data(all.tot.stan.data, k_model, data.means, data.sds)  
unscaled.pred.data.df <- unscale.data.frame(pred.data.df.cultivars, data.means, data.sds) # pred.data.df.cultivars is the same for every k-fold loop


# backtransform the data to undo the normalisation by e.g. log/sqrt functions to make more normally distributed ----
b.unscaled.fit.extracted <- backtransform.data(unscaled.fit.extracted)
b.D.data2.cultivars <- backtransform.data(D.data2.cultivars)
b.unscaled.pred.data.df <- backtransform.data(unscaled.pred.data.df)
b.unscaled.tot.stan.data <- backtransform.data(unscaled.tot.stan.data)
if (length(all.pert.pred.traits) > 0) {
  data.for.perturbed.plot <- get_perturbed_trait_datatable_GP(unscaled.pred.data.df, unscaled.fit.extracted, pert.trait)
  b.data.for.perturbed.plot <- backtransform.data(data.for.perturbed.plot)
}


# Plot the results --------

# setup out sub-directories
curr.fold.graph.dir <- paste0(graph.dir, 'k', '-averaged/')
dir.create(curr.fold.graph.dir)
curr.pert.plt.dir <- paste0(curr.fold.graph.dir, pert.trait, '/') # plots of effects of perturbations  
dir.create(curr.pert.plt.dir)


# plots using the non-backtransformed data ------------------------------------------------------------------------

# plot predictions averaged across all k-fold dags agains real data ----
pred.traits <- get_predicted_traits(k_model)
p <- make_posterior_scatter_plots(unscaled.fit.extracted, unscaled.tot.stan.data, all.pred.traits)
ggsave(paste0(curr.pert.plt.dir, 'avg_pred_vs_measured.pdf'), p, width=15, height=15)
# write the datatable of the info in the 'avg_pred_vs_measured.pdf' plot
avg_pred_vs_measured_dt <- convert_to_data_and_predictions_dt(unscaled.fit.extracted, 
                                                              unscaled.tot.stan.data,
                                                              all.pred.traits, 
                                                              D.data2.cultivars) 
saveRDS(avg_pred_vs_measured_dt, file=paste0(curr.pert.plt.dir, 'avg_pred_vs_measured.data.rds'))

# stan_data4plot <- df2stan_data(subset(D.data2.cultivars, select=data.columns)) # use the original unscaled data to make this.
# p <- make_pairwise_posterior_plots_compact(unscaled.fit.extracted, stan_data4plot, pred.traits, k_model)
# ggsave(paste0(curr.pert.plt.dir, 'pairwise_predictions.png'), p, width=40, height=40)
# p <- make_pairwise_posterior_plots_real_x_compact(unscaled.fit.extracted, stan_data4plot, pred.traits, k_model) # WORKS
# ggsave(paste0(curr.pert.plt.dir, 'pairwise_predictions_experimental_x.png'), p, width=40, height=40)

# plot predicted effect of perturbations on yield -----------
if (length(all.pert.pred.traits) > 0) {
  data.cols <- names(data.for.perturbed.plot)[names(data.for.perturbed.plot) != 'cultivar']
  p <- plot_perturbation_results(subset(data.for.perturbed.plot, select=data.cols), pert.trait, line.plots=TRUE)
  ggsave(paste0(pert.graph.dir, pert.trait, '_', 'k-avg', '_predicted_perturbation_effect.pdf'), p, width=15, height=15)
  gc()
}


# plot the same stuff, but after backtransformation to original trait observation scales is applied ---------------------------

# plot and table of predicted vs observed values ----------
pred.traits <- get_predicted_traits(k_model)
p <- make_posterior_scatter_plots(b.unscaled.fit.extracted, b.unscaled.tot.stan.data, all.pred.traits)
ggsave(paste0(curr.pert.plt.dir, 'avg_pred_vs_measured_back.pdf'), p, width=15, height=15)
avg_pred_vs_measured_dt <- convert_to_data_and_predictions_dt(b.unscaled.fit.extracted, 
                                                              b.unscaled.tot.stan.data,
                                                              all.pred.traits, 
                                                              b.D.data2.cultivars) 
saveRDS(avg_pred_vs_measured_dt, file=paste0(curr.pert.plt.dir, 'avg_pred_vs_measured.data_back.rds'))

# stan_data4plot <- df2stan_data(subset(b.D.data2.cultivars, select=data.columns)) # use the original unscaled data to make this.
# p <- make_pairwise_posterior_plots_compact(b.unscaled.fit.extracted, stan_data4plot, pred.traits, k_model)
# ggsave(paste0(curr.pert.plt.dir, 'pairwise_predictions_back.png'), p, width=40, height=40)
# p <- make_pairwise_posterior_plots_real_x_compact(b.unscaled.fit.extracted, stan_data4plot, pred.traits, k_model) # WORKS
# ggsave(paste0(curr.pert.plt.dir, 'pairwise_predictions_experimental_x_back.png'), p, width=40, height=40)
# gc()

# plot predicted effects of trait perturbations on yield (with backtransformed scales)
if (length(all.pert.pred.traits) > 0) {
  data.cols <- names(data.for.perturbed.plot)[names(data.for.perturbed.plot) != 'cultivar']
  p <- plot_perturbation_results(subset(b.data.for.perturbed.plot, select=data.cols), pert.trait, line.plots=TRUE)
  ggsave(paste0(pert.graph.dir, pert.trait, '_', 'k-avg', '_predicted_perturbation_effect_back.pdf'), p, width=15, height=15)
  gc()
}

