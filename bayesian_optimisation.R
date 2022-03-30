rm(list=ls())
library(here) # here_1.0.1 
library(data.table) # data.table_1.13.2
library(rstan) # rstan_2.21.2

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source('BO_functions.R')


# set RNG seed
set.seed(123)

#### setup paramters ---------------------------------------------------------
data.dir <- here('data', 'BO_data') #../final_traits_spring_GP_raw_imputed_v4_v2/'

kernel_function_file <- './GP_ARD_functions.stan' # location of Stan kernel function source
output.dir <- here('BO_output')
out.suffix <- '_PCA_4PC_constrained_0.001e' # suffix for output files

ecotype <- 'spring' # the ecotype data and model to be used ("spring" / "winter")
do.PCA <- TRUE # whether to transform basis to PCA prior to doing any model fitting, in order to impose correlation constraints
numPCs <- 4 # number of principle components to use during optimisation
epsilon <- 0.01 # epsilon value to be used in expected improvement calculation


# setup output directory -------
dir.create(output.dir)
model.dir <- paste0(output.dir, '/', ecotype, '_hyperparam_model/')
dir.create(model.dir)


# SPECIFY MODEL TO BE USED ---------------------------------------------
# specify model structure based on the avg DAG structure found
# BUT only interested in yield, and the traits directly connected to yield.
# nb if do PCA to deal with correlation, will overwrite the k_model later
# to use PCs rather than trait names directly.
if (ecotype=='spring') {
  data.path <- paste0(data.dir, '/spring_trait_data.rds')
  # spring model
  indep_traits = c('TGW', 'totSeedNum', 'totSeedArea','tenPodWeight','numPodsMain','numPods2ary')
  k_model = '[TGW][totSeedNum][totSeedArea][tenPodWeight][numPodsMain][numPods2ary][totSeedW|TGW:totSeedNum:totSeedArea:tenPodWeight:numPodsMain:numPods2ary]'
  ecotype <- 'spring'
  } else {
  if (ecotype=='winter') {
    data.path <- paste0(data.dir, '/winter_trait_data.rds')
    # winter model
    indep_traits = c('num2aryBranch','numPodsMain','numPods2ary','totSeedNum','TGW','tenPodWeight','beakLen','ovuleArea')
    k_model = '[num2aryBranch][numPodsMain][numPods2ary][totSeedNum][TGW][tenPodWeight][beakLen][ovuleArea][totSeedW|num2aryBranch:numPodsMain:numPods2ary:totSeedNum:TGW:tenPodWeight:beakLen:ovuleArea]'
  } else {
    stop("'ecotype' doesn't specify spring or winter model and data to be used!")
  }
}


# prepare the data ----------------------------------------------
# load the normalised trait data for the appropriate (spring/winter) cultivars
# produced by "id_trait_structure.R"
D.data2.cultivars <- readRDS(data.path)
D.data2.cultivars[, rep:=seq(1,nrow(.SD)), by=.(cultivar)] # add column to keep track of individual plants

# Scale the data for model fitting
data.columns <- names(D.data2.cultivars)[!(names(D.data2.cultivars) %in% c('cultivar', 'rep'))]
# keep means and sds so can unscale later
data.means <- apply(subset(D.data2.cultivars, select=data.columns), 2, mean)
data.sds <- apply(subset(D.data2.cultivars, select=data.columns), 2, sd)
D.data.sc <- data.table(data.frame(scale(data.frame(subset(D.data2.cultivars, select=data.columns)))))


# convert trait space to independent (not correlated) basis-------------
if (do.PCA) {
  X_orig <- D.data.sc[, ..indep_traits]
  # confirm that is already centred and scaled
  # apply(X, 2, FUN=mean)
  # apply(X, 2, FUN=sd)
  
  # perform the PCA
  Xpca <- prcomp(X_orig, center=F, scale=F)
  write.csv(summary(Xpca)$importance, 
         file=paste0(output.dir, '/', ecotype, out.suffix, '.csv'),
         row.names=T)
  
  # update the data and model to be the PCA data, and use the PC names
  Xprime <- data.frame(Xpca$x[,1:numPCs])
  D.data.sc <- cbind(D.data.sc[, c('totSeedW')], Xprime)
  PCstring1 <- paste0('[', paste(colnames(summary(Xpca)$importance)[1:numPCs], collapse=']['), ']')
  PCstring2 <- paste0('[totSeedW|', paste(colnames(summary(Xpca)$importance)[1:numPCs], collapse=':'), ']')
  k_model <- paste0(PCstring1, PCstring2)
}


# OPTIMIZATION OF HYPERPARAMS ----------------------------------------------------------------

# reformat the data for stan modelling
stan.data <- df2stan_data_GP(D.data.sc, k_model)
    
# where will save stan source file
write.path.opt <- paste0(model.dir, 'opt_model.stan')
# where will save the max likelihood hyperparameters
param.file <- paste0(model.dir, 'opt_hyperparams.csv')

# write the model source file
bnmodel2stanmodel_GP_opt(k_model, stan.data, kernel_function_file, write.path.opt)
# compile the optimization model
gp_opt <- stan_model(file=write.path.opt) #stan_model(file='GP_model_test_opt_ARD.stan')
# optimize the hyperparams
opt_fit <- optimizing(gp_opt, data=stan.data, seed=95848338, hessian=FALSE)
# store the optimal values
write_opt_parameters_to_file(opt_fit, param.file)
stan.data <- add_opt_parameters_to_stan_data(opt_fit, stan.data)
gc()
      

# q-EXPECTED IMPROVEMENT CONDITIONAL UPON OBSERVED DATA & ML HYPERPARAMS ----------------------------------------------

# where to find the Stan model for yield prediction. Nb spring and winter can share
# a model, becasue X dims are passed as paramaters.
pred.model.file <- paste0(data.dir, '/pred_model.stan')

# set up compiled stan fit model object to save compiling in the optimisation loop ----
# setup dummy data for prediction
X_pred = runif(length(indep_traits)) # generate dummy q point so stan data will have 
pred.data.df <- data.frame(t(X_pred))
if (!(do.PCA)) {
  names(pred.data.df) <- indep_traits
} else {
  names(pred.data.df) <- colnames(summary(Xpca)$importance)
}
pred.data.df <- pred.data.df[, colnames(stan.data$X_totSeedW)] # make sure X col order is consistent
stan.data$N_pred <- nrow(pred.data.df)
stan.data$X_pred <- as.matrix(pred.data.df)
# compile the model
fit <- stan(file=pred.model.file, data=stan.data, 
                  warmup = 0, iter = 1, chains = 1, refresh=0, algorithm='Fixed_param')

# carry out q-EI using constant_liar algo from 
# http://www.cs.ubc.ca/labs/beta/EARG/stack/2010_CI_Ginsbourger-ParallelKriging.pdf
# nb expected_improvement() defined in BO_functions.R
q.values <- constant_liar(expected_improvement, stan.data, fit, epsilon,
                          q=2,
                          n_restarts=25)

# convert returned points back from PCs to original basis  ------------------------------------
if (do.PCA) {
  # convert observed points
  stan.data$X_totSeedW <- as.matrix(X_orig)
  # convert suggested q-values
  i=1
  for (i in 1:length(q.values)) {
    q.values[[i]]$par <- q.values[[i]]$par %*% t(Xpca$rotation[, 1:numPCs])
    q.values[[i]]$par <- as.vector(q.values[[i]]$par) # convert from 1d matrix to vector to be consistent with non-pca result
  }
}


# unscale the q.values and stan.data ------------------------------------------------------
# get the means and the sds in the same order the X matrices in the stan data are.
covar.means <- data.means[colnames(stan.data$X_totSeedW)]
covar.sds <- data.sds[colnames(stan.data$X_totSeedW)]
totSeedW.mean <- data.means['totSeedW']
totSeedW.sds <- data.sds['totSeedW']

# get the unscaled q-value points and yields
i=1
for (i in 1:length(q.values)) {
  q.values[[i]]$par.unscaled <- unscale(q.values[[i]]$par, covar.means, covar.sds)
  q.values[[i]]$mu.unscaled <- unscale(q.values[[i]]$mu, totSeedW.mean, totSeedW.sds)
  q.values[[i]]$sigma.unscaled <- unscale(q.values[[i]]$sigma, totSeedW.mean, totSeedW.sds)
}
# unscale the observation data too
stan.data$X_totSeedW.unscaled <- matrix(0, 
                                        nrow=nrow(stan.data$X_totSeedW), 
                                        ncol=ncol(stan.data$X_totSeedW))
stan.data$totSeedW.unscaled <- rep(0, length(stan.data$totSeedW))
for (i in 1:nrow(stan.data$X_totSeedW)) {
  stan.data$X_totSeedW.unscaled[i, ] <- unscale(stan.data$X_totSeedW[i, ],
                                                covar.means, covar.sds)
  stan.data$totSeedW.unscaled[i] <- unscale(stan.data$totSeedW[i], totSeedW.mean, totSeedW.sds)
}

# Save the results -------------------------------------------------------------
saveRDS(stan.data, file=paste0(output.dir, '/', ecotype, '_stan.data', out.suffix, '.RDS'))
saveRDS(q.values, file=paste0(output.dir, '/', ecotype, '_q.values', out.suffix, '.RDS'))

  

