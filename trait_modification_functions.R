get_predicted_traits <- function(bnmodel) {
  # take model string, and return a list of the nodes which have parents
  modelList <- parseModelString(bnmodel)
  predicted.traits <- c()
  for (n in names(modelList)){
    if (!(anyNA(modelList[[n]]))){
      predicted.traits <- c(predicted.traits, n)
    }
  }
  return(predicted.traits)
}

parseModelString <- function(bnmodel) {
  # takes bn style model string, and converts it to a list of the parents of each node.
  childs <- strsplit(bnmodel, split='\\[')[[1]]
  modelList <- list()
  
  for (node in childs) {
    if (node == '') { # skip the first, empty one
      next
    }
    d <- substr(node, start=1, stop=nchar(node)-1) # drop closing bracket
    if (!(grepl('\\|', d))){ # if has no parents
      modelList[[d]] <- NA
    } else { # if has at least one parent
      child <- strsplit(d, split='\\|')[[1]][1]
      parent.string <- strsplit(d, split='\\|')[[1]][2]
      parents <- strsplit(parent.string, ':')[[1]] # names of parents of child as list
      modelList[[child]] <- parents
    }
  }
  return(modelList)
}

get_perturbed_independent_traits <- function(bnmodelString, perturbed.node) {
  # return vector of the names of the nodes which are NOT descended from "node" in the bnmodelString
  # either directly or indirectly
  
  independent.nodes <- c()
  modelList <- parseModelString(bnmodelString)
  
  for (curr.test.node in names(modelList)) {
    # check if perturbed is ancestor of curr.test.node
    if (!(is.node.indirect.parent(modelList, curr.test.node, perturbed.node)) &
        curr.test.node != perturbed.node) {
      independent.nodes <- c(independent.nodes, curr.test.node)
    }
  }
  return(independent.nodes)
}

is.node.indirect.parent <- function(modelList, node, candidate.parent) {
  # return TRUE / FALSE whether candidate.parent is parent of node in modelList DAG
  parents <- modelList[[node]]
  parents <- na.omit(parents)
  while (length(parents)!=0) {
    # pop first parent
    curr.parent <- parents[1]
    if (length(parents)==1){
      parents = c()
    } else {
      parents <- parents[2:length(parents)]
    }
    
    if (curr.parent==candidate.parent) {
      return(TRUE)
    }
    parents <- c(parents, na.omit(modelList[[curr.parent]]))
  }
  return(FALSE) # if get to here, never encountered candidate.parent in parents
}  

get_all_predicted_and_perturbed_traits <- function(results, PERT.TRAIT) {
  ##' return vectors of all the traits which have parents (and so are predictions of a model)
  ##' (all.pred.traits)
  ##' returns vector of traits which are predicted, AND descendents (direct or indirect) of the PERT.TRAIT
  ##' (all.pert.pred.traits)
  ##' in at least one of the k-fold bnlearn model strings in results
  
  i = 1
  all.pred.traits <- c()
  all.pert.pred.traits <- c()
  for (i in seq(from=1, to=5)) {
    # get the traits which are predicted in at least one model
    curr.pred.traits <- get_predicted_traits(results$model.strings[i])
    all.pred.traits <- c(all.pred.traits, curr.pred.traits)
    
    # and get the traits which are predicted under perturbation in at least one model
    perturbed.independent.traits <- get_perturbed_independent_traits(results$model.strings[i], PERT.TRAIT) 
    data.traits <- c(perturbed.independent.traits, PERT.TRAIT)
    curr.pert.pred.traits <- setdiff(curr.pred.traits, data.traits)
    all.pert.pred.traits <- c(all.pert.pred.traits, curr.pert.pred.traits)
  }
  all.pred.traits <- unique(all.pred.traits)
  all.pert.pred.traits <- unique(all.pert.pred.traits)
  
  return(list(all.pred.traits, all.pert.pred.traits))
}

get_perturbed_trait_datatable_GP <- function(perturbed.in.data, extractedfitObj, perturbed.trait.name) {
  # take the stan data, and the predictions for the traits if change the "perturbed" parameter, and combine them
  
  # get the perturbed value data
  #x.data <- data.table(perturbed.in.data[[c(perturbed.trait.name, 'cultivar')]])
  x.data <- subset(perturbed.in.data, select=c(perturbed.trait.name, 'cultivar', 'rep'))
  #names(x.data) <- c(perturbed.trait.name)
  
  # get the predictions in consistent format with the input data.
  pred.D <- extractedfitObj # list of array predictions for each trait in "generated quantities" section.
  pred.traits <- names(pred.D)[grepl('_pred', names(pred.D))] # get the predicted trait names
  pred.traits <- pred.traits[pred.traits != paste0(perturbed.trait.name, '_pred')] # drop the perturbed trait
  p.t <- pred.traits[1]
  predictions.list <- list()
  for (p.t in pred.traits) {
    curr.predictions <- pred.D[[p.t]] # get predictions for trait of interest (num. sampled points x num. datapoints to predict)
    curr.predictions <- data.table(t(curr.predictions)) # get in correct orientation for binding to input data
    
    curr.data <- cbind(x.data, curr.predictions)
    curr.data.m <- melt(curr.data, id.vars=c(perturbed.trait.name, 'cultivar', 'rep'))
    curr.data.m$variable <- NULL
    #names(curr.data.m) <- c(perturbed.trait.name, p.t)
    names(curr.data.m)[names(curr.data.m)=='value'] <- p.t
    curr.perturbed.dt <- curr.data.m[,1:3] # seperate the perturbed trait values, cultivar name, and rep number from the prediction
    curr.pred.dt <- curr.data.m[,4]
    predictions.list <- c(predictions.list, list(curr.pred.dt))
  }
  predictions.dt <- do.call('cbind', predictions.list) # bind all the predictions together...
  predictions.dt <- cbind(curr.perturbed.dt, predictions.dt) # bind the pertubed data and cultivars to the predictions too
  
  return(predictions.dt)
}

make_stan_data_NULL <- function(df) {
  out <- list()
  out[['N']] <- nrow(df)
  for (trait in names(all.data.df)){
    out[[trait]] <- df[[trait]]
  }
  return(out)
}

make_stan_model_NULL <- function(stan_data, write.path) {
  # writes stan model text path to make predictions for the distribution of each trait if 
  # had no parent.
  outStr <- ''
  outStr <- paste0(outStr, make_data_NULL(stan_data))
  outStr <- paste0(outStr, make_parameters_NULL(stan_data))
  outStr <- paste0(outStr, make_model_NULL(stan_data))
  outStr <- paste0(outStr, make_generated_NULL(stan_data))
  
  cat(file=write.path, outStr)
}

make_data_NULL <- function(stan_data) {
  outstr <- 'data {\n'
  outstr <- paste0(outstr, '\tint<lower=0> N;\n')
  for (name in names(stan_data)) {
    if (name != 'N'){
      outstr <- paste0(outstr, '\tvector[N] ', name, ';\n')
    }
  }
  outstr <- paste0(outstr, '}\n\n')
  return(outstr)
}

make_parameters_NULL <- function(stan_data) {
  outStr <- 'parameters {\n'
  for (name in names(stan_data)) {
    if (name != 'N') {
      outStr <- paste0(outStr, '\treal mu_', name, ';\n\treal<lower=0> sigma_', name, ';\n\n')
    }
  }
  outStr <- paste0(outStr, '}\n\n')
  return(outStr)
}

make_model_NULL <- function(stan_data) {
  outStr <- 'model {\n'
  for (name in names(stan_data)) {
    if (name != 'N') {
      outStr <- paste0(outStr, '\t', name, ' ~ normal(mu_', name, ', sigma_', name, ');\n')
    }
  }
  outStr <- paste0(outStr, '}\n\n')
}

make_generated_NULL <- function(stan_data) {
  outStr <- 'generated quantities {\n'
  for (name in names(stan_data)) {
    if (name != 'N') {
      outStr <- paste0(outStr, '\tvector[N] ', name, '_rep;\n')
    }
  }
  outStr <- paste0(outStr, '\n\tfor (n in 1:N) {\n')
  for (name in names(stan_data)) {
    if (name != 'N') {
      outStr <- paste0(outStr, '\t\t', name, '_rep[n] = normal_rng(mu_', name, ', sigma_', name, ');\n')
    }
  }
  outStr <- paste0(outStr, '\t}\n}\n')
}

downsample <- function(extracted.NULL.fitObj, numChains, num.iterations.per.chain) {
  num.samples <- num.iterations.per.chain * numChains
  sample.idx <- sample.int(length(extracted.NULL.fitObj[[1]]), num.samples)
  thinned.fitObj <- list()
  for (n in names(extracted.NULL.fitObj)) {
    orig <- extracted.NULL.fitObj[[n]]
    
    if (!('matrix' %in% class(orig))) { # if 1d array
      thinned <- orig[sample.idx]
    } else { # if 2d array
      thinned <- orig[sample.idx, ]
    }
    thinned.fitObj[[n]] <- thinned
  }
  return(thinned.fitObj)
}

make_stan_prediction_data_df_GP <- function(D.data2, t, num.perturbations) {
  # make new stan data list for GP model, for when trait (t) is perturbed.
  # for each real datapoint, generates num.perturbations number of values for the perturbed parameter:
  # D.data2 is the data (was df2stan_data_GP(D.data2, k_model), was used to fit the original model)
  # bn.model.altered is the perturbed DAG (parents of perturbed trait removed)
  # t is the perturbed trait ***
  # num.perturbations : how many perturbation values to generate 
  # fit : is the previously fitted model (has alpha, rho, eta, and sigma values)
  # numSamples : the number of samples to take for the parameter values (downsample from all the ones used in the fitting), if too big, R crashes.
  
  traits <- names(D.data2)
  N <- nrow(D.data2) # number of original datapoints
  
  # make perturbed data values
  t.range <- D.data2[[t]]
  t.delta <- max(t.range) - min(t.range)
  t.min <- min(t.range) - 0.3*t.delta
  t.max <- max(t.range) + 0.3*t.delta
  t.perturb.range <- seq(from=t.min, to=t.max, length.out=num.perturbations) # for each row, make 10 new datapoints
  t.perturb.values <- rep(t.perturb.range, N) # replicate the perturb.range for each row
  # combine the new perturbed values with the real data
  D.new <- subset(D.data2, select=names(D.data2)[names(D.data2) != t]) # drop t from D.data2
  D.new.repd <- D.new[rep(seq_len(nrow(D.new)), each = num.perturbations), ]
  D.new.repd[[t]] <- t.perturb.values
  # convert to stan data list
  
  return(D.new.repd)
}

combine_pred_data_with_plant_id <- function(pred.data.df, D.data.sc.cultivar) {
  # add cultivar and rep info to pred.data.df
  merge.traits <- names(pred.data.df)[names(pred.data.df)!=pert.trait]
  curr.D.data.sc.cultivar <- D.data.sc.cultivar
  curr.D.data.sc.cultivar[[pert.trait]] <- NULL
  pred.data.df$rowOrder <- 1:nrow(pred.data.df) # add column to reorder after merge to maintain order same as in the X.data
  pred.data.df.cultivars <- merge(pred.data.df, curr.D.data.sc.cultivar, by=merge.traits)
  pred.data.df.cultivars <- pred.data.df.cultivars[order(rowOrder),] # NEED pred.data.df.cultivars to be in same order as pred.data.df to be able to add cultivars to predictions!!!
  pred.data.df.cultivars$rowOrder <- NULL
  
  return(pred.data.df.cultivars)
}

make_stan_data_GP <- function(D.data, pred.data.df, modelString, seperate.predictions) {
  # take dataframe of observed datapoints, dataframe of datapoints want to make 
  # predictions for, and bnlearn style string describing the model. 
  # combine them into a list for passing data to stan model
  
  # D.data : dataframe of observed data for use fitting the GP model
  # pred.data.df : dataframe of datapoints want to make predictions for 
  # (include all the traits, no need to remove the ones want to predict).
  
  # reformat the training data df for the stan model to list
  s_data <- df2stan_data_GP(D.data, modelString) # make the data for fitting the model
  # names(s_data)
  # reformat the prediction data df to stan list
  if (seperate.predictions == F) { # if passing in prediction data as matrices of parents for each trait
    s_data_pred <- df2stan_data_GP(pred.data.df, modelString)
    #names(s_data_pred)
    pattern <- 'X_|^N$' # remove the "y" values from the predicted data 
    # (keep the ones with X_, or which just N)
    s_data_pred <- s_data_pred[grepl(pattern, names(s_data_pred))]
    # append "Pred" to all of the predictor data names
    for (i in seq(1, length(names(s_data_pred)))) {
      names(s_data_pred)[i] <- paste0(names(s_data_pred[i]), 'Pred')
    }
  } else if (seperate.predictions == T) { # if passing in prediction data as vectors for combo in stan
    s_data_pred <- list('NPred'=nrow(pred.data.df))
    for (trait in names(pred.data.df)) {
      s_data_pred[[paste0(trait, '_4pred')]] <- pred.data.df[[trait]]
    }
  }
  # combine training and predicting datapoints into single list.
  all_stan_data <- c(s_data, s_data_pred)
  return(all_stan_data)
}

df2stan_data_GP <- function(D.data, modelString) {
  # uses model structure to convert D.data to correct list() format to pass into Stan
  # for GP model (unlike linear), parents of each node need to be in a matrix! (as well as seperate vectors
  # for each node to be predicted).
  
  data_list <- list()
  
  modelList <- parseModelString(modelString)
  node <- names(modelList)[1]
  for (node in names(modelList)) {
    if (!(anyNA(modelList[[node]]))) { # if node has no parents, won't predict it
      paramName <- node # the parameter name for prediction
      Xname <- paste0('X_',node) # the parameter name for the predictor values 
      Xdim <- paste0('D_', node)
      paramValues <- D.data[[node]]
      X.df <- subset(D.data, select=c(modelList[[node]])) 
      setcolorder(X.df, sort(names(X.df))) # need the columns of the parents of the node in alphabetical order, as need to remake in stan, so need known order
      X.M <- as.matrix(X.df) # stan needs matrix
      data_list[[paramName]] <- paramValues
      data_list[[Xname]] <- X.M
      data_list[[Xdim]] <- ncol(X.df) # the number of parents of the current node
    }
  }
  # the number of observations (should be the same for all traits in model, so just use the last one)
  data_list[['N']] <- nrow(X.df)
  
  return(data_list)
}

# functions for writing Stan hyperparameter optimisation model file via bnmodel2stanmodel_GP_opt() -----
bnmodel2stanmodel_GP_opt <- function(bnmodel.string, stan_data, kernel_function_file, write.path) {
  ##' write the source code for the stan model used for hyperparameter optimisation
  
  function.section <- make_function_string(kernel_function_file)
  #cat(function.section)
  data.section <- make_data_string_GP_opt(stan_data)
  #cat(data.section)
  parameters.section <- make_parameter_string_GP_opt(bnmodel.string)
  #cat(parameters.section)
  model.section <- make_model_string_GP_opt(bnmodel.string)
  #cat(model.section)
  
  stan.model <- paste0(function.section, data.section, parameters.section,
                       model.section)
  
  # write to file
  cat(stan.model, file=write.path)
}

make_function_string <- function(kernel_function_file){
  functionString <- paste0(readLines(kernel_function_file), collapse="\n")
  functionString <- paste0(functionString, '\n\n')
  return(functionString)
}

make_data_string_GP_opt <- function(stan_data) {
  # make the data section for the stan file for parameter optimization for Maximum Marginal Likelihood
  # approach to hyperparameters
  #print('generating data section...')
  data_string <- "data{\n"
  data_string <- paste0(data_string, stan_indent(1), "// Define variables in data\n\t//X_node = parent data of node\n\t//D_node = num parents of node\n")
  data_string <- paste0(data_string, stan_data_declare_GP_opt(stan_data))
  data_string <- paste0(data_string, "}\n\n")
  return(data_string)
}

stan_indent <- function(n) {
  return(strrep("\t", n))
}

stan_data_declare_GP_opt <- function(stan_data) {
  # declare data for optimization file for GP Max Marginal Likelihood estimation of hyperparameters
  
  t <- names(stan_data)[1]
  paramString <- ""
  #t <- 'b_TSW_NFB'
  
  # want to order by length to ensure N, and N_samples written first
  trait.lengths <- c()
  trait.names <- c()
  for (t in names(stan_data)) {
    #print(t)
    trait.names <- c(trait.names, t)
    trait.lengths <- c(trait.lengths, length(stan_data[[t]]))
  }
  len.order <- order(trait.lengths)
  ordered.traits <- trait.names[len.order]
  
  for (t in ordered.traits) {
    # if is N / P (number of datapoints)
    #print(t)
    if (length(stan_data[[t]]) == 1) {
      if (t != 'NPred') { # don't want "NPred" param for this
        currString <- paste0("\tint<lower=1> ", t, ";\n") # if is N / D / N_samples
      } else { # skip if NPred
        next
      }
    } else if ('matrix' %in% class(stan_data[[t]])) {
      if (nrow(stan_data[[t]]==stan_data[['N']])) { # if is a predictor matrix (X_)
        node <- strsplit(t, '_')[[1]][2]
        node_dims <- paste0('D_',node)
        currString <- paste0('\tvector[',node_dims, '] ', t, '[N];\n') #    vector[D] X[N];
      } else if (nrow(stan_data[[t]]==stan_data[['NPred']])) { # if is a predictor matrix (X_) for perturbation
        node <- strsplit(t, '_')[[1]][2]
        node_dims <- paste('D_node')
        currString <- paste0('\tvector[',node_dims, '] ', t, '[N_samples];\n')
      }
    } else if (length(stan_data[[t]]) == stan_data[['N']]) { # if is a predicted variable
      currString <- paste0("\tvector[N] ", t, ";\n")
    } else if (length(stan_data[[t]]) == stan_data[['NPred']]) { # if is a predicted, perturbed variable
      next # skip if 4_pred
      #currString <- paste0("\tvector[N_samples] ", t, ";\n")
    } 

    paramString <- paste0(paramString, currString)
  }
  return(paramString)
}

make_parameter_string_GP_opt <- function(bnmodel.string) {
  modelList <- parseModelString(bnmodel.string)
  
  #print('generating parameters section...')
  paramString <- "parameters {\n"
  paramString <- paste0(paramString, stan_indent(1), "// Define parameters\n")
  
  for (node in names(modelList)) {
    if (!(anyNA(modelList[[node]]))) { # if parents exist. any.na equivalent to is.na, here, but without the length warning
      curr.str <- make_a_GP_opt_parameter_string(node, modelList)
      paramString <- paste0(paramString, curr.str, '\n')
    }
  }
  paramString <- paste0(paramString, '}\n\n')
  return(paramString)
}

make_a_GP_opt_parameter_string <- function(node, modelList) {
  outStr <- paste0('\tvector<lower=0>[D_',node,'] Rho_', node, '; //for automatic relevance determination, each parent gets own Rho\n') 
  outStr <- paste0(outStr, '\treal<lower=0> Alpha_', node, ';\n')
  outStr <- paste0(outStr, '\treal<lower=0> Sigma_', node, ';\n')
  #outStr <- paste0(outStr, '\tvector[N] Rho_', node, ';\n')
  return(outStr)
}

make_model_string_GP_opt <- function(bnmodel.string) {
  # make model section for Max Marginal Likelihood Estimation of hyperparameters
  pred.traits <- get_predicted_traits(bnmodel.string)
  outStr <- 'model {\n\t// Calculate Kernels\n'
  
  # make the kernel section
  n <- 'A'
  for (n in pred.traits) { # for each node with parents
    curr.str <- paste0('\tmatrix[N, N] K_', n, ' = rbf(X_', n, ', X_', n, ', Alpha_', n , ', Rho_', n, ')\n')
    curr.str <- paste0(curr.str, '\t\t\t + diag_matrix(rep_vector(square(Sigma_',n, '), N));\n')
    curr.str <- paste0(curr.str,  '\tmatrix[N, N] K_L_', n, ' = cholesky_decompose(K_', n, ');\n')
    
    outStr <- paste0(outStr, curr.str)
  }
  
  # make the priors section
  outStr <- paste0(outStr, '\n\t// Priors act as regularisation\n')
  for (n in pred.traits) {
    curr.str <- paste0('\tAlpha_', n, ' ~ normal(0,1);\n')
    curr.str <- paste0(curr.str, '\tRho_', n, ' ~ inv_gamma(5, 5);\n')
    curr.str <- paste0(curr.str, '\tSigma_', n, ' ~ normal(0, 1);\n')
    curr.str <- paste0(curr.str, '\n')
    
    outStr <- paste0(outStr, curr.str)
  }
  
  # make the likelihood section
  outStr <- paste0(outStr, '\n\t// Likelihood section\n')
  for (n in pred.traits) {
    outStr <- paste0(outStr, '\t', n, ' ~ multi_normal_cholesky(rep_vector(0,N), K_L_', n, ');\n')
  }
  outStr <- paste0(outStr, '}\n\n')
  
  return(outStr)
}

add_opt_parameters_to_stan_data <- function(opt_fit, tot.stan.data) {
  # add Max. Marginal Likelihood hyperparamters (alpha, rho, sigma) to stan data
  # to pass into the predictive model.
  
  n <- 'Rho_podLen[1]'
  for (n in names(opt_fit$par)) {
    if (!(grepl('Rho', n))) {
      p.name <- n
      p.val <- opt_fit$par[[n]]
    } else { # if is Rho_<trait>[i]
      p.name <- strsplit(n, '\\[')[[1]][1]
      p.val <- opt_fit$par[grepl(p.name, names(opt_fit$par))] # relies on being ordered rho[1] rho[2] etc in opt_fit
      p.val <- as.array(p.val)
    }
    print(paste0(p.name, ' : ', p.val))
    tot.stan.data[[p.name]] = p.val
  }
  return(tot.stan.data)
}

write_opt_parameters_to_file <- function(opt_fit, write.path) {
  # write max. marginal likelihood hyperparameters (alpha, rho, sigam) to write.path file
  out.str <- ''
  n <- 'Rho_podLen[1]'
  for (n in names(opt_fit$par)) {
    p.name <- n
    p.val <- opt_fit$par[[n]]
    out.str <- paste0(out.str, p.name, '\t', p.val, '\n')
  }
  
  cat(out.str, file=write.path)
  return()
}

# functions for writing Stan perturbation prediction model file via bnmodel2stanmodel_GP_pred_v2() -----
bnmodel2stanmodel_GP_pred_v2 <- function(bnmodel.string, perturbed.trait, stan_data, kernel_function_file, write.path) {
  # as bnmodel2stanmodel_GP_pred, except makes _rep and _pred predictions at the same time to try and avoid the weird problem 
  # where _rep is unrepresentative of _pred 
  
  print('***generating .stan file using bnmodel2stanmodel_GP_pred_v2() ***')
  function.section <- make_function_string(kernel_function_file)
  #cat(function.section)
  data.section <- make_data_string_GP_pred(stan_data)
  #cat(data.section)
  transformed.data.section <- make_transformed_data_string_GP()
  #cat(transformed.data.section)
  parameters.section <- 'parameters{\n}\n\n'
  #cat(parameters.section)
  model.section <- 'model{\n}\n\n'
  #cat(model.section)
  gen.section <- make_generative_string_GP_pred_v2(bnmodel.string, perturbed.trait, stan_data)
  #cat(gen.section)
  
  stan.model <- paste0(function.section, data.section, transformed.data.section, parameters.section,
                       model.section, gen.section)
  #cat(stan.model)
  # write to file
  cat(stan.model, file=write.path)
}

make_data_string_GP_pred <- function(stan_data) {
  #print('generating data section...')
  data_string <- "data{\n"
  data_string <- paste0(data_string, stan_data_declare_GP_pred(stan_data))
  data_string <- paste0(data_string, '}\n\n')
  #cat(data_string)
  return(data_string)
}

stan_data_declare_GP_pred <- function(stan_data) {
  t <- names(stan_data)[1]
  paramString <- ""
  #t <- 'b_TSW_NFB'
  
  # want to order by length to ensure N, and N_samples written first
  trait.lengths <- c()
  trait.names <- c()
  for (t in names(stan_data)) {
    #print(t)
    trait.names <- c(trait.names, t)
    trait.lengths <- c(trait.lengths, length(stan_data[[t]]))
  }
  len.order <- order(trait.lengths)
  ordered.traits <- trait.names[len.order]
  
  #t <- 'Rho_tenPodNumber'
  for (t in ordered.traits) {
    # if is N / P (number of datapoints)
    #print(t)
    if (grepl('Alpha_|Sigma_|Rho_', t)) { # if is a parameter
      if (grepl('Alpha_|Sigma_', t)) { # if is a real number
        currString <- paste0('\treal<lower=0> ', t, ';\n')
      } else if (grepl('Rho_', t)) {
        D.t <- paste0('D_', strsplit(t, '_')[[1]][2])
        currString <- paste0('\tvector<lower=0>[', D.t, '] ', t, ';\n')
      }
    } else if (length(stan_data[[t]]) == 1) {
      currString <- paste0("\tint<lower=1> ", t, ";\n") # if is N / D / N_samples
    } else if ('matrix' %in% class(stan_data[[t]])) {
      if (nrow(stan_data[[t]]==stan_data[['N']])) { # if is a predictor matrix (X_)
        node <- strsplit(t, '_')[[1]][2]
        node_dims <- paste0('D_',node)
        currString <- paste0('\tvector[',node_dims, '] ', t, '[N];\n') #    vector[D] X[N];
      } else if (nrow(stan_data[[t]]==stan_data[['NPred']])) { # if is a predictor matrix (X_) for perturbation
        node <- strsplit(t, '_')[[1]][2]
        node_dims <- paste('D_node')
        currString <- paste0('\tvector[',node_dims, '] ', t, '[N_samples];\n')
      }
    } else if (length(stan_data[[t]]) == stan_data[['N']]) { # if is a predicted variable
      currString <- paste0("\tvector[N] ", t, ";\n")
    } else if (length(stan_data[[t]]) == stan_data[['NPred']]) { # if is a predicted, perturbed variable
      currString <- paste0("\tvector[NPred] ", t, ";\n")
    } 
    # for the X matrices, need to get the correct D dimension name
    else if (length(stan_data[[t]]) == stan_data[['P']]) {
      currString <- paste0("\tvector[P] ", t, ";\n")
    } 
    paramString <- paste0(paramString, currString)
  }
  return(paramString)
}

make_transformed_data_string_GP <- function(){
  s <- 'transformed data {\n\treal delta = 1e-9; // needed to ensure K positive semidef\n}\n\n'
  return(s)
}

make_generative_string_GP_pred_v2 <- function(bnmodel.string, perturbed.trait, stan_data) {
  # make generative section for making _rep and _pred from the same model fitting function call 
  nodes <- get_predicted_traits(bnmodel.string)
  
  outStr <- 'generated quantities {\n'
  
  # declare _rep (predictions for observed data) outputs. (including for the perturbed trait)
  declareStr_rep <- ''
  for (n in nodes) {
    declareStr_rep <- paste0(declareStr_rep, '\tvector[N] ', n, '_rep;\n')
  }
  outStr <- paste0(outStr, declareStr_rep, '\n')
  
  # declare _pred (predictions for unobserved/perturbed data) outputs
  # don't predict the perturbed trait, or anything which is not a descendent of it!
  data_4pred_traits <- names(stan_data)[grepl('_4pred', names(stan_data))]
  declareStr_pred <- ''
  n <- 'C'
  for (n in nodes[nodes != perturbed.trait]) {
    if (!(any(grepl(paste0('^', n, '_'), data_4pred_traits)))) { # if n is not one of the _4pred traits, try and predict it
      declareStr_pred <- paste0(declareStr_pred, '\tvector[NPred] ', n, '_pred;\n')
    }
  }
  outStr <- paste0(outStr, declareStr_pred, '\n\n')
  #cat(outStr)
  
  # make the predictions for the _rep data points
  outStr <- paste0(outStr, '\t// Predictions for the observed data\n')
  rep_str <- make_rep_predictive_str_v2(bnmodel.string, data_4pred_traits)
  outStr <- paste0(outStr, rep_str, '\n\n')
  
  # make the predictions for the _pred (unobserved) data, using other predictions
  # need to get order to predict in (dependencies). 
  pred_str <- make_pred_predictive_str_v2(bnmodel.string, perturbed.trait, stan_data)
  outStr <- paste0(outStr, pred_str)
  
  outStr <- paste0(outStr, '}\n\n')
  
  return(outStr)
}

make_rep_predictive_str_v2 <- function(bnmodel.string, non_pred_traits) {
  # make the '_rep' section of the generative data section of the GP predictive model in stan
  # in v2, don't do this for the traits will also predict _pred for (only do it for the non-root ones which
  # are also not _pred predicted, i.e. for roots in non_pred_traits)
  
  
  
  predtraits <- get_predicted_traits(bnmodel.string) # the non-root nodes, which will want a _rep result for
  
  # make the set of traits need to declare and define _rep seperately to _pred bit
  repTraits <- c()
  for (t in predtraits) {
    if ((any(grepl(paste0('^', t, '_'), non_pred_traits)))) {
      repTraits <- c(repTraits, t)
    }
  }
  
  outStr <- ''
  
  # if any rep traits which need defining in this section
  if (length(repTraits != 0)) {
    
    outStr <- paste0(outStr, '\t{\n')
    # declare the f(x) intermediate output
    for (t in repTraits) {
      outStr <- paste0(outStr, '\t\tvector[N] f_', t, '_rep;\n')
    }
    outStr <- paste0(outStr, '\n')
    
    # predict the f(x) values
    for (t in repTraits) {
      curr.str <- paste0('\t\tf_',t,'_rep = gp_pred_rng(X_', t, ', ', t, ', X_', t, ',\n\t\t\t\t\t',
                         'Alpha_', t, ', Rho_', t, ', Sigma_', t, ', delta);\n')
      outStr <- paste0(outStr, curr.str)
    }
    
    # predict the y(x) values
    outStr <- paste0(outStr,'\n\t\tfor (n in 1:N) {\n')
    for (t in repTraits) {
      curr.str <- paste0('\t\t\t', t, '_rep[n] = normal_rng(f_', t,'_rep[n], Sigma_', t, ');\n')
      outStr <- paste0(outStr, curr.str)
    }
    outStr <- paste0(outStr, '\t\t}\n')
    outStr <- paste0(outStr, '\t}\n')
  }
  
  return(outStr)
}

make_pred_predictive_str_v2 <- function(bnmodel.string, perturbed.trait, stan_data) {
  
  # get dependency order of traits will predict values for
  outStr <- '\t//Predictions for the perturbed data\n'
  ordered_traits <- get_dependency_order(bnmodel.string) # all the traits in dependency order
  
  # only want to predict the ones which descendents of the perturbed trait
  data_4pred_traits <- names(stan_data)[grepl('_4pred', names(stan_data))] # the traits which passed in data for predicting perturbation for
  data_4pred_trait_names <- tstrsplit(data_4pred_traits, '_4')[[1]]
  ordered_pred_traits <- ordered_traits[!(ordered_traits %in% data_4pred_trait_names)]
  
  pred.trait <- ordered_pred_traits[1]
  for (pred.trait in ordered_pred_traits) {
    # make the string defining how to predict the current "pred.trait"
    curr.pred.string <- make_curr_pred_trait_v2(bnmodel.string, pred.trait, perturbed.trait, data_4pred_trait_names)
    outStr <- paste0(outStr, curr.pred.string)
  }
  #outStr <- paste0(outStr, '}\n\n')
  return(outStr)
}

get_dependency_order <- function(bnmodel.string) {
  # return list of predicted traits in the order (first to last), so roots are first, 
  #then stuff that only depends on them, then stuff that depends on that layer etc.
  # requires that model be DAG!! (but so does everything else :) )
  
  model.list <- parseModelString(bnmodel.string)
  #print(model.list)
  num.pred.traits <- length(model.list)
  inf.counter <- 0
  
  ordered.traits <- c()
  updated.model.list <- model.list
  while (num.pred.traits > length(ordered.traits)) { # while still not all resolved
    # sanity check not Inf loop
    inf.counter <- inf.counter + 1
    if (inf.counter > 100) {
      stop("get_dependency_order(): more than 100 iterations : check that specified modelString is a DAG!")
    }
    
    #print(ordered.traits)
    #print(model.list)
    model.list <- updated.model.list
    n <- names(model.list)[1]
    for (n in names(model.list)) {
      if (anyNA(model.list[[n]])) { # if has no unresolved dependencies
        ordered.traits <- c(ordered.traits, n) # append to output
        updated.model.list[[n]] <- NULL # delete from modellist
        # remove from other nodes' dependencies
        other.node <- names(updated.model.list)[2]
        for (other.node in names(updated.model.list)) {
          if (!(anyNA(updated.model.list[[other.node]]))) { # if other node has dependencies still
            # remove node (n) from any depenencies
            updated.model.list[[other.node]] <- updated.model.list[[other.node]][updated.model.list[[other.node]] != n]
            if (length(updated.model.list[[other.node]])==0) { # if empty, set to NA
              updated.model.list[[other.node]] <- NA
            }
          }
        }
      }
    }
  }
  return(ordered.traits)
}

make_curr_pred_trait_v2 <- function(bnmodel.string, n, perturbed.trait, data_traits) {
  # for a single trait (n), make the string for how to 
  # construct the necessary X to predict it etc.
  
  #root.traits <- get_root_traits(bnmodel.string) # if root, will use "<trait>_4pred", 
  # otherwise, will use <triat_pred" to fill in X_<n>
  
  
  
  # get parents of trait n
  model.list <- parseModelString(bnmodel.string)
  ordered.parents <- sort(model.list[[n]]) # must be alphabetical to be consistent with the X_<n>
  # in the real data
  
  outStr <- paste0('\t{\n\t\t// ', n, '\n')
  outStr <- paste0(outStr, '\t\tvector[NPred+N] f_', n, '_pred;\n') # f_ for current trait
  outStr <- paste0(outStr, '\t\tvector[NPred+N] p_', n, '_pred;\n') # all the predictions (_rep and _pred) for the current trait
  outStr <- paste0(outStr, '\t\tvector[D_', n, '] X_', n, 'Pred[NPred]; // declare X array for prediction\n') # the X matrix for the _pred values
  outStr <- paste0(outStr, '\t\tvector[D_', n, '] X_', n, 'All[NPred+N]; // declare X array for prediction (_pred and _rep)\n') # the X matrix for the _pred values and the _rep values
  
  outStr <- paste0(outStr, '\t\t// fill in X array\n')
  outStr <- paste0(outStr, '\t\tfor (n in 1:NPred) {\n')
  # loop to fill in X_<trait>
  for (i in seq(1, length(ordered.parents))){
    cp <- ordered.parents[i]
    if (cp %in% data.traits | cp == perturbed.trait) { # if is root, or is the perturbed trait. Then use the data values fed into the model
      currStr <- paste0('\t\t\tX_', n, 'Pred[n, ', i, '] = ', cp, '_4pred[n]; //parents ordered alphabetically\n')
    } else { # if current parent is not root trait, or the perturbed trait, use its model predicted value
      currStr <- paste0('\t\t\tX_', n, 'Pred[n, ', i, '] = ', cp, '_pred[n]; //parents ordered alphabetically\n')
    }
    outStr <- paste0(outStr, currStr)
  }
  outStr <- paste0(outStr, '\t\t}\n\n') # close the X_<n> filling in loop
  
  # concat the X for _rep to the start of it
  outStr <- paste0(outStr, '\t\tX_', n, 'All = append_array(X_', n, ', X_', n, 'Pred);\n\n')
  
  outStr <- paste0(outStr, '\t\tf_', n, '_pred = gp_pred_rng(X_', n, 'All', ',',n,', X_', n,',\n')
  outStr <- paste0(outStr, '\t\t\t\t\tAlpha_', n, ', Rho_', n, ', Sigma_', n, ', delta);\n\n')
  outStr <- paste0(outStr, '\t\tfor (n in 1:N+NPred) {\n')
  outStr <- paste0(outStr, '\t\t\tp_', n, '_pred[n] = normal_rng(f_', n, '_pred[n], Sigma_', n, ');\n')
  outStr <- paste0(outStr, '\t\t}\n\n')
  
  # split the _rep and the _pred output
  outStr <- paste0(outStr, '\t\t//Split the _rep and _pred predictions\n')
  outStr <- paste0(outStr, '\t\t', n, '_rep = p_', n, '_pred[1:N];\n')
  outStr <- paste0(outStr, '\t\t', n, '_pred = p_', n, '_pred[N+1:N+NPred];\n')
  
  outStr <- paste0(outStr, '\t}\n\n') # close the current section
  
  return(outStr)
}

write_convergence_check <- function(fitObj, write.file) {
  # print summary table (with Rhat values) to file. (for checking that chains well mixed, and neff
  # not too small)
  S <- data.table(data.frame(summary(fit)[['summary']]))
  names(S) <- c('mean', 'se_mean', 'sd', '2.5_perc', '25_perc', '50_perc', '75_perc', '97.5_perc', 'n_eff', 'Rhat')
  S <- data.frame(summary(fit)[['summary']])
  fwrite(S, write.file, row.names=TRUE) 
}

# Unscaling data functions -----
unscale.extracted.fit <- function(extractedfitObj, modelString, data.means, data.sds) {
  traits <- names(parseModelString(modelString))
  out.list <- list()
  for (t in traits) {
    for (i in names(extractedfitObj)) {
      if (grepl(t, i)) {
        transformed.fit <- (extractedfitObj[[i]] * data.sds[t]) + data.means[t]
        out.list[[i]] <- transformed.fit
      }
    }
  }
  out.list[['lp__']] <- extractedfitObj[['lp__']]
  return(out.list)
}

unscale.stan.data <- function(tot.stan.data, modelString, data.means, data.sds) {
  # reverse the "scaling" of the data performed prior to fitting and estimation
  # leave the hyperparameters alone, as rescaling seems misleading / not sure 
  # how would work
  traits <- names(parseModelString(modelString))
  t='podLen'
  out.stan.data <- list()
  for (t in traits) {
    for (i in names(tot.stan.data)) {
      if (grepl(t, i) & !(grepl('D_|Alpha_|Rho_|Sigma_', i))) { # if the tot.stan.data relates to the trait
        transformed.data <- (tot.stan.data[[i]] * data.sds[t]) + data.means[t]
        out.stan.data[[i]] <- transformed.data
      }
    }
  }
  return(out.stan.data)
}

unscale.data.frame <- function(scaled.data, data.means, data.sds) {
  # takes dataframe of scaled data, and unscales it using stored mean and standard deviation
  out.df <- scaled.data
  for (c in names(data.means)) {
    out.df[[c]] <- (out.df[[c]]*data.sds[c]) + data.means[c]
  }
  return(out.df)
}

# backtransform data functions -----
backtransform.data <- function(IN) {
  log.traits <- c('num2aryBranch', 'numFBranch', 'numFlowers',
                  'ovaryLen', 'tenPodArea', 'totSeedArea', 
                  'TGW')
  
  sqrt.traits <- c('height', 'numPods2ary', 'numPodsTot')
  
  perc.logit.traits <- c('percAborted2ary', 'percAbortedMain', 'percAbortedTot')
  
  logit.traits <- c('tenPodCompactness', 'totCompactness')
  
  inv.traits <- c('oilContent')
  
  #t <- names(IN)[1]
  #t <- 'X_TGW'
  for (t in names(IN)) {
    trait <- strsplit(t, '_')[[1]][1]
    if (trait == 'X') { # format of some is X_traitname, 
      trait <- strsplit(t, '_')[[1]][2]
    }
    
    if (trait %in% log.traits) {
      inv.trait <- inv.log(IN[[t]])
    } else if (trait %in% perc.logit.traits) {
      inv.trait <- inv.perc.logit(IN[[t]])
    } else if (trait %in% logit.traits) {
      inv.trait <- inv.logit(IN[[t]])
    } else if (trait %in% sqrt.traits) {
      inv.trait <- inv.sqrt(IN[[t]])
    } else if (trait %in% inv.traits) {
      inv.trait <- inv.inv(IN[[t]])
    }else {
      inv.trait <- IN[[t]]
    }
    
    IN[[t]] <- inv.trait
  }
  return(IN)
}

inv.log <- function(v) {
  exp(v)
}
inv.sqrt <- function(v) {
  v**2
}
inv.perc.logit <- function(v) {
  (exp(v) / (1+exp(v)))*100
}
inv.logit <- function(v) {
  (exp(v) / (1+exp(v)))
}
inv.inv <- function(v) {
  50 * exp(v) / (1 + exp(v))
}

# plotting & output functions -----
make_posterior_scatter_plots <- function(extractedfitObj, fitData, traitNameList) {
  # use bayesplot to
  # make a multiplot cowplot of the mean posterior value against 
  # the true value for each trait in traitNameList. 
  # requires that the posterior has the same name "_rep" as the trait 
  # data
  t <- 'totSeedW' #traitNameList[1]
  p.list <- list()
  for (t in sort(traitNameList)) {
    print(t)
    t_real <- fitData[t]
    #t_rep <- as.matrix(fitObj, pars = paste0(t, '_rep'))
    t_rep <- as.matrix(extractedfitObj[[paste0(t, '_rep')]])
    # calculate proportion of variance explained by model
    R <- calc.prop.of.variance(t_real[[1]], t_rep)
    
    p <- ppc_scatter_avg(y=t_real[[1]], yrep=t_rep)
    p <- p + ggtitle(paste0(t, ', ', R, '% explained')) +
      theme(plot.title = element_text(size=2))+
      theme_bw()
    p.list <- c(p.list, list(p))
  }
  ncol <- ceiling(sqrt(length(traitNameList)))
  nrow <- ceiling(length(traitNameList) / ncol)
  p.grid <- plot_grid(plotlist=p.list, ncol=ncol, nrow=nrow)
  
  # add common title
  title <- ggdraw()+draw_label(
    "Measured (y) vs mean predicted (Average y_rep) trait values",
    fontface='bold',
    x=0,
    hjust = 0)+
    theme(plot.margin = margin(0, 0, 0, 7))
  p.grid <- plot_grid(title, p.grid, ncol=1, rel_heights=c(0.1, 1))
  #p.grid
  return(p.grid)
}

calc.prop.of.variance <- function(y_obs, y_pred_matrix) {
  # see http://onlinestatbook.com/2/regression/partitioning.html
  # calculates proportion of variance explained as 1 - (SSY / SSE)
  # SSY = sum((y - y_mean)**2)
  # SSY = SSY' + SSE
  # SSE = sum((y - y_pred)**2)
  # SSY' = sum((y_pred - y_mean)**2)
  
  # proportion explained = SSY' / SSY 
  # proportion not explained = SSE / SSY
  
  y_mean <- mean(y_obs)
  SSY <- sum((y_obs - y_mean)**2)
  
  avg_y_pred <- apply(y_pred_matrix, 2, mean)
  SSE <- sum((y_obs - avg_y_pred)**2)
  
  R <- (1 - (SSE / SSY))*100
  return(round(R, digits=2))
}

convert_to_data_and_predictions_dt <- function(extractedfitObj, fitData, traitNameList, D.data2.cultivars) {
  # recombine the original trait data, 
  # the predicted _rep values
  # and the cultivar corresponding to each datapoint into a datatable for writing to .rds object
  # for use in GWAS
  
  # get the cultivar info for the data, relying on the order of the data being unchanged
  out.df <- data.table(data.frame('cultivar'=D.data2.cultivars$cultivar))
  
  t <- 'totSeedNum'
  for (t in sort(traitNameList)) {
    
    t_real <- fitData[[t]] # the unscaled data was fitting model to
    t_original <- D.data2.cultivars[, ..t][[1]] # the original trait data fed into the "simple_yield_cluster_stan.R" script
    
    # check that the order of the data hasn't been screwed up, or that 
    # the unscaling hasn't gone wrong somehow
    stopifnot(all(round(t_real)==round(t_original)))
    
    t_rep <- as.matrix(extractedfitObj[[paste0(t, '_rep')]])
    t_rep_mean <- apply(t_rep, 2, mean)
    
    #plot(t_rep_mean, t_real) # plot is same as the info in avg_pred_vs_measured.pdf (what expect)
    out.df[, eval(t):=t_real]
    out.df[, paste0(eval(t), '_hat'):= t_rep_mean]
  }
  return(out.df)
}

plot_perturbation_results <- function(perturbed.dt, perturbed.trait.name, line.plots=FALSE) {
  # perturbed.dt is data table of perturbed values fed into model, against model predicted values, using measured values of root nodes as inputs
  
  # get rid of columns not needed for this plot
  if ('cultivar' %in% names(perturbed.dt)) {
    perturbed.dt$cultivar <- NULL
  }
  if ('rep' %in% names(perturbed.dt)) {
    perturbed.dt$rep <- NULL
  }
  perturbed.dt <- unique(perturbed.dt)
  
  num.bins <- length(unique(perturbed.dt[[perturbed.trait.name]])) - 1
  
  dt.m <- melt(perturbed.dt, id.vars=perturbed.trait.name)
  
  dt.m$variable <- factor(dt.m$variable, levels=c(sort(unique(as.character(dt.m$variable)))))
  
  # p <- ggplot(dt.m, aes_string(x=perturbed.trait.name, y='value'))+
  #   geom_bin2d(bins=num.bins)+
  #   scale_fill_continuous(type='viridis')+
  #   facet_wrap(~variable, scales='free')+
  #   theme_bw()
  #seperate plots so can have seperate fills.
  p.list <- list()
  curr.trait <-  sort(unique(as.character(dt.m$variable)))[1]
  for (curr.trait in sort(unique(as.character(dt.m$variable)))) {
    curr.dt <- dt.m[dt.m$variable==curr.trait, ]
    curr.label <- strsplit(curr.trait, '_')[[1]][1]
    
    # make dataframe for calculating medians & quantiles (used in both plots)
    curr.dt.copy <- curr.dt
    
    if (line.plots==TRUE) { # plot quanitle lineplots of CIs of results rather than heatmap
      # calculate quantiles, at each value of the perturbed trait
      curr.dt.copy[, q0:=quantile(value, probs=0.0), by=.(get(perturbed.trait.name))]
      curr.dt.copy[, q10:=quantile(value, probs=0.1), by=.(get(perturbed.trait.name))]
      curr.dt.copy[, q20:=quantile(value, probs=0.2), by=.(get(perturbed.trait.name))]
      curr.dt.copy[, q30:=quantile(value, probs=0.3), by=.(get(perturbed.trait.name))]
      curr.dt.copy[, q40:=quantile(value, probs=0.4), by=.(get(perturbed.trait.name))]
      curr.dt.copy[, q50:=quantile(value, probs=0.5), by=.(get(perturbed.trait.name))]
      curr.dt.copy[, q60:=quantile(value, probs=0.6), by=.(get(perturbed.trait.name))]
      curr.dt.copy[, q70:=quantile(value, probs=0.7), by=.(get(perturbed.trait.name))]
      curr.dt.copy[, q80:=quantile(value, probs=0.8), by=.(get(perturbed.trait.name))]
      curr.dt.copy[, q90:=quantile(value, probs=0.9), by=.(get(perturbed.trait.name))]
      curr.dt.copy[, q100:=quantile(value, probs=1.0), by=.(get(perturbed.trait.name))]
      
      summary.curr.dt <- unique(subset(curr.dt.copy, select=c(perturbed.trait.name, 'q0','q10', 'q20', 'q30',
                                                              'q40', 'q50', 'q60', 'q70', 'q80', 'q90', 'q100') 
      ))
      
      p <- ggplot(curr.dt, aes_string(x=perturbed.trait.name))+
        #geom_bin2d(bins=num.bins)+
        #scale_fill_continuous(type='viridis')+
        # geom_point(data=summary.curr.dt, 
        #            aes_string(x=perturbed.trait.name, y='q50'),
        #            color='grey')+
        geom_line(data=summary.curr.dt, 
                  aes_string(x=perturbed.trait.name, y='q50'),
                  color='grey', size=1)+
        geom_line(data=summary.curr.dt, 
                  aes_string(x=perturbed.trait.name, y='q0'),
                  color='grey', size=0.5)+
        geom_line(data=summary.curr.dt, 
                  aes_string(x=perturbed.trait.name, y='q10'),
                  color='grey', size=0.5)+
        geom_line(data=summary.curr.dt, 
                  aes_string(x=perturbed.trait.name, y='q20'),
                  color='grey', size=0.5)+
        geom_line(data=summary.curr.dt, 
                  aes_string(x=perturbed.trait.name, y='q30'),
                  color='grey', size=0.5)+
        geom_line(data=summary.curr.dt, 
                  aes_string(x=perturbed.trait.name, y='q40'),
                  color='grey', size=0.5)+
        geom_line(data=summary.curr.dt, 
                  aes_string(x=perturbed.trait.name, y='q60'),
                  color='grey', size=0.5)+
        geom_line(data=summary.curr.dt, 
                  aes_string(x=perturbed.trait.name, y='q70'),
                  color='grey', size=0.5)+
        geom_line(data=summary.curr.dt, 
                  aes_string(x=perturbed.trait.name, y='q80'),
                  color='grey', size=0.5)+
        geom_line(data=summary.curr.dt, 
                  aes_string(x=perturbed.trait.name, y='q90'),
                  color='grey', size=0.5)+
        geom_line(data=summary.curr.dt, 
                  aes_string(x=perturbed.trait.name, y='q100'),
                  color='grey', size=0.5)+
        #facet_wrap(~variable, scales='free')+
        theme_bw()+
        theme(legend.position = 'none')+
        ylab(curr.label)+
        ggtitle(curr.label)
      
    } else { # plot heatmap
      curr.dt.copy[, q50:=quantile(value, probs=0.5), by=.(get(perturbed.trait.name))]
      summary.curr.dt <- unique(subset(curr.dt.copy, select=c(perturbed.trait.name, 'q50') ))
      
      p <- ggplot(curr.dt, aes_string(x=perturbed.trait.name, y='value'))+
        geom_bin2d(bins=num.bins)+
        scale_fill_continuous(type='viridis')+
        geom_point(data=summary.curr.dt, 
                   aes_string(x=perturbed.trait.name, y='q50'),
                   color='grey')+
        geom_line(data=summary.curr.dt, 
                  aes_string(x=perturbed.trait.name, y='q50'),
                  color='grey')+
        #facet_wrap(~variable, scales='free')+
        theme_bw()+
        theme(legend.position = 'none')+
        ylab(curr.label)+
        ggtitle(curr.label)
    }
    
    p.list[[curr.trait]] <- p
  }   
  p.grid <- plot_grid(plotlist=p.list)    
  #p.grid
  
  return(p.grid)
}

