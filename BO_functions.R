expected_improvement <- function(X_pred, fit, stan.data, epsilon) {
  # calculate ei at the point X, given observation data (stan.data), and compiled model (fit)
  # following http://krasserm.github.io/2018/03/21/bayesian-optimization/
  # epsilon is exploration exploitation tradeoff. Bigger leads to more exploration
  # ASSUMES ONLY ONE POINT X_PRED BEING CALCULATED!
  

  # add the point to calc ei for (X_pred) points to the stan.data 
  pred.data.df <- data.frame(t(X_pred))
  #names(pred.data.df) <- indep_traits
  #pred.data.df <- pred.data.df[, colnames(stan.data$X_totSeedW)] # make sure X order is consistent
  stan.data$N_pred <- nrow(pred.data.df)
  stan.data$X_pred <- as.matrix(pred.data.df)
  
  # by using a previously fit model object, avoids need to recompile.
  results <- stan(fit=fit, data=stan.data, warmup = 0, 
                  iter = 1, chains = 1, 
                  refresh=0, algorithm='Fixed_param')
  extractedfitObj <- extract(results)
  
  # calculate EI here, as can't figure out the functions stan.
  mu <- extractedfitObj$f_pred[1,1]
  sigma <- sqrt(extractedfitObj$cov_pred[1,1,1])
  mu_sample <- extractedfitObj$y[1,]
  
  # print('epsilon')
  # print(epsilon)
  
  if (sigma == 0) {
    ei = 0.0
  } else {
    mu_sample_opt = max(mu_sample)
    imp = mu - mu_sample_opt - epsilon
    Z = imp / sigma
    ei = imp * pnorm(Z) + sigma * dnorm(Z)
  }
  
  return(ei)
}

propose_location <- function(acquisition_fun, 
                             stan.data,
                             fit, 
                             epsilon,
                             n_restarts=25) {
  # proposes the next sampling point, by optimisation of the
  # acquisition function (e.g. neg_expected_improvement() ).
  # constrained to search within the observed min and max for each dimension.
  
  # returns:
  # -1: no optimisations converged successfully
  # OR
  # list of:
  # $par: X at optimal
  # $value: expected improvement at this point
  # $counts : ??
  # $convergence: 0 indicating convergence
  # $message: optim message
  # $times.found : the number of random restarts which found this solution.
  
  # calculate the lower and upper bounds to be within the observed range
  lower <- apply(stan.data$X_totSeedW, MARGIN=2, FUN=min)
  upper <- apply(stan.data$X_totSeedW, MARGIN=2, FUN=max)
  
  # optim minimizes by default, so multiply by -1.
  min_obj <- function(X) {
    return (-1 * acquisition_fun(X, fit, stan.data, epsilon))
  }
  
  # do the optimisations
  optims = list()
  for (n in 1:n_restarts) {
    print(paste0('random restart: ', n))
    init.vals <- rnorm(stan.data$D_totSeedW)
    
    # given bad starting params, sometimes gives error
    opt <- tryCatch( {
      optim(par=init.vals, min_obj, method='L-BFGS-B', lower=lower, upper=upper)},
      error=function(cond) {
        message(cond)
        return(list('par'=NA, 'value'=NA, 
                    'counts'=NA, 'convergence'=-1,
                    'message'=message(cond)))}
    )
    
    opt$value <- -1*opt$value # convert minimize to maximise
    optims[[n]] <- opt
  }
  
  # get the runs which actually converged
  convs <- list()
  count <- 1
  for (i in 1:n_restarts) {
    if (optims[[i]]$convergence==0) {
      convs[[count]] <- optims[[i]]
      count = count + 1
    }
  }
  
  if (length(convs) == 0) {
    return -1 # no converging optimisations.
  }
  
  # check that optimal value found multiple times - makes more believable
  # that found a decent maxima
  curr.opt.count <- 0
  curr.opt <- -999
  curr.opt.X <- NA
  for (i in 1:length(convs)) {
    # if found a new maximum
    if (convs[[i]]$value > curr.opt) {
      curr.opt <- convs[[i]]$value
      curr.opt.count <- 1
      curr.opt.X <- convs[[i]]$par
      return.idx <- i
    }
    if (isTRUE(all.equal(convs[[i]]$par, curr.opt.X,
                         tolerance = 0.01))) { # if found the optimal location again
      curr.opt.count <- curr.opt.count + 1
    }
  }
  
  # set up return
  convs[[return.idx]]$times.found <- curr.opt.count
  return(convs[[return.idx]])
}

constant_liar <- function(acquisition_fun, 
                          stan.data,
                          fit, 
                          epsilon,
                          q=5,
                          n_restarts=5) {
  # propose next 5 points to sample, following constant liar
  # algo. with constant liar value being max observed y value.
  # http://www.cs.ubc.ca/labs/beta/EARG/stack/2010_CI_Ginsbourger-ParallelKriging.pdf
  # section 4.2.2 
  
  orig.stan.data <- copy(stan.data)
  
  # calculate the constant liar value to use as the max() observed value
  L = max(stan.data$totSeedW)
  
  q.points <- list()
  for (curr_q in 1:q) {
    print(paste0('curr q: ', curr_q))
    
    # get the next proposed location
    next_loc <- propose_location(acquisition_fun, stan.data, fit, epsilon, 
                                 n_restarts=n_restarts)
    
    # update the observations data in line with the constant liar algo
    # as if the proposed value was observed to have y value of L
    stan.data$N <- stan.data$N + 1
    stan.data$totSeedW <- c(stan.data$totSeedW, L)
    stan.data$X_totSeedW <- rbind(stan.data$X_totSeedW, next_loc$par)
    
    q.points[[curr_q]] <- next_loc
  }
  
  # for theproposed points, calculate the mean expected y-value, 
  # and sd
  for (i in 1:length(q.points)) {
    X_pred = q.points[[i]]$par
    orig.stan.data$X_pred <- matrix(X_pred, nrow=1)
    orig.stan.data$N_pred <- 1
    
    results <- stan(fit=fit, data=orig.stan.data, warmup = 0, 
                    iter = 1, chains = 1, 
                    refresh=0, algorithm='Fixed_param')
    extractedfitObj <- extract(results)
    
    # calculate EI here, as can't figure out the functions stan.
    mu <- extractedfitObj$f_pred[1,1]
    sigma <- sqrt(extractedfitObj$cov_pred[1,1,1])
    
    q.points[[i]]$mu <- mu
    q.points[[i]]$sigma <- sigma
  }
  
  return(q.points)
}

unscale <- function(scaled.values, means.vec, sds.vec) {
  # inverse of R scale() function.
  return ((scaled.values * sds.vec) + means.vec)
}


# vector2matrices <- function(networkObj) {
#   # https://github.com/hredestig/pcaMethods/blob/master/R/vector2matrices.R
#   ##' Tranform the vectors of weights to matrix structure
#   ##' 
#   ##' networkObj should be generated as e.g. 
#   ##' P <- pcaMethods::pca(X, method='nlpca', nPcs=ncol(X))
#   ##' networkObj <- P@network
#   
#   netDim <- dim(networkObj@net)
#   posBegin <- 1
#   posEnd <- 0
#   result <- list()
#   i=1
#   for (i in 1:(netDim[2]-1)) {
#     wSize <- c(networkObj@net[i+1], networkObj@net[i]+1) # includes bias
#     posEnd <- posEnd + prod(wSize)
#     result[[i]] <-  
#       matrix(networkObj@weights$current()[posBegin:posEnd], wSize[1], wSize[2])
#     posBegin <- posEnd + 1
#   }
#   
#   if (posEnd > length(networkObj@weights$current())) {
#     stop("weight vector has too many elements\n")
#   }
#   return(result)
# }

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


# functions to write stan optimisation source files -------------------
bnmodel2stanmodel_GP_opt <- function(bnmodel.string, stan_data, kernel_function_file, write.path) {
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
    # for the X matrices, need to get the correct D dimension name
    #else if (length(stan_data[[t]]) == stan_data[['P']]) {
    #  currString <- paste0("\tvector[P] ", t, ";\n")
    #} 
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

# functions to store max likelihood hyperparams -----
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
    #print(paste0(p.name, ' : ', p.val))
    tot.stan.data[[p.name]] = p.val
  }
  return(tot.stan.data)
}

# functions to backtransfrom normalisation of data  --------
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

