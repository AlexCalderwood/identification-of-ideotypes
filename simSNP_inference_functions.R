build_corr.M <- function(num.causal.SNPs, num.non.causal.SNPs, c) {
  CC <- matrix(0, nrow=num.causal.SNPs, ncol=num.causal.SNPs)
  diag(CC) <- 1
  
  if (num.non.causal.SNPs == 0) {
    return(CC)
  }
  
  CNC <- matrix(c, nrow=num.causal.SNPs, ncol=num.non.causal.SNPs)
  NCC <- t(CNC)
  
  NCNC <- matrix(0, nrow=num.non.causal.SNPs, ncol=num.non.causal.SNPs)
  diag(NCNC) <- 1
  
  top <- cbind(CC, CNC)
  bottom <- cbind(NCC, NCNC)
  M <- rbind(top, bottom)
  
  return(M)
}

calculate_var_explained_old <- function(trait.value, SNP.table ) {
  # calculate the variance explained by each SNP as eta**2
  # http://davidmlane.com/hyperstat/B160638.html
  # varies depending on the number of causal SNPs, the size of B, and the trait noise!
  
  # trait.value = vector of the n.obs of the trait
  # SNPs = n.obs * n.SNPs matrix of 1/0 for the SNP complement of each observed plant
  
  global.mean <- mean(trait.value)
  sum.square.total <- sum((trait.value - global.mean)**2)
  
  # calculate the mean of the members of each SNP group
  #tmp <- SNP.table * as.vector(trait.value) 
  # for each SNP (col) calcualte the between group error (groups are SNP is off, and SNP is on)
  group.SSB <- apply(SNP.table, 2, calc.group.mean, trait.val=trait.value, global.mean=global.mean) 
  
  # calc. proportion of variance explained by each SNP
  prop.var.explained <- group.SSB / sum.square.total
  
  return(prop.var.explained)
}

calulate_variance_explained_by_SNPs <- function(sim.data, results) {
  # calculate the variance of each trait explained by each snp, following
  # http://davidmlane.com/hyperstat/B160638.html
  
  results$var.explained <- NA
  curr.trait <- 'parent.trait'
  curr.SNP <- 'p.SNP.1'
  for (curr.trait in unique(results$trait)) {
    for (curr.SNP in unique(results$SNP)) {
      
      trait.value <- sim.data[, curr.trait]
      snp.value <- sim.data[, curr.SNP]
      var.explained <- calculate_var_explained(trait.value, snp.value)
      
      results$var.explained[results$trait==curr.trait & 
                              results$SNP==curr.SNP] <- var.explained
      
    }
  }
  return(results)
}

calculate_var_explained <- function(trait.value, SNP.value) {
  # calculate the variance explained by each SNP as eta**2
  # http://davidmlane.com/hyperstat/B160638.html
  
  global.mean <- mean(trait.value)
  sum.square.total <- sum((trait.value - global.mean)**2)
  
  group.SSB <- calc.group.mean(SNP.value, trait.value, global.mean)
  
  prop.var.explained <- group.SSB / sum.square.total
  
  return(prop.var.explained)
}

calc.group.mean <- function(SNP.vec, trait.val, global.mean) {
  # should be run for each SNP.
  # for the groups SNP==1, and SNP==0, calculates the sum of square between groups, following
  # http://davidmlane.com/hyperstat/B86009.html#ssbg
  
  num.off <- sum(SNP.vec==0)
  num.on <- sum(SNP.vec==1)
  off.mean <- mean(trait.val[SNP.vec==0])
  on.mean <- mean(trait.val[SNP.vec==1])
  
  off.SSB <- num.off*((off.mean-global.mean)**2)
  on.SSB <- num.on*((on.mean-global.mean)**2)
  
  return(off.SSB + on.SSB)
}

calculate_precision_recall <- function(f.results, FDR.th) {
  # calculate precision = TP / (TP + FP)
  # for child trait, and child trait residuals
  # results : significance df
  # should be considered significant
  # FDR.th : the threshold value below which FDR rate is positive
  
  # not interested in the parent trait
  results <- f.results[f.results$trait %in% c('child.trait', 'child.trait.residuals'),]
  
  # vectors of results
  out.child.only <- c()
  out.trait <- c()
  out.FDR.th <- c()
  out.precision <- c()
  out.recall <- c()
  out.numParentSnps <- c()
  out.numChildSnps <- c()
  out.numDummySnps <- c()
  out.numParentSnpsIdd <- c()
  out.numChildSnpsIdd <- c()
  out.numDummySnpsIdd <- c()
  
  for (child.only in c(TRUE, FALSE)) {   # boolean for whether only the child.SNPs, or the child and parent.SNPs
    # define the real positives
    if (child.only==TRUE) { # if only consider the direct SNP as causal
      results$positive <- FALSE
      results$positive[grepl('c.SNP.', results$SNP)] <- TRUE
    } else { # if consider the SNPs causal to the child trait, or to the parent trait as positive
      results$positive <- FALSE
      results$positive[!(grepl('d.SNP.', results$SNP))] <- TRUE
    }
    
    # define the inferred positives
    results$inferred.positive <- results$pAdj <= FDR.th
    
    
    # calculate precision and recal for each trait
    results$true.positive <- (results$positive & results$inferred.positive)
    results$false.positive <- (results$inferred.positive & !(results$positive))
    
    for (curr.trait in unique(results$trait)) {
      curr.results <- results[results$trait==curr.trait,]
      
      if (sum(curr.results$inferred.positive)==0) {
        precision = 0
      } else {
        precision <- sum(curr.results$true.positive) / sum(curr.results$inferred.positive)
      }      
      recall <- sum(curr.results$true.positive) / sum(curr.results$positive)
      
      if(child.only==TRUE){
        child.only.str <- 'c.SNPs'
      } else {
        child.only.str <- 'p.SNPs+c.SNPs'
      }
      
      # count the true number of child and parent SNPs
      numChildSnps <- sum(grepl('c.SNP.', unique(curr.results$SNP)))
      numParentSnps <- sum(grepl('p.SNP.', unique(curr.results$SNP)))
      numDummySnps <- sum(grepl('d.SNP.', unique(curr.results$SNP)))
      # count the number of each of these identified as statistically associated with the child trait
      tmp <- unique(curr.results[, c('SNP', 'inferred.positive')])
      numChildSnpsIdd <- sum(grepl('c.SNP.', tmp$SNP) * tmp$inferred.positive)
      numParentSnpsIdd <- sum(grepl('p.SNP.', tmp$SNP) * tmp$inferred.positive)
      numDummySnpsIdd <- sum(grepl('d.SNP.', tmp$SNP) * tmp$inferred.positive)
      
      
      # save current results
      out.child.only <- c(out.child.only, child.only.str)
      out.trait <- c(out.trait, curr.trait)
      out.FDR.th <- c(out.FDR.th, FDR.th)
      out.precision <- c(out.precision, precision)
      out.recall <- c(out.recall, recall)
      out.numParentSnps <- c(out.numParentSnps, numParentSnps)
      out.numChildSnps <- c(out.numChildSnps, numChildSnps)
      out.numDummySnps <- c(out.numDummySnps, numDummySnps)
      out.numParentSnpsIdd <- c(out.numParentSnpsIdd, numParentSnpsIdd)
      out.numChildSnpsIdd <- c(out.numChildSnpsIdd, numChildSnpsIdd)
      out.numDummySnpsIdd <- c(out.numDummySnpsIdd, numDummySnpsIdd)
      
    } # loop over traits
  }# loop over child.only
  
  out.results <- data.frame('trait'=out.trait, 'SNPs.considered.causal'=out.child.only,
                            'FDR.th'=out.FDR.th, 'precision'=out.precision, 
                            'recall'=out.recall,
                            'num.p.SNPs'=out.numParentSnps,
                            'num.c.SNPs'=out.numChildSnps,
                            'num.d.SNPs'=out.numDummySnps,
                            'num.p.SNPs.idd'=out.numParentSnpsIdd,
                            'num.c.SNPs.idd'=out.numChildSnpsIdd,
                            'num.d.SNPs.idd'=out.numDummySnpsIdd)
                            
  return(out.results)
}

sample_SNPs <- function(n.obs, num.causal.SNPs, num.non.causal.SNPs, p) {
  # generate "SNPs" - bernoulli variables with desired correlation
  
  # n.obs = the number of observations (the number of plants)
  # num.causal.SNPs = the number of uncorrelated snps (we assume the 
  # causal SNPs are independent).
  # num.non.causal.SNPs = the number of SNPs which are correlated to the
  # "causal" ones (which are not correlated to each other)
  # p = the probability of a SNP being 1 - (in which case it contributes
  # to the phenotype).
  
  # c = the correlation between the non-causal, and teh causal SNPs (a
  # single, scalar value). - GET RID OF THIS AS A VARIABLE, JUST TREAT
  # ALL AS INDEPENDENT - DOESN'T WORK TOO CONTROLLABLY IF WANT CAUSAL 
  # TO BE INDEPENDENT, AND NON-CAUSAL DEPENDENT!!! - set all are independent!!
  
  # if only generating for 1 SNP:
  if (num.causal.SNPs==1) {
    b <- as.matrix(rbinom(n.obs, 1, p))
    return(b)
  }
  
  
  # setup params required for sampling
  n.rows <- n.obs
  d <- num.causal.SNPs + num.non.causal.SNPs
  p <- p
  mean.vec <- rep(p, d) # vec of means
  c <- 0.0 # wantet to try and vary this from 0, BUT doesn't work very well.
  # if c > ~0.25- easy to 
  # define an impossible matrix, and simulated values aren't that similar anyway...
  # so just leave unvaried @ 0 for now.
  corr.M <- build_corr.M(num.causal.SNPs, num.non.causal.SNPs, c)
  
  # generate binomial samples
  # my homemade copula method doesn't seem to work well!
  # MultiRNG uses some weird algo for binary observations
  # correlation amoung binomial are much more similar to desired!
  b <- MultiRNG::draw.correlated.binary(no.row=n.rows,
                                        d=d,
                                        prop.vec=mean.vec,
                                        corr.mat=corr.M)
  
  #cor(b)
  return(b)
}

simulate_data <- function(num.obs, num.parent.causal.SNPs, num.child.causal.SNPs,
                          num.non.causal.SNPs, p,
                          b, parent.sd, 
                          d, g, child.sd) {
  # p is probability SNP==1 (rather than 0)
  # b is coefficient of parent SNPs on parent
  # parent.df is standard deviation of parent trait noise after SNPs
  # d is coefficient of child SNPs on child trait
  # g is coefficient of parent trait on child trait
  # child.sd is sd of the noise of the child trait
  
  # generate the SNPs which control the PARENT TRAIT
  parent.trait.SNPs <- sample_SNPs(num.obs, num.parent.causal.SNPs, 0, p)
  #cor(parent.trait.SNPs)
  
  # generate not involved SNPs (potential false positives) - slow, and want all to be independent now anyway!, so this is a waste
  #dummy.SNPs <- sample_SNPs(num.obs, num.non.causal.SNPs, 0, p)
  #cor(dummy.SNPs)
  # make all non.causal independent
  dummy.SNPs <- matrix(0, nrow=num.obs, ncol=num.non.causal.SNPs)
  for (i in 1:num.non.causal.SNPs) {
    dummy.SNPs[, i] <- rbinom(num.obs, 1, p)
  }
  
  # generate the SNPs which control the CHILD TRAIT
  child.trait.SNPs <- sample_SNPs(num.obs, num.child.causal.SNPs, 0, p)
  
  # bind all the snps together
  all.SNPs <- cbind(parent.trait.SNPs, child.trait.SNPs, dummy.SNPs)
  
  # generate the parent trait (needed, as need to control for it!), using parent trait SNPs
  parent.trait <- generate.parent.trait(parent.trait.SNPs, b=b, parent.sd=parent.sd)
  #plot(x=parent.trait.SNPs, y=parent.trait)
  
  # generate the child trait using the parent (rather thant the Parent.SNPs directly)
  #d = 1 # the coefficients of the child trait SNPs
  #g = 1 # the coefficient of the parent trait - by keeping as 1 means than just can subtract rather than
  # having to infer...
  child.trait <- generate.child.trait(child.trait.SNPs, parent.trait, d, g, child.sd=child.sd) 
  
  # combine datas, and convert SNPs to factors
  sim.data <- combine_generated_data(parent.trait, child.trait, parent.trait.SNPs, 
                                     child.trait.SNPs, dummy.SNPs)
  rm(parent.trait, child.trait, parent.trait.SNPs, child.trait.SNPs, dummy.SNPs, 
     all.SNPs)
  
  # calculate "residual" child trait, and add to "sim.data" as a column
  sim.data <- calculate_child_residual(sim.data)
  
  return(sim.data)
}

generate.parent.trait <- function(parent.trait.SNPs, b, parent.sd) {
  # parent.trait.SNPs : matrix 0/1 for whether each SNP in each obs
  # B the scalar coefficient for how much each causal SNP affects the trait (B in linear model)
  # all are the same
  # parent.sd : sd of the (normally distributed) noise in the parent trait
  # returns a vec of the parent trait values
  
  num.causal.SNPs <- ncol(parent.trait.SNPs)
  num.obs <- nrow(parent.trait.SNPs)
  
  # coefficents for each SNP
  B <- as.vector(rep(b, num.causal.SNPs))
  # parent noise
  parent.noise <- rnorm(num.obs, mean=0, sd=parent.sd)
  # linear model to add together - can't make a 1x1 matrix multiply work
  if (num.causal.SNPs > 1) {
    parent <- (parent.trait.SNPs[, 1:num.causal.SNPs] %*% B) + parent.noise
  } else {
    parent <- (parent.trait.SNPs * B) + parent.noise
  }
  
  return(parent)
}

generate.child.trait <- function(child.trait.SNPs, parent.trait, d, g, child.sd=3) {
  # generate the child trait 
  # d = coefficient for child SNP effectiveness
  # g = coefficient for parent trait effect on child trait
  
  num.causal.SNPs <- ncol(child.trait.SNPs)
  num.obs <- nrow(child.trait.SNPs)
  
  # coefficents for each SNP
  D <- as.vector(rep(d, num.causal.SNPs))
  
  # parent noise
  child.noise <- rnorm(num.obs, mean=0, sd=child.sd)
  
  # linear model to add together - can't make a 1x1 matrix multiply work
  if (num.causal.SNPs > 1) {
    child <- (child.trait.SNPs %*% D) + g * parent.trait + child.noise
  } else {
    child <- (child.trait.SNPs * D) + g * parent.trait + child.noise
  }
  
  return(child)
}

combine_generated_data <- function(parent.trait, child.trait, parent.trait.SNPs, 
                                   child.trait.SNPs, dummy.SNPs) {
  
  sim.data <- data.frame('parent.trait'=parent.trait, 'child.trait'=child.trait)
  parent.snps <- data.frame(apply(parent.trait.SNPs, 2, as.factor), stringsAsFactors = TRUE)
  names(parent.snps) <- paste0('p.SNP.', 1:ncol(parent.snps))
  
  child.snps <- data.frame(apply(child.trait.SNPs, 2, as.factor), stringsAsFactors = TRUE)
  names(child.snps) <- paste0('c.SNP.', 1:ncol(child.snps))
  
  dummy.snps <- data.frame(apply(dummy.SNPs, 2, as.factor), stringsAsFactors = TRUE)
  names(dummy.snps) <- paste0('d.SNP.', 1:ncol(dummy.snps))
  
  sim.data <- cbind(sim.data, parent.snps, child.snps, dummy.snps)
  
  return(sim.data)
}

calculate_significance_of_SNPs <- function(sim.data) {
  # for each column with "trait" in the name, fit a 1 SNP model using each column with "SNP" in as a predictor
  # and save the P greater than value, and the traitwise BH adjusted p value 
  
  all.SNPs <- names(sim.data)[grepl('SNP', names(sim.data))]
  all.traits <- names(sim.data)[grepl('trait', names(sim.data))]
  
  out.trait <- c()
  out.SNP <- c()
  out.pVal <- c()
  out.pAdj <- c()
  # curr.trait <- 'parent.trait'
  # curr.SNP <- 'p.SNP.1'
  for (curr.trait in all.traits) {
    curr.pVals <- c() # collect per trait, as will apply FDR correction per trait
    for (curr.SNP in all.SNPs) {
      model.string <- paste0(curr.trait, ' ~ ', 'factor(', curr.SNP,')') # don't enforce a 0 intercept, as then forces ALL SNPs to be signif, as have to explaint why non-zero!!!
      model <- lm(as.formula(model.string), data=sim.data)
      s <- summary(model)
      p.val <- s$coefficients[2,4] # probability of the SNP involved
      
      out.trait <- c(out.trait, curr.trait)
      out.SNP <- c(out.SNP, curr.SNP)
      curr.pVals <- c(curr.pVals, p.val)
    }
    # for each trait, calculate the FDR adjusted p-value
    curr.pAdj <- p.adjust(curr.pVals, method='BH')
    
    out.pVal <- c(out.pVal, curr.pVals)
    out.pAdj <- c(out.pAdj, curr.pAdj)
  }
  
  results <- data.frame('trait'=out.trait, 'SNP'=out.SNP, 'pVal'=out.pVal, 'pAdj'=out.pAdj)
}

calculate_child_residual <- function(sim.data) {
  # fit model of child.trait ~ parent.trait, and use to calculate child trait residuals, and add to sim.data
  
  model <- lm(child.trait ~ parent.trait, data=sim.data)
  pred.child.trait <- predict(model, newdata=sim.data)
  
  sim.data$child.trait.pred <- pred.child.trait
  
  ggplot(sim.data, aes(x=parent.trait, y=child.trait))+
    geom_point(color='blue')+
    geom_point(aes(x=parent.trait, y=child.trait.pred), color='red')
  
  sim.data$child.trait.residuals <- sim.data$child.trait - sim.data$child.trait.pred
  sim.data$child.trait.pred <- NULL  
  
  return(sim.data)
}
