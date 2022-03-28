
read_blacklist <- function(blacklistFile) {
  # read blacklist file, and convert to bnlearn compatible format
  bl.vec <- scan(blacklistFile, what="", sep="\n")
  bl.list <- list()
  for (line in bl.vec) {
    curr.line <- strsplit(line, ',')
    bl.list <- c(bl.list, curr.line)
  }
  bl <- bnlearn::tiers2blacklist(bl.list)
  return(bl)
}

get_sensitivity_ranges <- function(DAG, data,  n.bootstraps) {
  # estimate coefficients relating traits using data. Calculates the mean estimate, using the full dataset, 
  # and the 95% confidence interval by bootstrapping the data.
  
  fitted <- bn.fit(DAG, data)
  # get sample values
  true_estimates <- get_sensitivity_point_estimate(fitted, data)
  
  # get bootstrapped values
  bootstrap.ests <- list()
  for (i in seq(1, n.bootstraps)) {
    # sample random rows from the data
    b.data <- data[sample(nrow(data), nrow(data), replace=TRUE), ]
    b.fitted <- bn.fit(DAG, b.data)
    curr.estimates <- get_sensitivity_point_estimate(b.fitted, b.data)
    names(curr.estimates)[names(curr.estimates)=='Values'] <- 'boot.Value'
    curr.estimates <- merge(true_estimates, curr.estimates, by=c('Y', 'X', 'Vars'))
    bootstrap.ests[[i]] <- curr.estimates
  }
  bootstrap.ests.df <- do.call('rbind', bootstrap.ests)
  bootstrap.ests.df$delta <- bootstrap.ests.df$boot.Value - bootstrap.ests.df$Values
  # get the d values for each parameter
  
  # confidence intervals - 95% CIs - 
  tmp <- quantile(c(0.01, 0.5, 0.99), probs=c(0.025))
  bootstrap.ests.df[, lower.delta:=quantile(delta, probs=0.025), by=.(Y,X,Vars)]
  bootstrap.ests.df[, upper.delta:=quantile(delta, probs=0.975), by=.(Y,X,Vars)]
  
  CI.df <- unique(bootstrap.ests.df[, c('Y', 'X', 'Vars', 'Values', 'lower.delta', 'upper.delta')])
  CI.df$Values.UL <- CI.df$Values + CI.df$upper.delta
  CI.df$Values.LL <- CI.df$Values + CI.df$lower.delta
  CI.df$lower.delta <- NULL 
  CI.df$upper.delta <- NULL
  return(CI.df)
}

get_sensitivity_point_estimate <- function(fitted, data) {
  # calculate conditional parameter coefficients for a single fitted model, 
  # uses the data it was fit with for calculating sds.
  
  fitted.rbmn <- bnfit2nbn(fitted)
  fitted.gema <- nbn2gema(fitted.rbmn)
  fitted.mn <- gema2mn(fitted.gema)
  
  #Y <- names(fitted)[1]
  # Y <- 'TSW'
  # X <- 'NP'
  Ys <- c()
  Xs <- c()
  values <- c()
  vars <- c()
  for (Y in names(fitted)) {
    for (X in names(fitted)) {
      if (X != Y) {
        L <- unlist(condi4joint(fitted.mn, par=Y, pour=X, x2=NULL))
        mean <- L[1]
        coefficient <- L[2] # coefficient per X unit
        sd <- L[3]
        
        #sd.tst <- sd(data$NP) 
        X.sd <- sd(data[, X, with=F][[1]])
        coeff.sd <- coefficient / X.sd # coefficient per standard deviation of X
        
        Ys <- c(Ys, c(Y,Y,Y,Y))
        Xs <- c(Xs, c(X,X,X,X))
        values <- c(values, c(mean, coefficient, sd, coeff.sd))
        vars <- c(vars, c('mean', 'coeff.per.unit', 'sd', 'coeff.per.sd'))
      }
    }
  }
  # make output table
  out.df <- data.table('Y'=Ys, 'X'=Xs, 'Values'=values, 'Vars'=vars)
  return(out.df)
}

make_test_predictions <- function(fitted, test.Data) {
  # take fitted model, test data and make predictions for 
  # values of test data. 
  # returns a list of data frames
  
  bnmodel <- modelstring(fitted)
  predicted.traits <- get_predicted_traits(bnmodel) # traits which aren't root nodes
  
  xval.list <- list()
  # test.trait = 'total.seed.weight'
  for (test.trait in predicted.traits) {
    y.hat <- predict(fitted, node=test.trait, data=test.Data)
    y <- test.Data[, test.trait, with=F]
    tmp.df <- data.table('y'=y, 'y.hat'=y.hat, 'trait'=test.trait, 'fold'=curr.fold)
    names(tmp.df) <- c('y', 'y.hat', 'y.trait', 'fold')
    tmp.df <- cbind(tmp.df, test.Data)
    xval.list <- c(xval.list, list(tmp.df))
  }
  return(xval.list)
}

# making nicer DAG plot...
# https://stackoverflow.com/questions/29282522/setting-edge-width-in-rgraphviz
setEdgeAttr <- function( graph, attribute, value, ID1, ID2) {
  
  idfunc <- function(x) paste0( sort(x), collapse="~" )
  all.ids <- sapply( AgEdge(graph), 
                     function(e) idfunc( c( attr(e, "head"), attr( e, "tail" ))))
  
  #sel.ids <- names(value)
  
  sel.ids <- apply(cbind( ID1, ID2 ), 1, idfunc )
  setdiff(sel.ids, all.ids)
  if(!all(sel.ids %in% all.ids)) stop( "only existing edges, please" )
  sel <- match( sel.ids, all.ids )
  
  for(i in 1:length(sel)) {
    attr( attr( graph, "AgEdge" )[[ sel[i] ]], attribute ) <- value[i]
  }
  
  return(graph)
}

setNodeAttr <- function( graph, attribute, value) {
  idfunc <- function(x) paste0( sort(x), collapse="~" )
  all.names <- sapply( AgNode(graph), 
                       function(e) idfunc( attr(e, "name")))
  sel.ids <- names(value)
  if(!all(sel.ids %in% all.names)) stop( "only existing nodes, please" )
  sel <- match( sel.ids, all.names )
  
  for(i in 1:length(sel)) {
    attr( attr( graph, "AgNode" )[[ sel[i] ]], attribute ) <- value[[i]]
  }
  
  return(graph)
}

setNodeLabelAttr <- function( graph, attribute, value) {
  
  idfunc <- function(x) paste0( sort(x), collapse="~" )
  
  all.names <- sapply( AgNode(graph), 
                       function(e) idfunc( attr(e, "name")))
  sel.ids <- names(value)
  if(!all(sel.ids %in% all.names)) stop( "only existing nodes, please" )
  sel <- match( sel.ids, all.names )
  
  for(i in 1:length(sel)) {
    attr( attr( attr( graph, "AgNode" )[[ sel[i] ]], 'txtLabel' ), attribute) <- value[[i]]
  }
  
  return(graph)
}

reformat.rag <- function(L, i, curr.strs, curr.rag, should.show.nodes) {
  
  # get node positions for avg.strength.plot
  
  # set colours
  edgeCols <- L[[i]]@renderInfo@edges$col # get the edges in the current plot
  eNames <- names(edgeCols)
  ID1 <- tstrsplit(names(edgeCols), '~')[[1]]
  ID2 <- tstrsplit(names(edgeCols), '~')[[2]]
  curr.rag <- setEdgeAttr(curr.rag, 'color', edgeCols, ID1, ID2)
  #plot(curr.rag)
  
  # set weights
  # calculate the edgeWidths
  curr.strs
  edge <- eNames[1]
  weights <- c()
  for (edge in eNames) {
    f <- strsplit(edge, '~')[[1]][1]
    t <- strsplit(edge, '~')[[1]][2]
    str <- curr.strs$strength[curr.strs$from==f & curr.strs$to==t]
    if (length(str) == 0) { # if edge not in this DAG
      str = 0
    }
    weights <- c(weights, str)
    #weights(edge=str)
  }
  weights <- (scale(weights, center=F))*0.5
  names(weights) <- eNames
  ID1 <- tstrsplit(names(weights), '~')[[1]]
  ID2 <- tstrsplit(names(weights), '~')[[2]]
  curr.rag <- setEdgeAttr(curr.rag, 'lwd', weights, ID1, ID2)
  #plot(curr.rag)
  
  # set edge line stype (lty)
  # calculate lty
  edgeLty <- c()
  for (e in eNames) {
    if (edgeCols[[e]] == 'black') {
      curr.lty <- 'solid'
    } else {
      curr.lty <- 'dotted'
    }
    edgeLty <- c(edgeLty, curr.lty)
  }
  names(edgeLty) <- eNames
  ID1 <- tstrsplit(names(edgeLty), '~')[[1]]
  ID2 <- tstrsplit(names(edgeLty), '~')[[2]]
  curr.rag <- setEdgeAttr(curr.rag, 'lty', edgeLty, ID1, ID2)
  
  # if doing the k1-5, just want shapes, as node labels will be too small
  if (!(should.show.nodes)) {
    # set fill color to black
    ID1 <- L[[i]]@nodes
    value <- rep('black', length(ID1))
    names(value) <- ID1  
    curr.rag <- setNodeAttr(curr.rag, 'fillcolor', value)
    #plot(curr.rag)
    
    # set font size
    ID1 <- L[[i]]@nodes
    value <- rep(0, length(ID1)) # fontsize
    names(value) <- ID1
    curr.rag <- setNodeLabelAttr(curr.rag, 'labelFontsize', value)
    
    # set shape
    ID1 <- L[[i]]@nodes
    value <- rep('ellipse', length(ID1))
    names(value) <- ID1  
    curr.rag <- setNodeAttr(curr.rag, 'shape', value)
    #plot(curr.rag)
    
    # set width
    ID1 <- L[[i]]@nodes
    value <- rep(25, length(ID1))
    names(value) <- ID1
    curr.rag <- setNodeAttr(curr.rag, 'rWidth', value)
    curr.rag <- setNodeAttr(curr.rag, 'lWidth', value)
  }
  # if doing the k-avg, just node labels
  if (should.show.nodes) {
    # set shape
    ID1 <- L[[i]]@nodes
    value <- rep('ellipse', length(ID1))
    names(value) <- ID1  
    curr.rag <- setNodeAttr(curr.rag, 'shape', value)
    #plot(curr.rag)
    
    # set width
    ID1 <- L[[i]]@nodes
    value <- rep(25, length(ID1))
    names(value) <- ID1
    curr.rag <- setNodeAttr(curr.rag, 'rWidth', value)
    curr.rag <- setNodeAttr(curr.rag, 'lWidth', value)
    
    # set height
    ID1 <- L[[i]]@nodes
    value <- rep(25, length(ID1))
    names(value) <- ID1
    curr.rag <- setNodeAttr(curr.rag, 'height', value)
    
    # set font size
    ID1 <- L[[i]]@nodes
    value <- rep(40, length(ID1)) # fontsize
    names(value) <- ID1
    curr.rag <- setNodeLabelAttr(curr.rag, 'labelFontsize', value)
    
  }
  
  return(curr.rag)
}

plot_raw_data <- function(D, traits.of.interest) {
  D.plot <- subset(D, select=traits.of.interest)
  p <- ggpairs(D.plot)
  p <- p+
    theme_bw()+
    theme(
      strip.text=element_text(size=6), 
      axis.text=element_text(size=6)
    )
  return(p)
}

plot_DAGs <- function(avg.dag.list, str.df.list, graph.dir, width, height) {
  
  pdf.path <- paste0(graph.dir, 'k-fold_DAGs.pdf')
  # make the big k-fold comparison plot
  pdf(file=pdf.path, width=width, height=height)
  par(mfrow=c(n.folds,1))
  L = graphviz.compare(avg.dag.list[[6]],
                       avg.dag.list[[1]],avg.dag.list[[2]], avg.dag.list[[3]],avg.dag.list[[4]],avg.dag.list[[5]],
                       shape='ellipse', 
                       main=c('k-averaged', 'k-1', 'k-2', 'k-3', 'k-4', 'k-5'))
  dev.off()
  
  
  # make the avg DAG have edge thickness == weight
  avg <- L[[1]]
  edgeCols <- avg@renderInfo@edges$col # get the edges in the current plot
  eNames <- names(edgeCols)
  
  curr.strs <- str.df.list[[6]] # get the averaged strengths
  edge <- eNames[1]
  weights <- c()
  for (edge in eNames) {
    f <- strsplit(edge, '~')[[1]][1]
    t <- strsplit(edge, '~')[[1]][2]
    str <- curr.strs$strength[curr.strs$from==f & curr.strs$to==t]
    if (length(str) == 0) { # if edge not in this DAG
      str = 0
    }
    weights <- c(weights, str)
    #weights(edge=str)
  }
  weights <- weights*7
  names(weights) <- eNames
  avg@renderInfo@edges$lwd <- weights # change the edgeWeights for the avg plot
  
  pdf(file=paste0(graph.dir, 'avg_DAG_strength.pdf'), width=10, height=10)
  par(mfrow=c(1,1))
  renderGraph(avg)
  dev.off()
  
  
  # pdf(paste0(graph.dir, 'nice_DAG_kfold.pdf'), width = width, height=(width/4))
  # par(mfrow=c(1,5))
  # for (i in 2:6) {
  #   curr.rag <- agopen(L[[i]], '')
  #   curr.str <- str.df.list[[i-1]] # this is ordered so the avg one is at position 6
  #   curr.rag <- reformat.rag(L, i, curr.str, curr.rag, should.show.nodes=FALSE)
  #   plot(curr.rag, main = paste0('k-', i-1))
  # }
  # dev.off()
  
  
}

plot_DAG_strength <- function(avg.dag.list, str.df.list, scores, graph.dir, width, height) {
  pdf.path <- paste0(graph.dir,'k-fold_DAGs_strength.pdf')
  pdf(file=pdf.path, width=width, height=height)
  par(mfrow=c(floor(sqrt(n.folds)),ceiling(sqrt(n.folds))))
  for (i in seq(1,n.folds+1)){
    strength.plot(avg.dag.list[[i]], str.df.list[[i]], shape='ellipse', render=T, main=paste0('SCORE:', round(scores[i])))
  }
  dev.off()
}

plot_prediction_vs_witheld_data <- function(xval.df) {
  p <- ggplot(xval.df, aes(x=y, y=y.hat, color=fold))+
    geom_point(size=1)+
    facet_wrap(~y.trait, scales='free')+
    theme_bw()+
    theme(strip.text=element_text(size=8))+
    geom_abline(intercept=0, slope=1, linetype='dashed')
  return(p)
}

plot_pairwise_predictions <- function(xval.df) {
  # makes pairwise plots of predicted and measured data for each pairwise combination of traits - want to check 
  # that not overvalueing e.g. increasing branch number is actually there's a peak value.
  # each prediction is made using all the traits except for the withheld one.
  # only do it for the traits actually making predictions for!
  
  # y vals are predicted and measured. X vals are measured only.
  
  i <- unique(xval.df$y.trait)[1]
  p.list <- list()
  # iterate over all the traits
  i <- 'pod_len'
  for (i in unique(xval.df$y.trait)) {
    curr.df <- xval.df[xval.df$y.trait==i, ]
    curr.df.m <- melt(curr.df, id.vars=c('y', 'y.hat', 'y.trait', 'fold'))
    
    names(curr.df.m) <- c('y', 'y.hat', 'y.trait', 'fold', 'x.name', 'x.true')
    curr.df.m.m <- melt(curr.df.m, id.vars=c('y.trait', 'fold', 'x.name', 'x.true'))
    
    p <- ggplot(curr.df.m.m, aes(x=x.true, y=value))+
      geom_point(aes(colour=variable), size=0.5, alpha=0.6)+
      facet_wrap(~x.name, nrow=1, scales='free')+
      ylab(i)+
      theme_bw()+
      theme(axis.text=element_text(size=6),
            axis.title.x=element_blank(),
            legend.position='top')
    p
    p.list <- c(p.list, list(p))
    
  }
  p.grid <- plot_grid(plotlist = p.list, ncol=1)
  return(p.grid) 
}

plot_trait_sensitivity <- function(sensitivity.df) {
  p.list.major <- list()
  curr.var <- 'sd'
  curr.t <- 'ASW'
  for (curr.var in unique(sensitivity.df$Vars)) {
    p.list.minor <- list()
    for (curr.t in unique(sensitivity.df$Y)) {
      curr.df <- sensitivity.df[sensitivity.df$Y==curr.t & sensitivity.df$Vars==curr.var, ]
      p <- ggplot(curr.df, aes(x=X, y=Values, color=fold))+
        geom_point()+
        geom_errorbar(aes(ymin=Values.LL, ymax=Values.UL), width=0.4)+
        theme_bw()+
        theme(legend.position = 'none', 
              axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=6), 
              axis.text.y = element_text(size=6),
              axis.title = element_text(size = 8))+
        xlab('')+
        ylab('')
      
      if(curr.var==unique(sensitivity.df$Vars[1])) {
        p <- p+ylab(curr.t)
      }
      if(curr.t == unique(sensitivity.df$Y)[1]){
        p <- p + ggtitle(curr.var)+theme(plot.title = element_text(size=8))
      }
      #p
      p.list.minor <- c(p.list.minor, list(p))
    }
    p.grid <- plot_grid(plotlist=p.list.minor, ncol=1)
    #p.grid
    p.list.major <- c(p.list.major, list(p.grid))
  }
  p.list.major[[1]]
  p.grid <- plot_grid(plotlist=p.list.major, nrow=1)
  return(p.grid)
}

calculate_average_str.df <- function(n.folds, results) {
  # calculate the edge weights for each link by weighted averaged over all the k-folds. 
  # edges are weighted according to the bic scores of the network they're found in.
  
  curr.fold = 1
  for (curr.fold in seq(1, n.folds)) {
    weight = abs(results$scores[curr.fold]) / sum(abs(results$scores)) # to weight each result by its relative score. (assumes that all will be negative)
    tmp.str.df <- results$str.df.list[[curr.fold]]
    tmp.str.df$strength <- tmp.str.df$strength*weight
    tmp.str.df$direction <- tmp.str.df$direction*weight
    if (curr.fold == 1) {
      all.str.df <- tmp.str.df
    } else {
      all.str.df <- rbind(all.str.df, tmp.str.df)
    }
  }
  # strength is probability that a link exists between the two nodes
  # direction is probability that is in that direction
  summed.strength <- aggregate(.~from+to, 
                               all.str.df, 
                               FUN=sum)

    for (row in 1:nrow(summed.strength)){
    all.str.df$from[row] <- summed.strength$from[row]
    all.str.df$to[row] <- summed.strength$to[row]
    all.str.df$direction[row] <- summed.strength$direction[row]
    all.str.df$strength[row] <- summed.strength$strength[row]
  }
  all.str.df <- all.str.df[1:nrow(summed.strength), ]
  
  return(all.str.df)
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

parseModelString <- function(bnmodel) {
  # takes bnlearn style model string, and converts it to a list of the parents of each node.
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
