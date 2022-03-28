rm(list=ls())
library(here) # here_1.0.1 
library(data.table) # data.table_1.13.2 
library(ggplot2) # ggplot2_3.3.2 
library(bnlearn) # bnlearn_4.5
library('Rgraphviz') # Rgraphviz_2.32.0
source('id_trait_structure_functions.R')



########## START OF MAIN ##########

# set RNG seed
set.seed(123456)

# Set parameters for testing -----
out.dir <- './id_trait_structure_output/'
data.file <- here('data', 'final_data_normed_cleaned_imputed.rds')
ecotype = 'Spring_OSR'
blacklistFile <- here('data', 'final_traits_blacklist.txt')
whitelistFile <- here('data', 'final_traits_whitelist.txt')
traits.of.interest <- c('height','num2aryBranch',
                        'numFlowers','timeToFlower',
                        'ovaryLen','styleLen','ovuleNum','ovuleArea','ovuleACov','gynLen',
                        'percAbortedMain','numPodsMain','podLen','beakLen',
                        'numPods2ary','percAborted2ary','timeToMature',
                        'tenPodNumber',
                        'tenPodArea','tenPodWeight',
                        'totSeedArea','totSeedNum','seedACov',
                        'totCompactness','tenPodCompactness',
                        'TGW',
                        'oilContent','totSeedW')


# set up variables and output directories -------
if (ecotype=='Spring_OSR') {
  ecotype <- 'Spring OSR'
} else if (ecotype=='Winter_OSR') {
  ecotype <- 'Winter OSR'
} else {
  ecotype <- ecotype
}
graph.dir <- here(out.dir, 'graphs/')
model.dir <- here(out.dir, 'data/')
dir.create(out.dir)
dir.create(graph.dir)
dir.create(model.dir)


# load the cleaned, imputed data -------
D <- readRDS(data.file)
D <- D[D$type==ecotype,] # cut down to "Spring OSR" or "Winter OSR
D <- subset(D, select=c(traits.of.interest, 'cultivar')) # cut down to just the trait data 
                                                         # and the cultivar info



# load whte/blacklists of required / forbidden links --------
bl <- read_blacklist(blacklistFile)
if (whitelistFile != 'NA') {
  print('using whitelist:')
  wl = fread(whitelistFile)
  print(wl)
  print('')
}


# Set up the k-fold data split -----
D.shuffled <- D[sample(nrow(D)), ] # shuffle rows for random subsetting
n.folds=5
D.data.cultivars <- subset(D.shuffled, select=c(traits.of.interest, 'cultivar'))
D.data2 <- subset(D.shuffled, select=c(traits.of.interest))
folds <- cut(seq(1, nrow(D.shuffled)), breaks=n.folds, labels=F)


# for each fold, learn DAG ------
results <- list(avg.dag.list=list(),
                str.df.list=list(),
                dag.plot.list=list(),
                xval.list=list(),
                scores=c(),
                sensitivity.list=list(),
                model.strings=c())
curr.fold <- 2
for (curr.fold in seq(1,n.folds)) {

  # split the data for the current fold
  testInds <- which(folds==curr.fold, arr.ind=T)
  test.Data <- D.data2[testInds, ]
  train.Data <- D.data2[-testInds, ]
  
  # learn DAG structure 
  if (exists('wl')) {
    str.df <- bnlearn::boot.strength(train.Data, R=500, algorithm='tabu', algorithm.args=list(blacklist=bl, whitelist=wl))
  } else {
    str.df <- bnlearn::boot.strength(train.Data, R=500, algorithm='tabu', algorithm.args=list(blacklist=bl))
  }
  
  keep.th <- attr(str.df, 'threshold') # calculated as Scutari M, Nagarajan R (2013).
  keep.th <- keep.th+0.05 
  avg.dag <- bnlearn::averaged.network(str.df, threshold=keep.th)
  
  # if has undirected arc (for this, both directions must be ok by blacklist), allocate arbitraty direction
  undirected.arcs <- undirected.arcs(avg.dag)
  if (nrow(undirected.arcs) != 0) {
    'has undirected arcs: randomly assigning direction'
    for (r in 1:nrow(undirected.arcs)) {
      avg.dag = set.arc(avg.dag, from=undirected.arcs[r, 1], to=undirected.arcs[r, 1])
    }
  }
  
  # store model for current fold data
  results$str.df.list[[curr.fold]] <- str.df
  results$avg.dag.list[[curr.fold]] <- avg.dag 
  results$scores <- c(results$scores, bnlearn::score(avg.dag, data=train.Data, type='bic-g'))
  results$model.strings <- c(results$model.strings, modelstring(avg.dag))
}

# Combine the models for each fold into an averaged model ------
# combine the str.df results together and average to make an "across the k-folds" averaged network
all.str.df <- calculate_average_str.df(n.folds, results)
all.dag <- averaged.network(all.str.df)

# save averaged across k-fold results to "results" (in the n.folds+1 position)
results$str.df.list[[n.folds+1]] <- all.str.df
results$avg.dag.list[[n.folds+1]] <- all.dag
results$scores <- c(results$scores, bnlearn::score(all.dag, data=D.data2, type='bic-g'))
results$model.strings <- c(results$model.strings, modelstring(avg.dag))

# for each fold found DAG, and for the averaged DAG found at the end, 
# fit model parameters, and make predictions about the with hold data ----
for (curr.fold in seq(1, n.folds)) { # just for the k-folds, NOT for the averaged DAG
  print(curr.fold)
  testInds <- which(folds==curr.fold, arr.ind=T)
  train.Data <- D.data2[-testInds, ]
  test.Data <- D.data2[testInds, ]
  avg.dag <- results$avg.dag.list[[curr.fold]]
  
  # fit model parameters, and predict performance on training data
  fitted <- bn.fit(avg.dag, data=train.Data)
  # predict new observations for each trait included
  results$xval.list <- c(results$xval.list, make_test_predictions(fitted, test.Data))
}

# convert lists to df as appropriate
results$xval.df <- do.call('rbind', results$xval.list) # predicted, and actual values for withheld data
results$xval.df$fold <- as.factor(results$xval.df$fold)
results$sensitivity.df <- do.call('rbind', results$sensitivity.list) # dataframe of the conditional distributions (Y|X)
results$sensitivity.df$fold <- as.factor(results$sensitivity.df$fold)


#### end of learning DAGs and making bnlearn predictions and sensitivitied---------
###################################################################################

# save bnlearn model results and the data used to generate it -------
saveRDS(results, file=paste0(model.dir, 'bnlearn_results.rds')) # the bnlearn resutls (stan needs the found model strings)
saveRDS(D.data.cultivars, file=paste0(model.dir, 'bnlearn_data.rds')) # the cut down data used for fitting the bnlearn models


# plot bnlearn results -------
# plot the dags found in each fold - for comparison
plot_DAGs(results$avg.dag.list, results$str.df.list, graph.dir, width=10, height=40)
# plot the strength plots for the dags
#plot_DAG_strength(results$avg.dag.list, results$str.df.list, results$scores, graph.dir, width=15, height=15)
# plot the test predicted values vs real values for each variable
p <- plot_prediction_vs_witheld_data(results$xval.df)
ggsave(paste0(graph.dir, 'trait_value_vs_predictions.pdf'), plot=p, width=10, height=9)


