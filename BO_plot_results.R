library(here) # here_1.0.1 
library(ggplot2) # ggplot2_3.3.2
library(data.table) # data.table_1.13.2
source('BO_functions.R')

ecotype <- 'spring'
out.suffix <-'_PCA_4PC_constrained_0.01e' # the suffix used in the files for the setup (whether PCA or not)
in.dir <- here('BO_output')
out.dir <- here('BO_output')

# load the data produced by bayesian_optimisation.R --------------
stan.data <- readRDS(paste0(in.dir, '/', 
                            ecotype, '_stan.data', out.suffix, '.RDS'))
q.values <- readRDS(paste0(in.dir,  '/', 
                           ecotype, '_q.values', out.suffix, '.RDS'))


# Reformat it for plotting ---------------------------------------
# q points to data frame 
trait.names <- colnames(stan.data$X_totSeedW)
q.X.list <- list()
q.y.list <- list()
i=1
for (i in 1:length(q.values)) {
  pos <- q.values[[i]]$par.unscaled
  curr.df.X <- data.frame('q'=i,
                        'position'=pos,
                        'traits'=names(pos))
  curr.df.y <- data.frame('q'=i,
                        'yield'=q.values[[i]]$mu.unscaled,
                        'sd'=q.values[[i]]$sigma.unscaled,
                        'EI'=q.values[[i]]$value,
                        'times.found'=q.values[[i]]$times.found)
  q.X.list[[i]] <- curr.df.X
  q.y.list[[i]] <- curr.df.y
}
q.X.df <- do.call('rbind', q.X.list)
q.X.df$q <- as.factor(q.X.df$q)
q.y.df <- do.call('rbind', q.y.list)
q.y.df$q <- as.factor(q.y.df$q)
rm(q.X.list, curr.df.X, curr.df.y, q.values, q.y.list)


# observations to data frame
obs.df <- data.frame(stan.data$X_totSeedW.unscaled)
names(obs.df) <- colnames(stan.data$X_totSeedW)  
obs.df$yield <- stan.data$totSeedW.unscaled



# backtransform the data -------------------------------------------
# data was transformed to make more normally distributed prior to any modelling
# (so backtransform back from from e.g. log scale back to observed scale)
obs.back.df <- backtransform.data(obs.df)
obs.back.df.m <- melt(data.table(obs.back.df), id.vars='yield')
names(obs.back.df.m) <- c('yield', 'traits', 'position')
rm(obs.back.df, obs.df)

q.X.df.c <- dcast(data.table(q.X.df), q ~ traits, value.var='position')
q.X.back.df.c <- backtransform.data(q.X.df.c)
q.X.back.df <- melt(q.X.back.df.c, id.vars='q')
names(q.X.back.df) <- c('q', 'traits', 'position')
rm(q.X.df, q.X.back.df.c, q.X.df.c)

# combo the q-point X values with their y values
q.df <- merge(q.X.back.df, q.y.df, by='q')
rm(q.X.back.df, q.y.df)

## fix the trait names to final names:
look <- data.frame('old'=c('TGW', 'beakLen', 'num2aryBranch', 'numPods2ary', 
                           'numPodsMain', 'ovuleArea', 'tenPodWeight', 'totSeedNum',
                           'totSeedArea'),
                   'new'=c('TGW', 'BeakLength', 'NumberSecondInfl', 'NumberPods S', 
                           'NumberPods M', 'OvuleArea', 'SeedWeight M', 'SeedNumber',
                           'SeedArea'))
q.df$traits <- as.character(q.df$traits)
q.df$traits <- look$new[match(q.df$traits, look$old)]

obs.back.df.m$traits <- as.character(obs.back.df.m$traits)
obs.back.df.m$traits <- look$new[match(obs.back.df.m$traits, look$old)]

# hacky fix for rounding errors giving observations of less than 0
obs.back.df.m[obs.back.df.m < 0] <- 0

# Make the plots ----------------------------------------------------  
ggplot(q.df, aes(x=position, y=yield, color=q))+
  geom_point(data=obs.back.df.m, aes(x=position, y=mean(yield)), shape='|', inherit.aes=F)+
  geom_hline(data=obs.back.df.m, aes(yintercept=min(yield)), size=0.1)+
  geom_hline(data=obs.back.df.m, aes(yintercept=max(yield)), size=0.1)+
  geom_hline(data=obs.back.df.m, aes(yintercept=mean(yield)), size=0.5)+
  geom_point(data=q.df, aes(x=position, y=yield, color=q), size=1.2, alpha=0.8)+
  facet_wrap(~traits, scales='free_x')+
  xlab('')+
  ylab('seedYield')+
  labs(color='exploration priority')+
  guides(color=guide_legend(ncol=10))+
  theme_bw()+
  theme(legend.position='top', 
        legend.margin = margin(t=0, r=0, b=-10, l=0, unit='pt'),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=10),
        legend.title=element_text(size=10))

if (ecotype=='spring') {
ggsave(paste0(out.dir, '/', ecotype, out.suffix, '.pdf'),
       width=6, height=3.5)
} 
if (ecotype=='winter') {
  ggsave(paste0(out.dir, '/', ecotype, out.suffix, '.pdf'),
         width=6, height=5)
}

                           