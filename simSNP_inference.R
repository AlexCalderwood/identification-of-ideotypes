library(here) # here_1.0.1
library(data.table) # data.table_1.13.2
library(ggplot2) # ggplot2_3.3.2
source('simSNP_inference_functions.R')


set.seed(100100100)

# set parameters for simulation
num.obs = 100 # number of observed "plants" to generate for each parameter set
num.parent.causal.SNPs = 5
num.child.causal.SNPs = 5
# num.non.causal.SNPs = 200 
p = 0.5 # probability snp = 1 (rather than 0)
b = 1 # coeffient of p.SNPs on parent trait
d = 1 # coefficient of c.SNPs on child trait
# g = 3 # coefficient of parent.trait on child trait
#parent.sd = 0.1 # sd of noise in generating parent trait
#child.sd = 0.1 # sd of noise in generating child trait
# shared.sd = 0.1 # sd or parent and child trait held to be equal
out.file <- 'simSNP_inference_out.csv'

all.results <- list()

# loop over varied parameters
count = 0
for (rep in c(1:10)) { # was 1:5
    for (g in c(0, 0.1, 0.5, 1, 3))  
      for (shared.sd in c(0.1, 0.5, 1)) {
        for (num.non.causal.SNPs in c(1000)) {
          
          count = count + 1
          print(paste0('loop: ', count))
          print(num.non.causal.SNPs)
          
          parent.sd <- shared.sd
          child.sd <- shared.sd
        
          # generate the simulated data
          sim.data <- simulate_data(num.obs, num.parent.causal.SNPs, num.child.causal.SNPs,
                                    num.non.causal.SNPs, p,
                                    b, parent.sd, 
                                    d, g, child.sd)
          
          
          # sanity plotting phenotypes and e.g. SNP
          # effect of child SNP
          # ggplot(sim.data, aes(x=c.SNP.1, y=parent.trait))+
          #   geom_point()
          # ggplot(sim.data, aes(x=c.SNP.1, y=child.trait))+
          #   geom_point()
          # ggplot(sim.data, aes(x=c.SNP.1, y=child.trait.residuals))+
          #   geom_point()
          # # effect of parent SNP
          # ggplot(sim.data, aes(x=p.SNP.1, y=parent.trait))+
          #   geom_point()
          # ggplot(sim.data, aes(x=p.SNP.1, y=child.trait))+
          #   geom_point()
          # ggplot(sim.data, aes(x=p.SNP.1, y=child.trait.residuals))+
          #   geom_point()
          
          # infer the probability each SNP is involved in each trait
          results <- calculate_significance_of_SNPs(sim.data)
          
          # calculate the variance explained by each SNP of each trait - useful summary 
          # of the parameters used to generate - for biology expect 10->20% for complex 
          # phenotypes
          results <- calulate_variance_explained_by_SNPs(sim.data, results)
          
          #calculate child trait variance explained by parent trait
          parent.trait.var <- (cor(sim.data$parent.trait, sim.data$child.trait))**2
          var.results <- results[results$trait=='child.trait',]
          mean.p.SNP.var <- mean(var.results$var.explained[grepl('p.SNP', var.results$SNP)])
          sd.p.SNP.var <- sd(var.results$var.explained[grepl('p.SNP', var.results$SNP)])
          mean.c.SNP.var <- mean(var.results$var.explained[grepl('c.SNP', var.results$SNP)])
          sd.c.SNP.var <- sd(var.results$var.explained[grepl('c.SNP', var.results$SNP)])
          
          # score inference performance
          # about whether you consider that the traits which act via the PARENT trait 
          # to act on the child trait or not.
          PR.results <- calculate_precision_recall(results, FDR.th=0.01)  
          
          # add parameter info to PR.results
          PR.results$rep <- rep
          PR.results$num.plants <- num.obs
          PR.results$num.p.SNPs <- num.parent.causal.SNPs
          PR.results$num.c.SNPs <- num.child.causal.SNPs
          PR.results$num.d.SNPs <- num.non.causal.SNPs
          PR.results$prob.SNP.on <- p
          PR.results$b <- b
          PR.results$parent.sd <- parent.sd
          PR.results$d <- d
          PR.results$g <- g
          PR.results$child.sd <- child.sd
          PR.results$parent.trait.var.explained <- parent.trait.var
          PR.results$p.SNP.var.explained <- mean.p.SNP.var
          PR.results$p.SNP.var.explained.sd <- sd.p.SNP.var
          PR.results$c.SNP.var.explained <- mean.c.SNP.var
          PR.results$c.SNP.var.explained.sd <- sd.c.SNP.var
        
          all.results <- c(all.results, list(PR.results)) 
        }
      }
    }

all.PR.results <- do.call('rbind', all.results)

write.csv(all.PR.results, file=out.file, row.names=F)
