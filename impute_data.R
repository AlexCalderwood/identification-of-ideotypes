library(here) # here_1.0.1
library(mice) # mice_3.11.0
library(data.table) # data.table_1.13.2
library(sjmisc) # sjmisc_2.8.6 

impute_missing_data <- function(D) {
  # impute missing values using mice, fit 5 imputed values for each missing value
  # see https://www.rdocumentation.org/packages/mice/versions/2.25/topics/mice.impute.pmm
  # may refuse to impute a column if is too highly correlated with others.
  mice_imputes <- mice::mice(D, m=5, method='pmm', maxit = 40)

  # merge the imputed values with the original data. 
  M <- sjmisc::merge_imputations(
    D,
    mice_imputes,
    summary=c('dens')
  )
  
  # define information columns (as opposed to trait measurement columns) 
  D2 <- M$data
  if ('KtrtID' %in% names(D)) {
    info.cols <- c('gHouse', 'R', 'C', 'Bnarb', 'BN', 'type', 'KtrtID', 'trtID', 'rep', 'abr', 'cultivar')
  } else {
    info.cols <- c('gHouse', 'R', 'C', 'Bnarb', 'BN', 'type', 'abr', 'cultivar')
  }
  info.df <- subset(D, select=info.cols)
  D2 <- cbind(info.df, D2)
  
  # check nothing lost somehow
  stopifnot(length(setdiff(names(D), names(D2)))==0)
  
  out <- list(D2, M$plot)
  return(out)
}

set.seed(123)


D <- readRDS(here('data', 'final_data_normed_cleaned.rds'))
O <- impute_missing_data(D)
D.imp <- O[[1]]
p <- O[[2]]

# output is:
# saveRDS(D.imp, file=here('final_data_normed_cleaned_imputed.rds'))