

rm(list=ls())

drive         = "G:\\My Drive\\"
work_dir      = paste0(drive,"SOP\\methylation\\gsms\\exposure\\papers\\charlie\\paper\\kfold_simulation\\")
source(paste0(work_dir,"functions_mrs_simulation.R")) 
library(MASS)
library(glmnet)

set.seed(1)
n_sample = 500
n_pred   = 1000
n_sims   = 10000

n_folds               = 10               # number of fold for EN analysis, "loo" for leave one out
alpha                 = 0                # alpha = 0 is ridge regression that will keep all variables in the mode and alpha is 1 is lasso select the most predictive sites (issue is that it may select different sites for different folds but you could run it once on all sample to get a single set but keep the cross validation as an estimate of well that set predicts).
family                = "gaussian"       # use "binomial" for logistic regression and "gaussian" for normal regression
n_effects             = 0
r                     = 0.0 # size of the effects
if (n_pred<n_effects) stop("Cannot have more effects than predictors")


analysis_label = paste0("scale_nSims",n_sims,"_nfolds",n_folds,"_alpha",alpha) 

for (i in 1:n_sims) { #   i = 1

  cat(paste0("Simulation: ", i, "\n"))
  
  outcome   = rnorm( n_sample )
  sigma     = diag(rep(1,n_pred),n_pred,n_pred)
  if ( n_effects==0 ) {
    pred_data = mvrnorm(n_sample,rep(0,n_pred),sigma) 
      } else {
       pred_data = mvrnorm(n_sample,rep(0,n_pred),sigma) 
       pred_data[,1:n_effects] = r*outcome + sqrt(1-r^2)*pred_data[,1:n_effects]
    }
  
 
  row.names(pred_data) = paste0("sample_",1:n_sample)
  if (i==1) results =               run_regression(outcome,pred_data,n_folds,alpha,family)  else 
            results = rbind(results,run_regression(outcome,pred_data,n_folds,alpha,family) )
  
}

summary( results )

t.test(results[,"estimate"], mu = 0, alternative = "greater")

write.csv( results,paste0(work_dir,"result_",analysis_label,".csv"),row.names=T,quote=F)

par(mar=c(6,6,2,2))
opar=par(ps=15)
hist( results[,"estimate"],freq = T,  breaks = 50,cex=5, main="", xlab = "Predictive power"  )


