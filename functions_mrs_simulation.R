

run_regression = function(outcome,pred_data,n_folds,alpha=0,family="gaussian") {

  n_sample       = length(outcome)
  results        = matrix(NA,3) 
  
  if (n_folds=="loo") {
    n_folds = n_sample
    foldid  = 1:n_folds
  } else {
    fold_size = rep(floor(n_sample/n_folds),n_folds)
    rest      = n_sample-sum(fold_size)
    if (rest>0) fold_size[1:rest] = fold_size[1:rest] + 1
    foldid    = sample(rep(1:n_folds,fold_size),n_sample) # do stratified assigment to folds
    #    foldid    = sample(1:n_folds,size=nrow(pred_data), replace=TRUE)  
    # print(table(foldid))
  }

  cv.glmnet.result = cv.glmnet(pred_data,outcome,parallel=FALSE,alpha=alpha,family=family)
  lambda_min       = cv.glmnet.result$lambda.min
 
  for (i in 1:n_folds) { # i=1
    
    if(i %% 100==0) cat(paste0("Fold=: ", i, "\n"))
    
    train_data = pred_data[ which(foldid != i),,drop=F]
    train_out  = outcome[which(foldid != i)]
    test_data  = pred_data[ which(foldid == i),,drop=F]
    test_out   = outcome[which(foldid == i)]
    
    glmnet.result    = glmnet(train_data,train_out,lambda=lambda_min,alpha=alpha,family=family)	
    temp             = predict.glmnet(glmnet.result,test_data,alpha=alpha,family=family)
    temp             = cbind(test_out,temp )
      
    # coef( cv.glmnet.result)
    colnames(temp)   = c("outcome","mrs")
    rownames(temp)   = row.names(pred_data)[which(foldid == i)]  
    
    temp[,'mrs']=scale(temp[,'mrs']) # to avoud fold effects
    
   
       
    if (i==1) mrs = temp else mrs = rbind(mrs,temp)
    
  } #   for (i in 1:n_folds)
  
  results[1:3] = unlist( cor.test(mrs[,"outcome"],mrs[,"mrs"],alternative="greater")[c("estimate","statistic","p.value")])

  results        = as.vector(results)
  names(results) = c("estimate","statistic","pval")
  
  results 
  
}
