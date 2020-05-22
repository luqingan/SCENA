
single_methods = function(bulk,sc,cor_kegg,cor_string,cor_biogrid,cor_new){
  bulk = bulk[match(rownames(sc),rownames(bulk)),]
  ### normalize
  sc_norm = (log2(sc+1))
  cor_x = cor(t(scalematrix(sc_norm)))
  cor_x[is.na(cor_x)] = 0
  ## calculate correlation
  cor_b = cor(t(scalematrix(bulk)))
  
  sc_normalized_NA = sc_norm
  sc_normalized_NA[sc_norm == 0] <- NA ## replace 0 with NA
  cor_s = cor(t(sc_normalized_NA),use = "pairwise.complete.obs")
  
  ########### direct complete apply to dat
  cor_comp = CovarianceWithMissing(t(sc_normalized_NA))
  cor_shrink = cov.shrink(t(sc_norm))
  cor_shrink = as.matrix(cov2cor(cor_shrink))
  ####### apply to cor_s
  ## match biogrid
  ind = match(rownames(sc), rownames(cor_biogrid))
  cor_biogrid = cor_biogrid[ind,ind]
  cor_biogrid[is.na(cor_biogrid)]=0
  ## match string
  ind = match(rownames(sc), rownames(cor_string))
  cor_string = cor_string[ind,ind]
  cor_string[is.na(cor_string)]=0
  ## match new 
  ind = match(rownames(sc), rownames(cor_new))
  cor_new = cor_new[ind,ind]
  cor_new[is.na(cor_new)]=0
  
  ## match kegg
  ind = match(rownames(sc), rownames(cor_kegg))
  cor_kegg = cor_kegg[ind,ind]
  cor_kegg[is.na(cor_kegg)]=0
  
  ## methods
  cor_bs = correlation_NA_to_aux(sc_norm,cor_b,cor_s)
  cor_alpha = correlation_weighted_average(sc_norm,cor_b,cor_s)
  
  cor_ks = correlation_NA_to_aux(sc_norm,cor_kegg,cor_s)
  cor_alpha_kegg = correlation_weighted_average(sc_norm,cor_kegg,cor_s)
  
  cor_ss = correlation_NA_to_aux(sc_norm,cor_string,cor_s)
  cor_alpha_string = correlation_weighted_average(sc_norm,cor_string,cor_s)
  
  cor_sbio = correlation_NA_to_aux(sc_norm,cor_biogrid,cor_s)
  cor_alpha_biogrid = correlation_weighted_average(sc_norm,cor_biogrid,cor_s)
  
  ######### other data , need to deal with unmatched genes
  cor_ns = correlation_NA_to_aux(sc_norm,cor_new,cor_s)
  cor_alpha_n = correlation_weighted_average(sc_norm,cor_new,cor_s)
  
  cor = list(cor_x,cor_b,cor_bs,cor_alpha, ## bulk
             cor_kegg,cor_ks,cor_alpha_kegg, ## kegg
             cor_string,cor_ss,cor_alpha_string,## string
             cor_biogrid,cor_sbio,cor_alpha_biogrid, ## biogrid
             cor_comp,cor_shrink,## direct
             cor_new,cor_ns,cor_alpha_n) ## new data
  ## 18
  
  # cor_pd = lapply(cor, function(x) cov2cor(as.matrix(nearPD(x)$mat)))
  return(cor)
}

ridge = function(mini_ref,B,test,type ='b'){
  cor_mini_ref = cor(t(scalematrix(log2(mini_ref+1))))
  res = list()
  for (b in c(1:B)){
    set.seed(347*d+i*72+b*2)
    mini_samp = downsample(mini_ref,100) ## use rate=100, lower percent of zeros
    
    saver_mini = saver(mini_samp, do.fast = TRUE,ncores=8)
    dr_mini = DrImpute(mini_samp)
    prime_mini <- PRIME(mini_samp)
    scRMD_mini <- rmd(mini_samp)$exprs
    single_mini = try_noimpute_methods(aux[[i]][[1]],mini_samp,cor_kegg,aux[[i]][[2]],aux[[i]][[3]],aux[[i]][[4]])
    
    cor_x = single_mini[[1]]
    cor_saver <- cor.genes(saver_mini)
    cor_drimpute = correlation_aux(dr_mini)
    cor_scRMD= correlation_aux(scRMD_mini)
    cor_prime= correlation_aux(prime_mini)
    
    single_mini = c(single_mini,list(cor_saver,cor_drimpute,cor_scRMD,cor_prime))
    
    train = lapply(single_mini,function(mat) mat[upper.tri(mat, diag = FALSE)])
    train = data.frame(matrix(unlist(train), ncol=length(single_mini), byrow=FALSE))
    ref_train = cor_mini_ref[upper.tri(cor_mini_ref, diag = FALSE)]
    train = data.frame(y = ref_train,train)
    X <- as.matrix(train[,-1])
    Y <- train$y
    para = NULL
    X_long = NULL
    if (type=='b'){
      fit.ridge <- cv.glmnet(X, Y, family="gaussian",type.measure="mse", alpha=0)
      para = rbind(para,as.vector(predict(fit.ridge,type="coef", s = fit.ridge$lambda.1se)))
    }
    if (type=='a'){
      X_long = rbind(X_long,X)
    }
  }
  
  if (type=='b'){
    para_final = colMeans(para)
    res$para = para_final
    test_x = cbind(rep(1,nrow(test)),test)
    cor_y = as.matrix(test_x)%*%as.matrix(para_final)
    cor_pred = vec2cor(cor_y)
    cor_pred = cov2cor(as.matrix(nearPD(cor_pred)$mat))
    res$y = cor_pred
    return(res)
  }
  
  if (type=='a'){
    y_long = rep(Y,B)
    fit.ridge <- cv.glmnet(X_long, y_long, family="gaussian",type.measure="mse", alpha=0)
    newX = model.matrix(~.,data=test)
    newX = newX[,-1]
    cor_y <- predict(fit.ridge, s=fit.ridge$lambda.1se, newx=newX)
    cor_pred = vec2cor(cor_y)
    cor_pred = cov2cor(as.matrix(nearPD(cor_pred)$mat))
    
    res$para = as.vector(predict(fit.ridge,type="coef", s = fit.ridge$lambda.1se))
    res$y = cor_pred
    return(res)
  }
}

ridge_real = function(mini_ref,B,test,type ='b'){
  cor_mini_ref = cor(t(scalematrix(log2(mini_ref+1))))
  res = list()
  for (b in c(1:B)){
    set.seed(347*d+i*72+b*2)
    mini_samp = downsample(mini_ref,100) ## use rate=100, lower percent of zeros
    
    saver_mini = saver(mini_samp, do.fast = TRUE,ncores=8)
    dr_mini = DrImpute(mini_samp)
    prime_mini <- PRIME(mini_samp)
    scRMD_mini <- rmd(mini_samp)$exprs
    single_mini = try_noimpute_methods(aux_real[[i]][[1]],mini_samp,cor_kegg,aux_real[[i]][[2]],aux_real[[i]][[3]],aux_real[[i]][[4]])
    
    cor_x = single_mini[[1]]
    cor_saver <- cor.genes(saver_mini)
    cor_drimpute = correlation_aux(dr_mini)
    cor_scRMD= correlation_aux(scRMD_mini)
    cor_prime= correlation_aux(prime_mini)
    
    single_mini = c(single_mini,list(cor_saver,cor_drimpute,cor_scRMD,cor_prime))
    
    train = lapply(single_mini,function(mat) mat[upper.tri(mat, diag = FALSE)])
    train = data.frame(matrix(unlist(train), ncol=length(single_mini), byrow=FALSE))
    ref_train = cor_mini_ref[upper.tri(cor_mini_ref, diag = FALSE)]
    train = data.frame(y = ref_train,train)
    X <- as.matrix(train[,-1])
    Y <- train$y
    para = NULL
    X_long = NULL
    if (type=='b'){
      fit.ridge <- cv.glmnet(X, Y, family="gaussian",type.measure="mse", alpha=0)
      para = rbind(para,as.vector(predict(fit.ridge,type="coef", s = fit.ridge$lambda.1se)))
    }
    if (type=='a'){
      X_long = rbind(X_long,X)
    }
  }
  
  if (type=='b'){
    para_final = colMeans(para)
    res$para = para_final
    test_x = cbind(rep(1,nrow(test)),test)
    cor_y = as.matrix(test_x)%*%as.matrix(para_final)
    cor_pred = vec2cor(cor_y)
    cor_pred = cov2cor(as.matrix(nearPD(cor_pred)$mat))
    res$y = cor_pred
    return(res)
  }
  
  if (type=='a'){
    y_long = rep(Y,B)
    fit.ridge <- cv.glmnet(X_long, y_long, family="gaussian",type.measure="mse", alpha=0)
    newX = model.matrix(~.,data=test)
    newX = newX[,-1]
    cor_y <- predict(fit.ridge, s=fit.ridge$lambda.1se, newx=newX)
    cor_pred = vec2cor(cor_y)
    cor_pred = cov2cor(as.matrix(nearPD(cor_pred)$mat))
    
    res$para = as.vector(predict(fit.ridge,type="coef", s = fit.ridge$lambda.1se))
    res$y = cor_pred
    return(res)
  }
}
