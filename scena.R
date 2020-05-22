library(SAVER) 
library(DrImpute)
library(scRMD)
library(PRIME)
library(mclust)

###### change link_data to interested one 
setwd("~/Dropbox/GQ - Applications/Single-cell/data")
#### load code
# source('install.R')
source('single_ensemble.R')
source('graphfun.R')
source('get_aux.R')

## load raw data 
load('chu.rdata')
load('chu_time.rdata')
load('darmanis.rdata')

## gene match 
bulk = readRDS("RNA_data_norm_hg19_all_name.rds") ## 41989*444
cor_kegg = readRDS('cor_kegg.rds')
## other sc for chu
yan = readRDS('yan.rds')
yann <- assay(yan, "logcounts")
## for darmanis 
lake = readRDS('lake.rds') 
lakee <- assay(lake, "normcounts")
############################

match_gene = function(sc,bulk,cor_kegg){
  genes = intersect(intersect(rownames(sc),rownames(bulk)),rownames(cor_kegg)) ## 6038
  return(sc[na.omit(match(genes,rownames(sc))),])
}

sc_chu = match_gene(chu$sc_cnt,bulk,cor_kegg)
sc_chu_time = match_gene(chu_time$sc_cnt,bulk,cor_kegg)
sc_darmanis = match_gene(darmanis$sc_cnt,bulk,cor_kegg)
# aux = readRDS('network.rds')
label = list(chu$sc_label,chu_time$sc_label,darmanis$sc_label)

###############functions 
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- rowSds(data, center = cm)
  csd[which(csd==0)]=0.001
  (data - cm) / csd
}

## get off-diagonal 
upper = function(mat){
  mat[upper.tri(mat, diag = FALSE)]
}

get_mse = function(cor_x,cor_ref){
  x = upper(cor_x)
  y = upper(cor_ref)
  return(mean((x-y)^2,na.rm = TRUE))
}

### ARI
get_percent = function(eigenvalue){
  percent = NULL
  for (p in c(1:length(eigenvalue))){
    percent[p] = sum(eigenvalue[1:p])/sum(eigenvalue)
  }
  return(percent)
}

get_range = function(acc_var,a,b){
  pcs = which(acc_var<b&acc_var>a)
  perc = acc_var[pcs]
  p = rbind(perc,pcs)
  rownames(p) = c('percent','Num_pc')
  return(p)
}

### clustering 
fit_cluster = function(mydata){
  d0 = dist(mydata,method = 'manhattan')
  fit <- hclust(d0,method = 'ward.D')
  return(fit)
}

## pca
get_pca = function(data,eigen_vector){
  pca = data.frame((as.matrix(t(data))%*%eigen_vector))
  return(pca)
}




build_ref = function(x1,label){
  # Filter out library size greater than 4500000
  x2 <- x1[, which(colSums(x1) <= quantile(colSums(x1),0.95))]
  label = label[which(colSums(x1) <= quantile(colSums(x1),0.95))]
  x3 <- x2[rowMeans(x2) >= quantile(rowMeans(x2),0.25), ]
  x4 <- x3[rowSums(x3 != 0) >= quantile(rowSums(x3 != 0),0.15), ]
  
  # build reference
  lib.size <- colSums(x4)
  non.zero.prop <- apply(x4, 1, function(x) sum(x != 0)/length(x))
  cells.filt <- which(lib.size > quantile(lib.size,0.05))
  genes.filt <- which(non.zero.prop > quantile(non.zero.prop,0.5))
  data.filt <- x4[genes.filt, cells.filt]
  label = label[cells.filt]
  datt = list(data.filt,label)
  return(datt)
}

ref_chu = build_ref(sc_chu,chu$sc_label)
ref_chu_time = build_ref(sc_chu_time,chu_time$sc_label)
ref_darmanis = build_ref(sc_darmanis,darmanis$sc_label)

ref = list(ref_chu[[1]],ref_chu_time[[1]],ref_darmanis[[1]])
label_ref = list(ref_chu[[2]],ref_chu_time[[2]],ref_darmanis[[2]])
# saveRDS(label_ref,file = 'label_ref.rds')


data_name = c('chu','chut','darmanis')
cor_ref <- vector("list", 3)
names(cor_ref) <- data_name
for (i in 1:3) { 
  ref_norm = log2(ref[[i]]+1)
  cor_ref[[i]] = cor(t(ref_norm))
}
####################################################

downsample = function(refe,rate){
  alpha <- rgamma(ncol(refe), 10, rate = rate) 
  data.samp <- t(apply(sweep(refe, 2, alpha, "*"), 1, function(x)
    rpois(length(x), x)))
  colnames(data.samp) = colnames(refe)
  return(data.samp)
}
model_name = c('Reference','X_s','SAVER','drImpute','scRMD','PRIME','SCENA_average','SCENA_ridge')
data_name = c('chu','chu_time',"darmanis")

correlation_aux = function(x_aux){
  cor_aux = cor(t(scalematrix(log2(x_aux+1))))
  cor_aux[is.na(cor_aux)]=0
  return(cor_aux)
}
ave = function(list_of_cor){
  x =Reduce('+',list_of_cor)/length(list_of_cor)
  x = as.matrix(data.frame(x[1:nrow(x),1:ncol(x)]))
  cor_pred = cov2cor(as.matrix(nearPD(x)$mat))## find nearest pd and change it to -1 1 
  return(cor_pred)
}
########### get aux 
aux_chu = list(get_bulk(ref[[1]],bulk),get_string(ref[[1]]),get_biogrid(ref[[1]]),get_new(ref[[1]],yann))
aux_chu_time = list(get_bulk(ref[[2]],bulk),get_string(ref[[2]]),get_biogrid(ref[[2]]),get_new(ref[[2]],yann))
aux_darmanis = list(get_bulk(ref[[3]],bulk),get_string(ref[[3]]),get_biogrid(ref[[3]]),get_new(ref[[3]],lakee))
aux = list(aux_chu,aux_chu_time,aux_darmanis)
#######conduct 30 downsampling 

## 
rate = c(3000,1000,1000)

mse <- vector("list", 3)
cmd <- vector("list", 3)
ARI_all <- vector("list", 3)
cor_noimpute = vector("list", 3)
B = 1 ## number of iteration for mini-downsample 
stan_x = NULL
for (d in c(1:10)){ ## downsample repetition
  for (i in c(1:3)){ ## data loop 
    print(c(i,d))
    r = rate[i]
    set.seed(i*21+36*d)
    ### generate rv from gamma with mean 0.1 for each sample 
    alpha <- rgamma(ncol(ref[[i]]), 10, rate = r) ## mean is 0.1 for each sample 
    ## for each gene, generate poisson with mean = alpha*entry 
    data.samp <- t(apply(sweep(ref[[i]], 2, alpha, "*"), 1, function(x)
      rpois(length(x), x)))
    colnames(data.samp) = colnames(ref[[i]])
    stan_x = scalematrix(log2(data.samp+1))
    
    ###### imputation 
    saver = saver(data.samp, do.fast = TRUE,ncores=8)
    dr = DrImpute(data.samp)
    prime <- PRIME(data.samp)
    scRMD <- rmd(data.samp)$exprs
    ############ do single methods to downsampled 
    
    single = single_methods(aux[[i]][[1]],data.samp,cor_kegg,aux[[i]][[2]],aux[[i]][[3]],aux[[i]][[4]])
    ################# have no imputation at all
    cor_x = single[[1]]
    cor_saver <- cor.genes(saver)
    cor_drimpute = correlation_aux(dr)
    cor_scRMD= correlation_aux(scRMD)
    cor_prime= correlation_aux(prime)
    ################# have imputation only in model stacking 
    single= c(single,list(cor_saver,cor_drimpute,cor_scRMD,cor_prime))
    cor_ave= ave(single)
    ##################################################
    test = lapply(single,function(mat) mat[upper.tri(mat, diag = FALSE)])
    test =data.frame(matrix(unlist(test), ncol=length(single), byrow=FALSE))
    ########## ridge model stacking (build mini-ref and mini-down to fit the model)
    ridge = ridge(mini_ref,B,test,'b')
    cor_ridge= ridge$y 
    cor_all = list(cor_ref[[i]],cor_x,cor_saver,cor_drimpute,cor_scRMD,cor_prime,cor_ave,cor_ridge)
    names(cor_all) = model_name
    cor_noimpute[[i]][[d]]=cor_all
    
    ############### for result 
    ## MSE
    mse[[i]][[d]] = sapply(cor_all, function(x) get_mse(x,cor_ref[[i]]))
    cmd[[i]][[d]] = sapply(cor_all, function(x) cmd(x,cor_ref[[i]]))    
    ## eigens 
    eigens = lapply(cor_all, function(x) eigen(x))
    eigen_vec = lapply(eigens, function(x) x$vectors)
    eigen_value = lapply(eigens, function(x) x$values)
    
    impo = lapply(eigen_value, function(x) get_percent(x))
    pc_range = lapply(impo, function(acc_var) get_range(acc_var,0.9,0.99)) 
    
    pca = lapply(eigen_vec, function(ev) get_pca(stan_x,ev))  
    
    ######### ari 
    ari = NULL
    pca_tran = lapply(pca, function(x) scalematrix((x)))
    start = sapply(pc_range, function(x) x[[2,1]])
    end = sapply(pc_range, function(x) x[[2,ncol(x)]])
    byy = round((end-start)/20)
    
    for (m in c(1:length(cor_all))){
      pcs = seq(start[m],end[m],by=byy[m])
      acc_cell = array(0,dim = c(length(pcs),3))
      for (pp in c(1:length(pcs))){
        pc = pcs[pp]
        fit_models= fit_cluster(pca_tran[[m]][,1:pc])
        K = length(unique(label_ref[[i]]))
        cls = cutree(fit_models, k=K)
        acc_cell[pp,3] = adjustedRandIndex(label_ref[[i]], cls)
      } ## end PCs clustering loop
      acc_cell[,1] = rep(model_name[m],pp)
      acc_cell[,2] = impo[[m]][pcs]
      ari = rbind(ari,acc_cell)
    }
    
    colnames(ari) = c('Method','percent','ARI')
    ARI_all[[i]][[d]] = ari
    saveRDS(cor_noimpute,file = ('cor.rds')) 
    saveRDS(mse,file = ('mse.rds')) 
    saveRDS(cmd,file = ('cmd.rds')) 
    saveRDS(ARI_all,file =('ari.rds')) 
  }
}
