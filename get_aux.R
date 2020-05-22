library(data.table)

library(biomaRt)
library(Matrix)

library(gage)
library('org.Hs.eg.db')

data('egSymb')
library(BiocManager)
library(TCGAbiolinks)
library(MASS)

#############3
######## kegg
# keggset = kegg.gsets(species = "hsa", id.type = "entrez", check.new=FALSE)
# kegg.gs.sym<-lapply(keggset$kg.sets, eg2sym)
# kegg = kegg.gs.sym[-keggset$dise.idx] ## pathways 
# 
# ## genes that are in kegg , 6860 genes, 239 paths 
# GENE = sort(unique(unlist(kegg)))
# GENE.PATH = sapply(GENE,function(gene) sapply(1:length(kegg), function(i) gene %in% kegg[[i]]))
# kegg_relation = t(GENE.PATH)%*%GENE.PATH
# 
# ##GET kegg correlation matrix 
# GENE.PATH = 1*(GENE.PATH)
# cor_kegg = cor(GENE.PATH)
# saveRDS(cor_kegg,file = 'cor_kegg.rds')


scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- rowSds(data, center = cm)
  csd[which(csd==0)]=0.001
  (data - cm) / csd
}


############### string 
library(STRINGdb)
string_db = STRINGdb$new( version="10", species=9606,score_threshold=0, input_directory="" )

get_string = function(dat){
  gene_of_interest = data.frame(gene = rownames(dat))
  string_id = string_db$map(gene_of_interest, "gene", removeUnmappedRows = TRUE )
  
  ## genes cannot be matched 
  xx = gene_of_interest[which(is.na(match(gene_of_interest[,1], string_id$gene))),1]
  
  ## build interaction matrix 
  interaction = string_db$get_interactions(string_id$STRING_id)[,c(1,2,16)]
  interaction$combined_score = interaction$combined_score/1000
  adj <- matrix(0, nrow(string_id), nrow(string_id))
  mat = as.matrix(interaction[,1:2])
  rownames(adj) = colnames(adj) = string_id$STRING_id
  
  adj[mat] <- interaction[,3]
  rownames(adj) = colnames(adj) = string_id$gene
  adj_t = t(adj)
  adj = adj_t+adj
  diag(adj)=1
  adj = data.frame(adj)
  nas = array(0,dim = c(length(xx),nrow(adj)))
  rownames(nas)  = xx
  adj = rbind(as.matrix(adj),as.matrix(nas))
  nas = array(0,dim = c(length(xx),nrow(adj)))
  rownames(nas)  = xx
  adj = cbind(adj,t(nas))
  return(adj)
}


########## biogrid 

get_biogrid = function(dat){
  names.genes.de <- rownames(dat)
  tmp.biogrid <- data.frame("Official.Symbol.Interactor.A" = names.genes.de,
                            "Official.Symbol.Interactor.B" = rev(names.genes.de))
  biogridd <- getAdjacencyBiogrid(tmp.biogrid, names.genes.de)
  bio = diag(colSums(biogridd))-biogridd
  biogrid = ginv(bio)
  biogrid = cov2cor(biogrid)
  colnames(biogrid) = colnames(bio)
  rownames(biogrid) = rownames(bio)
  return(biogrid)
}

######################################################
# dat = ref[[1]]
# new = yan
# dim(yan)
get_new = function(dat,new){
  # dat = sc_chu
  # new = yan
  genes = match(rownames(dat),rownames(new))
  nomatch = which(is.na(genes))
  yan_dat = new[na.omit(genes),]
  cor_ref = cor(t((log2(dat+1))))
  cor_yan = cor(t(scalematrix(yan_dat)))
  add = array(0,dim=c(length(nomatch),ncol(cor_yan)))
  rownames(add) = rownames(cor_ref)[nomatch]
  cor_y=rbind(cor_yan,add)
  add = array(0,dim=c(nrow(cor_y),length(nomatch)))
  colnames(add) = rownames(cor_ref)[nomatch]
  cor_y = cbind(cor_y,add)
  cor_y = cor_y[match(rownames(cor_ref),rownames(cor_y)),]
  cor_y = cor_y[,match(colnames(cor_ref),colnames(cor_y))]
  cor_y[nomatch,] = cor_ref[nomatch,]
  cor_y[,nomatch] = cor_ref[,nomatch]
  return(cor_y) 
}


get_bulk = function(dat,bulk){
  genes = intersect(intersect(rownames(dat),rownames(bulk)),rownames(cor_kegg)) ## 6038
  bulk = bulk[na.omit(match(genes,rownames(bulk))),]
  return(bulk)
}

