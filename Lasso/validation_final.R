library(stringr)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(plsRglm)
library(glmnet)
library(PRROC)
library(pROC)
library(data.table)
library(sva)
#library(DESeq2)
coding_gene_convert<-function(pre_list){
  pre_list=toupper(pre_list)
  pre_list[which(pre_list=='')]='-'
  protein_gene_dat=read.csv('protein-coding_gene.txt',sep='\t')
  protein_gene_dat$symbol=toupper(protein_gene_dat$symbol)
  converted_list=match(pre_list,protein_gene_dat$symbol)
  #converted_list[which(!is.na(converted_list))]=protein_gene_dat$symbol[converted_list[which(!is.na(converted_list))]]
  gene_set=str_split(toupper(protein_gene_dat$prev_symbol),'\\|')
  g <- rep(seq_along(gene_set), sapply(gene_set, length))
  converted_list[which(is.na(converted_list))]=g[match(pre_list[which(is.na(converted_list))], unlist(gene_set))]
  post_list=protein_gene_dat$symbol[converted_list]
  coding_gene_convert<-post_list
}
set.seed(2)
gene_list=c('PLCH1','SERPINB10','CCDC87','SPRY1','MT1G','NPAS1')
clinical_dat=fread('validation_sample.tsv')
exprs=fread('fpkm_raw_matrix_val.tsv',sep='\t')
gene_name=exprs$name
exprs=exprs[,-1]
exprs=data.frame(t(exprs))

colnames(exprs)=gene_name
exprs$id=row.names(exprs)
clinical_exprs_dat=merge(clinical_dat,exprs,by='id',all.x = F,all.y = F)
exprs$id=NULL
# batch=clinical_exprs_dat$batch
# mod <- model.matrix(~ 1, data = data.frame(batch = batch))
# expr_data_combat <- ComBat(dat = t(as.matrix(clinical_exprs_dat[,c(-1:-6)])), batch = batch, mod = NULL, par.prior = FALSE, prior.plots = FALSE)
# colnames(expr_data_combat)=clinical_exprs_dat$id
# write.csv(expr_data_combat,'fpkm_raw_matrix_val_rmbatch.csv')
exprs_rmbatch=as.data.frame(fread('fpkm_raw_matrix_val_rmbatch.csv',sep=','))
rownames(exprs_rmbatch) <- exprs_rmbatch$V1
exprs_rmbatch$V1=NULL
exprs_rmbatch=t(exprs_rmbatch)
exprs_rmbatch=data.frame(exprs_rmbatch)
# cor_list=rep(NA,ncol(exprs))
# for(i in 1:ncol(exprs)){
#   cor_list[i]=cor(exprs[,i],exprs_rmbatch[,i],method = 'spearman')
# }
# hist(cor_list)
exprs_rmbatch$id=row.names(exprs_rmbatch)
clinical_exprs_rmbatch_dat=merge(clinical_dat,exprs_rmbatch,by='id',all.x = F,all.y = F)
index0=which(clinical_exprs_rmbatch_dat$PNI==0)
index1=which(clinical_exprs_rmbatch_dat$PNI==1)
clinical_exprs_rmbatch_dat$PNI[index0]='Non-PNI'
clinical_exprs_rmbatch_dat$PNI[index1]='PNI'
clinical_exprs_rmbatch_dat$PNI=as.factor(clinical_exprs_rmbatch_dat$PNI)


temp1=apply(-0.1398607+exprs_rmbatch[c('SPRY1','MT1G','NPAS1')]*c(-0.0053208199961898,0.0100221845066824,0.0667469878828907),1,sum)
score = 1/(1+exp(-temp1))
clinical_exprs_rmbatch_dat$model1_score=score
pni=rep(0,nrow(clinical_exprs_rmbatch_dat))
pni[clinical_exprs_rmbatch_dat$PNI=='PNI']=1
pROC_obj <- roc.curve(weights.class0=pni,scores.class0=as.numeric(score),curve = T)
pROC_obj$auc=round(pROC_obj$auc,4)
print(plot(pROC_obj,color=3, cex.main =2,lwd=5))

p <- ggplot(clinical_exprs_rmbatch_dat,aes(x=PNI,y=model1_score,group=PNI,color=PNI)) + geom_boxplot(outlier.shape = NA) 
p <- p + geom_signif(comparisons = list(c('Non-PNI','PNI')),color='black')
p <- p + geom_jitter(size=2)
print(p)


