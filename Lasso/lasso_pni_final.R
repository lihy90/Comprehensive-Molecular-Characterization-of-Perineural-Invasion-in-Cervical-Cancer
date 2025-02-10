library(stringr)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(plsRglm)
library(glmnet)
library(PRROC)
library(pROC)

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
clinical_dat=read.csv('TableS1.csv')
exprs=read.csv('FPKM.45sample.tsv',sep='\t')
exprs=exprs[,-which(names(exprs) %in% c('P07','P20','P06','p45'))]
exprs$Gene=str_remove_all(exprs$Gene,'[^1-9A-Za-d]')
exprs$Gene=coding_gene_convert(exprs$Gene)
exprs=exprs[-which(is.na(exprs$Gene)),]
exprs=as.data.frame(exprs  %>% group_by(Gene)%>% mutate_all('max'))
exprs=exprs[-which(duplicated(exprs$Gene)),]
row.names(exprs)=exprs$Gene
exprs=exprs[,-1]
exprs=as.data.frame(t(exprs))
train_dat=merge(clinical_dat[,c(1,2,6)],exprs,by.x = 'Patient_ID',by.y = 'row.names')
train_dat=train_dat[,-1]
train_dat=train_dat[,c(1,2,which((sqrt(apply(train_dat[,c(-1,-2)],2,var))/apply(train_dat[,c(-1,-2)],2,mean))>1)+2)]
train_dat=train_dat[,c(1,2,which(apply(train_dat[,c(-1,-2)],2,mean)>1)+2)]
pred_score=c()
lamda_list=c()
for(j in 1:nrow(train_dat)){
  temp_dat=train_dat[-j,]
  
  pre_dat=train_dat[j,,drop=F]
  if(length(lamda_list)==0){
    lasso1=cv.glmnet(x =as.matrix(temp_dat[,-2]),y=temp_dat$PNI,nfolds=nrow(temp_dat),family = "binomial")
    lamda_list=lasso1$lambda
    predict_dat=data.frame(lamda=lamda_list)
  }else{
    lasso1=cv.glmnet(x =as.matrix(temp_dat[,-2]),y=temp_dat$PNI,family = "binomial",lambda = lamda_list)
  }
  
  pred_score1=-predict(lasso1,s=lasso1$lambda.min,newx=as.matrix(pre_dat[,-2]),type='response')+1
  pred_score=c(pred_score,pred_score1)
  print(j)
}
auc_list=c()
#png('lasso_auc_v5.png',width = 800,height =800)
pROC_obj <- roc.curve(weights.class0=train_dat$PNI,scores.class0=pred_score,curve = T)
pROC_obj$auc=round(pROC_obj$auc,4)
print(plot(pROC_obj,color=3, cex.main =2,lwd=5))

#dev.off()
result_dat=train_dat
result_dat$Status=ifelse(result_dat$PNI==1,'PNI','non-PNI')
result_dat$pred_score=pred_score
#png('lasso_boxplot_v5.png',width = 1200,height =800)
theme_set(
  theme_gray(base_size = 12)
)
p <- ggplot(result_dat,aes(x=Status,y=pred_score,group=Status,color=Status))+ geom_boxplot(outlier.shape = NA) +ylab('Predicted Score')+xlab('Status')+geom_signif(comparisons = list(c('PNI','non-PNI')),color='black')
p <- p + geom_jitter(size=2)

print(p)

#dev.off()
lasso1=cv.glmnet(x =as.matrix(train_dat[,-2]),y=train_dat$PNI,lambda = lamda_list,family = "binomial")
lasso1$lambda.min
write.csv(sort(lasso1$glmnet.fit$beta[,5][which(lasso1$glmnet.fit$beta[,5]!=0)]),'coef_lasso_final.csv')



exprs=read.csv('FPKM.45sample.tsv',sep='\t')
exprs=exprs[,-which(names(exprs) %in% c('P07','P20','P06','p45'))]
exprs$Gene=str_remove_all(exprs$Gene,'[^1-9A-Za-d]')
exprs$Gene=coding_gene_convert(exprs$Gene)
exprs=exprs[-which(is.na(exprs$Gene)),]
exprs=as.data.frame(exprs  %>% group_by(Gene)%>% mutate_all('max'))
exprs=exprs[-which(duplicated(exprs$Gene)),]
row.names(exprs)=exprs$Gene
exprs=exprs[,-1]
exprs=as.data.frame(t(exprs))
train_dat=merge(clinical_dat[,c(1,2,6)],exprs,by.x = 'Patient_ID',by.y = 'row.names')
train_dat$pred_score=pred_score

write.csv(train_dat[,c(1,2,3,17671)],'result_lasso_final.csv')