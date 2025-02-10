library(ggplot2)
train_dat=read.csv('test_dat_test.csv')


drug_dat=read.csv('drug_info.csv')
cell_dat=read.csv('sample_info.csv')
train_dat$name=drug_dat$drug_name[match(train_dat$name,drug_dat$id)]
#pdf('test_result.pdf')
cor1=cor(train_dat$value,train_dat$predicted_IC50,method = 'spearman')
cor1.pval=cor.test(train_dat$value,train_dat$predicted_IC50,method = 'spearman')$p.value
print(ggplot(train_dat,aes(x=value,y=predicted_IC50))+geom_point(size=0.0000000001)+xlab('IC50')+ylab('predicted IC50')+ggtitle(paste0("correlation = ",round(cor1,3),", p-value = ",formatC(cor1.pval,format = "e", digits = 2)))+geom_smooth(method = lm))
correlation_list=c()
p_value_list=c()
drug_list=c()
num_list=c()
for(each in unique(train_dat$name)){
  index1=which(train_dat$name==each)
  if(length(index1)>10){
    temp_dat=train_dat[index1,]
    cor1=cor(temp_dat$value,temp_dat$predicted_IC50,method = 'spearman')
    cor1.pval=cor.test(temp_dat$value,temp_dat$predicted_IC50,method = 'spearman',alternative ='greater' )$p.value
    print(ggplot(temp_dat,aes(x=value,y=predicted_IC50))+geom_point()+xlab('IC50')+ylab('predicted IC50')+ggtitle(paste0(each,", correlation = ",round(cor1,3),", p-value = ",formatC(cor1.pval,format = "e", digits = 2)))+geom_smooth(method = lm))
    correlation_list=c(correlation_list,cor1)
    p_value_list=c(p_value_list,cor1.pval)
    drug_list=c(drug_list,each)
    num_list=c(num_list,length(index1))
  }
}
result_dat=data.frame(cor=correlation_list,p_val=p_value_list,drug=drug_list,num=num_list)
result_dat=result_dat[order(result_dat$cor),]
result_dat$col1=ifelse(result_dat$p_val<0.05,'p_value<0.05','p_value>0.05')
result_dat$x=1:(dim(result_dat)[1])
print(ggplot(data = result_dat,aes(x=x,y=cor))+geom_point(aes(col=col1))+geom_line()+xlab('Drug')+ylab('Similarity Score'))
#dev.off()
write.csv(result_dat,'drug_cor_result.csv')