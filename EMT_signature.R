###############################################
# Topic: EMT_REO signature
# Author: Wangkai
##############################################

####使用emt2refergenes(CV<0.3)构建基因对标志####
rm(list=ls())
setwd('E:/colorectal_sur/GSE_TCGA/')

#####
#1.识别CV基因####
load('data/emt.Rdata')
myfun_div <- function(GSE,cli){
  gid=GSE[,1];exp=GSE[,-1]
  cli1=cli[grep('1',cli$stage),]
  exp1=exp[,match(cli1[,1],colnames(exp))]
  cli2=cli[grep('3|4',cli$stage),]
  exp2=exp[,match(cli2[,1],colnames(exp))]
  data=list(gid=gid,exp1=exp1,exp2=exp2)
  return(data)
}
myfun_cv <- function(x){
  cv=sd(x)/mean(x)
  return(cv)
}
#GSE39582
load('data/GSE39582_inf.Rdata')
data=myfun_div(GSE,cli)
gid1=data$gid;exp1=as.matrix(data$exp1);exp2=as.matrix(data$exp2)
exp=cbind(exp1,exp2)
cv1=apply(exp, 1, myfun_cv)
cvgenes1=gid1[which(cv1<0.3)]
#TCGA
load('data/TCGA_inf.Rdata')
data=myfun_div(fpkm,cli)
gid2=data$gid;exp1=as.matrix(data$exp1);exp2=as.matrix(data$exp2)
exp=cbind(exp1,exp2)
# log2exp=log2(exp+1)
log2exp=log2(exp)
log2exp[log2exp<0]=0
cv2=apply(log2exp, 1, myfun_cv)
cvgenes2=gid2[which(cv2<0.3)]
cvgenes=intersect(cvgenes1,cvgenes2)
gid=intersect(gid1,gid2)
emt2plat=intersect(emt,gid)#识别emt与检测平台交叠基因
CvMinusEmt2plat=setdiff(cvgenes,emt2plat)
refergenes=c(emt2plat,CvMinusEmt2plat)
save(emt2plat,refergenes,file='result/emt_refergenes.Rdata')
# NcvEmt=setdiff(emt2plat,cvgenes)
# save(NcvEmt,cvgenes,file='refergenes.Rdata')

#2.识别srpairs_emt2refergenes改进算法(fisher)####
# load('GSE_TCGA/emt_refergenes.Rdata')
load('refergenes.Rdata')
myfun_div <- function(GSE,cli){
  gid=GSE[,1];exp=GSE[,-1]
  cli1=cli[grep('1',cli$stage),]
  exp1=exp[,match(cli1[,1],colnames(exp))]
  cli2=cli[grep('3|4',cli$stage),]
  exp2=exp[,match(cli2[,1],colnames(exp))]
  data=list(gid=gid,exp1=exp1,exp2=exp2)
  return(data)
}

myfun_fisher <- function(x){
  ftest=fisher.test(matrix(x,2,2,byrow=T))
  fd=x[1]/(x[1]+x[2])-x[3]/(x[3]+x[4])
  fresult=c(fd,ftest$estimate,ftest$p.value)
  return(fresult)
}

myfun_srpairs <- function(GSE,cli,genes1,genes2){
  begin=Sys.time()
  data=myfun_div(GSE,cli)
  gid=data$gid;exp1=as.matrix(data$exp1);exp2=as.matrix(data$exp2)
  Lgenes1=length(genes1);Lgenes2=length(genes2)
  Lsam1=ncol(exp1);Lsam2=ncol(exp2)
  allpairs=matrix(NA,nrow=(((Lgenes2-Lgenes1)+(Lgenes2-1))*Lgenes1/2),ncol=5)
  m=1
  for(i in 1:(Lgenes1)){
    print(i)
    diff1=rep(1,times=(Lgenes2-i))%*%exp1[match(genes1[i],gid),,drop=F]-exp1[match(genes2[(i+1):Lgenes2],gid),,drop=F]
    diff2=rep(1,times=(Lgenes2-i))%*%exp2[match(genes1[i],gid),,drop=F]-exp2[match(genes2[(i+1):Lgenes2],gid),,drop=F]
    count1=rowSums(diff1>0);count2=rowSums(diff2>0)
    pairs=cbind(genes1[i],genes2[(i+1):Lgenes2])
    sep=cbind(count1,(Lsam1-count1),count2,(Lsam2-count2))
    index=which((sep[,1]/(sep[,1]+sep[,2]))>0.5&(sep[,4]/(sep[,3]+sep[,4]))>0.5
                |(sep[,2]/(sep[,1]+sep[,2]))>0.5&(sep[,3]/(sep[,3]+sep[,4]))>0.5)
    if(length(index)!=0){
      sep_index=sep[index,,drop=F]
      pairs_index=pairs[index,,drop=F]
      result <- t(apply(sep_index,1,myfun_fisher))
      fresult0=cbind(pairs_index,result)
      n=m+nrow(fresult0)-1
      allpairs[m:n,]=fresult0
      m=n+1
    }else next
  }
  allpairs=na.omit(allpairs)
  allpairs=cbind(allpairs,p.adjust(allpairs[,5],method='BH'))
  colnames(allpairs) <- c('Ga','Gb','fd','estimate','p','fdr')
  srpairs <- allpairs[allpairs[,6]<0.05,]
  srpairs <- rbind(srpairs[srpairs[,3]<0,c(2,1,3:6)],srpairs[srpairs[,3]>0,])
  srpairs[,3]=abs(srpairs[,3])
  srpairs=srpairs[order(srpairs[,3],decreasing = T),]
  #colnames(srpairs) <- c('Ga','Gb','fd','estimate','p','fdr')
  fresult=list(allpairs=allpairs,srpairs=srpairs)
  end=Sys.time()
  print(end-begin)
  return(fresult)
}
#GSE39582
load('data/GSE39582_inf.Rdata')
fresult <- myfun_srpairs(GSE,cli,emt2plat,refergenes)
save(fresult,file='srpairs/srpairs_GSE39582.Rdata')
#TCGA
load('data/TCGA_inf.Rdata')
fresult <- myfun_srpairs(fpkm,cli,emt2plat,refergenes)
save(fresult,file='srpairs/srpairs_TCGA.Rdata')

#3.基因对交叠####
load('srpairs/srpairs_GSE39582.Rdata')
srpairs=fresult$srpairs
str1=paste(srpairs[,1],srpairs[,2],sep=',')
load('srpairs/srpairs_TCGA.Rdata')
srpairs=fresult$srpairs
str2=paste(srpairs[,1],srpairs[,2],sep=',')
str22=paste(srpairs[,2],srpairs[,1],sep=',')
costr=intersect(str1,str2)
costr22=intersect(str1,str22)#costr22=1,交叠一致性99.9%
srpairs=matrix(unlist(strsplit(costr,split=',')),ncol=2,nrow=length(costr),byrow=T)
save(srpairs,file='srpairs/cosrpairs.Rdata')

#4.prognosis_related pairs####
load('data/GSE39582_inf.Rdata')#GSE24551
load('srpairs/cosrpairs.Rdata')
myfun_cox <- function(GSE,cli,pairs){
  gid=GSE[,1];exp=GSE[,-1]
  cli=cli[which(cli$stage=='2'),]
  #cli=cli[which(cli$ACT=='N'),]
  cli=cli[!is.na(cli$RFS_event),]
  exp=exp[,match(cli$sample,colnames(exp))]
  pair_cox=c()
  library(survival)
  for(i in 1:nrow(pairs)){
    print(i)
    diff=exp[match(pairs[i,1],gid),]-exp[match(pairs[i,2],gid),]
    sam1=diff[which(diff>=0)]
    sam2=diff[which(diff<0)]
    if(ncol(sam1)>0&ncol(sam2)>0){
      cli1=cli[match(colnames(sam1),cli[,1]),]
      cli2=cli[match(colnames(sam2),cli[,1]),]
      library(plyr)
      cli11=transform(cli1,class='0')#control
      cli22=transform(cli2,class='1')#case,,,hr=1:0
      cli_new=rbind(cli11,cli22)
      label=cli_new$class
      survi=Surv(cli_new$RFS_time,cli_new$RFS_event) #Surv(sur_time,sur_state)
      cox=coxph(survi~label+location+gender+age+MMR+CIMP+CIN+TP53+KRAS+BRAF,data=cli_new)#survi~label+location+gender+age+MMR+CIMP+CIN+TP53+KRAS+BRAF
      b=summary(cox)
      p = b$coefficients[1,5]
      cindex = b$concordance[1]
      hr = b$conf.int[1,1]
      lower.hr = b$conf.int[1,3]
      upper.hr = b$conf.int[1,4]
      pair0=t(c(pairs[i,1:2],p,cindex,hr,lower.hr,upper.hr))
      pair_cox=rbind(pair_cox,pair0)
    }
    else next
  }
  colnames(pair_cox)=c('Ga','Gb','p','cindex','hr','down_hr','up_hr')
  return(pair_cox)
}#p<0.05&hr>1
pairs_cox <- myfun_cox(GSE,cli,srpairs)
pairs_sur <- pairs_cox[which(pairs_cox[,3]<0.05&pairs_cox[,5]>1),]
save(pairs_cox,pairs_sur,file='result/pairs_multicox.Rdata')

#5.结合算法#####结合算法后发现naiveBayes的分类效果最好，后续使用nb算法####
#整理数据
myfun_div <- function(GSE,cli){
  gid=GSE[,1];exp=GSE[,-1]
  cli1=cli[grep('1',cli$stage),]
  exp1=exp[,match(cli1[,1],colnames(exp))]
  cli2=cli[grep('3|4',cli$stage),]
  exp2=exp[,match(cli2[,1],colnames(exp))]
  data=list(gid=gid,exp1=exp1,exp2=exp2)
  return(data)
}

myfun_diff <- function(GSE,cli,pairs){
  data=myfun_div(GSE,cli)
  gid=data$gid;exp1=data$exp1;exp2=data$exp2
  exp=cbind(exp1,exp2)
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  diffindex=as.data.frame(t(ifelse(diff>0,1,0)))
  colnames(diffindex)=c(paste('pair',1:ncol(diffindex),sep=''))
  risk=factor(c(rep('L',ncol(exp1)),rep('H',ncol(exp2))),levels=c('L','H'))#1--lowrisk;-1--highrisk
  diffdata=cbind(diffindex,risk)
  diff=list(diffindex=diffindex,diffdata=diffdata)
  return(diff)
}
load('GSE_TCGA/data/GSE39582_inf.Rdata')
diffresult1=myfun_diff(GSE,cli,pairs_sur)
load('GSE_TCGA/data/TCGA_inf.Rdata')
diffresult2=myfun_diff(fpkm,cli,pairs_sur)
diffindex=rbind(diffresult1$diffindex,diffresult2$diffindex)
diffdata=rbind(diffresult1$diffdata,diffresult2$diffdata)
save(diffindex,diffdata,file='GSE_TCGA/diffresult5175.Rdata')

#5_1.SVM####
#install.packages('e1071')
library(e1071)
load('GSE_TCGA/diffresult.Rdata')
#筛选最优model参数
tuned <- tune.svm(risk ~.,data=diffdata,gamma=10^(-6:-1),cost=10^(-2:2))
summary(tuned)
#构建model
svm.fit<-svm(risk~.,data=diffdata,gamma=0.001,cost=10)#kernel = "linear"
summary(svm.fit)
#验证model
svm.pre=predict(svm.fit,newdata=diffindex)
(index=table(actual=diffdata$risk,predict=svm.pre))
(consistency=sum(diag(index))/sum(index))
save(svm.fit,index,consistency,file='GSE_TCGA/1.svm/model.Rdata')
#5_2.Naive Bayes####
library(e1071)
load('GSE_TCGA/diffresult.Rdata')
nb.fit <- naiveBayes(risk~.,data=diffdata)#nb要求输入数据为data.frame,分类信息为因子数据
summary(nb.fit)
nb.pre <- predict(nb.fit,newdata =diffindex)#预测结果
#生成实际与预测交叉表和预测精度
(index <- table(actual=diffdata$risk,predict=nb.pre))
(consistency=sum(diag(index))/sum(index))
save(nb.fit,index,consistency,file='GSE_TCGA/2.naiveBayes/model5175.Rdata')

#5_3.decision tree####
# install.packages('rpart')
library(rpart)
load('GSE_TCGA/diffresult.Rdata')
# install.packages('maptree')
# library(maptree)
set.seed(1)
dt.fit<-rpart(risk~., data = diffdata,method='class',parms=list(split="information"))
#draw.tree(dt.fit)
printcp(dt.fit)#观察误差
pfit=prune(dt.fit,cp=dt.fit$cptable[which.min(dt.fit$cptable[,"xerror"]),"CP"])#选择xerror最小对应的CP值剪枝
dt.pre<-predict(pfit,newdata=diffindex,type="class")#利用预测集进行预测
(index <- table(actual=diffdata$risk,predict=dt.pre))
(consistency=sum(diag(index))/sum(index))
save(pfit,index,consistency,file='GSE_TCGA/3.decision tree/model.Rdata')

#验证模型model####
#GSE39582####
myfun_km_model <- function(GSE,cli,pairs,model){
  gid=GSE[,1];exp=GSE[,-1]
  # cli=cli[which(cli$ACT=='N'),]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  diffindex=as.data.frame(t(ifelse(diff>0,1,0)))
  colnames(diffindex)=c(paste('pair',1:ncol(diffindex),sep=''))
  model.pre <- predict(model,newdata=diffindex,type="class")
  (index=table(model.pre))
  sample1=rownames(diffindex)[which(model.pre=='L')]
  sample2=rownames(diffindex)[which(model.pre=='H')]
  sam1=cli[match(sample1,cli$sample),]
  sam2=cli[match(sample2,cli$sample),]
  library(plyr)
  sam11=transform(sam1,class='0')
  sam22=transform(sam2,class='1')
  cli_new=rbind(sam11,sam22)
  cli_new=cli_new[!is.na(cli_new$RFS_event),]
  label=cli_new$class
  library(survival)
  survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
  fit=survfit(survi~label)
  cox=coxph(survi~label)
  b=summary(cox)
  hr=round(b$conf.int[1],2)
  lower.hr=round(b$conf.int[3],2)
  upper.hr=round(b$conf.int[4],2)
  cindex=round(b$concordance[1],3)
  p=round(b$sctest[3],3)
  #pdf('GSE_TCGA/km_GSE39582_RFS.pdf')
  plot(fit,col=c('blue','red'),main='GSE39582',ylab='survival rate',
       xlab='survival months',lty=c(1,2),lwd=1.5,mark.time=T)
  legend('bottomright',pch=15:18,lty=c(1,2),merge=F,
         legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                  paste('highrisk',length(which(label==1)),sep=' : ')),
         col=c('blue','red'))
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n','HR=',signif(hr,4),
                        '\n','(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),bty='n',font=2)
  #dev.off()
}
load('GSE_TCGA/data/GSE39582_inf.Rdata')
load('GSE_TCGA/pairs_multicox.Rdata')
myfun_km_model(GSE,cli,pairs_sur,svm.fit)
myfun_km_model(GSE,cli,pairs_sur,nb.fit)
myfun_km_model(GSE,cli,pairs_sur,pfit)

#TCGA&others####
myfun_km_model <- function(GSE,cli,pairs,model,sur_status){
  gid=GSE[,1];exp=GSE[,-1]
  # cli=cli[which(cli$ACT=='N'),]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  diffindex=as.data.frame(t(ifelse(diff>0,1,0)))
  colnames(diffindex)=c(paste('pair',1:ncol(diffindex),sep=''))
  model.pre <- predict(model,newdata=diffindex,type="class")
  (index=table(model.pre))
  sample1=rownames(diffindex)[which(model.pre=='L')]
  sample2=rownames(diffindex)[which(model.pre=='H')]
  sam1=cli[match(sample1,cli$sample),]
  sam2=cli[match(sample2,cli$sample),]
  library(plyr)
  sam11=transform(sam1,class='0')
  sam22=transform(sam2,class='1')
  cli_new=rbind(sam11,sam22)
  if(sur_status==1){
    cli_new=cli_new[!is.na(cli_new$DFS_event),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
  }
  if(sur_status==2){
    cli_new=cli_new[!is.na(cli_new$OS_event),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$OS_time,cli_new$OS_event)
  }
  fit=survfit(survi~label)
  cox=coxph(survi~label)
  b=summary(cox)
  hr=round(b$conf.int[1],2)
  lower.hr=round(b$conf.int[3],2)
  upper.hr=round(b$conf.int[4],2)
  cindex=round(b$concordance[1],3)
  p=signif(b$sctest[3],3)
  #pdf('GSE_TCGA/km.pdf')
  plot(fit,col=c('blue','red'),main='GSE17536',ylab='survival rate',
       xlab='survival months',lty=c(1,2),lwd=1.5,mark.time=T)
  legend('bottomright',pch=15:18,lty=c(1,2),merge=F,
         legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                  paste('highrisk',length(which(label==1)),sep=' : ')),
         col=c('blue','red'))
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n','HR=',signif(hr,4),
                        '\n','(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),bty='n',font=2)
  #dev.off()
}
# load('GSE_TCGA/pairs_multicox.Rdata')
# load('GSE_TCGA/data/TCGA_inf.Rdata')
# myfun_km_model(fpkm,cli,pairs_sur,nb.fit,sur_status=2)
##注意更换图片的题目
myfun_km_model(GSE30378,cli30378,pairs,nb.fit,sur_status=1)


#6.散点图&ROC####
myfun_div <- function(GSE,cli){
  gid=GSE[,1];exp=GSE[,-1]
  cli1=cli[grep('1',cli$stage),]
  exp1=exp[,match(cli1[,1],colnames(exp))]
  cli2=cli[grep('3|4',cli$stage),]
  exp2=exp[,match(cli2[,1],colnames(exp))]
  data=list(gid=gid,exp1=exp1,exp2=exp2)
  return(data)
}

myfun_counts <- function(GSE,cli,pairs){
  data=myfun_div(GSE,cli)
  gid=data$gid;exp1=data$exp1;exp2=data$exp2
  diff1=exp1[match(pairs[,1],gid),]-exp1[match(pairs[,2],gid),]
  count1=colSums(diff1>0)
  count1=sort(count1,decreasing=T)
  diff2=exp2[match(pairs[,1],gid),]-exp2[match(pairs[,2],gid),]
  count2=colSums(diff2>0)
  count2=sort(count2,decreasing=T)
  col_label=c(rep('green',length(count1)),rep('red',length(count2)))
  plot(c(count1,count2),xlab='stage',ylab='counts',main='GSE39582',col=col_label)
  abline(h=nrow(pairs)/2,lwd=1,col='blue')
  sensitivity=mean(count1>=nrow(pairs)/2)
  specificity=mean(count2<nrow(pairs)/2)
  legend('topright',pch=c(21,21),legend=c(paste('sensitivity',signif(sensitivity,digits=3),sep=':'),
                                          paste('specificity',signif(specificity,digits=3),sep=':')),col=c('green','red'))
}
#pdf('VoteK/散点图.pdf')
myfun_counts(GSE,cli,pairs_sur)
#dev.off()

myfun_roc <- function(GSE,cli,pairs){
  data=myfun_div(GSE,cli)
  gid=data$gid;exp1=data$exp1;exp2=data$exp2
  exp = cbind(exp1,exp2)
  sam1 = rep(0,ncol(exp1))#0--control
  sam2 = rep(1,ncol(exp2))#1--case
  label = c(sam1,sam2)
  diff = exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  score = colMeans(diff>0)
  library(pROC)
  roc=roc(label,score,percent=F,ci=T,auc=T)
  plot(roc,type="o",col="cyan2",main='GSE39582',ylim=c(0:1))
  text(x=0.2,y=0.2,paste("AUC=",signif(roc$auc,4),seq=""),bty="n",font=2)
}
#pdf('VoteK/ROC.pdf')
myfun_roc(GSE,cli,pairs_sur)
#dev.off()

#7_1.确定K值投票####
myfun_div <- function(GSE,cli){
  gid=GSE[,1];exp=GSE[,-1]
  cli1=cli[grep('1',cli$stage),]
  exp1=exp[,match(cli1[,1],colnames(exp))]
  cli2=cli[grep('3|4',cli$stage),]
  exp2=exp[,match(cli2[,1],colnames(exp))]
  data=list(gid=gid,exp1=exp1,exp2=exp2)
  return(data)
}

myfun_cindex_k <- function(GSE,cli,pairs){
  gid=GSE[,1];exp=GSE[,-1]
  # cli=cli[which(cli$ACT=='N'),]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  result=c()
  for(k in 0:nrow(pairs)){
    print(k)
    label=rep(0,ncol(diff))
    label[count<k]=1
    cli_new=cbind(cli2,label)
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event) #Surv(sur_time,sur_state)
    cox=coxph(survi~label+location+gender+age+MMR+CIMP+CIN+TP53+KRAS+BRAF,data=cli_new)#survi~label+location+gender+age+MMR+CIMP+CIN+TP53+KRAS+BRAF
    b=summary(cox)
    cindex=b$concordance[1]
    p=b$coefficients[1,5]
    result=rbind(result,c(k,cindex,p))
  }
  colnames(result)=c('k','cindex','p')
  return(result)
}
VoteK=myfun_cindex_k(GSE,cli,pairs_sur)#k=25(>=25-L;<25-H)   #k=20
write.table(VoteK,'VoteK/VoteK96.txt',sep='\t',col.names = T,row.names = F)
# plot(x, y, main, xlab, ylab, xlim, ylim, axes)
# pdf('VoteK/VoteK96.pdf')
plot(VoteK[,1],VoteK[,2],
     xlab='K',ylab='C-index')
lines(VoteK[,1],VoteK[,2],col='red')

# dev.off()

#验证VoteK####
myfun_km_k <- function(GSE,cli,pairs,k=25,sur_status){
  gid=GSE[,1];exp=GSE[,-1]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  sam1=cli[match(names(count)[which(count>=k)],cli$sample),]
  sam2=cli[match(names(count)[which(count<k)],cli$sample),]
  library(plyr)
  sam11=transform(sam1,class='0')
  sam22=transform(sam2,class='1')
  cli_new=rbind(sam11,sam22)
  if(sur_status==1){
    cli_new=cli_new[!is.na(cli_new$DFS_event),]
    cli_new=cli_new[!is.na(cli_new$DFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
  }
  if(sur_status==2){
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    cli_new=cli_new[!is.na(cli_new$RFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
  }
  fit=survfit(survi~label,data=cli_new)
  sur_time=summary(fit)
  # save(sur_time,file='sur_time.Rdata')
  #surv_diff <- survdiff(survi~ label)
  # cox=coxph(survi~label)
  cox=coxph(survi~label,data=cli_new)#survi~label+location+gender+age+MMR+CIMP+CIN+TP53+KRAS+BRAF
  b=summary(cox)
  hr=round(b$conf.int[1],2)
  lower.hr=round(b$conf.int[3],2)
  upper.hr=round(b$conf.int[4],2)
  cindex=round(b$concordance[1],3)
  p=signif(b$sctest[3],3)
  pdf('GPS_GSE.pdf')
  plot(fit,col=c('blue','red'),main='GSE38832',ylab='survival rate',
       xlab='survival months',lty=c(1,2),lwd=1.5,mark.time=T)
  abline(v=60,lwd=1,lty=2,col='black')
  legend('bottomright',pch=15:18,lty=c(1,2),merge=F,
         legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                  paste('highrisk',length(which(label==1)),sep=' : ')),
         col=c('blue','red'))
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n',
                        'C-index=',cindex,'\n',
                        'HR=',signif(hr,4),'\n',
                        '(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),
       bty='n',font=2)
  dev.off()
} #sur_status:DFS-1,RFS-2
myfun_km_k(GSE,cli,pairs_sur,sur_status = 1)#ACT=N
myfun_km_k(GSE,cli,Unipairs,sur_status = 1)

#ggplot2
library("survminer")
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF")
)

#7_2.显著多/少投票####
myfun_km_sig <- function(GSE,cli,pairs,sur_status){
  gid=GSE[,1];exp=GSE[,-1]
  # cli=cli[which(cli$ACT=='N'),]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  sigindex1=qbinom(0.95,nrow(pairs),0.5)
  sigindex2=qbinom(0.05,nrow(pairs),0.5)
  sam1=cli[match(names(count)[which(count>=sigindex1)],cli$sample),]
  sam2=cli[match(names(count)[which(count<=sigindex2)],cli$sample),]
  sam3=cli2[match(setdiff(cli2$sample,union(sam1$sample,sam2$sample)),cli2$sample),]
  library(plyr)
  sam11=transform(sam1,class='0')
  sam22=transform(sam2,class='1')
  sam33=transform(sam3,class='2')
  cli_new=rbind(sam11,sam22,sam33)
  cli_new$class=as.numeric(as.character(cli_new$class))
  if(sur_status==1){
    cli_new=cli_new[!is.na(cli_new$DFS_event),]
    cli_new=cli_new[!is.na(cli_new$DFS_time),]
    label=cli_new$class
    library(survival)
    cli_new2=cli_new[grep('0|1',cli_new$class),]
    label2=cli_new2$class
    survi2=Surv(cli_new2$DFS_time,cli_new2$DFS_event)
    survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
  }
  if(sur_status==2){
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    cli_new=cli_new[!is.na(cli_new$RFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
    cli_new2=cli_new[grep('0|1',cli_new$class),]
    label2=cli_new2$class
    survi2=Surv(cli_new2$RFS_time,cli_new2$RFS_event)
  }
  fit=survfit(survi~label)
  fit2=survfit(survi2~label2)
  sur_time=summary(fit2)
  save(sur_time,file='sur_time.Rdata')
  cox=coxph(survi2~label2)
  b=summary(cox)
  hr=round(b$conf.int[1],2)
  lower.hr=round(b$conf.int[1,3],2)
  upper.hr=round(b$conf.int[1,4],2)
  cindex=round(b$concordance[1],3)
  p=signif(b$sctest[3],3)
  # pdf('VoteSignif_km.pdf')
  plot(fit,col=c('blue','red','gray25'),main='GSE',ylab='survival rate',
       xlab='survival months',lty=c(1,2,3),lwd=1.5,mark.time=T)
  abline(v=60,lwd=1,lty=2,col='black')
  legend('bottomright',pch=15:18,lty=c(1,2,3),merge=F,
         legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                  paste('highrisk',length(which(label==1)),sep=' : '),
                  paste('gray',length(which(label==2)),sep=' : ')),
         col=c('blue','red','gray25'))
  
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n','HR=',signif(hr,4),
                        '\n','(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),bty='n',font=2)
  # dev.off()
} #sur_status:DFS-1,RFS-2,OS-3
myfun_km_sig(GSE,cli,pairs_sur,sur_status = 2)#ACT=N
myfun_km_sig(GSE,cli,Unipairs,sur_status = 1)

#8.构建其他平台基因对标志####
pairs=pairs_sur;gid=GSE[,1]
index=c()
for(i in 1:nrow(pairs)){
  print(i)
  index0=sum(pairs[i,1:2]%in%gid)
  index=c(index,index0)
}
Unipairs=pairs[which(index==2),]
save(Unipairs,file='pairs_GPL96.Rdata')
#多因素cox####
GSE1=GSE;cli1=cli
GSE2=GSE;cli2=cli#ACT
GSE3=GSE;cli3=cli
GSE4=GSE;cli4=cli
covar1=intersect(colnames(cli1),colnames(cli2))
covar2=intersect(colnames(cli3),colnames(cli4))
covar=intersect(covar1,covar2)
cli11=cli1[,c(covar,'RFS_time','RFS_event')]
cli22=cli2[,c(covar,'RFS_time','RFS_event')]
cli33=cli3[,c(covar,'RFS_time','RFS_event')]
cli44=cli4[,c(covar,'DFS_time','DFS_event')]
library(plyr)
cli11=as.matrix(transform(cli11,batch=0))
cli22=as.matrix(transform(cli22,batch=1))
cli33=as.matrix(transform(cli33,batch=2))
cli44=as.matrix(transform(cli44,batch=3))
GSE=cbind(GSE1,GSE2[,-1],GSE3[,-1],GSE4[,-1])
cli=as.data.frame(rbind(cli11,cli22,cli33,cli44),stringsAsFactors = F)
# cli[which(cli$gender=='F'),'gender']=1
# cli[which(cli$gender=='M'),'gender']=2
cli_inf=as.data.frame(apply(cli[,-1], 2, FUN=as.numeric))
cli=transform(cli_inf,sample=cli$sample)
cli=cli[,c(5,1:4)]
cli$sample=as.character(cli$sample)
GSE=GSE[,c(1,match(cli$sample,colnames(GSE)))]
save(GSE,cli,file='GSE_union.Rdata')

myfun_km_k <- function(GSE,cli,pairs,k=25,sur_status){
  gid=GSE[,1];exp=GSE[,-1]
  # cli=cli[which(cli$ACT=='N'),]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  sam1=cli[match(names(count)[which(count>=k)],cli$sample),]
  sam2=cli[match(names(count)[which(count<k)],cli$sample),]
  library(plyr)
  sam11=transform(sam1,class='0')
  sam22=transform(sam2,class='1')
  cli_new=rbind(sam11,sam22)
  if(sur_status==1){
    cli_new=cli_new[!is.na(cli_new$DFS_event),]
    cli_new=cli_new[!is.na(cli_new$DFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
  }
  if(sur_status==2){
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    cli_new=cli_new[!is.na(cli_new$RFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
  }
  if(sur_status==3){
    cli_new=cli_new[!is.na(cli_new$OS_event),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$OS_time,cli_new$OS_event)
  }
  fit=survfit(survi~label,data=cli_new)
  cox=coxph(survi~label+batch,data=cli_new)#label+location+gender+age+batch
  b=summary(cox)
  hr=round(b$conf.int[1,1],2)
  lower.hr=round(b$conf.int[1,3],2)
  upper.hr=round(b$conf.int[1,4],2)
  cindex=round(b$concordance[1],3)
  p=signif(b$coefficients[1,5])
  pdf('VoteK_km_multicox.pdf')
  plot(fit,col=c('blue','red'),main='GSE',ylab='survival rate',
       xlab='survival months',lty=c(1,2),lwd=1.5,mark.time=T)
  legend('bottomright',pch=15:18,lty=c(1,2),merge=F,
         legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                  paste('highrisk',length(which(label==1)),sep=' : ')),
         col=c('blue','red'))
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n','HR=',signif(hr,4),
                        '\n','(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),bty='n',font=2)
  dev.off()
} #sur_status:DFS-1,RFS-2,OS-3
# myfun_km_k(GSE,cli,Unipairs,sur_status = 1)
myfun_km_k(GSE,cli,pairs_sur,sur_status = 2)

myfun_km_sig <- function(GSE,cli,pairs,sur_status){
  gid=GSE[,1];exp=GSE[,-1]
  # cli=cli[which(cli$ACT=='N'),]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  sigindex1=qbinom(0.95,nrow(pairs),0.5)
  sigindex2=qbinom(0.05,nrow(pairs),0.5)
  sam1=cli[match(names(count)[which(count>=sigindex1)],cli$sample),]
  sam2=cli[match(names(count)[which(count<=sigindex2)],cli$sample),]
  sam3=cli2[match(setdiff(cli2$sample,union(sam1$sample,sam2$sample)),cli2$sample),]
  library(plyr)
  sam11=transform(sam1,class='0')
  sam22=transform(sam2,class='1')
  sam33=transform(sam3,class='2')
  cli_new=rbind(sam11,sam22,sam33)
  if(sur_status==1){
    cli_new=cli_new[!is.na(cli_new$DFS_event),]
    cli_new=cli_new[!is.na(cli_new$DFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
  }
  if(sur_status==2){
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    cli_new=cli_new[!is.na(cli_new$RFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
  }
  if(sur_status==3){
    cli_new=cli_new[!is.na(cli_new$OS_event),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$OS_time,cli_new$OS_event)
  }
  fit=survfit(survi~label,data=cli_new)
  cox=coxph(survi~label+batch,data=cli_new)#label+location+gender+age+batch
  b=summary(cox)
  hr=round(b$conf.int[1,1],2)
  lower.hr=round(b$conf.int[1,3],2)
  upper.hr=round(b$conf.int[1,4],2)
  cindex=round(b$concordance[1],3)
  p=signif(b$coefficients[1,5])
  pdf('VoteSignif_km_multicox.pdf')
  plot(fit,col=c('blue','red','gray25'),main='GSE',ylab='survival rate',
       xlab='survival months',lty=c(1,2,3),lwd=1.5,mark.time=T)
  legend('bottomright',pch=15:18,lty=c(1,2,3),merge=F,
         legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                  paste('highrisk',length(which(label==1)),sep=' : '),
                  paste('gray',length(which(label==2)),sep=' : ')),
         col=c('blue','red','gray25'))
  
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n','HR=',signif(hr,4),
                        '\n','(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),bty='n',font=2)
  dev.off()
} #sur_status:DFS-1,RFS-2,OS-3
myfun_km_sig(GSE,cli,pairs_sur,sur_status = 2)

#9.与其他标志的对比####
rm(list=ls())
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/Compare/PRGPs/')
#9.1对比PRGPs####
RiskCluster=read.csv('RiskCluster.txt',sep='\t',header=T,fill=T,as.is=T)
RiskCluster=RiskCluster[,c(1,3,5,7,9,12,14,16)]
colnames(RiskCluster)=c('GSE39582','riskcluster1','TCGA','riskcluster2','GSE14333','riskcluster3','GSE33113','riskcluster4')
RiskCluster[,3]=gsub('-','_',RiskCluster[,3])
#save(RiskCluster,file='RiskCluster.Rdata')
myfun_km <- function(cli,riskcluster,sur_status){
  cli2=cli[grep('2',cli$stage),]
  risksam=riskcluster[match(intersect(cli2$sample,riskcluster[,1]),riskcluster[,1]),]
  save(risksam,file='risksam.Rdata')
  sam1=cli2[match(risksam[grep('Low Risk',risksam[,2]),1],cli2$sample),]
  sam2=cli2[match(risksam[grep('High Risk',risksam[,2]),1],cli2$sample),]
  #library(plyr)
  sam11=transform(sam1,class='0')
  sam22=transform(sam2,class='1')
  cli_new=rbind(sam11,sam22)
  if(sur_status==1){
    cli_new=cli_new[!is.na(cli_new$DFS_event),]
    cli_new=cli_new[!is.na(cli_new$DFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
  }
  if(sur_status==2){
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    cli_new=cli_new[!is.na(cli_new$RFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
  }
  fit=survfit(survi~label,data=cli_new)
  sur_time=summary(fit)
  cox=coxph(survi~label)
  b=summary(cox)
  hr=round(b$conf.int[1],2)
  lower.hr=round(b$conf.int[3],2)
  upper.hr=round(b$conf.int[4],2)
  cindex=round(b$concordance[1],3)
  p=signif(b$sctest[3],3)
  # pdf('PRGPs_GSE.pdf')
  plot(fit,col=c('blue','red'),main='GSE33113',ylab='survival rate',
       xlab='survival months',lty=c(1,2),lwd=1.5,mark.time=T)
  abline(v=60,lwd=1,lty=2,col='black')
  legend('bottomright',pch=15:18,lty=c(1,2),merge=F,
         legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                  paste('highrisk',length(which(label==1)),sep=' : ')),
         col=c('blue','red'))
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n','HR=',signif(hr,4),
                        '\n','(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),bty='n',font=2)
  dev.off()
}
riskcluster=RiskCluster[,7:8]
myfun_km(cli,riskcluster,sur_status=1)
#9.2对比ColoGuideEx,直接比较两套数据合并后结果，结果不好，弃用####
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/Compare/ColoGuideEx/')
sam=read.csv('Ex.txt',sep='\t',header = T,fill=T,as.is=T)
sam_v2=sam[grep('GSE14333|GSE17538',sam[,1]),]
save(sam_v2,file='sam_v2.Rdata')
sam_test=sam[grep('Test',sam$Sample.series),]
save(sam_test,file='sam_test.Rdata')
GSE1=GSE;cli1=cli
GSE2=GSE;cli2=cli
cli=rbind(cli1[,c(1,3,6,7)],cli2[,c(1,2,7,8)])
cli=cli[match(intersect(cli$sample,sam_v2[,2]),cli$sample),]
GSE=cbind(GSE1,GSE2[,-1])
GSE=GSE[,c(1,match(cli$sample,colnames(GSE)))]
save(GSE1,GSE2,cli1,cli2,file = 'GSE_test.Rdata')
#9.2对比ColoGuideEx,自行比对结果####
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/Compare/ColoGuideEx/')
cli12=cli1[grep('2',cli1$stage),]
stage2=cli1[match(intersect(cli12$sample,sam_test$GEOsample),cli1$sample),]
cli214=cli2[grep('1|4',cli2$sample),]
stage14=cli2[match(intersect(cli214$sample,sam_test$GEOsample),cli2$sample),]
GSE14=GSE2[,c(1,match(stage14$sample,colnames(GSE2)))]
GSE2=GSE1[,c(1,match(stage2$sample,colnames(GSE1)))]
# save(GSE2,stage2,GSE14,stage14,file='GSE_inf_test.Rdata')
load("D:/data/TCGA/ensg2gid.Rdata")
symbols=c('PIGR','CXCL13','MMP3','TUBA1B','CXCL10',
       'SESN1','AZGP1','KLK6','EPHA7','SEMA3A','DSC3','ENPP3','BNIP3')
geneid=ensg2gid[match(symbols,ensg2gid$symbol),2]
lgid=geneid[1:5]
hgid=geneid[6:13]

GSE=na.omit(GSE[match(GeneRank$gid,GSE[,1]),])
# x=GSE[,2]
myfun_Ex <- function(x){
  samrank=GSE[order(x,decreasing = F),1]
  lexp=samrank[1:floor(nrow(GSE)*0.2)]
  hexp=samrank[ceiling(0.8*nrow(GSE)):nrow(GSE)]
  score=sum(sum(lgid%in%lexp),sum(hgid%in%hexp))
  # if(score>=5) {'highrisk'}
  # else {'lowrisk'}
  ifelse(score>=5,'highrisk','lowrisk')
}
cluster=apply(GSE[,-1],2,myfun_Ex)
table(cluster)
#9.3对比IRGPI####
rm(list=ls())
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/Compare/IRGPI')
cluster=read.csv('cluster.txt',sep='\t',header=T,fill=T,as.is=T)
save(cluster,file='cluster.Rdata')
#km_gps#####
myfun_km_k <- function(GSE,cli,pairs,GEO,k=25,sur_status){
  gid=GSE[,1];exp=GSE[,-1]
  cli=cli[grep('2',cli$stage),]
  GEOsamples=cluster[grep(GEO,cluster[,2]),]
  cli=cli[match(intersect(cli$sample,GEOsamples[,1]),cli$sample),]
  exp=exp[,match(cli$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  sam1=cli[match(names(count)[which(count>=k)],cli$sample),]
  sam2=cli[match(names(count)[which(count<k)],cli$sample),]
  library(plyr)
  sam11=transform(sam1,class='0')
  sam22=transform(sam2,class='1')
  cli_new=rbind(sam11,sam22)
  if(sur_status==1){
    cli_new=cli_new[!is.na(cli_new$DFS_event),]
    cli_new=cli_new[!is.na(cli_new$DFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
  }
  if(sur_status==2){
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    cli_new=cli_new[!is.na(cli_new$RFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
  }
  fit=survfit(survi~label,data=cli_new)
  sur_time=summary(fit)
  # save(sur_time,file='sur_time.Rdata')
  #surv_diff <- survdiff(survi~ label)
  cox=coxph(survi~label)
  b=summary(cox)
  hr=round(b$conf.int[1],2)
  lower.hr=round(b$conf.int[3],2)
  upper.hr=round(b$conf.int[4],2)
  cindex=round(b$concordance[1],3)
  p=signif(b$sctest[3],3)
  # pdf('GPS_GSE.pdf')
  plot(fit,col=c('blue','red'),main='GSE',ylab='survival rate',
       xlab='survival months',lty=c(1,2),lwd=1.5,mark.time=T)
  abline(v=60,lwd=1,lty=2,col='black')
  legend('bottomright',pch=15:18,lty=c(1,2),merge=F,
         legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                  paste('highrisk',length(which(label==1)),sep=' : ')),
         col=c('blue','red'))
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n',
                        'C-index=',cindex,'\n',
                        'HR=',signif(hr,4),'\n',
                        '(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),
       bty='n',font=2)
  # dev.off()
} #sur_status:DFS-1,RFS-2
myfun_km_k(GSE,cli,pairs_sur,GEO='GSE',sur_status = 1)#ACT=N

#km_IRGPI#####
myfun_km_k <- function(GSE,cli,GEO,sur_status){
  gid=GSE[,1];exp=GSE[,-1]
  cli=cli[grep('2',cli$stage),]
  GEOsamples=cluster[grep(GEO,cluster[,2]),]
  cli=cli[match(intersect(cli$sample,GEOsamples[,1]),cli$sample),]
  samcluster=GEOsamples[match(intersect(cli$sample,GEOsamples[,1]),GEOsamples[,1]),]
  sam1=samcluster[grep('Low Immune Risk',samcluster$immune.risk.group),1]
  sam2=samcluster[grep('High Immune Risk',samcluster$immune.risk.group),1]
  sam11=transform(cli[match(sam1,cli$sample),],class='0')
  sam22=transform(cli[match(sam2,cli$sample),],class='1')
  cli_new=rbind(sam11,sam22)
  if(sur_status==1){
    cli_new=cli_new[!is.na(cli_new$DFS_event),]
    cli_new=cli_new[!is.na(cli_new$DFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
  }
  if(sur_status==2){
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    cli_new=cli_new[!is.na(cli_new$RFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
  }
  fit=survfit(survi~label,data=cli_new)
  sur_time=summary(fit)
  # save(sur_time,file='sur_time.Rdata')
  #surv_diff <- survdiff(survi~ label)
  cox=coxph(survi~label)
  b=summary(cox)
  hr=round(b$conf.int[1],2)
  lower.hr=round(b$conf.int[3],2)
  upper.hr=round(b$conf.int[4],2)
  cindex=round(b$concordance[1],3)
  p=signif(b$sctest[3],3)
  # pdf('IRGPI_GSE.pdf')
  plot(fit,col=c('blue','red'),main='GSE17536',ylab='survival rate',
       xlab='survival months',lty=c(1,2),lwd=1.5,mark.time=T)
  abline(v=60,lwd=1,lty=2,col='black')
  legend('bottomright',pch=15:18,lty=c(1,2),merge=F,
         legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                  paste('highrisk',length(which(label==1)),sep=' : ')),
         col=c('blue','red'))
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n',
                        'C-index=',cindex,'\n',
                        'HR=',signif(hr,4),'\n',
                        '(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),
       bty='n',font=2)
  # dev.off()
} #sur_status:DFS-1,RFS-2
myfun_km_k(GSE,cli,GEO='GSE17536',sur_status = 1)#ACT=N

#9.4对比oncotypeDx colon####
rm(list=ls())
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/Compare/OncotypeDx/')
load("D:/data/TCGA/ensg2gid.Rdata")
myfun_Dx <- function(x){
  gene1=c('BGN','FAP','INHBA')
  gene2=c('MKI67','MYC','MYBL2')
  gene3='GADD45B'
  refergene=c('ATP5E','GPX1','PGK1','VDAC2','UBB')
  gene1=ensg2gid[match(gene1,ensg2gid$symbol),2]
  gene2=ensg2gid[match(gene2,ensg2gid$symbol),2]
  gene3=ensg2gid[match(gene3,ensg2gid$symbol),2]
  refergene=ensg2gid[match(refergene,ensg2gid$symbol),2]
  refer=mean(x[match(refergene,GSE[,1])])
  RSu=0.15*mean(x[match(gene1,GSE[,1])])-0.3*mean(x[match(gene2,GSE[,1])])+0.15*(x[match(gene3,GSE[,1])])
  RS=44*(RSu+0.82)
  return(RS)
}
# RS=apply(exp,2,myfun_Dx)
# RSu = 0·15* mean(BGN,FAP,INHBA)-0·3*mean(MKI67,MYC,MYBL2)+ 0·15* GADD45B)
# This score was then rescaled RS=44*(RSu+0.82)

myfun_km_Dx <- function(GSE,cli,sur_status){
  gid=GSE[,1];exp=GSE[,-1]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  RS=apply(exp,2,myfun_Dx)
  sam1=cli[match(names(RS)[which(RS<30)],cli$sample),]
  sam2=cli[match(names(RS)[which(RS>=41)],cli$sample),]
  library(plyr)
  sam11=transform(sam1,class='0')
  sam22=transform(sam2,class='1')
  cli_new=rbind(sam11,sam22)
  if(sur_status==1){
    cli_new=cli_new[!is.na(cli_new$DFS_event),]
    cli_new=cli_new[!is.na(cli_new$DFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
  }
  if(sur_status==2){
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    cli_new=cli_new[!is.na(cli_new$RFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
  }
  fit=survfit(survi~label,data=cli_new)
  sur_time=summary(fit)
  # save(sur_time,file='sur_time.Rdata')
  #surv_diff <- survdiff(survi~ label)
  cox=coxph(survi~label)
  b=summary(cox)
  hr=round(b$conf.int[1],2)
  lower.hr=round(b$conf.int[3],2)
  upper.hr=round(b$conf.int[4],2)
  cindex=round(b$concordance[1],3)
  p=signif(b$sctest[3],3)
  # pdf('Dx_GSE.pdf')
  plot(fit,col=c('blue','red'),main='GSE',ylab='survival rate',
       xlab='survival months',lty=c(1,2),lwd=1.5,mark.time=T)
  abline(v=60,lwd=1,lty=2,col='black')
  legend('bottomright',pch=15:18,lty=c(1,2),merge=F,
         legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                  paste('highrisk',length(which(label==1)),sep=' : ')),
         col=c('blue','red'))
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n',
                        'C-index=',cindex,'\n',
                        'HR=',signif(hr,4),'\n',
                        '(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),
       bty='n',font=2)
  # dev.off()
} #sur_status:DFS-1,RFS-2
myfun_km_Dx(fpkm,cli,sur_status = 1)#ACT=N

#9.5对比13genes####
rm(list=ls())
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/Compare/13genes/')
Risk score = (0.040757*THBS2)+(0.077189*CAV2)+(0.047062*SCG2)+(0.252962*SLC6A1)+
  (0.306452*SAV1)+(-0.06446*MRPL35) +(0.076798*SEZ6L2) +(0.407432*ERO1A)+
  (0.22306*RAB3B) +(0.49833*OBSL1)+ (0.038808*CD109) + (0.015854*PTPN14) + 
  (-0.04333*LRPAP1)
#9.6对比3GPS####
rm(list=ls())
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/Compare/3GPS')
gps=read.table('3GPS.txt',sep='\t',header=T,fill=T,as.is=T)
load("D:/data/TCGA/ensg2gid.Rdata")
gps[,1]=ensg2gid[match(gps[,1],ensg2gid$symbol),2]
gps[,2]=ensg2gid[match(gps[,2],ensg2gid$symbol),2]
save(gps,file='gps.Rdata')
myfun_sur <- function(GSE,cli,pairs,sur_status){
  gid=GSE[,1];exp=GSE[,-1]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  exp=exp[,match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  sam1=cli[match(names(count)[which(count>=2)],cli$sample),]
  sam2=cli[match(names(count)[which(count<2)],cli$sample),]
  library(plyr)
  sam11=transform(sam1,class='0')
  sam22=transform(sam2,class='1')
  cli_new=rbind(sam11,sam22)
  if(sur_status==1){
    cli_new=cli_new[!is.na(cli_new$DFS_event),]
    cli_new=cli_new[!is.na(cli_new$DFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
  }
  if(sur_status==2){
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    cli_new=cli_new[!is.na(cli_new$RFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
  }
  fit=survfit(survi~label,data=cli_new)
  sur_time=summary(fit)
  cox=coxph(survi~label,data=cli_new)
  b=summary(cox)
  hr=round(b$conf.int[1],2)
  lower.hr=round(b$conf.int[3],2)
  upper.hr=round(b$conf.int[4],2)
  cindex=round(b$concordance[1],3)
  p=signif(b$sctest[3],3)
  pdf('GPS_GSE.pdf')
  plot(fit,col=c('blue','red'),main='GSE33113',ylab='Disease-free survival',
       xlab='Disease-free survival (months)',lty=c(1,1),lwd=2,mark.time=T)
  legend('bottomright',pch=15:18,lty=c(1,1),merge=F,
         legend=c(paste('Low-risk',length(which(label==0)),sep=' : '),
                  paste('High-risk',length(which(label==1)),sep=' : ')),
         col=c('blue','red'))
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n',
                        'C-index=',cindex,'\n',
                        'HR=',signif(hr,4),'\n',
                        '(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),
       bty='n',font=2)
  dev.off()
}
myfun_sur(GSE,cli,gps,sur_status=2)

#cindex条形图对比结果####
rm(list=ls())
# install.packages('reshape2')
# library(reshape2)
data=data.frame(class=c('51-GPS','PRGPI','Oncotype Dx'),
                GSE39582=c(0.627,0.692,0.516),
                GSE14333=c(0.639,0.779,0.494),
                GSE17538=c(0.71,0.661,0.494),
                GSE33113=c(0.656,0.608,0.583),
                TCGA=c(0.639,0.564,0.518))
mydata=as.matrix(data[,-1])
row.names(mydata)=data[,1]
# pdf('Three Signatures C-index.pdf')
barplot(mydata, col=colors()[c(23,89,12)], border='white', font.axis=2, beside=T, 
        legend=rownames(mydata), ylab='C-index', font.lab=2 , ylim = 0:1, xpd=F)
dev.off()
####################################
#######TCGA多组学验证结果###########
####################################
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/multi-omics/')
#####
#1.edgeR####
rm(list=ls())
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/enrichment/')
# source("http://bioconductor.org/biocLite.R")
# BiocManager::install("edgeR")
library(edgeR)
myfun_pre <- function(GSE,cli,pairs,k=25){
  gid=GSE[,1];exp=GSE[,-1]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  sam_l=cli[match(names(count)[which(count>=k)],cli$sample),]
  sam_h=cli[match(names(count)[which(count<k)],cli$sample),]
  cli_pre=list(sam_l=sam_l,sam_h=sam_h)
  return(cli_pre)
}
cli_pre=myfun_pre(fpkm,cli,pairs_sur)

myfun_edgeR <- function(cli_pre,counts){
  sam_l=cli_pre$sam_l;sam_h=cli_pre$sam_h
  gid=counts$gid
  exp_l=counts[,match(sam_l$sample,colnames(counts))]
  exp_h=counts[,match(sam_h$sample,colnames(counts))]
  # colnames(exp_l)=paste('L',1:ncol(exp_l),sep='')
  # colnames(exp_h)=paste('H',1:ncol(exp_h),sep='')
  exp=cbind(exp_l,exp_h)
  rownames(exp)=gid
  group=factor(c(rep('L',ncol(exp_l)),rep('H',ncol(exp_h))))#以首字母靠前的为参照
  data1=DGEList(counts=exp,genes=gid,group=group)#创建edgeR数据格式
  index=rowSums(cpm(data1)>1)>=(ncol(exp)/2)#过滤掉在至少50%的样本中<=1的基因
  data2=data1[index,]
  data2=calcNormFactors(data2)#标准化，默认为TMN
  data2=estimateCommonDisp(data2)
  data2=estimateTagwiseDisp(data2)
  et2=exactTest(data2)
  fdr=p.adjust(et2$table[,3],method="BH")
  et2$table[,1]=-et2$table[,1]#以低风险做参考
  gid=as.numeric(rownames(et2$table))
  edgeR_result=cbind(gid,et2$table,fdr)
  DEGs=edgeR_result[which(edgeR_result$fdr<0.05),]
  return(DEGs)
}
DEGs=myfun_edgeR(cli_pre,counts)
save(DEGs,file='DEGs_TCGA.Rdata')

#limma####
rm(list=ls())
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("limma")
myfun_pre <- function(GSE,cli,pairs,k=25){
  gid=GSE[,1];exp=GSE[,-1]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  sam_l=cli[match(names(count)[which(count>=k)],cli$sample),]
  sam_h=cli[match(names(count)[which(count<k)],cli$sample),]
  cli_pre=list(sam_l=sam_l,sam_h=sam_h)
  return(cli_pre)
}
cli_pre=myfun_pre(GSE,cli,pairs_sur)#Unipairs;cli=cli[grep('Primary tumor resection',cli$tissue),]
library(limma)
myfun_limma <- function(cli_pre,GSE){
  sam_l=cli_pre$sam_l;sam_h=cli_pre$sam_h
  gid=GSE[,1]
  exp_l=GSE[,match(sam_l$sample,colnames(GSE))]
  exp_h=GSE[,match(sam_h$sample,colnames(GSE))]
  exp=cbind(exp_l,exp_h)
  rownames(exp)=gid
  label=factor(c(rep('control',ncol(exp_l)),rep('case',ncol(exp_h))))
  design=model.matrix(~0+label)
  colnames(design)=c('control','case')
  fit=lmFit(exp,design)
  cont.matrix=makeContrasts('control-case',levels=design)#case相对于control
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  final=topTable(fit2,coef=1,number=nrow(exp),adjust.method='BH',sort.by='B',resort.by='M')
  limmaResult=cbind(as.numeric(rownames(final)),final)
  colnames(limmaResult)=c('gid','logFC','AveExpr','t','p','fdr','B')
  limmaResult=limmaResult[order(limmaResult$gid,decreasing = F),]
  return(limmaResult)
}#上下调方向指exp_h相对于exp_l
limmaResult <- myfun_limma(cli_pre,GSE)
DEGs=limmaResult[which(limmaResult$fdr<0.5),]
save(DEGs,file='DEGs_limma.Rdata')

#交叠一致性####
rm(list=ls())
setwd('../差异分析一致性结果/')
tcga=DEGs
cogenes=intersect(tcga$gid,DEGs$gid)
consistency=mean((tcga[match(cogenes,tcga$gid),2]*DEGs[match(cogenes,DEGs$gid),2])>0)
library (VennDiagram)
venn.diagram(x=list(TCGA=tcga$gid,GSE=DEGs$gid),
             filename = "TCGA_GSE.png",
             main=paste('Consistency = ',signif(consistency,3),sep=''),
             main.cex=0.7,
             height=450,width = 450,resolution =300,imagetype="png",col="transparent",
             fill=c('cornflowerblue','darkorchid1'),
             alpha = 0.50, cex=0.45, cat.cex=0.4,margin = 0.2)

#富集分析####
# source("http://www.bioconductor.org/biocLite.R")
# BiocManager::install("AnnotationHub")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("clusterProfiler")
# BiocManager::install("topGO")
# BiocManager::install("Rgraphviz")
library(AnnotationHub)
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
id=DEGs$gid
#KEGG分析#
ekk <- enrichKEGG(id,organism='hsa',pvalueCutoff=0.05)#KEGG富集分析
kegg=as.data.frame(ekk)
pdf('Bubble plot_Enrichment KEGG.pdf')
dotplot(ekk,font.size=12,showCategory=10,title='Bubble plot_Enrichment KEGG')#画气泡图
dev.off()
write.table(kegg,"KEGG-enrich.txt",sep='\t',row.names =F)
# browseKEGG(ekk,'hsa00190')	#显示通路图

#GO富集分析#
ego <- enrichGO(OrgDb="org.Hs.eg.db",gene=id,keyType = "ENTREZID",
                ont='BP',pvalueCutoff=0.05,readable=F) #GO富集分析;ont='BP/MF/CC'
go=as.data.frame(ego)
write.table(go,"GO-enrich_BP.txt",sep='\t',row.names =F)
pdf('Enrichment GO_BP.pdf')
dotplot(ego,showCategory=10,title="Enrichment GO Top10") #泡泡图
barplot(ego, showCategory=20,title="EnrichmentGO")  #柱状图
dev.off()
# plotGOgraph(ego) 	#GO图

#GSEA数据准备######
rm(list=ls())
myfun_pre <- function(GSE,cli,pairs,k=25){
  gid=GSE[,1];exp=GSE[,-1]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  sam_l=cli[match(names(count)[which(count>=k)],cli$sample),]
  sam_h=cli[match(names(count)[which(count<k)],cli$sample),]
  cli_pre=list(sam_l=sam_l,sam_h=sam_h)
  return(cli_pre)
}
cli_pre=myfun_pre(fpkm,cli,pairs_sur)

sam_l=cli_pre$sam_l;sam_h=cli_pre$sam_h
gid=fpkm$gid
exp_l=fpkm[,match(sam_l$sample,colnames(fpkm))]
exp_h=fpkm[,match(sam_h$sample,colnames(fpkm))]
exp=cbind(exp_l,exp_h)
data=cbind(gid,exp)
subdata=data[match(DEGs$gid,gid),]#导入DEGs
subdata1=na.omit(subdata)
a=ensg2gid[match(subdata1$gid,ensg2gid$entrez_id),3]#导入ensg2gid
sep=cbind(a,subdata1)
colnames(sep)[1:2]=c('symbol','gene')
write.table(sep,'DEGs_exp.txt',sep='\t',col.names = T,row.names = F)
label=c(rep(0,ncol(exp_l)),rep(1,ncol(exp_h)))
write.table(label,'label.txt',sep='\t',col.names = T,row.names = F)

#2.mutation#####
rm(list=ls())
load('cli_pre.Rdata')
mcoad=read.table('D:/data/TCGA/COAD/COAD_mutation/TCGA.COAD.maf',sep='\t',header=T,fill=T)
mread=read.table('D:/data/TCGA/READ/READ_mutation/TCGA.READ.maf',sep='\t',header=T,fill=T,quote = '')
m1=mcoad[,c(2,9,16)]
m2=mread[,c(2,9,16)]
mdata=rbind(m1,m2)
mdata[,2:3]=apply(mdata[,2:3],2,as.character)
mdata=mdata[-grep('Silent',mdata$Variant_Classification),]
mdata[,3]=substr(mdata$Tumor_Sample_Barcode,1,16)
mdata=mdata[grep('01A$',mdata[,3]),]
mdata[,3]=gsub('-','_',mdata[,3])
mdata[,3]=substr(mdata$Tumor_Sample_Barcode,1,12)
save(mdata,file='TCGA_mutation.Rdata')
#导入突变信息，分组信息，ensg2gid
gid=unique(mdata$Entrez_Gene_Id)
fresult=c()
for(i in 1:length(gid)){
  print(i)
  samples=unique(mdata[which(gid[i]==mdata$Entrez_Gene_Id),3])
  lM=sum(samples%in%sam_l$sample)
  hM=sum(samples%in%sam_h$sample)
  ftest=fisher.test(matrix(c(lM,(nrow(sam_l)-lM),hM,(nrow(sam_h)-hM)),2,2,byrow = T))
  fd=lM/nrow(sam_l)-hM/nrow(sam_h)
  fresult=rbind(fresult,c(gid[i],fd,ftest$estimate,ftest$p.value))
}
# fresult=cbind(fresult,p.adjust(fresult[,3],method = 'BH'))#因为突变频率较低，卡p值即可
fresult=cbind(fresult,ensg2gid[match(fresult[,1],ensg2gid$entrez_id),3])
fresult=fresult[,c(1,5,2:4)]
colnames(fresult)=c('gid','symbol','fd','estimate','p')#fd值为低风险相对于高风险突变率
fresult1=as.data.frame(fresult)
fresult1[,c(1,3:5)]=apply(fresult1[,c(1,3:5)],2,as.numeric)
fresult1[,2]=as.character(fresult1[,2])
fresult=fresult1[order(fresult1[,1],decreasing = F),]
SigMgenes=fresult[fresult$p<0.05,]
save(SigMgenes,fresult,file='mutation_fisher.Rdata')

#3.CNV####
rm(list=ls())
cnv=read.table('D:/data/TCGA/cnv/COADREAD-TP.CopyNumber/all_lesions.conf_99.txt',sep='\t',header=T,fill=T,as.is=T)
cnv1=cnv[1:72,c(2,10:625)]
colnames(cnv1)=substr(colnames(cnv1),1,16)
cnv2=cnv1[,grep('01A$',colnames(cnv1))]
colnames(cnv2)=gsub('[.]','_',colnames(cnv2))
colnames(cnv2)=substr(colnames(cnv2),1,12)
cnv=cbind(cnv1[,1],cnv2)
cnv[,1]=as.character(cnv[,1])
save(cnv,file='TCGA_cnv.Rdata')

load('cli_pre.Rdata')
lsam=intersect(colnames(cnv),sam_l$sample)
hsam=intersect(colnames(cnv),sam_h$sample)
fresult=c()
for(i in 1:nrow(cnv)){
  print(i)
  l0=sum(cnv[i,match(lsam,colnames(cnv))]==0)
  h0=sum(cnv[i,match(hsam,colnames(cnv))]==0)
  ftest=fisher.test(matrix(c(l0,(length(lsam)-l0),h0,(length(hsam)-h0)),2,2,byrow = F))
  fresult=rbind(fresult,c(cnv[i,1],ftest$estimate,ftest$p.value))
}
# fresult1=cbind(fresult,p.adjust(fresult[,3],method = 'BH'))
colnames(fresult)=c('segment','estimate','p')
fresult1=as.data.frame(fresult,stringsAsFactors = F)
fresult1[,c(2:3)]=apply(fresult1[,c(2:3)],2,as.numeric)
fresult1[,1]=gsub(' ','',fresult1[,1])
SigCNV=fresult1[fresult1$p<0.05,]
save(SigCNV,fresult,file='CNV_fisher.Rdata')

amp=read.table('D:/data/TCGA/cnv/COADREAD-TP.CopyNumber/amp_genes.conf_99.txt',sep='\t',header = F,fill=T,stringsAsFactors = F)
del=read.table('D:/data/TCGA/cnv/COADREAD-TP.CopyNumber/del_genes.conf_99.txt',sep='\t',header = F,fill=T,stringsAsFactors = F)
segAmp=amp[,match(SigCNV$segment[1:3],amp[1,])]
segDel=del[,match(SigCNV$segment[4:6],del[1,])]
colnames(segAmp)=segAmp[1,];segAmp=segAmp[-(1:4),]
colnames(segDel)=segDel[1,];segDel=segDel[-(1:4),]
save(segAmp,segDel,SigCNV,file='SigCNV.Rdata')

#4.methylation####
rm(list=ls())
data=read.table('D:/data/TCGA/SK/Methylation27/pcexp.txt',sep='\t',header=T,fill=T)
gid=as.character(data$probe)
setwd('D:/data/TCGA/READ/READ_methylation/')
library(data.table)
time_fread <- system.time(
  methydata <- fread('Illumina Human Methylation 450.merge.txt')
)
paste("数据的大小为：",format(object.size(methydata),units="auto"))## 数据的大小
methyread=methydata[match(gid,methydata$Tags),]
setwd('D:/data/TCGA/COAD/COAD_methylation/')
time_fread <- system.time(
  methydata <- fread('Illumina Human Methylation 450.merge.txt')
)
paste("数据的大小为：",format(object.size(methydata),units="auto"))## 数据的大小
methycoad=methydata[match(gid,methydata$Tags),]
methycrc=cbind(methycoad,methyread[,-1])
methycrc=na.omit(methycrc)
colnames(methycrc)=substr(colnames(methycrc),1,12)
colnames(methycrc)=gsub('-','_',colnames(methycrc))
save(methycrc,file='TCGA_methylation.Rdata')

myfun_limma <- function(cli_pre,GSE){
  gid=GSE[,1]
  exp_l=GSE[,match(intersect(sam_l$sample,colnames(GSE)),colnames(GSE))]
  exp_h=GSE[,match(intersect(sam_h$sample,colnames(GSE)),colnames(GSE))]
  exp=cbind(exp_l,exp_h)
  rownames(exp)=gid
  label=factor(c(rep('control',ncol(exp_l)),rep('case',ncol(exp_h))))
  design=model.matrix(~0+label)
  colnames(design)=c('control','case')
  fit=lmFit(exp,design)
  cont.matrix=makeContrasts('control-case',levels=design)#case相对于control
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  final=topTable(fit2,coef=1,number=nrow(exp),adjust.method='BH',sort.by='B',resort.by='M')
  limmaResult=cbind(as.character(rownames(final)),final)
  colnames(limmaResult)=c('gid','logFC','AveExpr','t','p','fdr','B')
  limmaResult=limmaResult[order(limmaResult$gid,decreasing = F),]
  return(limmaResult)
}#上下调方向指exp_h相对于exp_l
limmaResult <- myfun_limma(cli_pre,methycrc)
limmaResult$gid=as.character(limmaResult$gid)
DEprobe=limmaResult[which(limmaResult$fdr<0.05),]
save(limmaResult,DEprobe,file='methylation_limma.Rdata')

pid2gid=read.table('D:/data/TCGA/SK/Methylation27/Level_3/jhu-usc.edu_COAD.HumanMethylation27.1.lvl-3.TCGA-A6-2670-01A-02D-0820-05.txt',sep='\t',header = T,fill=T,stringsAsFactors = F)
colnames(pid2gid)=pid2gid[1,]
pid2gid=pid2gid[-1,c(1,3)]
DEprobe=cbind(pid2gid[match(DEprobe$gid,pid2gid$`Composite Element REF`),2],DEprobe)
colnames(DEprobe)[1]='symbol'
DEprobe[,1]=as.character(DEprobe[,1])
data1=DEprobe[grep('^$',DEprobe$symbol),]
data2=DEprobe[-grep('^$',DEprobe$symbol),]
sep1=table(data2$symbol)
data21=data2[match(names(sep1)[which(sep1==1)],data2$symbol),]
data22=data2[which(data2$symbol%in%names(sep1)[which(sep1!=1)]),]
myfun_all <- function(x){
  ifelse(all(x>0)|all(x<0),1,-1)
}
sep2=apply(data22[,3,drop=F],2,function(x) tapply(x,as.factor(data22$symbol),myfun_all))
# sep=apply(data22[,3,drop=F],2,function(x) tapply(x,as.factor(data22$symbol),prod))#prod表示连乘，使用错误
data220=data22[match(rownames(sep2)[which(sep2[,1]>0)],data22$symbol),]
DEGs=rbind(data1,data21,data220)
DEGs=cbind(ensg2gid[match(DEGs$symbol,ensg2gid$symbol),2],DEGs)#需要导入ensg2gene
colnames(DEGs)[c(1,3)]=c('gid','segment')
save(DEGs,file='SigMethylation.Rdata')


#多因素Cox风险分析####
gid=GSE[,1];exp=GSE[,-1];k=25
cli2=cli[grep('2',cli$stage),]
exp=exp[,match(cli2$sample,colnames(exp))]
diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
count=colSums(diff>0)
sam1=cli[match(names(count)[which(count>=k)],cli$sample),]
sam2=cli[match(names(count)[which(count<k)],cli$sample),]
library(plyr)
sam11=transform(sam1,class='0')
sam22=transform(sam2,class='1')
cli_new=rbind(sam11,sam22)
cli_new=cli_new[!is.na(cli_new$DFS_event),]
cli_new=cli_new[!is.na(cli_new$DFS_time),]
label=as.numeric(cli_new$class)
library(survival)
survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
# fit=survfit(survi~label,data=cli_new)
# sur_time=summary(fit)
#surv_diff <- survdiff(survi~ label)
cox=coxph(survi~label+location+gender+age+MMR+CIMP+CIN+TP53+KRAS+BRAF,data = cli_new)
b=summary(cox)
hr=round(b$conf.int[1],2)
lower.hr=round(b$conf.int[3],2)
upper.hr=round(b$conf.int[4],2)
cindex=round(b$concordance[1],3)
p=signif(b$sctest[3],3)

####################################
######TCGA后续CellLines验证实验#####
####################################
rm(list=ls())
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/CellLines/')
####
#1.TCGA_stage2高低风险样本间差异基因edgeR####
rm(list=ls())
# source("http://bioconductor.org/biocLite.R")
# BiocManager::install("edgeR")
library(edgeR)
load('data/TCGA_inf.Rdata')
load('pairs_multicox.Rdata')
myfun_pre <- function(GSE,cli,pairs,k=25){
  gid=GSE[,1];exp=GSE[,-1]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  sam_l=cli[match(names(count)[which(count>=k)],cli$sample),]
  sam_h=cli[match(names(count)[which(count<k)],cli$sample),]
  cli_pre=list(sam_l=sam_l,sam_h=sam_h)
  return(cli_pre)
}
cli_pre=myfun_pre(fpkm,cli,pairs_sur)

myfun_edgeR <- function(cli_pre,counts){
  sam_l=cli_pre$sam_l;sam_h=cli_pre$sam_h
  gid=counts$gid
  exp_l=counts[,match(sam_l$sample,colnames(counts))]
  exp_h=counts[,match(sam_h$sample,colnames(counts))]
  exp=cbind(exp_l,exp_h)
  rownames(exp)=gid
  group=factor(c(rep('L',ncol(exp_l)),rep('H',ncol(exp_h))))#以首字母靠前的为参照
  data1=DGEList(counts=exp,genes=gid,group=group)#创建edgeR数据格式
  index=rowSums(cpm(data1)>1)>=(ncol(exp)/2)#过滤掉在至少50%的样本中<=1的基因
  data2=data1[index,]
  data2=calcNormFactors(data2)#标准化，默认为TMN
  data2=estimateCommonDisp(data2)
  data2=estimateTagwiseDisp(data2)
  et2=exactTest(data2)
  fdr=p.adjust(et2$table[,3],method="BH")
  et2$table[,1]=-et2$table[,1]#以低风险做参考,本身以字母小的做参考
  gid=as.numeric(rownames(et2$table))
  edgeR_result=cbind(gid,et2$table,fdr)
  DEGs=edgeR_result[which(edgeR_result$fdr<0.05),]
  return(DEGs)
}
DEGs=myfun_edgeR(cli_pre,counts)
save(DEGs,file='DEGs_TCGA.Rdata')
#2.组织数据交叠一致性####
hexp2=DEGs[which(DEGs$logFC>0),]
hexpstage=DEGs[which(DEGs$logFC>0),]
cohexp=intersect(hexp2$gid,hexpstage$gid)
hexp=cbind(hexp2[match(cohexp,hexp2$gid),],hexpstage[match(cohexp,hexpstage$gid),])
#3.划分CellLines数据差异基因####
myfun_pre <- function(GSE,pairs,k=25){
  gid=GSE[,1];exp=GSE[,-1]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  sam_l=names(count)[which(count>=k)]
  sam_h=names(count)[which(count<k)]
  cli_pre=list(sam_l=sam_l,sam_h=sam_h)
  return(cli_pre)
}
cli_pre=myfun_pre(rpkm,pairs_sur)
# save(cli_pre,file='CellLines_pre.Rdata')
myfun_edgeR <- function(cli_pre,counts){
  sam_l=cli_pre$sam_l;sam_h=cli_pre$sam_h
  gid=counts[,1]
  exp_l=counts[,match(sam_l,colnames(counts))]
  exp_h=counts[,match(sam_h,colnames(counts))]
  exp=cbind(exp_l,exp_h)
  rownames(exp)=gid
  group=factor(c(rep('L',ncol(exp_l)),rep('H',ncol(exp_h))))#以首字母靠前的为参照
  data1=DGEList(counts=exp,genes=gid,group=group)#创建edgeR数据格式
  index=rowSums(cpm(data1)>1)>=(ncol(exp)/2)#过滤掉在至少50%的样本中<=1的基因
  data2=data1[index,]
  data2=calcNormFactors(data2)#标准化，默认为TMN
  data2=estimateCommonDisp(data2)
  data2=estimateTagwiseDisp(data2)
  et2=exactTest(data2)
  fdr=p.adjust(et2$table[,3],method="BH")
  et2$table[,1]=-et2$table[,1]#以低风险做参考,本身以字母小的做参考
  gid=as.numeric(rownames(et2$table))
  edgeR_result=cbind(gid,et2$table,fdr)
  DEGs=edgeR_result[which(edgeR_result$fdr<0.05),]
  return(DEGs)
}
DEGs=myfun_edgeR(cli_pre,counts)
save(DEGs,file='CellLines_DEGs.Rdata')
#4.基因一致性####
TargetGenes=read.csv('D:/data/DrugBank/drugid2geneid.txt',sep='\t',header=T,fill=T,as.is=T)
Tgenes=unique(TargetGenes$geneid)
Drugid=read.csv('D:/data/DrugBank/drug2drugid.txt',sep='\t',header=T,fill=T,as.is=T)
DEGs_CellLines=DEGs
DEGs_stage2=DEGs
DEGs_stage134=DEGs
coDEGs=intersect(DEGs_CellLines$gid,DEGs_stage2$gid)
index=DEGs_CellLines[match(coDEGs,DEGs_CellLines$gid),2]*DEGs_stage2[match(coDEGs,DEGs_stage2$gid),2]
sum(index>0)
CellLines=DEGs_CellLines[match(coDEGs[which(index>0)],DEGs_CellLines$gid),]
stage2=DEGs_stage2[match(coDEGs[which(index>0)],DEGs_stage2$gid),]
DEGs=cbind(CellLines,stage2)
hDEGs=DEGs[which(DEGs[,2]>1&DEGs[,7]>1),]#在细胞系和TCGA二期样本中均显著差异且高表达logFC>1
# save(hDEGs,file='hDEGs.Rdata')
hTgenes=hDEGs[match(intersect(hDEGs[,1],Tgenes),hDEGs[,1]),]
symbols=ensg2gid[match(hTgenes[,1],ensg2gid$entrez_id),3]
hTgenes=cbind(symbols,hTgenes)
drug=c()
index=c()
for(i in 1:nrow(hTgenes)){
  print(i)
  drugid0=TargetGenes[which(TargetGenes$geneid==hTgenes$gid[i]),1]
  drug0=Drugid[match(drugid0,Drugid$drugID),1]
  drug=c(drug,paste(drug0, collapse = ','))
  index=c(index,length(drug0))
}
hTgenes=cbind(hTgenes,drug,index)
hTgenes$symbols=as.character(hTgenes$symbols)
hTgenes$drug=as.character(hTgenes$drug)
sep=read.table('CellLines.txt',sep='\t',header=T,as.is=T)
hTgenes=cbind(hTgenes,sep[match(hTgenes$gid,sep$Entrezid),])
# save(hTgenes,file='hTgenes.Rdata')
write.table(hTgenes,'hTgenes.txt',sep='\t',col.names = T,row.names = F)

##聚类分析结果#########
rm(list=ls())
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/data/')
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/data/emt_refergenes.Rdata")
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/data/GSE39582_inf.Rdata")
# load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/data/TCGA_inf.Rdata")
# load('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/cli_pre.Rdata')
myfun_pre <- function(GSE,cli,pairs,k=25){
  cli2=cli[grep('2',cli$stage),]
  gid=GSE[,1];exp=GSE[,match(cli2$sample,colnames(GSE))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  sam_l=names(count)[which(count>=k)]
  sam_h=names(count)[which(count<k)]
  sam_l=cli[match(sam_l,cli$sample),]
  sam_h=cli[match(sam_h,cli$sample),]
  cli_pre=list(sam_l=sam_l,sam_h=sam_h)
  return(cli_pre)
}
cli_pre=myfun_pre(GSE,cli,pairs_sur,k=25)
sam_l=cli_pre$sam_l;sam_h=cli_pre$sam_h
# save(sam_l,sam_h,file='cli_pre_GSE39582.Rdata')

exp=fpkm[match(emt2plat,fpkm$gid),]
log2exp=log2(exp)
log2exp[log2exp<0]=0
exp=log2exp

# exp=GSE[match(emt2plat,GSE$gid),]
exp_l=exp[,match(sam_l$sample,colnames(exp))]
exp_h=exp[,match(sam_h$sample,colnames(exp))]
exp=as.matrix(cbind(exp_l,exp_h))
library(gplots)
# hr=hclust(as.dist(1-cor(t(exp),method='pearson')),method='complete')
hc=hclust(as.dist(1-cor(exp,method='pearson')),method='complete')
colcol=c(rep('green',ncol(exp_l)),rep('red',ncol(exp_h)))
# pdf('emtHeatmap_TCGA.pdf')
heatmap.2(exp,#Rowv=as.dendrogram(hr),
          Colv=as.dendrogram(hc),
          scale='none',trace="none",col=greenred,
          breaks=seq(min(exp),max(exp),length=100),
          ColSideColors=colcol,labCol=F,main='TCGA_stage2')#labCol=colnames(exp)
# dev.off()
a=cutree(hc,2)
cliL=cli[match(names(which(a==1)),cli$sample),]
cliH=cli[match(names(which(a==2)),cli$sample),]
library(plyr)
sam11=transform(cliL,class='0')
sam22=transform(cliH,class='1')
cli_new=rbind(sam11,sam22)
cli_new=cli_new[!is.na(cli_new$RFS_event),]
cli_new=cli_new[!is.na(cli_new$RFS_time),]
label=cli_new$class
library(survival)
survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
fit=survfit(survi~label,data=cli_new)
cox=coxph(survi~label)
b=summary(cox)
hr=round(b$conf.int[1],2)
lower.hr=round(b$conf.int[3],2)
upper.hr=round(b$conf.int[4],2)
cindex=round(b$concordance[1],3)
p=signif(b$sctest[3],3)
# pdf('emt_GSE39582.pdf')
plot(fit,col=c('blue','red'),main='GSE39582',ylab='survival rate',
     xlab='survival months',lty=c(1,2),lwd=1.5,mark.time=T)
abline(v=60,lwd=1,lty=2,col='black')
legend('bottomright',pch=15:18,lty=c(1,2),merge=F,
       legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                paste('highrisk',length(which(label==1)),sep=' : ')),
       col=c('blue','red'))
text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n',
                      'C-index=',cindex,'\n',
                      'HR=',signif(hr,4),'\n',
                      '(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),
     bty='n',font=2)
# dev.off()
Lcli=intersect(sam_l$sample,cliL$sample)#115
Hcli=intersect(sam_h$sample,cliH$sample)#67

#kmeanCluster####
rm(list=ls())
exp=GSE[match(emt2plat,GSE$gid),]
exp_l=exp[,match(sam_l$sample,colnames(exp))]
exp_h=exp[,match(sam_h$sample,colnames(exp))]
exp=as.matrix(cbind(exp_l,exp_h))
ann_col=data.frame(gpscluster=c(rep('L',ncol(exp_l)),rep('H',ncol(exp_h))))
rownames(ann_col)=colnames(exp)
#利用k-mean是进行聚类
texp=t(exp)
set.seed(123)
k_result <- kmeans(texp, 2, nstart = 25)
table(k_result$cluster)
kcluster=data.frame(kcluster=k_result$cluster)
data=cbind(kcluster,ann_col)
data1=data[grep('1',data[,1]),]
data2=data[grep('2',data[,1]),]
data=rbind(data1,data2)
exp_new=exp[,match(rownames(data),colnames(exp))]
ann_colors=list(kcluster=c('1'='#FFC0CB','2'='#ADD8E6'),
                gpscluster=c('L'='green','H'='red'))
                # label=c(LL='#FF1493',LH='#1E90FF',HH='#FFD700',HL='#00FF00'))
# pdf('k_means_GSE39582.pdf')
library(pheatmap)
pheatmap(exp_new,color=colorRampPalette(rev(c('red','black','green')))(100),
         clustering_distance_rows='correlation',
         clustering_distance_cols='correlation',
         cluster_rows =T, cluster_cols=F,
         clustering_method='complete',
         annotation_col =data,
         annotation_colors=ann_colors,
         border_color="grey",legend = T,
         show_colnames=F,show_rownames=F,
         main='GSE39582')
# dev.off()
sam1=cli[match(rownames(kcluster)[grep('1',kcluster$kcluster)],cli$sample),]
sam2=cli[match(rownames(kcluster)[grep('2',kcluster$kcluster)],cli$sample),]
sam11=transform(sam1,class='0')
sam22=transform(sam2,class='1')
cli_new=rbind(sam11,sam22)
cli_new=cli_new[!is.na(cli_new$RFS_event),]
cli_new=cli_new[!is.na(cli_new$RFS_time),]
label=cli_new$class
library(survival)
survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
fit=survfit(survi~label,data=cli_new)
cox=coxph(survi~label)
b=summary(cox)
hr=round(b$conf.int[1],2)
lower.hr=round(b$conf.int[3],2)
upper.hr=round(b$conf.int[4],2)
cindex=round(b$concordance[1],3)
p=signif(b$sctest[3],3)
# pdf('emt_GSE39582.pdf')
plot(fit,col=c('blue','red'),main='GSE39582',ylab='survival rate',
     xlab='survival months',lty=c(1,2),lwd=1.5,mark.time=T)
abline(v=60,lwd=1,lty=2,col='black')
legend('bottomright',pch=15:18,lty=c(1,2),merge=F,
       legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                paste('highrisk',length(which(label==1)),sep=' : ')),
       col=c('blue','red'))
text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n',
                      'C-index=',cindex,'\n',
                      'HR=',signif(hr,4),'\n',
                      '(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),
     bty='n',font=2)
# dev.off()


#################################################
########ManuscriptPlot###########################
#################################################
rm(list=ls())
setwd('C:/Users/17362/Desktop/参考文献/Manuscript/Figures_Tables/ManuscriptPlot')
######
#1.KM_GPS####
k=25;pairs=pairs_sur
gid=GSE[,1];exp=GSE[,-1]
cli2=cli[grep('2',cli$stage),]
exp=exp[,match(cli2$sample,colnames(exp))]
diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
count=colSums(diff>0)
sam1=cli[match(names(count)[which(count>=k)],cli$sample),]
sam2=cli[match(names(count)[which(count<k)],cli$sample),]
library(plyr)
sam11=transform(sam1,class='0')
sam22=transform(sam2,class='1')
cli_new=rbind(sam11,sam22)
if(sur_status==1){
  cli_new=cli_new[!is.na(cli_new$DFS_event),]
  cli_new=cli_new[!is.na(cli_new$DFS_time),]
  label=cli_new$class
  library(survival)
  survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
}
if(sur_status==2){
  cli_new=cli_new[!is.na(cli_new$RFS_event),]
  cli_new=cli_new[!is.na(cli_new$RFS_time),]
  label=cli_new$class
  library(survival)
  survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
}
fit=survfit(survi~label,data=cli_new)
sur_time=summary(fit)
cox=coxph(survi~label,data=cli_new)
b=summary(cox)
hr=round(b$conf.int[1],2)
lower.hr=round(b$conf.int[3],2)
upper.hr=round(b$conf.int[4],2)
cindex=round(b$concordance[1],3)
p=signif(b$sctest[3],3)
pdf('GPS_GSE.pdf')
plot(fit,col=c('blue','red'),main='GSE39582',ylab='Disease-free survival',
     xlab='Disease-free survival (months)',lty=c(1,1),lwd=2,mark.time=T)
legend('bottomright',pch=15:18,lty=c(1,1),merge=F,
       legend=c(paste('Low-risk',length(which(label==0)),sep=' : '),
                paste('High-risk',length(which(label==1)),sep=' : ')),
       col=c('blue','red'))
text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n',
                      'HR=',signif(hr,4),'\n',
                      '(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),
     bty='n',font=2)
dev.off()

#2.C-index确定K值####
myfun_cindex_k <- function(GSE,cli,pairs){
  gid=GSE[,1];exp=GSE[,-1]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  count=colSums(diff>0)
  result=c()
  for(k in 0:nrow(pairs)){
    print(k)
    label=rep(0,ncol(diff))
    label[count<k]=1
    cli_new=cbind(cli2,label)
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event) #Surv(sur_time,sur_state)
    cox=coxph(survi~label+location+gender+age+MMR+CIMP+CIN+TP53+KRAS+BRAF,data=cli_new)#survi~label+location+gender+age+MMR+CIMP+CIN+TP53+KRAS+BRAF
    b=summary(cox)
    cindex=b$concordance[1]
    p=b$coefficients[1,5]
    result=rbind(result,c(k,cindex,p))
  }
  colnames(result)=c('k','cindex','p')
  return(result)
}
VoteK=myfun_cindex_k(GSE,cli,pairs_sur)#k=25(>=25-L;<25-H)   #k=20
write.table(VoteK,'VoteK570.txt',sep='\t',col.names = T,row.names = F)
# pdf('VoteK570.pdf')
plot(VoteK[,1],VoteK[,2],
     xlab='K',ylab='C-index',type='p',pch= 2 )
lines(VoteK[,1],VoteK[,2],col='red',lwd=3)
abline(v=25,lty=2)
# dev.off()

#3.ROC_curve####
myfun_div <- function(GSE,cli){
  gid=GSE[,1];exp=GSE[,-1]
  cli1=cli[grep('1',cli$stage),]
  exp1=exp[,match(cli1[,1],colnames(exp))]
  cli2=cli[grep('3|4',cli$stage),]
  exp2=exp[,match(cli2[,1],colnames(exp))]
  data=list(gid=gid,exp1=exp1,exp2=exp2)
  return(data)
}
myfun_roc <- function(GSE,cli,pairs){
  data=myfun_div(GSE,cli)
  gid=data$gid;exp1=data$exp1;exp2=data$exp2
  exp = cbind(exp1,exp2)
  sam1 = rep(0,ncol(exp1))#0--control
  sam2 = rep(1,ncol(exp2))#1--case
  label = c(sam1,sam2)
  diff = exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
  score = colMeans(diff>0)
  library(pROC)
  roc=roc(label,score,percent=F,ci=T,auc=T)
  # plot(roc,type="o",col="red",main='GSE39582',ylim=c(0:1))
  plot(roc,type="p", pch=20, main='Discovery Cohort',ylim=c(0:1))
  lines(roc,col='red',lwd=2)
  text(x=0.3,y=0.2,paste("AUC=",signif(roc$auc,4),seq=""),bty="n",font=2,cex=2)
}
#pdf('ROC.pdf')
myfun_roc(GSE,cli,pairs_sur)
#dev.off()
#4_1.对比IRGPI#####
myfun_km_k <- function(GSE,cli,GEO,sur_status){
  gid=GSE[,1];exp=GSE[,-1]
  cli=cli[grep('2',cli$stage),]
  GEOsamples=cluster[grep(GEO,cluster[,2]),]
  cli=cli[match(intersect(cli$sample,GEOsamples[,1]),cli$sample),]
  samcluster=GEOsamples[match(intersect(cli$sample,GEOsamples[,1]),GEOsamples[,1]),]
  sam1=samcluster[grep('Low Immune Risk',samcluster$immune.risk.group),1]
  sam2=samcluster[grep('High Immune Risk',samcluster$immune.risk.group),1]
  sam11=transform(cli[match(sam1,cli$sample),],class='0')
  sam22=transform(cli[match(sam2,cli$sample),],class='1')
  cli_new=rbind(sam11,sam22)
  if(sur_status==1){
    cli_new=cli_new[!is.na(cli_new$DFS_event),]
    cli_new=cli_new[!is.na(cli_new$DFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
  }
  if(sur_status==2){
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    cli_new=cli_new[!is.na(cli_new$RFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
  }
  fit=survfit(survi~label,data=cli_new)
  sur_time=summary(fit)
  cox=coxph(survi~label)
  b=summary(cox)
  hr=round(b$conf.int[1],2)
  lower.hr=round(b$conf.int[3],2)
  upper.hr=round(b$conf.int[4],2)
  cindex=round(b$concordance[1],3)
  p=signif(b$sctest[3],3)
  # pdf('IRGPI_GSE.pdf')
  plot(fit,col=c('blue','red'),main='GSE',ylab='survival rate',
       xlab='survival months',lty=c(1,2),lwd=1.5,mark.time=T)
  abline(v=60,lwd=1,lty=2,col='black')
  legend('bottomright',pch=15:18,lty=c(1,2),merge=F,
         legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                  paste('highrisk',length(which(label==1)),sep=' : ')),
         col=c('blue','red'))
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n',
                        'C-index=',cindex,'\n',
                        'HR=',signif(hr,4),'\n',
                        '(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),
       bty='n',font=2)
  # dev.off()
} #sur_status:DFS-1,RFS-2
myfun_km_k(GSE,cli,GEO='GSE',sur_status = 1)#ACT=N

#4_2.对比oncotypeDx colon####
load("D:/data/TCGA/ensg2gid.Rdata")
myfun_Dx <- function(x){
  gene1=c('BGN','FAP','INHBA')
  gene2=c('MKI67','MYC','MYBL2')
  gene3='GADD45B'
  refergene=c('ATP5E','GPX1','PGK1','VDAC2','UBB')
  gene1=ensg2gid[match(gene1,ensg2gid$symbol),2]
  gene2=ensg2gid[match(gene2,ensg2gid$symbol),2]
  gene3=ensg2gid[match(gene3,ensg2gid$symbol),2]
  refergene=ensg2gid[match(refergene,ensg2gid$symbol),2]
  refer=mean(x[match(refergene,GSE[,1])])
  RSu=0.15*mean(x[match(gene1,GSE[,1])])-0.3*mean(x[match(gene2,GSE[,1])])+0.15*(x[match(gene3,GSE[,1])])
  RS=44*(RSu+0.82)
  return(RS)
}
# RS=apply(exp,2,myfun_Dx)
# RSu = 0·15* mean(BGN,FAP,INHBA)-0·3*mean(MKI67,MYC,MYBL2)+ 0·15* GADD45B)
# This score was then rescaled RS=44*(RSu+0.82)
myfun_km_Dx <- function(GSE,cli,sur_status){
  gid=GSE[,1];exp=GSE[,-1]
  cli2=cli[grep('2',cli$stage),]
  exp=exp[,match(cli2$sample,colnames(exp))]
  RS=apply(exp,2,myfun_Dx)
  sam1=cli[match(names(RS)[which(RS<30)],cli$sample),]
  sam2=cli[match(names(RS)[which(RS>=41)],cli$sample),]
  library(plyr)
  sam11=transform(sam1,class='0')
  sam22=transform(sam2,class='1')
  cli_new=rbind(sam11,sam22)
  if(sur_status==1){
    cli_new=cli_new[!is.na(cli_new$DFS_event),]
    cli_new=cli_new[!is.na(cli_new$DFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$DFS_time,cli_new$DFS_event)
  }
  if(sur_status==2){
    cli_new=cli_new[!is.na(cli_new$RFS_event),]
    cli_new=cli_new[!is.na(cli_new$RFS_time),]
    label=cli_new$class
    library(survival)
    survi=Surv(cli_new$RFS_time,cli_new$RFS_event)
  }
  fit=survfit(survi~label,data=cli_new)
  sur_time=summary(fit)
  cox=coxph(survi~label)
  b=summary(cox)
  hr=round(b$conf.int[1],2)
  lower.hr=round(b$conf.int[3],2)
  upper.hr=round(b$conf.int[4],2)
  cindex=round(b$concordance[1],3)
  p=signif(b$sctest[3],3)
  # pdf('Dx_TCGA.pdf')
  plot(fit,col=c('blue','red'),main='GSE',ylab='survival rate',
       xlab='survival months',lty=c(1,2),lwd=1.5,mark.time=T)
  abline(v=60,lwd=1,lty=2,col='black')
  legend('bottomright',pch=15:18,lty=c(1,2),merge=F,
         legend=c(paste('lowrisk',length(which(label==0)),sep=' : '),
                  paste('highrisk',length(which(label==1)),sep=' : ')),
         col=c('blue','red'))
  text(x=40,y=0.1,paste('log-rank P=',signif(p,digits=3),'\n',
                        'C-index=',cindex,'\n',
                        'HR=',signif(hr,4),'\n',
                        '(','95%CI,',signif(lower.hr,4),'-',signif(upper.hr,4),')'),
       bty='n',font=2)
  # dev.off()
} #sur_status:DFS-1,RFS-2
myfun_km_Dx(fpkm,cli,sur_status = 1)#ACT=N

#4_3.cindex条形图对比结果####
rm(list=ls())
setwd('C:/Users/17362/Desktop/参考文献/Manuscript/Figures_Tables/ManuscriptPlot')
x=matrix(c(1.5,2.5,3.5,5.5,6.5,7.5,9.5,10.5,11.5,13.5,14.5,15.5,17.5,18.5,19.5),
         ncol=5,nrow=3)
data=read.table('cindex.txt',sep='\t',header = T,fill=T,as.is=T)
cindex=as.matrix(data[1:3,-1])
se=as.matrix(data[5:7,-1])
pdf('C-index for Three Signatures.pdf')
barplot(cindex, col=colors()[c(23,89,12)], border='white', font.axis=2, beside=T, 
        legend=data[1:3,1], ylab='C-index', font.lab=2 , ylim = 0:1, xpd=F)
plot.error <- function(x, y, sd, len = 0.05, col = "black", horiz = FALSE) {
  if (!horiz) {
    arrows(x0 = x, y0 = y, x1 = x, y1 = y - sd/2, col = col, angle = 90, length = len)
    arrows(x0 = x, y0 = y, x1 = x, y1 = y + sd/2, col = col, angle = 90, length = len)
  } else {
    arrows(x0 = y, y0 = x, x1 = y - sd/2, y1 = x, col = col, angle = 90, length = len)
    arrows(x0 = y, y0 = x, x1 = y + sd/2, y1 = x, col = col, angle = 90, length = len)
  }
}
for (i in 1:3) plot.error(x[i, ], cindex[i, ], sd = se[i, ])
dev.off()

#5.保存cli信息&MultiCox####
rm(list=ls())
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/data/cli/')
# write.table(cli,'cli.txt',sep='\t',col.names = T,row.names = F)
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/pairs_multicox.Rdata")
cli=read.table('cli17538.txt',sep='\t',header = T,fill=T,as.is = T)
gid=GSE[,1];exp=GSE[,-1];k=25;pairs=pairs_sur
cli2=cli[grep('2',cli$stage),]
exp=exp[,match(cli2$sample,colnames(exp))]
diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
count=colSums(diff>0)
sam1=cli[match(names(count)[which(count>=k)],cli$sample),]
sam2=cli[match(names(count)[which(count<k)],cli$sample),]
library(plyr)
sam11=transform(sam1,class='0')
sam22=transform(sam2,class='1')
cli_new=rbind(sam11,sam22)
cli_new=cli_new[!is.na(cli_new$DFS_event),]
cli_new=cli_new[!is.na(cli_new$DFS_time),]
label=as.numeric(cli_new$class)
library(survival)
survi=Surv(cli_new$DFS_time,cli_new$DFS_event);colnames(cli_new)
cox=coxph(survi~label+gender+age+grade,data = cli_new)
summary(cox)

b=summary(cox)
hr=round(b$conf.int[1],2)
lower.hr=round(b$conf.int[3],2)
upper.hr=round(b$conf.int[4],2)
cindex=round(b$concordance[1],3)
p=signif(b$sctest[3],3)

#6.多组学分析热图####
#转变为突变离散谱数据
rm(list=ls())
load("D:/data/TCGA/TCGA_mutation.Rdata")
mdata[grep('3UTR',mdata$Variant_Classification),2]='3UTR'
mdata[grep("3'UTR",mdata$Variant_Classification),2]='3UTR'
mdata[grep('5UTR',mdata$Variant_Classification),2]='5UTR'
mdata[grep("5'UTR",mdata$Variant_Classification),2]='5UTR'
gene=unique(mdata$Entrez_Gene_Id)
sample=unique(mdata$Tumor_Sample_Barcode)
data=matrix(0,nrow=length(gene),ncol=length(sample))
sep=paste(mdata$Entrez_Gene_Id,mdata$Tumor_Sample_Barcode,sep=',')
for(i in 1:length(sample)){
  print(i)
  for(j in 1:length(gene)){
    index<-unique(mdata[which(paste(gene[j],sample[i],sep=',')==sep),2])
    if(length(index)==1){data[j,i]<-index}
    else if(length(index)>1){
      if("Missense_Mutation"%in%index){data[j,i]<-"Missense_Mutation"}
      else if("Frame_Shift_Del"%in%index){data[j,i]<-"Frame_Shift_Del"}
      else if("Intron"%in%index){data[j,i]<-"Intron"}
      else if("3UTR"%in%index){data[j,i]<-"3UTR"}
      else if("Nonsense_Mutation"%in%index){data[j,i]<-"Nonsense_Mutation"}
      else if("Frame_Shift_Ins"%in%index){data[j,i]<-"Frame_Shift_Ins"}
      else if("Splice_Region"%in%index){data[j,i]<-"Splice_Region"}
      else if("Splice_Site"%in%index){data[j,i]<-"Splice_Site"}
      else if("5UTR"%in%index){data[j,i]<-"5UTR"}
      else if("In_Frame_Del"%in%index){data[j,i]<-"In_Frame_Del"}
      else if("In_Frame_Ins"%in%index){data[j,i]<-"In_Frame_Ins"}
      else {data[j,i]<-"other"}
    }
  }
}
mut=cbind(gene,data)
colnames(mut)=c('gene',sample)
save(mut,file='mut.Rdata')
#6.1_ComplexHeatmap 绘制肿瘤突变分布图####
# source("http://bioconductor.org/biocLite.R")
# BiocManager::install('ComplexHeatmap')
# install.packages('circlize')
library(circlize)
library(grid)
library(ComplexHeatmap)
load("C:/Users/17362/Desktop/参考文献/Manuscript/Figures_Tables/ManuscriptPlot/mut.Rdata")
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/multi-omics/mutataion/mutation_fisher.Rdata")
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/cli_pre_TCGA.Rdata")
load("D:/data/TCGA/ensg2gid.Rdata")
lsam=intersect(sam_l$sample,colnames(mut));hsam=intersect(sam_h$sample,colnames(mut))
mdata=mut[match(SigMgenes$gid,mut[,1]),c(1,match(c(lsam,hsam),colnames(mut)))]#lrisk在前,hrisk在后
mdata1=mdata[order(rowSums(mdata!=0),decreasing = T),]
rownames(mdata1)=ensg2gid[match(mdata1[,1],ensg2gid$entrez_id),3]
mdata=mdata1[1:50,-1]
mdata[which(mdata=='3UTR',arr.ind = T)]='UTR3'
mdata[which(mdata=='5UTR',arr.ind = T)]='UTR5'
save(mdata,file='mdata.Rdata')
mdata[which(mdata==0,arr.ind = T)]=""
col = c("Missense_Mutation" = "red","Intron" = "violet","Frame_Shift_Del"="steelblue1","Nonsense_Mutation"="palegreen","Frame_Shift_Ins"="yellow3","UTR3"="darkorchid","Splice_Region"="tan1","UTR5"="palevioletred1","Splice_Site"="tan4","In_Frame_Del"="lightblue","In_Frame_Ins"="wheat")
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "gray85", col = NA))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  Intron = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Intron"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  UTR3 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["UTR3"], col = NA))
  },
  UTR5 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["UTR5"], col = NA))
  },
  Splice_Region = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Splice_Region"], col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  }
)
# pdf("ComplexHeatmapTop30.pdf")
muttype=c("Missense_Mutation","Intron","Frame_Shift_Del","Nonsense_Mutation","Frame_Shift_Ins","UTR3","Splice_Region","UTR5","Splice_Site","In_Frame_Del","In_Frame_Ins")
heatmap_legend_param = list(cluster_columns = TRUE,title = "Alternations", at = muttype,labels = muttype)
oncoPrint(mdata ,alter_fun = alter_fun, col = col, 
          column_order = colnames(data),
          remove_empty_columns = TRUE,
          show_pct = TRUE, pct_gp = gpar(fontsize = 8),row_names_gp = gpar(fontsize = 8),
          column_title = "Altofrequency mutation genes",
          heatmap_legend_param = heatmap_legend_param)
# dev.off()
#确定高低风险样本分界线
count=colSums(mdata!="")
index=count[-which(count==0)]
which(names(index)=='TCGA_A6_2675')#第43个为第一个高风险样本

#6.2_pheatmap 绘制CNVs分布图####
rm(list=ls())
setwd('C:/Users/17362/Desktop/参考文献/Manuscript/Figures_Tables/ManuscriptPlot/')
load("D:/data/TCGA/TCGA_cnv.Rdata")
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/multi-omics/cnv/CNV_fisher.Rdata")
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/cli_pre_TCGA.Rdata")
region=c(SigCNV$segment)
sam=c(sam_l$sample,sam_h$sample)
myfun_strsplit <- function(x){
  sep=unlist(strsplit(x,split = ' '))
  return(sep[1])
}
cnvdata=cnv[match(SigCNV$segment,sapply(cnv[, 1],FUN=myfun_strsplit)),
            match(intersect(sam,colnames(cnv)),colnames(cnv))]
library(circlize)
library(grid)
library(ComplexHeatmap)
mdata=apply(cnvdata,2,as.numeric)
rownames(mdata)=SigCNV$segment
rowSums(mdata)
mdata[which(mdata==0,arr.ind = T)]=""
col=c('2'='red','1'='violet')
# col = c("Missense_Mutation" = "red","Intron" = "violet","Frame_Shift_Del"="steelblue1","Nonsense_Mutation"="palegreen","Frame_Shift_Ins"="yellow3","UTR3"="darkorchid","Splice_Region"="tan1","UTR5"="palevioletred1","Splice_Site"="tan4","In_Frame_Del"="lightblue","In_Frame_Ins"="wheat")
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "gray85", col = NA))
  },
  '1' = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), 
              gp = gpar(fill = col["1"], col = NA))
  },
  '2' = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"),  h-unit(0.3, "mm"), 
              gp = gpar(fill = col["2"], col = NA))
  }
)
muttype=c("2","1")
heatmap_legend_param = list(cluster_columns = TRUE,title = "CNV level" , at = muttype,labels = muttype)

pdf("ComplexHeatmap_CNVs.pdf",width = 7, height = 1.8)
oncoPrint(mdata ,
          alter_fun = alter_fun, col = col, 
          column_order = colnames(mdata),
          remove_empty_columns = TRUE,
          show_pct = TRUE, pct_gp = gpar(fontsize = 7),
          row_names_gp = gpar(fontsize = 7),
          gap = unit(5, "mm"),
          heatmap_legend_param = heatmap_legend_param)
dev.off()
count=colSums(mdata!="")
index=count[-which(count==0)]
which(names(index)=='TCGA_5M_AATE')#第65个为第一个高风险样本
#7.SFRP4基因对表达变化####
rm(list=ls())
setwd('C:/Users/17362/Desktop/参考文献/Manuscript/Figures_Tables/ManuscriptPlot/')
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/pairs_multicox.Rdata")
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/data/TCGA_inf.Rdata")
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/cli_pre_TCGA.Rdata")
load("D:/data/TCGA/ensg2gid.Rdata")
sep=pairs_sur[which(pairs_sur=='6424',arr.ind = T)[1],1:2]#SFRP4基因对
ensg2gid[match(sep,ensg2gid$entrez_id),]
data=fpkm[match(sep,fpkm$gid),]
expl=data[,match(sam_l$sample,colnames(data))]
expl1=expl[,order(expl[1,],decreasing = F)]
exph=data[,match(sam_h$sample,colnames(data))]
exph1=exph[,order(exph[1,],decreasing = T)]
exp=cbind(expl1,exph1)
pdf('GPSEXP1.pdf')
plot(x=1:ncol(exp),y=exp[1,],pch=20,
     xlab='Samples',ylab='Expression level',
     col='blue',
     main='SFRP4')
points(x=1:ncol(expl),y=expl[2,],pch=14,col='seagreen')
points(x=(1:ncol(exph))+ncol(expl),y=exph[2,],pch=17,col='red')
# plot(x=1:ncol(exp),y=exp[1,],pch=20,
#      xlab='Samples',ylab='Expression level',
#      col='blue')
# for(i in 1:ncol(expl1)){
#   print(i)
#   if(expl1[1,i]>expl1[2,i]){
#     points(x=i,y=expl1[2,i],pch=17,col='seagreen')
#   }
#   else {points(x=i,y=expl1[2,i],pch=2)}
# }
# for(i in 1:ncol(exph1)){
#   print(i)
#   j=i+ncol(expl1)
#   if(exph1[1,i]<exph1[2,i]){
#     points(x=j,y=exph1[2,i],pch=17,col='red')
#   }
#   else {points(x=j,y=exph1[2,i],pch=2)}
# }
abline(v=(ncol(expl)+0.5),lwd=1)
legend('bottomright',pch=c(19,17),
       legend=c('ATP23','SFRP4'),
       col=c('blue','gray28'))
dev.off()

#7.共表达分析####
rm(list=ls())
setwd('E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/validation/Coexpression')
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/data/TCGA_inf.Rdata")
load("D:/data/TCGA/ensg2gid.Rdata")
gid=fpkm[,1];exp2=fpkm[,match(cli[grep('2',cli$stage),1],colnames(fpkm))]
SFRP4=ensg2gid[match('SFRP4',ensg2gid$symbol),2]#SFRP4的entrezid:6424
SFRP4exp=exp2[match('6424',fpkm$gid),,drop=F]
gid=gid[-which(gid=='6424')];exp=exp2[match(gid,fpkm[,1]),]
result=c()
for(i in 1:length(gid)){
  print(i)
  cortest=cor.test(as.numeric(SFRP4exp[1,]),as.numeric(exp[i,]),method='pearson')
  result=rbind(result,c(gid[i],cortest$estimate,cortest$p.value))
}
result=cbind(result,p.adjust(result[,3],method = 'BH'))
result=cbind(ensg2gid[match(result[,1],ensg2gid$entrez_id),3],result)
colnames(result)=c('symbol','gid','cor','p','fdr')
SigCorgene=result[which(result[,5]<0.05),]

#富集分析####
# source("http://www.bioconductor.org/biocLite.R")
# BiocManager::install("AnnotationHub") #下载安装数据包，缺少的数据包都可以这样安装
# BiocManager::install("org.Mmu.eg.db")
# BiocManager::install("clusterProfiler")
# BiocManager::install("topGO")
# BiocManager::install("Rgraphviz")
rm(list=ls())
library(AnnotationHub)  #library导入需要使用的数据包
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
id=SigCorgene[,2]
#KEGG分析#
ekk <- enrichKEGG(id,organism='hsa',pvalueCutoff=0.05)#KEGG富集分析
kegg=as.data.frame(ekk)
#pdf('Bubble plot_Enrichment KEGG.pdf')
dotplot(ekk,font.size=12,showCategory=10,title='Bubble plot_Enrichment KEGG-18_22')#画气泡图
#dev.off()
write.table(kegg,"KEGG-enrich.txt",sep='\t',row.names =F)
# browseKEGG(ekk,'mmu04961')  #显示通路图

#GO富集分析
ego <- enrichGO(OrgDb="org.Hs.eg.db",gene=id,keyType = "ENTREZID",
                ont="BP",pvalueCutoff=0.05,readable=F) #GO富集分析
dotplot(ego,showCategory=10,title="Enrichment GO Top10") #泡泡图
barplot(ego, showCategory=20,title="EnrichmentGO")  #柱状图
plotGOgraph(ego)  #GO图
gobp=as.data.frame(ego)
write.table(gobp,"GObp-enrich.txt",sep='\t',row.names =F)

sep=fpkm[match(id,fpkm$gid),]
# gid=sep[,1];data=sep[,-1]
explh=sep[,c(1,match(c(sam_l$sample,sam_h$sample),colnames(sep)))]
explh=cbind(ensg2gid[match(explh[,1],ensg2gid$entrez_id),3],explh)
write.table(explh,'explh.txt',sep='\t',col.names = T,row.names = F,quot=F)

#survAUC####
# install.packages('survivalROC')
rm(list=ls())
library(survivalROC)
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/data/GSE39582_inf.Rdata")
load("E:/colorectal_sur/10_emt2minusgenes/GSE_TCGA/pairs_multicox.Rdata")
pairs=pairs_sur;k=25;cutoff <- 12*5
gid=GSE[,1];exp=GSE[,-1]
cli2=cli[grep('2',cli$stage),]
exp=exp[,match(cli2$sample,colnames(exp))]
diff=exp[match(pairs[,1],gid),]-exp[match(pairs[,2],gid),]
count=colSums(diff>0)
sam1=cli[match(names(count)[which(count>=k)],cli$sample),]
sam2=cli[match(names(count)[which(count<k)],cli$sample),]
library(plyr)
sam11=transform(sam1,class='0')
sam22=transform(sam2,class='1')
cli_new=rbind(sam11,sam22)
cli_new=cli_new[!is.na(cli_new$RFS_event),]
nobs <- nrow(cli_new)
Mayo4.1= survivalROC(Stime=cli_new$RFS_time,##生存时间
                     status=cli_new$RFS_event,## 终止事件    
                     marker = as.numeric(cli_new$class), ## marker value    
                     predict.time = cutoff,## 预测时间截点
                     span = 0.25*nobs^(-0.20))##span,NNE法的namda
str(Mayo4.1)## list结构
## 绘图
plot(Mayo4.1$FP, Mayo4.1$TP, ## x=FP,y=TP
     type="l",col="red", ##线条设置
     xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.1$AUC,3)), ##连接
     ylab="TP",
     main="Mayoscore 4, Method = NNE \n  Year = 1")## \n换行符
abline(0,1,col="gray",lty=2)##线条颜色

