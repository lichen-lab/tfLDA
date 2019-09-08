**Reproduce all the figures for HeLaS3 in the manuscript**
----------------------
```
library(lda)
library(pheatmap)
library(RColorBrewer)
library(tfLDA)
data("HeLaS3")
obj=createtfLDAObj(documents=HeLaS3@documents,wordnames=HeLaS3@wordnames,samplenames=HeLaS3@samplenames)


#run LDA model
ntopic=12
obj=tfLDA(obj,topic=ntopic,iterations=1000,alpha=50)
plotLogLikelihood(obj)

topics=obj@models[[1]]$topics
document_sums=obj@models[[1]]$document_sums
topic_sums=obj@models[[1]]$topic_sums
assignments=obj@models[[1]]$assignments
log.likelihoods=obj@models[[1]]$log.likelihoods
samplenames=obj@samplenames
wordnames=obj@wordnames



# Figure1: the whole heatmap with all samples and TFs
theta=matrix(0,length(assignments),ntopic)
for(i in 1:length(assignments)){
  a=table(assignments[[i]])
  b=match(as.numeric(names(a)),0:(ntopic-1))
  theta[i,b]=a
  theta[i,]=theta[i,]/sum(theta[i,])
}
theta_mean=aggregate(x=theta,by=list(obj@samplenames), FUN="mean")
rownames(theta_mean)=theta_mean[,1]
theta_mean=theta_mean[,-1]
colnames(theta_mean)=paste("Topic",1:ntopic,sep="")
tfs=rownames(theta_mean)

b1=c("ELK4","ELK1")
b2=c("AP2gamma","AP2alpha")
b3=c("BAF170","BAF155")
b4=c("Max","cMyc")
b5=c("NFYA","NFYB")
b6=c("Rad21","SMC3")

p1=match(b1,tfs)
p2=match(b2,tfs)
p3=match(b3,tfs)
p4=match(b4,tfs)
p5=match(b5,tfs)
p6=match(b6,tfs)

df=rep(0,nrow(theta_mean))
df[p1]=1
df[p2]=2
df[p3]=3
df[p4]=4
df[p5]=5
df[p6]=6

labels=rep("others",nrow(theta_mean))
labels[p1]="block1"
labels[p2]="block2"
labels[p3]="block3"
labels[p4]="block4"
labels[p5]="block5"
labels[p6]="block6"

blist=plist=list()
blist[[1]]=b1
blist[[2]]=b2
blist[[3]]=b3
blist[[4]]=b4
blist[[5]]=b5
blist[[6]]=b6

plist[[1]]=p1
plist[[2]]=p2
plist[[3]]=p3
plist[[4]]=p4
plist[[5]]=p5
plist[[6]]=p6

annotation_row=data.frame(TFs = labels)
rownames(annotation_row)=rownames(theta_mean)

pdf(file='Helas3-heatmap.pdf')
pheatmap(theta_mean,fontsize_row =7,fontsize=8,annotation_row=annotation_row,main="Helas3",treeheight_row=0,treeheight_col=0,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100))
dev.off()



# Figure2: subset of whole heatmap
tfset=c("STAT3","STAT1","Rad21","SMC3","NFYA","NFYB","cMyc","Max","ELK1","ELK4")
topicset=paste("Topic",c(6,9,8,4,11,1),sep="")
rowid=match(tfset,rownames(theta_mean))
colid=match(topicset,colnames(theta_mean))
theta_mean2=theta_mean[rowid,colid]
dim(theta_mean2)

labels=rep("",length(tfset))
labels[1:2]="block1"
labels[3:4]="block2"
labels[5:6]="block3"
labels[7:8]="block4"
labels[9:10]="block5"

annotation_row=data.frame(TFs = labels)
rownames(annotation_row)=rownames(theta_mean2)

pdf(file='Helas3-subset-heatmap.pdf')
pheatmap(theta_mean2,cluster_rows=F,cluster_cols=F,annotation_row=annotation_row,main="Helas3",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100))
dev.off()


  

# Figure3: Helas3-topicmotif.pdf
mypar=function (a = 2, b = 2, brewer.n = 8, brewer.name = "Dark2", ...) {
  par(mar = c(2.5, 5, 1.6, 1.1), mgp = c(1.5, 0.5, 0))
  par(mfrow = c(a, b), ...)
}


obj=topMotif(obj,ntopic=ntopic,method='rank')
topic.word=obj@params$topic.word
prob.topic.word=obj@params$prob.topic.word
topicorder=order(obj@models[[1]]$topic_sums,decreasing = T)
topic.word=topic.word[topicorder]
prob.topic.word=prob.topic.word[topicorder]
cols=rainbow(ntopic)
cols=cols[topicorder]

pdf(file="Helas3-topicmotif.pdf",width=10,height=5)
mypar(3,5)
for(i in 1:ntopic){
  rbPal=colorRampPalette(c("white",cols[i]))
  barplot(rev(prob.topic.word[[i]]),col=rbPal(length(topic.word[[i]])),horiz=T,names.arg=rev(topic.word[[i]]),las=1,cex.names=0.8,main=paste("Topic",topicorder[i],sep=""),
          col.main=cols[i],col.lab=cols[i],font.axis =2,axes=T,xlim=c(0,0.3))
}
dev.off()



obj=topMotif(obj,ntopic=ntopic,method='gamma')
topic.word=obj@params$topic.word
prob.topic.word=obj@params$prob.topic.word
topicorder=order(obj@models[[1]]$topic_sums,decreasing = T)
topic.word=topic.word[topicorder]
prob.topic.word=prob.topic.word[topicorder]
cols=rainbow(ntopic)
cols=cols[topicorder]

pdf(file="Helas3-topicmotif2.pdf",width=10,height=5)
mypar(3,5)
for(i in 1:ntopic){
  rbPal=colorRampPalette(c("white",cols[i]))
  barplot(rev(prob.topic.word[[i]]),col=rbPal(length(topic.word[[i]])),horiz=T,names.arg=rev(topic.word[[i]]),las=1,cex.names=0.8,main=paste("Topic",topicorder[i],sep=""),
          col.main=cols[i],col.lab=cols[i],font.axis =2,axes=T,xlim=c(0,0.3))
}
dev.off()




```

