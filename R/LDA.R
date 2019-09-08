#setClass("tfLDA", slots=list(documents="list", wordnames="character", samplenames="character",models = "list",params="list"))

tfLDA<-methods::setClass(
  "tfLDA",
  slots=list(documents="list", wordnames="character", samplenames="character",models = "list",params="list")
)


createtfLDAObj <- function(documents,wordnames,samplenames){
  obj=new(Class = "tfLDA",documents=documents,
          wordnames=wordnames,samplenames=samplenames,models=list(),params=list())
  obj
}


addModels <- function(modelList){
  names(modelList) <- laply(1:length(modelList), function(x) sapply(modelList[x], function(y) nrow(y$topic_sums)))
  modelList <- modelList[as.character(sort(as.numeric(names(modelList))))]
  modelList <- modelList[unique(names(modelList))]
  modelList
}



tfLDA <- function(obj,topic=c(5, 10, 20),seed=1234,iterations = 500,burnin = NULL,alpha = 50,eta=0.1,alphaByTopic=TRUE) {
    
  print('Running tfLDA...')
  
  set.seed(seed)

  if (alphaByTopic==TRUE){
      models <- suppressWarnings(llply(.data=topic, .fun=function(t) 
          lda.collapsed.gibbs.sampler(obj@documents, t, obj@wordnames, num.iterations=iterations, alpha=alpha/t, eta, compute.log.likelihood = TRUE, burnin=burnin), .progress = progress_text(char = ".")))
  }else{
      models <- suppressWarnings(llply(.data=topic, .fun=function(t)
          lda.collapsed.gibbs.sampler(obj@documents, t, obj@wordnames, num.iterations=iterations, alpha=alpha, eta, compute.log.likelihood = TRUE, burnin=burnin), .progress = progress_text(char = ".")))
  }
  
  if (length(obj@models)>0){
      models <- addModels(c(obj@models, models))
  } else {
    names(models) <- laply(1:length(models), function(x) sapply(models[x], function(y) nrow(y$topic_sums)))
  }
  
  obj@models <- models
  obj
  
}



plotLogLikelihood <- function(obj){
  
  loglikelihoods <- sapply(seq_along(obj@models), FUN=function(i) obj@models[[i]]$log.likelihood[2,])
  iterations <- ncol(obj@models[[1]]$log.likelihood)
  topics <-  sapply(seq_along(obj@models), FUN=function(i) nrow(obj@models[[i]]$topics))
  colnames(loglikelihoods) <- paste(topics, 'topics')
  
  par(bty = 'n')
  matplot(1:iterations,loglikelihoods,  type = 'l', lty=1, lwd=4, col = 1:length(topics), xlab="Iteration", ylab="Log likelihood", main='')
  legend("bottomright", legend =topics, fill=1:length(topics))

}



hdpConvert<-function(documents,outfile){
  for(i in 1:length(documents)){
    b=table(documents[[i]][1,])
    a=as.numeric(names(b))
    N=length(a)
    d=paste(a,b,sep=":")
    d=paste(c(N,d),"",collapse=" ")
    write.table(d,file=outfile,quote=F,col.names=F,row.names=F,append=T)
  }
}



createDocuments<-function(motifdir,database=c('JASPAR','TRANSFAC')){
  
  if(database=='JASPAR'){
    data(JASPAR)
    #load('../JASPAR.rda')
    allmotifs=JASPAR
  }else if(database=='TRANSFAC'){
    data(TRANSFAC)
    #load('../TRANSFAC.rda')
    allmotifs=TRANSFAC
  }
  allmotif.names=names(allmotifs)
  files=dir(motifdir,full.name=TRUE)
  documents=NULL;index=1
  motifcounts=list()
  for(i in 1:length(files)){
    mat=read.table(files[i])
    motifcounts[[i]]=mat[,c(1,2,4,8)]
    colnames(motifcounts[[i]])=c("motifname","nhits","nhits.bg","nword")
    matchid=match(allmotif.names,mat[,1])
    matchcount=mat[matchid,][,8]
    matchlist=sapply(c(1:length(allmotif.names)),function(x) rep((x-1),matchcount[x]))
    mat1=matrix(c(unlist(matchlist),rep(1,sum(matchcount))),byrow=TRUE,nrow=2)
    storage.mode(mat1)="integer"
    documents[[index]]=mat1
    index=index+1
  }
  
  documents2=list()
  for(i in 1:length(documents)){
    b=table(documents[[i]][1,])
    a=as.numeric(names(b))
    m=rbind(as.integer(a),as.integer(b))
    rownames(m)=colnames(m)=NULL
    documents2[[i]]=m
  }
  
  names(documents)=names(documents2)=names(motifcounts)=basename(files)
  #unlist(strsplit(basename(files),split='\\.')[[1]][1])
  list(documents=documents,documents2=documents2,motifcounts=motifcounts)

}







topMotif<-function(obj,ntopic,method=c('gamma','rank'),rank=10,thrP=0.96){
  
  
  if(length(obj@models)==0){
    stop('tfLDA should run first...\n')
  }
  
  ntopic=as.character(ntopic)
  topics=obj@models[[ntopic]]$topics
  document_sums=obj@models[[ntopic]]$document_sums
  topic_sums=obj@models[[ntopic]]$topic_sums
  assignments=obj@models[[ntopic]]$assignments
  log.likelihoods=obj@models[[ntopic]]$log.likelihoods
  samplenames=obj@samplenames
  wordnames=obj@wordnames
  nword=length(wordnames)
  topic.word=top.topic.words(topics, num.words = nword)
  ntopic=as.numeric(ntopic)
  num.topic.word=prob.topic.word=matrix(double(ntopic*nword),ntopic,nword)
  for(i in 1:ntopic){
    num.topic.word[i,]=topics[i,][order(topics[i,],decreasing=T)]
    prob.topic.word[i,]=num.topic.word[i,]/topic_sums[i]
  }
  prob.topic.word=t(prob.topic.word)
  
  prob.topic.word.ls=topic.word.ls=list()
  for (i in 1:ncol(prob.topic.word)){
    if(method=='gamma'){
      distr=suppressWarnings(fitdist(prob.topic.word[,i],"gamma", method="mme"))
      cutoff=as.numeric(unlist(quantile(distr, probs = thrP))[1])
      topic.word.ls[[i]]=topic.word[which(prob.topic.word[,i] > cutoff), i, drop = FALSE]
      prob.topic.word.ls[[i]]=prob.topic.word[,i][which(prob.topic.word[,i] > cutoff)]
    }else if(method=='rank'){
      topic.word.ls[[i]]=topic.word[1:rank,i]
      prob.topic.word.ls[[i]]=prob.topic.word[1:rank,i]
    }
  }
  names(topic.word.ls)=names(prob.topic.word.ls)=paste("Topic",1:ncol(prob.topic.word),sep="")
  obj@params$topic.word=topic.word.ls
  obj@params$prob.topic.word=prob.topic.word.ls
  obj
  
}



