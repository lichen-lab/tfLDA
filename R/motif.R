
fmotifscan<-function(seqIP,seqCT,outfile,database=c('JASPAR','TRANSFAC'),score="80%",ncore=1){
  
  registerDoParallel(cores=ncore)
  if(database=='JASPAR'){
    #data(JASPAR)
    load('../JASPAR.rda')
    allmotifs=JASPAR
  }else if(database=='TRANSFAC'){
    #data(TRANSFAC)
    load('../TRANSFAC.rda')
    allmotifs=TRANSFAC
  }
  allmotifs.name=names(allmotifs)
  motif.num=length(allmotifs)
 
  
  seqIP=readDNAStringSet(seqIP, format="fasta")
  seqIPlist=as.list(seqIP)
  seqCT=readDNAStringSet(seqCT, format="fasta")
  seqCTlist=as.list(seqCT)
  nsample=length(seqIPlist)
  
  Z=foreach(i=1:length(allmotifs),.combine='rbind') %dopar%{
    motif.len=ncol(allmotifs[[i]])
    y<-foreach(j=1:nsample,.combine='cbind') %dopar%{
      seqIP=toString(seqIPlist[[j]])
      matchres.IP=suppressWarnings(matchPWM(allmotifs[[i]], DNAString(seqIP), min.score=score))
      seqCT=toString(seqCTlist[[j]])
      matchres.CT=suppressWarnings(matchPWM(allmotifs[[i]], DNAString(seqCT), min.score=score))
      matchnum.IP=nrow(as.matrix(ranges(matchres.IP)))
      matchnum.CT=nrow(as.matrix(ranges(matchres.CT)))
      c(matchnum.IP,matchnum.CT)
    }
    nIP=length(which(y[1,]>0))
    nCT=length(which(y[2,]>0))
    test.mat=c(nIP,nsample-nIP,nCT,nsample-nCT)
    teststat=fisher.test(matrix(test.mat,nrow=2,byrow=TRUE))
    nIPCT=ifelse(nIP>nCT,nIP-nCT,1)
    result=c(allmotifs.name[i],test.mat,teststat$p.value,teststat[[3]],nIPCT)
    result
  }
  write.table(as.data.frame(Z),file=outfile,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
  
}




fmotifscanbatch<-function(seqIPfile,seqCTfile,seqdir,motifdir,database=c('JASPAR','TRANSFAC'),score="80%",ncore=1){
  seqIPs=as.matrix(read.table(seqIPfile))
  seqCTs=as.matrix(read.table(seqCTfile))
  for(i in 1:length(seqIPs)){
    tempfile=strsplit(basename(seqIPs[i]),split=".fa.upper")[[1]][1]
    outfile=paste(motifdir,"/",tempfile,"-count",sep="")
    fmotifscan(file.path(seqdir,basename(seqIPs[i])),file.path(seqdir,basename(seqCTs[i])),outfile,database,score,ncore)
  }
}




fpeakscanbatch<-function(seqIPfile,seqCTfile,seqdir,peakfile,peakdir,database=c('JASPAR','TRANSFAC'),score="80%",ncore=1){
  seqIPs=as.matrix(read.table(seqIPfile))
  seqCTs=as.matrix(read.table(seqCTfile))
  peaks=as.matrix(read.table(peakfile))
  regions=NULL
  for(i in 1:length(peaks)){
    tempfile=strsplit(basename(seqIPs[i]),split=".fa.upper")[[1]][1]
    temp=fpeakscan(file.path(seqdir,basename(seqIPs[i])),file.path(seqdir,basename(seqCTs[i])),
              file.path(peakdir,basename(peaks[i])),database,score,ncore)
    temp=cbind(temp,rep(peaks[i],nrow(temp)))
    regions=rbind(regions,temp)
  }
  colnames(regions)=c('chr','start','end','motif','sample')
  regions
}
  

fpeakscan<-function(seqIP,seqCT,peak,database=c('JASPAR','TRANSFAC'),score="80%",ncore=1){
  
  registerDoParallel(cores=ncore)
  if(database=='JASPAR'){
    #data(JASPAR)
    load('../JASPAR.rda')
    allmotifs=JASPAR
  }else if(database=='TRANSFAC'){
    #data(TRANSFAC)
    load('../TRANSFAC.rda')
    allmotifs=TRANSFAC
  }
  allmotifs.name=names(allmotifs)
  motif.num=length(allmotifs)
  
  seqIP=readDNAStringSet(seqIP, format="fasta")
  seqIPlist=as.list(seqIP)
  seqCT=readDNAStringSet(seqCT, format="fasta")
  seqCTlist=as.list(seqCT)
  nsample=length(seqIPlist)
  regions=read.table(peak)


  Z=foreach(i=1:length(allmotifs),.combine='rbind') %dopar%{
    motif.len=ncol(allmotifs[[i]])
    y<-foreach(j=1:nsample,.combine='cbind') %dopar%{
      seqIP=toString(seqIPlist[[j]])
      matchres.IP=suppressWarnings(matchPWM(allmotifs[[i]], DNAString(seqIP), min.score=score))
      seqCT=toString(seqCTlist[[j]])
      matchres.CT=suppressWarnings(matchPWM(allmotifs[[i]], DNAString(seqCT), min.score=score))
      matchnum.IP=nrow(as.matrix(ranges(matchres.IP)))
      matchnum.CT=nrow(as.matrix(ranges(matchres.CT)))
      c(matchnum.IP,matchnum.CT)
    }
    id=which(y[1,]>0 & y[2,]==0)
    if(length(id)>0){
      res=cbind(regions[id,1:3],rep(allmotifs.name[i],length(id)))
    }else{
      res=NULL
    }
    result=res
  }
  

}










