topPeak<-function(rawpeakfile,peakdir='peaks',ntop=500,ext=500){
  files=as.matrix(read.table(rawpeakfile))
  for(i in 1:length(files)){
    tempmat=read.table(file.path(files[i]))
    tempmat.sort=tempmat[order(tempmat[,8],decreasing=TRUE),] #broadPeak column8 is -logpvalue #broadPeak column5 score
    if(nrow(tempmat.sort)>ntop){
      tempmat.sort=tempmat.sort[c(1:ntop),]
    }else if(nrow(tempmat.sort)<200){
      next
    }
    tempmat.mid=floor((tempmat.sort[,2]+tempmat.sort[,3])/2)
    tempmat.summit=data.frame(chr=tempmat.sort[,1],start=tempmat.mid-ext,end=tempmat.mid+ext,id=c(1:nrow(tempmat.sort)))
    tempmat.summit[,2]=str_trim(format(tempmat.summit[,2], scientific = FALSE),side="both")
    tempmat.summit[,3]=str_trim(format(tempmat.summit[,3], scientific = FALSE),side="both")
    tempfile.sort=paste(peakdir,"/",basename(files[i]),".topIP",sep="")
    write.table(tempmat.summit,file=tempfile.sort,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
  }
}



CTPeak<-function(topIPpeakfile,peakdir,ext=10^5){

  Hmat=as.matrix(seqlengths(Hsapiens))
  Hmat=data.frame(chr=rownames(Hmat),len=Hmat)
  files=as.matrix(read.table(topIPpeakfile))
  
  for(i in 1:length(files)){
    CTpeakfile=paste(strsplit(files[i],".topIP")[[1]],".topCT",sep="")
    mat=read.table(file.path(peakdir,files[i]))
    chr.order=match(mat[,1],Hmat[,1])
    len.order=Hmat[,2][chr.order]
    pos1=which(mat[,2]>ext)
    if(length(pos1)>0){
      mat[pos1,][,2]=mat[pos1,][,2]-ext
      mat[pos1,][,3]=mat[pos1,][,3]-ext
    }
    pos2=which(mat[,2]<=ext)
    pos3=which((mat[,3]+ext)<len.order)
    pos4=intersect(pos2,pos3)	
    if(length(pos4)>0){
      mat[pos4,][,2]=mat[pos4,][,2]+ext
      mat[pos4,][,3]=mat[pos4,][,3]+ext
    }
    mat=as.data.frame(mat)
    mat[,2]=str_trim(format(mat[,2], scientific = FALSE),side="both")
    mat[,3]=str_trim(format(mat[,3], scientific = FALSE),side="both")
    write.table(mat,file=paste(peakdir,"/",CTpeakfile,sep=""),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
  }
}



