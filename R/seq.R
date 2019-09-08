genSeq<-function(peakfile,peakdir,seqdir,genomeseq){
  files=as.matrix(read.table(peakfile))
  for(i in 1:length(files)){
    seqfile=paste(files[i],".fa",sep="")
    seqfile.upper=paste(seqfile,".upper",sep="")
    cmd=paste("bedtools getfasta -fi", genomeseq, " -bed ",file.path(peakdir,files[i])," -fo ",file.path(seqdir,seqfile))
    system(cmd)
    cmd=paste(" tr \"[:lower:]\" \"[:upper:]\" < ", file.path(seqdir,seqfile), " > ",file.path(seqdir,seqfile.upper))
    system(cmd)
  }
  cmd=paste("rm ",seqdir,"/*fa",sep="")
  system(cmd)
}










