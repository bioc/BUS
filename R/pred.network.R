"pred.network"<-function(pM,similarity,thresh){
if(is.numeric(pM)){
Ttestcorr<-function(pM,similarity,alpha){
  Ttest<-similarity
  w=as.vector(which(pM>=alpha))
  Ttest[w]<-0
  return(Ttest)
  }
M<-Ttestcorr(pM,similarity,thresh)
}
else{M<-NULL}
return(M)
}

