"gene.similarity"<-function(EXP,measure,net.trim,na.replica=50){
complete <-function(x,smooth){
  x.all <- x[!is.na(x)]
  n.missing <- sum(is.na(x))
  if(smooth){
    x[is.na(x)] <- as.numeric(sample(x.all,n.missing,replace=FALSE))+runif(n.missing,-0.001,0.001)
    }else{
    x[is.na(x)] <- sample(x.all,n.missing,replace=FALSE)
    }
  return(x)
  }
s.similarity=function(x,measure,net.trim){
  mat <- apply(x,1,FUN="complete",smooth=TRUE)
  if(measure=="corr"){
      out=cor(mat)
      }
      else{
       out=build.mim(discretize(mat))
      }
  if(net.trim=="mret") {
    out=mrnet(out)
    }
    else if((net.trim=="aracne")){
      out=aracne(out)
      }
      else if(net.trim=="clr"){
        out=clr(out)
        }
        else {
          if(length(unique(as.vector(out)))>1)  out=norm(out)
          }
  diag(out)=1
  return(out)
  }
if((measure!="corr")&&(measure!="MI")) stop("measure is not correct")
if((measure!="MI")&&(net.trim!="none")) stop("the trim method is only correct under mutual information metric") 
if((net.trim!="mrnet")&&(net.trim!="aracne")&&(net.trim!="clr")&&(net.trim!="none")) stop("wrong method for net trim")
if(any(is.na(EXP))){
    ut <- replicate(na.replica,s.similarity(EXP,measure,net.trim))
    out <- apply(array(ut,c(nrow(EXP),nrow(EXP),na.replica)),c(1,2),mean)
    return(out)
  }else{
    return(s.similarity(EXP,measure,net.trim))
  }
}