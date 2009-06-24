"gene.pvalue"<-function(EXP,measure,net.trim,n.replica=400){
perm<-function(data){
  return(t(apply(data,1,FUN="sample",replace=FALSE)))
  }
repli.matrix<-function(EXP,measure,net.trim,n.replica){
  ut <- replicate(n.replica,gene.similarity(perm(EXP),measure,net.trim))
  out <- array(ut,c(nrow(EXP),nrow(EXP),n.replica))
  }
beta.estimate<-function(data){
  me <- mean(as.vector(data))
  va <- var(as.vector(data))
  alpha <- me*(me*(1-me)/va-1)
  beta <- (1-me)*(me*(1-me)/va-1)
  return(list(shape1=alpha,shape2=beta))
  }
betatail<-function(data,x){
  if(x>max(data)){
    para <- beta.estimate(data)
    return(1-pbeta(x,shape1=para$shape1,shape2=para$shape2))
    }else{
    return(sum(data>x)/length(data))
    }
  }
if((measure!="MI")&&(measure!="corr")) stop("measure is not correct")
if(measure=="corr"){
    fi=function(i,EXP)
      {
      fj=function(j,i,EXP)
        {
        return(cor.test(EXP[i,],EXP[j,])$p.value)
        }
      return(apply(matrix(1:nrow(EXP),nc=1),1,fj,i,EXP))
      }
    out.single=t(apply(matrix(1:nrow(EXP),nc=1),1,fi,EXP))
  out.corrected=NULL
  }
  else{
    real <- gene.similarity(EXP,measure,net.trim)
    rep <- repli.matrix(EXP,measure,net.trim,n.replica)
    per=rep
    fi=function(i,per,real)
      {
      fj=function(j,i,per,real)
        {
        return(betatail(abs(c(per[i,j,],per[j,i,])),abs(real[i,j])))
        }
      return(apply(matrix(1:ncol(real),nc=1),1,fj,i,per,real))
      }
      out.single=t(apply(matrix(1:nrow(real),nc=1),1,fi,per,real))
        ## to delete the diagonal element in calculation of corrected p-values
    fi1=function(i,per,real)
      {
      fj1=function(j,i,per,real)
        {
        temp <- abs(c(per[i,,],per[,j,]))
        if(length(temp[temp!=1])==0) stop("corrected permutation method is not appropiate in this case")
       retn<- betatail(temp[temp!=1],abs(real[i,j]))
       if(is.na(retn)) retn=0
        return(retn)
        }
      return(apply(matrix(1:ncol(real),nc=1),1,fj1,i,per,real))
      }
     out.corrected=t(apply(matrix(1:nrow(real),nc=1),1,fi1,per,real))   
    }
 mm.corrected <- MM.correction(out.single)
 return(list(single.perm.p.value=out.single,multi.perm.p.value=out.corrected,MMcorrected.p.value=mm.corrected))
}


