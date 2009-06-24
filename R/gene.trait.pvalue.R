"gene.trait.pvalue"<-function(EXP,trait,measure,method.permut=2,n.replica=400){
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
perm.trait<-function(data){
    return(t(apply(data,1,FUN="sample",replace=FALSE)))
  }
repli.matrix.both<-function(EXP,trait,measure,n.replica){
    ut <- replicate(n.replica,gene.trait.similarity(EXP,perm.trait(trait),measure))
    out <- array(ut,c(nrow(EXP),nrow(trait),n.replica))
  }
ecd <- function(x,distro,use.beta){
  return(betatail(abs(distro),abs(x)))
}
if((measure!="MI")&&(measure!="corr")) stop("measure is not correct")
if((method.permut!=1)&&(method.permut!=2)&&(method.permut!=3)) stop("method.permut is out of bound")
if(measure=="corr"){
  out.single=matrix(0,nr=nrow(EXP),nc=nrow(trait))
  fi=function(i,EXP,trait)
    {
    fj=function(j,i,EXP,trait)
     {
     return(cor.test(EXP[i,],trait[j,])$p.value)
     }
    return(apply(matrix(1:nrow(trait),nc=1),1,fj,i,EXP,trait))
    }
  out.single=t(apply(matrix(1:nrow(EXP),nc=1),1,fi,EXP,trait))
  out.corrected=NULL
  }
  else{
    real <- gene.trait.similarity(EXP,trait,measure)
    perm<- repli.matrix.both(EXP,trait,measure,n.replica)
    correction=method.permut;use.beta=TRUE
  fi=function(i,real,perm,use.beta)
    {
    fj=function(j,i,real,perm,use.beta)
     {
     return(ecd(real[i,j],perm[i,j,],use.beta))
     }
    return(apply(matrix(1:ncol(real),nc=1),1,fj,i,real,perm,use.beta))
    }
  out.single=t(apply(matrix(1:nrow(real),nc=1),1,fi,real,perm,use.beta))
  if(correction==3){
     cum <- perm
     fi=function(i,real,cum,use.beta)
       {
        fj=function(j,i,real,cum,use.beta)
         {
         return(ecd(real[i,j],cum,use.beta))
         }
       return(apply(matrix(1:ncol(real),nc=1),1,fj,i,real,cum,use.beta))
      }
     out.corrected=t(apply(matrix(1:nrow(real),nc=1),1,fi,real,cum,use.beta))
    }
  if(correction==1){
       fi=function(i,real,perm,use.beta)
       {
        fj=function(j,i,real,cum,use.beta)
         {
         return(ecd(real[i,j],cum,use.beta))
         }
       cum=perm[i,,]
       return(apply(matrix(1:ncol(real),nc=1),1,fj,i,real,cum,use.beta))
      }
     out.corrected=t(apply(matrix(1:nrow(real),nc=1),1,fi,real,perm,use.beta))
    }
  if(correction==2){
       fj=function(j,real,perm,use.beta)
       {
        fi=function(i,j,real,cum,use.beta)
         {
         return(ecd(real[i,j],cum,use.beta))
         }
       cum=perm[,j,]
       return(apply(matrix(1:nrow(real),nc=1),1,fi,j,real,cum,use.beta))
       }
     out.corrected=apply(matrix(1:ncol(real),nc=1),1,fj,real,perm,use.beta)
    }
  }
out.MM.corrected <- MM.correction(out.single)
return(list(single.perm.p.value=out.single,multi.perm.p.value=out.corrected,MMcorrected.p.value=out.MM.corrected))
}

