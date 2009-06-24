"gene.trait.similarity"<-function(EXP,trait,measure,data.type=c("continuous","continuous"),na.replica=50){
if((data.type[1]!="continuous")&&(data.type[1]!="discrete")) stop("Wrong data type for gene expression data") 
if((data.type[2]!="continuous")&&(data.type[2]!="discrete")) stop("Wrong data type for trait data") 
complete<-function(x,smooth){
  x.all <- x[!is.na(x)]
  n.missing <- sum(is.na(x))
  if(smooth){
    x[is.na(x)] <- as.numeric(sample(x.all,n.missing,replace=FALSE))+runif(n.missing,-0.001,0.001)
  }else{
    x[is.na(x)] <- sample(x.all,n.missing,replace=FALSE)
  }
  return(x)
  }
corr.both<-function(mat,nat,data.type){
    mat <- apply(mat,1,FUN="complete",smooth=(data.type[1]=="continuous"))
    nat <- apply(nat,1,FUN="complete",smooth=(data.type[2]=="continuous"))
    return(cor(mat,nat))
  }
mi.both<-function(mat,nat,data.type){
    mat <- t(apply(mat,1,FUN="complete",smooth=(data.type[1]=="continuous")))
    nat <- t(apply(nat,1,FUN="complete",smooth=(data.type[2]=="continuous")))
    mi.core <- function(x,y,data.type){
      len <- length(x)
      n.bin <- round(sqrt(len))
      if(data.type[1]=="continuous"){
        xo <- ceiling((rank(x))/n.bin)
      }else{
        xo <- x
      }
      if(data.type[2]=="continuous"){
        yo <- ceiling((rank(y))/n.bin)
      }else{
        yo <- y
      }
      xt <- table(xo)/len
      yt <- table(yo)/len
      ## joint distribution for x and y
      xyt <- table(xo,yo)/len
      xyt <- xyt[xyt!=0]
      ## entropy for x and y
      x.h <- sum(xt*log(xt))
      y.h <- sum(yt*log(yt))
      ## joint entropy
      xy.h <- sum(xyt*log(xyt))
      ## standardized by entropy for y so that MI values are within [0,1]
      return((xy.h-x.h-y.h)/min(c(-y.h,-x.h)))
    }
    out.mi <- function(i,j,data.type){
    return(mi.core(mat[i,],nat[j,],data.type))
    }
    v <- Vectorize(out.mi,vectorize.args=c("i","j"))
    return(outer(1:nrow(mat),1:nrow(nat),FUN="v",data.type=data.type))
  }
b.similarity=function(mat,nat,measure,data.type){
  if(measure=="MI"){
       return(mi.both(mat,nat,data.type))
   }
    else{
    if(measure=="corr"){
      return(corr.both(mat,nat,data.type))}
      else{
         stop("the measure input is not available")
        }
    }
  }
if((measure!="MI")&&(measure!="corr")) stop("measure is not correct")
if(any(c(is.na(EXP),is.na(trait)))){
  ut <- replicate(na.replica,b.similarity(EXP,trait,measure,data.type))
  out <- apply(array(ut,c(nrow(EXP),nrow(trait),na.replica)),c(1,2),mean)
  return(out)
  }else{
    return(b.similarity(EXP,trait,measure,data.type))
  }
}