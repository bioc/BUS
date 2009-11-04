`gene.trait.similarity` <-
function(EXP,trait,measure,na.replica=50){

gene.trait.min<-function(mat,nat)
{
	n1 <- ncol(mat)
	n2 <- ncol(nat)
	N <- nrow(mat)
	
	mat[which(is.na(mat))] <- -2000000
	nat[which(is.na(nat))] <- -2000000
    mi.vector <- .Call( "MINempirical", mat, nat, N, n1, n2,DUP=FALSE)
    dim(mi.vector) <- c(n1,n2)
    mi.matrix <- as.matrix(mi.vector)
    return(mi.matrix)

}

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
	
    mat <- apply(mat,1,FUN="complete",smooth=(data.type[1]=="continuous"))
    nat <- apply(nat,1,FUN="complete",smooth=(data.type[2]=="continuous"))
    
    if(data.type[1]=="continuous")
      	mat<-as.matrix(discretize(mat))
    else
        mat<-as.matrix(mat)
    if(data.type[2]=="continuous")
      	nat <-as.matrix(discretize(nat))
    else
        nat<-as.matrix(nat)
 
    out.mi<-gene.trait.min(mat, nat)
    
    return(norm(out.mi))
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
  
# data type recognition
data.type.exp<-"continuous"
if(!is.integer(EXP))
data.type.exp<-"continuous"
else
data.type.exp<-"discrete"
if(!is.integer(trait))
data.type.trait<-"continuous"
else
data.type.trait<-"discrete"
data.type<-c(data.type.exp,data.type.trait)

gene.names<-rownames(EXP)
trait.names<-rownames(trait)

if((measure!="MI")&&(measure!="corr")) stop("measure is not correct")
if(any(c(is.na(EXP),is.na(trait))))
{
  ut <- replicate(na.replica,b.similarity(EXP,trait,measure,data.type))
  out <- apply(array(ut,c(nrow(EXP),nrow(trait),na.replica)),c(1,2),mean)

  }
  else
  {
  	out<-b.similarity(EXP,trait,measure,data.type)

    
  }
dimnames(out)<-list(gene.names,trait.names)
return(out)
}

