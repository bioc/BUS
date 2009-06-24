"MM.correction"<-function(data){
fisher.core<-function(x){
  ran <- order(-x)
  out <- numeric(length(x))
  for(i in 1:length(x)){
    fisher.inv <- -2*sum(log(x[ran[1:i]]))
    out[ran[i]] <- 1-pchisq(fisher.inv,df=2*i)
    }
  return(out)
  }
out <- matrix(rep(0,length(data)),dim(data))  
for(i in 1:nrow(data)){
  for(j in 1:ncol(data)){
    if(i!=j){
      extract <- matrix(rep(FALSE,length(data)),dim(data))
      extract[i,] <- TRUE
      extract[,j] <- TRUE
      diag(extract) <- FALSE
      extract[i,j] <- FALSE
      temp <- c(data[i,j],data[extract])
      out[i,j] <- fisher.core(temp)[1]
    }
   }
  }
return(out)
}

