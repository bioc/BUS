"BUS"<-function(EXP,trait=NULL,measure,method.permut=2,method.correct,n.replica=400,net.trim=NULL,thresh=NULL,nflag) {
if((nflag!=1)&&(nflag!=2)) stop("nflag is out of bound")
if(nflag==1){
  if((thresh>1)||(thresh<0)) stop("thresh is out of bound")
	if((net.trim!="mrnet")&&(net.trim!="aracne")&&(net.trim!="clr")&&(net.trim!="none")) stop("net.trim is out of bound")
  }
if((measure!="MI")&&(measure!="corr")) stop("measure is not correct")
if((method.correct!="MMcorrection")&&(method.correct!="permutation")&&(method.correct!="both")) stop("method.correct is not correct")
if(!is.matrix(EXP)) stop("EXP should be a matrix")
if(nflag==2){
  if(!is.numeric(trait)) stop("trait should be a numeric matrix in this case")
  if(!is.matrix(trait)) stop("trait should be a matrix")
  }
if(nflag==2){
if((method.permut!=1)&&(method.permut!=2)&&(method.permut!=3)) stop("method.permut is out of bound")
}   
if(nflag==1)
{
  MI=gene.similarity(EXP,measure,net.trim)
  rlist=gene.pvalue(EXP,measure,net.trim,n.replica)
  single.p<-rlist$single.perm.p.value
  fisher.p<-rlist$MMcorrected.p.value
  corrected.p<-rlist$multi.perm.p.value
  RE.fisher=RE.corr=NULL
  if(method.correct=="MMcorrection")
  {
    RE.fisher<-pred.network(fisher.p,MI,thresh)
    return(list(single.perm.p.value=single.p,MMcorrected.p.value=fisher.p,net.pred.MMcorrected=RE.fisher))
   }
  if(method.correct=="permutation")
  {
     RE.corr<-pred.network(corrected.p,MI,thresh)
     return(list(single.perm.p.value=single.p,multi.perm.p.value=corrected.p,net.pred.permut=RE.corr))
  }
  if(method.correct=="both")
  {
    RE.fisher<-pred.network(fisher.p,MI,thresh)
    RE.corr<-pred.network(corrected.p,MI,thresh)
    return(list(single.perm.p.value=single.p,multi.perm.p.value=corrected.p,MMcorrected.p.value=fisher.p,net.pred.permut=RE.corr,net.pred.MMcorrected=RE.fisher))
  }
}
else{
    n.trait=nrow(trait)
    MI=gene.trait.similarity(EXP,trait,measure)
    rlist=gene.trait.pvalue(EXP,trait,measure,method.permut,n.replica)
    single.p<-rlist$single.perm.p.value
    fisher.p<-rlist$MMcorrected.p.value
    corrected.p<-rlist$multi.perm.p.value
    if(method.correct=="MMcorrection")
     {
      return(list(single.perm.p.value=single.p,MMcorrected.p.value=fisher.p))
     }
    if(method.correct=="permutation")
     {
       return(list(single.perm.p.value=single.p,multi.perm.p.value=corrected.p))
     }
    if(method.correct=="both")
     {
      return(list(single.perm.p.value=single.p,multi.perm.p.value=corrected.p,MMcorrected.p.value=fisher.p))
     }
  }
}


