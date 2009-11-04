`BUS` <-
function(EXP,trait=NULL,measure,method.permut=2,n.replica=400,net.trim=NULL,thresh=NULL,nflag) {

if(nflag==2)
{
	message("This is a 'gene-trait' case. If a discrete trait data is tested, each element in 'trait' should be a really 'integer'.")
	}
if((nflag!=1)&&(nflag!=2)) stop("nflag is out of bound")
 
if(nflag==1){
  if((thresh>1)||(thresh<0)) stop("thresh is out of bound")
	if((net.trim!="mrnet")&&(net.trim!="aracne")&&(net.trim!="clr")&&(net.trim!="none")) stop("net.trim is out of bound")
  }
if((measure!="MI")&&(measure!="corr")) stop("measure is not correct")

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
  rlist=gene.pvalue(EXP,measure, net.trim,n.replica)
  single.p<-rlist$single.perm.p.value
  corrected.p<-rlist$multi.perm.p.value
 
  RE.corr<-pred.network(corrected.p,MI,thresh)
 
  return(list(similarity=MI,single.perm.p.value=single.p,multi.perm.p.value=corrected.p,net.pred.permut=RE.corr))

}
else
{
    n.trait=nrow(trait)
    MI=gene.trait.similarity(EXP,trait,measure)
    rlist=gene.trait.pvalue(EXP,trait,measure, method.permut,n.replica)
    single.p<-rlist$single.perm.p.value
    corrected.p<-rlist$multi.perm.p.value

    return(list(similarity=MI,single.perm.p.value=single.p,multi.perm.p.value=corrected.p))
 
  }

}

