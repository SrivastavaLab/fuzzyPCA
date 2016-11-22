## A is the abundance matrix and T is the trait matrix and dset is the data set

# Fuzzit produces for every bromeliad the traits weighted by the abundance of every species in a plant. 

Fuzzit<-function(A,T,dset)
{
  bromeliad.id<-colnames(A)
  a.fuz<-matrix(NA,nrow=dim(A)[2],ncol=dim(T)[2])	#bromeliad x  mean fuzzy coded trait matrix
  n.plant<-dim(A)[2]
  
  for(i in 1:dim(A)[2]){
    a.fuz[i,]<-colSums(A[,i]*T,na.rm=TRUE)/n.plant		# average [A]abundance weighted traits for a plant
  }
  
  colnames(a.fuz)<-colnames(T)
  a.fuz2<-data.frame(bromeliad.id,rep(dset,n.plant),a.fuz)
  names(a.fuz2)[2]<-"dset"
  return(a.fuz2)
}