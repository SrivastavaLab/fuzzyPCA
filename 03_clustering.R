library(NbClust)
library(RColorBrewer)

res_fpca <- readRDS('Data/PCA_results.rds')
fdat2 <- read.csv('Data/Abundance_weighted_traits.csv')

################################################################################
#cluster analysis
par(mfrow=c(1,1),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
hc1<-hclust(d = dist(res_fpca$li), method = "ward.D")
plot(hc1,cex=.6,labels=FALSE)
rect.hclust(hc1,k=4)
c1<-rect.hclust(hc1,k=4)

par(mfrow=c(1,1),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
hc1<-hclust(d = dist(res_fpca$c1), method = "ward.D")
plot(hc1,cex=.6,labels=rownames(res_fpca$c1))
rect.hclust(hc1,k=4)
c1<-rect.hclust(hc1,k=4)


################################## Using NbClust

## Coordinates of the bromeliads clustered 

set.seed(1)
res<-NbClust(res_fpca$li, distance = "euclidean", min.nc=2, max.nc=10,
             method = "ward.D", index = "all")
res$All.index
res$Best.nc
res$All.CriticalValues
res$Best.partition

postscript("Figures/Bromeliad_coordinates_cluster.ps",horizontal=TRUE)

par(mfrow=c(1,1),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
newdf<-data.frame(res_fpca$li,as.factor(fdat2$region),res$Best.partition)
names(newdf)<-c("Axis1","Axis2","region","clust")
mycol2<-brewer.pal(5,"Paired")
s.class(res_fpca$li,as.factor(res$Best.partition),cpoint=0,cstar=0,axesell=TRUE,col=mycol2)	
points(newdf$Axis1,newdf$Axis2,col=newdf$region)
legend(min(newdf$Axis1),min(newdf$Axis1),inset=c(-0.25,0),fill=mycol2,legend=newdf$region,title="Region", horiz=TRUE)

dev.off()



###### Coordinates of the traits clustered

res.t<-NbClust(res_fpca$co, distance = "euclidean", min.nc=2, max.nc=10,
               method = "ward.D", index = "all")
res.t$All.index
res.t$Best.nc
res.t$All.CriticalValues
res.t$Best.partition

postscript("Figures/Trait_coordinates_cluster.ps",horizontal=TRUE)

newdf.t<-data.frame(res_fpca$co,res.t$Best.partition)
names(newdf.t)<-c("CS1","CS2","clust")
mycol2<-brewer.pal(7,"Paired")
s.class(res_fpca$co,as.factor(res.t$Best.partition),cpoint=0,cstar=0,axesell=TRUE,col=mycol2)	
points(newdf.t$CS1,newdf.t$CS2,pch=1)
text(newdf.t$CS1,newdf.t$CS2,labels=rownames(newdf.t), cex= 0.7, pos=3)

dev.off()
