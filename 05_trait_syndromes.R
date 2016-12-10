library()

res.t <- readRDS( "Data/clustering_traits.rds")
res <- readRDS("Data/clustering_broms.rds")
fdat2 <- read.csv('Data/Abundance_weighted_traits.csv')
res_fpca <- readRDS('Data/PCA_results.rds')


##############################################################################
###traits syndromes (clusters of traits)

names(sort(res.t$Best.partition[res.t$Best.partition==1]))
names(sort(res.t$Best.partition[res.t$Best.partition==2]))
names(sort(res.t$Best.partition[res.t$Best.partition==3]))
names(sort(res.t$Best.partition[res.t$Best.partition==4]))
names(sort(res.t$Best.partition[res.t$Best.partition==5]))
names(sort(res.t$Best.partition[res.t$Best.partition==6]))
names(sort(res.t$Best.partition[res.t$Best.partition==7]))

tsyn1<-c("AS1", "AS2", "BS4", "FD3", "FG2", "FG3", "LO4", "LO6", "RE4", "RE7", "RM3", "CP2")
tsyn2<-c("AS3", "BS2", "FD4", "FG1", "FG4", "LO2", "LO3", "RE2", "MD8")
tsyn3<-c( "AS4", "BS3", "BS5", "FD1", "FD2", "FD5", "FD6", "FD7", "FD8", "FG5", "FG6", "LO1", "LO5",
          "LO7", "RE1", "RE3", "RE5", "RE6", "RE8", "RF1", "RF2", "RM2", "RM5", "MD1", "MD2", "MD4",
          "MD5", "BF1")
tsyn4<-c("BS1", "DM1", "RM1", "MD7", "BF4")
tsyn5<-c("DM2", "RM4")
tsyn6<-c("RF3", "MD3", "CP1", "BF3")
tsyn7<-c("RF4", "MD6", "CP3", "BF2")


##############################################################################
###traits of plant that cluster together
fdat2b<-fdat2
fdat2b[fdat2b==0]<-NA

fdat3<-split(fdat2b,sort(res$Best.partition))

a<-NULL
par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn1]
  #boxplot(log(temp))
  #abline(h=0,lty=8)
  a[[i]]<-boxplot(log(temp),plot=FALSE)[[1]]
  #mtext(paste("cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
#mtext("Trait syndrom 1",side=3,line=0.5,outer=TRUE)

par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn1]
  temp[,which(colSums(sign(a[[i]][c(2,4),]))==0)]<-1
  boxplot(log(temp),plot=TRUE)
  abline(h=0,lty=8)
  mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
mtext("log(AWT)",side=2,line=0.5,outer=TRUE,las=0)
mtext("Trait syndrom 1",side=3,line=0.5,outer=TRUE)

a<-NULL
par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn2]
  #boxplot(log(temp))
  #abline(h=0,lty=8)
  a[[i]]<-boxplot(log(temp),plot=FALSE)[[1]]
  #mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
#mtext("Trait syndrom 2",side=3,line=0.5,outer=TRUE)

par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn2]
  temp[,which(colSums(sign(a[[i]][c(2,4),]))==0)]<-1
  boxplot(log(temp),plot=TRUE)
  abline(h=0,lty=8)
  mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
mtext("log(AWT)",side=2,line=0.5,outer=TRUE,las=0)
mtext("Trait syndrom 2",side=3,line=0.5,outer=TRUE)

a<-NULL
par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn3]
  #boxplot(log(temp))
  #abline(h=0,lty=8)
  a[[i]]<-boxplot(log(temp),plot=FALSE)[[1]]
  #mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
#mtext("Trait syndrom 3",side=3,line=0.5,outer=TRUE)

par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn3]
  temp[,which(colSums(sign(a[[i]][c(2,4),]))==0)]<-1
  boxplot(log(temp),plot=TRUE)
  abline(h=0,lty=8)
  mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
mtext("log(AWT)",side=2,line=0.5,outer=TRUE,las=0)
mtext("Trait syndrom 3",side=3,line=0.5,outer=TRUE)

a<-NULL
par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn4]
  #boxplot(log(temp))
  #abline(h=0,lty=8)
  a[[i]]<-boxplot(log(temp),plot=FALSE)[[1]]
  #mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
#mtext("Trait syndrom 4",side=3,line=0.5,outer=TRUE)

par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn4]
  temp[,which(colSums(sign(a[[i]][c(2,4),]))==0)]<-1
  boxplot(log(temp),plot=TRUE)
  abline(h=0,lty=8)
  mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
mtext("log(AWT)",side=2,line=0.5,outer=TRUE,las=0)
mtext("Trait syndrom 4",side=3,line=0.5,outer=TRUE)

a<-NULL
par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn5]
  #boxplot(log(temp))
  #abline(h=0,lty=8)
  a[[i]]<-boxplot(log(temp),plot=FALSE)[[1]]
  #mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
#mtext("Trait syndrom 5",side=3,line=0.5,outer=TRUE)

par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn5]
  temp[,which(colSums(sign(a[[i]][c(2,4),]))==0)]<-1
  boxplot(log(temp),plot=TRUE)
  abline(h=0,lty=8)
  mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
mtext("log(AWT)",side=2,line=0.5,outer=TRUE,las=0)
mtext("Trait syndrom 5",side=3,line=0.5,outer=TRUE)

a<-NULL
par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn6]
  #boxplot(log(temp))
  #abline(h=0,lty=8)
  a[[i]]<-boxplot(log(temp),plot=FALSE)[[1]]
  #mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
#mtext("Trait syndrom 6",side=3,line=0.5,outer=TRUE)

par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn6]
  temp[,which(colSums(sign(a[[i]][c(2,4),]))==0)]<-1
  boxplot(log(temp),plot=TRUE)
  abline(h=0,lty=8)
  mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
mtext("log(AWT)",side=2,line=0.5,outer=TRUE,las=0)
mtext("Trait syndrom 6",side=3,line=0.5,outer=TRUE)

a<-NULL
par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn7]
  #boxplot(log(temp))
  #abline(h=0,lty=8)
  a[[i]]<-boxplot(log(temp),plot=FALSE)[[1]]
  #mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
#mtext("Trait syndrom 7",side=3,line=0.5,outer=TRUE)

par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn7]
  temp[,which(colSums(sign(a[[i]][c(2,4),]))==0)]<-1
  boxplot(log(temp),plot=TRUE)
  abline(h=0,lty=8)
  mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
mtext("log(AWT)",side=2,line=0.5,outer=TRUE,las=0)
mtext("Trait syndrom 7",side=3,line=0.5,outer=TRUE)



################################################################################
################################################################################
#top 10 functional traits based on loadings
res_fpca$co$dist<-sqrt((res_fpca$co$Comp1 - res_fpca$co$Comp2) ^ 2)		#Euclidian distance; length of arrow on fig 1?
res_fpca$co<-res_fpca$co[rev(order(res_fpca$co$dist)),]					#order; biggest loading to smallest 

T10<-rownames(res_fpca$co)[1:10]										#select top 10 loadings. 

# [1] "DM2" "DM1" "BF2" "RM1" "MD6" "CP1" "MD8" "BF3" "CP3" "RM4"



tsyn2<-c("MD8")
tsyn4<-c("DM1", "RM1")
tsyn5<-c("DM2", "RM4")
tsyn6<-c("CP1", "BF3")
tsyn7<-c("CP3","BF2","MD6")


## copy and paste the boxplots



a<-NULL
par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn2]
  #boxplot(log(temp))
  #abline(h=0,lty=8)
  a[[i]]<-boxplot(log(temp),plot=FALSE)[[1]]
  #mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
#mtext("Trait syndrom 2",side=3,line=0.5,outer=TRUE)

par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn2]
  temp[,which(colSums(sign(a[[i]][c(2,4),]))==0)]<-1
  boxplot(log(temp),plot=TRUE)
  abline(h=0,lty=8)
  mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
mtext("log(AWT)",side=2,line=0.5,outer=TRUE,las=0)
mtext("Trait syndrom 2",side=3,line=0.5,outer=TRUE)

a<-NULL
par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn4]
  #boxplot(log(temp))
  #abline(h=0,lty=8)
  a[[i]]<-boxplot(log(temp),plot=FALSE)[[1]]
  #mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
#mtext("Trait syndrom 4",side=3,line=0.5,outer=TRUE)

par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn4]
  temp[,which(colSums(sign(a[[i]][c(2,4),]))==0)]<-1
  boxplot(log(temp),plot=TRUE)
  abline(h=0,lty=8)
  mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
mtext("log(AWT)",side=2,line=0.5,outer=TRUE,las=0)
mtext("Trait syndrom 4",side=3,line=0.5,outer=TRUE)

a<-NULL
par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn5]
  #boxplot(log(temp))
  #abline(h=0,lty=8)
  a[[i]]<-boxplot(log(temp),plot=FALSE)[[1]]
  #mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
#mtext("Trait syndrom 5",side=3,line=0.5,outer=TRUE)

par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn5]
  temp[,which(colSums(sign(a[[i]][c(2,4),]))==0)]<-1
  boxplot(log(temp),plot=TRUE)
  abline(h=0,lty=8)
  mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
mtext("log(AWT)",side=2,line=0.5,outer=TRUE,las=0)
mtext("Trait syndrom 5",side=3,line=0.5,outer=TRUE)

a<-NULL
par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn6]
  #boxplot(log(temp))
  #abline(h=0,lty=8)
  a[[i]]<-boxplot(log(temp),plot=FALSE)[[1]]
  #mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
#mtext("Trait syndrom 6",side=3,line=0.5,outer=TRUE)

par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn6]
  temp[,which(colSums(sign(a[[i]][c(2,4),]))==0)]<-1
  boxplot(log(temp),plot=TRUE)
  abline(h=0,lty=8)
  mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
mtext("log(AWT)",side=2,line=0.5,outer=TRUE,las=0)
mtext("Trait syndrom 6",side=3,line=0.5,outer=TRUE)

a<-NULL
par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn7]
  #boxplot(log(temp))
  #abline(h=0,lty=8)
  a[[i]]<-boxplot(log(temp),plot=FALSE)[[1]]
  #mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
#mtext("Trait syndrom 7",side=3,line=0.5,outer=TRUE)

par(mfrow=c(3,2),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
for(i in 1:length(fdat3)){
  temp<-fdat3[[i]][,tsyn7]
  temp[,which(colSums(sign(a[[i]][c(2,4),]))==0)]<-1
  boxplot(log(temp),plot=TRUE)
  abline(h=0,lty=8)
  mtext(paste("plant cluster ",i,sep=''),side=3,line=0.5,outer=FALSE)
}
mtext("log(AWT)",side=2,line=0.5,outer=TRUE,las=0)
mtext("Trait syndrom 7",side=3,line=0.5,outer=TRUE)


