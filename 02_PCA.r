#####################################################################################################
##Functionaltrait analysis
#CESAB data from BWGdb [Bromeliad Working Group data base] June 18 2016
#co-PI Regis, Diane and Kurt
library(ade4)
library(vegan)
library(gam)
library(plyr)
library(dplyr)
library(RColorBrewer)

source('functions.R')

trait_matrix_by_dataset <- readRDS("Data/trait_matrix_by_dataset.rds")
abundance_matrix_by_dataset <- readRDS("Data/abundance_matrix_by_dataset.rds")

#produces trait_matrix_by_dataset and abundance_matrix_by_dataset
#test: names(trait_matrix_by_dataset)==names(trait_matrix_by_dataset)


#run Fuzzit by dataset 
#Fuzzit produces traits weighted by abundance for every bromeliad

nset<-length(names(trait_matrix_by_dataset))
fdat<-NULL


### Selecting the abundance and the trait matrix for each data set and using Fuzzit on those matrices 

for(i in 1:nset){ 
  dset<-as.numeric(names(abundance_matrix_by_dataset[i]))
  A<-abundance_matrix_by_dataset[[i]]
  T<-trait_matrix_by_dataset[[i]]

  fdat[[i]]<-Fuzzit(A,T,dset)

}

fdat2<-bind_rows(fdat)
fdat2$dset

####group dataset into "regions" 
#see: 01_datasets.csv

## we are adding a region based on the country to each of the datasets -> all of them are in a dataframe fdat2
## SCIENCE NOTE: review region definitions

##Add region by dataset ID 
regions <- read.csv('Data/Regions.csv', header = TRUE)
fdat2 <- left_join(fdat2,regions) %>% select(-dataset)


#BWGdb 
# We need to define groups for the traits. Within each trait there are a x number of groups
groups<-c(4,5,2,8,6,7,8,4,5,8,3,4)
names(groups)=c("AS","BS","DM","FD","FG","LO","RE","RF","RM","MD","CP","BF")##noms des traits

reg1<-as.numeric(as.factor(sort(unique(fdat2$region))))
names(reg1)<-sort(unique(fdat2$region))

#attention: no all-zero columns/rows, must delete!!
which(colSums(fdat2[,3:(dim(fdat2)[2]-1)])==0)
which(rowSums(fdat2[,3:(dim(fdat2)[2]-1)])==0) 

postscript("Figures/Ordination.ps",horizontal=TRUE)
par(mfrow=c(1,1),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)

#display.brewer.pal(9,"Paired")  #start with this and repeat for county but not within country
mycol=c("#A6CEE3","#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#E31A1C","#E31A1C","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6")

#Fuzzy Principal Correspondence analysis (FPCA)
data_fpca <- prep.fuzzy.var(fdat2[,3:(dim(fdat2)[2]-1)], groups) #code trait modalities 
res_fpca <- dudi.fpca(data_fpca, scannf = FALSE, nf = 2)
#FPCA plots
#scatter(res_fpca)#samples and traits
#s.arrow(res_fpca$l1)#samples
s.arrow(res_fpca$c1)#traits
summary(res_fpca)#Eigenvalues, projected inertia (by axis and cumulated)
s.class(res_fpca$li,as.factor(fdat2$region),cpoint=0,cstar=0,axesell=TRUE,col=mycol)	#show clusters by region
##legend("bottomleft", pch=16, col=brewer.pal(17,"Spectral"), legend=unique(as.factor(fdat3$region)))

dev.off()
  
saveRDS(res_fpca, "Data/PCA_results.rds")
write.csv(fdat2, 'Data/Abundance_weighted_traits.csv', row.names = FALSE)
