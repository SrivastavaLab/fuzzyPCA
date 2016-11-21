#####################################################################################################
##[f]unctional [t]rait analysis
#CESAB data from BWGdb [Bromeliad Working Group data base] June 18 2016
#co-PI Regis, Diane and Kurt
library(ade4)
library(vegan)
library(gam)
library(plyr)
library(dplyr)

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

#####################################################################################################
#use Melissa's code to extract data
#extracting_matrices.R
#last extraction 28 June 2016

#produces trait_matrix_by_dataset and abundance_matrix_by_dataset
#test: names(trait_matrix_by_dataset)==names(trait_matrix_by_dataset)


#run Fuzzit by dataset
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
fdat2$region<-"XX"

##currently not done by visit ID
#fdat2 <- fdat2 %>%
#  mutate(region = ifelse(visit_id %in%c(141,361,156,171),"PR_df", region))%>% # dwarf forest
#  mutate(region = ifelse(visit_id %in%c(151,166,356,181,136),"PR_pc", region))%>% # Palo colorado
#  mutate(region = ifelse(visit_id %in%c(146,161,351,176,131),"PR_t", region))%>% # tabunocco
#  mutate(region = ifelse(visit_id %in%c(451),"SABA_dry", region))%>% #Saba dry forest
#  mutate(region = ifelse(visit_id %in%c(121),"SABA_montange", region))%>% #Saba Montagne
#  mutate(region = ifelse(visit_id %in%c(126),"SABA_cloud", region))%>% #Saba cloud forest
#  mutate(region = ifelse(visit_id %in%c(116),"SABA_montange", region))%>% #Saba SC montagne
#  mutate(region = ifelse(visit_id %in%c(201),"DOM_cloud", region))%>% #Dominica cloud forest
#  mutate(region = ifelse(visit_id %in%c(196),"DOM_montange", region))%>% #Dominica montagne thicket
#  mutate(region = ifelse(visit_id %in%c(191),"DOM_subtropical", region))%>% #Dominica subtropical
#  mutate(region = ifelse(visit_id %in%c(446),"Sona_1000", region))%>% #Sonadora 1000
#  mutate(region = ifelse(visit_id %in%c(376),"Sona_400", region))%>% #Sonadora 400
#  mutate(region = ifelse(visit_id %in%c(391),"Sona_400", region))%>% #Sonadora 450
#  mutate(region = ifelse(visit_id %in%c(396),"Sona_500", region))%>% #Sonadora 500
#  mutate(region = ifelse(visit_id %in%c(401),"Sona_500", region))%>% #Sonadora 550
#  mutate(region = ifelse(visit_id %in%c(506),"Sona_600", region))%>% #Sonadora 600
#  mutate(region = ifelse(visit_id %in%c(411),"Sona_600", region))%>% #Sonadora 650
#  mutate(region = ifelse(visit_id %in%c(416),"Sona_700", region))%>% #Sonadora 700
#  mutate(region = ifelse(visit_id %in%c(421),"Sona_700", region))%>% #Sonadora 750
#  mutate(region = ifelse(visit_id %in%c(426),"Sona_800", region))%>% #Sonadora 800
#  mutate(region = ifelse(visit_id %in%c(431),"Sona_800", region))%>% #Sonadora 850 
#  mutate(region = ifelse(visit_id %in%c(436),"Sona_900", region))%>% #Sonadora 900
#  mutate(region = ifelse(visit_id %in%c(441),"Sona_900", region))  #Sonadora 950
  

## we are adding a region based on the country to each of the datasets -> all of them are in a dataframe fdat2
## SCIENCE NOTE: review region definitions


####NOTE!!! DO DATA FRAME DATA SET AND REGIONS AND THEN MERGE, REPLACE CODE BELOW 

##by data set ID
fdat2 <- fdat2 %>%
  mutate(region = ifelse(dset %in%c(6),"BR_car_closed", region))%>%   #Cardoso2008	Brazil closed habitat
  mutate(region = ifelse(dset %in%c(51),"CR", region))%>%   #Pitilla1997	Costa Rica
  mutate(region = ifelse(dset %in%c(56),"CR", region))%>%   #Pitilla2000	Costa Rica
  mutate(region = ifelse(dset %in%c(61),"CR", region))%>%   #Pitilla2002	Costa Rica
  mutate(region = ifelse(dset %in%c(66),"CR", region))%>%   #Pitilla2004	Costa Rica
  mutate(region = ifelse(dset %in%c(71),"CR", region))%>%   #Pitilla2010	Costa Rica
  mutate(region = ifelse(dset %in%c(76),"CO", region))%>%   #Guasca2001	Colombia
  mutate(region = ifelse(dset %in%c(81),"CO", region))%>%   #Sisga2000	Colombia
  mutate(region = ifelse(dset %in%c(86),"CO", region))%>%   #RioBlanco2012	Colombia
  mutate(region = ifelse(dset %in%c(91),"CO", region))%>%   #RioBlanco2014	Colombia
  mutate(region = ifelse(dset %in%c(96),"JAM", region))%>%   #DiscoveryBay1997	Jamaica
  mutate(region = ifelse(dset %in%c(101),"HON", region))%>%   #Cusuco2006	Honduras
  mutate(region = ifelse(dset %in%c(106),"HON", region))%>%   #Cusuco2007	Honduras
  mutate(region = ifelse(dset %in%c(111),"SABA", region))%>%   #Saba2009	Netherlands Antilles
  mutate(region = ifelse(dset %in%c(116),"PR", region))%>%   #ElVerde2010	Puerto Rico
  mutate(region = ifelse(dset %in%c(121),"PR", region))%>%   #ElVerde1993	Puerto Rico
  mutate(region = ifelse(dset %in%c(126),"PR", region))%>%   #ElVerde1994	Puerto Rico
  mutate(region = ifelse(dset %in%c(131),"PR", region))%>%   #ElVerde1997	Puerto Rico
  mutate(region = ifelse(dset %in%c(136),"DOM", region))%>%   #Dominica2002	Dominica
  mutate(region = ifelse(dset %in%c(141),"BR_mac", region))%>%   #Macae2008	Brazil
  mutate(region = ifelse(dset %in%c(146),"BR_car_open", region))%>%   #Cardoso2011	Brazil
  mutate(region = ifelse(dset %in%c(151),"BR_pic", region))%>%   #Picinguaba2011	Brazil
  mutate(region = ifelse(dset %in%c(156),"BR_jur", region))%>%   #Jureia2013	Brazil
  mutate(region = ifelse(dset %in%c(161),"BR_s_do", region))%>%   #SerraDoJapi2011	Brazil
  mutate(region = ifelse(dset %in%c(166),"AR", region))%>%   #LasGamas2010	Argentina
  mutate(region = ifelse(dset %in%c(171),"AR", region))%>%   #LasGamas2012	Argentina
  mutate(region = ifelse(dset %in%c(176),"BR_car", region))%>%   #IlhaCardoso2011	Brazil open and closed
  mutate(region = ifelse(dset %in%c(181),"AR", region))%>%   #LasGamas2013	Argentina
  mutate(region = ifelse(dset %in%c(186),"FG_ps", region))%>%   #PetitSaut2007	French Guiana
  mutate(region = ifelse(dset %in%c(191),"FG_Kaw", region))%>%   #Kaw2008	French Guiana
  mutate(region = ifelse(dset %in%c(196),"FG_ps", region))%>%   #PetitSaut2008	French Guiana
  mutate(region = ifelse(dset %in%c(201),"FG_nour", region))%>%   #Nouragues2006	French Guiana
  mutate(region = ifelse(dset %in%c(206),"FG_nour", region))%>%   #Nouragues2009	French Guiana
  mutate(region = ifelse(dset %in%c(211),"FG_sin", region))%>%   #Sinnamary2011	French Guiana
  mutate(region = ifelse(dset %in%c(216),"FG_ps", region))%>%   #PetitSaut2014	French Guiana
  mutate(region = ifelse(dset %in%c(221),"PR", region))%>%   #ElVerde1996	Puerto Rico
  mutate(region = ifelse(dset %in%c(231),"PR", region))   #Sonadora2004	Puerto Rico


#BWGdb 
# We need to define groups for the traits. Within each trait there are a x number of groups
groups<-c(4,5,2,8,6,7,8,4,5,8,3,4)
names(groups)=c("AS","BS","DM","FD","FG","LO","RE","RF","RM","MD","CP","BF")##noms des traits

reg1<-as.numeric(as.factor(sort(unique(fdat2$region))))
names(reg1)<-sort(unique(fdat2$region))

#attention: no all-zero columns/rows, must delete!!
which(colSums(fdat2[,3:(dim(fdat2)[2]-1)])==0)
which(rowSums(fdat2[,3:(dim(fdat2)[2]-1)])==0) 

postscript("plot.ft3.ps",horizontal=TRUE)
par(mfrow=c(1,1),omi=c(1,1,1,1),mar=c(5,6,3,1)*.5,las=1)
library(RColorBrewer)
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
  
  
  
  
