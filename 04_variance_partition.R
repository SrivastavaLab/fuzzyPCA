
res_fpca <- readRDS('PCA_results.rds')

######## Procrustes
#Partition the axis of the bromeliads res_fpca$li by region

vare.proc <-procrustes(res_fpca$li,as.numeric(as.factor(fdat2$region)))
vare.proc
summary(vare.proc)
plot(vare.proc)
plot(vare.proc, kind=2)
residuals(vare.proc)
protest(res_fpca$li,as.numeric(as.factor(fdat2$region)), scores = "sites", permutations = how(nperm = 999))


## How is the clustering related into region
