library(vegan)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

res_fpca <- readRDS('Data/PCA_results.rds')

fdat2 <- read.csv('Data/Abundance_weighted_traits.csv')


######## Procrustes
#Partition the axis of the bromeliads res_fpca$li by region

vare.proc <-procrustes(res_fpca$li,as.numeric(as.factor(fdat2$region)))
vare.proc
summary(vare.proc)
plot(vare.proc)
plot(vare.proc, kind=2)
residuals(vare.proc)
protest(res_fpca$li,as.numeric(as.factor(fdat2$region)), scores = "sites", permutations = how(nperm = 999))


## Turnover accross traits

betadiver_all <- fdat2 %>% select(-region) %>% betadiver()
png("Figures/BetaDiv/BetaAll.png")
plot(betadiver_all)
dev.off()

reg_num_j <- fdat2 %>% select(region) %>% unique() %>% mutate(reg_num = c(14, 5, 6, 4, 3, 2 , 1, 11, 15, 12, 10, 13, 16, 7, 8, 9, 6))

nested_sorted <- left_join(fdat2, reg_num_j) %>% arrange(reg_num) %>% select(-region, -reg_num) %>% nestednodf()
png("Figures/BetaDiv/NestedAll.png")
plot(nested_sorted)
dev.off()

part_beta_all <- fdat2 %>% select(-region) %>% nestedbetajac()
part_beta_all


##  Trait turnover based on region 

fdat_region <- split(fdat2, fdat2$region) %>% 
  map(~ select(.x, -region))

betadiver_region <- fdat_region %>% 
  map(~ betadiver(.x)) 

for(i in 1:length(betadiver_region)){
  tit <- betadiver_region[i] %>% names()
  png(paste0('Figures/BetaDiv/',tit, 'betadiv.png'))
  plot(betadiver_region[[i]], main = tit)
  dev.off()
}

for(i in 1:length(betadiver_region)){
  region_name <-betadiver_region[i] %>% names()
  nested <- fdat2 %>% filter(region == region_name) %>% select(-region) %>% nestednodf()
  png(paste0('Figures/BetaDiv/',region_name, 'nested.png'))
  plot(nested)
  dev.off()
}

beta_part <- list()

for(i in 1:length(betadiver_region)){
  region_name <-betadiver_region[i] %>% names()
  part <- fdat2 %>% filter(region == region_name) %>% select(-region) %>% nestedbetajac()
  beta_part[[i]] <- data.frame(region_name, part)
}

Beta_part_table <- beta_part %>% 
  map(~ rownames_to_column(.x)) %>% 
  map(~ spread(.x, key = rowname, value = part)) %>% 
  bind_rows() %>% 
  arrange(jaccard)

## Beta disper
betadisper_fdat <- betadisper(dist(fdat2 %>% select(-region)), fdat2$region, type = 'centroid')
permutest(betadisper_fdat, pairwise = TRUE, permutations = 99)

anova(betadisper_fdat)

(mod.HSD <- TukeyHSD(betadisper_fdat))
plot(mod.HSD)


## Change plot axis <- not working 
plot(betadisper_fdat, axes = c(3, 2), ellipe = TRUE)
jpeg('Figures/BetaDiv/BetaDisper.jpeg')
boxplot(betadisper_fdat)
dev.off()










## Turnover across species

abundance <- abundance_matrix_by_dataset %>% 
  map(~ as.data.frame(.x)) %>% bind_rows(.id = 'dset')

abundance[is.na(abundance)] <- 0

abundance <- abundance %>% mutate(dset = as.numeric(dset))

regions <- read.csv('Data/Regions.csv', header = TRUE)
abun2 <- left_join(abundance,regions) %>% select(-dataset)


betadiver_all_sp <- abun2 %>% select(-region) %>% betadiver()
plot(betadiver_all_sp)

reg_num_j <- abun2 %>% select(region) %>% unique() %>% mutate(reg_num = c(14, 5, 6, 4, 3, 2 , 1, 11, 15, 12, 10, 13, 16, 7, 8, 9, 6))

nested_sorted_sp <- left_join(abun2, reg_num_j) %>% arrange(reg_num) %>% select(-region, -reg_num) %>% nestednodf()
plot(nested_sorted_sp)

part_beta_all <- abun2 %>% select(-region) %>% nestedbetajac()
part_beta_all

