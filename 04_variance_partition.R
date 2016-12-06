library(vegan)
library(dplyr)

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
plot(betadiver_all)

reg_num_j <- fdat2 %>% select(region) %>% unique() %>% mutate(reg_num = c(14, 5, 6, 4, 3, 2 , 1, 11, 15, 12, 10, 13, 16, 7, 8, 9, 6))

nested_sorted <- left_join(fdat2, reg_num_j) %>% arrange(reg_num) %>% select(-region, -reg_num) %>% nestednodf()
plot(nested_sorted)

part_beta_all <- fdat2 %>% select(-region) %>% nestedbetajac()
part_beta_all

## Turnover across species
library(purrr)

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

