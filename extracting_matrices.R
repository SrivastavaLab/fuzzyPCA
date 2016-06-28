library(fwdata)
library(dplyr)
library(purrr)
library(tidyr)

fwdata:::fw_auth()

full <- fw_data()

full$abundance

full$traits

abundance_dataset_list <- full$abundance %>% 
  mutate(dataset_id = as.factor(dataset_id)) %>% 
  split(.$dataset_id) 

unique_species_id_list <- abundance_dataset_list %>% 
  map(~ .$species_id) %>% 
  map(~ unique(.)) 

trait_by_dataset <- unique_species_id_list %>% 
  map(~ data.frame(species_id = .)) %>% 
  map(~ left_join(., full$traits))

trait_matrix_by_dataset <- trait_by_dataset %>% 
  map(~ select(., species_id, matches("^[A-Z]{2}\\d"))) %>% 
  map(~ as.matrix(.))

abundance_matrix_by_dataset <- abundance_dataset_list %>% 
  map(~ select(., -dataset_id, -bwg_name)) %>% 
  map(~ spread(., brm, abd)) %>% 
  map(~ as.matrix(.))

trait_matrix_by_dataset[[1]][,1] == abundance_matrix_by_dataset[[1]][,1]

