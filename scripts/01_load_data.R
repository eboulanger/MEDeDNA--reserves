#'##############################################################################
#' 
#'  codes to load and prepare data to reproduce analyses and figures of the Boulanger et al. (under review) paper: 
#'  "Environmental DNA metabarcoding reveals and unpacks a biodiversity conservation paradox in Mediterranean marine reserves." 
#'  
#'##############################################################################

#### load all packages
# for wrangling
library(dplyr)
library(forcats)
library(stringr)
# for (a) 
library(rsq)
library(margins)
# for (b)
library(betapart)
library(reshape)
library(tidyselect)
library(vegan)
# for (c)
# for figures
library(ggplot2)
library(patchwork)
library(ggalt)
library(ggrepel)
library(grid)

#### load and wrangle data 

# 1. species detection ----
species_ab <- read.csv(here::here("data", "MEDeDNA_species_ab.csv")) 
# keep only species names as columns and replicates as rows, and make presence-absence
species <- species_ab %>% 
  `rownames<-`(.$Species) %>% 
  select(-Class, -Order, -Family, -Genus, -Species) %>% 
  t(.)
species[species>0] <- 1

# 2. sample metadata and explanatory variables ----
meta_env <- read.csv(here::here("data", "MEDeDNA_sample_metadata.csv")) %>% 
  mutate(Replicate = Code) %>% 
  mutate(Site = str_sub(Replicate, end = -3)) %>% 
  `rownames<-`(.$Code)

#    make sure species and site data are in the same order
rownames(species) == rownames(meta_env)
meta_env <- meta_env[rownames(species), ]

#    set the order of protection and region factor levels
meta_env <- meta_env %>% 
  mutate(Protection = fct_relevel(Protection, c("reserve", "outside5", "outside10"))) %>%   # set the correct protection level order
  mutate(Region = fct_inorder(Region)) # reorder factor levels in same order as regions appear in dataset

# 3. species traits ----
sp_traits <- read.csv(here::here("data","MEDeDNA_species_traits.csv")) %>% 
  `rownames<-`(.$Species)

# 4. extract species lists by trait ----

 # cryptobenthic species
crypto_sp <- species_ab %>% 
  filter(Family %in% c("Tripterygiidae","Grammatidae","Creediidae","Aploactinidae","Gobiidae","Chaenopsidae","Gobiesocidae",
                       "Labrisomidae","Pseudochromidae","Bythitidae","Plesiopidae","Dactyloscopidae","Blenniidae","Apogonidae",
                       "Callionymidae","Opistognathidae","Syngnathidae")) %>% # list of cryptobenthic families by Brandl et al 2018 Biological Reviews
  pull(Species) %>% as.character()

 # pelagic species
pelagic_sp <- sp_traits %>% 
  filter(Vertical_Distribution_adj == "Pelagic") %>% 
  pull(Species) %>% as.character()

# demersal species
demersal_sp <- sp_traits %>% 
  filter(Vertical_Distribution_adj == "Demersal") %>% 
  pull(Species) %>% as.character()

# benthic species
benthic_sp <- sp_traits %>% 
  filter(Vertical_Distribution_adj == "Benthic") %>% 
  pull(Species) %>% as.character()

 # detected species rarity
 # for each region, get species detected in only 1 or 2 samples
regions <- c("BAN", "CR", "CAL", "PQ", "CX", "CV")
rarity <- list()
for (i in regions) {
  subset <- select(as.data.frame(t(species)), starts_with(i)) # subset each region
  subset$rarity <- apply(subset, 1, sum)            # add rowsum to infer rarity
  rarity[[i]] <- rownames(subset)[subset$rarity > 0 & subset$rarity < 3] # extract species with rep-level occurrence < 3
}
rare_sp <- unlist(rarity, use.names = FALSE) %>% unique()

 # vulnerable species
vulnerable_sp <- sp_traits %>% 
  filter(Vulnerability > 70) %>% 
  pull(Species) %>% as.character()

