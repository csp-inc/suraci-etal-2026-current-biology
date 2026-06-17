## ---------------------------
##
## Script name: assess-covariate-soe.R
##
## Purpose of script: prep detection rate data and
## camera covariates and run scale of effect analysis
##
## Author: Justin Suraci
##
## Date last updated: 4/9/2025
##
## ---------------------------

library(tidyverse)
library(sf)
library(glmmTMB)
library(reshape2)
library(lubridate)
library(hms)
library(MuMIn)
library(corrplot)
source("ct-analysis/image-prep-utils.R")

#######################################
# STEP 1 - PREP DATA
#######################################
# Load full detection dataset
load("data/op-ct-all-detections-2025-02-14.rda", verbose = T)
indDet <- det30min$indDet

# Load camera data
load("data/all_projects_covariates_20250814.rda")

##################
##################
# Truncate detection dataset to each camera's study-specific activity period
detTrunc <- data.frame()
for(i in 1:nrow(cams)){
  c_temp <- cams[i,]
  d_temp <- indDet %>% filter(CamID == c_temp$CamID,
                              Date >= c_temp$Station_setup & Date <= c_temp$Station_takedown) %>% 
    arrange(timeStamp)
  detTrunc <- rbind(detTrunc, d_temp)
}

# Summarize total detections of relevant species
sp_counts <- detTrunc %>% group_by(CamID) %>% 
  summarize(deer_count = length(which(Species == 'Deer')),
            elk_count = length(which(Species == 'Elk')),
            mbeaver_count = length(which(Species == 'MountainBeaver')),
            puma_count = length(which(Species == 'Puma')),
            human_count = length(which(Species == 'Human')))

# Join species counts to camera dataset and calculate detection rate
cams <- cams %>% left_join(sp_counts, by = 'CamID')
cams <- cams %>% mutate(deer_rate = deer_count / (count/100),
                        elk_rate = elk_count / (count/100),
                        mbeaver_rate = mbeaver_count / (count/100),
                        puma_rate = puma_count / (count/100),
                        human_rate = human_count / (count/100))

#######################################
# STEP 2 - RUN SINGLE COV MODELS
#######################################

# Compile list of cov pairs to compare and run comparison
cov_list = list(c("elevation_100m","elevation_1km"),
                c("slope_100m","slope_1km"),
                c("aspect_100m","aspect_1km"),
                c("vrm_100m","vrm_1km"),
                c("edge_pcov_100m","edge_pcov_1km"),
                c("tree_pcov_100m","tree_pcov_1km"),
                c("shrub_pcov_100m","shrub_pcov_1km"),
                c("grass_pcov_100m","grass_pcov_1km"),
                c("impervious_100m","impervious_1km"),
                c("all_ag_pcov_100m","all_ag_pcov_1km"),
                c("developed_pcov_100m","developed_pcov_1km"),
                c("human_pop_den_100m","human_pop_den_1km"),
                c("nightlight_100m","nightlight_1km"),
                c("dist_edge_100m", "dist_edge_1km"),
                c("killsite_mod_100m","killsite_mod_1km"),
                c("movement_mod_100m","movement_mod_1km"),
                c("regen_dist_5yrs_pcov_1km","regen_dist_5yrs_pcov_100m"),
                c("regen_dist_15yrs_pcov_1km","regen_dist_15yrs_pcov_100m"))

# Function to compare single cov models
cmFit <- function(covs, data, species){
  cov1 <- data %>% pull(all_of(covs[1]))
  cov2 <- data %>% pull(all_of(covs[2]))
  sp_dat <- data %>% pull(paste0(species, "_count"))
  
  print(paste0("running ", covs[1]))
  mod1 <- glmmTMB(sp_dat ~ cov1 + offset(log(count)), data = data, family = poisson, ziformula = ~ 1)
  print(paste0("running ", covs[2]))
  mod2 <- glmmTMB(sp_dat ~ cov2 + offset(log(count)), data = data, family = poisson, ziformula = ~ 1)
  
  # collect AIC scores
  AIC1 <- AICc(mod1)
  AIC2 <- AICc(mod2)
  tab_temp <- data.frame(covs = covs, AIC = c(AIC1, AIC2))
  if(is.na(AIC1) & is.na(AIC2)==FALSE){
    keep_temp <- covs[2]
  } else if(is.na(AIC1) == FALSE & is.na(AIC2)){
    keep_temp <- covs[1]
  } else if(is.na(AIC1) & is.na(AIC2)){
    keep_temp <- "fit failed"
  } else if(AIC1 < AIC2){
    keep_temp <- covs[1]
  } else keep_temp <- covs[2]
  
  # keep_temp <- if(AIC1 < AIC2) covs[1] else covs[2]
  
  # Report results
  print(paste(covs[1], "AIC:", AIC1, sep = ""))
  print(paste(covs[2],"AIC:", AIC2, sep = ""))
  
  temp_list = list(tab_temp, keep_temp)
  return(temp_list)
  gc()
}

# Fit list of mod pairs and compile output
# Run for deer
mod_map_deer <- cov_list %>% purrr::map(cmFit, data = cams, species = 'deer')
cov_keep_deer <- mod_map_deer %>% purrr::map(2) %>% unlist()
cov_tab_deer <- mod_map_deer %>% purrr::map(1) %>% do.call(rbind, .)
write_csv(cov_tab_deer, "output/occ-model/cov-tab-deer-20250409.csv")

# Run for pumas
mod_map_puma <- cov_list %>% purrr::map(cmFit, data = cams, species = 'puma')
cov_keep_puma <- mod_map_puma %>% purrr::map(2) %>% unlist()
cov_tab_puma <- mod_map_puma %>% purrr::map(1) %>% do.call(rbind, .)
write_csv(cov_tab_puma, "output/occ-model/cov-tab-puma-20250409.csv")

# Run for elk
mod_map_elk <- cov_list %>% purrr::map(cmFit, data = cams, species = 'elk')
cov_keep_elk <- mod_map_elk %>% purrr::map(2) %>% unlist()
cov_tab_elk <- mod_map_elk %>% purrr::map(1) %>% do.call(rbind, .)
write_csv(cov_tab_elk, "output/occ-model/cov-tab-elk-20250409.csv")

# Run for mountain beaver
mod_map_mbeaver <- cov_list %>% purrr::map(cmFit, data = cams, species = 'mbeaver')
cov_keep_mbeaver <- mod_map_mbeaver %>% purrr::map(2) %>% unlist()
cov_tab_mbeaver <- mod_map_mbeaver %>% purrr::map(1) %>% do.call(rbind, .)
write_csv(cov_tab_mbeaver, "output/occ-model/cov-tab-mbeaver-20250409.csv")

#######################################
# STEP 3 - TEST FOR COV CORRELATIONS
#######################################
# Test for correlations
cam_sub<-cams %>% as.data.frame() %>% 
  select(c(cov_keep_deer, "road_dist", "riparian_dist", "developed_dist","ag_dist","FPA_TH_activity")) 
camCor <- cor(cam_sub[complete.cases(cam_sub),])
diag(camCor) <- NA
camCor %>% apply(2, function(x){ifelse(abs(x)<0.6, 0, x)})
corrplot(abs(camCor), method = "color", col = colorRampPalette(c("blue", "white", "red"))(200), 
         is.corr = FALSE, tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.8)

