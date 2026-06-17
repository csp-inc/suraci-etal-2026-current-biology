## ---------------------------
##
## Script name: 00-calc-psi-at-cam-locs.R
##
## Purpose of script: Use the top puma-deer
## CTO model and camera level covariates
## to calculate occupancy probability at
## each camera site. Calculations based on
## Rota et al. 2016
##
## Author: Justin Suraci
##
## ---------------------------
library(tidyverse)
library(sf)

# Load Camera Data
load("data/all_projects_covariates_20250814.rda", verbose = TRUE)

# Load model output
# USING PUMA-DEER CTO SM4DM2
fit <- readRDS("RD1DM2-puma-deer-fit.Rds")
results <- readRDS("RD1DM2-puma-deer-results.Rds")

# Prep covariates that went into the model
edge <- scale(cams$edge_pcov_1km)[,1]
slope <- scale(cams$slope_1km)[,1]
riparianDist <- scale(cams$riparian_dist)[,1]
devDist <- scale(cams$developed_dist)[,1]
roadDist <- scale(cams$road_dist)[,1]
roadDens <- scale(cams$road_density)[,1]
predMatrix <- data.frame(int = 1,
                         edge = edge, 
                         slope = slope, 
                         riparianDist = riparianDist,
                         devDist = devDist,
                         devDis_sq = devDist^2,
                         roadDist = roadDist,
                         roadDis_sq = roadDist^2,
                         roadDens = roadDens) %>% 
  as.matrix()

#Get the beta means and 95%CIs
beta <- fit$par
sigma <- solve(fit$hessian)

#Index the appropriate betas for the f params
f1_ind<-1:9
f2_ind<-10:18
f12_ind<-19:27

f1_b <- beta[f1_ind]
f2_b <- beta[f2_ind]
f12_b <- beta[f12_ind]

# Predict from model for each camera location
f1_pred <- sweep(predMatrix, 2, f1_b, "*") %>% 
  apply(1, sum)
f2_pred <- sweep(predMatrix, 2, f2_b, "*") %>% 
  apply(1, sum)
f12_pred <- sweep(predMatrix, 2, f12_b, "*") %>% 
  apply(1, sum)

# Calculate Psi based on Equation 2 in Rota et al. 2016 MEE
psi11 = exp(f1_pred + f2_pred + f12_pred) / 
  (1 + exp(f1_pred) + exp(f2_pred) + exp(f1_pred + f2_pred + f12_pred))
psi10 = exp(f1_pred) / 
  (1 + exp(f1_pred) + exp(f2_pred) + exp(f1_pred + f2_pred + f12_pred))
psi01 = exp(f2_pred) / 
  (1 + exp(f1_pred) + exp(f2_pred) + exp(f1_pred + f2_pred + f12_pred))
psi00 = 1 / (1 + exp(f1_pred) + exp(f2_pred) + exp(f1_pred + f2_pred + f12_pred))

# Calculate occupancy probability for deer and pumas
# based on equation 2 in Rota et al. 2016 MEE
deerOccProb = psi11 + psi10
pumaOccProb = psi11 + psi01

# Join back to camera dataset
cams$deerOccProb = deerOccProb
cams$pumaOccProb = pumaOccProb
cams$probUnoccupied = psi00
cams$deerOnlyProb = psi10
cams$pumaOnlyProb = psi01
cams$deerPumaProb = psi11

# Save it
save(cams, file = "data/cams-w-puma-deer-occ.rda")
st_write(cams, "data/cams-w-puma-deer-occ.gpkg",
         driver = "GPKG")
