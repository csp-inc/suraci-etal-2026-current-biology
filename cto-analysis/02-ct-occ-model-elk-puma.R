## ---------------------------
##
## Script name: ct-occ-model.R
##
## Purpose of script: this script formats input datasets for our continuous-time  
##                    occupancy model and then fits the model
##
## Author: Mae Lacey (adapted from example script shared/developed by Arielle Parsons)
## All C++ code for fitting CTO models was based on: 
## Kelner et al. 2021 https://link.springer.com/article/10.1007/S13253-021-00482-Y
##
## Date last updated: 9/8/2025
##
## Email contact: mae[at]csp-inc.org
##
## ---------------------------

library(dplyr)
library(sf)
library(parallel)
library(Rcpp)
library(stringr)
library(tidyverse)
library(reshape2)
library(lubridate)
library(hms)



# ------------------------------------------------------------------------------
# Load relevant datasets and functions
# ------------------------------------------------------------------------------

# load indDetection function
source("cto-analysis/image-prep-utils.R")

'%!in%' <- function(x,y)!('%in%'(x,y))

# load detections
load("data/op-ct-all-detections-2025-02-14-final.rda")

# load covariates
load("data/all_projects_covariates_20250814.rda") 
cams <- elk_cams # renaming to cams


# ------------------------------------------------------------------------------
# User inputs
# ------------------------------------------------------------------------------

# set time zone for datasets loaded above
timezone = "America/Los_Angeles"

# set det based on time thinning of choice
det_select <- det60min # using 60 min throughout moving forward
det <- det_select$indDet[order(det_select$indDet$Date, decreasing=FALSE),]

# set some parameters for naming
prey = "elk"
pred = "puma"
#modtype = "RD2DM2"
modtype = "null"


# ------------------------------------------------------------------------------
# Data prep
# ------------------------------------------------------------------------------

# ------------------------------- Detection data ------------------------------- 

# pulling middle 90 days' worth of data
cams$days_elapsed <- difftime(cams$Station_takedown, cams$Station_setup,units = "days") 
cams$days_mid_deploy <- round(cams$days_elapsed/2, digits = 0)
cams$mid_deploy_start <- cams$days_mid_deploy - 45 # use 15 for the middle 30 days, 45 for the middle 90 days
cams$mid_deploy_end <- cams$days_mid_deploy + 45 # use 15 for the middle 30 days, 45 for the middle 90 days
cams <- cams %>%
  mutate(alt_deploy_startdate = case_when(days_elapsed <= 90 ~ Station_setup, # use 30 for 30 days, 90 for 90 days
                                          days_elapsed > 90 ~ as.Date(Station_setup) + mid_deploy_start),
         alt_deploy_enddate = case_when(days_elapsed <= 90 ~ Station_takedown, 
                                        days_elapsed > 90 ~ as.Date(Station_setup) + mid_deploy_end)) %>%
  dplyr::select(-c(Station_setup, Station_takedown)) # removing original dates to rename alt date columns
# then rename
names(cams)[names(cams) == "alt_deploy_startdate"] <- "Station_setup"
names(cams)[names(cams) == "alt_deploy_enddate"] <- "Station_takedown"

# we're setting covs as the cams dataset because that contains all site-specific covariates
# currently contains raw covariates, but I'm scaling the data down below
covs <- cams

# filter detection dataset to match deployment timeframe from cams
det$Date <- as.Date(det$Date , format = "%m/%d/%y") # make sure date field is formatted properly
deploy_dates <- cams %>%
  dplyr::select(CamID, Station_setup, Station_takedown)
det_deploy <- merge(det, deploy_dates, by = "CamID", all = FALSE) # Note that some cameras fall out
# quick check to make sure we're filtering based on elk_cams and not retaining any other cameras
u1 <- unique(det_deploy$CamID)
u2 <- unique(cams$CamID)
setdiff(u1, u2) # values in det_deploy but not in cams
setdiff(u2, u1) # values in cams by not in det_deploy

# then actually filter dets based on deployment dates
det_deploy <- det_deploy %>% filter(Date >= Station_setup & Date <= Station_takedown)

# pull lat and long from cam dataset, add to detection dataset
cam_latlon <- as.data.frame(unique(covs[c("Project", "Station", "CamID", "year", "X", "Y")])) %>%
  dplyr::select(CamID, X, Y)
det <- left_join(det_deploy, cam_latlon, by=c("CamID"))

# check count of detections by species
table(det$Species)

# you need to do this again to go back from the numeric timestamp
det$Date=as.POSIXct(det$Date,format="%Y-%m-%d")
det$Date=lubridate::with_tz(det$Date, tzone = timezone)
det$Date=lubridate::date(det$Date)
time=det$Time
time=as.character(time)
time=gsub("[.]", ":", time)
dt=paste(det$Date, time, sep = " ")
det$timestamp=as.POSIXct(det$timeStamp, format="%Y-%m-%d %H:%M:%S",origin="01-01-1900", tz=timezone)

# removing those deployments with < 24 hours of data
det$Begin<-as.POSIXct(det$timestamp, format="%Y-%m-%d %H:%M:%S", origin="01-01-1900",
                       tz=timezone)

duration<-det%>%group_by(CamID)%>%summarize(duration=max(as.Date(Begin,format="%Y-%m-%d"))-min(as.Date(Begin,format="%Y-%m-%d")))

l<-duration[which(duration$duration<1),]

which(det$CamID %in% l$CamID)

dets<-det[which(det$CamID %!in% l$CamID),]

# check count of detections by species and project - final sample sizes
# then filter to only retain elk and puma and then get project n
dets_samp <- dets %>% filter(Species == "Elk" | Species == "Puma")
#dets_samp <- dets %>% filter(Species %in% c("Elk", "Puma")) # alt method
table(dets_samp$Species)
table(dets_samp$Project)
ggplot(as.data.frame(table(dets_samp$Project)), aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.7) +
  labs(title = "Sample Size by Project",
       x = "Project",
       y = "Count of Detections") +
  theme_classic() +
  scale_fill_brewer(palette = "Set2") 
#ggsave(paste0('summary-plots/sample-size-by-proj-', time_filt, '_filt.png'), height=4, width=6, dpi=300)


# unique deployments
deps <- unique(dets$CamID)

# date / time of each detection
dt <- as.POSIXlt(dets$Begin, tz = timezone, format = "%Y-%m-%d %H:%M:%S")

#Sort the detections

# blank lists for storing detection times
sp1 <- sp2 <- sp1_ind <- sp2_ind <- vector('list', length(deps))

# deployment length (in days)
dep_len <- numeric(length(deps))

dep_start <- as.POSIXct(rep(NA, length(deps)),
                        tz = timezone)

# identifying start and end time of each deployment
for(i in 1:length(deps)){
  
  # index of deployment i
  ind <- which(dets$CamID == deps[i])
  
  dep_start[i] <- min(dt[ind])
  
  # deployment length
  dep_len[i] <- as.numeric(difftime(max(dt[ind]), min(dt[ind]),
                                    units = 'days'))
  
  # any sp1 detected?
  if(any(dets$Species[ind] == 'Elk')){ # sp1 is the subordinate species
    
    # indices of detections
    sp1_ind[[i]] <- which(dets$Species[ind] == 'Elk')
    
    # time of sp1 detection, relative to camera setup
    sp1[[i]] <- as.numeric(difftime(dt[ind][sp1_ind[[i]]], min(dt[ind]),
                                    units = 'days'))
    
  }
  
  # any sp2 detected?
  if(any(dets$Species[ind] == 'Puma')){ # sp2 is the dominant species
    
    # indices of detections
    sp2_ind[[i]] <- which(dets$Species[ind] == 'Puma')
    
    # time of sp2 detection, relative to camera setup
    sp2[[i]] <- as.numeric(difftime(dt[ind][sp2_ind[[i]]], min(dt[ind]),
                                    units = 'days'))
    
  }
}

source('cto-utils/splt_lik.R')  # calculated likely of detection times
cl <- makeCluster(10)  # sending to 20 cores
clusterEvalQ(cl, source('cto-utils/splt_lik.R')) # loading functions
clusterExport(cl, c('sp1', 'sp2', 'dep_len')) # sending variables
sourceCpp("cto-utils/mmpp_covs.cpp")

# Some additional data required for cpp version
# this is also calculated for R version, 
# but inside the actual lik function
y1_i <- sapply(sp1, function(x) ifelse(is.null(x), 0, 1))
y2_i <- sapply(sp2, function(x) ifelse(is.null(x), 0, 1))

# Get time difference between detection times and interval boundaries
#   y = detection times at a camera site
#   J = total time length of survey at a camera site
#   inc = amount of time in each interval (in units of days)
get_yd <- function(y, J, inc=1){
  d <- seq(0, J, by=inc)
  if((J-d[length(d)]) > 0){
    d <- c(d, J)
  }
  
  if(is.null(y)){
    groups <- lapply(1:(length(d)-1), function(x) numeric(0))
  } else{
    print(d)
    groups <- split(y, cut(y, d))
  }
  
  groups2 <- lapply(1:length(groups), function(i){
    c(d[i], groups[[i]], d[i+1])
  })
  
  out <- lapply(groups2, function(x) diff(x))
  
}

# Get time differences yd for each species at each site
yd_sp1 <- lapply(1:length(sp1), function(i) get_yd(sp1[[i]], 
                                                   dep_len[i], inc=1/24))
yd_sp2 <- lapply(1:length(sp2), function(i) get_yd(sp2[[i]], 
                                                   dep_len[i], inc=1/24))

yd1 <- unlist(yd_sp1)
yd2 <- unlist(yd_sp2)

# --------------------------------- Covariates --------------------------------- 
#Match deployment id to site id
#covs$CamID<-covs$site
dep_to_site <- dets %>% group_by(CamID) %>% summarize()

site_covs <- data.frame(dep_to_site) %>%
  left_join(covs) %>%
  # rename and scale all relevant covs
  mutate(
    slope = scale(slope_1km),
    edge = scale(edge_pcov_1km),
    shrub = scale(shrub_pcov_100m),
    riparianDist = scale(riparian_dist),
    agDist = scale(ag_dist),
    devDist = scale(developed_dist),
    roadDist = scale(road_dist),
    nightlight = scale(nightlight_100m),
    regen = scale(regen_dist_15yrs_pcov_100m),
    killMod = scale(killsite_mod_1km),
    roadDens = scale(road_density),
    roadDistAll = scale(road_dist_all)
  ) %>% 
  dplyr::select(CamID, slope, edge, shrub, riparianDist, agDist, devDist, 
                roadDist, nightlight, regen, killMod, roadDens, roadDistAll) %>%
  as_tibble()

t_covs <- data.frame(dep_to_site) %>%
  left_join(covs) %>%
  # rename and scale all relevant covs
  mutate(
    dev = scale(developed_pcov_1km),
    roadDist = scale(road_dist),
    nightlight = scale(nightlight_100m)
    #killMod = scale(killsite_mod_1km)
  ) %>% 
  dplyr::select(CamID, dev, roadDist, nightlight #killMod
                ) %>%
  as_tibble()
#write.csv(t_covs, "t_covs_ct_occ_model.csv")
#saveRDS(t_covs, "t_covs_ct_occ_model.Rds")

# Make observation covariate (time of day)
# Make sequence of times for each deployment
names(yd_sp1) <- deps
ndet <- sapply(yd_sp1, length)

sec_in_inc <- 60*60 #seconds in each increment (1 hr)

time_list <- lapply(1:length(ndet), function(i){
  # Check if dep_start[i] is finite
  if(!is.finite(dep_start[i])) {
    return(NULL)  # Skip or return NULL for non-finite dep_start[i]
  }
  
  tseq <- as.POSIXlt(seq(dep_start[i], by=sec_in_inc, length.out=ndet[i]))
  tseq$hour + tseq$min/60 + tseq$sec/3600
})
time_vec <- unlist(time_list)

# Vectorize the time since the last detection of the dominant species
# on the same scale as time_list
time_list2 <- lapply(1:length(ndet), function(i){
  seq(0, dep_len[[i]], by=1/24)
})
time_vec2 <- unlist(time_list2)
length(time_vec)==length(time_vec2)

# other covs - Fourier series
obs_covs <- data.frame(deploy = rep(deps, ndet),
                       f1c = cos(pi*time_vec/12),
                       f2c = cos(2*pi*time_vec/12),
                       f1s = sin(pi*time_vec/12),
                       f2s = sin(2*pi*time_vec/12)
                       )
# bringing in site-specific covs that are going into the detection model
# now we need to join our t_covs to obs_covs to repeat the same value for each cam across all records (because these are annual covs, not hourly or daily) 
obs_covs <- merge(x = obs_covs, y = t_covs, by.x = "deploy", by.y = "CamID")

# Index to subset lambda values by site
lidx_i <- matrix(NA, nrow=length(yd_sp1), ncol=2)
idx <- 0
for (i in seq_along(yd_sp1)){
  lidx_i[i,1] <- idx
  lidx_i[i,2] <- idx + length(yd_sp1[[i]]) - 1
  idx <- idx + length(yd_sp1[[i]])
}

# Index to subset yd (y-d) values by site i and interval j
# yd is now a vector instead of a list of lists so the index is needed
maxj <- max(sapply(yd_sp1, length))
yd1_st_idx <- yd1_en_idx <- matrix(NA, nrow=length(yd_sp1), ncol=maxj)
idx <- 0
for (i in seq_along(yd_sp1)){
  yd_sub <- yd_sp1[[i]]
  for (j in seq_along(yd_sub)){
    yd1_st_idx[i,j] <- idx
    yd1_en_idx[i,j] <- idx + length(yd_sub[[j]]) - 1
    idx <- idx + length(yd_sub[[j]])
  }
}

yd2_st_idx <- yd2_en_idx <- matrix(NA, nrow=length(yd_sp2), ncol=maxj)
idx <- 0
for (i in seq_along(yd_sp2)){
  yd_sub <- yd_sp2[[i]]
  for (j in seq_along(yd_sub)){
    yd2_st_idx[i,j] <- idx
    yd2_en_idx[i,j] <- idx + length(yd_sub[[j]]) - 1
    idx <- idx + length(yd_sub[[j]])
  }
}


# ------------------------------------------------------------------------------
# Fit model
# ------------------------------------------------------------------------------
# first set which model we're using based on the model type specified above
if (modtype == "RD1DM2") {
  print("running road model 1 x detection model 2")
  # ----- STATE MODEL -----
  X_f1 <- X_f12 <- X_f2 <- model.matrix(~edge + slope + riparianDist + devDist + I(devDist^2) + roadDist + I(roadDist^2) + roadDens, site_covs)
  # ----- DETECTION MODEL -----
  X_lam1 <- X_lam2 <- X_lam3 <- model.matrix(~f1c + f2c + f1s + f2s + dev + roadDist, obs_covs)
} else if (modtype == "null") {
  print("running null model")
  # ----- STATE MODEL -----
  X_f1 <- X_f12 <- X_f2 <- model.matrix(~1, site_covs)
  # ----- DETECTION MODEL -----
  X_lam1 <- X_lam2 <- X_lam3 <- model.matrix(~f1c + f2c + f1s + f2s, obs_covs)
}

# index structure would only need to change here if we added another species, other model changes
pind <- matrix(NA, nrow=8, ncol=2)
pind[1,] <- c(0, 0+ncol(X_f1)-1)
pind[2,] <- c(pind[1,2]+1, pind[1,2]+1+ncol(X_f2)-1)
pind[3,] <- c(pind[2,2]+1, pind[2,2]+1+ncol(X_f12)-1)
pind[4,] <- c(pind[3,2]+1, pind[3,2]+2)
pind[5,] <- c(pind[4,2]+1, pind[4,2]+2)
pind[6,] <- c(pind[5,2]+1, pind[5,2]+1+ncol(X_lam1)-1)
pind[7,] <- c(pind[6,2]+1, pind[6,2]+1+ncol(X_lam2)-1)
pind[8,] <- c(pind[7,2]+1, pind[7,2]+1+ncol(X_lam3)-1)

# Speed may increase with openMP (w/ threads argument)
# depends on real cores not virtual cores - using 10 for now
starts <- optim(rep(0,max(pind)+1), mmpp_covs, method = 'SANN',
                control = list(maxit=1000, trace=1, REPORT =5),
                pind=pind, X_f1=X_f1, X_f2=X_f2, X_f12=X_f12,
                X_lam1=X_lam1, X_lam2=X_lam2, X_lam3=X_lam3, 
                yd1=yd1, yd2=yd2,
                lidx_i=lidx_i, yd1_st_idx=yd1_st_idx, yd1_en_idx=yd1_en_idx, 
                yd2_st_idx=yd2_st_idx, yd2_en_idx=yd2_en_idx, 
                y1_i=y1_i, y2_i=y2_i, threads=10)
starts$par
starts$value

# fit the model
fit <- optim(starts$par, mmpp_covs, method = 'L-BFGS-B', hessian=TRUE,
             control = list(maxit=2000, trace = 1, REPORT = 5),
             pind=pind, X_f1=X_f1, X_f2=X_f2, X_f12=X_f12,
             X_lam1=X_lam1, X_lam2=X_lam2, X_lam3=X_lam3, 
             yd1=yd1, yd2=yd2,
             lidx_i=lidx_i, yd1_st_idx=yd1_st_idx, yd1_en_idx=yd1_en_idx, 
             yd2_st_idx=yd2_st_idx, yd2_en_idx=yd2_en_idx, 
             y1_i=y1_i, y2_i=y2_i, threads=10)


# ------------------------------------------------------------------------------
# Export model outputs
# ------------------------------------------------------------------------------
# change wd for saving all model outputs
saveRDS(fit, paste0(modtype, "-", pred, "-", prey, "-fit.Rds"))

est <- fit$par
names(est) <- c(paste0("f1_",colnames(X_f1)), paste0("f2_",colnames(X_f2)),
                paste0("f12_",colnames(X_f12)),
                "log_mu1[1]","log_mu1[2]","log_mu2[1]","log_mu2[2]",
                paste0("loglam1_",colnames(X_lam1)),
                paste0("loglam2_",colnames(X_lam2)), 
                paste0("loglam3_",colnames(X_lam3)))
se <- sqrt(diag(abs(solve(fit$hessian))))

results <- data.frame(est=round(est, 3), se=round(se,3))
results$lower <- results$est - 1.96*results$se
results$upper <- results$est + 1.96*results$se
results

saveRDS(results, paste0(modtype, "-", pred, "-", prey, "-results.Rds"))
write.csv(results, paste0(modtype, "-", pred, "-", prey, "-results.csv"), row.names = T)

# AIC
(2*length(fit$par))-(2*-fit$value)
