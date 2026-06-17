## ---------------------------
##
## Script name: 03-wvc-vs-puma-kde.R
##
## Purpose of script: Corroborating analysis
## for effect of pumas on DVC count
## using data from collared pumas
##
## Author: Justin Suraci
##
## ---------------------------

library(tidyverse)
library(sf)
library(terra)
library(amt)
library(patchwork)
library(sfheaders)
library(data.table)
library(tictoc)
library(spdep)
library(adespatial)
library(glmmTMB)
library(nnet)
library(effects)
library(DHARMa)
library(reshape2)
library(MuMIn)

# ---------- COLLECT DATASETS -------------------
# Get AOI
puma_crs <- st_crs(26910)
aoi <- st_read("data/puma-refined-aoi.geojson") %>% st_transform(puma_crs)

# Get WVC datasets
load(file = "data/deer-vehicle-collision-data/wvc-analysis-dataset-20251109.rda", verbose = T)

# Get puma location data
load(file = 'data/puma-locs-full-age-class-20230321.rda', verbose = TRUE)
#subset to AOI
puma_sp <- st_as_sf(puma_ac, coords = c("UTM_X", "UTM_Y"), crs = puma_crs)
puma_in <- st_within(puma_sp, aoi, sparse = F) %>% as.vector()
puma_ac <- puma_ac[puma_in,]
# identify individuals with relocation time <= 2 hrs
reloc <- puma_ac %>% group_by(name) %>% 
  summarize(hr_median = (median(dt_hr))) %>%
  filter(hr_median <=2.5) %>% as.data.frame()

# Cut-off date
co_date <- as.POSIXct("2020-01-01", tz = "America/Los_Angeles")

# Process puma location dataset
puma_locs <- puma_ac %>%
  # Remove dependent young
  filter(age_class != "dep_young",
         name %in% reloc$name) %>% 
  # Create age + year + id column to break tracks up by age class and year (for yearly covs)
  mutate(yr = year(date_time),
         # age_ID = paste(name, age_class, yr, sep = "_"),
         age_ID = paste(name, age_class, sep = "_")) %>% 
  # Remove earlier tracks
  filter(date_time >= co_date)

# ---------- PREP KDE -------------------
### Get unique puma ids and fit a 95% kernel density estimated 'home range' to the data 
### Test of this with those pumas with at least a week of steps 
## Get list of unique ids 
pID <- unique(puma_locs$age_ID)
## Set up list to hold KDEs
pumaKDEs <- list()
# Template raster
slopeRast <- rast("data/ssf-cov-rasters/slope-90m.tif")

for(i in 1:length(pID)){
  # track progress
  print(paste("Running", pID[i], "......", i, "of", length(pID)))
  
  # Skip small # track
  n <- puma_locs %>% 
    dplyr::filter(puma_locs$age_ID == pID[i]) %>% nrow()
  print(n)
  if(n < 100) next

  # Make KDE
  kdeTemp <- puma_locs %>%
    dplyr::filter(age_ID == pID[i]) %>%
    make_track(.x = UTM_X, .y= UTM_Y, .t= date_time, crs=puma_crs) %>%
    # hr_kde(h = hr_kde_lscv(.)$h,
    #        trast = slopeRast,
    #        levels = 0.95)
    hr_kde(trast = slopeRast)
  outRast <- kdeTemp$ud
  names(outRast) <- pID[i]

  # Add track length as metadata tag
  mtag <- paste0("track_length=",n)
  metags(outRast) <- mtag

  # Save individual rasters
  outName = paste0("output/puma-kdes/kde-v2/",pID[i],"-kde-temp.tif")
  writeRaster(outRast, outName, overwrite = TRUE)
}

### Weight the rasters by the total number of steps in the trajectory
### Sum the weighted rasters to produce final output
# Get individual kde file names
kdeFiles <- list.files("output/puma-kdes/kde-v2")
# Make an all zero raster to add to
sumWeighted <- slopeRast
values(sumWeighted) <- 0
for(i in 1:length(kdeFiles)){
  # Read in individual KDE file and extract track_length metadata tag as weight
  inFile <- paste0("output/puma-kdes/kde-v2/",kdeFiles[i])
  KDE <- rast(inFile)
  ### Had to adjust this to properly access the metadata tag for track length 
  tags <- metags(KDE)
  weight <- tags$value[tags$name=='track_length'] %>% as.numeric()

  # Create weighted individual raster
  KDE_weighted <- KDE*weight
  
  # add to population raster
  sumWeighted <- sumWeighted + KDE_weighted
}

kde_op <- crop(sumWeighted, vect(aoi))   # clips to bounding box
kde_op <- mask(kde_op, vect(aoi))

#save
writeRaster(kde_op, "output/puma-kdes/weighted-population-kde.tif",
            overwrite = T)


#-----/////-----/////-----/////-----/////-----/////-----
#              GRID ANALYSIS - WVC COUNTS
#-----/////-----/////-----/////-----/////-----/////-----

# ---------- EXTRACT KDE OVER WVC GRID -------------------
kde_op <- rast("output/puma-kdes/weighted-population-kde.tif")

rescale01 <- function(r){
  mm <- minmax(r, compute = TRUE)
  rout <- (r - mm[1])/(mm[2] - mm[1])
  return(rout)
}
# Rescale and reproject kde raster
vecGrid <- vect(wvcGrid)
kde_proj <- kde_op %>% 
  rescale01() %>% 
  project(terra::crs(vecGrid))

# Extract mean pixel value in each wvc grid cell
pumaKDE_mean <- terra::extract(kde_proj, vecGrid, fun = mean, na.rm = TRUE)

# Add to dataset 
wvcGrid$pumaKDE_mean <- pumaKDE_mean$slope

# ---------- RUN ANALYSIS -------------------

# ACCOUNTING FOR SPATIAL STRUCTURE
# Prep spatial eigenvectors
coords <- wvcGrid %>% dplyr::select(x, y) %>% st_drop_geometry()
nb <- dnearneigh(coords, d1 = 0, d2 = 12000)
plot(nb, coords)
lw <- nb2listw(nb, style = "W")
mem <- scores.listw(lw, MEM.autocor = "positive")

# Identify SEV with strongest influence
sevDF <- data.frame()
for(i in 1:ncol(mem)){
  dat <- wvcGrid
  dat$sev <- mem[,i]
  temp <- glmmTMB(wvc_count_deer_5yrs ~ sev + offset(log(road_density)),
                  data = dat,
                  ziformula = ~ 1,
                  family = "poisson")
  pval_count <- summary(temp)$coefficients$cond["sev", "Pr(>|z|)"]
  outDF <- data.frame(sev = i, pval_count = pval_count)
  sevDF <- rbind(sevDF, outDF)
}
sevDF %>% arrange(pval_count)
# Add em to the dataset
wvcGrid$sev63 <- mem[,63]
wvcGrid$sev1 <- mem[,1]
wvcGrid$sev4 <- mem[,4]


# FIT AND CHECK MODEL ----------------------------------------------
countModKDE <- glmmTMB(wvc_count_deer_5yrs ~ 
                         pumaKDE_mean +   
                         devProp +
                         sev4,
                       ziformula = ~ 1,
                       data = wvcGrid, 
                       family = "poisson")

countModKDE_2 <- glmmTMB(wvc_count_deer_5yrs ~ 
                           pumaKDE_mean +   
                           deerOcc +
                           devProp +
                           sev4 +
                           pumaKDE_mean:deerOcc,
                       ziformula = ~ 1,
                       data = wvcGrid, 
                       family = "poisson")

summary(countModKDE)
countModKDESim <- simulateResiduals(fittedModel = countModKDE, plot = TRUE)
testUniformity(countModKDESim)
testDispersion(countModKDESim)
testZeroInflation(countModKDESim)
testSpatialAutocorrelation(countModKDESim, x = coords[,1], y = coords[,2])

# PLOT PRDICTIONS ------------------------------------------------
# prep plotting dataset
pdCountKDE <- data.frame(
  pumaKDE_mean = seq(min(wvcGrid$pumaKDE_mean), max(wvcGrid$pumaKDE_mean), length.out = 100),
  sev63 = mean(wvcGrid$sev63),
  sev4 = mean(wvcGrid$sev4),
  # devProp = mean(wvcGrid$devProp) # mean dev value
  devProp = 0.19 #dev value at max wvc count
)

predCountKDE <- predict(countModKDE, newdata = pdCountKDE, type = 'response', se.fit = TRUE)
predDatCountKDE <- data.frame(pumaKDE_mean = pdCountKDE$pumaKDE_mean, 
                              wvc = predCountKDE$fit,
                              se_upper = predCountKDE$fit + predCountKDE$se.fit,
                              se_lower = predCountKDE$fit - predCountKDE$se.fit)

# make plots
countPlotKDE <- ggplot(predDatCountKDE, aes(x = pumaKDE_mean, y = wvc)) +
  geom_ribbon(aes(ymin = se_lower, ymax = se_upper), fill = '#375cbf', alpha = 0.8) +
  geom_line() +
  theme_classic() + 
  ylab("Number of DVCs") +
  xlab("Puma UD") +
  theme(legend.position = "none")
pdf("output/occ-model/plots/raw/revision-plots/wvc-vs-kde-model.pdf",
    height = 3, width = 3)
countPlotKDE
dev.off()

wvcKDEraw <- ggplot(wvcGrid, aes(x = pumaKDE_mean, y = wvc_count_deer_5yrs)) +
  geom_point(shape = 21, size = 1.5, fill = '#375cbf', alpha = 0.75, stroke = 0.5) +
  theme_classic() + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5),
                     labels = c("0", "0.25", "0.5")) +
  scale_y_continuous(breaks = c(0, 10, 20)) +
  ylab("DVCs") +
  xlab("Puma UD") +
  theme(legend.position = "none") +
  theme(
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA)
  )
pdf("output/occ-model/plots/raw/revision-plots/wvc-vs-kde-inset.pdf",
    height = 2, width = 2)
wvcKDEraw
dev.off()

#-----/////-----/////-----/////-----/////-----/////-----
#           POINT ANALYSIS - WVC DIEL PATTERN
#-----/////-----/////-----/////-----/////-----/////-----

# PREPARE POINT DATA -----------------------------------
# Subset to WVC points across the study in the last 5 years
load(file = "data/olympic-peninsula-grid/op-rd-grid.rda", verbose = T)
rdGrid <- opGrid %>% filter(road_density > 0) 
secIn5Yrs <- 157788000
fiveYrsAgo <- max(wvc$date_time)-secIn5Yrs
wvc5yr <- wvc %>% filter(date_time > fiveYrsAgo) %>% st_filter(rdGrid)


# Add mean KDE extracted w/in 1km buffer of each point
wvcBuff1km <- wvc5yr %>% st_buffer(1000)
vecPoint1km <- vect(wvcBuff1km)
pumaKDE_point1km <- terra::extract(kde_proj, vecPoint1km, 
                                   fun = mean, na.rm = TRUE)
wvc5yr$pumaKDE_1km <- pumaKDE_point1km$slope

# Subset to just deer collisions and adjust variables
wvc5yrDeer <- wvc5yr %>% 
  filter(species %in% c("Deer")) %>%
  mutate(injuredBinary = ifelse(injured_dead > 0, 1, 0),
         time_of_day = factor(time_of_day, levels = c("Daylight","Dark","Crepuscular")),
         pumaKDE_LH = case_when(pumaKDE_1km <= quantile(pumaKDE_1km, probs = c(0.5)) ~ "Low",
                                pumaKDE_1km > quantile(pumaKDE_1km, probs = c(0.5)) ~ "High"
         ) %>% 
           factor(levels = c("Low", "High"))
  )


# FIT AND SUMMARIZE MODEL ------------------------------------
multModLH_kde <- nnet::multinom(time_of_day ~ pumaKDE_LH, data = wvc5yrDeer)
summary(multModLH_kde)
z <- summary(multModLH_kde)$coefficients / summary(multModLH_kde)$standard.errors
p <- 2 * (1 - pnorm(abs(z)))  # two-tailed test
p
multNull <- nnet::multinom(time_of_day ~ 1, data = wvc5yrDeer)
anova(multModLH_kde, multNull, test = "Chisq")

# PLOT MULT MOD -------------------------------------
# Get predictions and error using the effects package
multEff <- allEffects(multModLH_kde)
multProb <- multEff$pumaKDE_LH$prob
multSE <- multEff$pumaKDE_LH$se.prob
multLowCI <- multEff$pumaKDE_LH$lower.prob
multUpCI <- multEff$pumaKDE_LH$upper.prob

# Prep the prediction dataset for plotting
multPDat <- rbind(melt(multProb), melt(multSE), melt(multLowCI), melt(multUpCI)) %>% 
  mutate(pumaOcc = case_when(Var1 == 1 ~ "Low",
                             Var1 == 2 ~ "High") %>% 
           factor(levels = c("Low", "High")),
         TOD = str_split(Var2, pattern = "prob.") %>% 
           lapply(function(x) x[[2]]) %>% 
           unlist() %>% 
           factor(levels = c("Daylight","Dark","Crepuscular")),
         pred = str_split(Var2, pattern = "prob.") %>% 
           lapply(function(x) x[[1]]) %>% 
           unlist() %>% 
           paste("prob", sep = "")) %>% 
  dplyr::select(pumaOcc, TOD, pred, value) %>% 
  pivot_wider(names_from = pred)

# Make plot
multPlotLH <- ggplot(multPDat, aes(x = TOD, y = prob, fill = pumaOcc)) + 
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = prob - se.prob, ymax = prob + se.prob), 
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("#a2c4eb","#375cbf")) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("P(WVC occurred during diel period)") +
  xlab("Diel period")
pdf("output/occ-model/plots/raw/revision-plots/diel-vs-kde.pdf",
    height = 3, width = 3)
multPlotLH
dev.off()

