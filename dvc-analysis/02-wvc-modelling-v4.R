## ---------------------------
##
## Script name: 02-wvc-modelling.R
##
## Purpose of script: Fit, check, and plot
## predictions from puma-DVC models
##
## Author: Justin Suraci
##
## ---------------------------

library(tidyverse)
library(sf)
library(spdep)
library(adespatial)
library(glmmTMB)
library(nnet)
library(effects)
library(DHARMa)
library(reshape2)
library(patchwork)
library(MuMIn)
library(performance)

# Load relevant datasets -------------------------------------------------
load(file = "data/deer-vehicle-collision-data/wvc-analysis-dataset-20251109.rda", verbose = T)
wvcGrid <- wvcGrid %>% mutate(wvc_binDeer = ifelse(wvc_count_deer > 0, 1, 0),
                              wvc_binDeer_5yrs = ifelse(wvc_count_deer_5yrs > 0, 1, 0),
                              wvc_binDeer_10yrs = ifelse(wvc_count_deer_10yrs > 0, 1, 0))
# Subset to WVCs across the entire peninusal in the last 5 years
load(file = "data/olympic-peninsula-grid/op-rd-grid.rda", verbose = T)
rdGrid <- opGrid %>% filter(road_density > 0) 
secIn5Yrs <- 157788000
fiveYrsAgo <- max(wvc$date_time)-secIn5Yrs
wvc5yr <- wvc %>% filter(date_time > fiveYrsAgo) %>% st_filter(rdGrid)

#-------------------------------------------------------------------------------
# WVC COUNT GRID-BASED MODELS --------------------------------------------------
#-------------------------------------------------------------------------------
# ACCOUNTING FOR SPATIAL STRUCTURE
# Prep spatial eigenvectors
coords <- wvcGrid %>% select(x, y) %>% st_drop_geometry()
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
  temp2 <- glmmTMB(injured_dead_5yrs ~ sev + offset(log(road_density)),
                  data = dat,
                  ziformula = ~ deerOcc,
                  family = "poisson")
  pval_count <- summary(temp)$coefficients$cond["sev", "Pr(>|z|)"]
  pval_injury <- summary(temp2)$coefficients$cond["sev", "Pr(>|z|)"]
  outDF <- data.frame(sev = i, pval_count = pval_count, pval_injury = pval_injury)
  sevDF <- rbind(sevDF, outDF)
}
sevDF %>% arrange(pval_count)
sevDF %>% arrange(pval_injury)
# Add em to the dataset
wvcGrid$sev1 <- mem[,1]


# FIT AND CHECK MODEL ----------------------------------------------
countMod <- glmmTMB(wvc_count_deer_5yrs ~ 
                         pumaOcc + 
                         deerOcc + 
                         devProp +
                         sev1 +
                         pumaOcc:deerOcc, 
                       ziformula = ~ 1,
                       data = wvcGrid, 
                       family = "poisson")

summary(countMod)
countModSim <- simulateResiduals(fittedModel = countMod, plot = TRUE)
testUniformity(countModSim)
testDispersion(countModSim)
testZeroInflation(countModSim)
testSpatialAutocorrelation(countModSim, x = coords[,1], y = coords[,2])


# PLOT PRDICTIONS ------------------------------------------------
# prep plotting dataset
pdCountLow <- data.frame(
  # pumaOcc = seq(0.06, 0.99, length.out = 100),
                         pumaOcc = seq(min(wvcGrid$pumaOcc), max(wvcGrid$pumaOcc), length.out = 100),
                         sev63 = mean(wvcGrid$sev63),
                         sev1 = mean(wvcGrid$sev1),
                         # devProp = mean(wvcGrid$devProp),
                         devProp = 0.19, #dev value at max wvc count
                         road_density = mean(wvcGrid$road_density),
                         deerOcc = quantile(wvcGrid$deerOcc, probs = c(0.25)))
pdCountHigh <- data.frame(
  # pumaOcc = seq(0.06, 0.99, length.out = 100),
                          pumaOcc = seq(min(wvcGrid$pumaOcc), max(wvcGrid$pumaOcc), length.out = 100),
                          sev63 = mean(wvcGrid$sev63),
                          sev1 = mean(wvcGrid$sev1),
                          # devProp = mean(wvcGrid$devProp),
                          devProp = 0.19,
                          road_density = mean(wvcGrid$road_density),
                          deerOcc = quantile(wvcGrid$deerOcc, probs = c(0.75)))
predCountLow <- predict(countMod, newdata = pdCountLow, type = 'response', se.fit = TRUE)
predCountHigh <- predict(countMod, newdata = pdCountHigh, type = 'response', se.fit = TRUE)
predDatCountLow <- data.frame(pumaOcc = pdCountLow$pumaOcc, 
                              wvc = predCountLow$fit,
                              se_upper = predCountLow$fit + predCountLow$se.fit,
                              se_lower = predCountLow$fit - predCountLow$se.fit, 
                              deerCat = 'Low')
predDatCountHigh <- data.frame(pumaOcc = pdCountHigh$pumaOcc, 
                               wvc = predCountHigh$fit,
                               se_upper = predCountHigh$fit + predCountHigh$se.fit,
                               se_lower = predCountHigh$fit - predCountHigh$se.fit, 
                               deerCat = "High")
predDatCount <- rbind(predDatCountLow, predDatCountHigh)

# make plots
countPlot <- ggplot(predDatCount, aes(x = pumaOcc, y = wvc, fill = deerCat, lty = deerCat)) +
  geom_ribbon(aes(ymin = se_lower, ymax = se_upper), alpha = 0.8) +
  geom_line() +
  scale_fill_manual(values = c('#7f32a8', '#ce93ed')) +
  theme_classic() + 
  ylab("Number of DVCs") +
  xlab("Puma occupancy probability") +
  labs(fill = "Deer occupancy", lty = "Deer occupancy") +
  theme(legend.position = "none")
  

#-------------------------------------------------------------------------------
# BINARY WVC LOGISTIC GRID-BASED MODELS ----------------------------------------
#-------------------------------------------------------------------------------
# DEAL W/ SPATIAL STRUCTURE/AUTOCORRELATION
# Prep spatial eigenvectors
coords <- wvcGrid %>% select(x, y) %>% st_drop_geometry()
nb <- dnearneigh(coords, d1 = 0, d2 = 12000)
plot(nb, coords)
lw <- nb2listw(nb, style = "W")
mem <- scores.listw(lw, MEM.autocor = "positive")

# Identify SEV with strongest influence
sevDF2 <- data.frame()
for(i in 1:ncol(mem)){
  dat <- wvcGrid
  dat$sev <- mem[,i]
  temp <- glm(wvc_binDeer_5yrs ~ sev, 
                  data = dat, 
                  family = "binomial")
  pval <- summary(temp)$coefficients["sev", "Pr(>|z|)"]
  outDF <- data.frame(sev = i, pval = pval)
  sevDF2 <- rbind(sevDF2, outDF)
}
sevDF2 %>% arrange(pval)
# Add em to the dataset
wvcGrid$sev2_2 <- mem[,2]
wvcGrid$sev2_66 <- mem[,66]

# FIT AND CHECK MODEL ----------------------------------------------
# Fit model
binMod <- glm(wvc_binDeer_5yrs ~ pumaOcc +
                I(pumaOcc^2) +
                devProp,
              data = wvcGrid, family = "binomial")
summary(binMod)
# Check binMod
binModSim <- simulateResiduals(fittedModel = binMod, plot = TRUE)
testUniformity(binModSim)
testDispersion(binModSim)
testZeroInflation(binModSim)
testSpatialAutocorrelation(binModSim, x = coords[,1], y = coords[,2])

# PLOT WVC PRES/ABS LOGISTIC GRID-BASED MODELS --------------------------------------
# prep plotting dataset
pdBin <- data.frame(pumaOcc = seq(min(wvcGrid$pumaOcc), max(wvcGrid$pumaOcc), length.out = 100),
                    # sev2_2 = mean(wvcGrid$sev2_2),
                    road_density = mean(wvcGrid$road_density),
                    # devProp = mean(wvcGrid$devProp)
                    devProp = 0.19 #dev value at max wvc count
                    )
predBin <- predict(binMod, newdata = pdBin, type = 'response', se.fit = TRUE)
predDat <- data.frame(pumaOcc = pdBin$pumaOcc, 
                      wvc = predBin$fit,
                      se_upper = predBin$fit + predBin$se.fit,
                      se_lower = predBin$fit - predBin$se.fit)
wvcMax <- predDat$pumaOcc[predDat$wvc == max(predDat$wvc)]

# make plots
binPlot <- predDat %>% ggplot(aes(x = pumaOcc, y = wvc)) +
  geom_ribbon(aes(ymin = se_lower, ymax = se_upper), fill = "#fcd7d7") +
  geom_line() +
  # geom_vline(xintercept = 0.55, linetype = 2) +
  # geom_vline(xintercept = 0.72, linetype = 2) +
  theme_classic() + 
  ylab("P(DVC occurrence)") +
  xlab("Puma occupancy probability")

#-------------------------------------------------------------------------------
# MULTINOMIAL MODEL FOR WVC TOD ------------------------------------------------
#-------------------------------------------------------------------------------
# PREP THE DATA ----------------------------------------------
wvcCam5yrs <- wvc5yr %>% filter(nn_cam_dist <= 2500) %>% 
  filter(species %in% c("Deer")) %>%
  mutate(time_of_day = factor(time_of_day, levels = c("Daylight","Dark","Crepuscular")),
         daytime = ifelse(time_of_day == "Daylight", 1, 0),
         daytime2 = ifelse(sunIllum == 0, 0, 1),
         pumaOccLH = case_when(pumaOccProb <= 0.5 ~ "Low",
                                pumaOccProb > 0.5 ~ "High",
                               ) %>%
           factor(levels = c("Low", "High")),
         injuredBinary = ifelse(injured_dead > 0, 1, 0)
         )

# FIT AND SUMMARIZE MODEL ----------------------------------------------
multModLH <- nnet::multinom(time_of_day ~ pumaOccLH, data = wvcCam5yrs)
summary(multModLH)
z <- summary(multModLH)$coefficients / summary(multModLH)$standard.errors
p <- 2 * (1 - pnorm(abs(z)))  # two-tailed test
p
multNull <- nnet::multinom(time_of_day ~ 1, data = wvcCam5yrs)
anova(multModLH, multNull, test = "Chisq")

# PLOT MULT MOD -------------------------------------
# Get predictions and error using the effects package
multEff <- allEffects(multModLH)
multProb <- multEff$pumaOccLH$prob
multSE <- multEff$pumaOccLH$se.prob
multLowCI <- multEff$pumaOccLH$lower.prob
multUpCI <- multEff$pumaOccLH$upper.prob

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
  select(pumaOcc, TOD, pred, value) %>% 
  pivot_wider(names_from = pred)

# Make plot
multPlotLH <- ggplot(multPDat, aes(x = TOD, y = prob, fill = pumaOcc)) + 
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = prob - se.prob, ymax = prob + se.prob), 
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("#f7c6ad","#db5916")) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("P(WVC occurred during diel period)") +
  xlab("Diel period")

# ----------------------------------------------------------------
# SAVE OUTPUTS ---------------------------------------------------
# ----------------------------------------------------------------
save(countMod, binMod, multModLH,
     file = "output/wvc-models.rda")