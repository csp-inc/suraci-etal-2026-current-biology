# Suraci et al. 2026. Large carnivores shape road safety through effects on prey space use. Current Biology

### Organization: Conservation Science Partners, Inc.

### Contact: Dr. Justin Suraci ([justin\@csp-inc.org](mailto:justin@csp-inc.org))

This repo contains all code required to replication the analyses described in the above publication.

#### *Abstract*

Quantifying the ecosystem services provided by large carnivores through their top-down effects on prey is an important strategy for promoting human-carnivore coexistence. Such “predation services” may include reductions in wildlife-vehicle collisions mediated by changes in the abundance or distribution of ungulate prey. However, empirical demonstrations are few and present limited opportunity to understand the mechanisms by which carnivores may influence road safety through their impacts on prey space use. We leveraged a large camera trap dataset and spatial data on wildlife-vehicle collisions to examine how pumas shape the spatio-temporal habitat use of ungulate prey and in turn affect the probability and frequency of collisions across the Olympic Peninsula of Washington, USA. Occupancy models revealed that, where pumas were present, deer were more diurnally active, were 15% less likely to occur in high-road density areas, and were 86% more likely to occupy remote habitat away from development (and associated traffic). These predator-related changes in deer habitat use were associated with a 67% reduction in the probability of a deer-vehicle collision and a 76% reduction in the total expected number of collisions over a five-year period where puma occupancy probability was high. Puma collar data corroborated these findings and further indicated a shift in deer-vehicle collision timing towards daylight hours where puma space use was high, corresponding to observed changes in deer diel activity. Thus, even where predator and prey populations are well established, pumas may substantially reduce collision risk by altering deer behavior, benefiting society through increased road safety.

#### *Analysis Workflow*

This analysis consists of two primary component: 
1. Modeling the effects of puma presence on deer and elk spatiotemporal habitat use using Continuous Time Occupancy (CTO) models. Scripts in the [cto-analysis](cto-analysis/) directory:
    - Prepare and examine [model covariates](cto-analysis/00-assess-covariate-soe.R)
    - Fit CTO models for puma interactions with [deer](cto-analysis/01-ct-occ-model-deer-puma.R) and [elk](cto-analysis/01-ct-occ-model-elk-puma.R)
2. Quantifying the effects of puma presence and activity on deer-vehicle collisions. Scripts in the [dvc-analysis](dvc-analysis/) directory:
    - Calculate [puma and deer occupancy probability](dvc-analysis/00-calc-psi-at-cam-locs-v2.R) (from CTO model results) at locations across the study area
    - Prepare [data on deer-vehicle collisions](dvc-analysis/01-prep-wvc-dataset-v2.R)
    - Run [primary](dvc-analysis/02-wvc-modelling-v4.R) and [corroborative](dvc-analysis/03-wvc-vs-puma-kde.R) models of puma effects on DVC probability and counts  

