# Creating a moonscape soil for Eden landing
library(tidyverse)

# 1. If rCTM is loaded and in the memory, forget rCTM
if ("rCTM" %in% (.packages())){
  detach("package:rCTM", unload=TRUE) 
}

# 2. If remotes is not already installed, install it
if (! ("remotes" %in% installed.packages())) {
  install.packages("remotes")
}


# 4. Load most current version into memory
library(rCTM)

library(gridExtra)

rootDepthMax = 30
omPackingDensity = 0.085
mineralPackingDensity = 1.99
rootPackingDensity = omPackingDensity

moonscapeCohorts <- data.frame(age=rep(0, round(rootDepthMax+0.6)), # set all ages to 0
                      fast_OM=rep(0, round(rootDepthMax+0.6)), # Note that fast, slow, and root mass initial conditions are all set to 0
                      slow_OM=rep(0, round(rootDepthMax+0.6)),
                      respired_OM=rep(0, round(rootDepthMax+0.6)),
                      mineral=rep(mineralPackingDensity,round(rootDepthMax+0.6)), # Mineral mass is 100% so equal to the packing density
                      root_mass=rep(0,round(rootDepthMax+0.6)),
                      layer_top=0:((round(rootDepthMax+0.6)-1)), # 1 cm depth increments
                      layer_bottom=1:round(rootDepthMax+0.6)) %>% 
  dplyr::mutate(cumCohortVol = cumsum(layer_bottom-layer_top))

edenLandingMem <- rCTM::runCohortMem(startYear = 2010,
                                     relSeaLevelRiseInit = 0.24, 
                                     relSeaLevelRiseTotal = 40, 
                                     initElv = 165, 
                                     meanSeaLevel = 110, 
                                     meanSeaLevelDatum = 98.1, 
                                     meanHighWaterDatum = 194.7, 
                                     meanHighHighWaterDatum = 213.3, 
                                     meanHighHighWaterSpringDatum = NA, 
                                     suspendedSediment = .000045, 
                                     lunarNodalAmp = 2.5, 
                                     lunarNodalPhase = 2011.181, 
                                     nFloods = 705.79, 
                                     bMax = 0.09, 
                                     zVegMin = 80, 
                                     zVegMax = 350, 
                                     zVegPeak = NA, 
                                     plantElevationType = "orthometric",
                                     rootToShoot = 2, 
                                     rootTurnover = 0.5, 
                                     abovegroundTurnover = NA, 
                                     speciesCode = NA, 
                                     rootDepthMax = 30, 
                                     shape = "linear", 
                                     omDecayRate = 0.2, 
                                     recalcitrantFrac = 0.35, 
                                     captureRate = 0.275, 
                                     omPackingDensity = 0.085, 
                                     mineralPackingDensity = 1.99, 
                                     initialCohorts = moonscapeCohorts, 
                                     uplandCohorts = NA, 
                                     supertidalCohorts = NA, 
                                     supertidalSedimentInput = NA
                                     )

plot(edenLandingMem$annualTimeSteps$year,
     edenLandingMem$annualTimeSteps$meanSeaLevel, 
     type="l")

plot(edenLandingMem$annualTimeSteps$year,
     edenLandingMem$annualTimeSteps$surfaceElevation, type="l")

animateCohorts(cohorts = edenLandingMem$cohorts, scenario = edenLandingMem$annualTimeSteps,
               filename = "EdenLandingMem.gif")

core <- simulateSoilCore(cohorts=edenLandingMem$cohorts, coreYear=2020)

accretionPlot <-  ggplot(data = core, aes(x=layer_bottom, y=(layer_bottom-layer_top)/input_yrs)) +
  geom_point() +
  geom_line() +
  xlab("Depth (cm)") +
  ylab("Accretion Rate (cm/yr)") +
  scale_x_reverse() +
  coord_flip()

loiPlot <- ggplot(data = core, aes(x=layer_bottom, y=om_fraction)) +
  geom_point() +
  geom_line() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  ylab("Organic Matter (fraction)") +
  scale_x_reverse() +
  coord_flip()

bdPlot <- ggplot(data = core, aes(x=layer_bottom, y=dry_bulk_density)) +
  geom_point() +
  geom_line() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  ylab(expression(paste("Bulk Density (g cm"^"-3",")",sep=""))) +
  scale_x_reverse() +
  coord_flip()

grid.arrange(accretionPlot, loiPlot, bdPlot, nrow=1)

ggplot(edenLandingMem$annualTimeSteps, aes(x = year, y = cFlux*10000)) +
  geom_line()

ggplot(edenLandingMem$annualTimeSteps, aes(x = year, y = cSequestration*10000)) +
  geom_line()
