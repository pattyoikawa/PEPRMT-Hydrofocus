#MEM inputs for Eden Landing site
# created: 9/15/2020 
# by: Michael Najarro

# purpose: to simply bring your the required 26 inputs
# needed for MEM (June 2020 version) to work
# excludes all but 1 optional variable.

# note: this function requires no inputs.

var_input <- function() {
  startYear<<-2018
  endYear<<-startYear+99
  meanSeaLevel<<-98.1
  meanSeaLevelDatum<<- 0#meanSeaLevel
  relSeaLevelRiseInit<<-0.9
  relSeaLevelRiseTotal<<-250
  rslrTotal<<-250
  initElv<<-198.1
  meanTidalHeight<<-194.7
  meanHighWater<<-194.7
  meanHighWaterDatum<<- #213.7 see noaa tide gauge or schile paper??
  meanHighHighWater<<-213.7
  meanHighHighWaterSpring<<-287
  suspendedSediment<<-.000025
  lunarNodalAmp<<-2.5
  bMax<<-.048
  zVegMin<<-100#80
  zVegMax<<-350#160
  zVegPeak<<-180#130
  plantElevationType<<-"orthometric"
  rootToShoot<<-2
  rootTurnover<<-0.5
  rootDepthMax<<-30
  shape<<-"linear"
  omDecayRate<<-0.2
  recalcitrantFrac<<-0.2
  settlingVelocity<<-2.8
  omPackingDensity<<-0.082
  rootPackingDensity<<-omPackingDensity
  mineralPackingDensity<<-2.43
  coreYear<<-2050
  
  # #initial elevation value~~
  # core <- read.csv(file="./data/CCRCN_cores.csv")
  # rr_core <- core %>%
  #   filter(site_id == "Rush_Ranch") %>%
  #   separate(core_date, into = c("year", "month", "day"), sep = "-")
  # 
  #   initElv <<- rr_core %>%
  #   summarise(initelv = mean(core_elevation)) %>%
  #   pull()
  # 
  # initElv <<- round(initElv *100, digits = 4)
  # #~~~~
  
}