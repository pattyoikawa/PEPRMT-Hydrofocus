#####################
# PEPRMT-DAMMM RECO #
#####################
#Written by Patty Oikawa
#patty.oikawa@gmail.com

#modified by Michael Najarro on
#September 17th 2020

#This model predicts Ecosystem Respiration using a
# Dual Arrhenius Michaelis-Menten Approach.

#Note:
#Some exogenous variables are not used. However it is easier 
# to keep this same input variable stucture for all PEPRMT codes.


PEPRMT_final_sys_CO2_Reco <- function(data,
                                      theta,
                                      wetland_type) {
  library(tidyverse)
  library(magrittr)
  
  #load(file=data)
  xdata <- data
  #theta <- theta <- c(-1, 0, -1, 0, 59280)
   
  #Exogenous Variables
   Time_2 <- xdata[,1] #day of year variable continuous
   DOY_disc_2 <- xdata[,2] #day of year variable discontinuous, starts over again at beginning of ea year
   TA_2 <- xdata[,3] #daily ave (C)
   WT_2 <- xdata[,4] #water table height (cm)
   #PAR_2 <- xdata[,5] #photosynthetically active radiation (umol m-2 d-1)
   #LAI_2 <- xdata[,6] #Leaf area index
   GPP_2 <- xdata[,7] #gpp modeled--g C m-2 d-1 (negative for uptake)
   GI_2 <- xdata[,8]#xdata[,8] #greeness index from Phenocam or Landsat etc
   Season_drop_2 <- xdata[,9]#xdata[,9] #Season variable that is set to 1 in winter (DOY 1-88, 336-365), 2 pre-spring (DOY 89-175), 3 spring (DOY 176-205), 4 summer (DOY 206-265), 5 fall (DOY 266-335)
   #wc_90CI_2 <- xdata[,10] # 90 percent confidence interval around NEE determined by gapfilling error and random error
   wetland_age_2 <- xdata[,11] #age in years of wetland; only important if wetland is >4yrs old, otherwise just set to 10
   #FPAR <- xdata[,12] #If using LAI data, set FPAR variable to 1's, if using a greeness index set FPAR to 0's
   
  # ## In current version the -1 below is commented out!
  # #Static C allocation theme
   NPPsum_avail_2 <- c(GPP_2)#*-1) #g C m-2 day-1 change to + numbers & give Reco access to all GPP
  
  #SET UP Reco
  # C_allocation <- theta[1]
  alpha1 <- 3e3 #g C m-2 d-1;--SET AS CONSTANT; Initially 3e3 umol CO2 m-2 s-1=3e3*10^-6*44.01*0.273*60*60*24=3e3 g C m-2 d-1
  ea1 <- (theta[1]+18)*1000 # 5/9/15 initial=19 Joule mol-1, 22.5 for SOM C pool
  km1 <- theta[2]+1.8e3 # g C m-3; ORIGINALLY 1.5e8 umol C m-3 soil = 1.2e-5 g C cm-3 #smaller? Km for SOM C pool
  #km1_gCcm3 <- km1*1e-6*12*1e-6
  alpha2 <- 3e3 #g C m-2 d-1 --SET AS CONSTANT
  ea2 <- (theta[3]+17.5)*1000 #5/9/15 initial 17 convert to Joule, 20.7
  km2 <- theta[4]+38 # g C m-3; ORIGINALLY 3.17e6umol C m-3 soil; 1e-6 umol C m-3 soil = 1.2e-5 g C cm-3 #larger? Km for Ps C pool
  #km2_gCcm3 <- km2*1e-6*12*1e-6
  
  #initialize C pools
  #C1_init <- theta[5] #total C avail in g C m-3; for peat soils BD~0.38g soil cm-3 (Rongzhong Davis team)
  #, with 15.6# C runs around 5.9e4 g C m-3
  #C2_init <- 0 #in g C m-3
  
  C1_init <- theta[5] # all decomposed soil organic matter in top meter from MEM in gC m-3 
  # inclusive of current year
  C2_init <- 0 #in g C d-1 labile/GPP pool starts as 0       
  
  #Reco inhibited when WT high
  a1 <- 0.00033
  a2 <- 0.0014
  a3 <- 0.75
  
  #Time Invariant
  R <- 8.314 #J K-1 mol-1
  RT <- R* (TA_2 + 274.15) #T in Kelvin-all units cancel out
  Vmax1 <- alpha1* exp(-ea1/RT) #g C m-2 d-1 SOM
  Vmax2 <- alpha2* exp(-ea2/RT) #g C m-2 d-1 labile
  
  #preallocating space
  S1sol <- vector("numeric",length(Time_2))
  S2sol <- vector("numeric",length(Time_2))
  R1 <- vector("numeric",length(Time_2))
  R2 <- vector("numeric",length(Time_2))
  S1 <- vector("numeric",length(Time_2))
  S2 <- vector("numeric",length(Time_2))
  percent_reduction <- vector("numeric",length(Time_2))
  percent_enhancement <- vector("numeric",length(Time_2))
  Reco_1 <- vector("numeric",length(Time_2))
  Reco_full <- vector("numeric",length(Time_2))
  #NPP_1 <- vector("numeric",length(Time_2))
  #NPP_full <- vector("numeric",length(Time_2))
  Ps <- vector("numeric",length(Time_2))
  C2in <- vector("numeric",length(Time_2))
  percent_available <- vector("numeric",length(Time_2))
  
  for (t in 1:length(Time_2)) { 
    #C allocation
    C2in[t] <- NPPsum_avail_2[t] # gC m-2 d-1
    #     GPPsum_avail[t] <- NPPsum_avail_2[t]*0.5 #only 50% of GPP is available to methanogens
    #     C2in_ch4[t] <- GPPsum_avail[t]
    
    #if (t == 1 | site_change[t]>0 #if beginning of model or switch sites, start C1 pool over
    if (t == 1) {
      S1[t] <- C1_init #substrate avail NOT affected by water avail-- SOM pool
      S2[t] <- C2_init + C2in[t]  # Ps C pool-- some initial Ps C lingering in soil + day 1 GPPavail
    } else {
      S1[t] <- S1sol[t-1]
      S2[t] <- S2sol[t-1] #substrate availability based on Ps on time step previous
    }
    
    #Empirical factor for increased availability of SOC during the first 3 yrs following restoration
    # insert logic tree here for peatland or wetland switch!!
    #just commenting this out for now, not sure if I can keep the peatland stuff in here as theta 5 has changed meanings in diff models
   
    
    # if(wetland_type == "Peatland" &wetland_age_2[t]<1){
     # percent_available[t] <- 0.6
      #  } else {
      #percent_available[t] <- 0.2 #only 20# of this pool is available
    #} else {}
    
    #S1[t] = S1[t]*percent_available[t]  #SOM pool
    
    
    #following Davidson and using multiple eq for diff substrate pools
    R1[t] <- Vmax1[t] * S1[t]/(km1 + S1[t]) #g C m2 d-1 Reaction velocity
    R2[t] <- Vmax2[t] * S2[t] /(km2 + S2[t]) #g C m2 d-1   
    if (R1[t]<0) {R1[t]=0}
    
    if (R2[t]<0) {R2[t]=0}
    
    #Reco is reduced by 25% when WT is at or above soil surface
    #--McNicol Silver 2015
    #   a1 <- 0.00033
    #   a2 <- 0.0014
    #   a3 <- 0.75
    #   WT_ex <- c(-30, -20, -10, 0)
    #   percent_red_ex <- c(1, 0.85, 0.77, 0.75)
    
    #   plot(percent_red_ex ~ WT_ex)
    
    # insert logic tree here for peatland or wetland switch!!
    #Also commenting out for now as we won't be able to run the same PEPRMTs for tidal and peatlands
    #if(wetland_type == "Peatland"){
     # percent_reduction[t] <- (a1*WT_2[t]^2) - (a2*WT_2[t]) + a3
    #  if (WT_2[t]>5) {percent_reduction[t] <- 0.75}
     # if (percent_reduction[t]>1.25) {percent_reduction[t] <- 1.25}
    #  if (percent_reduction[t]<0.75) {percent_reduction[t] <- 0.75}
      
    #  R1[t] <- R1[t]*percent_reduction[t] #g C m2 d-1  Reaction velocity
    #  R2[t] <- R2[t]*percent_reduction[t] #g C m2 d-1 
    #} else {}
    
    #Empirical factor for elevated Reco during the first 3 yrs following restoration
    if (wetland_age_2[t]<4) {
      percent_enhancement[t] <- 1.2
    } else {
      percent_enhancement[t] <- 1
    }
    
    R1[t] = R1[t]*percent_enhancement[t] #umol m2 sec Reaction velocity
    R2[t] = R2[t]*percent_enhancement[t] #umol m2 sec
    
    if (t==1) {
      S1sol[t] = C1_init - (R1[t]) #accounts for depletion of C sources in soil due to Reco and methane production
      S2sol[t] = (C2_init+C2in[t]) - (R2[t])
    } else {
      S1sol[t] = S1sol[t-1] - (R1[t])
      S2sol[t] = (S2sol[t-1]+C2in[t])- (R2[t])
    }
    
    if (S1sol[t]<0) {S1sol[t] <- 0}
    if (S2sol[t]<0) {S2sol[t] <- 0}
    
    ###########EDITED OUT IN CURRENT VERSION############
    #in autumn time or season 6, labile PS C pool empties into SOM
    # if (Season_drop_2[t]>5) {
    #   S1sol[t] = S1sol[t]+(0.2*S2sol[t]) #move part of labile C into SOM pool--mimicing plant matter dying
    #   S2sol[t] = S2sol[t]-(0.2*S2sol[t])
    # }
    # 
    # #in winter time or season 1, labile PS C pool empties into SOM
    # if (Season_drop_2[t]<2) {
    #   S1sol[t] = S1sol[t]+(S2sol[t]) #move entire labile C into SOM pool--mimicing plant matter dying
    #   S2sol[t] = 0
    # }
    ########################################
    
    Reco_1[t] <- R1[t] + R2[t] 
    Reco_full[t] <- (R1[t]) + (R2[t]) #umol m2 d-1
  }
  
  NEE_mod <- GPP_2 + Reco_1 #umol m-2 d-1
  
  Reco_output <- cbind(Reco_full, NEE_mod, S1, S2) %>%
    as.data.frame(.)
  #ss1 <- Reco_output$S1
  #ss2 <- Reco_output$S2
  
  #CH4_input <<- read.csv(file="./Data_sets/EL_csv") %>%
  #  select(-X) %>%
  #  mutate(S1 = ss1, S2 = ss2)
  
  
  #save(Reco_output,file="./Data_sets/CH4_input.Rdata")
  return(Reco_output)
}

