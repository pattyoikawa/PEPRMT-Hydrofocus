########################################################
### PEPRMT-DAMM CH4 model V7.2.2: Sulfate Inhibition ### 
########################################################
#Based on Davidson 2012 DAMM model

#Written by Patty Oikawa, Michael Najarro
#patty.oikawa@gmail.com, mnajarro3@horizon.csueastbay.edu

#About the model:
#1. Originally PEPRMT Ch4 was parameterized for restored wetlands
#   in the Sacramento-San Joaquin River Delta; i.e. freshwater.
#2. All PEPRMT models use the same input structure (data) for CH4 models
#   however not all models use all variables in the structure;
#3. All variables are at the daily time step.
#4. You need to run the Reco model first to get the SOC pools.


# What is new in V7.2.2 model?
#1. 7.2.2 implements Debjani Sihi's Sulfate inihibition DAMM model
#   by adding additional Mikalis-Menten Equation for Oxygen.
#   - Each MM equation is included in Patty's DAMM equations M1 & M2
#     (previously just M1).
#   - The values for the half saturation constants, salinity
#     unit conversions, and concentrations are found within the code.

#2. The downscaled, daily GPP is fed to both PEPRMT Reco and Methane,
#   not observed (tower/eddy covariance) GPP.

#3. Removal of the for-loop that controls the Sihi microsite sampling (see 5)
#   - Removes the multiple iteration and averaging of PEPRMT-CH4 results.
#   - Removes the input parameter x, which was the # of times you wish
#     to iterate the entire methane model, to simulate microsites.

#   - We decided that Sihi microsatellite approach to sampling and applying PEPRMT CH4
#     is no longer necessary.
#   - By coercion, this removes the multiple iteration loop that controls the iteration
#     of the entire model and generation fo carbon sampling from normal pools. (i.e. the
#     the number of times you want to run the methane loop on the data.)
#     Therefore, both the normal PDF sampling, the outside loop, and all code
#     used to average multiple runs and combine them into a sinle list are commented out.


#Inputs: 
#1. Theta: a vector of 7 important inputs to PEPRMT CH4.
#            these should come from FME.

#2. Data: a data frame containing 9 variables as decribed at start of function.
 
#3. Wetland_type
# - "wetland_type == 1" correpsonds to a "peatland wetland", which activates
#     the % reduction loop. 
# - "wetland_type == any other integer" corresponds to a "coastal_wetland" and
#     turns off the % reduction loop.

#Outputs:
# a data frame that containing
# 1. pulse_emission_total: total amount of methane, sum of plant + water flux.
# 2. Plant_flux_net: net maount of methane released from plants.
# 3. Hydro_flux: net amount of methane that bubbles from soil through water.
# 4. M1: amount of methane from pool 1, the labile pool
# 5. M2: amount of methane from pool 2, the SOC pool
# 6. trans2:...

#units of output
# 1. pulse_emission_total: g C methane/ m^-2 day^-1
# 2. Plant_flux_net: " 
# 3. Hydro_flux: "
# 4. M1: g C/ m^-2 day^-1
# 5. M2: "
# 6. trans2:...

library(tidyverse)
library(magrittr)


CH4_daily_step <- function(theta,
                           data,
                           #x,
                           wetland_type){
  #Constants
  R = 8.314 #J K-1 mol-1
  
  #Exogenous Variables
  Time_2 =data[,1] # day of year
  DOY_disc_2=data[,2] #day of year that starts over every year
  TA_2 = data[,3] #Air temperature- measured (C)
  WT_2 = data[,4] #water table height (cm) equals 0 when water table at soil surface
  #PAR_2 = data(:,5)#Photosynthetically active radiation (umol m-2 d-1)
  #LAI_2 = data(:,6)#Leaf area index
  GPP_2 = data[,5] #Modeled GPP - use output from PEPRMT-GPP (gC m-2 day-1) where - = uptake
  #GI_2=data(:,8)#Greeness index daily average
  #Season_drop_2=data(:,9)#not used here
  #wc_90CI_2=data(:,10)
  wetland_age_2=data[,6]
  #FPAR=data(:,12)#not used here
  Sal <- data[,7]
  S1_2=data[,8]#Modeled SOC pool - use output from PEPRMT Reco model -cummulative grams m-2
  S2_2=data[,9]#Modeled labile C pool - use output from PEPRMT Reco model -cummulative grams m-2
  
  
  WT_2_adj=(WT_2/100)+1 #makes a new variable where wt=1 at soil surface
  
  
  #CH4 PARAMETERS
  #SOC pool
  M_alpha1 = 9e10# gC m-3 d-1 originally 9e10 umol m-3 s-1
  M_ea1 = (theta[1]+65.2)*1000# parameter in kJ mol-1 multiplied by 1000 = J mol-1
  M_km1 =theta[2]+2.2e-10#g C m-3 originally 1.8e-5 umol m-3 
  #Labile C pool
  M_alpha2 = 9e11# gC m-3 d-1 originally 9e11 umol m-3 s-1
  M_ea2 = (theta[3]+68)*1000
  M_km2 =theta[4]+1.5e-11#g C m-3 originally 1.3e-6 umol m-3 
  
  #CH4 oxidation parameters
  M_alpha3 = 9e10# originally 9e10 umol m-3 s-1=7.8e15 umol m-3 d-1
  M_ea3 = (theta[5]+59.6)*1000
  M_km3 =theta[6]+253e-5
  
  #Salinity sulfate parameters
  #kISo4: inhibition coefficient of oxygen for Sulfate production
  #conc_so4: concentration of sulfate in sea water at 35 psu, 20C, according to Michael E.Q. Pilson: An Introduction to the Chemistry of the Sea, 2nd edition
  kiSO4 = theta[7] #old: 87.99  mg L^-1; Jeong et al. 2009
  conc_so4 = 2.712 # grams sulfate / kilgram solution (I interpret as sea water)
  
  #Oxygen mikalis-menten equation parameters
  #KIO2: Inhibition coefficient of oxygen on CH4 production;
  # units: micromol O2/ L water, converted to mg/L water.
  # source: debjani Sihi paper, 2019
  #conc_O2: concentration of available oxygen in SF bay water
  #         sit2 29, near san west san mateo bridge.
  #         measure = average calculated oxygen con. for 2016.
  #         no measures for 2010 exist.
  #  units: mg/ L bay water.
  #  Source: Water quality measurements in San Francisco Bay by the U.S. Geological Survey,1969–2015, Tara S. Schraga1& James E. Cloern, 2017 
  #  https://www.nature.com/articles/sdata201798.pdf
  KIO2 = (0.43/10^-6) * 15.99943
  conc_O2 = 7.429435
  
  #parameter for hydrodynamic flux
  k=0.04 #gas transfer velocity(m day-1)
  #parameter for plant-mediated transport
  Vtrans=0.24#gas transfer velocity through plants(m d-1)
  Oxi_factor=0.35#percent oxidized during transport
  
  #empirical factors for inhibition of ch4 production when WT drops
  beta1=0.48
  beta2=-0.18
  beta3=0.0042
  #empirical factors for decaying inhibition of CH4 production following first
  #flooding of wetland
  zeta1=5.1e-6 #7.4e-6
  zeta2=0.00058
  zeta3=0.11
  
  GPP_2=GPP_2*-1 
  
  GPPmax=max(GPP_2) #parameter for plant-mediated transport
  
  #Time Invariant
  RT = R * (TA_2 + 274.15)#T in Kelvin-all units cancel out
  M_Vmax1 = M_alpha1 * exp(-M_ea1/RT)#g C m-2 d-1 
  M_Vmax2 = M_alpha2 * exp(-M_ea2/RT)#gC m-2 d-1 
  M_Vmax3 = M_alpha3 * exp(-M_ea3/RT)#gC m-2 d-1 
  
  #preallocating space
  S1sol = vector('numeric',length(Time_2))
  S2sol = vector('numeric',length(Time_2))
  
  M1 = vector('numeric', length(Time_2))
  M2 = vector('numeric', length(Time_2))
  M1_full = vector('numeric', length(Time_2))
  M2_full = vector('numeric', length(Time_2))
  M_full = vector('numeric', length(Time_2))
  M_percent_reduction= vector('numeric', length(Time_2))
  M_percent_reduction_2= vector('numeric', length(Time_2))
  
  CH4water= vector('numeric', length(Time_2))
  Hydro_flux= vector('numeric', length(Time_2))
  Plant_flux= vector('numeric', length(Time_2))
  Plant_flux_net= vector('numeric', length(Time_2))
  CH4water_store= vector('numeric', length(Time_2))
  CH4water_0= vector('numeric', length(Time_2))
  Oxi_full= vector('numeric', length(Time_2))
  R_Oxi= vector('numeric', length(Time_2))
  CH4water_0_2= vector('numeric', length(Time_2))
  #Vtrans=zeros(1,length(Time_2))
  #Oxi_factor=zeros(1,length(Time_2))
  trans2=vector('numeric', length(Time_2))
  S1= vector('numeric', length(Time_2))
  S2= vector('numeric', length(Time_2))
  
  #---CREATE A SPACE FOR TO COLLECT ITERATIVE RESULTS---#
  outcome_lst <- vector('list', nrow(data))   
  # outcome_m1 <- vector('list', nrow(data))   
  # outcome_m2 <- vector('list', nrow(data))   
  # outcome_hydro <- vector('list', nrow(data))   
  # outcome_plant <- vector('list', nrow(data))
  # outcome_methane_water <- vector('list', nrow(data))
  # outcome_roxi <- vector('list', nrow(data))
  # 
  # #------ITERATIVE LOOP 10K TIMES---------#
  #for(i in 1:x){
  #   
  #   #Recreate S1 & S2 from random sample
  #   # from Normal distb. of Carbon each iteration. 
  #   S1 <- rnorm(length(Time_2), mean = mean(S1_2), sd = mean(S1_2) + .06*mean(S1_2))
  #   S2 <- rnorm(length(Time_2), mean = mean(S2_2), sd = mean(S2_2) + .06*mean(S2_2))   
  
  #--METHANE TRANSPORT ACROSS DATA---
  for(t in 1:length(Time_2)) {
    
    #parameter for plant-mediated transport--function of GPP
    trans2[t]=((GPP_2[t]+(GPPmax))/GPPmax)-1
    if (trans2[t]<0) {trans2[t]=0}
    if (trans2[t]>1) {trans2[t]=1}
    
    #calculate your sea water density values
    #pw: density of standard mean ocean water (Michael E.Q. Pilson,"An Intro to the Chemistry of the Sea)
    #     from chapter 3.6, Salinity and Densiy; units in kg / m^3
    pw = 999.842594 +(6.793952*(10^-2)*TA_2[t])-(9.095290*(10^-3)*(TA_2[t])^2) + (1.001685*(10^-4)*(TA_2[t])^3) - (1.120083*(10^-6)*(TA_2[t]^4)) + (6.536332*(10^-9)*(TA_2[t]^5))
    
    p1 = ((8.24993*10^-1) - ((4.0894*10^-3)*TA_2[t]) + ((7.6438*10^-5)*TA_2[t]^2) - ((8.2467*10^-7)*TA_2[t]^3) +((5.3875*10^-9)*TA_2[t]^4)*Sal[t])
    p2 = (-5.72466*10^-3) +(1.0227*10^-4)*TA_2[t] -(((1.6546*10^-6)*TA_2[t]^2)*Sal[t]^(3/2)) + (4.8314*10^-4* Sal[t]^2)
    seawater_density = pw + p1 + p2
    
    #calculate your ratio of sulfate to salinity
    # full unit conversion is: g sulfate/kg seawater * kg seawater/g salt
    sul_ratio = conc_so4/Sal[t] # g sulfate/g salinity 
    
    # now calculate the available sulfate concentration, which 
    # converts the salinity measured to sulfate.
    # full unit conversion:
    # g salt/kg seawater * kg seawater/m^3 seawater * g sulfate/g salt * 1000 mg sulfate/ g sulfate * m^3 seawater/1000L seawter * 1000L seawter/ 1L seawater 
    conc_so4AV = Sal[t] * seawater_density * sul_ratio * 1000 *1/1000 *1000 # mg sulfate/L seawater
    
    # BELOW IS THE OLD WAY OF CALCULATING THE AVAILABLE SULFATE CONCENTRATION!
    #convert your salinity value to a concentration of sulfate.
    #conc_so4Av: concentration of available sulfate at enzyme reaction site; mg sulfate/L seawater
    # unit conversion: start with PSU = g salt/kg sea water
    # g salt/kg seawater *1kg seawater/1000g seawter = g salt/g seawater
    # g salt/g seawater * density seawater(g seawater/ ml seawter) = g salt/ml seawater
    # g salt/ml seawater *  %sulfate in sea salt * 1000 mg sulfate/1 g sulfate * 1000 ml seawater/1L seawater  = mg sulfate/L seawater
    #conc_so4AV = Sal[t] * (1/1000) * 1.02813 * 0.08 * 1000^2
    #....08 is too small; available s04 must be a smaller number...
    
    #get the concentration of oxygen into micrograms per liter.
    
    
    #following Davidson and using multiple eq for diff substrate pool
    #implement Debjani Sihi's Salinity model
    M1[t] = M_Vmax1[t] * (S1_2[t] /(M_km1 + S1_2[t])) * (kiSO4 /(kiSO4 + conc_so4AV)) * (KIO2/(KIO2 + conc_O2)) #gC m2 d-1 Reaction velocity
    
    M2[t] = M_Vmax2[t] * S2_2[t] /(M_km2 + S2_2[t])  * (kiSO4 /(kiSO4 + conc_so4AV)) * (KIO2/(KIO2 + conc_O2))   #gC m2 d-1
    
    if (M1[t]<0) { M1[t]=0}
    if (M2[t]<0) { M2[t]=0}
    
    # # Empirical eq Oikawa for CH4 inhibition when WT falls below soil  -----
    # #surface--if WT below soil surface any time in 10days previous--CH4
    # #production reduced
    
    if(wetland_type == 1){
      if (t<=20 & WT_2_adj[t]<0.9){ # original=1
        M_percent_reduction[t]=(beta1*WT_2_adj[t]^2)+(beta2*WT_2_adj[t])+beta3
      } else {M_percent_reduction[t]=1}
      
      if (t>20){ Sel=WT_2_adj[(t-19):t] } else {Sel=5}
      
      if (t>20 & min(Sel, na.rm=T)<0.9) { 
        M_percent_reduction[t]=(beta1*WT_2_adj[t]^2)+(beta2*WT_2_adj[t])+beta3
      } else if (t>20){ M_percent_reduction[t]=1}
      
      if (WT_2_adj[t]<0){  M_percent_reduction[t]=0}
      
      if (M_percent_reduction[t]<0){  M_percent_reduction[t]=0}
      if (M_percent_reduction[t]>1){  M_percent_reduction[t]=1}
      
      #Empirical eq Oikawa for CH4 inhibition following restoration
      if (wetland_age_2[t]<2){
        M_percent_reduction_2[t]=(zeta1*DOY_disc_2[t]^2)+(zeta2*DOY_disc_2[t])+zeta3
      } else { M_percent_reduction_2[t]=1 }
      
      if (M_percent_reduction_2[t]>1) { M_percent_reduction_2[t]=1 }
      if (M_percent_reduction[t] < 0.75){M_percent_reduction[t] = 0.75}
      
      M1[t] = M1[t]*M_percent_reduction[t]  #gC m2 d Reaction velocity
      M2[t] = M2[t]*M_percent_reduction[t]  #gC m2 d
      M1[t] = M1[t]*M_percent_reduction_2[t] #gC m2 d Reaction velocity
      M2[t] = M2[t]*M_percent_reduction_2[t] #gC m2 d
    } else {}
    
    #S1sol and S2sol are the new SOC and labile pools adjusted for C lost throug CH4
    S1sol[t] = S1[t] - (M1[t]) #accounts for depletion of C sources in soil due to Reco and methane production
    S2sol[t] = S2[t] - (M2[t])
    
    if (S1sol[t]<0) { S1sol[t]=0}
    if (S2sol[t]<0) { S2sol[t]=0}
    
    #fSal <- Sal_slope*Sal[t] + Sal_int
    #if (Sal[t] >= 12.90005) {fSal <- Sal_slope1*Sal[t] + Sal_int1} #Sets fSal for salinities greater or equal to 12.9ppt
    #if (Sal[t] < 12.90005) {fSal <- Sal_slope2*Sal[t] + Sal_int2} #Sets fSal for salinities less than 12.9ppt
    
    M_full[t]=(M1[t]+M2[t])#*fSal #total CH4 produced at this time step in gC m-3 soil day-1
    
    #---COMPUTE CH4 TRANSPORT---
    #make sure WT_2_adj is never negative
    if (WT_2_adj[t]<0) { WT_2_adj[t]=0 }
    
    #now start ch4 transport loop
    if (t==1){
      
      # WT_2_adj = 1 at soil surface
      if (WT_2_adj[t]>1) {
        
        #Methane oxidation is zero when WT above surface
        R_Oxi[t] = 0
        Oxi_full[t]=0#umol m-3 d-1
        
        #This assumes you start out the year with no CH4 in water
        #where CH4water_0 is the initial concentration of CH4 in water
        CH4water_0[t]=0 # umol m-3
        CH4water_0_2[t]=0 # CH4 produced in previous time step
        #Only modeling CH4 dissolved in the water that goes down 1 m3 into
        #soil
        CH4water[t]= ((M_full[t]*1)+(CH4water_0[t]*WT_2_adj[t]))/ WT_2_adj[t]#gC per m^3 - concentration in water doesn't change
        #multiply CH4 water (mol m-3) by water wt (m) to get to mol m-2, then add
        #to M_full (umol m-2 d-1)
        #based on the concentrations in the soil and water, you get hydro and
        #plant-mediated fluxes
        Hydro_flux[t]=k*CH4water[t]  #umol m-3 30min-1 Hydrodynamic flux Poindexter
        Plant_flux[t]=(Vtrans*CH4water[t])*trans2[t]  #umol m-3 30min-1 Plant mediated transport Tian 2001--which just uses CH4 in soil
        Plant_flux_net[t]=Plant_flux[t]*Oxi_factor  #umol m-3 30min-1 Plant mediated transport after oxidation
        
        #subtract the moles of methane lost from the pools (soil and
        #water) to the atm
        CH4water_store[t]=CH4water[t]-Hydro_flux[t]-Plant_flux[t]  #umol m-3 stored in the system
        
        #If you start out the year with no water above the surface      
      } else {
        
        # Gives CH4 concentration in water which is now less that 1 m33
        CH4water_0[t]=(M_full[t]*1)+(0.00001*WT_2_adj[t])/WT_2_adj[t]#umol per m^3 - concentration in soil and water are the same
        
        # Methane oxidation turns on when WT falls below the surface
        R_Oxi[t] = M_Vmax3[t] * CH4water_0[t] /(M_km3 + CH4water_0[t])  #umol m2 day-1 Reaction velocity
        Oxi_full[t]=R_Oxi[t] # gC m-3 d-1
        CH4water[t]=CH4water_0[t]-Oxi_full[t] # now you have less ch4
        if (CH4water[t]<0) {  CH4water[t]=0 }
        
        CH4water_0_2[t]=0
        
        #this hydroflux uses a 2nd k parameter which basically inhibits
        #diffusive flux when there is no water above the soil surface
        Hydro_flux[t]=k*CH4water[t]  #gC m-3 d-1 Hydrodynamic flux Poindexter
        Plant_flux[t]=(Vtrans*CH4water[t])*trans2[t]  #umol m-2 30min-1 Plant mediated transport Tian 2001--which just uses CH4 in soil
        Plant_flux_net[t]=Plant_flux[t]*Oxi_factor  #umol m-2 30min-1 Plant mediated transport after oxidation
        
        CH4water_store[t]=CH4water[t]-Plant_flux[t]-Hydro_flux[t]  #umol m-3 stored in the system
      }
      
      # If you have water, then the CH4 should mix between the 2 layers and concentrations should be the same in water and soil
    } else {
      
      if (WT_2_adj[t]>1) {
        
        #account for any changes in concentration of CH4 due to any change in
        #WT_2_adj height
        CH4water_0[t]=(CH4water_store[t-1] * WT_2_adj[t-1]) / WT_2_adj[t]
        
        #Methane oxidation is zero when WT above surface
        R_Oxi[t] = 0
        Oxi_full[t]=0 #umol m-3 d-1
        CH4water_0_2[t]=0
        #Now add the new CH4 produced today to the soil and let it increase in concentration
        CH4water[t]= ((M_full[t]*1)+(CH4water_0[t]*WT_2_adj[t]))/ WT_2_adj[t]#umol per m^3 - concentration in water doesn't change
        
        #again compute fluxes
        Hydro_flux[t]=k*CH4water[t]  # umol m-2 30min-1 Hydrodynamic flux Poindexter
        Plant_flux[t]=(Vtrans*CH4water[t])*trans2[t]  # umol m-2 30min-1 Plant mediated transport Tian 2001--which just uses CH4 in soil
        Plant_flux_net[t]=Plant_flux[t]*Oxi_factor  # umol m-2 30min-1 Plant mediated transport after oxidation
        
        #subtract the moles of methane lost from the pools (soil and #water) to the atm
        CH4water_store[t]=CH4water[t]-Plant_flux[t]-Hydro_flux[t] #umol m-3 stored in the system
        
      } else {
        #if you don't have WT_2_adj above the soil, then all the CH4 in the water goes to the soil
        #First, account for increased concentration of CH4 in soil now that WT_2_adj has dropped
        CH4water_0[t]=(CH4water_store[t-1] * WT_2_adj[t-1]) / WT_2_adj[t]
        
        #now add new CH4 to this new concentrated pool of CH4 in the soil
        CH4water_0_2[t]= ((M_full[t]*1)+(CH4water_0[t]*WT_2_adj[t]))/ WT_2_adj[t]#gc per m^3 - concentration in water doesn't change
        
        #Methane oxidation turns on when WT falls below the surface
        R_Oxi[t] = M_Vmax3[t] * CH4water_0_2[t] /(M_km3 + CH4water_0_2[t])  #gC m2 d Reaction velocity
        Oxi_full[t]=R_Oxi[t]#gC m-3 day-1
        CH4water[t]=CH4water_0_2[t]-Oxi_full[t]#now you have less ch4
        if (CH4water[t]<0) {CH4water[t]=0}
        
        
        Hydro_flux[t]=k*CH4water[t]  #umol m-3 30min-1 Hydrodynamic flux Poindexter
        Plant_flux[t]=(Vtrans*CH4water[t])*trans2[t]  #umol m-2 30min-1 Plant mediated transport Tian 2001--which just uses CH4 in soil
        Plant_flux_net[t]=Plant_flux[t]*Oxi_factor  #umol m-2 30min-1 Plant mediated transport after oxidation
        
        CH4water_store[t]=CH4water[t]-Plant_flux[t]-Hydro_flux[t] #umol m-3 stored in the system
      }
    }
  }
  
  pulse_emission_total = Plant_flux_net+Hydro_flux#*fSal  #gC CH4 m-2 day-1 total CH4 flux to atm
  
  return(data.frame(cbind(pulse_emission_total,
                          Plant_flux_net,
                          Hydro_flux,
                          M1,
                          M2,
                          trans2)))
  
  #   # take avg of pulse_emission_total, Plant_flux_net, Hydro_flux,M1,M2, trans2 and store in a vector  
  #   outcome_lst[[i]] <- c(pulse_emission_total) 
  #   outcome_m1[[i]] <- c(M1)
  #   outcome_m2[[i]] <- c(M2)
  #   outcome_plant[[i]] <- c(Plant_flux_net)
  #   outcome_hydro[[i]] <- c(Hydro_flux)
  #   outcome_methane_water[[i]] <- c(CH4water)
  #   outcome_roxi[[i]] <- c(R_Oxi)
  # #}
  # 
  # #merge your iterative pulse emissions
  # avg_ch4_outputs <- do.call('rbind', outcome_lst)
  # avg_ch4_outputs2 <- avg_ch4_outputs %>%
  #   as.data.frame() %>%
  #   summarise_all(mean)
  # 
  # u <- as.numeric(avg_ch4_outputs2[1,])
  # uu <- data.frame(cbind(c(1:length(u)), u))
  # 
  # #merge your M1s
  # avg_m1_op <- do.call('rbind', outcome_m1)
  # avg_m1_op <- avg_m1_op %>%
  #   as.data.frame() %>%
  #   summarise_all(mean)
  # 
  # y <- as.numeric(avg_m1_op[1,])
  # #colnames(y) <- c("M1")
  # 
  # #merge your M2s
  # avg_m2_op <- do.call('rbind', outcome_m2)
  # avg_m2_op <- avg_m2_op %>%
  #   as.data.frame() %>%
  #   summarise_all(mean)
  # 
  # v <- as.numeric(avg_m2_op[1,])
  # #colnames(v) <- c("M2")
  # 
  # #merge your hydro
  # avg_hydro <- do.call('rbind', outcome_hydro)
  # avg_hydro <- avg_hydro %>%
  #   as.data.frame() %>%
  #   summarise_all(mean)
  # 
  # w <- as.numeric(avg_hydro[1,])
  # #colnames(w) <- c("hydro_flux")
  # 
  # #merge your plant flux
  # avg_plant <- do.call('rbind', outcome_plant)
  # avg_plant <- avg_plant %>%
  #   as.data.frame() %>%
  #   summarise_all(mean)
  # 
  # s <- as.numeric(avg_plant[1,])
  # #colnames(s) <- c("Plant_flux")
  # 
  # #merge your methane water
  # avg_ch4water <- do.call('rbind', outcome_methane_water)
  # avg_ch4water <- avg_ch4water %>%
  #   as.data.frame() %>%
  #   summarise_all(mean)
  # 
  # g <- as.numeric(avg_ch4water[1,])
  # 
  # #merge your roxi
  # avg_roxi <- do.call('rbind', outcome_roxi)
  # avg_roxi <- avg_roxi %>%
  #   as.data.frame() %>%
  #   summarise_all(mean)
  # 
  # h <- as.numeric(avg_roxi[1,])
  # 
  # ch4_output<- as.data.frame(cbind(uu,y,v,w,s,g,h))
  # colnames(ch4_output) <- c("index",
  #                           "avg_Pulse_Emission_Total",
  #                           "M1",
  #                           "M2",
  #                           "hydro_flux",
  #                           "Plant_flux",
  #                           "ch4_water",
  #                           "r_oxi")
  # return(ch4_output)
}