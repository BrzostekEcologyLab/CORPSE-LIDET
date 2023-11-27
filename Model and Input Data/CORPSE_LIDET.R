### This code models C loss during decomposition of leaf litters in the LIDET experiment with the CORPSE model

# the code has options to use either baseline, LIDET, or any of the other 8 best parameter sets identified in the 
# modified Monte Carlo analysis for litter decomposition in the litter and litterbag layers

# V2
# C and N data get exported as model output 

# NOTE: This code has various options for litter decomposition parameter sets, 
#       the desired one needs to be uncommented to run (and all others commented)

# Remove all functions, clear memory
rm(list=ls(all=TRUE)) 

# Load Packages
library(readr)
library(dplyr)

# set the working directory
setwd() # add the filepath for the working directory


# Set the parameters for the run
# params data table below is for litterbag and litter with LIDET or baseline parameters, or other 8 sets
# params_bulk data table has only the baseline parameters for the bulk and rhizosphere

# Vector of sites
sites<- c("AND","BNZ","BSF","CDR","CPR","HBR","HFR","JUN","KBS","KNZ","NWT","OLY","SEV","SMR","UFL","VCR")

# Define BASELINE Parameters 
paramset <- "Baseline" # used to name files at end of run
Vmaxreffast <- 33 
Vmaxrefslow <- 0.6 
Easlow <- 30E3 
kCfast <- 0.01 
kCslow <- 0.01 
kCnecro <- 0.01 
eupslow <- 0.1 

# Define new LIDET Parameters 
# paramset <- "LIDET" # (this is set 108, "LIDET" used to name files at end of run)
# Vmaxreffast <- 2.5
# Vmaxrefslow <- 0.35
# Easlow <- 32.5e3
# kCfast <- 0.007
# kCslow <-  0.007
# kCnecro <- 0.009
# eupslow <- 0.1

# Set 13
# paramset <- "Set13" # used to name files at end of run
# Vmaxreffast <- 3.25
# Vmaxrefslow <- 0.35
# Easlow <- 25e3
# kCfast <- 0.015
# kCslow <-  0.015
# kCnecro <- 0.01
# eupslow <- 0.02
  
# Set 31
# paramset <- "Set31" # used to name files at end of run
# Vmaxreffast <- 2
# Vmaxrefslow <- 0.3
# Easlow <- 35e3
# kCfast <- 0.01
# kCslow <-  0.011
# kCnecro <- 0.011
# eupslow <- 0.04

# Set 64 
# paramset <- "Set64" # used to name files at end of run
# Vmaxreffast <- 2.75
# Vmaxrefslow <- 0.35
# Easlow <- 35e3
# kCfast <- 0.009
# kCslow <-  0.011
# kCnecro <- 0.011
# eupslow <- 0.001

# Set 75
# paramset <- "Set75" # used to name files at end of run
# Vmaxreffast <- 3.25
# Vmaxrefslow <- 0.35
# Easlow <- 25e3
# kCfast <- 0.011
# kCslow <-  0.005
# kCnecro <- 0.01
# eupslow <- 0.02

# Set 113
# paramset <- "Set113" # used to name files at end of run
# Vmaxreffast <- 2 
# Vmaxrefslow <- 0.35 
# Easlow <- 25E3 
# eupslow <- 0.001 
# kCfast <- 0.007 
# kCslow <- 0.005 
# kCnecro <- 0.007 


# Set 120
# paramset <- "Set120" # used to name files at end of run
# Vmaxreffast <- 2.5
# Vmaxrefslow <- 0.35
# Easlow <- 25e3
# kCfast <- 0.015
# kCslow <-  0.005
# kCnecro <- 0.01
# eupslow <- 0.01

# Set 122
# paramset <- "Set122" # used to name files at end of run
# Vmaxreffast <- 3
# Vmaxrefslow <- 0.3
# Easlow <- 25e3
# kCfast <- 0.007
# kCslow <- 0.007
# kCnecro <- 0.011
# eupslow <- 0.1

# Set 134
# paramset <- "Set134" # used to name files at end of run
# Vmaxreffast <- 2
# Vmaxrefslow <- 0.2
# Easlow <- 25e3
# kCfast <- 0.01
# kCslow <-  0.015
# kCnecro <- 0.009
# eupslow <- 0.02


# set the number of sites for the loop
nsites <- as.numeric(length(sites))


# CORPSE functions -----------------------------------------

  #Input Parameters
  ## Data frame with definitions of expected parameters 
  params_definitions <- data.frame( 
    "Vmaxref" = as.character ('Relative maximum enzymatic decom rates (Fast, Slow, Necro)'), 
    "Ea" = as.character ('Activation engery (Fast, Slow, Necro)'),
    "kC" = as.character ('Michealis-Menton parameter (Fast, Slow, Necro)'),
    "gas_diffusion_exp" = as.character ('Determines suppression of decomp at high soil moisture'), 
    "minMicrobeC" = as.character ('Minimum microbial biomass (fraction of total C)'), 
    "Tmic"= as.character ('Microbial lifetime at 20C (years)'),
    "et" = as.character ('Fraction of turnover not converted to CO2'),
    "eup" = as.character ('Carbon uptake efficiency (Fast, Slow, Necro)'),
    "nup" = as.character ('Nitrogen uptake efficiency (Fast, Slow, Necro)'),
    "tProtected" = as.character ('Protected C turnover time (years)'), 
    "frac_N_turnover_min" = as.character ('Fraction of microbial biomass N turnover that is mineralized'), 
    "protection_rate" = as.character ('Protected carbon formation rate (year-1)'), 
    "CN_Microbe" = as.character ('C:N ratio of microbial biomass for AM sites'), 
    "max_immobilization_rate" = as.character ('Maximum N immobilization rate (fraction per day)'),
    "substrate_diffusion_exp" = as.character ('Determines suppression of decomp at low soil moisture'),
    "frac_turnover_slow" = as.character('Fraction of microbial biomass N turnover that goes to slow pool'),
    "new_resp_units" = as.character ('If TRUE, Vmaxref has units of 1/years and assumes optimal soil moisture has a relative rate of 1.0'),
    "iN_loss_rate" = as.character ('Loss rate of inorganic N pool (year-1) > 1 because it takes much less than a year for it to be removed'))
  
  
  ## Data frame with CORPSE parameters including new LIDET litter decomposition parameters
  params <- data.frame( "Vmaxref_Fast"= Vmaxreffast, 
                        "Vmaxref_Slow"=Vmaxrefslow , 
                        "Vmaxref_Necro"= 22, 
                        "Ea_Fast"= 5e3 , 
                        "Ea_Slow"= Easlow, 
                        "Ea_Necro"= 3e3, 
                        "kC_Fast"= kCfast, 
                        "kC_Slow"= kCslow, 
                        "kC_Necro"= kCnecro, 
                        "gas_diffusion_exp"= 0.6, 
                        "minMicrobeC"= 1e-3, 
                        "Tmic"= 0.25, 
                        "et"= 0.6, 
                        "eup_Fast"=  0.6, 
                        "eup_Slow"= eupslow, 
                        "eup_Necro"= 0.6, 
                        "tProtected"= 100.0,  
                        "frac_N_turnover_min"= 0.2,
                        "protection_rate_Fast"= 0, 
                        "protection_rate_Slow"= 0,  
                        "protection_rate_Necro"= 0, 
                        "nup_Fast"= 0.3, 
                        "nup_Slow"= 0.3,  
                        "nup_Necro"= 0.3, 
                        "CN_Microbe"= 6.1, 
                        "max_immobilization_rate"= 3.65,  
                        "substrate_diffusion_exp"= 1.5, 
                        "new_resp_units"= TRUE, 
                        "iN_loss_rate"= 5.0,  
                        "frac_turnover_slow"= 0.2)
  
  ## Data frame of CORPSE parameters with baseline litter decomposition parameters
  params_bulk <- data.frame( "Vmaxref_Fast"= 33, 
                             "Vmaxref_Slow"=  0.6, 
                             "Vmaxref_Necro"= 22, 
                             "Ea_Fast"= 5e3 , 
                             "Ea_Slow"= 30e3, 
                             "Ea_Necro"= 3e3, 
                             "kC_Fast"= 0.01, 
                             "kC_Slow"= 0.01, 
                             "kC_Necro"= 0.01, 
                             "gas_diffusion_exp"= 0.6, 
                             "minMicrobeC"= 1e-3, 
                             "Tmic"= 0.25, 
                             "et"= 0.6, 
                             "eup_Fast"=  0.6, 
                             "eup_Slow"=0.1 , 
                             "eup_Necro"= 0.6, 
                             "tProtected"= 100.0,  
                             "frac_N_turnover_min"= 0.2,
                             "protection_rate_Fast"= 0.6, 
                             "protection_rate_Slow"= 0.001,  
                             "protection_rate_Necro"= 4, 
                             "nup_Fast"= 0.3, 
                             "nup_Slow"= 0.3,  
                             "nup_Necro"= 0.3, 
                             "CN_Microbe"= 6.1, 
                             "max_immobilization_rate"= 3.65,  
                             "substrate_diffusion_exp"= 1.5, 
                             "new_resp_units"= TRUE, 
                             "iN_loss_rate"= 5.0,  
                             "frac_turnover_slow"= 0.2)
  
  
  paramvals<-c(params$Vmaxref_Fast, 
               params$Vmaxref_Slow,
               params$Ea_Slow,
               params$eup_Slow,
               params$kC_Fast,
               params$kC_Slow,
               params$kC_Necro)


  ## CORPSE leaf and root litter parameters
  litter_transfer_to_soil <- 0.1/365
  rhizoC_flux <- 0.0354/365 # based on Fahey et al 2005 estimate of NPP in late 1990s: 354 g/m2/yr (0.354 kg/m2/yr), rhizo C flux is 10% of NPP, added uniformly over the year
  rhizo_frac <- 0.15
  chem_types<-c('Fast','Slow','Necro')
  claymod <- 1.0  #Scalar that modifies the ability of clays to sorb and protect C 
  inorg_Ndep <- 0.01 # inorganic N deposition
  
  ##Protection rate parameter based on percent clay content from Mayes et al (2012) Table 1
  prot_clay <- function (claypercent, slope=0.4833, intercept=2.3282, BD=1.15, porosity=0.4) {
    prot<-1.0*(10**(slope*log10(claypercent)+intercept)*BD*1e-6)
    return (prot)
  }
  
  ##Calculate rates of change for ALL CORPSE pools 
  ##T = Temperature (Kelvin)
  ##Theta = soil water content (fraction of saturation)
  
  ##Next 3 functions are needed to run model code
  ##Functions are ordered to allow code to run
  
  ##Function to calculate Vmax of microbial decomposition 
  ##Vmax function, normalized to Tref=293.15 (T is in Kelvin)
  Vmax <- function (T, params, Tref = 293.15, Rugas = 8.314472) {
    Fast <- params$Vmaxref_Fast*exp(-params$Ea_Fast*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))
    Slow <- params$Vmaxref_Slow*exp(-params$Ea_Slow*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))
    Necro <- params$Vmaxref_Necro*exp(-params$Ea_Necro*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))
    Vmax<-data.frame(Fast,Slow,Necro)
    return(Vmax)
  }
  
  ##Function to sum carbon types  
  sumCtypes <- function(SOM, prefix, suffix='C') {
    out <- SOM[[paste(prefix, chem_types[1], suffix, sep='')]]
    if (length(chem_types)>1){
      for (i in 2:length(chem_types)) {
        out <- out+SOM[[paste(prefix,chem_types[i],suffix,sep='')]]
        
      }
    }
    return(out)
  }
  
  
  ##Function to calculate Decomposition rate 
  decompRate <- function(SOM, T, theta, params){
    theta[theta<0] <- 0.0
    theta[theta>1.0] <- 1.0
    vmax <- Vmax(T,params)
    decompRate <- data.frame(matrix(nrow=nspp,ncol=0))
    for (ctypes in chem_types) {
      drate <- vmax[[ctypes]]*theta*(SOM[[paste('u',ctypes,'C',sep='')]]*SOM$livingMicrobeC/((sumCtypes(SOM,'u'))*params[[paste('kC',ctypes,sep='_')]]+SOM$livingMicrobeC))
      zeroMB = SOM$livingMicrobeC==0.0
      drate[zeroMB] <- 0.0
      decompRate[paste(ctypes,'C',sep='')] <- drate
      NC_ratio <- SOM[[paste('u',ctypes,'N',sep='')]]/SOM[[paste('u',ctypes,'C', sep='')]]
      NC_ratio[SOM[[paste('u',ctypes,'C', sep='')]]==0.0]<-0.0
      decompRate[paste(ctypes,'N',sep='')] <- drate*NC_ratio
    }
    return (decompRate)
  }
  
# CORPSE code  
  CORPSE_lidet <- function(SOM,T,theta,params,claymod=1.0) {
    ##Calculate maximum potential C decomposition rate
    decomp<-decompRate(SOM,T,theta,params)
    
    ##Microbial Turnover 
    microbeTurnover<-(SOM$livingMicrobeC-params$minMicrobeC*sumCtypes(SOM,'u','C'))/params$Tmic ##kg/m2/yr
    microbeTurnover[microbeTurnover<0.0]<-0.0
    
    ##Calculating maintenance respiration and creating zeros matrix for overflow_resp
    maintenance_resp <- microbeTurnover*(1.0-params$et)
    overflow_resp <- maintenance_resp*0
    
    ##Calculating fraction dead microbes in C and N
    deadmic_C_production <- microbeTurnover*params$et ##actual fraction of microbial turnover
    deadmic_N_production<-(microbeTurnover*params$et)/(params$CN_Microbe)
    
    ##C and N available for microbial growth
    carbon_supply<-vector(mode='numeric', length=nspp)
    nitrogen_supply<-vector(mode='numeric', length=nspp)
    for (ctypes in chem_types) {
      carbon_supply<-carbon_supply+decomp[[paste(ctypes,'C',sep='')]]*params[[paste('eup',ctypes,sep='_')]]
      nitrogen_supply<-nitrogen_supply+decomp[[paste(ctypes,'N',sep='')]]*params[[paste('nup',ctypes,sep='_')]]
    }
    
    IMM_N_max <- params$max_immobilization_rate*365*SOM$inorganicN/(SOM$inorganicN+params$max_immobilization_rate)
    
    dmicrobeC <- vector(mode='numeric', length=nspp)
    dmicrobeN <- vector(mode='numeric', length=nspp)
    CN_imbalance_term <- vector(mode='numeric', length=nspp)
    
    ##Growth is nitrogen limited, with not enough mineral N to support it with max immobilization
    loc_Nlim <- (carbon_supply - maintenance_resp)>((nitrogen_supply + IMM_N_max) * params$CN_Microbe)
    
    CN_imbalance_term[loc_Nlim] <- (-IMM_N_max[loc_Nlim])
    
    dmicrobeC[loc_Nlim] <- ((nitrogen_supply[loc_Nlim] + IMM_N_max[loc_Nlim]) * params$CN_Microbe - microbeTurnover[loc_Nlim] * params$et)
    
    dmicrobeN[loc_Nlim] <- (nitrogen_supply[loc_Nlim] + IMM_N_max[loc_Nlim] - microbeTurnover[loc_Nlim] * params$et/params$CN_Microbe)
    
    overflow_resp[loc_Nlim] <- carbon_supply[loc_Nlim] - maintenance_resp[loc_Nlim] - (nitrogen_supply[loc_Nlim] + IMM_N_max[loc_Nlim]) * params$CN_Microbe
    
    ##Growth must be supported by immobilization of some mineral nitrogen, but is ultimately carbon limited
    loc_immob <- (carbon_supply - maintenance_resp >= nitrogen_supply*params$CN_Microbe) & (carbon_supply - maintenance_resp < (nitrogen_supply+IMM_N_max)*params$CN_Microbe)
    
    CN_imbalance_term[loc_immob] <- (-((carbon_supply[loc_immob] - maintenance_resp[loc_immob])/params$CN_Microbe - nitrogen_supply[loc_immob]))
    
    dmicrobeC[loc_immob] <- (carbon_supply[loc_immob] - microbeTurnover[loc_immob])
    
    dmicrobeN[loc_immob] <- ((carbon_supply[loc_immob]-maintenance_resp[loc_immob])/params$CN_Microbe - microbeTurnover[loc_immob]*params$et/params$CN_Microbe)
    
    ##Growth is carbon limited and extra N is mineralized
    loc_Clim <- !(loc_Nlim | loc_immob)
    
    dmicrobeC[loc_Clim] <- (carbon_supply[loc_Clim] - microbeTurnover[loc_Clim]) 
    
    dmicrobeN[loc_Clim] <- ((carbon_supply[loc_Clim] - maintenance_resp[loc_Clim])/params$CN_Microbe - microbeTurnover[loc_Clim]*params$et/params$CN_Microbe)
    
    CN_imbalance_term[loc_Clim] <- nitrogen_supply[loc_Clim] - (carbon_supply[loc_Clim]-maintenance_resp[loc_Clim])/params$CN_Microbe;
    
    ##CO2 production and cumulative CO2 produced by cohort
    CO2prod <- maintenance_resp + overflow_resp
    for (ctypes in chem_types) {
      CO2prod <- CO2prod + decomp[[paste(ctypes,'C',sep='')]]*(1.0-params[[paste('eup',ctypes,sep='_')]])
    }
    
    ##Update protected carbon 
    protectedturnover <- 
      data.frame(matrix(nrow=nspp,ncol=0))
    
    protectedprod <- data.frame(matrix(nrow=nspp,ncol=0))
    
    for (ctypes in chem_types) {
      protectedturnover[paste(ctypes,'C',sep='')] <- SOM[paste('p',ctypes,'C',sep='')]/params$tProtected
      
      protectedprod[paste(ctypes,'C',sep='')] <- SOM[paste('u',ctypes,'C',sep='')]*params[[paste('protection_rate',ctypes,sep='_')]]*claymod 
      
      protectedturnover[paste(ctypes,'N',sep='')] <- SOM[paste('p',ctypes,'N',sep='')]/params$tProtected
      
      protectedprod[paste(ctypes,'N',sep='')] <- SOM[paste('u',ctypes,'N',sep='')]*params[[paste('protection_rate',ctypes,sep='_')]]*claymod 
    }
    
    derivs <- SOM
    
    derivs$livingMicrobeC <- dmicrobeC
    derivs$livingMicrobeN <- dmicrobeN
    derivs$CO2 <- CO2prod
    derivs$inorganicN <- CN_imbalance_term
    
    for (ctypes in chem_types) {
      derivs[[paste('u',ctypes,'C',sep='')]] <- (-decomp[[paste(ctypes,'C', sep='')]]) + protectedturnover[[paste(ctypes,'C', sep='')]] - protectedprod[[paste(ctypes,'C', sep='')]]
      
      derivs[[paste('p',ctypes,'C',sep='')]] <- protectedprod[[paste(ctypes,'C', sep='')]] - protectedturnover[[paste(ctypes,'C', sep='')]]
      
      derivs[[paste('u',ctypes,'N',sep='')]] <- (-decomp[[paste(ctypes,'N',sep='')]]) + protectedturnover[[paste(ctypes,'N', sep='')]] - protectedprod[[paste(ctypes,'N', sep='')]]
      
      derivs[[paste('p',ctypes,'N',sep='')]] <- protectedprod[[paste(ctypes,'N', sep='')]] - protectedturnover[[paste(ctypes,'N', sep='')]]
    }
    
    derivs['uNecroC'] <- derivs['uNecroC'] + deadmic_C_production * (1.0-params$frac_turnover_slow)
    
    derivs['uSlowC'] <- derivs['uSlowC'] + deadmic_C_production * params$frac_turnover_slow
    
    turnover_N_min <- deadmic_N_production * params$frac_N_turnover_min
    
    turnover_N_slow <- deadmic_N_production * params$frac_turnover_slow
    
    non_assimilated_N <- vector(mode = 'numeric', length = 1)
    
    for (ctypes in chem_types) { 
      non_assimilated_N <- non_assimilated_N + decomp[[paste(ctypes,'N',sep='')]] * (1.0 - params[[paste('nup',ctypes,sep='_')]])
    }
    
    derivs['uNecroN'] <- derivs['uNecroN'] + deadmic_N_production - turnover_N_min - turnover_N_slow
    
    derivs['uSlowN'] <- derivs['uSlowN'] + turnover_N_slow
    
    derivs['inorganicN'] <- derivs['inorganicN'] + turnover_N_min + non_assimilated_N
    
    return(derivs)
  }
  
  
  
  ##{Main Code Function, for multiple model runs}
  
  ##Function Definition: 
  ##Inputs have been separated out 
  CORPSE_loop <- function(nyears,nspp,litter_added,CORPSE_full_spinup_litter,CORPSE_full_spinup_bulk,CORPSE_full_spinup_rhizo,litter_production,froot_turnover,root_prod,litter_CN,root_CN,soil_T,soil_VWC){
    
    #Function Variables:
    timestep <- nyears*365
    CORPSEstep <- 1/365 ## this is daily values
    claymod <- 1.0  #Scalar that modifies the ability of clays to sorb and protect C 
    inorg_Ndep <- 0.01 # inorganic N deposition. 
    
    #Input files: these are read in below under "Running the model"
    # adding 100g in litterbag
    litterbag_int <- litter_added
    #Initial Data
    litter_int <- (CORPSE_full_spinup_litter)
    bulk_int <- (CORPSE_full_spinup_bulk)
    rhizo_int <- (CORPSE_full_spinup_rhizo)
    
    # litter production- amount added is in g
    litter_production <- (litter_production)

    colnames(litter_production)<-1:nspp
    froot_turnover<-litter_production
    root_prod<-litter_production
    
    # Load litter quality data (LIDET)
    # CN ratio
    CN_ratio_litter <- as.matrix((litter_CN))
    CN_ratio_root <- as.matrix((root_CN))
    
    ##soil abiotic conditions (temperature, water content)
    soilT <- as.matrix((soil_T))
    soilVWC <- as.matrix((soil_VWC))
    
    # calculate fast fraction (0-1) based on 100 g of litter added
    litter_fastfrac <- as.matrix(litter_int[,1]*10)
    
    ## Make an empty data frame with column names for pools (unprotected and protected) and chem_types (Fast, Slow, Necro) for each soil compartment
    litter <- 
      data.frame("uFastC"= numeric (nspp),"uSlowC"= numeric(nspp), "uNecroC"= numeric(nspp), "pFastC"= numeric(nspp), "pSlowC"= numeric(nspp),"pNecroC"= numeric(nspp), "livingMicrobeC"= numeric(nspp),"uFastN"= numeric(nspp), "uSlowN"= numeric(nspp), "uNecroN"= numeric(nspp), "pFastN"= numeric(nspp),"pSlowN"= numeric(nspp), "pNecroN"= numeric(nspp),"inorganicN"= numeric(nspp), "CO2"= numeric(nspp), "livingMicrobeN"= numeric(nspp)
      )
    
    litterbag <- 
      data.frame("uFastC"= numeric (nspp),"uSlowC"= numeric(nspp), "uNecroC"= numeric(nspp), "pFastC"= numeric(nspp), "pSlowC"= numeric(nspp),"pNecroC"= numeric(nspp), "livingMicrobeC"= numeric(nspp),"uFastN"= numeric(nspp), "uSlowN"= numeric(nspp), "uNecroN"= numeric(nspp), "pFastN"= numeric(nspp),"pSlowN"= numeric(nspp), "pNecroN"= numeric(nspp),"inorganicN"= numeric(nspp), "CO2"= numeric(nspp), "livingMicrobeN"= numeric(nspp)
      )
    
    bulk <- 
      data.frame("uFastC"= numeric (nspp),"uSlowC"= numeric(nspp), "uNecroC"= numeric(nspp), "pFastC"= numeric(nspp), "pSlowC"= numeric(nspp),"pNecroC"= numeric(nspp), "livingMicrobeC"= numeric(nspp),"uFastN"= numeric(nspp), "uSlowN"= numeric(nspp), "uNecroN"= numeric(nspp), "pFastN"= numeric(nspp),"pSlowN"= numeric(nspp), "pNecroN"= numeric(nspp),"inorganicN"= numeric(nspp), "CO2"= numeric(nspp), "livingMicrobeN"= numeric(nspp)
      )
    
    rhizo <- 
      data.frame("uFastC"= numeric (nspp),"uSlowC"= numeric(nspp), "uNecroC"= numeric(nspp), "pFastC"= numeric(nspp), "pSlowC"= numeric(nspp),"pNecroC"= numeric(nspp), "livingMicrobeC"= numeric(nspp),"uFastN"= numeric(nspp), "uSlowN"= numeric(nspp), "uNecroN"= numeric(nspp), "pFastN"= numeric(nspp),"pSlowN"= numeric(nspp), "pNecroN"= numeric(nspp),"inorganicN"= numeric(nspp), "CO2"= numeric(nspp), "livingMicrobeN"= numeric(nspp)
      )
    
    # load initial data into dataframes made in step above
    litter[1:nspp,1:16] <- litter_int[1:nspp,1:16] 
    litter$livingMicrobeN <- litter$livingMicrobeC/params$CN_Microbe
    
    litterbag[1:nspp,1:16] <- litterbag_int[1:nspp,1:16] 
    litterbag$livingMicrobeN <- litterbag$livingMicrobeC/params$CN_Microbe
    
    bulk[1:nspp,1:16] <- bulk_int[1:nspp,1:16] 
    bulk$livingMicrobeN <- bulk$livingMicrobeC/params_bulk$CN_Microbe
    
    rhizo[1:nspp,1:16] <- rhizo_int[1:nspp,1:16] 
    rhizo$livingMicrobeN <- rhizo$livingMicrobeC/params$CN_Microbe
    
    # Set up matrices for model output
    
    ##Set up empty matrices to hold CORPSE model outputs (as they change over time)
    
    litter_uFastC <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_uSlowC <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_uNecroC <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_pFastC <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_pSlowC <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_pNecroC <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_uFastN <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_uSlowN <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_uNecroN <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_pFastN <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_pSlowN <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_pNecroN <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_livingMicrobeC <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_inorganicN <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_CO2 <- matrix(NA, nrow=timestep, ncol=nspp)
    litter_livingMicrobeN <- matrix(NA,nrow=timestep, ncol=nspp)
    
    litterbag_uFastC <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_uSlowC <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_uNecroC <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_pFastC <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_pSlowC <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_pNecroC <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_uFastN <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_uSlowN <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_uNecroN <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_pFastN <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_pSlowN <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_pNecroN <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_livingMicrobeC <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_inorganicN <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_CO2 <- matrix(NA, nrow=timestep, ncol=nspp)
    litterbag_livingMicrobeN <- matrix(NA,nrow=timestep, ncol=nspp)
    
    bulk_uFastC <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_uSlowC <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_uNecroC <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_pFastC <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_pSlowC <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_pNecroC <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_uFastN <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_uSlowN <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_uNecroN <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_pFastN <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_pSlowN <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_pNecroN <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_livingMicrobeC <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_inorganicN <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_CO2 <- matrix(NA, nrow=timestep, ncol=nspp)
    bulk_livingMicrobeN <- matrix(NA,nrow=timestep, ncol=nspp)
    
    rhizo_uFastC <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_uSlowC <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_uNecroC <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_pFastC <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_pSlowC <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_pNecroC <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_uFastN <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_uSlowN <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_uNecroN <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_pFastN <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_pSlowN <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_pNecroN <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_livingMicrobeC <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_inorganicN <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_CO2 <- matrix(NA, nrow=timestep, ncol=nspp)
    rhizo_livingMicrobeN <- matrix(NA,nrow=timestep, ncol=nspp)
    
    
    # initializing shared inorganic N pool
    shared_inorganicN <- numeric(nspp)
    
    
    ## Call CORPSE code
    
    for (i in 1:timestep) {
      ##Start CORPSE model main loop
      
      #Parse step into DOY 
      ##this takes the timestep and makes it into DOY
      k<-(i%%365)
      if(k==0){print(k)}
      if(k==0){k<-365}
      
      ##Get Temperature and theta (soil moisture) values for this time point
      T_step <- soilT[k,1]+273.15
      theta_step <- soilVWC[k,1]
      
      ##Running the CORPSE function for LIDET data
      # Litter Layer
      litter$inorganicN <- shared_inorganicN
      results_litter <- CORPSE_lidet(litter, T_step, theta_step, params, claymod) 
      shared_inorganicN <- shared_inorganicN + CORPSEstep*results_litter$inorganicN
      
      # Litterbag layer
      litterbag$inorganicN <- shared_inorganicN
      results_litterbag <- CORPSE_lidet(litterbag, T_step, theta_step, params, claymod) 
      shared_inorganicN <- shared_inorganicN + CORPSEstep*results_litterbag$inorganicN
      
      # Bulk
      bulk$inorganicN <- shared_inorganicN
      results_bulk <- CORPSE_lidet(bulk, T_step, theta_step, params_bulk, claymod) 
      shared_inorganicN <- shared_inorganicN + CORPSEstep*results_bulk$inorganicN
      
      # Rhizo
      rhizo$inorganicN <- shared_inorganicN
      results_rhizo <- CORPSE_lidet(rhizo, T_step, theta_step, params_bulk, claymod) 
      shared_inorganicN <- shared_inorganicN + CORPSEstep*results_rhizo$inorganicN
      
      
      ## Update the pools in SOM by add derivs*dt (length of time step) to each SOM pool.  
      ## This  converts the units from mass per year to mass per the selected time step
      shared_inorganicN <- shared_inorganicN + CORPSEstep *(inorg_Ndep - params$iN_loss_rate*shared_inorganicN)
      
      litter <- litter + results_litter*CORPSEstep  
      litter$inorganicN <- shared_inorganicN 
      
      litterbag <- litterbag + results_litterbag*CORPSEstep  
      litterbag$inorganicN <- shared_inorganicN
      
      bulk <- bulk + results_bulk*CORPSEstep 
      bulk$inorganicN <- shared_inorganicN
      
      rhizo <- rhizo + results_rhizo*CORPSEstep 
      rhizo$inorganicN <- shared_inorganicN
      
      ##Save all pool values in output matrices
      litter_uFastC[i,]<-litter$uFastC
      litter_uSlowC[i,]<-litter$uSlowC
      litter_uNecroC[i,]<-litter$uNecroC
      litter_pFastC[i,]<-litter$pFastC
      litter_pSlowC[i,]<-litter$pSlowC
      litter_pNecroC[i,]<-litter$pNecroC
      litter_uFastN[i,]<-litter$uFastN
      litter_uSlowN[i,]<-litter$uSlowN
      litter_uNecroN[i,]<-litter$uNecroN
      litter_pFastN[i,]<-litter$pFastN
      litter_pSlowN[i,]<-litter$pSlowN
      litter_pNecroN[i,]<-litter$pNecroN
      litter_livingMicrobeC[i,]<-litter$livingMicrobeC
      litter_inorganicN[i,]<-litter$inorganicN
      litter_CO2[i,]<-litter$CO2
      litter_livingMicrobeN[i,]<-litter$livingMicrobeN
      
      litterbag_uFastC[i,]<-litterbag$uFastC
      litterbag_uSlowC[i,]<-litterbag$uSlowC
      litterbag_uNecroC[i,]<-litterbag$uNecroC
      litterbag_pFastC[i,]<-litterbag$pFastC
      litterbag_pSlowC[i,]<-litterbag$pSlowC
      litterbag_pNecroC[i,]<-litterbag$pNecroC
      litterbag_uFastN[i,]<-litterbag$uFastN
      litterbag_uSlowN[i,]<-litterbag$uSlowN
      litterbag_uNecroN[i,]<-litterbag$uNecroN
      litterbag_pFastN[i,]<-litterbag$pFastN
      litterbag_pSlowN[i,]<-litterbag$pSlowN
      litterbag_pNecroN[i,]<-litterbag$pNecroN
      litterbag_livingMicrobeC[i,]<-litterbag$livingMicrobeC
      litterbag_inorganicN[i,]<-litterbag$inorganicN
      litterbag_CO2[i,]<-litterbag$CO2
      litterbag_livingMicrobeN[i,]<-litterbag$livingMicrobeN
      
      bulk_uFastC[i,]<-bulk$uFastC
      bulk_uSlowC[i,]<-bulk$uSlowC
      bulk_uNecroC[i,]<-bulk$uNecroC
      bulk_pFastC[i,]<-bulk$pFastC
      bulk_pSlowC[i,]<-bulk$pSlowC
      bulk_pNecroC[i,]<-bulk$pNecroC
      bulk_uFastN[i,]<-bulk$uFastN
      bulk_uSlowN[i,]<-bulk$uSlowN
      bulk_uNecroN[i,]<-bulk$uNecroN
      bulk_pFastN[i,]<-bulk$pFastN
      bulk_pSlowN[i,]<-bulk$pSlowN
      bulk_pNecroN[i,]<-bulk$pNecroN
      bulk_livingMicrobeC[i,]<-bulk$livingMicrobeC
      bulk_inorganicN[i,]<-bulk$inorganicN
      bulk_CO2[i,]<-bulk$CO2
      bulk_livingMicrobeN[i,]<-bulk$livingMicrobeN
      
      rhizo_uFastC[i,]<-rhizo$uFastC
      rhizo_uSlowC[i,]<-rhizo$uSlowC
      rhizo_uNecroC[i,]<-rhizo$uNecroC
      rhizo_pFastC[i,]<-rhizo$pFastC
      rhizo_pSlowC[i,]<-rhizo$pSlowC
      rhizo_pNecroC[i,]<-rhizo$pNecroC
      rhizo_uFastN[i,]<-rhizo$uFastN
      rhizo_uSlowN[i,]<-rhizo$uSlowN
      rhizo_uNecroN[i,]<-rhizo$uNecroN
      rhizo_pFastN[i,]<-rhizo$pFastN
      rhizo_pSlowN[i,]<-rhizo$pSlowN
      rhizo_pNecroN[i,]<-rhizo$pNecroN
      rhizo_livingMicrobeC[i,]<-rhizo$livingMicrobeC
      rhizo_inorganicN[i,]<-rhizo$inorganicN
      rhizo_CO2[i,]<-rhizo$CO2
      rhizo_livingMicrobeN[i,]<-rhizo$livingMicrobeN
      
      
      ##t function transposes data frame
      
      ##Add leaf litter inputs to litter layer
      for (j in 1:nspp){
        litter$uFastC[j] <- litter$uFastC[j] + (t((litter_production[k,j])*litter_fastfrac[j,]))
        
        litter$uSlowC[j] <- litter$uSlowC[j] + (t((litter_production[k,j])*(1-litter_fastfrac[j,])))
        
        litter$uFastN[j] <- litter$uFastN[j]+(t(((litter_production[k,j])*litter_fastfrac[j,])/CN_ratio_litter[j,]))
        
        litter$uSlowN[j] <- litter$uSlowN[j] + (t(((litter_production[k,j])*(1-litter_fastfrac[j,]))/CN_ratio_litter[j,]))
        
        ##Add root litter inputs to the rhizosphere and bulk soil layers
        rhizo$uFastC[j] <- rhizo$uFastC[j] + (t(froot_turnover[k,j]) * litter_fastfrac[j,]*rhizo_frac)
        
        rhizo$uSlowC[j] <- rhizo$uSlowC[j] + (t(froot_turnover[k,j]) * (1-litter_fastfrac[j,])*rhizo_frac)
        
        rhizo$uFastN[j]<-rhizo$uFastN[j]+(t(froot_turnover[k,j])*litter_fastfrac[j,]*rhizo_frac/CN_ratio_root[j,])
        
        rhizo$uSlowN[j] <- rhizo$uSlowN[j] + (t(froot_turnover[k,j]) * (1-litter_fastfrac[j,])*rhizo_frac/CN_ratio_root[j,])
        
        bulk$uFastC[j] <- bulk$uFastC[j] + (t(froot_turnover[k,j])*litter_fastfrac[j,]*(1-rhizo_frac))
        
        bulk$uSlowC[j] <- bulk$uSlowC[j] + (t(froot_turnover[k,j]) * (1-litter_fastfrac[j,])*(1-rhizo_frac))
        
        bulk$uFastN[j] <- bulk$uFastN[j] + (t(froot_turnover[k,j]) * litter_fastfrac[j,]*(1-rhizo_frac)/CN_ratio_root[j,])
        
        bulk$uSlowN[j] <- bulk$uSlowN[j] + (t(froot_turnover[k,j]) * (1-litter_fastfrac[j,])*(1-rhizo_frac)/CN_ratio_root[j,])
      }
      
      ##Transfer a portion of the litter layer to the bulk and rhizosphere each time step
      newsoil <- litter * litter_transfer_to_soil
      rhizo <- rhizo + (newsoil * (rhizo_frac))
      bulk <- bulk + (newsoil * (1-rhizo_frac))
      litter <- litter - (litter * litter_transfer_to_soil)
      
      
      ##Add in non mycorrhizal C flux to rhizosphere
      rhizo$uFastC <- rhizo$uFastC + rhizoC_flux
      
      
    }
    
    ##This is a list of the information we want to get from this function
    Exported<-list(litter,litterbag,bulk,rhizo,litterbag_uFastC,litterbag_uSlowC,litterbag_uNecroC,litterbag_uFastN,litterbag_uSlowN,litterbag_uNecroN,litterbag_livingMicrobeC,litterbag_inorganicN,litterbag_CO2, litterbag_livingMicrobeN)
    
    ##These are the names for each item in the list to be called later
    names(Exported)<-c("litter","litterbag","bulk","rhizo","litterbag_uFastC","litterbag_uSlowC","litterbag_uNecroC","litterbag_uFastN","litterbag_uSlowN","litterbag_uNecroN","litterbag_livingMicrobeC","litterbag_inorganicN","litterbag_CO2", "litterbag_livingMicrobeN")
    
    return(Exported)
  }
  
  

# Running the Model -------------------------------------------------------

  
  ##INPUT DATA
  ##Reads in input data, runs the model, and writes export data to files for each site
  for(i in 1:nsites){

    # need to have folder "results_Baseline (or results_LIDET, etc) in working directory
    dir <- "results_Baseline/"
    
    # Read in input data each site. 
    # The file names for the sites must be the same except replacing site code for this to run. 
    CORPSE_full_spinup_litter <- read.csv(paste("Input Files/",sites[i],"Litterbag/CORPSE_full_spinup_litter.csv",sep=""))
    CORPSE_full_spinup_bulk <- read.csv(paste("Input Files/",sites[i],"Litterbag/CORPSE_full_spinup_bulk.csv",sep=""))
    CORPSE_full_spinup_rhizo <- read.csv(paste("Input Files/",sites[i],"Litterbag/CORPSE_full_spinup_rhizo.csv",sep=""))
    litter_production <- read.csv(paste("Input Files/",sites[i],"Litterbag/litter_production.csv",sep=""))
    froot_turnover <- read.csv(paste("Input Files/",sites[i],"Litterbag/litter_production.csv",sep=""))
    root_prod <- read.csv(paste("Input Files/",sites[i],"Litterbag/litter_production.csv",sep=""))
    litter_CN <- read.csv(paste("Input Files/",sites[i],"Litterbag/litter ",sites[i]," CN.csv",sep=""))
    root_CN <- litter_CN
    soil_T <- read.csv(paste("Input Files/",sites[i],"Litterbag/soilT ",sites[i]," DOY274start.csv", sep=""))
    soil_VWC <- read.csv(paste("Input Files/",sites[i],"Litterbag/soilVWC ",sites[i]," DOY274start.csv", sep=""))

    # Read in litterbag with 6 spp, adds 100 g of litter, same for all sites
    litter_added <- read.csv(paste("Input Files/","litterbag_init_100g_6spp.csv", sep=""))
    
    # Table of species at each site with number of years of data for each one
    # Runs each site for length of LIDET data years - i.e. HBR runs for 6 years while HFR runs for 10 years
    sitedata <- read.csv(paste("Input Files/",sites[i],"Litterbag/",sites[i],".csv",sep=""))
    
    # remove NAs in sitedata
    sitedata <- na.omit(sitedata)
    
    # Set up data to only run 6 common species done at most sites
    
    # make sure the number of rows in the initial files match the number of species
    CORPSE_full_spinup_bulk <- CORPSE_full_spinup_bulk[1:nrow(sitedata),]
    CORPSE_full_spinup_litter <- CORPSE_full_spinup_litter[1:nrow(sitedata),]
    CORPSE_full_spinup_rhizo <- CORPSE_full_spinup_rhizo[1:nrow(sitedata),]
    
    # add the column of species from "sitedata" to litter, rhizo, and bulk initial values,
    # and litter CN and root CN, then use that to filter rows for only 6 common spp
    CORPSE_full_spinup_bulk$spp <- sitedata$SPP 
    CORPSE_full_spinup_litter$spp <- sitedata$SPP 
    CORPSE_full_spinup_rhizo$spp <- sitedata$SPP 
    litter_CN$spp <- sitedata$SPP
    root_CN$spp <- sitedata$SPP
    
    # filter data for only 6 common species, and remove "spp" column
    CORPSE_full_spinup_bulk <- 
      CORPSE_full_spinup_bulk %>% 
      filter(spp == "ACSA" | spp == "DRGL" | spp == "PIRE" | spp == "QUPR" | spp == "THPL"| spp == "TRAE") %>% 
      select(-spp)
    
    CORPSE_full_spinup_litter <-
      CORPSE_full_spinup_litter %>% 
      filter(spp == "ACSA" | spp == "DRGL" | spp == "PIRE" | spp == "QUPR" | spp == "THPL"| spp == "TRAE") %>% 
      select(-spp)
    
    CORPSE_full_spinup_rhizo <- 
      CORPSE_full_spinup_rhizo %>% 
      filter(spp == "ACSA" | spp == "DRGL" | spp == "PIRE" | spp == "QUPR" | spp == "THPL"| spp == "TRAE") %>% 
      select(-spp)
    
    litter_CN <- 
      litter_CN %>% 
      filter(spp == "ACSA" | spp == "DRGL" | spp == "PIRE" | spp == "QUPR" | spp == "THPL"| spp == "TRAE") %>% 
      select(-spp)
    
    root_CN <- 
      root_CN %>% 
      filter(spp == "ACSA" | spp == "DRGL" | spp == "PIRE" | spp == "QUPR" | spp == "THPL"| spp == "TRAE") %>% 
      select(-spp)
    

    # filter sitedata for 6 common species that we present in the manuscript
    sitedata <-
      sitedata %>% 
      filter(SPP == "ACSA" | SPP == "DRGL" | SPP == "PIRE" | SPP == "QUPR" | SPP == "THPL"| SPP == "TRAE")
    
    # Because some sites are missing some of the 6 common spp, filter the litterbag for the spp present
    litter_added <- 
      litter_added %>% 
      mutate(spp = c("ACSA", "DRGL", "PIRE", "QUPR", "THPL", "TRAE")) %>% 
      filter(spp %in% sitedata$SPP ) %>% 
      select(-spp)
    
    # set columns for data equal to nspp in sitedata: froot_turnover, litter_production, root_prod
    # all columns have the same values so just remove excess columns
    froot_turnover <- froot_turnover[,1:nrow(sitedata)]
    litter_production <- litter_production[,1:nrow(sitedata)]
    root_prod <- root_prod[,1:nrow(sitedata)]
    
    
    ##Runs to max nyears of record for any lidet leaf at site i
    nyears <- max(sitedata$TIMES,na.rm=TRUE)
    
    ##Defines nspp for each site: how many leaves have yearly data from sitedata
    nspp <- as.numeric(length(which(!is.na(sitedata$TIMES))))
    
    
    ##Running the model
    MainFunction <- CORPSE_loop(nyears, nspp, litter_added, 
                                CORPSE_full_spinup_litter, CORPSE_full_spinup_bulk, 
                                CORPSE_full_spinup_rhizo, litter_production, 
                                froot_turnover, root_prod, litter_CN, root_CN, 
                                soil_T, soil_VWC)
    
    print(sites[i])
    
    
    ##Model Output
    
    # Read in LIDET observations to compare with model output
    dat_lidet_leaves <- read.csv("Input Files/LIDET_leaves.csv")
    #Lidet observed data for each site
    site_lidet <- subset(dat_lidet_leaves, site == sites[i])
    #Add a column for DOY
    site_lidet$DOY <- site_lidet$time*365
    
  
    ### Timeseries data for litterbag layer (unprotected pools only)
    # add spp as colnames
    colnames(MainFunction$litterbag_uFastC) <- t(sitedata[,1])
    colnames(MainFunction$litterbag_uSlowC) <- t(sitedata[,1])
    colnames(MainFunction$litterbag_uNecroC) <- t(sitedata[,1])
    colnames(MainFunction$litterbag_uFastN) <- t(sitedata[,1])
    colnames(MainFunction$litterbag_uSlowN) <- t(sitedata[,1])
    colnames(MainFunction$litterbag_uNecroN) <- t(sitedata[,1])
    colnames(MainFunction$litterbag_livingMicrobeC) <- t(sitedata[,1])
    colnames(MainFunction$litterbag_inorganicN) <- t(sitedata[,1])
    colnames(MainFunction$litterbag_CO2) <- t(sitedata[,1])
    colnames(MainFunction$litterbag_livingMicrobeN) <- t(sitedata[,1])
    
    
    # litterbag_uFastC
    write.csv(MainFunction$litterbag_uFastC, paste(dir, sites[i], "_" , paramset, "_", "litterbag_uFastC.csv",sep=""), row.names = FALSE)

    # litterbag_uSlowC
    write.csv(MainFunction$litterbag_uSlowC, paste(dir, sites[i], "_" , paramset, "_", "litterbag_uSlowC.csv",sep=""), row.names = FALSE)

    # litterbag_uNecroC
    write.csv(MainFunction$litterbag_uNecroC, paste(dir, sites[i], "_" , paramset, "_", "litterbag_uNecroC.csv",sep=""), row.names = FALSE)

    # litterbag_uFastN
    write.csv(MainFunction$litterbag_uFastN, paste(dir, sites[i], "_" , paramset, "_", "litterbag_uFastN.csv",sep=""), row.names = FALSE)

    # litterbag_uSlowN
    write.csv(MainFunction$litterbag_uSlowN, paste(dir, sites[i], "_" , paramset, "_", "litterbag_uSlowN.csv",sep=""), row.names = FALSE)

    # litterbag_uNecroN
    write.csv(MainFunction$litterbag_uNecroN, paste(dir, sites[i], "_" , paramset, "_", "litterbag_uNecroN.csv",sep=""), row.names = FALSE)

    # litterbag_livingMicrobeC
    write.csv(MainFunction$litterbag_livingMicrobeC, paste(dir, sites[i], "_" , paramset, "_", "litterbag_livingMicrobeC.csv",sep=""), row.names = FALSE)

    # litterbag_inorganicN
    write.csv(MainFunction$litterbag_inorganicN, paste(dir, sites[i], "_" , paramset, "_", "litterbag_inorganicN.csv",sep=""), row.names = FALSE)

    # litterbag_CO2
    write.csv(MainFunction$litterbag_CO2, paste(dir, sites[i], "_" , paramset, "_", "litterbag_CO2.csv",sep=""), row.names = FALSE)

    # litterbag_livingMicrobeN
    write.csv(MainFunction$litterbag_livingMicrobeN, paste(dir, sites[i], "_" , paramset, "_", "litterbag_livingMicrobeN.csv",sep=""), row.names = FALSE)

    
    #calculate total model C remaining in litterbag
    C_mod_kgm2_daily <- MainFunction$litterbag_uFastC + MainFunction$litterbag_uSlowC + 
      MainFunction$litterbag_uNecroC + MainFunction$litterbag_livingMicrobeC
    
    #sum initial litterbag C pools to add as first row to C_mod_kgm2_daily
    C_top <- t(litter_added$uFastC + litter_added$uSLowC + litter_added$livingMicrobeC)
    
    # add initial C values as first row in C_mod_kgm2_daily
    C_mod_kgm2_daily <- rbind(C_top, C_mod_kgm2_daily)
    
    #add spp as colnames
    colnames(C_mod_kgm2_daily) <- t(sitedata[,1])
    
    # calculate percent C remaining at each time step based on initial C in leaves 
    # initial C is very slightly greater than 1 because of microbial biomass C
    # (because model output is in kg/m2 and we are analyzing percent remaining as in the LIDET field data)
    # make a matrix to fill
    C_mod_pct_daily <- matrix(NA, nrow = nrow(C_mod_kgm2_daily), ncol = nspp)
    # name columns for each species
    colnames(C_mod_pct_daily) <- t(sitedata[,1])
    # fill in the matrix with percent C remaining at each time point
    for (l in 1:nrow(C_mod_pct_daily)) {
      C_mod_pct_daily[l,] <- C_mod_kgm2_daily[l,]/C_mod_kgm2_daily[1,]
    }
    
    #make a matrix with the day of each LIDET sampling time for each species (cumulative day for entire LIDET experiment)
    #matrix will be used to pull out model data from the same cumulative day as the field data
    modeldoy <- matrix(NA, max(site_lidet$DOY/365), nrow(sitedata))
    for (m in 1:nrow(sitedata)){
      spp1 <- subset(site_lidet, spp == as.character(sitedata[m,1]))
      modeldoy[spp1$time,m] <- spp1$DOY
    }
    
    #YEARLY data: Make a matrix of C remaining data (kg/m2) for each year of LIDET sampling
    C_mod_kgm2_yearly <- matrix(NA, max(site_lidet$DOY/365), nrow(sitedata))
    #pull out the day of model data from the daily data table to match the yearly LIDET data sampling day 
    for (m in 1:round(max(site_lidet$DOY)/365)){
      for (n in 1:nrow(sitedata)){
        C_mod_kgm2_yearly[m,n] <- C_mod_kgm2_daily[floor(modeldoy[m,n]),n]
      }    }
    # name columns for each species
    colnames(C_mod_kgm2_yearly) <- t(sitedata[,1])
    
    
    #YEARLY data: Make a matrix of C remaining data (% remaining) for each year of LIDET sampling
    C_mod_pct_yearly <- matrix(NA, max(site_lidet$DOY/365), nrow(sitedata))
    #pull out the day of model data from the daily data table to match the yearly LIDET data sampling day 
    for (m in 1:round(max(site_lidet$DOY)/365)){
      for (n in 1:nrow(sitedata)){
        C_mod_pct_yearly[m,n] <- C_mod_pct_daily[floor(modeldoy[m,n]),n]
      }    }
    # name columns for each species
    colnames(C_mod_pct_yearly) <- t(sitedata[,1])
    
    
    #calculate total model N remaining
    N_mod_kgm2_daily <- MainFunction$litterbag_uFastN + MainFunction$litterbag_uSlowN + 
      MainFunction$litterbag_uNecroN + MainFunction$litterbag_livingMicrobeN
    
    #sum initial litterbag N pools to add as first row to N_mod_kgm2_daily, transpose to row
    N_top <- t(litter_added$uFastN + litter_added$uSlowN + litter_added$livingMicrobeN)
    
    # add initial N values as first row in N_mod_kgm2_daily
    N_mod_kgm2_daily <- rbind(N_top, N_mod_kgm2_daily)
    
    #add spp as colnames
    colnames(N_mod_kgm2_daily) <- t(sitedata[,1])
    
    # calculate percent N remaining at each time step based on initial N in leaves 
    # initial N was calculated based on C:N of leaves so it varies widely. Initial total C was always 1.0001
    # (because model output is in kg/m2 and we are analyzing percent remaining as in the LIDET field data)
    # make a matrix to fill
    N_mod_pct_daily <- matrix(NA, nrow = nrow(N_mod_kgm2_daily), ncol = nspp)
    # name columns for each species
    colnames(N_mod_pct_daily) <- t(sitedata[,1])
    # fill in the matrix with percent N remaining at each time point
    for (l in 1:nrow(N_mod_pct_daily)) {
      N_mod_pct_daily[l,] <- N_mod_kgm2_daily[l,]/N_mod_kgm2_daily[1,]
    }
    
    #YEARLY data: Make a matrix of N remaining data (kg/m2) for each year of LIDET sampling
    #use the "modeldoy" matrix that was made for the C data above to pull out model data at each point in lidet data
    #modeldoy is used to pull out model data from the same cumulative day as the field data
    N_mod_kgm2_yearly <- matrix(NA, max(site_lidet$DOY/365), nrow(sitedata))
    # name columns for each species
    colnames(N_mod_kgm2_yearly) <- t(sitedata[,1])
    #pull out the day of model data from the daily data table to match the yearly LIDET data sampling day 
    for (m in 1:round(max(site_lidet$DOY)/365)){
      for (n in 1:nrow(sitedata)){
        N_mod_kgm2_yearly[m,n] <- N_mod_kgm2_daily[floor(modeldoy[m,n]),n]
      }    }
    
    #YEARLY data: Make a matrix of N remaining data (% remaining) for each year of LIDET sampling
    N_mod_pct_yearly <- matrix(NA, max(site_lidet$DOY/365), nrow(sitedata))
    # name columns for each species
    colnames(N_mod_pct_yearly) <- t(sitedata[,1])
    #pull out the day of model data from the daily data table to match the yearly LIDET data sampling day 
    for (m in 1:round(max(site_lidet$DOY)/365)){
      for (n in 1:nrow(sitedata)){
        N_mod_pct_yearly[m,n] <- N_mod_pct_daily[floor(modeldoy[m,n]),n]
      }    }
    
    
    # Save C Data
    # model yearly C data in kg/m2
    write.csv(C_mod_kgm2_yearly,
             paste(dir, sites [i], "_", paramset, "_C_mod_kgm2_yearly.csv", sep=""))
  
    # model yearly C data in % remaining
    write.csv(C_mod_pct_yearly,
              paste(dir, sites[i], "_", paramset, "_C_mod_pct_yearly.csv", sep=""))
    
    # model daily C data in kg/m2
    #add a day of year column, nyr column, and site column 
    exportC_day_kgm2 <- as.data.frame(C_mod_kgm2_daily)
    exportC_day_kgm2$cumday <- 0:(nyears*365)
    exportC_day_kgm2$nyr <- c(0, rep(1:nyears, each = 365))
    exportC_day_kgm2$doy <- c(0, rep_len(1:365, length.out=nrow(C_mod_kgm2_daily)-1))
    exportC_day_kgm2$site <- sites[i]
    write.csv(exportC_day_kgm2, 
              paste(dir, sites[i], "_", paramset, "_C_mod_kgm2_daily.csv",sep=""),
              row.names = FALSE) #daily model data
    
    # model daily C data in % remaining
    #add a day of year column, nyr column, and site column 
    exportC_day_pct <- as.data.frame(C_mod_pct_daily)
    exportC_day_pct$cumday <- 0:(nyears*365)
    exportC_day_pct$nyr <- c(0, rep(1:nyears, each = 365))
    exportC_day_pct$doy <- c(0, rep_len(1:365, length.out=nrow(C_mod_pct_daily)-1))
    exportC_day_pct$site <- sites[i]
    write.csv(exportC_day_pct, 
              paste(dir, sites[i], "_", paramset, "_C_mod_pct_daily.csv", sep=""),
              row.names = FALSE) #daily model data
  
    
    # Save N Data
    # model yearly N data in kg/m2
    write.csv(N_mod_kgm2_yearly,
              paste(dir, sites[i], "_", paramset, "_N_mod_kgm2_yearly.csv",sep=""))
    
    # model yearly N data in % remaining
    write.csv(N_mod_pct_yearly,
              paste(dir, sites[i], "_", paramset, "_N_mod_pct_yearly.csv",sep=""))
    
    # model daily N data in kg/m2
    #add a day of year column, nyr column, and site column 
    exportN_day_kgm2 <- as.data.frame(N_mod_kgm2_daily)
    exportN_day_kgm2$cumday <- 0:(nyears*365)
    exportN_day_kgm2$nyr <- c(0, rep(1:nyears, each = 365))
    exportN_day_kgm2$doy <- c(0, rep_len(1:365, length.out=nrow(N_mod_kgm2_daily)-1))
    exportN_day_kgm2$site <- sites[i]
    write.csv(exportN_day_kgm2, 
              paste(dir, sites[i], "_", paramset, "_N_mod_kgm2_daily.csv", sep = ""),
              row.names = FALSE) #daily model data
    
    # model daily N data in % remaining
    #add a day of year column, nyr column, and site column 
    exportN_day_pct <- as.data.frame(N_mod_pct_daily)
    exportN_day_pct$cumday <- 0:(nyears*365)
    exportN_day_pct$nyr <- c(0, rep(1:nyears, each = 365))
    exportN_day_pct$doy <- c(0, rep_len(1:365, length.out=nrow(N_mod_pct_daily)-1))
    exportN_day_pct$site <- sites[i]
    write.csv(exportN_day_pct, 
              paste(dir, sites[i], "_", paramset, "_N_mod_pct_daily.csv", sep = ""),
              row.names = FALSE) #daily model data
    
    
  }
  