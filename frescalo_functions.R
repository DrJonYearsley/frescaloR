##########################
# Function definitions for frescalo_parallel R Script
#
# This script is sourced in the frescalo_parallel.R script

min.fun <- function(alpha,fij,Phi) {
  # The function to minimise to fit Frescalo model to the data 
  # Page 5 column 2 in Hill (2011)
  F=1-(1-fij)^alpha
  F[abs(fij-1)<1e-10] = 1
  return( sum(F^2)/sum(F) - Phi)
}

speciesListFun = function(spList, species) {
  # A function that returns the presences and absecences of species per region
  # spList should be a data frame with columns location, species and time  

  # Create list of all locations
  locationList = unique(spList$location)
  
  # Initialize output list
  out = lapply(c(1:length(locationList)), function(i) {i})
  
  for (i in 1:length(locationList)) {
    spListSub = subset(spList, location==locationList[i])
    out[[i]] = as.integer(species %in% spListSub$species)
  }
  return(out)
}

frescalo = function(data_in, speciesList, spLocations, speciesNames, Phi=0.74, R_star=0.27, missing.data = 2) {
  # This function calculates the sampling effort multiplier, alpha, 
  # that would equate sampling effort across all regions
  # This is the main method of Frescalo.
  
  # The script returns a list containing two data frames: frescalo.out, freq.out
  # frescalo.out contains the recorder-effort multiplier for each focal location, and other related information
  # freq.out contains the local frequency of each species in each nehgbourhood of a focal region, along with other related information
  # Species identified as a benchmark in freq.out (benchmark=1) can be used to assess trends
  
  # Argument missing.data defines what to do if data in the neighbourhood are missing.
  # missing.data = 1 (default) do not proceed with the calculation 
  # missing.data = 2 set missing species data to all absences and proceed
  
  # R_star defines the corrected frequency threshold for benchmark species 
  # (default = 0.27 means the top 27% of species)
  
  location1List = unique(data_in$location1)
  nSpecies = sapply(speciesList[match(location1List,spLocations,nomatch=length(spLocations))], FUN=sum, simplify=T)
  out = data.frame(location=location1List, 
                   nSpecies=nSpecies, 
                   phi_in=NA, alpha=NA, phi_out=Phi, spnum_in=NA, spnum_out=NA, iter=NA)
  
  freq.out = data.frame(location=c(), species=c(), pres=c(), freq=c(), rank=c(), benchmark=c())
  
  for (f in 1:length(location1List)) {
    focal = location1List[f]
    focal_d = subset(data_in, location1==location1List[f])
    
    # Identify species in neighbourhood of focal region
    # If data are mssing assign it to the empty species list (last element of speciesList)
    neighbourhood = match(focal_d$location2, spLocations, nomatch=length(speciesList))
    speciesRegional = speciesList[neighbourhood]
    
    missingData = neighbourhood==length(speciesList)
    if (any(missingData) & missing.data==1) {
      # Species data missing from a location in the neighbourhood. Ignore this focal location
      warning(paste('Removing location ',focal,'. Missing data for locations in the neighbourhood.', sep=''))
      
      out$alpha[f] = NA
      out$iter[f] = NA
      out$phi_in[f]= NA 
    } else {
      # Calculate weights of locations in the neighbourhood
      weights = focal_d$w/sum(focal_d$w)
      
      # Create weighted neighbourhood frequencies (checked against Frescalo)
      frequency = Reduce('+',Map('*',as.list(weights), speciesRegional))  
      
      phi_in = sum(frequency^2) / sum(frequency)
      
      # Calculate the multiplier (alpha) that equalises recording effort 
      alpha.min = 1  # Minimum alpha ( =1 means no correction required)
      alpha.max = 5
      # Increase alpha.max until min.fun() becomes positive (i.e. ensure there is a zero)
      while (min.fun(alpha.max, frequency, Phi)<0) { alpha.max = alpha.max + 5}  
      while (min.fun(alpha.min, frequency, Phi)>0) { alpha.min = alpha.min/2}  
      
      # Find sampling-effort multiplier
      sol=uniroot(min.fun,interval=c(alpha.min,alpha.max), tol=0.001, frequency, Phi)
      out$alpha[f] = sol$root
      out$phi_in[f] = phi_in
      out$iter[f] = sol$iter
      # Expected species richness before recorder effort correct
      out$spnum_in[f] = sum(frequency)   
      # Expected species richness after recorder effort correct
      out$spnum_out[f] = sum(1-exp(sol$root*log(1-frequency)))
      
      # Create the data frame with the local species frequencies after correction
      freq.order = order(frequency, decreasing=T)
      focal.ind = match(focal,spLocations,nomatch=length(spLocations))
      
      # Pick out benchmark species (assumes that species are ordered by rank)
      benchmarkSpecies = rep(0, times=length(frequency))
      benchmarkSpecies[1:floor(R_star*length(frequency))] = 1
      
      freq.out = rbind(freq.out, data.frame(location=focal, 
                                            species=speciesNames[freq.order], 
                                            pres=speciesList[[focal.ind]][freq.order],
                                            freq=frequency[freq.order],
                                            rank=c(1:length(speciesList[[focal.ind]])),
                                            benchmark=benchmarkSpecies))
    }
  }
  return(list(frescalo.out=out, freq.out=freq.out))
}


trend = function(s_data, freq.out) {
  # Function to calculate the scaling factor for each species in one year bin
  
  timeBin = unique(s_data$time)
  if (length(timeBin)>1) {
    warning('More than one time bin supplied to trend() function')
  }
  
  locationList = unique(freq.out$location)
  spList = sort(unique(s_data$species))
  
  s_it = rep(0, times=length(locationList))
  # Calculate the proportion of benchmark species in each hectad (for this time bin)
  for (f in 1:length(locationList)) {
    focal_s = subset(s_data, location==locationList[f],select='species')
    focal_bench = subset(freq.out, location==locationList[f] & benchmark==1, select='species')
    
    s_it[f] = sum(focal_bench$species %in% focal_s$species) / length(focal_bench$species)
  }
  
  x = rep(NA, times=length(spList))   # Vector to contain the time factors
  for (s in 1:length(spList)) {
    focal_f = subset(freq.out, species==spList[s])
    
    # Rescale frequencies by effort for all hectads
    sf = focal_f$freq[match(locationList,focal_f$location)]*s_it
    
    # Calculate Q_ijt   
    # P_ijt = 1-exp(-Q_ijt x_jt) where x_jt is the time factor for species j at time t
    # and P_ijt is prob of observing species j in hectad i at time t
    Q_ijt = rep(-log(1-0.98), each=length(locationList))
    Q_ijt[sf<0.98] = -log(1-sf[sf<0.98])
    
    # Select a x max that ensire a sign change in min_trend_fun
    xMax = 5
    while(min_trend_fun(xMax, Q_ijt)<=0) {xMax = xMax+5}
    sol = uniroot(min_trend_fun,interval=c(0,xMax), tol=0.001, Q_ijt)
    x[s] = sol$root
  }
  
  return(data.frame(species=spList, time=timeBin, tFactor=x))
}


cfun = function(...) {
  # Bespoke function to combine the output from the frescalo() function
  input_list <- list(...)
  return(list(frescalo.out=Reduce('rbind',Map(function(x){x[[1]]},input_list)), 
              freq.out=Reduce('rbind',Map(function(x){x[[2]]},input_list))))
}


min_trend_fun = function(x,Q) {
  # A function to solve sum_i P_ijt = sum_i Q_ijt  in the Frescalo trend analysis Hill (2012)
  # P_ijt = 1-exp(-Q_ijt*x)
  
  minusQ = (1-Q)
  minusQ[minusQ<1e-12] = 1e-12
  return(sum(Q*x+log(minusQ)))
}
#######################################################
