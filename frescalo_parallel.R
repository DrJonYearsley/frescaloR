# A script to calculate the alpha values for Frescalo 
# The estimation of trends in species occurence has not been implemented
#
# This has been checked against Frescalo output and it agrees
#
# Written by Jon Yearsley (jon.yearsley@ucd.ie) 2nd Sept 2016

#setwd('/home/jon/MEGA/Jack/Frescalo')
rm(list=ls())  # Remove all variables from the memory

require(foreach, quietly=T)
require(doParallel, quietly=T)

# Set the size of the cluster (number of nodes). 
# cores=NULL will automatically pick the number of cores per node 
# (by defauult this is half the number of available cores)
cl <- makeCluster(1)
registerDoParallel(cl, cores=NULL)


# weight.file = './bsbi_2009_weight.txt'
# species.file = './bsbi_2009_sample.txt'

weight.file = './TestData/weights.txt'
species.file = './TestData/Test.txt'

outputFilename = './frescalo_out.txt'

Phi = 0.74    # The standardisation for recorder effort

# Import data
d = read.table(weight.file, header=F, col.names=c('Target','Hectad','w','w1','w2','nsim','ndist'), stringsAsFactors=F)
s = read.table(species.file, header=F, col.name=c('Hectad','Species','Year'), stringsAsFactors=F)


##########################
# Function definitions

min.fun <- function(alpha,fij,Phi) {
  # The function to minimise to fit Frescalo model to the data 
  # Page 5 column 2 in Hill (2011)
  F=1-(1-fij)^alpha
  F[abs(fij-1)<1e-10] = 1
  return( sum(F^2)/sum(F) - Phi)
}

frescalo = function(focal_d, speciesList, spLocations, Phi, missing = 1) {
# This function calculates the sampling effort multiplier, alpha, 
# that would equate sampling effort across all regions
# This is the main method of Frescalo.
  
  # Argument missing defines what to do if data in the neighbourhood are missing.
  # missing = 1 (default) do not proceed with the calculation 
  # missing = 2 set missing species data to all absences and proceed
  
  focal = focal_d$Target[1]
  
  # Identify species in neighbourhood of focal region
  # If data are mssing assign it to the empty species list (last element of speciesList)
  neighbourhood = match(focal_d$Hectad, spLocations, nomatch=length(speciesList))
  speciesRegional = speciesList[neighbourhood]
  
  missingData = neighbourhood==length(speciesList)
  if (any(missingData) & missing==1) {
    # Species data missing from a location inthe neighbourhood. Ignore this focal location
    warning(paste('Removing location ',focal,'. Missing data for locations in the neighbourhood.', sep=''))

    sol = list(root = NA, iter=NA)
    phi_in = NA
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
  }
  return(data.frame(location=focal, alpha=sol$root, phi_in=phi_in, iter=sol$iter))
}

#######################################################

grouping = as.factor(s$Hectad)
spLocations = levels(grouping)

species = as.character(unique(s$Species))   # Create list of unique species

# Limit focal regions to ones where we have species data
#dSub = subset(d, Target %in% spLocations)
dSub = d[,1:3]
targets = as.character(unique(dSub$Target))    # Create list of unique focal regions

# For each region record presence/absence of each species (a_ij in Hill 2011)
sSplit = split(s, grouping)  # Split species data up into hectads

idx = iter(sSplit)
speciesList <- foreach(spList = idx, .inorder=T) %do% {
  as.integer(species %in% spList$Species)
}
# Add an additional species list where everything is absent 
speciesList[[length(speciesList)+1]] = rep(0, times=length(species))

# For each focal regional calculate the sampling effort multipler
dSplit = split(dSub[,1:3], as.factor(dSub$Target))  # Split data up into focal regions
idx2 = iter(dSplit)
frescalo.out <- foreach(focal_data = idx2, .inorder=T, .combine='rbind') %dopar% {
  frescalo(focal_data, speciesList, spLocations, Phi, missing=2)
}

# Write the output to a text file
write.table(format(frescalo.out, digits=4,zero.print=T, width=10), file=outputFilename, col.names=T, row.names=F, quote=F, sep=' ')

stopCluster(cl)
gc()
