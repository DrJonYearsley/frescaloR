# A script to calculate the alpha values for Frescalo 
# The estimation of trends in species occurence has not been implemented
#
# This has been checked against Frescalo output and it agrees
#
# Written by Jon Yearsley (jon.yearsley@ucd.ie) 2nd Sept 2016

# Modified:
# 6/9/2016 (JY): Data sent to cluster nodes in chunks (variable=chunkSize)
# 8/9/2012 (JY): Save species frequency file and include expected species richness in output

#setwd('/home/jon/MEGA/Jack/Frescalo')
rm(list=ls())  # Remove all variables from the memory

require(foreach, quietly=T)
require(doParallel, quietly=T)
source('frescalo_functions.R')

# Set the size of the cluster (number of nodes). 
# cores=NULL will automatically pick the number of cores per node 
# (by defauult this is half the number of available cores)
cl <- makeCluster(2)
registerDoParallel(cl, cores=NULL)

R_star = 0.27
trend_analysis = FALSE

# weight.file = './bsbi_hectads_2000_2006_2010_weights.txt'
# species.file = './bsbi_hectads_2000_2006_2010_sample.txt'

weight.file = './TestData/weights.txt'
species.file = './TestData/Test.txt'

#outputFilename = './bsbi_hectads_frescalo_out.txt'
outputPrefix = './test'

Phi = 0.74    # The standardisation for recorder effort
chunkSize = 5 # Number of hectads to pass to each compute node

# Import data
d = read.table(weight.file, header=F, col.names=c('location1','location2','w','w1','w2','nsim','ndist'), stringsAsFactors=F)
s = read.table(species.file, header=F, col.names=c('location','species','Year'), stringsAsFactors=F)


##############################################################
spLocations = unique(s$location)
speciesNames = as.character(unique(s$species))   # Create list of unique species
# The variable species could be made numerical to make it more efficient

# For each region record presence/absence of each species (a_ij in Hill 2011)
locationGroups = as.factor(rep(c(1:ceiling(length(spLocations)/chunkSize)),each=chunkSize))
sSplit = split(s, locationGroups[match(s$location, spLocations)])  # Split species data up into hectads

idx = iter(sSplit)
speciesList <- foreach(spList = idx, .inorder=T, .combine='c') %dopar% {
  speciesListFun(spList, speciesNames)
}
# Add an additional species list where everything is absent 
speciesList[[length(speciesList)+1]] = rep(0, times=length(speciesNames))
spLocations = c(spLocations,'null_location')

#################################################################
# For each focal regional calculate the sampling effort multipler
dSub = d[,1:3]
location1List = as.character(unique(dSub$location1))    # Create list of unique focal regions
location1Groups = as.factor(rep(c(1:ceiling(length(location1List)/chunkSize)),each=chunkSize))

dSplit = split(dSub, location1Groups[match(dSub$location1, location1List)])  # Split data up into focal regions
idx2 = iter(dSplit)
output <- foreach(focalData = idx2, .inorder=T, .combine='cfun', .multicombine=TRUE) %dopar% {
  frescalo(focalData, speciesList, spLocations, speciesNames, Phi, R_star=0.27, missing=2)
}

#################################################################
# Trend analysis (still to be debugged)
if (trend_analysis) {
  # Do the Frescalo trend analysis if there are more than 1 year bins (use same location groups as sSplit)
  sSplit2 = split(s, as.factor(s$time))  # Split species data up into year bins
  idx3 = iter(sSplit2)
  trend.out <- foreach(sData=idx3, .inorder=T, .combine='rbind') %do% {
    trend(sData, output$freq.out)
  }
}

################################################################
# Write the output to a text file
write.table(format(output$frescalo.out, digits=4,zero.print=T, width=10), 
            file=paste(outputPrefix,'_frescalo_out.txt',sep=''), col.names=T, row.names=F, quote=F, sep=' ')

write.table(format(output$freq.out, digits=4, zero.print=T, width=10, scientific=F, justify='left'), 
            file=paste(outputPrefix,'_frescalo_freq.txt',sep=''), col.names=T, row.names=F, quote=F, sep=' ')

stopCluster(cl)
gc()

source('../R_Scripts/os2eastnorth.R')

en=os2eastnorth(GR=output$frescalo.out$location, hectad=T)$en
output$frescalo.out$x = en[,1]
output$frescalo.out$y = en[,2]

library(ggplot2)

ggplot(data=output$frescalo.out, aes(x=x, y=y, colour=alpha)) + geom_point()
