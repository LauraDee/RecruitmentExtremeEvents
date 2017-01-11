### Derive temperature time series by LME
# Laura Dee May 20, 2016 (cleaned)
## Per year and LME, this code computes mean SST, min SST, within-year SST SD, max SST, within-year SST CV, Within year mean SST minus long term mean SST 
# Within year mean SST minus within year min SST, Within year max SST minus within year mean SST 
# using NOAA Optimum Interpolation (OI) SST V2 data 1982-2015 (monthly mean SST by 1x1 degree cell)
## NOAA data from: http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html

rm(list=ls())
#install.packages(c("ncdf", "rgdal","bmp", "lubridate"))
install.packages("ncdf4")
#library(ncdf)
library(ncdf4)
library(bmp)
library(rgdal)
library(lubridate)
require(ggplot2)

#setwd("~/Desktop/temperature data/")
setwd("~/Dropbox/Temp var and recruitment/code")

## load files
lmeCSV = read.csv("LMES.csv", header=T)
lmes = lmeCSV[,2]
#read in lme mask from Malin as 180x360 matrix. indices: [lat,long]. Stored N-->S, W-->E
# so lmetxt[1,1] is northwest corner of map, lmetxt[180,360] is SE corner
lmetxt = as.matrix(read.csv("lmemask_2015-04-08.txt", sep=" ", header=F))
##  Inport list of LMEs: Malin's LMEs and SAUP LMEs
AllLMEs <- read.csv("MalinANDsaupLMEs.csv")

lats = seq(89.5,-89.5,-1) 
lons = seq(0.5,359.5,1)
nlon = length(lons) 
nlat = length(lats)

## Monthly mean SST by 1x1 degree cell with NOAA Optimum Interpolation (OI) SST V2 data
# for ncdf4 package 
nc = nc_open("sst.mnmean1982_2015.nc")

### Long-term means
#nc = nc_open("sst.ltm.1961-1990.nc")

ncdates = nc$dim$time$vals
ncdates = as.Date(ncdates,origin = '1800-1-1')

latIndex = 1
lonIndex = 1
timeIndex = 1

#to get Malin's mask ready for use to clip SST data, we have to transform it
#1. transform from 0 to 360E vs going from 180W to 180E
malinMask = cbind(lmetxt[,181:360],lmetxt[,1:180]) #sst data goes from 0 to 360E
#2. make longitude first index
malinMask = t(malinMask) #sst data has longitude as first index rather than latitude
## take a look at the LMEs and also ones Malin created:
lmeColors = sample(colours(), 70)
image(malinMask[,180:1],col=lmeColors)

####### this is just a mask with the LMEs from SAUP ########
##### this specifies which cells belong to which LME to extract data from:
#       lmeTIF = readGDAL("lmes1deg.tif")
#       rawLMEMask = matrix(lmeTIF@data$band1, nrow=180, ncol=360, byrow=TRUE)
#       rawLMEMask[rawLMEMask==0] = NA
#       rawLMEMask = cbind(rawLMEMask[,181:360],rawLMEMask[,1:180]) #sst data goes from 0 to 360E
#       rawLMEMask = t(rawLMEMask) #sst data has longitude as first index rather than latitude
#       image(rawLMEMask[,180:1], col=lmeColors)


## Create matrices to store processed data by LME and time
meanSSTs = matrix(NA, nrow=length(lmes), ncol=length(ncdates))
minSSTs = matrix(NA, nrow=length(lmes), ncol=length(ncdates))
maxSSTs = matrix(NA, nrow=length(lmes), ncol=length(ncdates))
sdSSTs = matrix(NA, nrow=length(lmes), ncol=length(ncdates))


# New loop! 
# Compute a set of SST metrics per LME. These will be computed temporally per
# cell first, then spatially averaged, vs old way of temporal summary 
# statistics computed on spatial averages (we do it both ways at the end)

#First, get all the data into a usable 3D array: [lon, lat, time]
#Note we only do ncdates 2:397 since that corresponds to 1982-2014 inclusive (starting at lon=1, lat=1, Jan 1982 
# (so say 2 because the dataset starts in Dec 1981, and we want the 2nd entry, Jan 1982. ))

## for the ncdf4 package:
sstall = ncvar_get(nc, varid = 'sst', start = c(1,1,2),
                   count = c(nlon,nlat,396))   
    # with ncdf package: 
      # sstall = get.var.ncdf(nc, varid = 'sst', start = c(1,1,2),count = c(nlon,nlat,396))   

#convert to Kelvin
sstallKelvin = sstall + 273.15

#store in 4D array so that within year computations are simpler.
#Dimensions are now lon, lat, year, month
NYEARS = 33 # 1982-2014
sstallKelvinByMonth = array(dim = c(360,180,NYEARS,12))
for(year in 1:NYEARS) {
  sstallKelvinByMonth[,,year,] = sstallKelvin[,,((year-1)*12+1):(year*12)]
}

#compute summary stats per cell
#1. Within year CV of SST - this will still be what we interact with FD
#2. Within year mean SST minus long term mean SST 
#3. Within year mean SST minus within year min SST
#4. Within year max SST minus within year mean SST 
#5. Within year sd SST 

#### need to store monthly temps by year and LME so then I can do histogram of the first 10 years VS last 10 years 

#convenience function for CV
coefVar = function(x) {
  return(sd(x)/mean(x))
}

#Compute within-year stats for SST per cell
#Dimensions are 1) lon, 2) lat, 3) year, 4) month  -- hold 1,2,3 constant and look at CV, mean, max, etc across MONTHS
cvstats = apply(sstallKelvinByMonth, c(1,2,3), FUN=coefVar)
meanstats = apply(sstallKelvinByMonth, c(1,2,3), FUN=mean)
maxstats = apply(sstallKelvinByMonth, c(1,2,3), FUN=max)
minstats = apply(sstallKelvinByMonth, c(1,2,3), FUN=min)
sdstats = apply(sstallKelvinByMonth, c(1,2,3), FUN=sd)

#Compute long term mean SST per cell
longtermmeans = apply(sstallKelvin, c(1,2), FUN=mean)


#############################################################################
####Compute per-LME temperature metric -- masking temp metrics by LME #######
############################################################################
tempDistData = NULL

getMaskedData = function(data, mask) {
  return(data*mask)
}

getMaskedDemeanedData = function(data, mask, means) {
  return((data-means)*mask)
}

###### MODIFY THIS PART TO DECIDE WHICH LMEs TO USE FOR THE MASK #################
#LD: need to modify our list of LMEs to use to include Malin's synthetic LMEs
## import file with SAUP lmes and Malin's added ones 
AllLMEs <- read.csv("MalinANDsaupLMEs.csv")
head(AllLMEs)
lmes <- AllLMEs[,2]

lmesToUse = lmes # to look at all lmes
nLMEs = length(lmesToUse)

maskToUse = malinMask

##############################################
#### Computing other SST metrics ############# DID THIS ON MAY 20, 2016
##############################################

maskToUse=malinMask
# to use Malin's list of LMEs
lmes <- AllLMEs[,2]
lmesToUse = lmes # to look at all lmes
nLMEs = length(lmesToUse)
#lmesToUse = c(23) # to look only at lme 23 or a particular LME

#STATS OF INTEREST PER CELL:
#1. Within year CV 
#2. Within year mean SST minus long term mean SST
meanDeviations = array(dim=c(360,180,NYEARS))
for(year in 1:25) {
  meanDeviations[,,year] =  meanstats[,,year] - longtermmeans
}
#3. within year mean minus within year min
meanMinusMin = meanstats - minstats
#4. within year max minus within year mean
maxMinusMean = maxstats - meanstats

# Now the per-cell stats we care about are in variables
# cvstats, meanDeviations, meanMinusMin, maxMinusMin.
# We want to go through each LME, mask the cells to those in that
# LME, and compute spatial averages of each per-cell stat across
# all cells in the LME, so we get one number for each stat for
# each LME. So we should get four stats per LME x year.

# Go through each LME, mask relevant stats, take spatial means,
# and store results in a data frame
allLMEsSSTData = NULL
allLMEsMonthlyTemps = NULL

# Compute a mean of the values in data according to the
# mask defined by mask. The mask should have 1s where values
# in data are valid, and NA otherwise.
maskedMean <- function(data, mask) {
  return(mean(data*mask, na.rm=T))
}

maskedMAx <- function(data, mask) {
  return(max(data*mask, na.rm=))
}

allMonthlyMeansDF = NULL

#process LMEs one by one
for(lmeIndex in 1:length(lmes)) {
  #generate an LME mask
  lmeId = as.numeric(lmes[lmeIndex])
  mask = maskToUse
  mask[mask!=lmeId] = NA
  mask[mask==lmeId] = 1
  
  # spatial mean first, then annual mean
  # hold year and month (dims 3 and 4) fixed, average (masked) across space
  monthlyMeanForLME = apply(sstallKelvinByMonth, c(3,4), maskedMean, mask)
  monthlyMeanDF = as.data.frame(monthlyMeanForLME)
  names(monthlyMeanDF) = paste("temp", 1:12, sep=".")
  monthlyMeanDF$LME = lmeId
  monthlyMeanDFLong = reshape(monthlyMeanDF, dir="long", varying = 1:12)
  names(monthlyMeanDFLong)[2] = "month"
  names(monthlyMeanDFLong)[4] = "year"
  monthlyMeanDFLong$year = monthlyMeanDFLong$year + 1981
  
  # tack LME-specific monthly mean data onto running data frame
  allMonthlyMeansDF = rbind(allMonthlyMeansDF, monthlyMeanDFLong)
  
  # average the monthly mean across months within a year, holding year (dim 1) fixed
  annualMeanForLME = apply(monthlyMeanForLME,1,mean)
  
  #take spatial means per stat across all cells in that LME
  maskedAnnualCV = apply(cvstats, 3, maskedMean, mask)
  maskedAnnualMeanDeviation = apply(meanDeviations, 3, maskedMean, mask)
  maskedAnnualMeanMinusMin = apply(meanMinusMin, 3, maskedMean, mask)
  maskedAnnualMaxMinusMean = apply(maxMinusMean, 3, maskedMean, mask)
  maskedAnnualSD = apply(sdstats, 3, maskedMean, mask)
  maskedAnnualMin = apply(minstats, 3, maskedMean, mask)
  maskedAnnualMax = apply(maxstats, 3, maskedMean, mask)
  
  ####### Find the annual max and mean per LME per year
  maskedAnnualMax = apply(maxstats, 3, maskedMean, mask)
  maskedAnnualMean = apply(meanstats, 3, maskedMean, mask)
  
  #put the results in a data frame, with four data points per year x LME
  lmeDF = data.frame(lmeId = lmeId, 
                     year = 1982:2014,
                     meanSST =  maskedAnnualMean,
                     meanSSTCV = maskedAnnualCV, 
                     maxSST = maskedAnnualMax,
                     sstMeanDev = maskedAnnualMeanDeviation,
                     minSST = maskedAnnualMin,
                     maxSST = maskedAnnualMax,
                     sstMeanMinusMin = maskedAnnualMeanMinusMin,
                     sstMaxMinusMean = maskedAnnualMaxMinusMean,
                     sstSD = maskedAnnualSD)
  
  #combine with other LME info
  allLMEsSSTData = rbind(allLMEsSSTData, lmeDF)
  
}

##data aggregrated by time (by cell), then space
write.csv(allLMEsSSTData, "sstTSspatialmeanofcells.csv", row.names=F)
## data aggregated first by space (monthly per LME), then taking annual average of months
write.csv(allMonthlyMeansDF, "sstByLMEandMonth.csv", row.names = F)


#############################################################################################################
#### compute other metrics from the spatially average, then time average data (dataframe = allMonthlyMeansDF) ##
#############################################################################################################

aggregate(temp$allMonthlyMeansDF, by = c(LME$allMonthlyMeansDF, year$allMonthlyMeansDF), FUN = min, data = allMonthlyMeansDF)
head(allMonthlyMeansDF)

ggplot(data = allMonthlyMeansDF, aes(x = month, y = temp, by = LME)) + 
  geom_point() + 
  facet_wrap(~ LME, scale = "free")

