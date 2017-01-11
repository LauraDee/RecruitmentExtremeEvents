## Process RAM data to derive Recruitment time series 
# Laura Dee Jan 3, 2017
## Matching RAM data to temp data at LME scale

## load packages for analyses ####
require(sandwich)
require(car)
require(plm)
require(ggplot2)
source("stdErrHelperFunctions.R")
## clustered robust standard errors, clustered on stock id
source("stdErrHelperFunctionsClusterStockid.R")  ###**LD TO DOUBLE CHECK THIS **

######### Read in different dataset files ############
setwd("~/Dropbox/Extreme Events and Fisheries/temperature data")

# full RAM dataset v 2.5
RamStocksData <- read.csv("RAMtimeseries.csv")
nrow(unique(RamStocksData[,c("assessid","stockid","year")]))
dim(RamStocksData)

# Matching RAM stocks to LMEs 
RamStocksByLME <- read.csv("RAMstocks2012_matchedtoLME.csv")

# Malin's life history database (ref for this??)
malinLHfile <- read.csv("SummarySheetFinalMissing_2015-04-21_forLauraDee.csv")

## temperature EE bin data masked by Malin's LMEs
#**LD note to self: confirm these files and the code that produced them are most recent ** #
nonDetrendedBelow10 = read.csv("malinLME_NonDT_EEbelow10.csv", header=T)
nonDetrendedAbove90 = read.csv("malinLME_NonDT_EEabove90.csv", header=T)
detrendedBelow10 = read.csv("malinLMEs_EEBelow10_1982_2014.csv", header=T)
detrendedAbove90 = read.csv("malinLMEs_EEOver90_1982_2014.csv", header=T)

## catch shares database with matches
#RamEDFmatching <- read.csv("GFR_CS_allmatches.csv")
#ramToEDF <- merge(RamEDFmatching, RamStocksData, by.x = "IdOrig", by.y = "assessid" )

##########################################################################################
### Prep data for matching ############################################################
##########################################################################################

## Recruitment is traditionally reported already lagged to the year of fertilization so
# Rt corresponds to SSBt 

## filter to assessments with Rt and SSBt 
RamStocksData = RamStocksData[!is.na(RamStocksData$R) & 
                                              !is.na(RamStocksData$SSB),]

#find how many years a stock has valid R and SSB estimates
nYearsWithRAndSSBPerStock = aggregate(year ~ stockid,
                                         data=RamStocksData,
                                         FUN=length)
names(nYearsWithRAndSSBPerStock)[2] = "nYears"

#find out how many stocks have at least 10 years with necessary data
nrow(nYearsWithRAndSSBPerStock[nYearsWithRAndSSBPerStock$nYears>=10,])
# 309 stocks 

#merge that info back in
RamStocksData = merge(RamStocksData, nYearsWithRAndSSBPerStock,
                             by=c("stockid"))

############################################################################################
# Merge Malin's Life History file with the RAM dataset (filtered by stocks with R and SSB) #
# & Matching RAM stocks to LMEs 
############################################################################################
MalinMatches <- merge(RamStocksData, malinLHfile, by = "stockid", all.x = TRUE, all.y = FALSE)

#check how many of these stocks that match btwn Malin's LH file and RAM have LMEs assigned # 
unmatchedStocks = MalinMatches[is.na(MalinMatches$lme_number),]
matchedStocks = MalinMatches[!is.na(MalinMatches$lme_number),]
length(unique(matchedStocks$stockid)) # 184 matched stocks (with LME)
length(unique(unmatchedStocks$stockid)) # 127 unmatched stocks (with LME)

## write a file with this list
write.csv(matchedStocks, "matchedMalinStocksWithRAndSSBJan3.csv", row.names=F)
write.csv(unmatchedStocks, "unmatchedMalinStocksWithRAndSSBJan3.csv", row.names=F)

###################################################################################################
### Use the file that Steve Miller and Laura Dee created in 2012 to match RAM stocks to #######
## LME to find more matches (called "RAMstocks2012_matchedtoLME.csv") #########################
#   called RamStocksByLME <- read.csv("RAMstocks2012_matchedtoLME.csv") #######################
  
#*****Laura Note to self **** To Do: either fix dataset to include more LME info in malinLHfile, or try to match against our 
    # LME->stock mapping from FD fisheries, and then fall back on manual matching. 
    # below proceeds as if matchedStocks is the final dataset to use.
    # #Note that some stocks in this file appear in multiple LMEs because the stock range
    # #spans more than one LME

RamStocksByLME$lmeId <- as.factor(RamStocksByLME$lmeId)
##** what lmeID codes are these? SAUP or FishBase? *** 

# Remove stocks that are listed in multiple LMEs
# First, count them.
nLMEsPerStock = aggregate(lmeId ~ STOCKID,
                          data=RamStocksByLME,
                          FUN=length)
names(nLMEsPerStock)[2] = "nLMEs"
  # alternate way of counting stocks:
       # nLMEsPerStock2 = data.frame(table(RamStocksByLME$STOCKID))
      # names(nLMEsPerStock2) = c("STOCKID", "nLMEs")

# Second, remove stocks with nLMEs > 1. i.e. keep only rows where the column with nLMEs==1
nLMEsPerStock = nLMEsPerStock[nLMEsPerStock$nLMEs==1,]

# Third, use the nLMEsPerStock DF to filter RamStocksByLME to only include stocks with a unique LME
RamStocksWithUniqueLME = merge(RamStocksByLME, nLMEsPerStock,
                               by="STOCKID",
                               all.x=F, all.y=F) #these are defaults, but just to be explicit: only keep rows in both datasets that match
nrow(RamStocksWithUniqueLME)

# Now, use this secondary LME <-> stock mapping to augment the matches we found using Malin's LME file
overallMatches = merge(MalinMatches, RamStocksWithUniqueLME[,c("STOCKID","lmeId")],
                       by.x=("stockid"),
                       by.y=("STOCKID"),
                       all.x=T,
                       all.y=F)
nrow(overallMatches) #should be 12779  
#LD to check: ****** It is 13077 now maybe due to the extra years by not needing lagged SSB data *****

#### Create LME master column: use Malin's if present. Otherwise, use LME info from other file. ####
overallMatches$lmeMatched = overallMatches$lme_number
overallMatches$lmeMatched[is.na(overallMatches$lmeMatched)] = overallMatches$lmeId[is.na(overallMatches$lmeMatched)]

# See what that bought us in terms of stocks that are matched to an LME
length(unique(overallMatches$stockid[!is.na(overallMatches$lme_number)]))  # 184
length(unique(overallMatches$stockid[!is.na(overallMatches$lmeId)]))       # 202
length(unique(overallMatches$stockid[!is.na(overallMatches$lmeMatched)]))  # 232

# how many stocks don't currently have LME matches but do have at least 10 years w/necessary R & SSB data?
length(unique(overallMatches$stockid[is.na(overallMatches$lmeMatched) & overallMatches$nYears>=10])) #78 out of 79 with no LME match

#write out results #writing out results for all stocks (even ones with <10years of R timeseries)
unmatchedOverallStocks = overallMatches[is.na(overallMatches$lmeMatched),]
matchedOverallStocks = overallMatches[!is.na(overallMatches$lmeMatched),]
write.csv(unmatchedOverallStocks, "unmatchedOverallStocksWithRAndSSBJan3.csv", row.names=F)
write.csv(matchedOverallStocks, "matchedOverallStocksWithRAndSSBJan3.csv", row.names=F)

#write out results after filtering to stocks with at least 10 years of data
write.csv(unmatchedOverallStocks[unmatchedOverallStocks$nYears>=10,], "unmatchedOverallStocksWith10yearsRAndSSBJan3.csv", row.names=F)
write.csv(matchedOverallStocks[matchedOverallStocks$nYears>=10,], "matchedOverallStocksWith10yearsRAndSSBJan3.csv", row.names=F)


###############################################################################################
##### Merge RAM data (stocks with LME matches) with binned temp data by LME & year  #############
##############################################################################################
modelData = merge(matchedOverallStocks, nonDetrendedBelow10,
                  by.x=c("lme_number", "year"),
                  by.y=c("lmeId", "year"))

modelData = merge(modelData, nonDetrendedAbove90,
                  by.x=c("lme_number", "year"),
                  by.y=c("lmeId", "year"))

modelData = merge(modelData, detrendedBelow10,
                  by.x=c("lme_number", "year"),
                  by.y=c("lmeId", "year"))

modelData = merge(modelData, detrendedAbove90,
                  by.x=c("lme_number", "year"),
                  by.y=c("lmeId", "year"))

write.csv(modelData, "matchedOverallStocksWith10yearsRAndSSBandTEMPJan3.csv")

########################################################################################
### RUN PRELIM MODELS ##################################################################
#######################################################################################
modelData$lme_number = as.factor(modelData$lme_number)
modelData$lmeMatched = as.factor(modelData$lmeMatched)

### *** NEED TO DETERMINE WHAT RECRUITMENT RESPONSE VARIABLES TO EXAMINE ****
modelData$outcome = modelData$R/modelData$SSB

# ### Look at subset of the data for  "California Current"
# CClmedata = modelData[modelData$lmeId==3,]
# CClmedata = modelData[modelData$lmeMatched==3,]


# mod1 = lm(outcome ~ stockid + lme_number:fracExtremeMonthsBelow10NonDetrended +
#             lme_number:fracExtremeMonthsAbove90NonDetrended + lme_number:fracExtremeMonthsBelow10 + 
#             stockid:fracExtremeMonthsAbove90,
#           data = modelData)
# 
# summary(mod1)


model1 = lm(R ~ stockid + SSB + stockid:fracExtremeMonthsBelow10NonDetrended +
              stockid:fracExtremeMonthsAbove90NonDetrended + stockid:fracExtremeMonthsBelow10 + 
              stockid:fracExtremeMonthsAbove90 + fecundity + eggdiamave + K + Troph + MaxLengthTL,
            data = modelData)


model2 = lm(R ~ lmeMatched + SSB + stockid:fracExtremeMonthsBelow10NonDetrended +
              stockid:fracExtremeMonthsAbove90NonDetrended + stockid:fracExtremeMonthsBelow10 + 
              stockid:fracExtremeMonthsAbove90 + fecundity + eggdiamave + K + Troph + MaxLengthTL,
            data = modelData)
summary(model2)

model3 = lm(R ~ stockid + SSB + stockid:fracExtremeMonthsBelow10NonDetrended +
              stockid:fracExtremeMonthsAbove90NonDetrended + stockid:fracExtremeMonthsBelow10 + 
              stockid:fracExtremeMonthsAbove90 ,
                            data = modelData)
summary(model3)


model3 = lm(R ~ lmeId + SSB + K + Troph + eggdiamave + MaxLengthTL + fecundity +
              stockid:fracExtremeMonthsAbove90 ,
            data = modelData)
summary(model3)

model3 = lm(R ~ lmeId + SSB + K + Troph + eggdiamave + MaxLengthTL + fecundity +
              stockid:fracExtremeMonthsAbove90NonDetrended,
            data = modelData)
summary(model3)


