#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# File name: fxModelInputs.R                                                  #
# Contact: chris.c.stewart@kp.org, michael.l.jackson@kp.org                   #
# Version Date 08/20/21                                                       #
#_____________________________________________________________________________#
# This program contains a series of custom functions for importing data used  #
# by the NmA simulation programs. These are ultimately downloaded from the    #
# Vaccine Impact Modeling Consortium at: https://montagu.vaccineimpact.org/   #
#_____________________________________________________________________________#
# Functions in this program:                                                  #
# (1) GetMontaguDemogData - downloads data from Montagu via API.              #
# (2) GetDemographicParameters                                                #
# (3) checkVIMCdates                                                          #
# (4) GetPopAgeDist                                                           #
# (5) GetVaccScenario                                                         #
# (6) GetDiseaseStateDist                                                     #
# (7) GetWAIFWmatrix                                                          #
# (8) GetModelParams                                                          #
# (9) GetLifeEx                                                               #
# (10) GetModelPct                                                            #
# (11) combineOutputFiles                                                     #
# (12) InitializePopulation                                                   #
#_____________________________________________________________________________#
# Parameters:                                                                 #
# Start and end dates of simulation - to be specified in calling program      #
# 3-letter country code                                                       #       
#	directory containing downloaded file                                        #
#	threshold - number of years to fill in if data is missing, defaults to 1    #
#  (Note: this is because although total pop goes to 2100, cbr, cdr, imr may  #
#   end at 2099- I just copy the nearest year's value)                        #
#_____________________________________________________________________________#
# Purpose: return a vector or dataframe of parameters for:                    #
#       -distribution for initializing population                             #
#       -demographics: birth, death, infant mortality rate for simulation     #
#       -vaccination scenarios                                                #
#_____________________________________________________________________________#
# These functions were written by Chris Stewart, with some editing by Chloe   #
# Krakhauer and Mike Jackson                                                  #
#_____________________________________________________________________________#


### (1) Download data from Montagu API ########################################
# Not currently used.                                                         #
# PURPOSE: GET DATA FROM API at "montagu.vaccineimpact.org"                   #
# This function does not work from KPWA network but tested from my laptop     #
# 11/16/18, "downloads the 5 data files we need by reading into dataframe     # 
#    then writing csv to the directory specified.  Reconstructs near-real     #
# filenames(compared to manual download),will be recognized by downstream fxn #
# requires R package "montagu", install with:                                 #
# install.packages("drat") # if needed                                        #
# drat:::add("vimc")                                                          #
# install.packages("montagu")                                                 #
#_____________________________________________________________________________#

GetMontaguDemogData<-function(username=NULL, password=NULL, touchstone="201710gavi-5", 
                               destpath=NULL){

  svr <- montagu_server(name="production", hostname="montagu.vaccineimpact.org", 
                        username=username, password=password)
  montagu_server_global_default_set(svr)
  tchlist<-montagu_touchstones_list(svr)
  
  if (touchstone %in% as.vector(tchlist$id)) {
    dlist<-montagu_demographics_list(touchstone_id = touchstone)
    demogidlist<-as.vector(dlist$id)
    needed_data<-c("cbr", "cdr", "unwpp_imr", "qq_pop", "tot_pop")
    
    for (i in 1:length(needed_data)) { 
      if (needed_data[i] %in% demogidlist) {
        dat<-montagu::montagu_demographic_data(type=needed_data[i], touchstone_id=tch)
        datrow<-dlist[dlist$id==needed_data[i],]
        #filename like: 201710gavi-5_dds-201710_unwpp_imr_both.csv
        # [touchstone id + source + demog id + _both]
        #but it doesn't matter, my function lookst for demog id in name
        
        if (length(destpath)>0) {
          filename<-paste0(datadir, touchstone, datrow$source, datrow$id, "_both.csv")
          write.csv(dat, filename)
          print(filename)
        }
      }
    }
  } 
  else {print("Touchstone not found.  Touchstone list:")
    print(tchlist)
  }
}

### (2) GetDemographicParameters ##############################################
# This function returns a dataset with a row for each year of simulation,     #
# with total pop, death rate, birth rate, and infant mortality rate           # 
# ASSUMES FILENAMES CONTAIN: "tot_pop_both", "cdr_both", "cbr_both","imr_both"#
# Chloe 7/12/19: now assume death rates and numbers of deaths from files      #
# called "p_dying_both" and "births" respectively.                            #
# Mike 12/10/19: Re-write so that death probabilities are estimated from      #
# annual age-specific UN interpolated population sizes, rather than probs     #
# from Montagu. For the 2019 model year, keep the current approach of using   #
# average death probs in 5-year age buckets, to keep things consistent with   #
# the MenA_OneSim.R program. Future work may move this to single-year buckets.#
#_____________________________________________________________________________#

GetDemographicParameters <- function(path, mycountry, start, end, fillThreshold=1) {
  
  setwd(path)
  
  #### (1) Get total population size
  totpop <- GetFilename(path, "tot_pop_both")
  if (is.character(totpop)==FALSE) { stop(mymsg) }
  dfpop <- read.csv(totpop)
  
  if (CheckDemogFileStructure(mycountry=mycountry, mydf=dfpop, dfdesc="tot_pop_both")==FALSE) { stop (filemsg)}
  ctrypop <- dfpop[dfpop$country_code==mycountry, c("country_code", "year", "value")]
  ctrypopfull <- checkVIMCdates(mydata=ctrypop, startyr=year(start), endyr=year(end), threshold=fillThreshold)
  
  if (is.data.frame(ctrypopfull)==FALSE) { stop(paste(datemsg, " tot_pop_both")) }
  ctrypopfull %>% group_by(country_code) %>% summarize(min(year), max(year))
  
  #### (2) Get birthrate and number of births each year
  cbr <- GetFilename(path, "cbr_both")
  births <- GetFilename(path, "births")
  if (is.character(cbr)==FALSE) { stop(mymsg) }
  if (is.character(births)==FALSE) { stop(mymsg) }
  
  dfbirth <- read.csv(cbr)
  numbirth <- read.csv(births)
  
  if (CheckDemogFileStructure(mycountry=mycountry, mydf=dfbirth, dfdesc="cbr_both")==FALSE) { stop (filemsg)}
  if (CheckDemogFileStructure(mycountry=mycountry, mydf=numbirth, dfdesc="births")==FALSE) { stop (filemsg)}
  ctrybirth <- dfbirth[dfbirth$country_code==mycountry, c("country_code", "year", "value")]
  numbirth_ctry <- numbirth[numbirth$country_code==mycountry, c("country_code", "year", "value")]
  ctrybirthfull <- checkVIMCdates(mydata=ctrybirth, startyr=year(start), endyr=year(end), threshold=fillThreshold)
  numbirthfull <- checkVIMCdates(mydata=numbirth_ctry, startyr=year(start), endyr=year(end), threshold=fillThreshold)
  
  if (is.data.frame(ctrybirthfull)==FALSE) { stop(paste(datemsg, " cbr_both")) }
  if (is.data.frame(numbirthfull)==FALSE) { stop(paste(datemsg, " births")) }
  ctrybirthfull %>% group_by(country_code) %>% summarize(min(year), max(year))
  numbirthfull %>% group_by(country_code) %>% summarize(min(year), max(year))
 
  build0 <- merge(x=ctrypopfull, y=ctrybirthfull, by=c("country_code", "year"), all=TRUE)
  colnames(build0)[colnames(build0)=="value.x"] <- "totalpop"
  colnames(build0)[colnames(build0)=="value.y"] <- "birthrate"
  build1 <- merge(x=build0, y=numbirthfull, by=c("country_code", "year"), all=TRUE)
  colnames(build1)[colnames(build1)=="value"] <- "births"
  
  #### (3) Calculate yearly probability of death by age group
  dr <- GetFilename(path, "int_pop_both")
  if (is.character(dr)==FALSE){ stop(mymsg) }
  int.pop.all <- read.csv(dr)
  if (CheckDemogFileStructure(mycountry=mycountry, mydf=int.pop.all, dfdesc="int_pop_both")==FALSE) { 
    stop (filemsg)
  }
  int.pop <- int.pop.all[int.pop.all$country_code==mycountry & int.pop.all$year>=1950, 
                         c("country_code", "year", "age_from", "age_to", "value")]
  
  # Error checking: should have 15251 rows (101 age groups x 151 years)
  if(dim(int.pop)[1] != 15251){
    stop("Incorrect number of entries in int_pop_both")
  }
  int.pop <- int.pop[order(int.pop$year, int.pop$age_from),]

  # For each age group, calculate deaths each year
  # First, ignore the age 100-120 (will set death rate to 1.0) 
  # and ignore year 2100 (since we don't have 2101 data for comparison)
  int.pop$next_value <- ifelse(int.pop$year==2100, NA,
                               ifelse(int.pop$age_to==120, NA, int.pop$value[103:15251]))
  int.pop$death.prob <- (int.pop$value - int.pop$next_value)/int.pop$value
  
  # Set death.prob to 1.0 for age 100-120
  int.pop$death.prob <- ifelse(int.pop$age_to==120, 1, int.pop$death.prob)
  
  # For 2100, use death probabilities from 2099
  int.pop$death.prob[int.pop$year==2100] <- int.pop$death.prob[int.pop$year==2099]
  
  # Fix death.prob where value is zero (set death prob to zero)
  int.pop$death.prob[int.pop$value==0] <- 0
  
  # In a handful of cases (at the upper age ranges) population size increases from
  # year to year. Rather than assume a negative death prob (which might break the
  # simulation in an unpredictable way) set death prob to zero for these.
  int.pop$death.prob[int.pop$death.prob < 0] <- 0
  
  # For every 5-year age bucket (except infant mortality) assume midpoint is death prob
  int.pop.small <- int.pop[int.pop$age_from %in% c(0, 2, seq(7, 82, by=5)),]
  int.pop.small$name <- ifelse(int.pop.small$age_from==0, "imr",
                               paste("dr", int.pop.small$age_from-2, int.pop.small$age_from+2, sep=""))
  int.pop.small$name <- ifelse(int.pop.small$name=="dr04", "dr14", int.pop.small$name)

  build2a <- spread(int.pop.small[,c("year", "death.prob", "name")], name, death.prob)
  build2a <- build2a[, c(1, 19, 3, 13, 2, 4:12, 14:18)]
  
  build2 <- merge(x=build1, y=build2a, by="year", all=TRUE)
  
  return(build2)
}

### (3) checkVIMCdates ########################################################
# Function to check dates.                                                    #

checkVIMCdates <- function(mydata, startyr, endyr, threshold=1) {
  # Assume data has variables country and year
  # Will fill in up to a threshold (default = 1 year) with values from nearest year
  datemsg<<-""
  datesum <- mydata %>% dplyr::group_by(country_code) %>% 
    dplyr::summarize(minyear=min(year), maxyear=max(year))
  if (datesum$minyear > startyr) {
    prediff <- datesum$minyear-startyr
    if (prediff <= threshold) {
      #fix by filling
      prefill <- mydata[mydata$year==datesum$minyear,]
      for (i in 1:prediff) {
        prefill$year <- datesum$minyear-i
        mydata <- rbind(mydata, prefill)
      }
    } else {
    # over threshold message
      datemsg<<- paste0("There is a ", prediff, "-year gap in the demographic data compared to your simulation begin date.  You can increase the threshold parameter to fill in data from the closest year.")
      return(FALSE)
    }
  }
  if (datesum$maxyear < endyr) {
    if (datesum$minyear-startyr <= threshold) {
      postdiff <- endyr-datesum$maxyear
      if (postdiff <= threshold) {
        postfill <- mydata[mydata$year==datesum$maxyear,]
        for (i in 1:postdiff) {
          postfill$year <- datesum$maxyear+i
          mydata <- rbind(mydata, postfill)
        }
      } else {
        # over threshold message
        datemsg<<-paste0("There is a ", postdiff, "-year gap in the demographic data compared to your simulation end date.  You can increase the threshold parameter to fill in data from the closest year.")
        return(FALSE)
      } 
    }
  }
  return(mydata)
}

### (4) GetPopAgeDist #########################################################
# This function returns a vector with 7 values, one for each 5-year age       #
#  band up to 30, calculated from quinquennial file, for the year closest to  #
#  specified start of simulation. Added names=age band 11/15/18               #
#  ASSUMES FILENAME CONTAINS "qq_pop_both"                                    #
#  Called by InitializePopulation.R                                           #
#_____________________________________________________________________________#

GetPopAgeDist<-function(path, mycountry, start) {
  disterr<<-""
  setwd(path)
  qqfile <- GetFilename(path, "qq_pop_both")
  if (is.character(qqfile)==FALSE) { stop(mymsg) }
  qqdf <- read.csv(qqfile)
  if (CheckDemogFileStructure(mycountry=mycountry, mydf=qqdf, dfdesc="qq_pop_both")==FALSE) { stop (filemsg)}
  #record for every 5th year - want the closest to start or before start?
  mround <- function(x,base){ 
    base*round(x/base) 
  }
  popyr <- mround(year(start), 5)
  qqkeep <- qqdf[qqdf$country_code==mycountry & qqdf$year==popyr,]
 
  if (nrow(qqkeep) > 0) {
    if (!(DemogNumVarExists("age_from", qqkeep) & DemogNumVarExists("age_to", qqkeep))) { 
      disterr<<- "Non-numeric data in age band variables"
      return(FALSE) 
    }
    
    if (!(DemogNumVarExists("age_to", qqkeep))) { return(FALSE) }
    #check min and max age band - now they all go up to 120 but lets be generous and say 90
    if (min(qqkeep$age_from > 0 || max(qqkeep$age_to) < 90) == FALSE) {
      qqkeep$ageband <- paste0("Age_", qqkeep$age_from, "_",qqkeep$age_to)
      bands <- qqkeep %>% group_by(country_code, ageband) %>% 
        summarize(tot=sum(value),minage=min(age_from)) 
      
      totpop <- qqkeep %>% group_by(country_code) %>% summarize(totpop=sum(value))
      if (totpop[,2] > 0) {
        numsall <- merge(x=bands, y=totpop, by="country_code")
        numsall$fraction <- numsall$tot/numsall$totpop
        agedist <- numsall[order(numsall$minage),]
        dist <- agedist$fraction
        names(dist) <- agedist$ageband
        return(dist)
      } else {
        disterr<<-"Population value is zero.  Please check qqpop file."
        return(FALSE)
      }
    } else { 
      disterr<<-"Incomplete age bands for this country and year" 
      return(FALSE)
    }
  } else {
    disterr<<-"No age distribution input found for this country and year"
    return(FALSE)
    }
}


### (5) GetVaccScenario #######################################################
# This function returns a dataset with a row for each year of simulation,     #
# with DosesCampaign, CoverRoutine, and AgeLimCampaign                        # 
#   NOTE: ASSUMES FILENAMES CONTAIN: "mena-routine" and "mena-campaign"       #
#    NOTE: NOT CURRENTLY LIMITED TO YEARS OF SIM, could use GetVIMCdates?     #
#     12/6/18 copying destring from taRifx to deal with "<NA>" in vacc files  #
#     11Dec2019: Added sub-scenario option to select between more options     #
#_____________________________________________________________________________#
destring <- function(x,keep="0-9.-") {
  return( as.numeric(gsub(paste("[^",keep,"]+",sep=""),"",x)) )
}

GetVaccScenario <- function(mycountry, scenario, sub.scenario, directory) { #sub.scenario allows for selection between bestcase and default vaccination scenario files
  vaccmsg <<- ""
  setwd(directory)
  if (scenario=="none") sub.scenario <- "NA"  #can't have a sub-scenario when there's no vaccinations
  if (scenario=="routine" | scenario=="both") {
    filename <- GetFilename(directory, "mena-routine", sub.scenario)
  }
  if (scenario=="campaign") {
    filename <- GetFilename(directory, "mena-campaign", sub.scenario)
  }
  if (scenario=="booster") {
    filename <- GetFilename(directory, "mena-booster", sub.scenario)
  }
  if (is.character(filename)==FALSE) { stop(mymsg) }
  
  dfvacc <- read.csv(filename, stringsAsFactors = FALSE)
  if (IsCountryAndColAvailable(country_code=mycountry,mydf=dfvacc, forVacc=1)==FALSE) { stop(countrymsg) }
  
  # Target and year validated above
  if (scenario %in% c("routine", "both", "booster")){
      if (!(DemogNumVarExists("coverage", dfvacc))) { 
        vaccmsg<<-"coverage variable missing from vaccination file"
        return(FALSE) 
        }
  }
  
  # Drop the records for ACWXY vaccination (these are handled separately)
  dfvacc <- dfvacc[dfvacc$vaccine=="MenA",]
  
  # Select needed variables
  ctryvacc <- dfvacc[dfvacc$country_code==mycountry, 
                     c("country_code", "year", "activity_type", "target" , "coverage",
                       "age_first", "age_last")]
  colnames(ctryvacc)[colnames(ctryvacc)=="coverage"] <-"CoverRoutine"
  
  # Convert ages to months (only needed for campaigns)
  ctryvacc$age_first <- ctryvacc$age_first*12
  ctryvacc$age_last <- ctryvacc$age_last*12 + 11
  
  ##target has "<NA>" where activity type = routine, hosing conversion
  #still getting coercion warning
  #getting this even though not strictly required by routine option
  ctryvacc$DosesCampaign <- destring(ctryvacc$target)
  newdf <- subset(ctryvacc, select=-c(target))
  
  return(newdf)
}


### (6) GetDiseaseStateDist ###################################################
# Function GetDiseaseStateDist, called by InitializePopulation.R              #
# Reads dist_both.csv, which is supplied with scripts; format should not vary #
#_____________________________________________________________________________#
GetDiseaseStateDist <- function(directory, region) {
  setwd(directory)
  dxfile <- GetFilename(directory, "dist_both.csv")
  if (is.character(dxfile)==FALSE) { 
    stop(mymsg) 
    print("File [dist_both.csv] is packaged with the R scripts and should be in the same directory.")
  }
  dist <- read.csv(dxfile, stringsAsFactors = TRUE)
  distcol <- ifelse(region=='hyper', 4, 3)
  statefract <- as.vector(dist[,distcol]) # fraction of each disease state in each of 7 population groups
  return(statefract)
}

### (7) GetWAIFWmatrix #################################################################
# Get WAIFWmatrix: constructs expanded WAIFW matrix from supplied parameter data frame;# 
#______________________________________________________________________________________#

GetWAIFWmatrix <- function(params) {
  
  Dwaifw <- matrix(0, nrow=5, ncol=5)
  Dwaifw[1,] <- rep(params$bd1, 5)
  Dwaifw[2,] <- rep(params$bd2, 5)
  Dwaifw[3,] <- c(rep(params$bd3, 2), params$bd6, rep(params$bd3, 2))
  Dwaifw[4,] <- c(rep(params$bd4, 3), params$bd7, rep(params$bd4, 1))
  Dwaifw[5,] <- rep(params$bd5, 5)
  Rwaifw <- Dwaifw * params$br
  Rwaifw.expanded <- Rwaifw[rep(seq_len(nrow(Rwaifw)), times=c(60, 60, 60, 60, 1201)),] #Replicates the first 4 rows 60 times, the last row 1201 times.  One row per month, so ages 0-4, 5-9, 10-14, 15-19, 20-120
  Dwaifw.expanded <- Dwaifw[rep(seq_len(nrow(Dwaifw)), times=c(60, 60, 60, 60, 1201)),]
  wboth <- array(data=NA, dim=c(1441, 5, 2))
  wboth[,,1] <- Rwaifw.expanded
  wboth[,,2] <- Dwaifw.expanded
  dimnames(wboth)[[3]]<-c("rainy", "dry")
  return(wboth)
}



### (8) GetModelParams ########################################################
# GetModelParams: Reads posterior_parameters.csv, which is supplied with      #
# scripts. Subset to hyper or non-hyper region.                               #
# Goal is to prevent hard-coding of parameters in model                       #
#_____________________________________________________________________________#

GetModelParams<-function(path=scripts.dir, region.val, 
                         scaling=c(1, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75)) {
  
  # Inputs: 
  # path        - character scaler with directory where parameter file exists
  # region.val  - character scaler indicating region type (hyper vs not)
  # scaling     - numeric vector to scale hyper bd1-bd7 down to not-hyper

  if (!(region.val %in% c("hyper", "not_hyper"))){
    mymsg <- "region must be hyper or not_hyper"
    stop(mymsg)
  }
  
  setwd(path)
  param.file <- GetFilename(path, "posterior_parameters.csv")
  if (is.character(param.file)==FALSE) { 
    stop(mymsg) 
    print("File [posterior_parameters.csv] is packaged with the R scripts and should be in the same directory.")
  }
  
  params <- read.csv(param.file, stringsAsFactors = FALSE)  #data frame
  if (region.val=="not_hyper"){
    params$bd1 <- params$bd1 * scaling[1]
    params$bd2 <- params$bd2 * scaling[2]
    params$bd3 <- params$bd3 * scaling[3]
    params$bd4 <- params$bd4 * scaling[4]
    params$bd5 <- params$bd5 * scaling[5]
    params$bd6 <- params$bd6 * scaling[6]
    params$bd7 <- params$bd7 * scaling[7]
  }
  return(params)
}


### (9) GetLifeExp ############################################################
# GetLifeExp: reads the life_ex_both file and creates a data.frame of life    #
# expectancy by age and year.                                                 #
#_____________________________________________________________________________#

GetLifeExp<-function(path, mycountry.s=mycountry) {
  
  # Inputs:
  # path        - character scalar with directory where life expectancy data exist
  # mycountry.s - character scalar with code for current country
  
  # Output:
  # data.frame with life expectancy by age and year
  
  setwd(path)
  
  lifex <- GetFilename(path, "life_ex_both")
  if (is.character(lifex)==FALSE) {stop(mymsg)}
  
  lifex.df <- read.csv(lifex)
  if (CheckDemogFileStructure(mycountry=mycountry.s, mydf=lifex.df, dfdesc="life_ex_both")==FALSE) { stop (filemsg)}

  lifex.df <- lifex.df[lifex.df$country_code==mycountry.s, 
                       c("age_from", "age_to", "year", "value")]
  
  # Expand to individual ages instead of age groups
  # First expand for the "missing" ages (e.g. 2-4 from the 1-4 group)
  # Then add the "present" ages (values of age_from)
  lifex.df$count <- lifex.df$age_to - lifex.df$age_from
  
  le.df <- data.frame(year=rep(lifex.df$year, times=lifex.df$count),
                      value=rep(lifex.df$value, times=lifex.df$count),
                      age_strt=rep(lifex.df$age_from, times=lifex.df$count),
                      counter=sequence(lifex.df$count))
  le.df$age_from <- le.df$age_strt + le.df$counter
  le.df <- rbind(lifex.df[, c("year", "value", "age_from")],
                 le.df[, c("year", "value", "age_from")])
  le.df <- le.df[order(le.df$year, le.df$age_from),]
  
  # Expand to every year instead of every 5 years
  results <- data.frame(year=rep(le.df$year, times=5),
                        value=rep(le.df$value, times=5),
                        age_from=rep(le.df$age_from, times=5),
                        counter=c(rep(0, times=length(le.df$year)),
                                  rep(1, times=length(le.df$year)),
                                  rep(2, times=length(le.df$year)),
                                  rep(3, times=length(le.df$year)),
                                  rep(4, times=length(le.df$year))))
  results$year <- results$year + results$counter
  
  results <- results[order(results$year, results$age_from), 
                     c("year", "age_from", "value")]
  names(results) <- c("year", "AgeInYears", "Life.Ex")
  
  # Duplicate 2099 to create 2100
  res2 <- results[results$year==2099,]
  res2$year <- 2100
  results <- rbind(results, res2)
  
  return(results)
    
}

### (10) Get proportion to model ##############################################
# For many countries, the target population for vaccination is only a portion #
# of the full country. This function imports data on the percent of the pop   #
# to use in the model. 

GetModelPct <- function(path=input.dir, mycountry.s=mycountry){
  # Inputs
  # path        - character scalar with directory where modeled percent data exist
  # mycountry.s - character scalar with code for current country

  # Outputs: scaler containing proportion of population to model  
  
  popmod.df <- read.csv(paste(path, "percent_pop_modeled.csv", sep="/"),
                        stringsAsFactors = FALSE)

  return(popmod.df$pct_pop_modeled[popmod.df$country_code==mycountry.s])  
}

### (11) Combine central estimate files #######################################
# For a defined vaccination program and sub-program, read in output files     #
# from all countries, format per VIMC specifications, and write the combined  #
# file to the output directory.                                               #       

combineOutputFiles <- function(path=output.dir, vacc_program="none",
                               vacc_subprogram="default", deliv.path,
                               touchstone){

  # Inputs
  # path            - character scalar indicating location of output files
  # vacc_program    - character scalar for program, as "none", "campaign", "both"
  # vacc_subprogram - character scalar for sub-program, as "default" or "bestcase
  # deliv.path      - character scalar for deliverables folder for compiled results
  
  # Output: writes a .csv file
  
  
  # (A) Check inputs
  if (!dir.exists(path)) {
    stop(paste(path, "is not a valid directory", sep=" "))
  }
  if (!(vacc_program %in% c("none", "campaign", "routine", "both", "booster"))){
    stop("vacc_progam must be one of 'none', 'campaign', 'routine', 'booster' or 'both'")
  }
  if (!(vacc_subprogram %in% c("default", "bestcase"))){
    stop("vacc_subprogram must be one of 'default', 'bestcase'")
  }
  if (!dir.exists(deliv.path)){
    stop(paste(path, "is not a valid directory", sep=" "))
  }
  
  # (B) Set up some needed data.frames and vectors
  names.df <- data.frame(
    country_code=c("BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "ERI", "ETH", "GHA", "GIN", "GMB",
                     "GNB", "KEN", "MLI", "MRT", "NER", "NGA", "RWA", "SDN", "SEN", "SSD", "TCD", "TGO",
                     "TZA", "UGA"),
    country=c("Burundi", "Benin", "Burkina Faso", "Central African Republic",
              "Cote d'Ivoire", "Cameroon", "Congo, the Democratic Republic of the",
              "Eritrea", "Ethiopia", "Ghana", "Guinea", "Gambia", "Guinea-Bissau", "Kenya",                                
              "Mali", "Mauritania", "Niger", "Nigeria", "Rwanda", "Sudan", "Senegal",                              
              "South Sudan", "Chad", "Togo", "Tanzania, United Republic of", "Uganda" ))
  
  oldnames=c("AgeInYears", "Cases", "Deaths", "DALYs", "cohortsize")
  newnames=c("age", "cases", "deaths", "dalys", "cohort_size")
  keepnames=c("disease",	"year", "age", "country", "country_name", "cohort_size",
              "cases", "dalys", "deaths", "cases_cwyx", "deaths_cwyx", "dalys_cwyx")


  output.df <- data.frame(disease=character(0), year=numeric(0), age=numeric(0),
                          country=character(0), country_name=character(0),
                          cohort_size=numeric(0), cases=numeric(0),
                          dalys=numeric(0), deaths=numeric(0), cases_cwyx=numeric(0),
                          deaths_cwyx=numeric(0), dalys_cwyx=numeric(0), 
                          stringsAsFactors = FALSE)

  # (C) Get list of all relevant output files and cycle through each
  # Read the file in, verify size, rename variables, and bind rows to output
  files.v <- list.files(path, pattern=paste(vacc_program, vacc_subprogram, sep="_"))
  files.v <- files.v[grep("^(?!PSA)", files.v, perl=TRUE)]
  
  for (f in 1:length(files.v)){
    country_code <- substr(files.v[f], 1, 3)
    result.df <- read.csv(paste(output.dir, files.v[f], sep="/"))
    
    size <- length(result.df$year)
    if (size != 7171){
      stop(paste("Input file is", size, "records and should be 7171", sep=" "))
    }
    
    result.df$disease <- rep("MenA", times=size)
    result.df$country <- rep(country_code, times=size)
    result.df$country_name <- as.character(rep(names.df$country[names.df$country_code==country_code], times=size))
    
    result.df <- result.df %>%
      rename_at(vars(oldnames), ~newnames)
    
    output.df <- bind_rows(output.df, result.df[, keepnames])
  }

  if(vacc_program=="none"){
    filename <- paste(deliv.path, "/mena-no-vaccination-", touchstone, ".MenA_KPW-Jackson.csv", sep="")
  } else {
   filename <- paste(deliv.path, "/mena-", vacc_program, "-", vacc_subprogram, "-", 
                    touchstone, ".MenA_KPW-Jackson.csv", sep="")
  }
  
  # (D) Output to the appropriate directory, with error-handling
  outval <- tryCatch({
      write.csv(x=output.df, file=filename, row.names=FALSE)
    }, warning=function(cond){
      message("Trying to output gave a warning:")
      message(cond)
      return(NA)
    }, error=function(cond){
      message("Trying to output gave an error:")
      message(cond)
      return(NA)
    }
  )
  
  if (is.null(outval)){print(paste("Output written to:", filename, sep=" "))}
  
}

#_____________________________________________________________________________#
# INITIALIZE POPULATION Input datasets: PopInputAgeDist.csv, dist_both.csv    #
#_____________________________________________________________________________#
# Parameters:                                                                 #
# Start and end dates of simulation - to be specified in calling program      #
# Startpop from country params                                                #
#	country name to get age distributions from popinputagedist.csv              #
#	region = hyper/not hyper to get age-specific disease state distribution.    #
#_____________________________________________________________________________#
# Purpose: fill first "slice" of the simulation matrix with the population at #   
# the starting point, divided by age and disease state. Matrix is 360 monthly #
# ages by 8 disease states.                                                   #
#_____________________________________________________________________________#
# Created as function 3/5/18, by Chris Stewart chris.c.stewart@kp.org         #
# Changes: remove loop, use vector multiplication instead,                    #
# 30 age group only needs one value in input vectors                          #    
# Change 10/23 EJ: removed the popsize argument, as it isn't used?            #
#_____________________________________________________________________________#
InitializePopulation<-function(scriptdir, inputdir, start, end, country="ETH", region="not_hyper", startSize) {
  initmsg<<-""
  #create the matrix
  dur<- ceiling(difftime(end, start, units = "weeks"))
  # poparray <- array(data=0, dim=c(361, 9, dur+1)) # Dimensions are [age groups, states, time]
  # Chloe edit 3/22/19: allow for monthly ages through age 120, including age=0 months.
  # Later on, for sake of comparison with Montagu estimates, will only look through age 70.
  # Note this could make code extremely (and unecessarily) slow...maybe just through age 70 if that's the case?
  # EJ edit: split Va into Vs and Vc
  poparray <- array(data=0, dim=c(1441, 10, dur+1)) # Dimensions are [age groups, states, time]
  dimnames(poparray)[[2]] <- c("Ns","Nc","Ls","Lc","Hs","Hc","Di","Vs", "Vc", "Inc")
  #parameters for initializing
  mypop<-GetPopAgeDist(path=inputdir, mycountry=country, start=start) 
  # Chloe edit 3/22/19: Now, this should start with pop size of age bands between 0-4 and 100-120.
  # Chloe edit 3/29/10: The error message below was originally intended for 7 age bands; 
  # I'm leaving this unedited for now, since I'm not sure what final groups of age bands we will use. 
  if (length(mypop) < 7) {
    initmsg<<-"Incomplete age distribution information.  Please check the downloaded file [touchstone]_qq_pop_both.csv."
    return(FALSE)
  }
  ##get age-specific proportions of each disease state into vectors: 7 ages X 7 disease states
  ## Chloe 3/22/19: now, we have 21 ages and 7 disease states.
  ## The line below gets fraction of each disease state in each of 7 population groups (5-year age bands up to age 30).
  ## To avoid having too many edits in too many functions, I decided to leave the "GetDiseaseStateDist" function alone, and as before,
  # assume fractions in each disease state in the 30-year age band match those in the older age bands, i.e. can just repeat the
  # final 7 values as needed here to provide proportion estimates up to age 120...so, if limiting to age 70 above, will need to edit this.
  # EJ 10/23: Change to 9 disease states, allowing for initial vaccination status.  dist_both.csv has those values in it.
  statefract<-GetDiseaseStateDist(directory=scriptdir, region=myregion)
  if (length(statefract) < 63) {  
    initmsg<<-"Incomplete disease state distribution information.  Please check the file dist_both.csv provided with the scripts."
    return(FALSE)
  }
  #expand age group fraction as vector to match pop matrix dimension 1
  # Chloe edit 3/22/19: changing this to match new dims; see note above.
  # agefract <- c(rep(as.numeric(mypop[1:6]), each=60), as.numeric(mypop[7]))
  # Chloe 3/29: In new age bands, all are 5-year except for last, which is 120 years.
  agefract <- c(rep(as.numeric(mypop[1:20]), each=60), rep(as.numeric(mypop[21]), each=(12*20)+1))
  # chunks <- c(rep(60, each=360), 1)
  # Chloe: Looks as though the starting assumption is that population age bands are evenly divided 
  # among each included age.
  chunks <- c(rep(60, each=1200), rep((12*20)+1,each=(12*20)+1))
  
  # Exanding this to more age brackets.
  # statemx<-rbind(
  #  matrix(rep(statefract[1:7], each=60), nrow=60),
  #  matrix(rep(statefract[8:14], each=60), nrow=60),
  #  matrix(rep(statefract[15:21], each=60), nrow=60),
  #  matrix(rep(statefract[22:28], each=60), nrow=60),
  #  matrix(rep(statefract[29:35], each=60), nrow=60),
  #  matrix(rep(statefract[36:42], each=60), nrow=60),
  #  matrix(rep(statefract[43:49], each=1), nrow=1)
  # )
  
  # Assuming last 7 values of statefract (originally just for 30-34 age band) apply to all older ages (which was effectively 
  # implied before).
  # Separate matrix for each age chunk,
  # separate column for disease state, separate row for each month of age.
  # EJ 10/23: shifting from groups of 7 disease states to 9 disease states
  statemx<-rbind(
    matrix(rep(statefract[1:9], each=60), nrow=60),
    matrix(rep(statefract[10:18], each=60), nrow=60),
    matrix(rep(statefract[19:27], each=60), nrow=60),
    matrix(rep(statefract[28:36], each=60), nrow=60),
    matrix(rep(statefract[37:45], each=60), nrow=60),
    matrix(rep(statefract[46:54], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=60), nrow=60),
    matrix(rep(statefract[55:63], each=(12*20)+1), nrow=(12*20)+1)
  )
  #initialize - 3rd dimension stays at 1
  poparray[, "Ns", 1]<-(startSize*agefract/chunks)*statemx[,1]
  poparray[, "Nc", 1]<-(startSize*agefract/chunks)*statemx[,2]
  poparray[, "Ls", 1]<-(startSize*agefract/chunks)*statemx[,3]
  poparray[, "Lc", 1]<-(startSize*agefract/chunks)*statemx[,4]
  poparray[, "Hs", 1]<-(startSize*agefract/chunks)*statemx[,5]
  poparray[, "Hc", 1]<-(startSize*agefract/chunks)*statemx[,6]
  poparray[, "Di", 1]<-(startSize*agefract/chunks)*statemx[,7]
  poparray[, "Vs", 1]<-(startSize*agefract/chunks)*statemx[,8]
  poparray[, "Vc", 1]<-(startSize*agefract/chunks)*statemx[,9]
  #leave Inc 0
  return(poparray)
}
