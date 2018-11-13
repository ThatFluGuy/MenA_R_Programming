

#### Program information ######################################################
# Source file name: ModelInputUtilities.R                                   #
#_______________________________________ ______________________________________#
# Input datasets: specify folder containing downloads from                    #
#   https://montagu.vaccineimpact.org/
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
#       -vaccination scenarios
#_____________________________________________________________________________#
# Created  11/8 /18, by Chris Stewart chris.c.stewart@kp.org                  #                                   
#_____________________________________________________________________________#


#_____________________________________________________________________________#
# This function returns a dataset with a row for each year of simulation,     #
# with total pop, death rate, birth rate, and infant mortality rate           # 
#       ASSUMES FILENAMES CONTAIN: "cdr_both", "cbr_both" and "imr_both"      #
#_____________________________________________________________________________#

GetDemographicParameters<-function(directory, mycountry, start, end, fillThreshold=1) {
  library(dplyr)
  library(lubridate)
  setwd("\\\\HOME/stewcc1/MenAModel/R_programming")
  source("ModelInputUtilities.R")
  setwd(path)
  flist<-list.files(path)
  totpop<-flist[grepl("tot_pop_both",flist)==TRUE]
  #check length of totpop to validate filename
  dfpop<-read.csv(totpop[1])
  #check nrow(ctrypop to validate country)
  #validate country with total pop file;
  ctrypop<-dfpop[dfpop$country_code==mycountry, c("country_code", "country", "year","gender", "value")]
  if (nrow(ctrypop) > 0) {
    ctrypopfull<-checkVIMCdates(mydata=ctrypop, startyr=year(start), endyr=year(end), threshold=1)
    ctrypopfull%>%group_by(country)%>%summarize(min(year), max(year))
  } else {
    #message: country not found, exit function
  }
  
  cbr<-flist[grepl("cbr_both",flist)==TRUE]
  dfbirth<-read.csv(cbr[1])
  ctrybirth<-dfbirth[dfbirth$country_code==mycountry, c("country_code", "country", "year","gender", "value")]
  if (nrow(ctrybirth) > 0) {
    ctrybirthfull<-checkVIMCdates(mydata=ctrybirth, startyr=year(start), endyr=year(end), threshold=1)
    ctrybirthfull%>%group_by(country)%>%summarize(min(year), max(year))
  } else {
    #message: country not found, exit function
  }
  
  build1<-merge(x=ctrypopfull, y=ctrybirthfull, by=c("country_code", "country", "gender", "year"), all=TRUE)
  colnames(build1)[colnames(build1)=="value.x"] <- "totalpop"
  colnames(build1)[colnames(build1)=="value.y"] <- "birthrate"
  build1$births<-build1$totalpop*build1$birthrate
  
  #keep empty age_from and age_to from imr file to preserve format
  imr<-flist[grepl("imr_both",flist)==TRUE]
  dfim<-read.csv(imr[1])
  ctryimr<-dfim[dfim$country_code==mycountry, c("country_code", "country", "age_from", "age_to", "year","gender", "value")]
  if (nrow(ctryimr) > 0) {
    ctryimrfull<-checkVIMCdates(mydata=ctryimr, startyr=year(start), endyr=year(end), threshold=1)
    ctryimrfull%>%group_by(country)%>%summarize(min(year), max(year))
  } else {
    #message: country not found, exit function
  }
  
  build2<-merge(x=build1, y=ctryimrfull, by=c("country_code", "country", "gender", "year"), all=TRUE)
  cdr<-flist[grepl("cdr_both",flist)==TRUE]
  dfcdr<-read.csv(cdr[1])
  ctrycdr<-dfcdr[dfcdr$country_code==mycountry, c("country_code", "country", "year","gender", "value")]
  if (nrow(ctrycdr) > 0) {
  ctrycdrfull<-checkVIMCdates(mydata=ctrycdr, startyr=year(start), endyr=year(end), threshold=1)
  ctrycdrfull%>%group_by(country)%>%summarize(min(year), max(year))
  } else {
    #message: country not found, exit function
  }
  
  build3<-merge(x=build2, y=ctrycdrfull, by=c("country_code", "country", "gender", "year"), all=TRUE)
  colnames(build3)[colnames(build3)=="value.x"] <-"imr"
  colnames(build3)[colnames(build3)=="value.y"] <-"v"
  
  return(build3)
  
}

checkVIMCdates<-function(mydata, startyr, endyr, threshold=1) {
  #assume data has variables country and year
  #will fill in up to a threshold (default = 1 year) with values from nearest year
  library(dplyr)
  datesum<-mydata%>%dplyr::group_by(country)%>%dplyr::summarize(minyear=min(year), maxyear=max(year))
  if (datesum$minyear > startyr) {
    prediff<-datesum$minyear-startyr
    if (prediff <= threshold) {
      #fix by filling
      prefill<-mydata[mydata$year==datesum$minyear,]
      for (i in 1:prediff) {
        prefill$year<-datesum$minyear-i
        mydata<-rbind(mydata, prefill)
      }
    }
    # over threshold message
  }
  if (datesum$maxyear < endyr) {
    if (datesum$minyear-startyr <= threshold) {
      postdiff<-endyr-datesum$maxyear
      if (postdiff <= threshold) {
        postfill<-mydata[mydata$year==datesum$maxyear,]
        for (i in 1:postdiff) {
          postfill$year<-datesum$maxyear+i
          mydata<-rbind(mydata, postfill)
        }
      } 
    }
  }
  return(mydata)
}

#_____________________________________________________________________________#
# This function returns a vector with 7 values, one for each 5-year age       #
#  band up to 30, calculated from quinquennial file, for the year closest to  #
#  specified start of simulation.                                             #
#  ASSUMES FILENAME CONTAINS "qq_pop_both"                                    #
#_____________________________________________________________________________#
GetPopAgeDist<-function(mycountry, start, directory) {
  setwd(path)
  flist<-list.files(path)
  qqfile<-flist[grepl("qq_pop_both",flist)==TRUE]
  qqdf<-read.csv(qqfile)
  #record for every 5th year - want the closest to start or before start?
  mround <- function(x,base){ 
    base*round(x/base) 
  }
  popyr<-mround(year(start), 5)
  qqkeep<-qqdf[qqdf$country_code==mycountry & qqdf$year==popyr,]
  qqkeep$ageband<-ifelse(qqkeep$age_from < 30, paste0("Age_", qqkeep$age_from, "_",qqkeep$age_to), "Age_30")
  bands<-qqkeep%>%group_by(country_code, ageband)%>%summarize(tot=sum(value),minage=min(age_from)) 
  totpop<-qqkeep%>%group_by(country_code)%>%summarize(totpop=sum(value))
  numsall<-merge(x=bands, y=totpop, by="country_code")
  numsall$fraction<-numsall$tot/numsall$totpop
  agedist <- numsall[order(numsall$minage),]
  return(agedist$fraction)
}


#_____________________________________________________________________________#
# This function returns a dataset with a row for each year of simulation,     #
# with DosesCampaign, CoverRoutine, and AgeLimCampaign                        # 
#   NOTE: ASSUMES FILENAMES CONTAIN: "mena-routine" and "mena-routine"        #
#    NOTE: NOT CURRENTLY LIMITED TO YEARS OF SIM                              #
#_____________________________________________________________________________#
GetVaccScenario<-function(mycountry, scenario, directory) {
  setwd(directory)
  flist<-list.files(directory)
  if (scenario=="routine") {
    filename<-flist[grepl("mena-routine",flist)==TRUE]
  }
  if (scenario=="campaign") {
    filename<-flist[grepl("mena-campaign",flist)==TRUE]
  }
  #   if (scenario=="none") {
  #    filename<-flist[grepl("mena-no-vacc",flist)==TRUE] #empty file, data types do not match others
  #  }

  dfvacc<-read.csv(filename[1], stringsAsFactors = FALSE)
  ctryvacc<-dfvacc[dfvacc$country_code==mycountry, c("country_code", "country", "year","vaccine", "activity_type", "age_last", "target" , "coverage")]
  colnames(ctryvacc)[colnames(ctryvacc)=="coverage"] <-"CoverRoutine"
  colnames(ctryvacc)[colnames(ctryvacc)=="age_last"] <-"AgeLimCampaign"
  ##target has "<NA>" as character, hosing conversion
  ctryvacc$DosesCampaign<-as.numeric(ctryvacc$target)
  newdf<-subset(ctryvacc, select=-c(target))
  #better way to make dataset with same structure if we need it:
  #if (scenario=="none") {
   # newdf<-newdf[newf$country_code=="XXX",]
  #}
}