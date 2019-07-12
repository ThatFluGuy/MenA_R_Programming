
#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: ModelInputUtilities.R                                     #
# Contact: chris.c.stewart@kp.org, michael.l.jackson@kp.org                   #
# Version Date 12/13/18                                                       #
#_______________________________________ _____________________________________#
# Input datasets: specify folder containing downloads from                    #
#   https://montagu.vaccineimpact.org/                                        #
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


#_____________________________________________________________________________#
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

# Chloe 7/12: note that death rates by age (TDB) and numbers of births will need to be added to this function (or something like it) to be used
# in the future; existence of that data in input folder currently assumed in functions below.

GetMontaguDemogData<-function( username=NULL, password=NULL, touchstone="201710gavi-5", destpath=NULL) {
  #GET DATA FROM APIdoes not work from KPWA network but tested from my laptop
  #destpath is where you want to put the data - should be same as path in GetDemographicParameters()
  #username and password for montagu site
  svr<-montagu_server(name="production", hostname="montagu.vaccineimpact.org", username=username, password=password)
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

#_____________________________________________________________________________#
# This function returns a dataset with a row for each year of simulation,     #
# with total pop, death rate, birth rate, and infant mortality rate           # 
# ASSUMES FILENAMES CONTAIN: "tot_pop_both", "cdr_both", "cbr_both","imr_both"#
# Chloe 7/12/19: now assume death rates and numbers of deaths from files
# called "p_dying_both" and "births" respectively.
#_____________________________________________________________________________#

# Chloe 7/12/19: kept crude birth rate data in here, but did use another file with number of births calculated.

GetDemographicParameters<-function(path, mycountry, start, end, fillThreshold=1) {
  setwd(path)
  totpop<-GetFilename(path, "tot_pop_both")
  if (is.character(totpop)==FALSE) { stop(mymsg) }
  dfpop<-read.csv(totpop)
  if (CheckDemogFileStructure(mycountry=mycountry, mydf=dfpop, dfdesc="tot_pop_both")==FALSE) { stop (filemsg)}
  ctrypop<-dfpop[dfpop$country_code==mycountry, c("country_code", "year", "value")]
  ctrypopfull<-checkVIMCdates(mydata=ctrypop, startyr=year(start), endyr=year(end), threshold=fillThreshold)
  if (is.data.frame(ctrypopfull)==FALSE) { stop(paste(datemsg, " tot_pop_both")) }
  ctrypopfull%>%group_by(country_code)%>%summarize(min(year), max(year))
  
  cbr<-GetFilename(path, "cbr_both")
  births <- GetFilename(path, "births")
  if (is.character(cbr)==FALSE) { stop(mymsg) }
  if (is.character(births)==FALSE) { stop(mymsg) }
  dfbirth<-read.csv(cbr)
  numbirth<-read.csv(births)
  if (CheckDemogFileStructure(mycountry=mycountry, mydf=dfbirth, dfdesc="cbr_both")==FALSE) { stop (filemsg)}
  if (CheckDemogFileStructure(mycountry=mycountry, mydf=numbirth, dfdesc="births")==FALSE) { stop (filemsg)}
  ctrybirth<-dfbirth[dfbirth$country_code==mycountry, c("country_code", "year", "value")]
  numbirth_ctry<-numbirth[numbirth$country_code==mycountry, c("country_code", "year", "value")]
  ctrybirthfull<-checkVIMCdates(mydata=ctrybirth, startyr=year(start), endyr=year(end), threshold=fillThreshold)
  numbirthfull<-checkVIMCdates(mydata=numbirth_ctry, startyr=year(start), endyr=year(end), threshold=fillThreshold)
  if (is.data.frame(ctrybirthfull)==FALSE) { stop(paste(datemsg, " cbr_both")) }
  if (is.data.frame(numbirthfull)==FALSE) { stop(paste(datemsg, " births")) }
  ctrybirthfull%>%group_by(country_code)%>%summarize(min(year), max(year))
  numbirthfull%>%group_by(country_code)%>%summarize(min(year), max(year))
 
  build0<-merge(x=ctrypopfull, y=ctrybirthfull, by=c("country_code", "year"), all=TRUE)
  colnames(build0)[colnames(build0)=="value.x"] <- "totalpop"
  colnames(build0)[colnames(build0)=="value.y"] <- "birthrate"
  build1 <- merge(x=build0, y=numbirthfull, by=c("country_code", "year"), all=TRUE)
  colnames(build1)[colnames(build1)=="value"] <- "births"
  
  # build1$births<-build1$totalpop*build1$birthrate
  # Chloe 7/12: replaced this crude calculation with WPP estimates, since they use more precise birth rates and further smoothing.
  
  # keep empty age_from and age_to from imr file to preserve format
  # currently infant mortality rates only here.
  imr<-GetFilename(path, "imr_both")
  if (is.character(imr)==FALSE) { stop(mymsg) }
  dfim<-read.csv(imr)
  if (CheckDemogFileStructure(mycountry=mycountry, mydf=dfim, dfdesc="imr_both")==FALSE) { stop (filemsg)}
  ctryimr<-dfim[dfim$country_code==mycountry, c("country_code", "year", "value")]
  ctryimrfull<-checkVIMCdates(mydata=ctryimr, startyr=year(start), endyr=year(end), threshold=fillThreshold)
  if (is.data.frame(ctryimrfull)==FALSE) { stop(paste(datemsg, " imr_both")) }
  ctryimrfull%>%group_by(country_code)%>%summarize(min(year), max(year))

  # build2<-merge(x=build1, y=ctryimrfull, by=c("country_code", "year"), all=TRUE)
  
  # currently other death rate here.
  # cdr<-GetFilename(path, "cdr_both")
  # if (is.character(cdr)==FALSE) { stop(mymsg) }
  # dfcdr<-read.csv(cdr)
  # if (CheckDemogFileStructure(mycountry=mycountry, mydf=dfcdr, dfdesc="cdr_both")==FALSE) { stop (filemsg)}
  # ctrycdr<-dfcdr[dfcdr$country_code==mycountry, c("country_code", "year", "value")]
  # ctrycdrfull<-checkVIMCdates(mydata=ctrycdr, startyr=year(start), endyr=year(end), threshold=fillThreshold)
  # if (is.data.frame(ctrycdrfull)==FALSE) { stop(paste(datemsg, " cdr_both")) }
  # ctrycdrfull%>%group_by(country_code)%>%summarize(min(year), max(year))
  
  # build3<-merge(x=build2, y=ctrycdrfull, by=c("country_code", "year"), all=TRUE)
  # colnames(build3)[colnames(build3)=="value.x"] <-"imr"
  # colnames(build3)[colnames(build3)=="value.y"] <-"v"
  
  # Chloe 5/22: now, getting all death rates across multiple age groups from single source.
  dr <- GetFilename(path, "p_dying_both")
  if (is.character(dr)==FALSE) { stop(mymsg) }
  dfcdr<-read.csv(dr)
  if (CheckDemogFileStructure(mycountry=mycountry, mydf=dfcdr, dfdesc="p_dying_both")==FALSE) { stop (filemsg)}
  ctrycdr<-dfcdr[dfcdr$country_code==mycountry, c("country_code", "age_to","year","value")]
  
  # Chloe 5/22: now, need to switch this to wide format for each year.
  ctrycdr.wide <- dcast(ctrycdr,country_code+year~age_to,value.var="value")
  curr.names <- as.numeric(colnames(ctrycdr.wide)[3:ncol(ctrycdr.wide)])
  for (i in 4:ncol(ctrycdr.wide)){
    colnames(ctrycdr.wide)[i] <- paste0("dr",curr.names[i-3]+1,curr.names[i-2])
  }
  colnames(ctrycdr.wide)[3] <- "imr"
  
  # Chloe 5/22: Since these death rates are only available every 5 years, I am assuming they are consistent across each 5-year span.
  # However, there is more detailed IMR data available: will retrieve and add later.
  ctrycdr.new <- ctrycdr.wide %>% slice(rep(1:n(), each = 5))
  for (i in 1:nrow(ctrycdr.wide)){
    new.years <- c(ctrycdr.wide$year[i]:(ctrycdr.wide$year[i]+4))
    first <- which(colnames(ctrycdr.new)=="year")
    ind1 <- ((i-1)*5)+1
    ind2 <- (((i-1)*5)+5)
    ctrycdr.new[ind1:ind2,first] <- new.years
  }
  # Replacing infant morality rates from more general file with ones from more-detailed file.
  # Applicable for the first year of age.
  # Note the final column only goes through age 84, so I'll assume the same death rate >84 yrs for now.
  ctryimrfull.part <- ctryimrfull[,-which(colnames(ctryimrfull)=="country_code")]
  ctrycdr.det <- merge(ctrycdr.new,ctryimrfull.part,by.x="year")
  ctrycdr.det$imr <- ctrycdr.det$value
  ctrycdr.det <- ctrycdr.det[,-which(colnames(ctrycdr.det)=="value")]
  
  ctrycdrfull<-checkVIMCdates(mydata=ctrycdr.det, startyr=year(start), endyr=year(end), threshold=fillThreshold)
  if (is.data.frame(ctrycdrfull)==FALSE) { stop(paste(datemsg, " p_dying_both")) }
  ctrycdrfull%>%group_by(country_code)%>%summarize(min(year), max(year))
  
  # Chloe 5/22: Making final data frame.
  build3<-merge(x=build1, y=ctrycdrfull, by=c("country_code", "year"), all=TRUE)
  
  return(build3)
  
}

checkVIMCdates<-function(mydata, startyr, endyr, threshold=1) {
  #assume data has variables country and year
  #will fill in up to a threshold (default = 1 year) with values from nearest year
  datemsg<<-""
  datesum<-mydata%>%dplyr::group_by(country_code)%>%dplyr::summarize(minyear=min(year), maxyear=max(year))
  if (datesum$minyear > startyr) {
    prediff<-datesum$minyear-startyr
    if (prediff <= threshold) {
      #fix by filling
      prefill<-mydata[mydata$year==datesum$minyear,]
      for (i in 1:prediff) {
        prefill$year<-datesum$minyear-i
        mydata<-rbind(mydata, prefill)
      }
    } else {
    # over threshold message
      datemsg<<- paste0("There is a ", prediff, "-year gap in the demographic data compared to your simulation begin date.  You can increase the threshold parameter to fill in data from the closest year.")
      return(FALSE)
    }
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
      } else {
        # over threshold message
        datemsg<<-paste0("There is a ", postdiff, "-year gap in the demographic data compared to your simulation end date.  You can increase the threshold parameter to fill in data from the closest year.")
        return(FALSE)
      } 
    }
  }
  return(mydata)
}

#_____________________________________________________________________________#
# This function returns a vector with 7 values, one for each 5-year age       #
#  band up to 30, calculated from quinquennial file, for the year closest to  #
#  specified start of simulation. Added names=age band 11/15/18               #
#  ASSUMES FILENAME CONTAINS "qq_pop_both"                                    #
#  Called by InitializePopulation.R                                           #
#_____________________________________________________________________________#
GetPopAgeDist<-function(path, mycountry, start) {
  disterr<<-""
  setwd(path)
  qqfile<-GetFilename(path, "qq_pop_both")
  if (is.character(qqfile)==FALSE) { stop(mymsg) }
  qqdf<-read.csv(qqfile)
  if (CheckDemogFileStructure(mycountry=mycountry, mydf=qqdf, dfdesc="qq_pop_both")==FALSE) { stop (filemsg)}
  #record for every 5th year - want the closest to start or before start?
  mround <- function(x,base){ 
    base*round(x/base) 
  }
  popyr<-mround(year(start), 5)
  qqkeep<-qqdf[qqdf$country_code==mycountry & qqdf$year==popyr,]
 
  if (nrow(qqkeep) > 0) {
    if (!(DemogNumVarExists("age_from", qqkeep) & DemogNumVarExists("age_to", qqkeep))) { 
      disterr<<- "Non-numeric data in age band variables"
      return(FALSE) 
    }
    if (!(DemogNumVarExists("age_to", qqkeep))) { return(FALSE) }
    #check min and max age band - now they all go up to 120 but lets be generous and say 90
    if (min(qqkeep$age_from > 0 || max(qqkeep$age_to) < 90) == FALSE) {
      # Chloe edit 3/22/19: Removing this line so as to not lump all 30+ groups together and replacing with next one.
      # qqkeep$ageband<-ifelse(qqkeep$age_from < 30, paste0("Age_", qqkeep$age_from, "_",qqkeep$age_to), "Age_30")
      qqkeep$ageband<-paste0("Age_", qqkeep$age_from, "_",qqkeep$age_to)
      bands<-qqkeep%>%group_by(country_code, ageband)%>%summarize(tot=sum(value),minage=min(age_from)) 
      totpop<-qqkeep%>%group_by(country_code)%>%summarize(totpop=sum(value))
      if (totpop[,2] > 0) {
        numsall<-merge(x=bands, y=totpop, by="country_code")
        numsall$fraction<-numsall$tot/numsall$totpop
        agedist <- numsall[order(numsall$minage),]
        dist<-agedist$fraction
        names(dist)<-agedist$ageband
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


#_____________________________________________________________________________#
# This function returns a dataset with a row for each year of simulation,     #
# with DosesCampaign, CoverRoutine, and AgeLimCampaign                        # 
#   NOTE: ASSUMES FILENAMES CONTAIN: "mena-routine" and "mena-campaign"       #
#    NOTE: NOT CURRENTLY LIMITED TO YEARS OF SIM, could use GetVIMCdates?     #
#     12/6/18 copying destring from taRifx to deal with "<NA>" in vacc files  #
#_____________________________________________________________________________#
destring <- function(x,keep="0-9.-") {
  return( as.numeric(gsub(paste("[^",keep,"]+",sep=""),"",x)) )
}
GetVaccScenario<-function(mycountry, scenario, directory) {
  vaccmsg<<-""
  setwd(directory)
  if (scenario=="routine" | scenario=="both") {
    filename<-GetFilename(directory, "mena-routine")
  }
  if (scenario=="campaign") {
    filename<-GetFilename(directory, "mena-campaign")
  }
  if (is.character(filename)==FALSE) { stop(mymsg) }
  dfvacc<-read.csv(filename, stringsAsFactors = FALSE)
  if (IsCountryAndColAvailable(country_code=mycountry,mydf=dfvacc, forVacc=1)==FALSE) { stop(countrymsg) }
  #target and year validated above.  Do we need AgeLimCampaign? No its not used.
  if (scenario=="routine" || scenario=="both") {
      if (!(DemogNumVarExists("coverage", dfvacc))) { 
        vaccmsg<<-"coverage variable missing from vaccination file"
        return(FALSE) 
        }
  }
  ctryvacc<-dfvacc[dfvacc$country_code==mycountry, c("country_code", "year", "activity_type", "target" , "coverage")]
  colnames(ctryvacc)[colnames(ctryvacc)=="coverage"] <-"CoverRoutine"
  ##target has "<NA>" where activity type = routine, hosing conversion
  #still getting coercion warning
  #getting this even though not strictly required by routine option
  ctryvacc$DosesCampaign<-destring(ctryvacc$target)
  newdf<-subset(ctryvacc, select=-c(target))
}


#_____________________________________________________________________________#
# Function GetDiseaseStateDist, called by InitializePopulation.R              #
# Reads dist_both.csv, which is supplied with scripts; format should not vary #
#_____________________________________________________________________________#
GetDiseaseStateDist<-function(directory, region) {
  setwd(directory)
  dxfile<-GetFilename(directory, "dist_both.csv")
  if (is.character(dxfile)==FALSE) { 
    stop(mymsg) 
    print("File [dist_both.csv] is packaged with the R scripts and should be in the same directory.")
  }
  dist<-read.csv(dxfile, stringsAsFactors = TRUE)
  distcol<-ifelse(region=='hyper', 4, 3)
  statefract<-as.vector(dist[,distcol]) # fraction of each disease state in each of 7 population groups
  return(statefract)
}

#_____________________________________________________________________________#                                                                  #
# Get WAIFWmatrix: Reads waifw_both.csv, which is supplied with scripts;      # 
# format should not vary                                                      #
#	expandWaifw : called by GetWAIFWmatrix: Expand WAIFW matrices               #
#      to match monthly population length, output is used by MenA_OneSim      #
#_____________________________________________________________________________#
# Chloe edit 3/29: need to expand further to account for higher ages part of sim now.
# Chloe edit 5/5: left all these new additions to the GetWAIFmatrix() function below.
expandWaifw<-function(waifw){
  # repeat what was originally columns :
  #b[x,1] 60x; b[x,2] 96x; b[1,3] 84x; b[1,4] 120x
  #needs to go to 361 - add extra line at end for last big bucket
  # Chloe edit 3/29: see note above.
  rbind ( 
    matrix(data=waifw[c(1,5,9,13)], nrow=60, ncol=4, byrow=TRUE),
    matrix(data=waifw[c(2,6,10,14)], nrow=96, ncol=4, byrow=TRUE),
    matrix(data=waifw[c(3,7,11,15)], nrow=84, ncol=4, byrow=TRUE),
    matrix(data=waifw[c(4,8,12,16)], nrow=121, ncol=4, byrow=TRUE)
  )
  
}

GetWAIFWmatrix<-function(path, region) {
  setwd(path)
  waifwfile<-GetFilename(path, "WAIFW_both.csv")
  if (is.character(waifwfile)==FALSE) { 
    stop(mymsg) 
    print("File [waifw_both.csv] is packaged with the R scripts and should be in the same directory.")
  }
  waifwin<-read.csv(waifwfile, stringsAsFactors = FALSE)  #vector
  Rwaifw<-waifwin[waifwin$region==region & waifwin$season=='rainy', 4]
  Dwaifw<-waifwin[waifwin$region==region & waifwin$season=='dry', 4]
  wboth<-array(c(expandWaifw(waifw=Rwaifw), expandWaifw(waifw=Dwaifw)), dim=c(361,4,2))
  # Chloe edit 3/29: Need to change dims to suit larger matrix.
  # Used to have one row per month up until age 30, then the same
  # value for ages 30 through 100; now need separate row for each month of age from 30 to the end.
  # i.e. repeating final row values.
  add.rainy <- array(rep(wboth[361,,1],each=(1441-361)),dim=c(1441-361,4,1))
  add.dry <- array(rep(wboth[361,,2],each=(1441-361)),dim=c(1441-361,4,1))
  rainy.new <- rbind(wboth[,,1],add.rainy[,,1])
  dry.new <- rbind(wboth[,,2],add.dry[,,1])
  wboth <- array(NA,dim=c(1441,4,2))
  wboth[,,1] <- rainy.new
  wboth[,,2] <- dry.new
  dimnames(wboth)[[3]]<-c("rainy", "dry")
  return(wboth)
}







