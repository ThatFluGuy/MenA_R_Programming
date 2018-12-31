#_____________________________________________________________________________#
# Functions called for validating inputs, called by MenA_VaccSims.R and       #
# ModelInputUtilities.R                                                       #
#_____________________________________________________________________________#
DemogNumVarExists<-function(varname, mydf) {
  colix<-grep(varname, colnames(mydf))
  if (length(colix > 0)) {
    if (is.numeric(mydf[,colix])) {
      return(TRUE)
    } else {
      dpmessage<<-paste0(varname, " in demographic parameters is not numeric.")
      return(FALSE)
    }
  } else {
    dpmessage<<-paste0(varname, " variable missing from demographic parameters.")
    return(FALSE)
  }
  return(TRUE)
}

ValidateFilename<-function(mypath, myfile) {
  #input dir has already been validated
  mymsg=""
  if (grepl(".csv", myfile)==FALSE) {
    mymsg<<-paste0(myfile, " is not a .csv file")
    return(FALSE)
  }
  if (file.exists(paste0(mypath, "/", myfile))==FALSE){
    mymsg<<-paste0(myfile, " not found in input directory")
    return(FALSE)
  }
  return(TRUE)
}

#_____________________________________________________________________________#
# This function checks for the country code in the data, and checks for the   #
# other reqd variable names (year and value, or year and target if forVacc)   #
#_____________________________________________________________________________#
IsCountryAndColsAvailable<-function(country_code, mydf, forVacc=0) {
  #make sure requested country is in data, and required columns are present
  countrymsg<<-""
  if ("country_code" %in% colnames(mydf)==FALSE) {
    countrymsg<<-"country_code variable not found in input data."
    return(FALSE)
  }
  if (country_code %in% levels(mydf$country_code) || country_code %in% as.vector(mydf$country_code)) {
    return(TRUE)
  } else {
    countrymsg<<-paste0(country_code, " not available in input data.")
    return(FALSE)
  }
  if (forVacc==0) {
    if ("year" %in% colnames(mydf) & "value" %in% colnames(mydf)) {
      return(TRUE)
    } else {
      countrymsg<<-"Required columns [year] and [value] are missing from input data."
      return(FALSE)
    }
  } else { #check target instead (its not numeric so don't use other column validation fxn)
    if ("year" %in% colnames(mydf) & "target" %in% colnames(mydf)) {
      return(TRUE)
    } else {
      countrymsg<<-"Required columns [year] and [target] are missing from input data."
      return(FALSE)
    }
  }
}

CheckDemogParameters<-function(params) {
  dpmessage<<-""
  if (is.data.frame(params)) {
    if (nrow(params)==0) {
      dpmessage<<-"Demographic parameter set is empty"
      return(FALSE)
    }
    if (!(DemogNumVarExists("totalpop", params))) { return(FALSE) }
    if (!(DemogNumVarExists("year", params))) { return(FALSE) }
    if (!(DemogNumVarExists("births", params))) { return(FALSE) }
    if (!(DemogNumVarExists("imr", params))) { return(FALSE) }
    if (!(DemogNumVarExists("v", params))) { return(FALSE) }
    
    mns<-colMeans(params[, c("totalpop", "year","births", "imr", "v" )])
    if (is.na(mns[1]) || mns[1]==0) {
      dpmessage<<-"Values for totalpop are missing or zero."
      return(FALSE)
    }
    if (is.na(mns[2]) || mns[2]==0) {
      dpmessage<<-"Values for year are missing or zero."
      return(FALSE)
    }
    if (is.na(mns[3]) || mns[3]==0) {
      dpmessage<<-"Values for births are missing or zero."
      return(FALSE)
    }
    if (is.na(mns[4]) || mns[4]==0) {
      dpmessage<<-"Values for infant mortality rate (imr) are missing or zero."
      return(FALSE)
    }
    if (is.na(mns[5]) || mns[5]==0) {
      dpmessage<<-"Values for death rate (v) are missing or zero."
      return(FALSE)
    }
  }
  return(TRUE)
}

CheckSetParameters<-function(setparams) {
  spmessage<<-""
  if (!(setparams$myregion %in% c("hyper", "not_hyper"))) {
    spmessage<<-"myregion must be set to either hyper or not_hyper."
    return(FALSE)
  }
  if (!(setparams$vacc_program %in% c("campaign", "routine", "both", "none"))) {
    spmessage<<-"Valid values for vacc_program are campaign, routine, both, or none."
    return(FALSE)
  }
  #start and end must be dates, end>start
  s<-try(as.Date(setparams$start, format="%Y-%m-%d"))
  e<-try(as.Date(setparams$end, format="%Y-%m-%d"))
  if (class(s) == "try-error" || is.na(s) || class(e) == "try-error" || is.na(e)) { 
    spmessage<<-"Start and end must be valid dates"
    return(FALSE)
  }
  else {
    if (as.Date(setparams$end) <= as.Date(setparams$start)) {
      spmessage<<-"End date must be greater than start date"
      return(FALSE)
      }
  }
  if (PSA != FALSE) {
    PSA<-FALSE
    spmessage<<-"Sorry, only PSA=FALSE is available at this time."
  }
  x<-suppressWarnings(try(as.numeric(setparams$phi)))
  y<-suppressWarnings(try(as.numeric(setparams$sd)))
  z<-suppressWarnings(try(as.numeric(setparams$nSims)))
  if (class(x) == "try-error" || is.na(x) || class(y) == "try-error" || is.na(y) || class(z) == "try-error" || is.na(z)) {
    spmessage<<-"phi, sd, and nSims must all be numeric."
    return(FALSE)
  }
  else {
    if (phi < 0 || phi > 1) {
      spmessage<<-"phi must be between 0 and 1"
      return(FALSE)
    }
    #warn if very large?
    if (nSims <=0) {
      spmessage<<-"nSims must be greater than 0"
      return(FALSE)
    }
  }
  if (dir.exists(inputdir)==FALSE) {
    spmessage<<-"inputdir is not a valid directory"
    return(FALSE)
  }
  if (dir.exists(outputdir)==FALSE) {
    spmessage<<-"outputdir is not a valid directory"
    return(FALSE)
  }
  return(TRUE)
}
  
