#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: MenA_paramCheck.R                                         #
# Contact: chris.c.stewart@kp.org, michael.l.jackson@kp.org                   #
# Version Date 12/31/2018                                                     #
#_____________________________________________________________________________#
# Functions called for validating inputs, called by MenA_VaccSims.R and       #
# ModelInputUtilities.R                                                       #
# ALL RETURN TRUE OR (FALSE + MESSAGE) except GetFilename, which returns a    #
#   filename or (FALSE + MESSAGE)                                             #
# Contents:                                                                   #
# -DemogNumVarExists: takes variable name and dataframe, checks to see if     #
#   variable exists in dataframe and that it is numeric                       #
# -GetFilename: takes directory path and pattern (part of file name.)         #
#   Makes sure a file is found and can be read with read.csv                  #
# -IsCountryAndColAvailable: takes country code and dataframe, makes sure df  #
#   contains country_code variable and data for requested country is present  #
#   if forVacc=1 is specified, also checks for target variable, which is not  #
#   currently numeric in source data                                          #
# -CheckDemogFileStructure: takes country, dataframe, and description         #
#   (which file, for error message only.)  Calls IsCountryandColAvailable     # 
#    makes sure numeric year and value variables are present.                 #
# -CheckDemogParamaters: input is output of GetDemogParameters (dataframe)    #
#   makes sure all variables are present and means are non-zero               # 
# -CheckSetParamaters: input is a list of parameters set by user in calling   #
#   script. Validates input: some are limited to specific values, others      #
#   must be numeric, directories must be valid path.                          #
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

IsCountryAndColAvailable<-function(country_code, mydf, forVacc=0) {
  #make sure requested country is in data, and required columns are present
  countrymsg<<-""
  if ("country_code" %in% colnames(mydf)==FALSE) {
    countrymsg<<-"country_code variable not found in input data:"
    return(FALSE)
  }
  if (!(country_code %in% levels(mydf$country_code) || country_code %in% as.vector(mydf$country_code))){
    countrymsg<<-paste0(country_code, " not available in input data:")
    return(FALSE)
  }
  if (forVacc==1) { #check target instead (its not numeric so don't use other column validation fxn)
    if ("target" %in% colnames(mydf)) {
      return(TRUE)
    } else {
      countrymsg<<-"Required column [target] is missing from input data:"
      return(FALSE)
    }
  }
  return(TRUE)
}

GetFilename<-function(path, pattern, subpattern="NA") {
  mymsg<<-""
  setwd(path)
  flist<-list.files(path)
  myfile<-flist[grepl(pattern,flist)==TRUE]
  if (length(myfile)==0) { 
    mymsg<<-paste("No file found for ", pattern)
    return(FALSE)
  }
  if (subpattern != "NA") {  #code to allow selection between bestcase and default vaccination scenarios
    myfile<-myfile[grepl(subpattern,myfile)==TRUE]
  }
  e<-try(read.csv(myfile[1]))
  if (class(e) == "try-error" || is.na(e)) { 
    mymsg<<-paste0(myfile[1], " is not a .csv file")
    return(FALSE)
  }
  return(myfile[1])
}

CheckDemogFileStructure<-function(mycountry, mydf, dfdesc){
  filemsg<<-""
  #just a wrapper/shortcut for 3 function calls
  if (IsCountryAndColAvailable(mycountry, mydf)==FALSE) {
    filemsg<<-paste(countrymsg, dfdesc)
    return(FALSE)
  }
  if (!(DemogNumVarExists("year", mydf))) { 
    filemsg<<-paste("Numeric variable [year] does not exist in file: ", dfdesc)
    return(FALSE)
  }
  if (!(DemogNumVarExists("value", mydf))) { 
    filemsg<<-paste("Numeric variable [value] does not exist in file: ", dfdesc) 
    return(FALSE)
  }
  return(TRUE)
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
    # if (!(DemogNumVarExists("v", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr14", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr59", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr1014", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr1519", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr2024", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr2529", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr3034", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr3539", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr4044", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr4549", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr5054", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr5559", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr6064", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr6569", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr7074", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr7579", params))) { return(FALSE) }
    if (!(DemogNumVarExists("dr8084", params))) { return(FALSE) }
                          
    
    # mns<-colMeans(params[, c("totalpop", "year","births", "imr", "v" )])
    mns<-colMeans(params[, c("totalpop", "year","births", "imr",  "dr14", "dr59", "dr1014","dr1519",
                             "dr2024" ,"dr2529" ,"dr3034" , "dr3539","dr4044", "dr4549","dr5054" ,
                             "dr5559" ,"dr6064" ,"dr6569" , "dr7074","dr7579","dr8084")])
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
      dpmessage<<-"Death rate for ages 1-4 are missing or zero."
      return(FALSE)
    }
    if (is.na(mns[6]) || mns[6]==0) {
      dpmessage<<-"Death rate for ages 5-9 are missing or zero."
      return(FALSE)
    }
    if (is.na(mns[7]) || mns[7]==0) {
      dpmessage<<-"Death rate for ages 10-14 are missing or zero."
      return(FALSE)
    }
    if (is.na(mns[8]) || mns[8]==0) {
      dpmessage<<-"Death rate for ages 15-19 are missing or zero."
      return(FALSE)
    }
    if (is.na(mns[9]) || mns[9]==0) {
      dpmessage<<-"Death rate for ages 20-24 are missing or zero."
      return(FALSE)
    }
  }
  if (is.na(mns[10]) || mns[10]==0) {
    dpmessage<<-"Death rate for ages 25-30 are missing or zero."
    return(FALSE)
  }
  if (is.na(mns[11]) || mns[11]==0) {
    dpmessage<<-"Death rate for ages 31-34 are missing or zero."
    return(FALSE)
  }
  if (is.na(mns[12]) || mns[12]==0) {
    dpmessage<<-"Death rate for ages 35-39 are missing or zero."
    return(FALSE)
  }
  if (is.na(mns[13]) || mns[13]==0) {
    dpmessage<<-"Death rate for ages 40-44 are missing or zero."
    return(FALSE)
  }
  if (is.na(mns[14]) || mns[14]==0) {
    dpmessage<<-"Death rate for ages 45-49 are missing or zero."
    return(FALSE)
  }
  if (is.na(mns[15]) || mns[15]==0) {
    dpmessage<<-"Death rate for ages 50-54 are missing or zero."
    return(FALSE)
  }
  if (is.na(mns[16]) || mns[16]==0) {
    dpmessage<<-"Death rate for ages 55-59 are missing or zero."
    return(FALSE)
  }
  if (is.na(mns[17]) || mns[17]==0) {
    dpmessage<<-"Death rate for ages 60-64 are missing or zero."
    return(FALSE)
  }
  if (is.na(mns[18]) || mns[18]==0) {
    dpmessage<<-"Death rate for ages 65-69 are missing or zero."
    return(FALSE)
  }
  if (is.na(mns[19]) || mns[19]==0) {
    dpmessage<<-"Death rate for ages 70-74 are missing or zero."
    return(FALSE)
  }
  if (is.na(mns[20]) || mns[20]==0) {
    dpmessage<<-"Death rate for ages 75-79 are missing or zero."
    return(FALSE)
  }
  if (is.na(mns[21]) || mns[21]==0) {
    dpmessage<<-"Death rate for ages 80-120 are missing or zero."
    return(FALSE)
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
  
