#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: MenA_CompileOutputs.R                                     #
# Version Date 01/13/2020                                                     #
#_____________________________________________________________________________#
# Purpose: Take simulation output files and convert to VIMC upload format.    #
# For central estimates, combine files across countries. For PSA, drop the    #
# row numbers and rename.                                                     #
#_____________________________________________________________________________#
# Inputs: Individual simulation results for each country/scenario.            #
#_____________________________________________________________________________#
# Outputs: For central estimates, combined .csv files across countries. For   #
# PSA, renamed files.                                                         #
#_____________________________________________________________________________#
# Steps in this program:                                                      #
#_____________________________________________________________________________#
# Author: Mike Jackson;  michael.l.jackson@kp.org                             #
#_____________________________________________________________________________#
# FLAG FOR FUTURE WORK: Fix MenA_VaccSims.R so that row.names are not output  #
# for PSA, and so that names are in VIMC suggested format. Then the PSA part  #
# of this program will not be necessary.                                      #
#_____________________________________________________________________________#


### (1) Set up options ########################################################

output.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Analysis/Simulation results"
deliv.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Deliverables/Deliverables 2020"

### (2) Compile central estimates #############################################

for (s in 1:9){
  combineOutputFiles(path=output.dir, vacc_program=paste0("scenario", s), deliv.path=deliv.dir)
  
}


### (3) Clean up PSA files ####################################################
# Need to rename according to VIMC's preferred convention. Also delete the    #
# column "x" that contains row numbers.                                       #

var.names <- c("disease", "run_id", "year", "age", "country", "country_name",
               "cohort_size",	"cases", "dalys", "deaths")

# (A) Get names of PSA files in output.dir
files.v <- list.files(output.dir, pattern="PSA_")

# (B) Pull out country and vaccine program codes
country_code.v <- sapply(strsplit(files.v, split="_"), '[', 2)
vacc_program.v <- sapply(strsplit(files.v, split="_"), '[', 3)
vacc_subprogram.v <- sapply(strsplit(files.v, split="_"), '[', 4)
vacc.v <- paste(vacc_program.v, vacc_subprogram.v, sep="_")

# (C) Index scenario files by country and create new names
country_num.v <- match(country_code.v, unique(country_code.v))

new.files.v <- paste("stochastic_burden_est_MenA_KPWA_", vacc.v, "_", 
                     country_num.v, ".csv", sep="")
  
for (f in 1:length(files.v)){
  print(paste("Processing file #", f, files.v[f], sep=" "))
  temp.df <- read.csv(paste(output.dir, files.v[f], sep="/"),
                      stringsAsFactors = FALSE)
  write.csv(x=temp.df[, names(temp.df) %in% var.names], 
            file=paste(deliv.dir, new.files.v[f], sep="/"),
            row.names=FALSE)
  rm(temp.df)
  gc()
}


