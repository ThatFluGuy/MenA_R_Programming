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
# (1) Set up options, get functions                                           #
# (2) Compile central estimates                                               #
#_____________________________________________________________________________#
# Author: Mike Jackson;  michael.l.jackson@kp.org                             #
#_____________________________________________________________________________#

### (1) Set up options, get functions #########################################

library(dplyr)

output.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Analysis/Simulation results"
deliv.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Deliverables/Deliverables 2021"
script.dir <- "C:/Users/O992928/documents/GAVI MenA predictions/R_programming"

setwd(script.dir)
source("fxModelInputs.R")

### (2) Compile central estimates #############################################
combineOutputFiles(path=output.dir, vacc_program="none", 
                   vacc_subprogram="default", deliv.path=deliv.dir, touchstone="202108")

combineOutputFiles(path=output.dir, vacc_program="routine", 
                   vacc_subprogram="default", deliv.path=deliv.dir, touchstone="202108")

combineOutputFiles(path=output.dir, vacc_program="booster", 
                   vacc_subprogram="default", deliv.path=deliv.dir, touchstone="202108")






