#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: Create_Parameter_Sets.R                                   #
# Version Date 12/23/2019                                                     #
#_____________________________________________________________________________#
# PUrpose: VIMC requires that PSA model runs use a minimum of 200 iterations  #
# and with reproducable parameters. This program takes the posterior          #
# distributions of the parameter estimates from Jackson PLoS One 2017,        #
# expands to 200 sets via sampling with replacement to make up the difference,#
# and adds randomly sampled values for case fatality ratio (CFR) and for VE.  #
# CFR sampling is as in the SAS model, random draws from binomial distrib     #
# based on the sample size of the relevant study(ies). At present, VE is      #
# assumed to range from ~0.85 to ~0.95, mean 0.9. Using binomial with draw of #
# 500 approximates this.                                                      #
#_____________________________________________________________________________#
# Input: Script builds some inputs from files downloaded                      #
# from vaccineimpact.org/montagu.  Location of the files is an input parameter#
#_____________________________________________________________________________#
# Output: csv file with defined parameter sets.                               #
#_____________________________________________________________________________#
# Created 23 December 2019, by Mike Jackson.                                  #
#_____________________________________________________________________________#


### (1) Import and expand the posterior parameters ############################

script.dir <- "C:/Users/jackml4/documents/Link_to_H_drive/GAVI MenA predictions/R_programming"

parm.df <- read.csv(paste(script.dir, "/posterior_parameters_SAS_original.csv", sep=""), 
                    stringsAsFactors = FALSE)

# Randomly identify rows to resample to bring the parameter sets to 200
sample.v <- sample(1:length(parm.df[,1]), 200-length(parm.df[,1]), replace=TRUE)
parm.df2 <- parm.df[sample.v,]
parm.df3 <- rbind(parm.df, parm.df2)

parm.df3$run_id <- 1:200


### (2) Create additional parameters ##########################################

extra.df <- data.frame(run_id=1:200, cfr1=numeric(200), cfr2=numeric(200), 
                     cfr3=numeric(200), cfr4=numeric(200), cfr5=numeric(200), 
                     cfr6=numeric(200), ve=numeric(200))

set.seed(2019)
extra.df$cfr1 <- rbinom(200, 220, 0.1057)/200
extra.df$cfr2 <- rbinom(200, 880, 0.096)/880
extra.df$cfr3 <- rbinom(200, 1100, 0.0887)/1100
extra.df$cfr4 <- rbinom(200, 1100, 0.0861)/1100
extra.df$cfr5 <- rbinom(200, 1100, 0.0789)/1100
extra.df$cfr6 <- rbinom(200, 6600, 0.1217)/6600

extra.df$ve <- rbinom(200, 500, 0.9)/500

### (3) Merge, arrange, and output ############################################

parm.df4 <- merge(parm.df3[,!(names(parm.df3)=="ve")], extra.df, by="run_id")

write.csv(parm.df4, paste(script.dir, "/posterior_parameters.csv", sep=""))