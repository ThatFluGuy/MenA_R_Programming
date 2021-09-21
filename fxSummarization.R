#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: fxSummarization.R                                         #
# Contact: chris.c.stewart@kp.org, michael.l.jackson@kp.org                   #
# Version Date 08/20/21                                                       #
#_____________________________________________________________________________#
# This program contains a series of custom functions for summarizing the      #
# output of multiple simulation runs to get means and population sizes.       #
#_____________________________________________________________________________#
# Functions in this program:                                                  #
# 1) getCohortSize: calulate total poplation in each year and year of age     #
# 2) summarizeOneSim: Collapse results of one simulation into a data frame    #
#   with Cases, Deaths, Dalys by year and year of age.                        #
#	3) SummarizeForOutput : takes a list of products of summarizeOneSim and     #  
#   calculates mean over all simulations; writes output file                  #
#_____________________________________________________________________________#
# These functions were written by Chris Stewart, with some editing by Chloe   #
# Krakhauer and Mike Jackson                                                  #
#_____________________________________________________________________________#

### (1) Function to calculate cohort size #####################################
getCohortSize <- function(poparray) { 
  
  # Inputs: poparray - array with population by age, model compartment, and week
  # Outputs: data.frame with cohort sizes
  
  # Use the midpoint of the year (last week of June) to get population
  # To speed up processing, first restrict the full cohort array to this week
  dates <- as.Date(as.numeric(dimnames(poparray)[[3]]), origin="1970-01-01")
  pop.a <- poparray[, , month(dates)==6 & day(dates)>=24 & year(dates)>=2000]
  
  # Sum the population across model states and convert to tall
  cohort <- apply(pop.a[ , 1:9, ], c(1,3), sum) # only a little slow / EJ edit: Va -> Vs, Vc.  Column 10 is Inc(idence), not used here.
  cohortlong <- melt(cohort)
  cohortlong$RealDate <- as.Date(cohortlong[,2], origin="1970-01-01")
  cohortlong$year <- year(cohortlong$RealDate)
  cohortlong$AgeInYears <- floor((cohortlong[,1]-1)/12)
  
  # Group all age 70+ together
  cohortlong$AgeInYears <- ifelse(cohortlong$AgeInYears > 70, 70, cohortlong$AgeInYears)

  # Sum by year of age
  cohortsizes <- cohortlong %>%
    group_by(year, AgeInYears) %>%
    summarize(cohortsize=sum(value))
  
  return(cohortsizes)
}


### (2) Count cases, deaths, and DALYs ########################################
summarizeOneSim <- function(poparray, sim.number=n, cfr.v=cfr, le.df=my.lifex) { 
  
  # Inputs:
  # poparray - array with population by age, model compartment, and week
  # sim.n    - scaler indicating simulation iteration
  # cfr.v    - vector of age-specific case fatality ratios
  # le.df    - data.frame with age- and year-specific life expectancy
  
  # Output: a data.frame with year, age, simulation, cases, deaths, and DALYs
  
  #summarize incident cases by year and year of age, calculate deaths and DALYs
  inclong <- melt(poparray[,"Inc",], varname=c("AgeWeek", "TimeWeek"))
  inclong$RealDate <- as.Date(inclong[,2], origin="1970-01-01")
  #summarize incident cases by year, group all 70+ together
  inclong$year <- year(inclong$RealDate)
  inclong$AgeInYears <- floor((inclong[,1]-1)/12)
  inclong$AgeInYears <- ifelse(inclong$AgeInYears > 70, 70, inclong$AgeInYears) 
  # Sum by age and year; much faster than aggregate()
  res <- inclong %>% 
    filter(year>=2000) %>% group_by(year, AgeInYears) %>% 
    summarize(Cases=sum(value))

  # Vector of case fatality ratios
  cfr.age <- c(cfr.v[1], rep(cfr.v[2], 4), rep(cfr.v[3], 5), rep(cfr.v[4], 5),
               rep(cfr.v[5], 5), rep(cfr.v[6], 51))
  
  res$Deaths <- res$Cases*cfr.age

  # Merge in life expectancy and calculate DALYs
  results <- merge(res, le.df, by=c("year", "AgeInYears"))
  results$DALYs <- results$Deaths*results$Life.Ex +
    (results$Cases - results$Deaths) * (0.26 * 0.072) * results$Life.Ex
  
  #results$DALYs<-results$Deaths*(70-results$AgeInYears) + (results$Cases-results$Deaths)* (0.26*0.072)*(70-results$AgeInYears)
  results$simulation <- sim.number
  return(results[order(results$year, results$AgeInYears), names(results) != "Life.Ex"])
}

### (3) Average results across sims ###########################################

summarizeForOutput <- function(results_list, cohort, write, filename, cwxy=cwxy.l[[1]]) {
  if (length(results_list) > 1) {
    allsims <- rbindlist(results_list)
    simmeans <- allsims %>% 
      select(year, AgeInYears, Cases, Deaths, DALYs) %>% 
      group_by(year, AgeInYears) %>% 
      summarize_all(.funs=(mean))
    simcount <- allsims %>% 
      group_by(year, AgeInYears) %>%
      summarize(simulations=max(simulation))
    simsummary <-merge(x=simmeans, y=simcount, by = c("year", "AgeInYears"))
  } else {
    allsims<-results_list
    simsummary<-allsims[, c("year", "AgeInYears", "Cases", "Deaths", "DALYs")]
    simsummary$simulations<-1
  }
  
  # Merge in cohort size and CWXY results
  finalsummary <- merge(x = simsummary, y = cohort, by = c("year", "AgeInYears"))
  finalsummary <- left_join(x=finalsummary, y=cwxy, by = c("year", "AgeInYears")) %>%
    arrange(year, AgeInYears)
  
  if (write==TRUE){
    write.csv(finalsummary, filename, row.names=FALSE)
    print(paste("Output written to", filename))
  }
  return(finalsummary)
}