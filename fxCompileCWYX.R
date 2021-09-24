#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: fxCompileCWYX.R                                           #
# Version Date 09/16/2021                                                     #
#_____________________________________________________________________________#
# Purpose: This program adds predicted cases/deaths/DALYs due to serogroup    #
# cwyx meningococcal disease to the output from the serogroup A simulations.  #
# The program [program name] creates the expected number of cases in each     #
# country/year/age group combination.                                         #
# This program takes those case counts, splits by year of age, estimates      #
# deaths and DALYs, and merges with the NmA output.                           #
#_____________________________________________________________________________#
# Input: MenCWYX_Sims_ccc.csv, where ccc is the country abbreviation.         #
#_____________________________________________________________________________#
# Output: A list with averaged CWYX cases and with individual model runs.     #
#_____________________________________________________________________________#
# Author: Mike Jackson                                                        #
#_____________________________________________________________________________#

compileCWYX <- function(cwyx.country, cwyx.cohortSize=cohortSize, 
                        cwyx.lifeEx=my.lifex, cwyx.scenario=vacc_program,
                        cwyx.sub.scenario=vacc_subprogram,
                        cwyx.cfr=paramfixed %>% select(run_id, cfr1, cfr2, cfr3, cfr4, cfr5, cfr6),
                        cwyx.path="G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Data/ACWXY",
                        cwyx.vacc.path=input.dir){
  
  ### (1) Import data, convert to cases by year of age ##########################
  # If using PSA, do this separately for each iteration. If not, average across #
  # iterations.                                                                 #
  # Within the age blocks given by the cwyx simulation output, distribute the   #
  # case counts according to population size.                                   #
  
  cwyx.df <- read.csv(paste0(cwyx.path, "/Mencwyx_Sims_", cwyx.country, ".csv"), 
                        stringsAsFactors = FALSE)
    
    
  # Convert cwyx.df from wide to tall and add minimim age for each group.
  # Warning is expected due to some age groups (AgeLT5 eg) not being convertable
  # to numeric when taking substring(x, 4, 5). But these are already assigned to
  # numeric earlier in the ifelse.
  cwyx.tall <- pivot_longer(cwyx.df, cols=starts_with("Age"), names_to="AgeGrp",
                            values_to="total.cases")
  cwyx.tall$min.age <- ifelse(cwyx.tall$AgeGrp=="AgeLT5", 0, 
                              ifelse(cwyx.tall$AgeGrp=="Age5t9", 5, 
                                     as.numeric(substr(cwyx.tall$AgeGrp, 4,5))))
  
  # Create matching size distribution for cohortSize
  cwyx.cohortSize$min.age <- ifelse(cwyx.cohortSize$AgeInYears<20, 
                               (cwyx.cohortSize$AgeInYears %/% 5)*5,
                               ifelse(cwyx.cohortSize$AgeInYears>=50, 50, 
                                      (cwyx.cohortSize$AgeInYears %/% 10)*10))
  
  # Get total population in each cohort group
  cohortGroup <- aggregate(cohortsize~year+min.age, data=cwyx.cohortSize, FUN=sum)
  names(cohortGroup)[3] <- "total.pop"
  
  # Merge group populations with cohort sizes
  cohortGroup <- left_join(cwyx.cohortSize, cohortGroup, by=c("year", "min.age"))
  cohortGroup$percent <- cohortGroup$cohortsize/cohortGroup$total.pop
  
  # Merge cases into cohortGroup
  cohortGroup <- full_join(cohortGroup, cwyx.tall %>% select(-AgeGrp), 
                           by=c("year", "min.age"))
  cohortGroup$Cases <- cohortGroup$total.cases * cohortGroup$percent
  
  # Keep only variables that are needed
  cohortGroup <- cohortGroup %>% select(index, year, AgeInYears, Cases) %>%
    arrange(index, year, AgeInYears)
  
  
  ### (2) Implement vaccination, if any #######################################
  # Assume protection is fully effective for 10 years and then absent         #
  
  if (cwyx.scenario=="booster"){
    filename <- GetFilename(cwyx.vacc.path, "mena-booster", cwyx.sub.scenario)
    if (is.character(filename)==FALSE) { stop(mymsg) }
    
    dfvacc <- read.csv(filename, stringsAsFactors = FALSE)
    if (IsCountryAndColAvailable(country_code=cwyx.country, mydf=dfvacc, forVacc=1)==FALSE) {
      stop(countrymsg) 
    }
    
    # Restrict to Acwyx entries
    dfvacc <- dfvacc[dfvacc$vaccine!="MenA",]
    
    # Restrict to modeled country
    dfvacc <- dfvacc[dfvacc$country_code==cwyx.country,]
    
    # Expand to coverage by year and age group
    cwyx.vacc <- data.frame(year=numeric(0), AgeInYears=numeric(0), coverage=numeric(0))
    for(y in unique(dfvacc$year)){
      # Create a dataframe for the modeled year
      temp.df <- dfvacc %>% filter(year==y) %>% 
        select(year, age_first, age_last, coverage)
      age.length <- temp.df$age_last - temp.df$age_first + 1
      temp.df2 <- data.frame(year = rep(y, times=age.length), 
                             AgeInYears = seq(temp.df$age_first, temp.df$age_last, 1),
                             coverage = rep(temp.df$coverage, times=age.length))
      
      cwyx.vacc <- rbind(cwyx.vacc, temp.df2)
      # Make the next nine years
      for (x in 1:9){
        temp.df2[,1:2] <- temp.df2[,1:2] + 1
        cwyx.vacc <- rbind(cwyx.vacc, temp.df2)
      }
      rm(temp.df, temp.df2)
    }
    rm(filename, dfvacc)
    
    # Convert coverage to effective coverage, assuming 95% VE
    cwyx.vacc$coverage <- cwyx.vacc$coverage * 0.95
    
    # Merge into full data
    cohortGroup <- left_join(cohortGroup, cwyx.vacc, by=c("year", "AgeInYears"))
    cohortGroup$coverage[is.na(cohortGroup$coverage)==TRUE] <- 0
    cohortGroup$Cases <- cohortGroup$Cases * (1 - cohortGroup$coverage)
    
    cohortGroup <- cohortGroup %>% select(-coverage)
  }
  
  
  ### (3) Add CFR and life expectancy to get deaths and DALYS ###################
  
  cwyx.cfr <- cwyx.cfr %>% rename(index=run_id)
  
  # Calculate deaths 
  deaths.df <- left_join(cohortGroup, cwyx.cfr, by="index")
  deaths.df$dth.rate <- ifelse(deaths.df$AgeInYears==0, deaths.df$cfr1,
                             ifelse(deaths.df$AgeInYears<5, deaths.df$cfr2,
                                    ifelse(deaths.df$AgeInYears<10, deaths.df$cfr3,
                                           ifelse(deaths.df$AgeInYears<15, deaths.df$cfr4,
                                                  ifelse(deaths.df$AgeInYears<20, deaths.df$cfr5, 
                                                         deaths.df$cfr6)))))
  deaths.df$Deaths <- deaths.df$dth.rate * deaths.df$Cases
  
  deaths.df <- deaths.df %>% select(index, year, AgeInYears, Cases, Deaths)
  
  # Calculate DALYs
  daly.df <- left_join(deaths.df, cwyx.lifeEx, by=c("year", "AgeInYears"))
  daly.df$DALYs <- daly.df$Deaths * daly.df$Life.Ex +
    (daly.df$Cases - daly.df$Deaths) * (0.26 * 0.072) * daly.df$Life.Ex
  
  
  ### (4) Compile output ######################################################
  
  # Aggregated across iterations
  final.df1 <- aggregate(cbind(Cases, Deaths, DALYs)~year+AgeInYears, data=daly.df,
                          FUN=mean)
  final.df1 <- final.df1 %>% rename(cases_cwyx=Cases, deaths_cwyx=Deaths,
                                    dalys_cwyx=DALYs)
  # Not aggregated
  final.df2 <- daly.df %>% select(index, year, AgeInYears, Cases, Deaths, DALYs) %>%
    arrange(index, year, AgeInYears) 
  final.df2 <- final.df2 %>% rename(cases_cwyx=Cases, deaths_cwyx=Deaths,
                                    dalys_cwyx=DALYs, run_id=index, age=AgeInYears) %>% 
    as.data.frame()

  return(list(final.df1, final.df2)) 
}
