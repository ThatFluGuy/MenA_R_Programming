#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: Compile_MenCWXY.R                                         #
# Version Date 09/16/2021                                                     #
#_____________________________________________________________________________#
# Purpose: This program adds predicted cases/deaths/DALYs due to serogroup    #
# CWXY meningococcal disease to the output from the serogroup A simulations.  #
# The program [program name] creates the expected number of cases in each     #
# country/year/age group combination.                                         #
# This program takes those case counts, splits by year of age, estimates      #
# deaths and DALYs, and merges with the NmA output.                           #
#_____________________________________________________________________________#
# Input: MenCWXY_Sims_ccc.csv, where ccc is the country abbreviation.         #
#_____________________________________________________________________________#
# Output:                                                                     #
#_____________________________________________________________________________#
# Steps in this program:                                                      #
# (1) Set up libraries                                                        #
# (4) Run simulations                                                         #
#_____________________________________________________________________________#
# Author: Mike Jackson                                                        #
#_____________________________________________________________________________#

compileCWXY <- function(cwxy.country, cwxy.cohortSize=cohortSize, 
                        cwxy.lifeEx=my.lifex,
                        cwxy.cfr=paramfixed %>% select(run_id, cfr1, cfr2, cfr3, cfr4, cfr5, cfr6),
                        cwxy.path="G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Data/ACWXY"){
  
  ### (1) Import data, convert to cases by year of age ##########################
  # If using PSA, do this separately for each iteration. If not, average across #
  # iterations.                                                                 #
  # Within the age blocks given by the CWXY simulation output, distribute the   #
  # case counts according to population size.                                   #
  
  cwxy.df <- read.csv(paste0(path.cwxy, "/MenCWXY_Sims_", country, ".csv"), 
                        stringsAsFactors = FALSE)
    
    
  # Convert cwxy.df from wide to tall and add minimim age for each group.
  # Warning is expected due to some age groups (AgeLT5 eg) not being convertable
  # to numeric when taking substring(x, 4, 5). But these are already assigned to
  # numeric earlier in the ifelse.
  cwxy.tall <- pivot_longer(cwxy.df, cols=starts_with("Age"), names_to="AgeGrp",
                            values_to="total.cases")
  cwxy.tall$min.age <- ifelse(cwxy.tall$AgeGrp=="AgeLT5", 0, 
                              ifelse(cwxy.tall$AgeGrp=="Age5t9", 5, 
                                     as.numeric(substr(cwxy.tall$AgeGrp, 4,5))))
  
  # Create matching size distribution for cohortSize
  cwxy.cohortSize$min.age <- ifelse(cwxy.cohortSize$AgeInYears<20, 
                               (cwxy.cohortSize$AgeInYears %/% 5)*5,
                               ifelse(cwxy.cohortSize$AgeInYears>=50, 50, 
                                      (cwxy.cohortSize$AgeInYears %/% 10)*10))
  
  # Get total population in each cohort group
  cohortGroup <- aggregate(cohortsize~year+min.age, data=cwxy.cohortSize, FUN=sum)
  names(cohortGroup)[3] <- "total.pop"
  
  # Merge group populations with cohort sizes
  cohortGroup <- left_join(cwxy.cohortSize, cohortGroup, by=c("year", "min.age"))
  cohortGroup$percent <- cohortGroup$cohortsize/cohortGroup$total.pop
  
  # Merge cases into cohortGroup
  cohortGroup <- full_join(cohortGroup, cwxy.tall %>% select(-AgeGrp), 
                           by=c("year", "min.age"))
  cohortGroup$Cases <- cohortGroup$total.cases * cohortGroup$percent
  
  # Keep only variables that are needed
  cohortGroup <- cohortGroup %>% select(index, year, AgeInYears, Cases) %>%
    arrange(index, year, AgeInYears)
  
  ### (2) Add CFR and life expectancy to get deaths and DALYS ###################
  
  cwxy.cfr <- cwxy.cfr %>% rename(index=run_id)
  
  # Calculate deaths 
  deaths.df <- left_join(cohortGroup, cwxy.cfr, by="index")
  deaths.df$dth.rate <- ifelse(deaths.df$AgeInYears==0, deaths.df$cfr1,
                             ifelse(deaths.df$AgeInYears<5, deaths.df$cfr2,
                                    ifelse(deaths.df$AgeInYears<10, deaths.df$cfr3,
                                           ifelse(deaths.df$AgeInYears<15, deaths.df$cfr4,
                                                  ifelse(deaths.df$AgeInYears<20, deaths.df$cfr5, 
                                                         deaths.df$cfr6)))))
  deaths.df$Deaths <- deaths.df$dth.rate * deaths.df$Cases
  
  deaths.df <- deaths.df %>% select(index, year, AgeInYears, Cases, Deaths)
  
  # Calculate DALYs
  daly.df <- left_join(deaths.df, cwxy.lifeEx, by=c("year", "AgeInYears"))
  daly.df$DALYs <- daly.df$Deaths * daly.df$Life.Ex +
    (daly.df$Cases - daly.df$Deaths) * (0.26 * 0.072) * daly.df$Life.Ex
  
  
  ### (3) Compile output ######################################################
  
  # Aggregated across iterations
  final.df1 <- aggregate(cbind(Cases, Deaths, DALYs)~year+AgeInYears, data=daly.df,
                          FUN=mean)
  final.df1 <- final.df1 %>% rename(cases_cwxy=Cases, deaths_cwxy=Deaths,
                                    DALYs_cwxy=DALYs)
  # Not aggregated
  final.df2 <- daly.df %>% select(index, year, AgeInYears, Cases, Deaths, DALYs) %>%
    arrange(index, year, AgeInYears) 
  final.df2 <- final.df2 %>% rename(cases_cwxy=Cases, deaths_cwxy=Deaths,
                                    DALYs_cwxy=DALYs, run_id=index, age=AgeInYears) %>% 
    as.data.frame()

  return(list(final.df1, final.df2)) 
}
