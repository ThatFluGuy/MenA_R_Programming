### Program information #######################################################
# Package: MenA_VaccSims                                                      #
# Program name: MenA_OutputDataChecks.R                                       #
# Version date: 14 January 2020                                               #
#_____________________________________________________________________________#
#  Input datasets: CCC_scenario_subscenario_date.csv, where CCC is a country  #
# abbreviation, scenario is the vaccine scenario, subscenario is the vaccine  #
# subscenario, and date is the date of the simulation.                        #
#_____________________________________________________________________________#
# Program description:                                                        #
# This program runs some basic checks on simulation output from the VIMC      #
# forecast outputs. Steps in this program:                                    #
# 1) Set up directories and packages                                          #
# 2) Import central data                                                      #
# 3) Plot deaths by scenario                                                      #
# 4) Plot deaths by age and year                                              #
# 5) Plot incidence by year and age                                           #
# 6) Vaccine impact                                                           #
#_____________________________________________________________________________#


### (1) Import packages, set up directories and datasets ######################
library(tibble)
library(haven) # Used to convert a SAS .sas7bdat file to an R tibble
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(dplyr)

fill.colors <- brewer.pal(11, name='Spectral')[c(1:4,8:11)]

output.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Analysis/Simulation results/"

type.list <- c("both_bestcase", "both_default", "campaign_bestcase", "campaign_default", 
               "none_default")


### (2) Import simulation data ################################################

# Create empty tibble as starting point
toPlot <- tibble(country=character(0), type=character(0), year=numeric(0), AgeInYears=numeric(0), 
                 Cases=numeric(0), Deaths=numeric(0), DALYs=numeric(0), cohortsize=numeric(0))

# Get list of central estimate files
files.v <- list.files(output.dir, include.dirs=FALSE, pattern="_") # Use pattern to exclude sub-dirs
files.v <- files.v[grep("^(?!PSA)", files.v, perl=TRUE)]           # Drop PSA files

# Extract vaccine program/subprogram (type) from filename
countries.v <- sapply(strsplit(files.v, "_"), "[", 1)
types.v <- paste(sapply(strsplit(files.v, "_"), "[", 2), 
                 sapply(strsplit(files.v, "_"), "[", 3), sep="_")

for (f in 1:length(files.v)){
  in.data <- read.csv(paste(output.dir, files.v[f], sep="/"), stringsAsFactors = FALSE)
  in.data$country <- rep(countries.v[f], times=dim(in.data)[1])
  in.data$type <- rep(types.v[f], times=dim(in.data)[1])
  
  toPlot <- bind_rows(toPlot, in.data[, names(toPlot)])
}

# Check for correct number of rows
# Should be 101 years * 71 age groups * x countries * 3 scenarios
length(toPlot$country)==(101*71*26*5)

# Add plotting variables
toPlot$AgeGroup <- cut(toPlot$AgeInYears, breaks=c(-1,4,9,14,19,24,29,49,70), 
                       labels=c('<5y','05-9y', '10-14y', '15-19y', '20-24y', '25-29y', '30-49y', '50-70y'))

toPlot$DeathRate <- 100000*toPlot$Deaths/toPlot$cohortsize
toPlot$IncRate <- 100000*toPlot$Cases/toPlot$cohortsize

### (3) Death rate by scenario ################################################

collapse <- toPlot %>%
  group_by(country, type, year) %>%
  summarize(cohortsize=sum(cohortsize), Deaths=sum(Deaths))
collapse$DthRate <- 100000*collapse$Deaths/collapse$cohortsize

p0 <- ggplot(data=collapse, aes(x=year, y=DthRate, color=type)) + theme_bw() +
  facet_grid(country~., scale="free_y") + theme(panel.grid=element_blank())

p0 + geom_line() +
  scale_color_manual(values=c("red", "firebrick4", "blue", "navyblue", "grey70"))

# Export as .pdf, 7"wide x 25"tall


### (4) Plot death rate by year and age #######################################

# The 'forcats::fct_rev' puts the youngest age groups at the bottom
g1 <- ggplot(data=toPlot, aes(x=year, y=DeathRate, fill=forcats::fct_rev(AgeGroup))) +
  theme_minimal() + facet_grid(country~type, scale="free_y") +
  guides(fill=guide_legend(title="Age group")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

g1 + geom_bar(stat='identity') +
  scale_fill_manual(values=fill.colors) +
  ylab("Death rate per 100,000")

# Aggregate across countries to match Montagu data check
agg.country <- toPlot %>%
  group_by(year, AgeInYears, type) %>%
  summarize(Deaths=sum(Deaths))

p1a <- ggplot(data=agg.country, aes(x=year, y=Deaths, fill=AgeInYears)) + theme_minimal() +
  theme(panel.grid=element_blank(), legend.position="top") + facet_grid(type~., scale="free_y")

p1a + geom_bar(stat='identity') +
  scale_fill_viridis(option="magma")


### (5) Plot Incidence by age and year ########################################

g2 <- ggplot(data=toPlot, aes(x=AgeInYears, y=IncRate, fill=year)) +
  theme_minimal() + facet_grid(country ~ type, scale="free_y") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

g2 + geom_bar(stat='identity') +
  scale_fill_viridis()


### (6) Vaccine impact ########################################################
# Since there are trivial differences between the default and bestcase vax    #
# coverage for the 2019 touchstone, just use default for this.                #

# Split by scenario
vax.none <- toPlot[toPlot$type=='none_default',][,c(1,3,4,5,6)]
vax.campaign <- toPlot[toPlot$type=='campaign_default',][,c(1,3,4,5,6)]
vax.both <- toPlot[toPlot$type=='both_default',][,c(1,3,4,5,6)]

names(vax.none)[4:5] <- c("Inc_none", "Dth_none")
names(vax.campaign)[4:5] <- c("Inc_camp", "Dth_camp")
names(vax.both)[4:5] <- c("Inc_both", "Dth_both")

# Merge to wide
impact <- vax.none %>%
  left_join(vax.campaign, by=c("country", "year", "AgeInYears")) %>%
  left_join(vax.both, by=c("country", "year", "AgeInYears"))


# Calculate impacts
impact$camp_dth <- impact$Dth_none - impact$Dth_camp
impact$both_dth <- impact$Dth_none - impact$Dth_both
impact$camp_inc <- impact$Inc_none - impact$Inc_camp
impact$both_inc <- impact$Inc_none - impact$Inc_both

impact$AgeGroup <- cut(impact$AgeInYears, breaks=c(-1,4,9,14,19,24,29,49,70), 
                       labels=c('<5y','05-9y', '10-14y', '15-19y', '20-24y', '25-29y', '30-49y', '50-70y'))


# Impact of campaigns, by year and age group
format(aggregate(cbind(camp_dth, both_dth, camp_inc, both_inc)~country, data=impact, FUN=sum),
       big.mark=",", scientific=F, digits=0)

g3 <- ggplot(data=impact, aes(x=year, y=camp_dth, fill=forcats::fct_rev(AgeGroup))) +
  theme_minimal() + facet_grid(~country) +
  guides(fill=guide_legend(title = "Age Group")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g3 + geom_bar(stat='identity') +
  scale_fill_manual(values=fill.colors)

g4 <- ggplot(data=impact, aes(x=year, y=both_dth, fill=forcats::fct_rev(AgeGroup))) +
  theme_minimal() + facet_grid(~country) +
  guides(fill=guide_legend(title = "Age Group")) +
  theme(panel.grid = element_blank())

g4 + geom_bar(stat='identity') +
  scale_fill_manual(values=fill.colors)





