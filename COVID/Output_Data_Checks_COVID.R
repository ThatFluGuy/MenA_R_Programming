### Program information #######################################################
# Program name: Output_Data_Checks.R                                          #
#_____________________________________________________________________________#
# Program location:                                                           #
# G:\CTRHS\Modeling_Infections\GAVI MenA predictions\Programming\Impact estimation
#_____________________________________________________________________________#
#  Input datasets: central_ccccc_tttt.sas7bdat, where ccccc is a country and  #
#  tttt is a vaccination scenario                                             #
#_____________________________________________________________________________#
# Program description:                                                        #
# This program runs some basic checks on simulation output from the GAVI      #
# forecast outputs. Steps in this program:                                    #
# 1) Set up directories                                                       #
# 2) Import data                                                              #
# 3) Plot deaths by year                                                      #
# 4) Plot deaths by age and year                                              #
# 5) Plot incidence by year and age                                           #
# 6) Vaccine impact                                                           #
# 7) PSA variability                                                          #
# 8) Compare submission .csv with .sas7bdat                                   #
#_____________________________________________________________________________#
# Version 1.0. Created 8th November 2017, by Mike Jackson.                    #
#_____________________________________________________________________________#



### (1) Import packages, set up directories and datasets ######################
library(tibble)
library(haven) # Used to convert a SAS .sas7bdat file to an R tibble
library(ggplot2)
library(RColorBrewer)
library(viridis)

fill.colors <- brewer.pal(11, name='Spectral')[c(1:4,8:11)]

#sas.dir <- "C:/Users/jackml4/Documents/Link_to_H_Drive/GAVI MenA predictions/Analysis/Simulation results/"
output.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Analysis/Simulation results/"

country.list <- c("bfa", "ner", "nga", "tcd")

### ANALYSES FOR CENTRAL ESTIMATES ############################################
### (2) Import simulation data ################################################

# Create empty tibble as starting point
toPlot <- tibble(year=numeric(0), AgeInYears=numeric(0), Cases=numeric(0), Deaths=numeric(0),
                 DALYs=numeric(0), simulations=numeric(0), cohortsize=numeric(0),
                 country=character(0), scenario=character(0))

files.v <- list.files(path=output.dir, pattern="scenario")
psa.v <- grep("PSA", files.v)
files.v <- files.v[-psa.v]

for (f in 1:length(files.v)){
  in.data <- read.csv(paste0(output.dir, files.v[f]), stringsAsFactors = FALSE)
  in.data$country <- substr(files.v[f], 1, 3)
  in.data$scenario <- strsplit(files.v[f], "_")[[1]][2]
  toPlot <- rbind(toPlot, in.data[, !(names(in.data) %in% "X")])
}

# Check for correct number of rows
# Should be 101 years * 71 age groups * x countries * 9 scenarios
length(toPlot$country)==(101*71*4*9)

# Add plotting variables
toPlot$AgeGroup <- cut(toPlot$AgeInYears, breaks=c(-1,4,9,14,19,24,29,49,70), 
                       labels=c('<5y','05-9y', '10-14y', '15-19y', '20-24y', '25-29y', '30-49y', '50-70y'))

toPlot$DeathRate <- 100000*toPlot$Deaths/toPlot$cohortsize
toPlot$IncRate <- 100000*toPlot$Cases/toPlot$cohortsize

### (3) Death rate by scenario ################################################

collapse <- aggregate(cbind(cohortsize, Deaths)~country+scenario+year, data=toPlot, FUN=sum)
collapse$DthRate <- 100000*collapse$Deaths/collapse$cohortsize

p0 <- ggplot(data=collapse, aes(x=year, y=DthRate, color=scenario)) + theme_bw() +
  facet_grid(country~., scale="free_y") + theme(panel.grid=element_blank())

p0 + geom_line()

# Export as .pdf, 7"wide x 25"tall


### (4) Plot death rate by year and age #######################################

# The 'forcats::fct_rev' puts the youngest age groups at the bottom
g1 <- ggplot(data=toPlot, aes(x=year, y=DeathRate, fill=forcats::fct_rev(AgeGroup))) +
  theme_minimal() + facet_grid(country~scenario, scale="free_y") +
  guides(fill=guide_legend(title="Age group")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

g1 + geom_bar(stat='identity') +
  scale_fill_manual(values=fill.colors) +
  ylab("Death rate per 100,000")

# Aggregate across countries to match Montagu data check
agg.country <- aggregate(Deaths~year+AgeInYears+scenario, data=toPlot, FUN=sum)

p1a <- ggplot(data=agg.country, aes(x=year, y=Deaths, fill=AgeInYears)) + theme_minimal() +
  theme(panel.grid=element_blank(), legend.position="top") + 
  facet_grid(scenario~., scale="free_y")

p1a + geom_bar(stat='identity') +
  scale_fill_viridis(option="magma")


### (5) Plot Incidence by age and year ########################################

g2 <- ggplot(data=toPlot, aes(x=Age, y=IncRate, fill=Year)) +
  theme_minimal() + facet_grid(country ~ type, scale="free_y") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

g2 + geom_bar(stat='identity') +
  scale_fill_viridis()


### (6) Vaccine impact ########################################################

# Split by scenario
vax.none <- toPlot[toPlot$type=='none',][,c(1,3,4,6:7)]
vax.campaign <- toPlot[toPlot$type=='campaign',][,c(1,3,4,6:7)]
vax.both <- toPlot[toPlot$type=='both',][,c(1,3,4,6:7)]

names(vax.none)[4:5] <- c("Dth_none", "Inc_none")
names(vax.campaign)[4:5] <- c("Dth_camp", "Inc_camp")
names(vax.both)[4:5] <- c("Dth_both", "Inc_both")

# Merge to wide
impact <- merge(vax.none, vax.campaign, by=c("country","Year","Age"))
impact <- as.tibble(merge(impact, vax.both, by=c("country","Year","Age")))

# Calculate impacts
impact$camp_dth <- impact$Dth_none - impact$Dth_camp
impact$both_dth <- impact$Dth_none - impact$Dth_both
impact$camp_inc <- impact$Inc_none - impact$Inc_camp
impact$both_inc <- impact$Inc_none - impact$Inc_both

impact$AgeGroup <- cut(impact$Age, breaks=c(-1,4,9,14,19,24,29,49,70), 
                       labels=c('<5y','05-9y', '10-14y', '15-19y', '20-24y', '25-29y', '30-49y', '50-70y'))


# Total impact
format(aggregate(cbind(camp_dth, both_dth, camp_inc, both_inc)~country, data=impact, FUN=sum),
       big.mark=",", scientific=F, digits=0)

g3 <- ggplot(data=impact, aes(x=Year, y=camp_dth, fill=forcats::fct_rev(AgeGroup))) +
  theme_minimal() + facet_grid(~country) +
  guides(fill=guide_legend(title = "Age Group")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g3 + geom_bar(stat='identity') +
  scale_fill_manual(values=fill.colors)

### (7) Compare .csv with .sas7bdat ###########################################

csv.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Deliverables/Deliverables 2017/"

none.csv <- read.csv(file=paste(csv.dir, "central_none.csv", sep=""), stringsAsFactors = F)
camp.csv <- read.csv(file=paste(csv.dir, "central_campaign.csv", sep=""), stringsAsFactors = F)
both.csv <- read.csv(file=paste(csv.dir, "central_both.csv", sep=""), stringsAsFactors = F)

names(none.csv) <- c("disease", "Year", "Age", "country", "country_name",
                     "cohort_size_csv", "Deaths_csv", "Cases_csv",
                     "DALYs_csv")
names(camp.csv) <- c("disease", "Year", "Age", "country", "country_name",
                     "cohort_size_csv", "Deaths_csv", "Cases_csv",
                     "DALYs_csv")
names(both.csv) <- c("disease", "Year", "Age", "country", "country_name",
                     "cohort_size_csv", "Deaths_csv", "Cases_csv",
                     "DALYs_csv")

none.check <- merge(toPlot[toPlot$type=="none", c(1,3:8)], none.csv, 
                    by=c("country", "Year", "Age"))
camp.check <- merge(toPlot[toPlot$type=="campaign", c(1,3:8)], camp.csv, 
                    by=c("country", "Year", "Age"))
both.check <- merge(toPlot[toPlot$type=="both", c(1,3:8)], both.csv, 
                    by=c("country", "Year", "Age"))

none.check$Deaths_diff <- none.check$Deaths - none.check$Deaths_csv
none.check$Cases_diff <- none.check$Cases - none.check$Cases_csv
none.check$DALYs_diff <- none.check$DALYs - none.check$DALYs_csv
none.check$cohort_diff <- none.check$cohort_size - none.check$cohort_size_csv
range(none.check$Deaths_diff)
range(none.check$Cases_diff)
range(none.check$DALYs_diff)
range(none.check$cohort_diff) # Slightly different; due to limits of variable size?


camp.check$Deaths_diff <- camp.check$Deaths - camp.check$Deaths_csv
camp.check$Cases_diff <- camp.check$Cases - camp.check$Cases_csv
camp.check$DALYs_diff <- camp.check$DALYs - camp.check$DALYs_csv
camp.check$cohort_diff <- camp.check$cohort_size - camp.check$cohort_size_csv
range(camp.check$Deaths_diff)
range(camp.check$Cases_diff)
range(camp.check$DALYs_diff)
range(camp.check$cohort_diff) # Slightly different; due to limits of variable size?


both.check$Deaths_diff <- both.check$Deaths - both.check$Deaths_csv
both.check$Cases_diff <- both.check$Cases - both.check$Cases_csv
both.check$DALYs_diff <- both.check$DALYs - both.check$DALYs_csv
both.check$cohort_diff <- both.check$cohort_size - both.check$cohort_size_csv
range(both.check$Deaths_diff)
range(both.check$Cases_diff)
range(both.check$DALYs_diff)
range(both.check$cohort_diff) # Slightly different; due to limits of variable size?


### (8) ANALYSES FOR PSA DATASETS #############################################
# Items here:                                                                 #
# a) Compare PSA cohort to central cohort                                     #
# b) Compare PSA mean cases/deaths/DALYs to central estimates                 #
# c) Display PSA variability                                                  #

for (country in country.list){
  for (type in type.list){
    print(paste(country, type, sep=" "))
    psa.df <- read_sas(paste(sas.dir, "psa_", country, "_", type, ".sas7bdat", sep=""))
    
    # a) Compare pop size from one PSA run to central
    comp.df1 <- toPlot[toPlot$country==toupper(country) & toPlot$type==type, 
                       names(toPlot) %in% c("Year", "Age", "cohort_size")]
    names(comp.df1)[3] <- "cohort_central"
    
    comp.df2 <- psa.df[psa.df$Loop==1, names(psa.df) %in% c("Year", "Age", "cohort_size")]
    names(comp.df2)[3] <- "cohort_psa"
    
    comp.df <- merge(comp.df1, comp.df2, by=c("Year", "Age"))
    print("Range of demographic differences:")
    print(range(comp.df$cohort_central - comp.df$cohort_psa))
    
    # b) Compare PSA means to central estimates
    psa.sum.df <- aggregate(cbind(Cases, Deaths, DALYs)~Year, data=psa.df, FUN=mean)
    central.sum.df <- aggregate(cbind(Cases, Deaths, DALYs)~Year,
                                data=toPlot[toPlot$country==toupper(country)&toPlot$type==type,],
                                FUN=mean)
    psa.sum.df <- psa.sum.df[order(psa.sum.df$Year),]
    central.sum.df <- central.sum.df[order(central.sum.df$Year),]
    
    ylims <- c(max(psa.sum.df$Cases, central.sum.df$Cases),
               max(psa.sum.df$Deaths, central.sum.df$Deaths),
               max(psa.sum.df$DALYs, central.sum.df$DALYs))
    
    plot(x=psa.sum.df$Year, psa.sum.df$Cases, main=paste(country, type, "Cases", sep=" "),
         type='l', col="red", ylim=c(0, ylims[1]))
    lines(x=central.sum.df$Year, central.sum.df$Cases, col="blue")

    plot(x=psa.sum.df$Year, psa.sum.df$Deaths, main=paste(country, type, "Deaths", sep=" "),
         type='l', col="red", ylim=c(0, ylims[2]))
    lines(x=central.sum.df$Year, central.sum.df$Deaths, col="blue")
    
    plot(x=psa.sum.df$Year, psa.sum.df$DALYs, main=paste(country, type, "DALYs", sep=" "),
         type='l', col="red", ylim=c(0, ylims[3]))
    lines(x=central.sum.df$Year, central.sum.df$DALYs, col="blue")
    
    # c) PSA variability
    psa.sd.df <- aggregate(Cases~Year, data=psa.df, FUN=sd)
    names(psa.sd.df)[2] <- "Cases_sd"
    
    psaPlot <- merge(psa.sum.df, psa.sd.df, by="Year")
    psaPlot$lcl <- psaPlot$Cases - psaPlot$Cases_sd
    psaPlot$ucl <- psaPlot$Cases + psaPlot$Cases_sd
    
    g4 <- ggplot(data=psaPlot, aes(x=Year, y=Cases)) + theme_bw() +
      theme(panel.grid = element_blank())
    
    g4 + geom_point() + 
      ggtitle(paste(country, type, "Cases and std dev", sep=" ")) +
      geom_errorbar(data=psaPlot, aes(x=Year, ymin=lcl, ymax=ucl))
    
  }
  # Pause execution to view graphs
  invisible(readline(prompt="Press [enter] to continue"))
  
}




