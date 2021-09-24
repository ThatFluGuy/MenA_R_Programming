#Provide solid documentation about datasets and sources, and non-a calculations (CWXY sum)

# To do: Percents in age distribution only add up to 0.998.
# To do: Account for undersurveillance? Already done in the .csv file?

library(dplyr)
glm.nb <- MASS::glm.nb
rnegbin <- MASS::rnegbin

cwyx.path <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Data/ACWXY"
pop.path <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Data/GAVI inputs/201910gavi_v4"

### (1) Import the raw data and fit the model #################################
# Exclude countries with <8 annual datapoints from the model.                 #
# Drop ETH, GMB, GIN, MRT, SSD, SDN, UGA on this basis.                       #

nona <- read.csv(paste0(cwyx.path, "/WHO Bulletin data.csv"),
                 stringsAsFactors = FALSE)

nona2 <- nona %>% filter(!(Country %in% c("ETH", "GMB", "GIN", "MRT", "SSD", "SDN", "UGA")))

nona2$Model_Country <- factor(nona2$Country)

nb.nona <- glm.nb(Num_nonA ~ Model_Country + offset(log(Pop)), data=nona2)

#ss <- subset(nona, Country=="Niger")

#plot(ss$Year, ss$Num_nonA)

#nona.2020 <- subset(nona, year==2020)

#testing <- predict(nb.nona, newdata=nona.2020)

#sims <- rnegbin(exp(testing), theta=nb.nona$theta)

#sims.df <- data.frame(Country=nona.2020$Country, total=sims, AgeLT5 = sims*0.123, Age5t9 = sims*0.149, Age10t14 = sims*0.185, Age15t19=sims*0.145, Age20t29=sims*0.163, Age30t39=sims*0.110, Age40t49=sims*0.061, Age50p=sims*0.062)

### (2) Get country- and year-specific predicted case counts ##################
# Import the population demographic data provided by Gavi for this.           #
# Link countries to the appropriate modeled country based on being hyper or   #
# non-hyperendemic.                                                           #
# The following countries also do not have any observations in the data:      #
# BDI, ERI, GNB, KEN, RWA, TZA                                                #
# Hyper-endemic countries are: BFA, ETH, MLI, NER, NGA, SDN, TCD.             #
# Assume that Benin is representative of the non-hyperendemic countries, and  #
# that Burkina is representative of hyperendemic countries.                   #

pop.df <- read.csv(paste0(pop.path, "/201910gavi-4_dds-201910_2_tot_pop_both.csv"),
                   stringsAsFactors = FALSE)

pop.df2 <- pop.df %>% select(country_code, year, value) %>% 
  rename(Country=country_code, Pop=value) %>%
  filter(year >= 2000,
         Country %in% c("BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "ERI", 
                        "ETH", "GHA", "GIN", "GMB", "GNB", "KEN", "MLI", "MRT", 
                        "NER", "NGA", "RWA", "SDN", "SEN", "SSD", "TCD", "TGO",
                        "TZA", "UGA"))

c.nonh <- c("BDI", "ERI", "GIN", "GMB", "GNB", "KEN", "MRT", "RWA", "SSD", "TZA", "UGA")
c.h <- c("ETH", "SDN")
         
pop.df2$Model_Country <- ifelse(pop.df2$Country %in% c.nonh, "BEN",
                                ifelse(pop.df2$Country %in% c.h, "BFA", pop.df2$Country))

pop.df2$Model_Country <- factor(pop.df2$Model_Country)

# Create a data.frame with 200 replicates of the population data
reptimes <- 200 
idx <- rep(1:nrow(pop.df2), reptimes) 
pred.df <- pop.df2[idx, ] %>% arrange(Country, year)
pred.df$index <- rep(1:200, times=max(idx))

# Get the predicted outcome for each country/year combination
pred.df$predict <- predict(nb.nona, newdata=pred.df)

# Random draw from negative binomial distribution for each country/year
pred.df$cases <- rnegbin(exp(pred.df$predict), theta=nb.nona$theta)

# Divide into age groups
pred.df$AgeLT5 <- pred.df$cases * 0.123
pred.df$Age5t9 <- pred.df$cases * 0.149
pred.df$Age10t14 <- pred.df$cases * 0.185
pred.df$Age15t19 <- pred.df$cases * 0.145
pred.df$Age20t29 <- pred.df$cases * 0.163
pred.df$Age30t39 <- pred.df$cases * 0.110
pred.df$Age40t49 <- pred.df$cases * 0.061
pred.df$Age50p <- pred.df$cases * 0.062

output.df <- pred.df %>%
  select(index, Country, year, cases, AgeLT5, Age5t9, Age10t14, Age15t19, Age20t29,
         Age30t39, Age40t49, Age50p)

### (3) Export ################################################################
# Export a separate .csv file for each modeled country.                       #

for (c in unique(output.df$Country)){
  temp.df <- output.df %>% filter(Country==c)
  write.csv(temp.df, paste0(cwyx.path, "/MenCWYX_Sims_", c, ".csv"), row.names=FALSE)
}







#Variances are larger than means on average, so zero-inflated Poisson isn't applicable
#round(cbind(tapply(nona$Num_nonA, nona$Country, mean, na.rm=TRUE), tapply(nona$Num_nonA, nona$Country, var, na.rm=TRUE)), 2)

