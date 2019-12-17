



#### Program information ######################################################
# Source file name: Testing_Sim_Output.R                                      #
# Package: MenA_VaccSims                                                      #
# Version Date 11/08/2019                                                     #
#___ Program Description _____________________________________________________#
# Runs simulations and tests results for expected incidence/prevalence/pop    #
#_____________________________________________________________________________#
# Evaluations to make:
# (1) Overall incidence - frequency and size of outbreaks
# (2) Incidence by age during epidemic and endemic years
# (3) Colonization by age during epidemic and endemic years
# (4) Population size and age distribution
# (5) Vaccination - expected impact of campaign only, campaign + routine


library(ggplot2)

#inputdir<-"C:/Users/jackml4/Documents/Link_to_H_Drive/GAVI MenA predictions/Data/GAVI inputs/201810synthetic_downloaded_2019"
#outputdir <- "C:/Users/jackml4/Documents/Link_to_H_Drive/GAVI MenA predictions/Scratch/Temporary output directory"
#script.dir <- "C:/Users/jackml4/Documents/Link_to_H_Drive/GAVI MenA predictions/R_programming"

inputdir<-"G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Data/GAVI inputs/201810synthetic_downloaded_2019"
outputdir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Scratch/Temporary output directory"
script.dir <- "H:/Git/MenA R"



start <- as.Date("1950-01-01")

### (1) Evaluation functions ##################################################

yearlyIR <- function(input=finalpop){
  # Returns a vector of annual incidence per 100,000
  cases <- apply(X=input[, "Inc", ], MARGIN=2, FUN=sum)
  pop <- numeric(dim(input)[3])
  for (t in 1:length(pop)){
    pop[t] <- sum(input[ , , t])
  }
  
  years <- ceiling(length(cases))
  
  out.v <- numeric(floor(length(cases)/52.178))
  
  for (y in 1:length(out.v)){
    start.d <- floor((52.178)*(y-1)) + 1
    stop.d <- floor(y*52.178)
    
    out.v[y] <- 100000 * sum(cases[start.d:stop.d]) / mean(pop[start.d:stop.d])
  }
  
  return(out.v)
}


yearlyCases <- function(input=finalpop){
  # Returns an array of annual cases and population size by age group
  out.a <- array(data=NA, dim=c(dim(input)[1], 2, floor(dim(input)[3])/52.178))

  dimnames(out.a)[[2]] <- c("Cases", "Pop")
  
  for (y in 1:dim(out.a)[3]){
    start.d <- floor((52.178)*(y-1)) + 1
    stop.d <- floor(y*52.178)
    
    out.a[, "Cases", y] <- apply(X=input[,"Inc",start.d:stop.d], MARGIN=1, FUN=sum)
    out.a[, "Pop", y] <- apply(X=input[, , start.d+26], MARGIN=1, FUN=sum)
    
  }

  return(out.a)
  
}


carrPrev <- function(input=finalpop, start.dt=start){
  
  # Calculates the age-specific prevalence of colonization during the dry (week 13)
  # and the rainy season (week 39) each year, in defined age-groups
  week.count <- as.numeric(floor((end - start.dt)/7)) + 1
  
  week.a <- array(data=NA, dim=c(2, week.count))
  week.a[1,] <- week(seq(start.dt, end, 7))==13
  week.a[2,] <- week(seq(start.dt, end, 7))==39

  # Number colonized in each age group at each time point
  age.col <- array(data=NA, dim=c(1441, week.count))
  # Pop size at each time point
  age.pop <- array(data=NA, dim=c(1441, week.count))
  
  # Prevalence array
  prev.a <- array(data=NA, dim=c(2, 9, floor(week.count/52.178))) # (season, ages, year)
  dimnames(prev.a)[[1]] <- c("Dry", "Rainy")
  dimnames(prev.a)[[2]] <- c("A_Lt_1", "A_01_04", "A_05_09", "A_10_14", "A_15_19",
                             "A_20_24", "A_25_29", "A_Ge_30", "A_All")

  for (a in 1:1441){
    age.col[a,] <- apply(X=input[a, c("Nc", "Lc", "Hc"), 1:week.count], MARGIN=2, FUN=sum)
    age.pop[a,] <- apply(X=input[a, , 1:week.count], MARGIN=2, FUN=sum)
  }
  
  for (s in 1:2){
    prev.a[s, 1, ] <- apply(X=age.col[1:11, week.a[s,]], MARGIN=2, FUN=sum) / 
      apply(X=age.pop[1:11, week.a[s,]], MARGIN=2, FUN=sum)  
    prev.a[s, 2, ] <- apply(X=age.col[12:59, week.a[s,]], MARGIN=2, FUN=sum) / 
      apply(X=age.pop[12:59, week.a[s,]], MARGIN=2, FUN=sum)
    prev.a[s, 3, ] <- apply(X=age.col[60:119, week.a[s,]], MARGIN=2, FUN=sum) / 
      apply(X=age.pop[60:119, week.a[s,]], MARGIN=2, FUN=sum)
    prev.a[s, 4, ] <- apply(X=age.col[120:179, week.a[s,]], MARGIN=2, FUN=sum) / 
      apply(X=age.pop[120:179, week.a[s,]], MARGIN=2, FUN=sum)  
    prev.a[s, 5, ] <- apply(X=age.col[180:239, week.a[s,]], MARGIN=2, FUN=sum) / 
      apply(X=age.pop[180:239, week.a[s,]], MARGIN=2, FUN=sum)
    prev.a[s, 6, ] <- apply(X=age.col[240:299, week.a[s,]], MARGIN=2, FUN=sum) / 
      apply(X=age.pop[240:299, week.a[s,]], MARGIN=2, FUN=sum)  
    prev.a[s, 7, ] <- apply(X=age.col[300:359, week.a[s,]], MARGIN=2, FUN=sum) / 
      apply(X=age.pop[300:359, week.a[s,]], MARGIN=2, FUN=sum)  
    prev.a[s, 8, ] <- apply(X=age.col[360:1441, week.a[s,]], MARGIN=2, FUN=sum) / 
      apply(X=age.pop[360:1441, week.a[s,]], MARGIN=2, FUN=sum)
    prev.a[s, 9, ] <- apply(X=age.col[, week.a[s,]], MARGIN=2, FUN=sum) /
      apply(X=age.pop[, week.a[s,]], MARGIN=2, FUN=sum)
  }  
  
  return(prev.a)
}
  

### (2) Overall incidence #####################################################
# Check to make sure the simulation is creating epidemics of the correct size #
# and frequency.                                                              #

paramfixed.all <- GetModelParams(path=inputdir, region.val=myregion)

incid.l <- list()
casepop.l <- list()
prev.l <- list()





for (i in 1:117){
  print(i)
  paramfixed.row <- paramfixed.all[i,]
  
  finalpop<-MenASimulation(startdt=start, enddt=end, fp=paramfixed.row, initpop=initpop, vacc_program=vacc_program,
                           countryparams=myparams, region=myregion, country=mycountry, inputdir=inputdir)

  incid.l[[i]] <- yearlyIR()[51:151]
  casepop.l[[i]] <- yearlyCases()[, , 51:151]
  prev.l[[i]] <- carrPrev(input=finalpop, start.dt=as.Date("2000-01-01"))

  plot(1:101, incid.l[[i]], pch=16, xlab="Year", ylab="Incidence per 100,000")
  lines(1:101, incid.l[[i]])
  
  colors.v <- c("red", "orange", "yellow", "green", "blue", "violet", "purple")
  if (i==117){
    plot(1:101, prev.l[[i]][1,1,], pch=16, col="grey", ylim=c(0, max(prev.l[[i]])), 
         xlab="Year", ylab="Prev", main="Dry season prevalence")
    lines(1:101, prev.l[[i]][1,1,], col="grey")
    for (a in 1:7){
      points(1:101, prev.l[[i]][1, a+1, ], pch=16, col=colors.v[a])
      lines(1:101, prev.l[[i]][1, a+1, ], col=colors.v[a])
    }
  }
}
  
  
### (3) Incidence by age and epidemic/endemic type ############################
# Compare with observed data on incidence by age group and epidemic size.     #

sim.df <- data.frame(iter=numeric(0), agegrp=numeric(0), year=numeric(0),
                     incid=numeric(0), epi.type=character(0))

for (i in 1:117){
  casepop.i <- casepop.l[[i]]
  
  epi.type <- ifelse(incid.l[[i]] < 20, "a.minor",
                     ifelse(incid.l[[i]] < 100, "b.middle", "c.major"))

  add.df <- data.frame(iter=numeric(606), agegrp=numeric(606), year=numeric(606),
                       incid=numeric(606), epi.type=character(606))
  
  add.df$iter <- rep(i, times=606)
  add.df$agegrp <- rep(1:6, times=101)[order(rep(1:6, times=101))]
  add.df$year <- rep(1:101, times=6)
  add.df$epi.type <- rep(epi.type, times=6)
  
  add.df$incid[1:101] <- 100000 * apply(X=casepop.i[1:59, "Cases", ], MARGIN=2, FUN=sum) /
    apply(X=casepop.i[1:59, "Pop", ], MARGIN=2, FUN=sum)
  add.df$incid[102:202] <- 100000 * apply(X=casepop.i[60:119, "Cases", ], MARGIN=2, FUN=sum) /
    apply(X=casepop.i[60:119, "Pop", ], MARGIN=2, FUN=sum)
  add.df$incid[203:303] <- 100000 * apply(X=casepop.i[120:179, "Cases", ], MARGIN=2, FUN=sum) /
    apply(X=casepop.i[120:179, "Pop", ], MARGIN=2, FUN=sum)
  add.df$incid[304:404] <- 100000 * apply(X=casepop.i[180:239, "Cases", ], MARGIN=2, FUN=sum) /
    apply(X=casepop.i[180:239, "Pop", ], MARGIN=2, FUN=sum)
  add.df$incid[405:505] <- 100000 * apply(X=casepop.i[240:359, "Cases", ], MARGIN=2, FUN=sum) /
    apply(X=casepop.i[240:359, "Pop", ], MARGIN=2, FUN=sum)
  add.df$incid[506:606] <- 100000 * apply(X=casepop.i[1:359, "Cases", ], MARGIN=2, FUN=sum) /
    apply(X=casepop.i[1:359, "Pop", ], MARGIN=2, FUN=sum)
  
  sim.df <- rbind(sim.df, add.df)
  
}  

sim.mean.df <- aggregate(incid~agegrp+epi.type, data=sim.df, FUN=mean)
names(sim.mean.df)[3] <- "value"
sim.mean.df$type <- "Simulated"
  
obs.df <- data.frame(epi.type=rep(c("a.minor", "b.middle", "c.major"), times=6),
                      agegrp=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),
                      type="Observed")
obs.df$value <- c(11.7, 70, 152.3,
                      16, 96, 208.8,
                      14, 83, 180.8,
                      11.3, 68, 147.8,
                      4.3, 25.5, 55,
                      11.3, NA, 146.3)

sas.df <- obs.df
sas.df$type <- "SAS model"
sas.df$value <- c(10.0, 74.1, 194.5,
                  14.7, 109.0, 293.8,
                  13.6, 99.7, 269.7,
                  7.6, 56.5, 163.8,
                  5.2, 39.4, 118.1,
                  6.4, 48.2, 140.7)


plot.df <- rbind(obs.df, sim.mean.df, sas.df)
color1 <- rgb(26, 102, 123, maxColorValue=255)
color2 <- rgb(199, 30, 29, maxColorValue=255)

p1 <- ggplot(data=plot.df, aes(x=agegrp, y=value, fill=type)) + theme_bw() +
  theme(panel.grid=element_blank()) + facet_grid(rows=vars(epi.type), scales="free")

p1 + geom_col(position="dodge") +
  scale_fill_manual(values=c("grey30", color2, color1)) +
  scale_x_continuous(limits=c(0.5, 6.5), breaks=1:6, 
                     labels=c("<5", "5-9", "10-14", "15-19", "20-29", "All <30")) +
  xlab("Age group") + ylab("Incidence per 100,000")

rm(sim.df, agg.df, casepop.a, epi.type, case.agegrp, case.lt30)

### (3) Prevalence by age, season, and epidemic ###############################

sim.df <- data.frame(iter=numeric(0), agegrp=numeric(0), year=numeric(0),
                     value=numeric(0), epi.type=character(0))

for (i in 1:117){
  prev.i <- prev.l[[i]]
  
  epi.type <- ifelse(prev.i[1, "A_All", ] <= 0.025, "b.minor", "c.major")
  
  # Lazy coding - not efficient processing
  for (y in 1:dim(prev.i)[3]){
    for (a in 1:8){
      add.df <- data.frame(iter=i, agegrp=a, year=y, 
                           value = prev.i[2, a, y], epi.type="a.rainy")
      sim.df <- rbind(sim.df, add.df)
      add.df <- data.frame(iter=i, agegrp=a, year=y, 
                           value = prev.i[1, a, y], epi.type=epi.type[y])
      sim.df <- rbind(sim.df, add.df)
    }
    
  }
}

sim.mean.df <- aggregate(value~epi.type+agegrp, data=sim.df, FUN=mean)
sim.mean.df$type <- "c. Sim mean"

sim.sd.df <- aggregate(value~epi.type+agegrp, data=sim.df, FUN=stats::sd)
names(sim.sd.df)[3] <- "stdev"
sim.sd.df <- merge(sim.mean.df, sim.sd.df, by=c("agegrp", "epi.type"))
sim.sd.df$lower <- sim.sd.df$value - sim.sd.df$stdev
sim.sd.df$upper <- sim.sd.df$value + sim.sd.df$stdev


obs.df <- data.frame(epi.type=rep(c("a.rainy", "b.minor", "c.major"), times=8),
                     agegrp=rep(1:8, times=3)[order(rep(1:8, times=3))],
                     type="a. Observed")
obs.df$value <- c(0.0026, 0.0039, 0.0290,
                  0.0035, 0.0054, 0.0398,
                  0.0058, 0.0087, 0.0643,
                  0.0068, 0.0104, 0.0768,
                  0.0060, 0.0091, 0.0674,
                  0.0046, 0.0072, 0.0535,
                  0.0040, 0.0062, 0.0458,
                  0.0030, 0.0046, 0.0338)

sas.df <- obs.df
sas.df$type <- "b. SAS model"
sas.df$value <- c(0.002, 0.004, 0.033,
                  0.002, 0.004, 0.036,
                  0.003, 0.008, 0.066,
                  0.005, 0.011, 0.086,
                  0.004, 0.010, 0.074,
                  0.002, 0.005, 0.044,
                  0.002, 0.005, 0.042,
                  0.002, 0.005, 0.042)
  

plot.df <- rbind(obs.df, sas.df, sim.mean.df)
color1 <- rgb(26, 102, 123, maxColorValue=255)
color2 <- rgb(199, 30, 29, alpha=5, maxColorValue=255)
color3 <- rgb(199, 0, 50, maxColorValue = 255)

p1 <- ggplot(data=plot.df, aes(x=agegrp, y=value, fill=type, color=type)) + theme_bw() +
  theme(panel.grid=element_blank()) + facet_grid(rows=vars(epi.type), scales="free")

p1 + geom_point(size=3) +
  geom_errorbar(data=sim.sd.df, aes(x=agegrp, ymin=lower, ymax=upper), color=color3, width=0) +
  scale_color_manual(values=c("grey30", color1, color3))

