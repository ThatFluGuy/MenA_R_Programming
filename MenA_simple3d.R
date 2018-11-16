begin<-Sys.time()
#first performance benchmarking: 1/26/18 simuation 70 sec
#after vectorization 2/8/18 its 3 sec!!!
#SAS is < 10 sec on server, a little over 30 locally
#after summarization re-write, 1 sim with summarization = 8 sec, 10 sims=12 sec, 100 sims=56 sec
library(lubridate)
library(dplyr)
library(data.table) 
#directory containing inputs
setwd("\\\\HOME/stewcc1/MenAModel/Rdata")
mycountry <- "ETH"
start <- as.Date("2001-01-01")
end <- as.Date("2100-12-31")
myregion <- "not_hyper"
PSA <- FALSE
Vaccination<-FALSE
program <- "none" ## "campaign" or "routine" or "both" or "none"
phi<-0.2
nSims<-10

#country-specific parameters
params<-read.csv("country_params.csv")
myparams <-params[params$country==mycountry & params$year>=year(start) & params$year<=year(end) ,]

if (Vaccination!=FALSE) {
  vaccdf <- read.csv("VaxCover_ETH.csv")
  myvacc <- vaccdf[vaccdf$country=='ETH',]
  #make as vector of years where nothing happens (empty except for campaign only) for efficiency
  if (program=="campaign") {
    nodoses<-as.vector(myvacc[is.na(myvacc$DosesCampaign) | myvacc$DosesCampaign==0,2])
  }
}

#fixed parameters
rc <- 0.23334 # Recovery from carriage, scaled to weeks (12.175/52.177) (30 days);
rd <- 0.70002 # Recovery from disease, scaled to weeks (36.5249/52.177) (10 days);
ld = 0.90 # Reduction in risk of invasive disease associated with low antibody;
hd = 1.00 # Reduction in risk of invasive disease associated with high antibody;
lc = 0.25 # Reduction in risk of carriage associated with low antibody;
hc = 0.75 #  Reduction in risk of carriage associated with high antibody;
foii = 0.00000096 # Constant force of infection from immigration (0.00005/52.177);
#disease model for rainy and dry
dd <-c(0.0019, -0.0000002)
dr <-c(0.0018, -0.00000021)

#Make WAIFW matrices, expanded:
expandWaifw<-function(waifw){
  # repeat what was originally columns :
  #b[x,1] 60x; b[x,2] 96x; b[1,3] 84x; b[1,4] 120x
  #needs to go to 361 - add extra line at end for last big bucket
  rbind ( 
    matrix(data=waifw[c(1,5,9,13)], nrow=60, ncol=4, byrow=TRUE),
    matrix(data=waifw[c(2,6,10,14)], nrow=96, ncol=4, byrow=TRUE),
    matrix(data=waifw[c(3,7,11,15)], nrow=84, ncol=4, byrow=TRUE),
    matrix(data=waifw[c(4,8,12,16)], nrow=121, ncol=4, byrow=TRUE)
  )
}

waifwin<-read.csv("WAIFW_both.csv", stringsAsFactors = FALSE)  #vector
Rwaifw<-waifwin[waifwin$region==myregion & waifwin$season=='rainy', 4]
br<-expandWaifw(waifw=Rwaifw)
Dwaifw<-waifwin[waifwin$region==myregion & waifwin$season=='dry', 4]
bd<-expandWaifw(waifw=Dwaifw)

#Parameters for and initializing population; the 3D matrix is the population array with depth time
startSize <- params[params$country==mycountry & params$year==year(start)-1, 5]
inpop<-read.csv("PopInputAgeDist.csv")
mypop<-inpop[inpop$Country=='ETH',]
#age-specific death rates (for PSA=no, other option not implemented yet)
cfr <- c(0.106, 0.096, 0.089, 0.086, 0.079, 0.122) #used AFTER simulation macro in SAS
dur<- ceiling(difftime(end, start, units = "weeks"))
#add a dxstate to array for incident cases
pop <- array(data=0, dim=c(361, 9, dur+1)) # Dimensions are [age groups, states, time]
dimnames(pop)[[2]] <- c("Ns","Nc","Ls","Lc","Hs","Hc","Di","Va", "Inc")
##get age-specific proportions of each disease state into vectors
dist<-read.csv("dist_both.csv", stringsAsFactors = TRUE)
distcol<-ifelse(myregion=='hyper', 4, 3)
#expand age group fraction as vector to match pop matrix dimension 1
agefract <- c(rep(mypop[,c(2:8)], each=60))
statefract<-as.vector(dist[,distcol]) # fraction of each disease state in each of 7 population groups
statemx<-rbind(
  matrix(rep(statefract[1:7], each=60), nrow=60),
  matrix(rep(statefract[8:14], each=60), nrow=60),
  matrix(rep(statefract[15:21], each=60), nrow=60),
  matrix(rep(statefract[22:28], each=60), nrow=60),
  matrix(rep(statefract[29:35], each=60), nrow=60),
  matrix(rep(statefract[36:42], each=60), nrow=60),
  matrix(rep(statefract[43:49], each=60), nrow=60)
)
#initialize - 3rd dimension stays at 1
chunks <- 60
for (i in 1:361) {
  if (i==361) {chunks<-1}
  pop[i, "Ns", 1]<-(startSize*agefract[[i]]/chunks)*statemx[i,1]
  pop[i, "Nc", 1]<-(startSize*agefract[[i]]/chunks)*statemx[i,2]
  pop[i, "Ls", 1]<-(startSize*agefract[[i]]/chunks)*statemx[i,3]
  pop[i, "Lc", 1]<-(startSize*agefract[[i]]/chunks)*statemx[i,4]
  pop[i, "Hs", 1]<-(startSize*agefract[[i]]/chunks)*statemx[i,5]
  pop[i, "Hc", 1]<-(startSize*agefract[[i]]/chunks)*statemx[i,6]
  pop[i, "Di", 1]<-(startSize*agefract[[i]]/chunks)*statemx[i,7]
  #leave Vacc, Inc 0
}
#startpop<-t(pop[,,1]) 
GetInfectiveRatio<-function(inpop){
  #returns a four-element vector containing the ratio of infectives in 4 age groups (0-<5, 5-<13, 14-<20, 20+)
  #sum total pop 
  mosums<-colSums(inpop)
  N1<-(sum(mosums[1:60]))
  N2<-(sum(mosums[61:156]))
  N3<-(sum(mosums[157:240]))
  N4<-(sum(mosums[241:361]))
  #sum infectious pop (hard code rows, be careful if changing the matrix)
  infpop<-inpop[c(2,4,6,7),]
  infmosums<-colSums(infpop)
  I1<-(sum(infmosums[1:60]))
  I2<-(sum(infmosums[61:156]))
  I3<-(sum(infmosums[157:240]))
  I4<-(sum(infmosums[241:361]))
  infRatio<-c(I1/N1, I2/N2, I3/N3, I4/N4)
  return(infRatio)
}
vaccinate<-function(popslice, vlookup, type, mydate) {
  eligibles<-c("Ns","Nc","Ls","Lc","Hs","Hc")
  if ((type=="campaign" | type=="both") & (month(mydate)==10)) {
    #get parameters
    cDoses <- vlookup[vlookup$year==year(mydate),3]
    #zero-length cDoses (not NA, apparently) is blowing things up
    if (length(cDoses)> 0) {
      if (!is.na(cDoses)) {
        #change ages i=12 to 359, here 13 to 360
        eligN <- sum(popslice[13:360, eligibles]) 
        pcNLH <- ifelse(cDoses<=eligN, cDoses/eligN, 1 )
        #why do the above if were not going to move anybody?
        popslice[13:360,"Va"] <- popslice[13:360,"Va"] + pcNLH * rowSums(popslice[13:360,eligibles])
        popslice[13:360,eligibles] <- (1-pcNLH) * popslice[13:360, eligibles]
      }
    }
  }
  if (type=="routine" | type=="both") {
    pr <- vlookup[vlookup$year==year(mydate),4]
    if (length(pr)> 0) {
      if (!is.na(pr)) {
        #routine vaccinations, when subjects turn 9 months- thats 10 here
        popslice[10,"Va"]<- popslice[10,"Va"] + (pr*sum(popslice[10,eligibles]))
        popslice[10,eligibles]<- (1-pr)*popslice[10,eligibles]
      }
    }
  }
  return(popslice)
}

theDate <- start
births <- myparams[myparams$year==year(theDate), 6]/52.1775
imr <- myparams[myparams$year==year(theDate), 10]/(1000*52.1775)
v <- myparams[myparams$year==year(theDate), 11] / (1000*52.1775)
deathvec<-c(rep(imr,12), rep(v, 349))
#waning age-group vectors 0-5mo, 6mo-2y, 3-10y, 11+y
wanev <- c(rep(1-imr, 7), rep(0.000172, 17), rep(0.000096, 107), rep(0.000307, 230))  #waning from vacc to hi ab scaled to wks-confirm wv(1) = NA
waneh <- c(rep(0.01092, 6), rep(0.00654, 18), rep(0.00527, 107), rep(0.00096, 230)) #waning from high to low ab scaled to weeks
wanel <- c(rep(0.00970,6), rep(0.00487, 18), rep(0.00364, 107), rep(0.00057, 230)) #waning from low 
j <- 2 #initialize time interval counter - the first position is filled with starting pop
iterDate<-start-7 #this will make aging happen the first iteration, like SAS - also supply date for initial pop
#dimnames(pop)[[3]][[1]] <- as.character(iterDate)
emptypop<-pop
LastMonth <- month(iterDate)
iterYear<-year(start-7)
random<-runif(1)
sto <- 1 + (phi * cos(random*pi))
#sto<-1 #print(paste0("new stochastic parameter ", sto))
my_data <- list()
for (n in 1:nSims) {
  while (theDate <= end) {
    #store the date in a vector to label iterations later
    #try naming them now, then maybe the following won't be necessary / doesn't seem to work
    #dimnames(pop)[[3]][[j]] <- as.character(theDate)
    iterDate<-c(iterDate, theDate)
    iterYear<-c(iterYear, year(theDate))
    if (day(theDate) <= 7 & month(theDate) == 9)  #yearly update to stochastic param
    {
      random<-runif(1)
      sto <- 1 + (phi * cos(random*pi))
      #sto<-1 #print(paste0("new stochastic parameter ", sto))
    }
    #yearly birth, death and imr rates, only needs to when year changes
    if (month(theDate)== 1 & LastMonth==12)
    {
      births <-myparams[myparams$year==year(theDate), 6]/52.1775
      imr <- myparams[myparams$year==year(theDate), 10]/(1000*52.1775)
      v <- myparams[myparams$year==year(theDate), 11] / (1000*52.1775)
      #new imr, need to update wanev and death
      wanev <- c(rep(1-imr, 7), wanev[8:361]) 
      deathvec<-c(rep(imr,12), rep(v, 349))
    }
    
    #  waifw matrix depends on rainy (Mar-Aug) or dry (Sep-Feb) season
    if (month(theDate)!=LastMonth)
    {
      if (month(theDate) %in% c(3,4,5,6,7,8)) { wmx <- br
      } else { wmx <-  bd
      }
      # Determine age-specific per-carrier risk of invasive disease;
      #like infection, this is also seasonal but peaks later, hence different months (rainy here = 6-10 June-Oct)
      iagevec<-seq(0, 360)
      if (month(theDate) %in% c(1,2,3,4,5,11,12)) {sigma = sigmavec<-dd[1] + dd[2] *iagevec
      } else {sigmavec<-dr[1] + dr[2] *iagevec
      }
    } #end of updates conditional on month
    
    #infectious ratios by pop age group (for calculating force) this happens every time point
    IR<-GetInfectiveRatio(t(pop[,,j-1]))
    # Includes force of infection from outside the population (foii, was for in SAS code)
    forcevec <- (sto*wmx[,1]*IR[1]) + (sto*wmx[,2]*IR[2]) + (sto*wmx[,3]*IR[3])+ (sto*wmx[,4]*IR[4]) + foii
    
    #vectorized (1-361) calculations instead of loop
    pop[,"Ns",j] <- wanel*pop[,"Ls",j-1] - (deathvec+forcevec)*pop[,"Ns",j-1] + pop[,"Ns",j-1]
    #births - this order is the same as sas (could just add births to line above)
    pop[1,"Ns",j] <- pop[1,"Ns",j] + births
    pop[,"Nc",j] <- forcevec*pop[,"Ns",j-1] - (deathvec+rc+sigmavec) * pop[,"Nc",j-1] + pop[,"Nc",j-1]
    pop[,"Ls",j] <- waneh*pop[,"Hs",j-1] + rc*pop[,"Nc",j-1] - (deathvec+wanel+(1-lc)*forcevec) * pop[,"Ls",j-1] + pop[,"Ls",j-1]
    pop[,"Lc",j] <- (1-lc)*forcevec*pop[,"Ls",j-1] - (deathvec+rc+(1-ld)*sigmavec) * pop[,"Lc",j-1]  + pop[,"Lc",j-1]
    pop[,"Hs",j] <- rc*(pop[,"Lc",j-1] + pop[,"Hc",j-1]) + rd*pop[,"Di",j-1]  + wanev*pop[,"Va",j-1] - (deathvec+waneh+(1-hc)*forcevec) * pop[,"Hs",j-1] + pop[,"Hs",j-1]
    pop[,"Hc",j] <- (1-hc)*forcevec*pop[,"Hs",j-1] - (deathvec+rc-(1-hd)*sigmavec)*pop[,"Hc",j-1] + pop[,"Hc",j-1]
    pop[,"Di",j] <- sigmavec*(pop[,"Nc",j-1] + (1-ld)*pop[,"Lc",j-1] + (1-hd)*pop[,"Hc",j-1]) - (deathvec+rd)*pop[,"Di",j-1] + pop[,"Di",j-1]
    pop[,"Va",j] <- -(deathvec+wanev)* pop[,"Va",j-1] + pop[,"Va",j-1]
    # Count incident cases of invasive disease NOTE THESE ARE before incrementing = use old pop
    pop[,"Inc",j] <- sigmavec*(pop[,"Nc",j-1] + (1-ld)*pop[,"Lc",j-1] + (1-hd)*pop[,"Hc",j-1])
    #Inc[,j] <- sigmavec*(pop[,"Nc",j-1] + (1-ld)*pop[,"Lc",j-1] + (1-hd)*pop[,"Hc",j-1])
    
    #aging handled monthly ** THIS WORKS **aging incident cases too
    if (month(theDate) != LastMonth){
      for (x in 1:9) {
        #if (x==1) { beforeaging<-pop[,x,j]}
        #save old 361 spot to add to last pot later
        save361<-pop[[361,x,j]]
        #vector to shift: pop[m,x,j]
        pop[,x,j]<-(c(0,  pop[,x,j])[1 : length(pop[,x,j])])
        #move last monthly group into 30+ pot
        pop[[361,x,j]]<-pop[[361,x,j]] + save361
      } #end of x loop
      #Vaccinate here, monthly
      if (Vaccination!=FALSE) {
        #there are some cases where theres nothing to do in many years
        if  (program=="campaign") { 
          if (year(theDate) %in% nodoses==FALSE) {
            pop[,,j] <- vaccinate(pop[,,j], myvacc, program, theDate)
            #print("vaccinating")
          }
        }
        else {
          pop[,,j] <- vaccinate(pop[,,j], myvacc, program, theDate)
        }
      }
    } #end of aging IF
    
    #}
    # if (year(theDate) > 2020) {break}
    LastMonth <- month(theDate)
    theDate<-theDate+7
    j <- j + 1
  } #end of date loop
  head(pop[,"Inc", 200:225])
  
  #str(pop)  #dims age*group*time
  # 1/22/18 add summarization code
  #need to be able to shut it off if breaking above or it will blow up
  summarizeme<-1
  if (summarizeme > 0) {
    #label matrix iteration dimension with dates
    dimnames(pop)[[3]] <- as.Date(iterDate)
    inclong<-(melt(pop[,"Inc",]))
    inclong$RealDate<-as.Date(inclong[,2], origin="1970-01-01")
    #summarize incident cases by year
    inclong$IterYear<-year(inclong$RealDate)
    inclong$AgeInYears<-floor((inclong[,1]-1)/12)
    res<-inclong%>%filter(IterYear>2000)%>%dplyr::group_by(IterYear,AgeInYears)%>%dplyr::summarize(Cases=sum(value))
    #split into under or over 30
    #over 30 generate a data frame with age 30-70, and cases sumingroup/41, for each yearly row
    over30<-res[res$AgeInYears==30,]
    over30expanded<-cbind(rep(over30$IterYear, 41),rep(30:70, each=100),rep(over30$Cases/41, 41))
    over30df<-as.data.frame(over30expanded)
    colnames(over30df)<-c("IterYear", "AgeInYears", "Cases")
    over30df$cfr<-cfr[6]
    #age-vector for death calc (under 30 only, over 30 is a constant)
    cfrVec<-c(rep(cfr[1],1), rep(cfr[2], 4),rep(cfr[3], 5), rep(cfr[4], 5), rep(cfr[5], 5), rep(cfr[6], 10))
    under30df<-cbind(as.data.frame(res[res$AgeInYears<30,]), "cfr"=rep(cfrVec, 100))
    results<-rbind(under30df, over30df)
    results$Deaths<-results$Cases*results$cfr
    results$DALYs<-results$Deaths*(70-results$AgeInYears) + (results$Cases-results$Deaths)* (0.26*0.072)*(70-results$AgeInYears)
    # res[k,"DALYs",] <- res[k,"death",] * (71-k) + ( res[k,"cases",]-res[k,"death",])* (0.26*0.072) * (71-k)
    results$simulation <- n
    my_data[[n]] <-results
    #cohort size - sum 2nd dimension except Inc / only need to do once - first simulation
    if (n==1) {
      cohort<-apply(pop[,1:8,], c(1,3), sum) # only a little slow
      cohortlong<-melt(cohort)
      cohortlong$RealDate<-as.Date(cohortlong[,2], origin="1970-01-01")
      cohortlong$IterYear<-year(cohortlong$RealDate)
      cohortlong$AgeInYears<-floor((cohortlong[,1]-1)/12)
      #pick last day of july per year for calculating cohort size  #WHERE MONTH(date) = 7 AND DAY(date) >= 25;
      cohortsample<-cohortlong[month(cohortlong$RealDate)==7 & day(cohortlong$RealDate)>24,]
      cohortsizes<-cohortsample%>%dplyr::group_by(IterYear, AgeInYears)%>%dplyr::summarize(cohortsize=sum(value))
      totalPop<-cohortsizes%>%dplyr::group_by(IterYear)%>%dplyr::summarize(tot=sum(cohortsize))
    }
  } #end of conditional summarization
  pop<-emptypop
  theDate <- start
  births <- myparams[myparams$year==year(theDate), 6]/52.1775
  imr <- myparams[myparams$year==year(theDate), 10]/(1000*52.1775)
  v <- myparams[myparams$year==year(theDate), 11] / (1000*52.1775)
  deathvec<-c(rep(imr,12), rep(v, 349))
  #only wanev depends on date via imr
  wanev <- c(rep(1-imr, 7), rep(0.000172, 17), rep(0.000096, 107), rep(0.000307, 230))  #waning from vacc to hi ab 
  j <- 2 #initialize time interval counter - the first position is filled with starting pop
  iterDate<-start-7 #this will make aging happen the first iteration, like SAS - also supply date for initial pop
  LastMonth <- month(iterDate)
  iterYear<-year(start-7)
  random<-runif(1)
  sto <- 1 + (phi * cos(random*pi))
} #end of looping nSims simulations
allsims <- rbindlist(my_data)
str(allsims)
if (nSims>1) {
  simsummary <-allsims%>%dplyr::select(IterYear, AgeInYears, Cases, Deaths, DALYs)%>%dplyr::group_by(IterYear, AgeInYears)%>%dplyr::summarize_all(.funs=(mean))
} else {
  simsummary<-allsims[, c("IterYear", "AgeInYears", "Cases", "Deaths", "DALYs")]
}
str(simsummary)
str(cohortsizes)
finalsummary<-merge(x = simsummary, y = cohortsizes, by = c("IterYear", "AgeInYears"))  
str(finalsummary)
head(finalsummary)
filename = paste0(mycountry, "_", program, "_", Sys.Date(), ".csv")
outfile<-paste0("\\\\HOME/stewcc1/MenAModel/Rdata/", filename)
write.csv(finalsummary, outfile)
print(begin)
#print(endsim)
print(Sys.time())
