#### Program information ######################################################
# Source file name: Initialize_population.R                                   #
#_____________________________________________________________________________#
# Input datasets: PopInputAgeDist.csv, dist_both.csv                          #
#_____________________________________________________________________________#
# Parameters:                                                                 #
# Start and end dates of simulation - to be specified in calling program      #
# Startpop from country params                                                #
#	country name to get age distributions from popinputagedist.csv              #
#	region = hyper or not hyper to get age-specific disease state distribution. #
#_____________________________________________________________________________#
# Purpose: fill first "slice" of the simulation matrix with the population at #   
# the starting point, divided by age and disease state. Matrix is 360 monthly #
# ages by 8 disease states.                                                   #
#_____________________________________________________________________________#
# Created as function 3/5/18, by Chris Stewart stewart.c@ghc.org              #
# Changes: remove loop, use vector multiplication instead,                    #
# 30 age group only needs one value in input vectors                          #           
#_____________________________________________________________________________#

InitializePopulation<-function(scriptdir, inputdir, start, end, popsize, country="ETH", region="not_hyper") {
  #create the matrix
  dur<- ceiling(difftime(end, start, units = "weeks"))
  poparray <- array(data=0, dim=c(361, 9, dur+1)) # Dimensions are [age groups, states, time]
  dimnames(poparray)[[2]] <- c("Ns","Nc","Ls","Lc","Hs","Hc","Di","Va", "Inc")
  #parameters for initializing
  mypop<-GetPopAgeDist(path=inputdir, mycountry=country, start=start) 
  #setwd("\\\\HOME/stewcc1/MenAModel/data/ModelInputs")
  #inpop<-read.csv("PopInputAgeDist.csv")
  #mypop1<-inpop[inpop$Country==country,]
  ##get age-specific proportions of each disease state into vectors: 7 ages X 7 disease states
  statefract<-GetDiseaseStateDist(directory=scriptdir, region=myregion)
  #dist<-read.csv("dist_both.csv", stringsAsFactors = TRUE)
  #distcol<-ifelse(region=='hyper', 4, 3)
  #statefract<-as.vector(dist[,distcol]) # fraction of each disease state in each of 7 population groups
  #expand age group fraction as vector to match pop matrix dimension 1
  agefract <- c(rep(as.numeric(mypop[1:6]), each=60), as.numeric(mypop[8]))
  chunks <- c(rep(60, each=360), 1)
 
  statemx<-rbind(
    matrix(rep(statefract[1:7], each=60), nrow=60),
    matrix(rep(statefract[8:14], each=60), nrow=60),
    matrix(rep(statefract[15:21], each=60), nrow=60),
    matrix(rep(statefract[22:28], each=60), nrow=60),
    matrix(rep(statefract[29:35], each=60), nrow=60),
    matrix(rep(statefract[36:42], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=1), nrow=1)
  )
  #initialize - 3rd dimension stays at 1
    poparray[, "Ns", 1]<-(startSize*agefract/chunks)*statemx[,1]
    poparray[, "Nc", 1]<-(startSize*agefract/chunks)*statemx[,2]
    poparray[, "Ls", 1]<-(startSize*agefract/chunks)*statemx[,3]
    poparray[, "Lc", 1]<-(startSize*agefract/chunks)*statemx[,4]
    poparray[, "Hs", 1]<-(startSize*agefract/chunks)*statemx[,5]
    poparray[, "Hc", 1]<-(startSize*agefract/chunks)*statemx[,6]
    poparray[, "Di", 1]<-(startSize*agefract/chunks)*statemx[,7]
    #leave Vacc, Inc 0
    return(poparray)
}