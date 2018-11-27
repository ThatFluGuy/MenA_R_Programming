#### Program information ######################################################
# Source file name: MenA_OneSim.R                                             #
#_____________________________________________________________________________#
# Input datasets: none                                                        #
#_____________________________________________________________________________#
# Parameters:                                                                 #
# Start and end dates of simulation - to be specified in calling program      #
# population array to fill                                                    #
#	data frame containing country parameters                                    #
#	WAIFW matrix (3d with dry vs rainy?)                                        #
# dd1&2, dr1&2 per-carrier risk of invasive disease (age-specific)            #
#_____________________________________________________________________________#
# Purpose: One simulation - iterate weekly from start date to end date        #   
# Move people per WAIFW weekly, age monthly, vaccinate monthly as speficiednnn#
# fill pop matrix slices 2 thru n                                             #
#_____________________________________________________________________________#
# Created as function 3/6/18, by Chris Stewart stewart.c@ghc.org              #
# Changes:                                                                    #
#   WAIFW and dxrisk as small 3D matrices                                     #           
#_____________________________________________________________________________#

MenASimulation<-function(startdt, enddt, pop, fixedparams, countryparams, WAIFWmx, dxr) {
  #setup before loop
  theDate <- start
  births <- countryparams[countryparams$year==year(theDate), 6]/52.1775
  imr <- countryparams[countryparams$year==year(theDate), 10]/(1000*52.1775)
  v <- countryparams[countryparams$year==year(theDate), 11] / (1000*52.1775)
  deathvec<-c(rep(imr,12), rep(v, 349))
  #waning age-group vectors 0-5mo, 6mo-2y, 3-10y, 11+y
  wanev <- c(rep(1-imr, 7), rep(0.000172, 17), rep(0.000096, 107), rep(0.000307, 230))  #waning from vacc to hi ab scaled to wks-confirm wv(1) = NA
  waneh <- c(rep(0.01092, 6), rep(0.00654, 18), rep(0.00527, 107), rep(0.00096, 230)) #waning from high to low ab scaled to weeks
  wanel <- c(rep(0.00970,6), rep(0.00487, 18), rep(0.00364, 107), rep(0.00057, 230)) #waning from low 
  j <- 2 #initialize time interval counter - the first position is filled with starting pop
  iterDate<-start-7 #this will make aging happen the first iteration, like SAS - also supply date for initial pop
  LastMonth <- month(iterDate)
  iterYear<-year(start-7)
  random<-runif(1)
  sto <- 1 + (phi * cos(random*pi))
  #date loop
  while (theDate <= enddt) {
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
      births <-countryparams[countryparams$year==year(theDate), "births"]/52.1775
      imr <- countryparams[countryparams$year==year(theDate), "imr"]/(1000*52.1775)
      v <- countryparams[countryparams$year==year(theDate), "v"] / (1000*52.1775)
      #new imr, need to update wanev and death
      wanev <- c(rep(1-imr, 7), wanev[8:361]) 
      deathvec<-c(rep(imr,12), rep(v, 349))
    }
    
    #  waifw matrix depends on rainy (Mar-Aug) or dry (Sep-Feb) season
    if (month(theDate)!=LastMonth)
    {
      if (month(theDate) %in% c(3,4,5,6,7,8)) { wmx <- WAIFWmx[,,"rainy"]
      } else { wmx <- WAIFWmx[,,"dry"]
      }
      # Determine age-specific per-carrier risk of invasive disease;
      #like infection, this is also seasonal but peaks later, hence different months (rainy here = 6-10 June-Oct)
      iagevec<-seq(0, 360)
      if (month(theDate) %in% c(1,2,3,4,5,11,12)) {sigma = sigmavec<-dxr["dry",1] + dxr["dry",2] *iagevec
      } else {sigmavec<-dxr["rainy",1] + dxr["rainy",2] *iagevec
      }
    } #end of updates conditional on month
    
    #infectious ratios by pop age group (for calculating force) this happens every time point
    IR<-GetInfectiveRatio(t(pop[,,j-1]))
    # Includes force of infection from outside the population (foii, was for in SAS code)
    pf<-fixedparams
    forcevec <- (sto*wmx[,1]*IR[1]) + (sto*wmx[,2]*IR[2]) + (sto*wmx[,3]*IR[3])+ (sto*wmx[,4]*IR[4]) + pf["foii"]
    
    #vectorized (1-361) calculations instead of loop
    pop[,"Ns",j] <- wanel*pop[,"Ls",j-1] - (deathvec+forcevec)*pop[,"Ns",j-1] + pop[,"Ns",j-1]
    #births - this order is the same as sas (could just add births to line above)
    pop[1,"Ns",j] <- pop[1,"Ns",j] + births
    pop[,"Nc",j] <- forcevec*pop[,"Ns",j-1] - (deathvec+pf["rc"]+sigmavec) * pop[,"Nc",j-1] + pop[,"Nc",j-1]
    pop[,"Ls",j] <- waneh*pop[,"Hs",j-1] + pf["rc"]*pop[,"Nc",j-1] - (deathvec+wanel+(1-pf["lc"])*forcevec) * pop[,"Ls",j-1] + pop[,"Ls",j-1]
    pop[,"Lc",j] <- (1-pf["lc"])*forcevec*pop[,"Ls",j-1] - (deathvec+pf["rc"]+(1-pf["ld"])*sigmavec) * pop[,"Lc",j-1]  + pop[,"Lc",j-1]
    pop[,"Hs",j] <- pf["rc"]*(pop[,"Lc",j-1] + pop[,"Hc",j-1]) + pf["rd"]*pop[,"Di",j-1]  + wanev*pop[,"Va",j-1] - (deathvec+waneh+(1-pf["hc"])*forcevec) * pop[,"Hs",j-1] + pop[,"Hs",j-1]
    pop[,"Hc",j] <- (1-pf["hc"])*forcevec*pop[,"Hs",j-1] - (deathvec+pf["rc"]-(1-pf["hd"])*sigmavec)*pop[,"Hc",j-1] + pop[,"Hc",j-1]
    pop[,"Di",j] <- sigmavec*(pop[,"Nc",j-1] + (1-pf["ld"])*pop[,"Lc",j-1] + (1-pf["hd"])*pop[,"Hc",j-1]) - (deathvec+pf["rd"])*pop[,"Di",j-1] + pop[,"Di",j-1]
    pop[,"Va",j] <- -(deathvec+wanev)* pop[,"Va",j-1] + pop[,"Va",j-1]
    # Count incident cases of invasive disease NOTE THESE ARE before incrementing = use old pop
    pop[,"Inc",j] <- sigmavec*(pop[,"Nc",j-1] + (1-pf["ld"])*pop[,"Lc",j-1] + (1-pf["hd"])*pop[,"Hc",j-1])
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
    # if (year(theDate) > 2020) {break}
    LastMonth <- month(theDate)
    theDate<-theDate+7
    j <- j + 1
  } #end of date loop
  #label matrix iteration dimension with dates
  dimnames(pop)[[3]] <- as.Date(iterDate)
  return(pop)
}