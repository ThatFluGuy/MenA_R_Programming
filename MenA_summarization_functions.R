#### Program information ######################################################
# Source file name: MenA_summarization_functions.R                            #
#_____________________________________________________________________________#
# Input datasets: none                                                        #
#_____________________________________________________________________________#
# Functions called by MenA_VaccSim                                            #
# Contents:                                                                   #
# getCohortSize: calulate total poplation in each year and year of age        #
# summarizeOneSim: Collapse results of one simulation into a data frame       #
#   with Cases, Deaths, Dalys by year and year of age.                        #
#	SummarizeForOutput : takes a list of products of summarizeOneSim and        #  
#   calculates mean over all simulations; writes output file                  #
#_____________________________________________________________________________#
#_____________________________________________________________________________#
# Created as functiosn 3/8/18, by Chris Stewart stewart.c@ghc.org             #
# Changes:                                                                    #
#_____________________________________________________________________________#

getCohortSize<-function(poparray) {
  #cohort size - sum 2nd dimension except Inc / only need to do once - first simulation
  #modify 3/14/18 - need to split 30+ pot into ages 30-70
    cohort<-apply(poparray[,1:8,], c(1,3), sum) # only a little slow
    cohortlong<-melt(cohort)
    cohortlong$RealDate<-as.Date(cohortlong[,2], origin="1970-01-01")
    cohortlong$IterYear<-year(cohortlong$RealDate)
    cohortlong$AgeInYears<-floor((cohortlong[,1]-1)/12)
    #pick last day of july per year for calculating cohort size  #WHERE MONTH(date) = 7 AND DAY(date) >= 25;
    cohortsample<-cohortlong[month(cohortlong$RealDate)==7 & day(cohortlong$RealDate)>24,]
    coh_under30<-cohortsample[cohortsample$AgeInYears!=30,]
    coh_over30<-cohortsample[cohortsample$AgeInYears==30,]
    coh_over30exp<-cbind(rep(coh_over30$IterYear, 41),rep(30:70, each=100),rep(coh_over30$value/41, 41))
    coh_over30expdf<-as.data.frame(coh_over30exp)
    colnames(coh_over30expdf)<-c("IterYear", "AgeInYears", "cohortsize")
    sum_under30<-coh_under30%>%group_by(IterYear, AgeInYears)%>%summarize(cohortsize=sum(value))
    cohortsizes<-rbind(as.data.frame(sum_under30), coh_over30expdf)
    return(cohortsizes)
}

summarizeOneSim<-function(poparray, n, cfr) {
  #summarize incident cases by year and year of age, calculate deaths and DALYs
  inclong<-(melt(poparray[,"Inc",]))
  inclong$RealDate<-as.Date(inclong[,2], origin="1970-01-01")
  #summarize incident cases by year
  inclong$IterYear<-year(inclong$RealDate)
  inclong$AgeInYears<-floor((inclong[,1]-1)/12)
  res<-inclong%>%filter(IterYear>2000)%>%group_by(IterYear,AgeInYears)%>%summarize(Cases=sum(value))
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
  return(results)
}

summarizeForOutput<-function(results_list, cohort, write, filename) {
  if (length(results_list) > 1) {
    allsims <- rbindlist(results_list)
    simsummary <-allsims%>%select(IterYear, AgeInYears, Cases, Deaths, DALYs)%>%group_by(IterYear, AgeInYears)%>%summarize_all(.funs=(mean))
  } else {
    allsims<-results_list
    simsummary<-allsims[, c("IterYear", "AgeInYears", "Cases", "Deaths", "DALYs")]
  }
  #need to join for cohort size
  finalsummary<-merge(x = simsummary, y = cohort, by = c("IterYear", "AgeInYears")) 
  if (write==TRUE){
    write.csv(finalsummary, filename)
    print(paste("Output written to", filename))
  }
  return(finalsummary)
}