

#### Program information ######################################################
# Source file name: MenA_plottingFunctions.R                                  #
#_______________________________________ _____________________________________#
# Input datasets: specify folder containing downloads from                    #
#   the output of MenA_VaccSims.R                                             #
#_____________________________________________________________________________#
# Parameters:                                                                 #
# Specify directory containing data to be plotted                             #
# 3-letter country code                                                       #       
#_____________________________________________________________________________#
# Purpose: Save a graph (various flavors) to the same directory as input data #
#_____________________________________________________________________________#
# Created  12/7/18, by Chris Stewart chris.c.stewart@kp.org                   #                                   
#_____________________________________________________________________________#

#_____________________________________________________________________________#
# Function: Get the most recent file with name matching a pattern             #
#_____________________________________________________________________________#
# Purpose: Given a directory and a pattern, return the name of the most recent# 
# file.  All plotting functions will use the most recent file                 #
#_____________________________________________________________________________#

GetMostRecentFile<-function(directory, mypattern) {
    filedf <- file.info(list.files(directory, pattern=mypattern, full.names = T))
    return(rownames(filedf)[which.max(filedf$mtime)])
}

PlotAllScenarios<-function(datadir, mycountry, saveit) {
  library(dplyr)
  library(data.table)
  library(ggplot2)
  
  allScenarios<-c("none","campaign","routine","both")
  #get the most recent file for every scenario available
  plot_files <- list()
  for (i in 1:length(allScenarios)) {
    patt<-paste0(mycountry, "_", allScenarios[i])
    myfile<-GetMostRecentFile(directory=datadir, mypattern=patt)
    plot_files[i]<-ifelse(length(myfile) > 0, myfile, NA)
  }
  mydata<-list()
  simcount<-list()
  for (j in 1:length(plot_files)){
    if (!is.na(plot_files[[j]])) {
      df<-read.csv(plot_files[[j]])
      vacc<-allScenarios[[j]]
      df$vaccination<-vacc
      mydata[[vacc]]<-df
    }
  }
  plotdata<-rbindlist(mydata)
  simcount<-plotdata%>%summarize(sim1=min(simulations), sim2=max(simulations))
  simstring<-ifelse(simcount$sim1==simcount$sim2, simcount$sim1, paste(simcount$sim1, simcount$sim2, sep="-"))
  colnames(plotdata)[colnames(plotdata)=="IterYear"] <- "year"
  YrSum<-plotdata%>%group_by(year,vaccination)%>%summarize(sumCases=sum(Cases))
  YrCum<-YrSum%>%group_by(year,vaccination)%>%mutate(cumCases=cumsum(sumCases))
  stitle<-paste("MenA modeling for", mycountry, ": ", simstring, " simulations" )
  plot1<-ggplot(data=YrCum, aes(x=year, y=cumCases, color=vaccination)) + geom_line() + ggtitle(stitle)
  if (saveit==TRUE) {
    filename<-paste0(datadir, "/", mycountry, "scenarios.png")
    ggsave(filename, plot = plot1)
  }
  return(plot1)
}


#tiff('\\\\home/stewcc1/MH_VDW/ICD10/papers/psychosis/Figures/Figure4ab.tiff', units="in", width=5, height=5, res=300)
#grid.arrange(plot4a, plot4b, ncol=1)
#dev.off()
PlotTenSims<-function(filename, stitle) {
  firstten<-read.csv(filename)
  firstten$factorsim<-as.factor(firstten$simulation)
  firstten<-firstten[, 3:5]
  firstten$program<-'R'
  colnames(firstten)[colnames(firstten)=="IterYear"] <- "Year"
  tenplot<-ggplot(data=firstten, aes(x=Year, y=sumCases, color=factorsim)) + geom_line() + ggtitle(stitle)
  tenplot
}
GridofTenSims<-function(datadir, mycountry, saveit) {
  #kind of messy but you can tell if it is working
  library(ggplot2)
  library(gridExtra)
  allScenarios<-c("none","campaign","routine","both")
  #get the most recent file for every scenario available
  plist<- list()
  for (k in 1:length(allScenarios)) {
    pattern<-patt<-paste0(mycountry, "_tensims_", allScenarios[k])
    myfile<-GetMostRecentFile(directory=datadir, mypattern=pattern)
    if (length(myfile) > 0) {
      mytitle<-paste0("1st 10 simulations for ", mycountry, "/",  allScenarios[k])
      detplot<-PlotTenSims(myfile, mytitle)
      vacc<-allScenarios[[k]]
      plist[[vacc]]<-detplot
    }
  }
  n <- length(plist)
  if (n > 0) {
    nCol <- floor(sqrt(n))
    if (saveit==TRUE) {
      fname<-paste0(datadir, "/Ten_Sim_Grid_for", mycountry, ".png")
      png(fname)
      do.call("grid.arrange", c(plist, ncol=nCol))
      dev.off()
    }
    else { thegrid<-do.call("grid.arrange", c(plist, ncol=nCol)) }
  }
}