### Compare scenarios in the "constant force of infection" model ##############

library(ggplot2)

output.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Analysis/ConstantFOI"

novax.raw <- read.csv(paste0(output.dir, "/ConstantFOI-none.csv"), stringsAsFactors = FALSE)
campaign.raw <- read.csv(paste0(output.dir, "/ConstantFOI-campaign.csv"), stringsAsFactors = FALSE)
routine.raw <- read.csv(paste0(output.dir, "/ConstantFOI-routine.csv"), stringsAsFactors = FALSE)

# Combine datasets and collapse by scenario and year

novax.raw$scenario <- "none"
campaign.raw$scenario <- "campaign"
routine.raw$scenario <- "routine"

combined.raw <- rbind(novax.raw, campaign.raw, routine.raw)


combined.df <- aggregate(cbind(Cases, Pop)~scenario+year, data=combined.raw, FUN=sum)

# Plot 

p1 <- ggplot(data=combined.df, aes(x=year, y=Cases, color=scenario)) + theme_bw() +
  theme(panel.grid = element_blank())

p1 + geom_line() +
  scale_color_manual(values=c("red", "black", "green")) +
  scale_y_continuous(limits=c(0, 6500))