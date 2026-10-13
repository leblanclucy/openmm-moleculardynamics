# graphs RMSF (residue on x axis vs RMSF on y axis) showing residue flexibility over a molecular dynamics trajectory
# uses multiple replicates and displays standard deviation to show spread in data
install.packages("ggplot2")
install.packages("tidyverse")
library(ggplot2)
library(tidyverse)

setwd('path')
getwd()

# load your RMSF data and set column names. i had 7 columns, the first being Residue ID, the others being replicates of my two conditions.
df <- read.csv("rmsf.csv")
colnames(df) <- c("Residue_ID", "condition1_rep1", "condition2_rep2")

# make your data from into long format where every row has info on condition, replicate, and RMSF
long <- df |> pivot_longer(-Residue_ID, names_to=c("condition", "replicate"), names_sep="_", values_to="rmsf")

# calculate mean and standard deviation so that you can graph both
# i floored the lower bound of standard deviation at zero since you cannot have negative RMSF
summ <- long |> group_by(condition, Residue_ID) |> summarize(mean=mean(rmsf), sd=sd(rmsf), n=n(), .groups="drop") |> mutate(lower = pmax(mean-sd, 0), upper = mean+sd)

# graph the mean as lines (geom_line) and the standard deviation as shading in the background (geom_ribbon) indicating the spread of the data
# alter the transparency (alpha) and shading color (fill) as desired. overlapping standard deviations are darker grey.
# manually set line colors as desired
ggplot(summ, aes(x=Residue_ID, y=mean))+geom_ribbon(aes(ymin = lower, ymax = upper, group = condition), fill = "grey30", alpha=0.45)+geom_line(aes(color=condition), linewidth=0.75)+scale_color_manual(values = c(condition1="firebrick", condition2="black"))+labs(x = "Residue Index", y="RMSF (Å)", color=NULL)+theme_classic()`

# save as eps in your working directory
# keep hxw consistent if making multiple graphs 
# use this device (cairo_ps) or else transparency won't save.
ggsave("rmsf.eps", height=5.21, width=7.91, dpi=300, device=cairo_ps) 
