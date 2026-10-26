# set working directory and confirm
setwd('path/to/data/')
getwd()

# import dependencies, install if needed
library(ggplot2)
library(tidyverse)

# import csv containing data, where column names are replicates for each enzyme, plus the first column is time in ns
dfr <- read.csv("data.csv")

# set column names, for however many enzymes you have
colnames(dfr) <- c("Time_ns", "enzyme1_seed13")

# use pivot_longer to reformat dataframe
longr <- dfr |> pivot_longer(-Time_ns, names_to=c("condition", "replicate"), names_sep="_", values_to="Rg")

# use summarize to get mean and standard deviation
summr <- longr |> group_by(condition, Time_ns) |> summarize(mean=mean(Rg), sd=sd(Rg), n=n(), .groups="drop") |> mutate(lower = pmax(mean-sd), upper = mean+sd)

# plot using ggplot, graphing the means +/- the standard deviations, color coding both by enzyme
# may want to set a y min that's non zero depending on how tightly clustered your data points are around a certain y value over the whole trajectory
ggplot(summr, aes(x=Time_ns, y=mean))+geom_ribbon(aes(ymin = lower, ymax = upper, group = condition), fill = condition, alpha=0.35)+geom_line(aes(color=condition), linewidth=0.75)
+scale_color_manual(values = c(enzyme1="blue4", enzyme2="green3", enzyme3="orange1", enzyme4="mediumpurple4"))+scale_fill_manual(values = c(enzyme1="blue4", enzyme2="green3", enzyme3="orange1", enzyme4="mediumpurple4"))
+labs(x = "Time (ns)", y="Rg (Å)", color=NULL)+theme_classic()+coord_cartesian(ylim = c(15,25))

# export, keep the same hxw for all graphs for each figure, use cairo_ps to save with transparency
ggsave("rg_v1.eps", height=5.21, width=7.91, dpi=300, device=cairo_ps)
