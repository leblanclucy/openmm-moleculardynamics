# set working directory and confirm
setwd('path/to/data/')
getwd()

# import dependencies, install if needed using install.packages("tidyverse") and install.packages("ggplot2")
library(ggplot2)
library(tidyverse)

# import csv containing data, where column names are replicates for each enzyme, plus the first column is time in ns
# pref using single quotes because " can give an error if curly
dfrmsd1 <- read.csv('rmsd_data.csv')

# set column names, for however many enzymes you have
# make sure that your time starts at 0 ns!
colnames(dfrmsd1) <- c('Time_ns', 'enzyme1_seed13')

# use pivot_longer to reformat dataframe
long1 <- dfrmsd1 |> pivot_longer(-Time_ns, names_to=c('condition', 'replicate'), names_sep='_', values_to='RMSD')

# use summarize to get mean only
summ1 <- long1 |> group_by(condition, Time_ns) |> summarize(mean=mean(RMSD), n=n(), .groups="drop")

# plot using ggplot, graphing the means, color coding by whole protein vs. catalytic domain
# replicates are shown separately as slightly thinner lines with lower opacity and different line style
# adjust ylim as needed, for RMSD it should start at 0 because that's basically the starting point
ggplot()+geom_line(data=long1, aes(x = Time_ns, y = RMSD, color = condition, linetype = replicate), linewidth = 0.65, alpha = 0.45)
+geom_line(data=summ1, aes(x = Time_ns, y = mean, color = condition), linewidth = 0.75)+scale_color_manual(values = c(whole='black', cat='maroon')) 
+ scale_linetype_manual(values = c(seed13='solid', seed14='dotted', seed15='longdash'))
+labs(x = 'Time(ns)', y='RMSD (Å)', color = 'Condition', linetype = 'Replicate')+theme_classic()+coord_cartesian(ylim = c(0,16))

# export, keep the same hxw for all graphs for each figure, use cairo_ps to save with transparency
ggsave('rmsd_v1.eps', height=5.21, width=7.91, dpi=300, device=cairo_ps)
