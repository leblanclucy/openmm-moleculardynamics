# set working directory and confirm
setwd('path/to/data/')
getwd()

# import dependencies, install if needed
library(ggplot2)
library(tidyverse)

# import csv containing data, where column names are secondary structures, and the first column has all replicates of each enzyme
dfd <- read.csv('data.csv')

# set column names for all secondary structure elements
colnames(dfd) <- c('replicate', 'ahelix', 'bbridge', 'bstrand', '310helix', 'pihelix', 'turn', 'bend', 'loop')

# use pivot_longer to reformat dataframe
# mutate here allows you to preserve the order of variables in the order given rather than R reordering alphabetically
longd <- dfd |> pivot_longer(-replicate, names_to = 'element', values_to = 'percent') |> separate(replicate, into = c('enzyme', 'replicate'), sep = '_') |> mutate(element = factor(element, levels = unique(element)), enzyme = factor(enzyme, levels = unique(enzyme)))

# use summarize to get mean and standard deviation
# i used mutate here to precalculate the upper and lower bounds but you can also do it in geom_errorbar
summd <- longd |> group_by(enzyme, element) |> summarize(mean=mean(percent), sd=sd(percent), n=n(), .groups='drop') |> mutate(lower = mean-sd, upper = mean+sd)

# plot using ggplot, graphing the means +/- the standard deviations, color coding both by enzyme, grouping by secondary structure element
# individual data points AND error bars are both graphed which is the standard for showing reproducibility
ggplot(summd, aes(x = element, y = mean, fill = enzyme)) + geom_col(position = position_dodge(width=0.9), color = 'black', linewidth = 0.5) 
+ geom_errorbar(aes(ymin=lower, ymax=upper), position = position_dodge(width=0.9), width = 0.5) 
+ geom_point(data = longd, aes(x = element, y = percent, fill = enzyme), color = 'grey20', position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), size = 1.5, show.legend = FALSE) 
+ scale_fill_manual(values = c(dm='blue4', ns='green3', cf='orange1', da='mediumorchid1')) 
+ labs(x = 'Secondary Structure Element', y = 'Percent (%)', fill = 'Enzyme') + theme_classic() + coord_cartesian(ylim = c(0,40))
             
# export, keep the same hxw for all graphs for each figure, use cairo_ps to save with transparency
ggsave('dssp_v1.eps', height=6, width=13, dpi=300, device=cairo_ps)
