# set working directory and confirm
setwd('path/to/data/')
getwd()

# import dependencies, install if needed
library(ggplot2)
library(tidyverse)

# import csv containing data, where column names are secondary structures, and the first column has all replicates of each enzyme
dfm <- read.csv('data.csv')

# set column names for both mm/gbsa and mm/pbsa since i calculated both, plus replicate indicating enzyme/seed
# not needed if csv already has appropriate column names
colnames(dfm) <- c('replicate', 'mmgbsa', 'mmpbsa')

# use pivot_longer to reformat dataframe
# mutate here allows you to preserve the order of variables in the order given rather than R reordering alphabetically
# binding_energy indicates whether it's gbsa or pbsa
longm <- dfm |> pivot_longer(-replicate, names_to = 'binding_energy', values_to = 'affinity') |> separate(replicate, into = c('enzyme', 'replicate'), sep = '_') |> mutate(binding_energy = factor(binding_energy, levels = unique(binding_energy)), enzyme = factor(enzyme, levels = unique(enzyme)))

# use summarize to get mean and standard deviation
# i used mutate here to precalculate the upper and lower bounds but you can also do it in geom_errorbar
summm <- longm |> group_by(enzyme, binding_energy) |> summarize(mean=mean(affinity), sd=sd(affinity), n=n(), .groups='drop') |> mutate(lower = mean-sd, upper = mean+sd)

# plot using ggplot, graphing the means +/- the standard deviations, color coding both by enzyme, grouping by either MM/GBSA or MM/PBSA
# individual data points AND error bars are both graphed which is the standard for showing reproducibility
ggplot(summm, aes(x = binding_energy, y = mean, fill = enzyme)) + geom_col(position = position_dodge(width=0.9), color = 'black', linewidth = 0.5) 
+ geom_errorbar(aes(ymin=lower, ymax=upper), position = position_dodge(width=0.9), width = 0.5) 
+ geom_point(data = longm, aes(x = binding_energy, y = affinity, fill = enzyme), color = 'grey20', position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.9), size = 3, show.legend = FALSE) 
+ scale_fill_manual(values = c(dm='dodgerblue4', ns='green3', cf='orange1', da='mediumorchid1')) + labs(x = 'Binding Affinity Method', y = 'Affinity (kcal mol-1)', fill = 'Enzyme') + theme_classic() + coord_cartesian(ylim = c(-100,0))
             
# export, keep the same hxw for all graphs for each figure, use cairo_ps to save with transparency
ggsave('mmgpbsa_bar_v1.eps', height=10, width=12, dpi=300, device=cairo_ps)
