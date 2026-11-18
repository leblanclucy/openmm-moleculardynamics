# load dependencies and read csv file
library(ggplot2)
library(tidyverse)
setwd('path')
getwd()
dfh <- read.csv('hbonds.csv')

# set order ahead of time
enzyme_order <- c('DmChi1', 'NsChi1', 'CfChi1', 'DaChi1')

# reshape using pivot_longer, mapping each value to its enzyme and residue
longh <- dfh |> pivot_longer(-replicate, names_to = c('enzyme', 'residue'), names_sep='_', values_to = 'occupancy') |> drop_na(occupancy) |> mutate(enzyme = factor(enzyme, levels=enzyme_order))

# mean and standard dev per bond, sorted by mean occupancy, going to test displaying with sd error bars first then dropping them if needed
summh <- longh |> group_by(enzyme, residue) |> summarize(mean = mean(occupancy), sd = sd(occupancy), .groups = 'drop') |> group_by(enzyme) |> mutate(slot=rank(mean, ties.method = 'first')) |> ungroup()

# create data points for geom_point ahead of time (makes the ggplot code block cleaner)
points <- left_join(longh, select(summh, enzyme, residue, slot), by = c('enzyme', 'residue'))

# set position dodge ahead of time to also make ggplot code block cleaner
dodge = position_dodge(width=0.9)

# plot as one large horizontal bar chart, ordered by enzyme, but not colored by enzyme (i will just label each part of the graph later)
p <- ggplot(summh, aes(x = enzyme, y = mean, group = slot)) + geom_col(fill='lightcyan1', color='black', width = 0.7, position = dodge, linewidth = 0.5) + geom_errorbar(aes(ymin = pmax(mean-sd,0), ymax=mean+sd), width=0.3, position = dodge) 
+ geom_point(data=points, aes(y = occupancy), color = 'grey20', position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.9), size = 1.5) + geom_text(aes(y = 0, label = residue), position=dodge, hjust = 1.25, size = 3) 
+ labs(x = NULL, y = 'HB Occupancy (%)') + theme_classic() + coord_flip() + scale_x_discrete(limits=rev) + scale_y_continuous(expand = expansion(mult = c(0.3, 0)), breaks = scales::breaks_width(25), limits = c(0,120))

# export and save as eps (vector-based)
ggsave(plot = p, 'hbonds_v1.eps', height=10, width=6, dpi=300, device=cairo_ps)
