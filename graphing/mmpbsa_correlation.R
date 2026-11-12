# set working directory and confirm
setwd('path/to/data/')
getwd()

# import dependencies, install if needed
library(ggplot2)
library(tidyverse)

# import csv containing data, where the first column is x, and the second is y
dfx <- read.csv('mmgbpsa_corr.csv')

# graph it using ggplot, adding a linear trend and confidence interval
# uses y ~ x as the regression method but you can do polynomial regression
# example: formula <- y ~ poly(x, 2, raw=TRUE), then instead of method=lm use method=formula
# se=TRUE means show a 95% confidence interval around the regression line
# also set x/y limits and breaks
p <- ggplot(dfx, aes(x=GBSA, y=PBSA)) + geom_smooth(method=lm, color='slateblue4', se=TRUE, alpha=0.5) + geom_point(size = 3) + theme_classic() 
+ coord_cartesian(xlim = c(-90,-25), ylim = c(-90,-25)) + stat_poly_eq(method=lm, geom='text', use_label('eq', 'r2'), label.x=-60, label.y=-40) 
+ scale_x_continuous(breaks = scales::breaks_width(15)) + scale_y_continuous(breaks = scales::breaks_width(15))

"""
for polynomial regression
ggplot(dfx, aes(x=GBSA, y=PBSA)) + geom_smooth(method=lm, formula=formula, color='slateblue4', se=TRUE, alpha=0.5) + geom_point(size = 3)
+ theme_classic() + coord_cartesian(xlim = c(-90,-20), ylim = c(-90,-20)) + stat_poly_eq(method=lm, formula=formula, geom='text', use_label('eq', 'r2'), label.x=-60, label.y=-40)
"""

# export, keep the same hxw for all graphs for each figure, use cairo_ps to save with transparency
ggsave(p, 'mmgpbsa_corr_v1.eps', height=10, width=10, dpi=300, device=cairo_ps)
