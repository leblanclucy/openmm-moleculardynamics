# load dependencies and read csv file
library(ggplot2)
library(tidyverse)
setwd('path')
getwd()
dfh <- read.csv('hbonds.csv')

# reshape using pivot_longer
longh <- dfh |> pivot_longer(-replicate, names_to = c('enzyme', 'residue'), names_sep='_', values_to = 'occupancy') |> drop_na(occupancy) |> mutate(residue = residue, enzyme = factor(enzyme))

# calculate mean and sd
summh <- longh |> group_by(enzyme, residue) |> summarize
