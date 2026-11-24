# after MM/GBSA (or MM/PBSA) per residue decomposition, this plots the individual binding energy of each residue as a line graph so you can identify hotspots for binding

# load packages
library(ggplot2)
library(tidyverse)

# read the csv where the first column is your residue numbers and each column after that is the delta H of each individual residue
dfg5 <- read.csv('decomp_os.csv')

# only needed if you don't have column names
colnames(dfg5) <- c('Residue', 'seed13', 'seed14', 'seed15', 'seed16')

# pivot longer
longg5 <- dfg5 |> pivot_longer(-Residue, names_to=c('replicate'), values_to='affinity')

# plot replicates separately using ggplot
# the more replicates you have, the more the lines overlap and so it can help to adjust the line width of the final replicates' lines in illustrator
# be very mindful of scale_y_continuous, make sure you aren't excluding any data
ggplot(longg5, aes(x=Residue, y=affinity, color=replicate))+geom_line(linewidth=0.75)+scale_color_manual(values=c(seed13='black',seed14='firebrick2',seed15='royalblue1',seed16='gray40'))+labs(x = "Residue Index", y="Affinity (kcal mol-1)", color=NULL)+theme_classic()+scale_y_continuous(limits=c(-10,3),breaks=scales::breaks_width(3))+scale_x_continuous(limits=c(0,320),breaks=scales::breaks_width(40))

# export the plot, use same dimensions for all line graphs
ggsave('decomp_os.eps', height=10, width=12, dpi=300, device=cairo_ps) 
