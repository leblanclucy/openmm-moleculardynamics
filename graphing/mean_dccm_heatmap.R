# purpose: look at which regions of the protein move with one another or not (correlation vs. anticorrelation)
# creates a dynamic cross correlation matrix and plots it as a heatmap
# takes the mean of replicates to make one averaged dccm

# load packages
library(bio3d)
library(tidyverse)
library(RColorBrewer)

# note this assumes your replicates live in different directories but their trajectory files etc have the same name
# create list of directories corresponding to each replicate
dirs <- c('rep1directory', 'rep2directory', 'rep3directory')

# loop through directories
# read the dcd file, pdb file, select alpha carbons, fit to trajectory, compute dccm on final 2500 frames
# store all of that per replicate in one variable, cij_list
cij_list <- lapply(dirs, function(d) {
  dcd <- read.dcd(file.path(d, "trajectory_centered.dcd"))
  pdb <- read.pdb(file.path(d, "solvated_system.pdb"))
  ca.inds <- atom.select(pdb, elety = "CA")
  xyz <- fit.xyz(pdb$xyz, dcd, fixed.inds = ca.inds$xyz, mobile.inds = ca.inds$xyz)
  dccm(tail(xyz, 2500)[, ca.inds$xyz])
})

# compute mean by averaging the matrices by element
cij_mean <- Reduce(`+`, cij_list) / length(cij_list)

# get residue numbers and label matrix
# reload pdb, ca.inds, and res variables because they are lost after the loop is done
setwd('path')
pdb <- read.pdb('solvated_system.pdb')
ca.inds <- atom.select(pdb, elety = "CA")
res <- pdb$atom$resno[ca.inds$atom]
dimnames(cij_mean) <- list(res, res)

# convert matrix into a long data frame, move the row names into a column, stack remaining columns into res_j and correlation
# one row of the result corresponds to one cell of the original matrix (res_i, res_j, correlation value)
cij_long <- as.data.frame(cij_mean) |> rownames_to_column('res_i') |> pivot_longer(-res_i, names_to = 'res_j', values_to= 'correlation') |> mutate(res_i = as.numeric(res_i), res_j = as.numeric(res_j))

# plot the data, coloring the correlations such that red is +1, blue is -1, and pale yellow is 0
# set manual axis limits and also manual breaks every 50 residues
# change the axis limits and breaks based on the size of the protein
# coord_fixed makes it square shaped even if you stretch out the window
p <- ggplot(cij_long, aes(x=res_i, y=res_j, fill = correlation)) + geom_tile() + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")), limits = c(-1, 1)) 
+ coord_fixed(xlim=c(0, 300), ylim=c(0,300), expand = FALSE) + labs (x = "Residue ID", y = "Residue ID") + theme_classic() + scale_x_continuous(breaks = scales::breaks_width(50)) 
+ scale_y_continuous(breaks = scales::breaks_width(50))

# save the plot
ggsave('mean_dccm_v1.eps', p, height=10, width=12, dpi=300, device=cairo_ps)
