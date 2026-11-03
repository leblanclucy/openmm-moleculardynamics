# purpose: look at which regions of the protein move with one another or not (correlation vs. anticorrelation)
# creates a dynamic cross correlation matrix and plots it as a heatmap

# install bio3D if you need to
install.packages('bio3d')

# load packages
library(bio3d)
library(tidyverse)
library(RColorBrewer)

# set path
setwd('path')

# load in trajectory and pdb file
# loading dcd file can take a while
# might also work with .cif but haven't tested it
dcd <- read.dcd('trajectory_centered.dcd')
pdb <- read.pdb('solvated_system.pdb')

# select alpha carbons and fit to trajectory
ca.inds <- atom.select(pdb, elety = 'CA')
xyz <- fit.xyz(fixed = pdb$xyz, mobile = dcd, fixed.inds = ca.inds$xyz, mobile.inds = ca.inds$xyz)

# compute dccm for the last 2500 frames
cij <- dccm(tail(xyz, 2500)[, ca.inds$xyz])

# get residue numbers and label matrix
res <- pdb$atom$resno[ca.inds$atom]
dimnames(cij) <- list(res, res)

# convert matrix into a long data frame, move the row names into a column, stack remaining columns into res_j and correlation
# one row of the result corresponds to one cell of the original matrix (res_i, res_j, correlation value)
cij_long <- as.data.frame(cij) |> rownames_to_column('res_i') |> pivot_longer(-res_i, names_to = 'res_j', values_to= 'correlation') |> mutate(res_i = as.numeric(res_i), res_j = as.numeric(res_j))

# plot the data, coloring the correlations such that red is +1, blue is -1, and pale yellow is 0
# set manual axis limits and also manual breaks every 50 residues
# coord_fixed makes it square shaped even if you stretch out the window
p <- ggplot(cij_long, aes(x=res_i, y=res_j, fill = correlation)) + geom_tile() + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")), limits = c(-1, 1)) 
+ coord_fixed(xlim=c(0, 300), ylim=c(0,300), expand = FALSE) + labs (x = "Residue ID", y = "Residue ID") + theme_classic() + scale_x_continuous(breaks = scales::breaks_width(50)) 
+ scale_y_continuous(breaks = scales::breaks_width(50))

# save the plot
ggsave('dccm_v1.eps', p, height=10, width=12, dpi=300, device=cairo_ps)
