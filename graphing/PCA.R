# does principal component analysis on MD trajectories using Bio3D
# pools replicates, last 50 ns of each trajectory
# uses catalytic domain (CD) only, change for each enzyme
library(bio3d)
library(ggplot2)

setwd('path')
getwd()

dirs_cf <- c('list', 'directories', 'here')

# define extent of whatever part of the protein you're interested in and load a reference pdb and list of alpha carbons in that domain
catdomain <- 57:300
ref_pdb <- read.pdb(file.path(dirs_cf[1], "solvated_system.pdb"))
ref_ca  <- atom.select(ref_pdb, "calpha", resno = catdomain)
head(ref_pdb$atom[ref_ca$atom, c("resno", "resid", "elety")])

# read dcd files in the directories supplied, only using the last 50 ns (for me, it was the last 2500 frames)
xyz_list <- lapply(dirs_cf, function(d) {
  dcd <- read.dcd(file.path(d, "trajectory_centered.dcd"))
  pdb <- read.pdb(file.path(d, "solvated_system.pdb"))
  ca  <- atom.select(pdb, "calpha", resno = catdomain)
  xyz <- fit.xyz(fixed = ref_pdb$xyz, mobile = dcd,
                 fixed.inds = ref_ca$xyz, mobile.inds = ca$xyz)
  tail(xyz, 2500)[, ca$xyz]
})

# verify only the final X number of frames have been loaded in
sapply(xyz_list, nrow)

# combines replicates into one matrix
xyz_all <- do.call(rbind, xyz_list)

# optional, but this keeps track of which replicate is which in the matrix
rep_id  <- rep(seq_along(xyz_list), times = sapply(xyz_list, nrow))

# computes PCA on the matrix containing all replicates
pc <- pca.xyz(xyz_all)

# stores the scores to be plotted by ggplot2 or other tools
scores <- data.frame(PC1 = pc$z[, 1], PC2 = pc$z[, 2], PC3 = pc$z[, 3], replicate = factor(rep_id))

# calculate the variance explained by each principal component
var_exp <- round(100 * pc$L / sum(pc$L), 1)

# graphs the scree plot as well as the first few principal components against one another
plot(pc)

# plot using ggplot2, coloring by replicate
ggplot(scores, aes(x = PC1, y = PC2, color = replicate)) + stat_density_2d(alpha = 0.45, linewidth = 0.5) + geom_point(alpha = 0.25, size = 1.2) + geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) + geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) + labs(x = paste0("PC1 (", var_exp[1], "% variance)"), y = paste0("PC2 (", var_exp[2], "% variance)"), color = "Replicate") + theme_classic()

# to check equilibration of large scale motions over time, plot pc1 vs pc2 with frame as the color, where later frames are yellow and initial frames are dark blue
scores$frame <- unlist(lapply(sapply(xyz_list, nrow), seq_len))
ggplot(scores, aes(x = PC1, y = PC2, color = frame)) + geom_point(alpha = 0.5, size = 1.2) + facet_wrap(~ replicate, labeller = label_both) + scale_color_viridis_c(option = "plasma", name = "Frame") + labs(x = paste0("PC1 (", var_exp[1], "% variance)"), y = paste0("PC2 (", var_exp[2], "% variance)")) + theme_bw()

# if you want to graph in 3D, need to use a different tool like plotly. this makes an interactive 3d plot that you can rotate in rstudio
plot_ly(data = scores, x = ~PC1, y = ~PC2, z = ~PC3, color = ~replicate, type = "scatter3d", mode = "markers", marker = list(size = 2.5, opacity = 0.65)) %>% layout(scene = list(xaxis = list(title = paste0("PC1 (", var_exp[1], "%)")), yaxis = list(title = paste0("PC2 (", var_exp[2], "%)")), zaxis = list(title = paste0("PC3 (", var_exp[3], "%)"))))
