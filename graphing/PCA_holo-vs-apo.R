# does principal component analysis on MD trajectories using Bio3D
# pools replicates, last 50 ns of each trajectory
# uses catalytic domain (CD) only, change for each enzyme
# compares holo vs apo
library(bio3d)
library(ggplot2)

setwd('path')
getwd()

dirs_holo <- c('list', 'directories', 'here')
dirs_apo <- c('list', 'directories', 'here')

# define extent of whatever part of the protein you're interested in and load a reference pdb and list of alpha carbons in that domain
catdomain <- 57:300
ref_pdb <- read.pdb(file.path(dirs_cf[1], 'solvated_system.pdb'))
ref_ca  <- atom.select(ref_pdb, 'protein', elety = 'CA', resno = catdomain)
head(ref_pdb$atom[ref_ca$atom, c('resno', 'resid', 'elety')])

# define a function that will take the dcd and pdb of each replicate then fit it to the same reference
load_group <- function(dir_vec, state_label, n_last = 2500) {
  xyz_list <- lapply(dir_vec, function(d) {
    dcd <- read.dcd(file.path(d, 'trajectory_centered.dcd'))
    pdb <- read.pdb(file.path(d, 'solvated_system.pdb'))
    ca  <- atom.select(pdb, 'protein', elety = 'CA', resno = catdomain)
    dcd_sub <- tail(dcd, n_last)
    xyz <- fit.xyz(fixed = ref_pdb$xyz, mobile = dcd_sub,
                   fixed.inds = ref_ca$xyz, mobile.inds = ca$xyz)
    xyz[, ca$xyz]
  })
  n_frames <- sapply(xyz_list, nrow)
  meta <- data.frame(
    state     = state_label,
    replicate = factor(rep(seq_along(xyz_list), times = n_frames)),
    frame     = unlist(lapply(n_frames, seq_len))
  )
  list(xyz = do.call(rbind, xyz_list), meta = meta)
}

# read dcd files in the directories supplied, only using the last 50 ns (for me, it was the last 2500 frames)
apo_data <- load_group(dirs_apo, 'Apo',  n_last = 2500)
holo_data  <- load_group(dirs_holo, 'Holo', n_last = 2500)

# pools all replicates for combined PCA so that holo and apo can be directly compared
xyz_combined <- rbind(apo_data$xyz, holo_data$xyz)

# add metadata which is whether a replicate is holo or apo
meta_all <- rbind(apo_data$meta, holo_data$meta)

# compute pca and load scores into a data frame
pc_apovsholo <- pca.xyz(xyz_combined)
scores_apovsholo <- cbind(data.frame(PC1 = pc_apovsholo$z[, 1], PC2 = pc_apovsholo$z[, 2], PC3 = pc_apovsholo$z[, 3]), meta_all)

# calculate variance accounted for by each principal component
var_exp_apovsholo <- round(100 * pc_apovsholo$L / sum(pc_apovsholo$L), 1)

# graphs the scree plot as well as the first few principal components against one another
# just to check if your PCA actually worked before committing to making and polishing the full figure
plot(pc_apovsholo)

# plot using ggplot2, coloring by apo vs. holo
ggplot(scores_apovsholo, aes(x = PC1, y = PC2, color = state)) + stat_density_2d(alpha = 0.8, linewidth = 0.5) + geom_point(alpha = 0.9, size = 1.2) + geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey70', linewidth = 0.4) + geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey70', linewidth = 0.4) + scale_color_manual(values = c('Apo'='grey20', 'Holo'='firebrick')) + labs(x = paste0('PC1 (', var_exp_apovsholo[1], '% variance)'), y = paste0('PC2 (', var_exp_apovsholo[2], '% variance)'), color = 'Condition') + theme_classic()

# export before you forget
ggsave('pca_cf_v1.eps', height=6, width=7, dpi=300, device=cairo_ps)
