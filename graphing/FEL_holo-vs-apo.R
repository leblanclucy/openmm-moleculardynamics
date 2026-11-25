# this computes the free energy landscape for two conditions, holoenzyme and apoenzyme
# this script assumes you already did PCA using the PCA_holo-vs-apo.R script
# you need to do PCA before FEL

library(MASS)

# define variables based on temperature used in the simulation plus the Boltzmann constant given units of kcal/(mol K)
T_kelvin <- 300
kB <- 0.001987204
kBT <- kB * T_kelvin

# define shared grid boundary applied to both apo and holo
pad_x <- diff(range(scores_apovsholo$PC1)) * 0.1
pad_y <- diff(range(scores_apovsholo$PC2)) * 0.1
shared_lims <- c(
  min(scores_apovsholo$PC1) - pad_x, max(scores_apovsholo$PC1) + pad_x,
  min(scores_apovsholo$PC2) - pad_y, max(scores_apovsholo$PC2) + pad_y
)

# function to calculate delta G
# set +5 kcal/mol as the upper limit as beyond that it's very rare, highly unfavorable conformations
calc_fel <- function(sub_df, state_name, max_dG = 5.0) {
  kde <- kde2d(x = sub_df$PC1, y = sub_df$PC2, n = 300, lims = shared_lims)
  grid <- expand.grid(PC1 = kde$x, PC2 = kde$y)
  grid$density <- as.vector(kde$z)
  grid$dG <- -kBT * log(grid$density / max(grid$density))
  grid$dG_capped <- ifelse(grid$dG > max_dG, NA, grid$dG)
  grid$state <- state_name
  grid
}

# calculate delta G for apo and holo
fel_apo  <- calc_fel(scores_apovsholo[scores_apovsholo$state == "Apo", ],  "Apo",  max_dG = 5.0)
fel_holo <- calc_fel(scores_apovsholo[scores_apovsholo$state == "Holo", ], "Holo", max_dG = 5.0)
fel_all  <- rbind(fel_apo, fel_holo)

# total variance accounted for by each principal component
var_exp_apovsholo <- round(100 * pc_apovsholo$L / sum(pc_apovsholo$L), 1)

# plot free energy landscape in 2D
# conformational basin representing the global minimum = very dark blue, more and more unfavorable = closer to green, yellow, then red
ggplot(fel_all, aes(x = PC1, y = PC2, z = dG_capped)) + geom_raster(aes(fill = dG_capped), interpolate = TRUE) + geom_contour(color = "white", alpha = 0.45, linewidth = 0.35, bins = 10) + facet_wrap(~ state) + scale_fill_viridis_c(option = "turbo", direction = 1, na.value  = "white", name = expression(Delta*G~"(kcal/mol)")) + labs(x = paste0("PC1 (", var_exp_apovsholo[1], "% variance)"), y = paste0("PC2 (", var_exp_apovsholo[2], "% variance)")) + theme_classic()

# export
ggsave('fel_2d_cf_v1.eps', height=6, width=13, dpi=300, device=cairo_ps)

# plot free energy landscape in 3D using plot3D which is easier than using plotly if you want to compare apo vs holo simultaneously
library(plot3D)

# first define a function that calculates delta G for holo and apo so it can be represented as a third axis
# then graph using persp3D, setting the same angles and axis limits so that holo and apo are viewed from the same perspective
plot_state_3d <- function(state_name, panel_title, show_colkey = FALSE) {
  sub_df <- scores_apovsholo[scores_apovsholo$state == state_name, ]
  kde <- kde2d(x = sub_df$PC1, y = sub_df$PC2, n = 200, lims = shared_lims)
  
  dG_mat <- -kBT * log(kde$z / max(kde$z))
  dG_mat[dG_mat > 5.0] <- NA
  
  persp3D(
    x = kde$x,
    y = kde$y,
    z = dG_mat,
    col = viridisLite::turbo(100),
    clim = c(0, 5),
    xlim = shared_lims[1:2],
    ylim = shared_lims[3:4],
    zlim = c(0, 5),
    theta = 40,
    phi = 30,
    expand = 0.6,
    border = NA,
    shade = 0.25,
    ticktype = "detailed",
    colkey = if (show_colkey) list(length = 0.5, width = 0.8) else FALSE,
    xlab = paste0("PC1 (", var_exp_apovsholo[1], "%)"),
    ylab = paste0("PC2 (", var_exp_apovsholo[2], "%)"),
    zlab = "Delta G (kcal/mol)",
    clab = c("Delta G", "(kcal/mol)"),
    main = panel_title
  )
}

# unlike ggplot it won't show in rstudio so you can export as png to quickly take a look
png("fel_3d_cf_v1.png", width = 3600, height = 1800, res = 300)
# par sets graphical parameters for plotting. mfrow divides the canvas into a 1 row by 2 column grid so the plots show side by side
# mar shrinks the margins around each panel. the defaults are bigger, so this makes it have less white space
par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))
# calls the function defined earlier, calculating delta G for apo, naming the plot apo, and then not showing the key on apo since you only need it once per panel
plot_state_3d("Apo",  "Apo",  show_colkey = FALSE)
# calls the function for the holo data
plot_state_3d("Holo", "Holo", show_colkey = TRUE)
# executes the graphing instructions. must do this or else the plot won't export at all
dev.off()

# and exporting for publication you need a vector image, so export as eps. same parameters as above for png
cairo_ps("fel_3d_cf_v1.eps", width = 12, height = 6, fallback_resolution = 600)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))
plot_state_3d("Apo",  "Apo",  show_colkey = FALSE)
plot_state_3d("Holo", "Holo", show_colkey = TRUE)
dev.off()
