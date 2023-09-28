###For Figure 9.9 in Chapter09
library(vagan)
library(rgl)
phyto_metadata <- readRDS("phyto_metadata.obj")
species_ryuko_data <- readRDS("phyto_ryuko_data.obj")
metadata_ecoplate <- readRDS("metadata_ecopl.obj")
summary_ecoplate <- readRDS("summary_ecopl.obj")
species_richness <- apply(species_ryuko_data > 0, 1, sum)
total_abundance <- apply(species_ryuko_data, 1, sum)

####multiple regression####
#one case
plot3d(
  x = phyto_metadata$temp, y = total_abundance, z = species_ryuko_data$`Anabaena affinis`,
  xlab = "", ylab = "", zlab = "",
  type = "s",
  size = 3,
  col = 9
)
fit_model <- lm(species_ryuko_data$`Anabaena affinis` ~ phyto_metadata$temp + total_abundance)
plane_coeff <- coef(fit_model)
planes3d(plane_coeff[2], plane_coeff[3], -1, plane_coeff[1], col = "blue", alpha = 0.5)

#another case
plot3d(
  x = phyto_metadata$temp, y = total_abundance, z = species_ryuko_data$`Microcystis aeruginosa`,
  xlab = "", ylab = "", zlab = "",
  type = "s",
  size = 3,
  col = 2
)
fit_model <- lm(species_ryuko_data$`Microcystis aeruginosa` ~ phyto_metadata$temp + total_abundance)
plane_coeff <- coef(fit_model)
planes3d(plane_coeff[2], plane_coeff[3], -1, plane_coeff[1], col = "blue", alpha = 0.5)


####PCA, PCoA####
plot3d(
  x = summary_ecoplate$s04, y = summary_ecoplate$s05, z = summary_ecoplate$s06,
  xlab = "", ylab = "", zlab = "",
  type = "s",
  size = 3,
  col = 9
)
model03 <- summary(capscale(summary_ecoplate[,4:6] ~ 1, distance="euclidean"))
#PC1-PC2 plane can be defined by the orthogonal vector (i.e, PC3)
PC1_coeff <- model03$species[1:3,1]
PC2_coeff <- model03$species[1:3,2]
PC3_coeff <- model03$species[1:3,3]
new_origin <- c(mean(summary_ecoplate$s04), mean(summary_ecoplate$s05), mean(summary_ecoplate$s06))
intercept_coeff <- (PC3_coeff %*% new_origin)*-1.0
planes3d(PC3_coeff[1], PC3_coeff[2], PC3_coeff[3], intercept_coeff, col = "blue", alpha = 0.5)
arrow3d(
  p0 = new_origin,
  p1 = new_origin + PC1_coeff*0.6,
  s = 1/4,
  type = "rotation"
)
arrow3d(
  p0 = new_origin,
  p1 = new_origin + PC2_coeff*0.6,
  s = 1/4,
  type = "rotation"
)
####Reduction to 2D####
plot(
  model03$sites[, 2] ~ model03$sites[, 1],
  pch = 1,
  cex = 3,
  xlim = c(-3,3), ylim = c(-3,3),
  xlab = "RDA1",
  ylab = "RDA2",
  asp = 1.0
)

plot3d(
  x = summary_ecoplate$s07, y = summary_ecoplate$s08, z = summary_ecoplate$s09,
  xlab = "", ylab = "", zlab = "",
  type = "s",
  size = 3,
  col = 9
)
model04 <- summary(capscale(summary_ecoplate[,7:9] ~ 1, distance="euclidean"))
#PC1-PC2 plane can be defined by the orthogonal vector (i.e, PC3)
PC1_coeff <- model03$species[1:3,1]
PC2_coeff <- model03$species[1:3,2]
PC3_coeff <- model03$species[1:3,3]
new_origin <- c(mean(summary_ecoplate$s07), mean(summary_ecoplate$s08), mean(summary_ecoplate$s09))
intercept_coeff <- (PC3_coeff %*% new_origin)*-1.0
planes3d(PC3_coeff[1], PC3_coeff[2], PC3_coeff[3], intercept_coeff, col = "blue", alpha = 0.5)
arrow3d(
  p0 = new_origin,
  p1 = new_origin + PC1_coeff*0.4,
  s = 1/4,
  type = "rotation"
)
arrow3d(
  p0 = new_origin,
  p1 = new_origin + PC2_coeff*0.6,
  s = 1/4,
  type = "rotation"
)
arrow3d(
  p0 = new_origin,
  p1 = new_origin + PC3_coeff*0.6,
  s = 1/4,
  type = "rotation"
)
