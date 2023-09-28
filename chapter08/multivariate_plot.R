###For chapter 08
library(vegan)
phyto_metadata <- readRDS("phyto_metadata.obj")
species_ryuko_data <- readRDS("phyto_ryuko_data.obj")
metadata_ecoplate <- readRDS("metadata_ecopl.obj")
summary_ecoplate <- readRDS("summary_ecopl.obj")
species_richness <- apply(species_ryuko_data > 0, 1, sum)
total_abundance <- apply(species_ryuko_data, 1, sum)

substrate_name <- c("Pyruvic-Acid-Methyl-Ester", "Tween-40", "Tween-80","alpha-Cyclodextrin", "Glycogen", "D-Cellobiose","alpha-D-Lactose", "beta-Methyl-D-Glucoside", "D-Xylose", "i-Erythritol", "D-Mannitol","N-Acetyl-D-Glucosamine", "D-Glucosaminic-Acid", "Glucose-1-Phosphate","alpha-Glycerol-Phosphate","D-Galactonic-Acid-gamma-Lactone", "D-Galacturonic-Acid", "2-Hydroxy-Benzoic-Acid", "4-Hydroxy-Benzoic-Acid", "gamma-Hydroxybutyric-Acid", "Itaconic-Acid", "alpha-Ketobutyric-Acid", "D-Malic-Acid", "L-Arginine", "L-Asparagine", "L-Phenylalanine", "L-Serine", "L-Threonine", "Glycyl-L-Glutamic-Acid", "Phenylethyl-amine", "Putrescine")

substrate_jpn <- c("Pyruvic-Acid-Methyl-Ester", "Tween-40", "Tween-80","alpha-Cyclodextrin", "Glycogen", "D-Cellobiose","alpha-D-Lactose", "beta-Methyl-D-Glucoside", "D-Xylose", "i-Erythritol", "D-Mannitol","N-Acetyl-D-Glucosamine", "D-Glucosaminic-Acid", "Glucose-1-Phosphate","alpha-Glycerol-Phosphate","D-Galactonic-Acid-gamma-Lactone", "D-Galacturonic-Acid", "2-Hydroxy-Benzoic-Acid", "4-Hydroxy-Benzoic-Acid", "gamma-Aminobutyric-acid", "Itaconic-Acid", "alpha-Ketobutyric-Acid", "D-Malic-Acid", "L-Arginine", "L-Asparagine", "L-Phenylalanine", "L-Serine", "L-Threonine", "Glycyl-L-Glutamic-Acid", "Phenylethyl-amine", "Putrescine")

plot(
  s05 ~ s04, data = summary_ecoplate, 
  type = "p", 
  cex = 3,
  pch = c(1,5)[as.factor(metadata_ecoplate$treatment)],
  xlab = substrate_name[4],
  ylab = substrate_name[5],
  xlim = c(0, 4),
  ylim = c(0, 4),
  asp = 1.0
)
model01 <- lm(s05 ~ s04, data = summary_ecoplate)
abline(model01, lty = 5)
model02 <- summary(capscale(summary_ecoplate[,4:5] ~ 1, distance = "euclidean"))
slope_PC1 <- model02$species[2,1]/model02$species[1,1]
intercept_PC1 <- mean(summary_ecoplate$s05) - slope_PC1*mean(summary_ecoplate$s04)
abline(intercept_PC1, slope_PC1)
par(new = T)
plot(
  mean(summary_ecoplate$s04), mean(summary_ecoplate$s05),
  cex = 3,
  pch = 3,
  xlim = c(0, 4),
  ylim = c(0, 4),
  xlab = "", ylab = "",
  asp = 1.0
)

y_dummy <- rep(0, length(summary_ecoplate[,1]))
plot(
  y_dummy ~ model02$sites[, 1],
  pch = c(1,5)[as.factor(metadata_ecoplate$treatment)],
  cex = 3,
  xlim = c(-2,2), ylim = c(-0.5,0.5),
  xlab = "New axis",
  ylab = ""
)
par(new = T)
plot(
  0.0, 0.0,
  cex = 3,
  pch = 3,
  xlim = c(-2,2), ylim = c(-0.5,0.5),
  xlab = "", ylab = ""
)

####Plot 3D data####
library(rgl)
plot3d(
  x = summary_ecoplate$s04, y = summary_ecoplate$s05, z = summary_ecoplate$s06,
  xlab = substrate_name[4], ylab = substrate_name[5], zlab = substrate_name[6],
  type = "s",
  size = 3,
  col = c(2,9)[as.factor(metadata_ecoplate$treatment)]
)
model03 <- summary(capscale(summary_ecoplate[,4:6] ~ 1, distance = "euclidean"))
#PC1-PC2 plane can be defined by the orthogonal vector (i.e, PC3)
PC1_coeff <- model03$species[1:3,1]
PC2_coeff <- model03$species[1:3,2]
PC3_coeff <- model03$species[1:3,3]
new_origin <- c(mean(summary_ecoplate$s04), mean(summary_ecoplate$s05), mean(summary_ecoplate$s06))
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
  p1 = new_origin + PC2_coeff*0.4,
  s = 1/4,
  type = "rotation"
)
####Reduction to 2D####
plot(
  model03$sites[, 2] ~ model03$sites[, 1],
  pch = c(1,5)[as.factor(metadata_ecoplate$treatment)],
  cex = 3,
  xlim = c(-3,3), ylim = c(-3,3),
  xlab = "New axis 1",
  ylab = "New axis 2",
  asp = 1.0
)
####Reduction to 1D####
y_dummy <- rep(0, length(summary_ecoplate[,1]))
plot(
  y_dummy ~ model03$sites[, 1],
  pch = c(1,5)[as.factor(metadata_ecoplate$treatment)],
  cex = 3,
  xlim = c(-3,3), ylim = c(-0.5,0.5),
  xlab = "New axis 1",
  ylab = ""
)

####PCA for 3D data with 3 PC axes####
plot3d(
  x = summary_ecoplate$s04, y = summary_ecoplate$s05, z = summary_ecoplate$s06,
  xlab = substrate_name[4], ylab = substrate_name[5], zlab = substrate_name[6],
  type = "s",
  size = 3,
  col = c(2,9)[as.factor(metadata_ecoplate$treatment)]
)
model03 <- summary(capscale(summary_ecoplate[,4:6] ~ 1, distance = "euclidean"))
PC1_coeff <- model03$species[1:3,1]
PC2_coeff <- model03$species[1:3,2]
PC3_coeff <- model03$species[1:3,3]
new_origin <- c(mean(summary_ecoplate$s04), mean(summary_ecoplate$s05), mean(summary_ecoplate$s06))
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
  p1 = new_origin + PC2_coeff*0.4,
  s = 1/4,
  type = "rotation"
)
arrow3d(
  p0 = new_origin,
  p1 = new_origin + PC3_coeff*0.4,
  s = 1/4,
  type = "rotation"
)

####PCA for ecoplate#####
PCA_model01 <- summary(capscale(summary_ecoplate ~ 1, distance = "euclidean"))
PCA_model01
PCA_model01$sites

PC1_01 <- PCA_model01$sites[,1]
PC2_01 <- PCA_model01$sites[,2]
plot(
  PC2_01 ~ PC1_01,
  cex = 3, pch = c(1,5)[as.factor(metadata_ecoplate$treatment)],
  xlab = "PC1 (21.1 %)", ylab = "PC2 (13.8 %) "
)

####PCA for phytoplankton####
PCA_model02 <- summary(capscale(species_ryuko_data ~ 1, distance = "euclidean"))
PCA_model02
PCA_model02$sites

PC1_02 <- PCA_model02$sites[,1]
PC2_02 <- PCA_model02$sites[,2]
plot(
  PC2_02 ~ PC1_02,
  cex = 3, pch = as.numeric(as.factor(phyto_metadata$month)),
  xlab = "PC1 (62.0%)", ylab = "PC2 (28.3 %) "
)

####PCA for airquality data####
data("airquality")
summary(airquality)
air_data <- na.omit(airquality)

PCA_model03 <- summary(capscale(air_data[,1:4] ~ 1, distance = "euclidean"))
PCA_model03

PC1_03 <- PCA_model03$sites[,1]
PC2_03 <- PCA_model03$sites[,2]
plot(
  PC2_03 ~ PC1_03,
  cex = 0.5, pch = air_data$Month,
  xlab = "PC1 (89.0 %)", ylab = "PC2 (10.5 %) ",
  asp = 1
)
text(PC1_03 + 0.5, PC2_03, labels = rownames(air_data), cex = 0.8)

air_data2 <- scale(air_data)
summary(air_data2[,1:4])

PCA_model04 <- summary(capscale(air_data2[,1:4] ~ 1, distance = "euclidean"))
PCA_model04

PC1_04 <- PCA_model04$sites[,1]
PC2_04 <- PCA_model04$sites[,2]
plot(
  PC2_04 ~ PC1_04,
  cex = 0.5, pch = air_data$Month,
  xlab = "PC1 (59.0 %)", ylab = "PC2 (22.4 %) ",  
  asp = 1
)
text(PC1_04+0.05, PC2_04, labels = rownames(air_data), cex = 0.8)

#Large variations between variables with the same unit
summary(species_ryuko_data$`Fragilaria crotonensis`)
summary(species_ryuko_data$`Micrasterias hardyi`)

####Standard way of PCA####
PCA_s01 <- prcomp(air_data[,1:4], scale. = F)
PCA_s02 <- prcomp(air_data[,1:4], scale. = T)
summary(PCA_s01)
biplot(PCA_s01)
biplot(PCA_s02)



####Ecological Dissimilarity#####
sp1 <- c(3.0, 0.0, 0.1, 6.0)
sp2 <- c(3.5, 0.5, 0.0, 6.0)
comm <- data.frame(sp.1 = sp1, sp.2 = sp2)
rownames(comm) <- c("A", "B", "C", "D")

plot(
  sp.2 ~ sp.1, data = comm,
  cex = 3.0,
  xlim = c(0,7), ylim = c(0,7),
  asp = 1.0
)
text(comm$sp.1, comm$sp.2, labels = rownames(comm), cex = 0.8)

#examples of distance
vegdist(comm, method = "euclidean")  #Euclidean 
vegdist(comm, method = "jaccard", binary = TRUE) #Jaccard
vegdist(comm, method = "bray") #Bray-Curtis

####PCoA for comm####
PCoA_comm_BC <- summary(capscale(comm ~ 1, distance = "bray"))
PCoA_comm_BC
PCoA1_comm_BC <- PCoA_comm_BC$sites[,1]
PCoA2_comm_BC <- PCoA_comm_BC$sites[,2]
plot(
  PCoA2_comm_BC ~ PCoA1_comm_BC,
  cex = 3,
  xlab = "PCoA1 (54.5 %)", ylab = "PCoA2 (41.8 %) ",
  asp = 1
)
text(PCoA1_comm_BC, PCoA2_comm_BC, labels = rownames(comm), cex = 0.8)


####PCoA for phytoplankton####
#With Jaccard
species_ryuko_data_b <- species_ryuko_data #copy
species_ryuko_data_b[species_ryuko_data_b > 0] <- 1 #binalization
PCoA_ryuko_J <- summary(capscale(species_ryuko_data_b ~ 1, distance = "jaccard"))
PCoA_ryuko_J

PCoA1_ryuko_J <- PCoA_ryuko_J$sites[,1]
PCoA2_ryuko_J <- PCoA_ryuko_J$sites[,2]
plot(
  PCoA2_ryuko_J ~ PCoA1_ryuko_J,
  cex = 3, pch = as.numeric(as.factor(phyto_metadata$month)),
  xlab = "PCoA1 (19.7 %)", ylab = "PCoA2 (12.3 %) ",
  asp = 1,
  main = "With Jaccard"
)

#With Bray-Curtis
PCoA_ryuko_BC <- summary(capscale(species_ryuko_data ~ 1, distance = "bray"))
PCoA_ryuko_BC

PCoA1_ryuko_BC <- PCoA_ryuko_BC$sites[,1]
PCoA2_ryuko_BC <- PCoA_ryuko_BC$sites[,2]
plot(
  PCoA2_ryuko_BC ~ PCoA1_ryuko_BC,
  cex = 3, pch = as.numeric(as.factor(phyto_metadata$month)),
  xlab = "PCoA1 (25.4 %)", ylab = "PCoA2 (14.7 %) ",
  asp = 1,
  main = "With Bray-Curtis"
)

#Hellinger
species_ryuko_data_H <- decostand(species_ryuko_data, method = "hellinger")
PCoA_ryuko_H <- summary(capscale(species_ryuko_data_H ~ 1, distance = "euclidean"))
PCoA_ryuko_H

PCoA1_ryuko_H <- PCoA_ryuko_H$sites[,1]
PCoA2_ryuko_H <- PCoA_ryuko_H$sites[,2]
plot(
  PCoA2_ryuko_H ~ PCoA1_ryuko_H,
  cex = 3, pch = as.numeric(as.factor(phyto_metadata$month)),
  xlab = "PCoA1 (22.4 %)", ylab = "PCoA2 (14. %) ",
  asp = 1,
  main = "With Hellinger"
)

####PCoA for ecoplate####
#With Jaccard
summary_ecoplate_b <- summary_ecoplate #copy
minimum_strength <- 0.2
summary_ecoplate_b[summary_ecoplate_b < minimum_strength] <- 0 #binalization
summary_ecoplate_b[summary_ecoplate_b > 0] <- 1 #binalization
PCoA_ecoplate_J <- summary(capscale(summary_ecoplate_b ~ 1, distance = "jaccard"))
PCoA_ecoplate_J 

PCoA1_ecoplate_J <- PCoA_ecoplate_J$sites[,1]
PCoA2_ecoplate_J<- PCoA_ecoplate_J$sites[,2]
plot(
  PCoA2_ecoplate_J ~ PCoA1_ecoplate_J ,
  cex = 3, pch = c(1,5)[as.factor(metadata_ecoplate$treatment)],
  xlab = "PCoA1 (51.9 %)", ylab = "PCoA2 (30.2 %) ",
  asp = 1,
  main = "With Jaccard"
)

####BOX10 Inconsistent behavior of visualization####
summary(prcomp(comm, scale. = T))
summary(prcomp(comm, scale. = F))

biplot(prcomp(comm, scale. = T))
biplot(prcomp(comm, scale. = F))
ordiplot(capscale(comm ~ 1, distance = "euclidean"), type = "text")
ordiplot(capscale(scale(comm) ~ 1, distance = "euclidean"), type = "text")
ordiplot(prcomp(comm, scale. = T), type = "text")
ordiplot(prcomp(comm, scale. = F), type="text")
ordiplot(prcomp(scale(comm), scale. = F), type = "text")

###Non-hierarchical clustering####
#as hypothetical example
plot(
  PCoA2_ryuko_BC ~ PCoA1_ryuko_BC,
  cex = 3, 
  asp = 1,
)

####Hierarchical clustering####
species_b.d <- vegdist(species_ryuko_data, method = "bray")
hclust_model <- hclust(species_b.d, method = "ward.D2")
plot(
  hclust_model,
  hang = -1,
  main = "phytoplankton composition with Bray-Curtis",
  label = phyto_metadata$YYMMDD
)

plot(
  PCoA2_ryuko_BC ~ PCoA1_ryuko_BC,
  cex = 2, pch = as.numeric(as.factor(phyto_metadata$month)),
  xlab = "PCoA1 (25.4 %)", ylab = "PCoA2 (14.7 %) ",
  asp = 1,
  main = "With Bray-Curtis"
)
text(PCoA1_ryuko_BC + 0.25, PCoA2_ryuko_BC, label = phyto_metadata$YYMMDD, cex = 0.8)
