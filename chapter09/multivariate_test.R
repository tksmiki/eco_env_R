###For chapter 09
####Preparation#####
library(vegan)
phyto_metadata <- readRDS("phyto_metadata.obj")
species_ryuko_data <- readRDS("phyto_ryuko_data.obj")
metadata_ecoplate <- readRDS("metadata_ecopl.obj")
summary_ecoplate <- readRDS("summary_ecopl.obj")
species_richness <- apply(species_ryuko_data > 0, 1, sum)
total_abundance <- apply(species_ryuko_data, 1, sum)

substrate_name <- c("Pyruvic-Acid-Methyl-Ester", "Tween-40", "Tween-80","alpha-Cyclodextrin", "Glycogen", "D-Cellobiose","alpha-D-Lactose", "beta-Methyl-D-Glucoside", "D-Xylose", "i-Erythritol", "D-Mannitol","N-Acetyl-D-Glucosamine", "D-Glucosaminic-Acid", "Glucose-1-Phosphate","alpha-Glycerol-Phosphate","D-Galactonic-Acid-gamma-Lactone", "D-Galacturonic-Acid", "2-Hydroxy-Benzoic-Acid", "4-Hydroxy-Benzoic-Acid", "gamma-Hydroxybutyric-Acid", "Itaconic-Acid", "alpha-Ketobutyric-Acid", "D-Malic-Acid", "L-Arginine", "L-Asparagine", "L-Phenylalanine", "L-Serine", "L-Threonine", "Glycyl-L-Glutamic-Acid", "Phenylethyl-amine", "Putrescine")

substrate_jpn <- c("Pyruvic-Acid-Methyl-Ester", "Tween-40", "Tween-80","alpha-Cyclodextrin", "Glycogen", "D-Cellobiose","alpha-D-Lactose", "beta-Methyl-D-Glucoside", "D-Xylose", "i-Erythritol", "D-Mannitol","N-Acetyl-D-Glucosamine", "D-Glucosaminic-Acid", "Glucose-1-Phosphate","alpha-Glycerol-Phosphate","D-Galactonic-Acid-gamma-Lactone", "D-Galacturonic-Acid", "2-Hydroxy-Benzoic-Acid", "4-Hydroxy-Benzoic-Acid", "gamma-Aminobutyric-acid", "Itaconic-Acid", "alpha-Ketobutyric-Acid", "D-Malic-Acid", "L-Arginine", "L-Asparagine", "L-Phenylalanine", "L-Serine", "L-Threonine", "Glycyl-L-Glutamic-Acid", "Phenylethyl-amine", "Putrescine")

####PCoA plot for ecoplate####
#With Bray-Curtis
PCoA_ecoplate_BC <- summary(capscale(summary_ecoplate ~ 1, distance = "bray"))
PCoA_ecoplate_BC

PCoA1_ecoplate_BC <- PCoA_ecoplate_BC$sites[,1]
PCoA2_ecoplate_BC <- PCoA_ecoplate_BC$sites[,2]
plot(
  PCoA2_ecoplate_BC ~ PCoA1_ecoplate_BC,
  cex = 3, pch = as.numeric(as.factor(metadata_ecoplate$treatment)),
  xlab = "PCoA1 (24.6 %)", ylab = "PCoA2 (12.4 %)",
  asp = 1,
  main = "Ecoplate With Bray-Curtis"
)

####Calculate F-values####
ecoplate_BC.d <- vegdist(summary_ecoplate, method = "bray") #Bray-Curtis
ecoplate_BC_F <- adonis(ecoplate_BC.d ~ metadata_ecoplate$treatment)$aov.tab$F.Model[1]
ecoplate_BC_F
####Example of permutation####
set.seed(1235) #fix the random seed
leng <- length(metadata_ecoplate$treatment)
for(i in 1:3) {
  perm_treatment <- sample(metadata_ecoplate$treatment, leng, replace = FALSE) #shuffling the index
  F_perm <- round(adonis(ecoplate_BC.d ~ perm_treatment)$aov.tab$F.Model[1],4)
  plot(
   PCoA2_ecoplate_BC ~ PCoA1_ecoplate_BC,
   cex = 3, pch = as.numeric(as.factor(perm_treatment)),
   xlab = "PCoA1 (24.6 %)", ylab = "PCoA2 (12.4 %)",
   asp = 1,
   main = paste("Permutation trial-", i, ": F value = ", F_perm, sep = "")
  )
}
####1000 permutations
test_for_Fperm <- adonis(ecoplate_BC.d ~ metadata_ecoplate$treatment, perm = 999)
hist(test_for_Fperm$f.perms, main = "Frequency of permutational F-values")
sum(test_for_Fperm$f.perms >= ecoplate_BC_F) # fraction with which permutational F is equal or greater than the observed F value. 
(sum(test_for_Fperm$f.perms >= ecoplate_BC_F) + 1) / (999 + 1)

####PERMANOVA for ecoplate data####
ecoplate_BC.d <- vegdist(summary_ecoplate, method = "bray") #Bray-Curtis
ecoplate_BC_permanova <- adonis(ecoplate_BC.d ~ metadata_ecoplate$treatment, perm = 999)
ecoplate_BC_permanova$aov.tab
####PERMDISP for ecoplate data####
ecoplate_BC.d <- vegdist(summary_ecoplate, method = "bray") #Bray-Curtis
ecoplate_BC_var <- betadisper(ecoplate_BC.d, metadata_ecoplate$treatment)
permutest(ecoplate_BC_var, perm = 999)
#Visualization
boxplot(
  ecoplate_BC_var$distances ~ metadata_ecoplate$treatment,
  outline = FALSE,
  col = "white"
)
stripchart(
  ecoplate_BC_var$distances ~ metadata_ecoplate$treatment,
  method = "stack",
  pch = c(1, 2),
  cex = 3,
  vertical = TRUE,
  add = TRUE
)

####PERMDISP for phytoplankton####
ryuko_BC.d <- vegdist(species_ryuko_data[c(1, 2, 13, 14, 15),], method = "bray") #Bray-Curtis
ryuko_BC_var <- betadisper(ryuko_BC.d, phyto_metadata$month[c(1, 2, 13, 14, 15)])
permutest(ryuko_BC_var)
#Visualization
boxplot(
  ryuko_BC_var$distances ~ phyto_metadata$month[c(1, 2, 13, 14, 15)],
  outline = FALSE,
  col = "white",
  xlab = "Month", ylab = "distance to center"
)
stripchart(
  ryuko_BC_var$distances ~ phyto_metadata$month[c(1, 2, 13, 14, 15)],
  method = "stack",
  pch = c(1, 2),
  cex = 3,
  vertical = TRUE,
  add = TRUE
)


####RDA for phytoplankton
phyto.rda <- rda(species_ryuko_data ~ phyto_metadata$temp + total_abundance + species_richness)
summary(phyto.rda)
plot(phyto.rda, scaling = 1, type = "text")
plot(phyto.rda, scaling = 2, type = "text")
permutest(phyto.rda, by = "terms", perm = 1999)

####CAP/dbRDA for phytoplankton#### 
ryuko_BC.d <- vegdist(species_ryuko_data, method = "bray") #Bray-Curtis
temperature <- phyto_metadata$temp

ryuko_CAP <- capscale(ryuko_BC.d ~ temperature + total_abundance + species_richness)
ryuko_dbRDA <- dbrda(ryuko_BC.d ~ temperature + total_abundance + species_richness)

summary(ryuko_CAP)
summary(ryuko_dbRDA)

plot(ryuko_CAP, scaling = 2, type = "text")
plot(ryuko_dbRDA, scaling = 2, type = "text")

permutest(ryuko_CAP, by = "terms", perm = 1999)
permutest(ryuko_dbRDA, by = "terms", perm = 1999)

month <- phyto_metadata$month
ryuko_CAP2 <- capscale(ryuko_BC.d ~ temperature + total_abundance + species_richness + month)
permutest(ryuko_CAP2, by = "terms", perm = 1999)

####PERMANOVA, CAP, and dbRDA####
ryuko_J.d <- vegdist(species_ryuko_data, method = "jaccard", binary = TRUE) #jaccard
adonis(ryuko_J.d  ~ month)$aov.tab
permutest(capscale(ryuko_J.d ~ month), by = "terms")
permutest(dbrda(ryuko_J.d ~ month), by = "terms")
