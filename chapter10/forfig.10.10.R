library(vegan)
phyto_metadata <- readRDS("phyto_metadata.obj")
species_ryuko_data <- readRDS("phyto_ryuko_data.obj")
species_richness <- apply(species_ryuko_data > 0, 1, sum)
total_abundance <- apply(species_ryuko_data, 1, sum)

#target matrix
sample_ID <- c("S01", "S02","S03","S04","S05","S06","S07","S08","S09","S10","S11","S12","S13","S14", "S15")
rownames(species_ryuko_data) <- sample_ID

PCoA_ryuko_BC <- capscale(species_ryuko_data[c(1,5,7,8,13),] ~ 1, distance = "bray")
PCoA1_ryuko_BC <- summary(PCoA_ryuko_BC)$sites[,1]
PCoA2_ryuko_BC <- summary(PCoA_ryuko_BC)$sites[,2]
plot(
  PCoA2_ryuko_BC ~ PCoA1_ryuko_BC,
  cex = 3, pch = 1,
  xlab = "", ylab = "",
  asp = 1,
  main = "With Bray-Curtis"
)
text(PCoA1_ryuko_BC + 0.2, PCoA2_ryuko_BC, label = sample_ID[c(1,5,7,8,13)], cex = 1.2)

#explanatory matrix
env_ryuko <- data.frame(SR = species_richness, TA = total_abundance, temp = phyto_metadata$temp)
env_ryuko_E <- capscale(env_ryuko[c(1,5,7,8,13),] ~ 1, distance = "euclidean")
env1_ryuko_E <- summary(env_ryuko_E)$sites[,1]
env2_ryuko_E <- summary(env_ryuko_E)$sites[,2]
plot(
  env2_ryuko_E ~ env1_ryuko_E,
  cex = 3, pch = 2,
  xlab = "", ylab = "",
  asp = 1,
  main = "With Bray-Curtis"
)
text(env1_ryuko_E + 3.0, env2_ryuko_E, label = sample_ID[c(1,5,7,8,13)], cex = 1.2)

test_pro <- procrustes(X = env_ryuko_E, Y = PCoA_ryuko_BC, symmteric = FALSE, scores = "sites", choice = c(1:2))

plot(test_pro, k = 1, type = "text")
