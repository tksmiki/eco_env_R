library(vegan)
library(ggplot2)
phyto_metadata <- readRDS("phyto_metadata.obj")
species_ryuko_data <- readRDS("phyto_ryuko_data.obj")
species_richness <- apply(species_ryuko_data > 0, 1, sum)
total_abundance <- apply(species_ryuko_data, 1, sum)

####Example 01: 2D scatter plot and regression  line####
#Classical plot (base)
model01 <- lm(species_richness ~ phyto_metadata$temp)
plot(
  species_richness ~ phyto_metadata$temp, 
  type = "p", 
  cex = 3,
  xlab = "temperature",
  ylab = "species richness"
)
abline(model01,col = 4)

#ggplot
#Combine metadata and numeric data into a single dataframe
lm_forggplot2 <- data.frame(SR = species_richness, TA = total_abundance, phyto_metadata)
#simple code
ggplot(lm_forggplot2, aes(x = temp, y = SR)) + 
  geom_point(size = 5, aes(shape = as.factor(year)),alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, colour = "blue", size = 1) + 
  xlab("temperature") +
  ylab("species richness")

#alternative way
lm_ggplot <- ggplot(lm_forggplot2, aes(x = temp, y = SR))
lm_ggplot <- lm_ggplot + geom_point(size = 5, aes(shape = as.factor(year)),alpha = 0.5)
lm_ggplot <- lm_ggplot + geom_smooth(method = "lm", se = FALSE, colour = "blue", size = 1) 
lm_ggplot <- lm_ggplot + xlab("temperature") + ylab("species richness")
print(lm_ggplot)


####Exmple 02: PCoA plot####
#Classical plot (base)
PCoA_ryuko_BC <- summary(capscale(species_ryuko_data ~ 1, distance = "bray"))
PCoA1_ryuko_BC <- PCoA_ryuko_BC$sites[,1]
PCoA2_ryuko_BC <- PCoA_ryuko_BC$sites[,2]
plot(
  PCoA2_ryuko_BC ~ PCoA1_ryuko_BC,
  cex = 3, pch = as.numeric(as.factor(phyto_metadata$month)),
  xlab = "PCoA1 (25.4 %)", ylab = "PCoA2 (14.7 %) ",
  asp = 1,
  main = "Phyto With Bray-Curtis"
)

#ggplot
#Combine metadata and numeric data of PCoA axes into a single dataframe
ryuko_forggplot2 <- data.frame(lm_forggplot2, PCoA1 = PCoA1_ryuko_BC, PCoA2 = PCoA2_ryuko_BC)

#ggplot2 for 2D scatter plot
g_PCoA <- ggplot(ryuko_forggplot2, aes(x = PCoA1, y = PCoA2))
g_PCoA <- g_PCoA + geom_point(aes(colour = month), size = 8, alpha = 0.5)
g_PCoA <- g_PCoA + labs(x = "PCoA (25.4 %)", y = "PCoA2 (14.7 %)") + ggtitle("Phyto With Bray-Curtis")
print(g_PCoA)

#Color should be carefully selected
g_PCoA2 <- ggplot(ryuko_forggplot2, aes(x = PCoA1,y = PCoA2, color = month))
g_PCoA2 <- g_PCoA2 + scale_colour_manual(values = c("#ff4b00", "#4dc4ff", "#f6aa00", "#804000"))
g_PCoA2 <- g_PCoA2 + geom_point(size = 8, aes(shape = month), alpha = 1)
g_PCoA2 <- g_PCoA2 + labs(x = "PCoA (25.4 %)", y = "PCoA2 (14.7 %)") + ggtitle("Phyto With Bray-Curtis")
g_PCoA2 <- g_PCoA2 + theme(panel.background = element_rect(fill = "transparent", colour = "black"),panel.grid = element_blank())
print(g_PCoA2)


