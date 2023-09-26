###For chapter 07
phyto_metadata <- readRDS("phyto_metadata.obj")
species_ryuko_data <- readRDS("phyto_ryuko_data.obj")
metadata_ecoplate <- readRDS("metadata_ecopl.obj")
summary_ecoplate <- readRDS("summary_ecopl.obj")
species_richness <- apply(species_ryuko_data > 0, 1, sum)
total_abundance <- apply(species_ryuko_data, 1, sum)

###multiple regression
plot(species_richness ~ phyto_metadata$temp, cex = 2)
plot(species_richness ~ total_abundance, cex = 2)

library(scatterplot3d)
scatterplot3d(
  x = phyto_metadata$temp, y = total_abundance, z = species_richness,
  xlab = "temperature", ylab = "abundance", zlab = "richness",
  angle = 30
)

library(rgl)
plot3d(
  x = phyto_metadata$temp, y = total_abundance, z = species_richness,
  xlab = "temperature", ylab = "abundance", zlab = "richness",
  type = "s",
  size = 3,
  col = c(2,9)[as.factor(phyto_metadata$year)]
)

model_temp <- lm(species_richness ~ phyto_metadata$temp)
model_abundance <- lm(species_richness ~ total_abundance)
plot(species_richness ~ phyto_metadata$temp, cex = 2)
abline(model_temp)
plot(species_richness ~ total_abundance, cex = 2)
abline(model_abundance)

plot3d(
  x = phyto_metadata$temp, y = total_abundance, z = species_richness,
  xlab = "temperature", ylab = "abundance", zlab = "richness",
  type = "s",
  size = 3,
  col = c(2,9)[as.factor(phyto_metadata$year)]
)
fit_model <- lm(species_richness ~ phyto_metadata$temp + total_abundance)
plane_coeff <- coef(fit_model)
planes3d(plane_coeff[2], plane_coeff[3], -1, plane_coeff[1], col = "blue", alpha = 0.5)

multi_model <- lm(species_richness ~ phyto_metadata$temp + total_abundance)
summary(multi_model)


#interaction
#ref https://momonoki2017.blogspot.com/2018/08/r019.html
data("ToothGrowth")
sub_TG <- subset(ToothGrowth, dose > 0.5)
boxplot(
  len ~ dose + supp, data = sub_TG, 
  outline = F,
  col = c("white", "white", "gray", "gray")
)
stripchart(
  len ~ dose + supp, data = sub_TG,
  method = "stack",
  pch = c(1,1,2,2),
  cex = 3,
  vertical = TRUE,
  add = TRUE
)
interaction.plot(
  x.factor = sub_TG$dose,
  response = sub_TG$len,
  trace.factor = sub_TG$supp,
  xlab = "dose level", ylab = "growth"
)

model_tooth01 <- lm(len ~ dose + supp, data = sub_TG)
anova(model_tooth01)

model_tooth02 <- lm(len ~ dose + supp + dose:supp, data = sub_TG)
anova(model_tooth02)

####General Linear Model01####
summary(lm(species_richness ~ phyto_metadata$month))

####General Linear Model02####
data("ToothGrowth")
class(ToothGrowth$supp)
class(ToothGrowth$dose)

interaction.plot(
  x.factor = ToothGrowth$dose,
  response = ToothGrowth$len,
  trace.factor = ToothGrowth$supp,
  xlab = "dose level", ylab = "growth"
)

model_tooth03 <- lm(len ~ dose + supp + dose:supp, data = ToothGrowth)
summary(model_tooth03)

####General Linear Model03####
summary(lm(species_richness ~ as.factor(phyto_metadata$year)))


####Model Selection####
#2018 vs 2019
boxplot(
  species_richness ~ phyto_metadata$year, 
  outline = F,
  col = c("white", "gray")
)
stripchart(
  species_richness ~ phyto_metadata$year,
  method = "stack",
  pch = c(1,2),
  cex = 3,
  vertical = TRUE,
  add = TRUE
)
model_temp <- lm(species_richness ~ phyto_metadata$temp)
model_abundance <- lm(species_richness ~ total_abundance)
plot(species_richness ~ phyto_metadata$temp, cex = 2)
abline(model_temp)
plot(species_richness ~ total_abundance, cex = 2)
abline(model_abundance)

selected_model <- step(lm(species_richness ~ phyto_metadata$temp + total_abundance + phyto_metadata$year))
summary(selected_model)



