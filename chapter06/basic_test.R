###For chapter 06
phyto_metadata <- readRDS("phyto_metadata.obj")
species_ryuko_data <- readRDS("phyto_ryuko_data.obj")
metadata_ecoplate <- readRDS("metadata_ecopl.obj")
summary_ecoplate <- readRDS("summary_ecopl.obj")
#from basic_graphics.R, calculating species richness
species_richness <- apply(species_ryuko_data > 0, 1, sum)

subset(species_richness, phyto_metadata$year == "2018")
subset(species_ryuko_data, phyto_metadata$year == "2018")

species_richness2018 <- subset(species_richness, phyto_metadata$year == "2018")
species_richness2019 <- subset(species_richness, phyto_metadata$year == "2019")

summary(species_richness2018)
summary(species_richness2019)

tapply(species_richness, phyto_metadata$year, summary)
tapply(species_richness, phyto_metadata$year, mean)
tapply(species_richness, phyto_metadata$year, var)

t.test(species_richness ~ phyto_metadata$year, var.equal = F)

total_abundance <- apply(species_ryuko_data,1,sum)
boxplot(
  total_abundance ~ phyto_metadata$year,
  outline = TRUE,
  col = "white"
)
stripchart(
  total_abundance ~ phyto_metadata$year,
  method = "stack",
  pch = c(1,2),
  cex = 3,
  vertical = TRUE,
  add = TRUE
)

tapply(total_abundance, phyto_metadata$year, mean)
t.test(total_abundance ~ phyto_metadata$year, var.equal = F)

boxplot(species_richness ~ as.factor(phyto_metadata$month), 
     outline = FALSE,
     col = "white",
     xlab = "month",
     ylab = "species richness",
     ylim = c(0,40)
)
stripchart(
  species_richness ~ as.factor(phyto_metadata$month),
  method = "stack",
  cex = 2,
  vertical = TRUE,
  add = TRUE
)

anova(lm(species_richness ~ as.factor(phyto_metadata$month)))

data("InsectSprays")
View(InsectSprays)
boxplot(count ~ spray, data = InsectSprays,
        outline = FALSE,
        col = "white",
        xlab = "spray type",
        ylab = "insect count"
)
stripchart(
  count ~ as.factor(spray), data = InsectSprays,
  method = "stack",
  cex = 2,
  vertical = TRUE,
  add = TRUE
)

anova(lm(count ~ spray, data = InsectSprays))

plot(
  species_richness~phyto_metadata$temp, 
  type="p", 
  cex = 3,
  xlab = "temperature",
  ylab = "species richness"
)

cor(species_richness, phyto_metadata$temp)
cor(summary_ecoplate[,1:5])
cor(species_richness, phyto_metadata$temp, method = "spearman")

cor.test(species_richness, phyto_metadata$temp)

model01 <- lm(species_richness ~ phyto_metadata$temp)
summary(model01)

plot(
  species_richness ~ phyto_metadata$temp, 
  type = "p", 
  cex = 3,
  xlab = "temperature",
  ylab = "species richness"
)
abline(model01,col = 4)


cor(phyto_metadata$temp, species_richness)
cor(phyto_metadata$temp*2, species_richness)
cor(phyto_metadata$temp, species_richness*10)

temp2 <- 2*phyto_metadata$temp
species_richness10 <- 10*species_richness
lm(species_richness ~ phyto_metadata$temp)
lm(species_richness ~ temp2)
lm(species_richness10 ~ phyto_metadata$temp)

cor(phyto_metadata$temp, species_richness)
cor(species_richness, phyto_metadata$temp)

model02 <- lm(phyto_metadata$temp ~ species_richness)
model02$coefficients
# phyto_metadata$temp = 5.9579 + 0.4186*species_richness
# species_richness = -0.9579/0.4186 + (1.0/0.4186)*phyto_metadata$temp
plot(
  species_richness ~ phyto_metadata$temp, 
  type="p", 
  cex = 3,
  xlab = "temperature",
  ylab = "species richness"
)
abline(model01,col = 4)
abline(c(-model02$coefficients[1]/model02$coefficients[2], 1.0/model02$coefficients[2]),lty=2)

data(anscombe)
model_a01 <- lm(y1~x1, data = anscombe)
model_a02 <- lm(y2~x2, data = anscombe)
model_a03 <- lm(y3~x3, data = anscombe)
model_a04 <- lm(y4~x4, data = anscombe)
summary(model_a01)
summary(model_a02)
summary(model_a03)
summary(model_a04)
plot(
  y1~x1, data = anscombe, 
  type = "p", 
  cex = 3,
  xlim = c(0, 20),
  ylim = c(0, 15),
  xlab = "",
  ylab = ""
)
abline(model_a01)
plot(
  y2~x2, data = anscombe, 
  type = "p", 
  cex = 3,
  xlim = c(0, 20),
  ylim = c(0, 15),
  xlab = "",
  ylab = ""
)
abline(model_a02)
plot(
  y3~x3, data = anscombe, 
  type = "p", 
  cex = 3,
  xlim = c(0, 20),
  ylim = c(0, 15),
  xlab = "",
  ylab = ""
)
abline(model_a03)
plot(
  y4~x4, data = anscombe, 
  type = "p", 
  cex = 3,
  xlim = c(0, 20),
  ylim = c(0, 15),
  xlab = "",
  ylab = ""
)
abline(model_a04)

data_s <- read.csv("DatasaurusDozen.csv", header=T)
data_s01 <- subset(data_s, dataset == "circle")[2:3]
data_s02 <- subset(data_s, dataset == "dino")[2:3]
cor(data_s01)
cor(data_s02)
