###For chapter 06
phyto_metadata <- readRDS("phyto_metadata.obj")
species_ryuko_data <- readRDS("phyto_ryuko_data.obj")
metadata_ecoplate <- readRDS("metadata_ecopl.obj")
summary_ecoplate <- readRDS("summary_ecopl.obj")

species_richness <- c()
for(j in 1: length(species_ryuko_data[,1])) {
  species_richness[j] <- sum(species_ryuko_data[j,] > 0)
}

species_richness <- apply(species_ryuko_data > 0, 1, sum)

plot(species_richness ~ as.numeric(phyto_metadata$year),
     type = "p", 
     cex = 3
     )

plot(species_richness ~ as.numeric(phyto_metadata$year), type = "p", cex = 3)

boxplot(
  species_richness ~ phyto_metadata$year,
  outline = TRUE,
  col = "white",
  xlab = "year 2018 vs 2019",
  ylab = "phytoplankton richness"
)

boxplot(
  species_richness ~ phyto_metadata$year,
  outline = TRUE,
  col = "white"
)
stripchart(
  species_richness ~ phyto_metadata$year,
  method = "jitter",
  pch = c(1,2),
  cex = 3,
  vertical = TRUE,
  add = TRUE
)

plot(
  species_richness ~ phyto_metadata$temp, 
  type = "p", 
  cex = 3,
  xlab = "temperature",
  ylab = "species richness"
)

plot(summary_ecoplate[,1:5], cex = 2)

