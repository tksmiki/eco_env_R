###For Chapter 4: Sections 4.1 & 4.3 
phyto_metadata <- read.csv("metadata_phyto.csv", header = T) #load file list
levels(as.factor(phyto_metadata$YYMMDD))

####load each data###############
raw_data <- list() #prepare empty list
no_sample <- length(phyto_metadata$file)
for(i in 1:no_sample) {
  print(i)
  raw_data[[i]] <- read.csv(as.character(phyto_metadata$file[i]), header = T)
}


###For Chapter 5: Section 5.1#######
phyto_metadata <- read.csv("metadata_phyto_e1.csv", header = T) #load file list

####load each data###############
raw_data <- list() #prepare empty list
no_sample <- length(phyto_metadata$file)
for(i in 1:no_sample) {
  print(i)
  raw_data[[i]] <- read.csv(as.character(phyto_metadata$file[i]), header = T)
}

raw_data[[5]] <- raw_data[[5]][-37, -4]

######For Section 5.2.1############### 
#change colnames
for(i in 1:length(phyto_metadata$file)) {
  colnames(raw_data[[i]])[3] <- as.character(phyto_metadata$file[i])
}
raw_data[[6]]$species[16] <- "Microcystis viridis"

merged_data <- merge(raw_data[[1]][,c(-1)], raw_data[[2]][,c(-1)], all = T)
#looping merge
for(i in 3:length(phyto_metadata$file)) {
  merged_data <- merge(merged_data, raw_data[[i]][,c(-1)], all = T)
}

View(merged_data)

#convert all NA to zero
merged_data[is.na(merged_data)] <- 0
#copy dataframe to new dataframe
species_ryuko_data <- merged_data
#copy species name to row names
rownames(species_ryuko_data) <- species_ryuko_data$species
#remove the redundant info
species_ryuko_data <- species_ryuko_data[,-1]

#transpose
species_ryuko_data <- as.data.frame(t(species_ryuko_data))

View(species_ryuko_data)
#check class
class(species_ryuko_data[,2])

