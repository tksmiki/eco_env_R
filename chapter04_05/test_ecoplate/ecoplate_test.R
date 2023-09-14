#Date, treatment, and file information from Xitou data sets
metadata_ecoplate <- read.csv("format_xitou_pattern.csv", header = T)

dat_list <- list()
no_sample <- length(metadata_ecoplate$data_file)
for(j in 1:no_sample){
  print(j)
  file_path <- paste("./text_file/", metadata_ecoplate$data_file[j], sep = "")
  dat_list[[j]] <- read.table(file_path, skip=5)
}

####Function to load ecoplate data
#Parameter list
#relative_path: relative path of data folder from the folder where the R script is saved
#file_list: the vector of file names that we intend to load 
#no_skip: the number of skip rows in ecoplate text file

load_ecoplate_data <- function(relative_path, file_list, no_skip = 5)    
{
  data_list <- list()   
  for(i in 1:length(file_list)) {
    file_name <- paste(relative_path, file_list[i], sep = "")
    data_list[[i]] <- read.table(file_name, skip = no_skip)  
  }
  return(data_list)  #output (return value) of this function
}

load_ecoplate_data(relative_path = "./text_file/", file_list = metadata_ecoplate$data_file, no_skip = 5)
data_ecoplate <- load_ecoplate_data(relative_path = "./text_file/", file_list = metadata_ecoplate$data_file, no_skip = 5)

load_ecoplate_data2 <-  function(relative_path, file_list, no_skip = 5)    
{
  data_list <-list() 
  for(i in 1:length(file_list)) {
    file_name <- paste(relative_path, file_list[i], sep="")
    e <- try(read.table(file_name, skip=no_skip), silent = FALSE)   #error management
    if(class(e) == "try-error") next  
    else data_list[[i]] <- read.table(file_name, skip=no_skip)
  }
  return(data_list)  #output (return value) of this function
}

metadata_ecoplate_e1 <- read.csv("format_xitou_pattern_e1.csv", header=T)
data_ecoplate2 <- load_ecoplate_data2(relative_path = "./text_file/", file_list = metadata_ecoplate_e1$data_file, no_skip = 5)

####Function to calculate the averages and standardization by control values
ave_ecoplate <- function(data_f){
  data_ave1<-(data_f$X1+data_f$X5+data_f$X9)/3.0  #take average
  data_ave2<-(data_f$X2+data_f$X6+data_f$X10)/3.0
  data_ave3<-(data_f$X3+data_f$X7+data_f$X11)/3.0
  data_ave4<-(data_f$X4+data_f$X8+data_f$X12)/3.0
  data_sum<-append(append(append(data_ave1,data_ave2),data_ave3),data_ave4)  
  data_sum_nor<-data_sum - data_sum[1] #normalizing by water well
  return(data_sum_nor) #output
}

stat_summary_ecoplate <- function(data_f, sample_name, variable_name)
{
  data_summary <- ave_ecoplate(data_f[[1]])
  for(i in 2:length(data_f)) {
    if(is.null(data_f[[i]])) next;  #error management, skipping the non-measured dates
      data_summary <- rbind.data.frame(data_summary, ave_ecoplate(data_f[[i]]))
    }#end of for i
  data_summary <- data_summary[,-1]
  data_summary[data_summary < 0] <- 0
  colnames(data_summary) <- variable_name
  rownames(data_summary) <- sample_name
  return(data_summary)
}

substrate_name <- c("s01", "s02","s03","s04","s05","s06","s07","s08","s09","s10","s11", "s12","s13","s14","s15","s16","s17","s18","s19","s20","s21", "s22","s23","s24","s25","s26","s27","s28","s29","s30", "S31")
summary_ecoplate <- stat_summary_ecoplate(data_f = data_ecoplate, sample_name = metadata_ecoplate$sample, variable_name = substrate_name)

