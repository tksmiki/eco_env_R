#From Chapter03 3.1.2
#simple generation of vectors
test_v <- c(1.0, 2.0, 2.5)
test2_v <- c("red", "yellow", "green")
class(test_v)
class(test2_v)
test3_v <- c(2.0, 3.0, 4.5)

test_v %*% test3_v

df1 <- data.frame(length = test_v, weight = test3_v)
rownames(df1) <- test2_v
df1
mat1 <- as.matrix(df1)
mat1
class(mat1)


test_list <- list()
test_list[[1]] <- df1
test_list[[2]] <- mat1
test_list[[3]] <- c(1,2,3,4)

test_list[[3]]
test_list[[1]]$weight

test_data01 <- read.csv("./data_sample2.csv", header = T)
test_data02 <- read.table("./20210810_04_day14.txt", skip = 13, header = T)
