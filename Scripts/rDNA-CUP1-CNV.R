library(zoo)



files <- list.files(pattern = "forLOH$")
for (file in files) {
  LOH1 <- read.delim(file, header = FALSE, quote = "", stringsAsFactors = FALSE)
  LOH <- subset(LOH1, LOH1[,3] != LOH1[,6], select = c(1, 2, 3, 4, 5, 6))
  n <- nrow(LOH)
  w303 <- c()
  yjm <- c()
  readepth <- mean(LOH[,4])
  
  for (i in 1:n) {
    w303[i] <- length(unlist(strsplit(LOH[i, 5], split = ",|\\.", fixed = FALSE))) - 1
    yjm[i] <- length(unlist(strsplit(toupper(LOH[i, 5]), split = LOH[i, 6], fixed = FALSE))) - 1
  }
  
  mapWJLOH <- cbind.data.frame(LOH[, 1], as.numeric(LOH[, 2]), w303, yjm)
  CUP1 <- subset(mapWJLOH, mapWJLOH[, 1] == "chrVIII" & (mapWJLOH[, 2] > 211729 & mapWJLOH[, 2] < 216235), select = c(1, 2, 3, 4))
  rDNA <- subset(mapWJLOH, mapWJLOH[, 1] == "chrXII" & (mapWJLOH[, 2] > 451575 & mapWJLOH[, 2] < 468812), select = c(1, 2, 3, 4))
  
  Wcup1 <- mean(CUP1[, 3])
  Ycup1 <- mean(CUP1[, 4])
  CUP1COV <- cbind(Wcup1, Ycup1)
  
  WrDNA <- mean(rDNA[, 3])
  YrDNA <- mean(rDNA[, 4])
  rDNACOV <- cbind(sample = sub("forLOH$", "", file), readepth, WrDNA, YrDNA)
  
  rDNACUP1 <- cbind(rDNACOV, CUP1COV)
  
  output_file <- paste0(sub("forLOH$", "", file), "CUP1-rDNA")
  
  write.table(rDNACUP1, file = output_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}
