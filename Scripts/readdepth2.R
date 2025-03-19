library(zoo)


chrom.length <- read.delim("/home/dqlab/chromlength.txt", header = FALSE, sep = '', stringsAsFactors = FALSE)


txt_files <- list.files(pattern = "_500readdepth\\.txt$")
w29 <- read.delim('/home/dqlab/JSC25.txt', header = FALSE, quote = "", stringsAsFactors = FALSE)
for (txt_file in txt_files) {
  sample_name <- sub("_500readdepth\\.txt$", "", txt_file)

  LOH1 <- read.delim(txt_file, header = FALSE, quote = "", stringsAsFactors = FALSE)
  
  min_rows <- min(nrow(LOH1), nrow(w29))
  LOH1 <- LOH1[1:min_rows, ]
  w29 <- w29[1:min_rows, ]
  LOH1 <- cbind(LOH1, w29[, 4])
  readepth <- mean(LOH1[, 4], na.rm = TRUE)
  readepthw29 <- mean(LOH1[, 5], na.rm = TRUE)
  
  pdf(paste(sample_name, '_500readdepth.pdf', sep = ""), width = 11, height = 8)
  par(mfrow = c(4, 1), mar = c(2, 2, 2, 2))
  n_axis <- 10
  
  for (a in 1:16) {
    chrom_plot <- subset(LOH1,
                         LOH1[, 1] == paste('chr', as.roman(a), sep = ''),
                         select = c(2, 3, 4, 5))
    if (nrow(chrom_plot) == 0) {
      next
    }
    X <- chrom_plot[, 1]
    
    plot(X, (chrom_plot[, 3] / chrom_plot[, 4]) / (readepth / readepthw29),
         type = 'p',
         ylim = c(-1, 2),
         col = 'red',
         main = chrom.length[a, 1])
    
    abline(h = seq(-1, 2, by = 0.5), col = "lightgray", lty = 2)
    axis(1, pretty(range(X), n_axis))
    box()
  }
  
  dev.off()
}
