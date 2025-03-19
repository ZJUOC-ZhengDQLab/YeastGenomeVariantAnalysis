library(zoo)

# Set the working directory
setwd(".")

# Get files ending with 'forLOH' in the current directory
files <- list.files(pattern = "forLOH$")

for (file in files) {
    tryCatch({
        # Read the file
        LOH1 <- read.delim(file, header = FALSE, quote = "", stringsAsFactors = FALSE)
        print(paste("Reading data from file:", file))

        # Remove the suffix from the filename
        base_name <- sub("forLOH$", "", file)

        # Process LOH1 data...
        LOH <- subset(LOH1, LOH1[, 3] != LOH1[, 6], select = c(1, 2, 3, 4, 5, 6))
        n <- nrow(LOH)
        print(paste("Number of rows after filtering:", n))

        readepth <- mean(LOH[, 4])
        
        w303 <- numeric(n)
        yjm <- numeric(n)

        for (i in 1:n) {
            w303[i] <- (length(unlist(strsplit(LOH[i, 5], split = ",|\\.", fixed = FALSE))) - 1) / readepth
            yjm[i] <- (length(unlist(strsplit(toupper(LOH[i, 5]), split = LOH[i, 6], fixed = FALSE))) - 1) / readepth
        }

        mapWJLOH <- cbind.data.frame(LOH[, 1], as.numeric(LOH[, 2]), w303, yjm)

        # Process specific chromosomal regions
        for (i in 1:n) {
            if (mapWJLOH[i, 1] == "chrVIII" & (mapWJLOH[i, 2] > 211729 & mapWJLOH[i, 2] < 216235)) {
                mapWJLOH[i, 3] <- as.numeric(mapWJLOH[i, 3] / 6)
                mapWJLOH[i, 4] <- as.numeric(mapWJLOH[i, 4] / 3)
            }
            if (mapWJLOH[i, 1] == "chrXII" & (mapWJLOH[i, 2] > 451575 & mapWJLOH[i, 2] < 468812)) {
                mapWJLOH[i, 3] <- as.numeric(mapWJLOH[i, 3] / 60)
                mapWJLOH[i, 4] <- as.numeric(mapWJLOH[i, 4] / 30)
            }
        }

        # Output the result file
        output_file <- paste(base_name, "data", sep = "")
        print(paste("Output file:", output_file))
        write.table(mapWJLOH, file = output_file, row.names = FALSE, col.names = FALSE)

        # Read the chromosome length file
        chrom.length <- read.table("/home/dqlab/chromlength.txt", header = FALSE, sep = ' ', col.names = c('X', 'Length'))

        # Generate PDF plots
        pdf_file <- paste(base_name, "LOH.pdf", sep = "")
        print(paste("PDF file:", pdf_file))
        pdf(pdf_file, width = 11, height = 8)
        par(mfrow = c(4, 1), mar = c(2, 2, 2, 2))

        for (a in 1:16) {
            chrom_plot <- subset(mapWJLOH, mapWJLOH[, 1] == paste('chr', as.roman(a), sep = ''), select = c(2, 3, 4))
            X <- chrom_plot[, 1]
            w303mean <- chrom_plot[, 2]
            yjmmean <- chrom_plot[, 3]

            if (length(X) == 0 || length(w303mean) == 0 || length(yjmmean) == 0) {
                cat("Empty data for chromosome", as.roman(a), "- Skipping plotting for this chromosome.\n")
                next
            }

            plot(X, w303mean, type = 'p', ylim = c(-1, 2), col = 'red', main = paste('Chr', as.roman(a)))
            points(X, yjmmean, type = 'p', col = 'blue')
            abline(h = seq(-1, 2, by = 0.5), col = "lightgray", lty = 2)
            axis(1)
            box()
        }
        dev.off()

    }, error = function(e) {
        cat("An error occurred:", conditionMessage(e), "\n")
        next
    })
}