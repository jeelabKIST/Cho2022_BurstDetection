# Install Libraries
if (!require(R.matlab)) install.packages("R.matlab")
if (!require(dunn.test)) install.packages("dunn.test")

# Import Libraries
library(R.matlab)
library(dunn.test)

# Read Data
data_dir <- "/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/Figure5"
data_bdur <- readMat(paste(data_dir, "burst_durations.mat", sep = "/"))
data_bsnr <- readMat(paste(data_dir, "burst_snrs.mat", sep = "/"))

burst_durations <- data_bdur$burst.durations
burst_snrs <- data_bsnr$burst.snrs

burst_durations <- list(burst_durations[[1]][[1]],
                        burst_durations[[2]][[1]],
                        burst_durations[[3]][[1]])
burst_snrs <- list(burst_snrs[[1]][[1]],
                   burst_snrs[[2]][[1]],
                   burst_snrs[[3]][[1]])

# Kruskal-Wallis Test
kw_result_dur <- kruskal.test(burst_durations)
kw_result_snr <- kruskal.test(burst_snrs)

# Post-hoc Dunn's Test
dunn_result_dur <- dunn.test(burst_durations, g = NA, method = "bonferroni",
                    kw = TRUE, label = TRUE, wrap = FALSE, table = TRUE,
                    list = FALSE, rmc = FALSE, alpha = 0.05, altp = FALSE)
dunn_result_snr <- dunn.test(burst_snrs, g = NA, method = "bonferroni",
                    kw = TRUE, label = TRUE, wrap = FALSE, table = TRUE,
                    list = FALSE, rmc = FALSE, alpha = 0.05, altp = FALSE)