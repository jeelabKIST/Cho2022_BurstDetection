# Install Libraries
if (!require(R.matlab)) install.packages("R.matlab")
if (!require(dunn.test)) install.packages("dunn.test")

# Import Libraries
library(R.matlab)
library(dunn.test)

# Read Data
data_dir <- "/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/Figure6"
data <- readMat(paste(data_dir, "ranked_mapping_theta.mat", sep = "/"))

ranked_mapping <- data$RANKED.MAPPING
ranked_dc <- ranked_mapping[[5]]

# Function to Remove NaN Values from a Vector
remove_nan <- function(x) {
  x[!is.nan(x)]
}

# Detection Confidence Values
group1 <- remove_nan(as.numeric(ranked_dc[[1]][[1]]))
group2 <- remove_nan(as.numeric(ranked_dc[[2]][[1]]))
group3 <- remove_nan(as.numeric(ranked_dc[[3]][[1]]))
group4 <- remove_nan(as.numeric(ranked_dc[[4]][[1]]))
group5 <- remove_nan(as.numeric(ranked_dc[[5]][[1]]))

ranked_dc <- list(group1, group2, group3, group4, group5)

# Kruskal-Wallis Test
kw_result <- kruskal.test(ranked_dc, na.action = na.omit)

# Post-hoc Dunn's Test
dunn_result <- dunn.test(ranked_dc, g = NA, method = "bonferroni",
                    kw = TRUE, label = TRUE, wrap = FALSE, table = TRUE,
                    list = FALSE, rmc = FALSE, alpha = 0.05, altp = FALSE)