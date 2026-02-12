# Load the required library
library(dplyr)
library(readr)
library(here)

#Set local directory
local_dir <- here("Fig_S2A/Stat_analysis")

# Read input CSV file located in same folder as script
data <- read.csv(file.path(local_dir, "LADL-poly-overlap.csv"))
data <- na.omit(data)

# Create an empty matrix to store absolute differences
diff_table <- matrix(NA, nrow = nrow(data), ncol = nrow(data))

# Generate differences table
# Iterate through all possible pairs and calculate the absolute difference
for (i in 1:nrow(data)) {
  for (j in 1:nrow(data)) {
    diff_table[i, j] <- abs(data$p[i] - data$p[j])
  }
}
# Create a new data frame with the results
diff_table_df <- as.data.frame(diff_table)
row.names(diff_table_df) <- data$Condition
colnames(diff_table_df) <- data$Condition

#Generate critical range table
# Chi square critical value
confidence <- 0.95 ####p-value=0.05
DF <- nrow(data)-1
chisqcv <- qchisq(confidence,DF)
# Create an empty matrix to store the calculated values
cr_table <- matrix(NA, nrow = nrow(data), ncol = nrow(data))
# Iterate through all possible pairs and calculate the values
for (i in 1:nrow(data)) {
  for (j in 1:nrow(data)) {
    pi <- data$p[i]
    pj <- data$p[j]
    Ni <- data$N[i]
    Nj <- data$N[j]
    cr_table[i, j] <- sqrt(chisqcv) * sqrt((pi*(1-pi)/Ni) + (pj*(1-pj)/Nj))
  }
}
# Create a new data frame with the results
cr_table_df <- as.data.frame(cr_table)
row.names(cr_table_df) <- data$Condition
colnames(cr_table_df) <- data$Condition

#Compare differences table with critical range table
# Create an empty matrix to store the comparison results
comparison_matrix <- matrix(NA, nrow = nrow(diff_table_df), ncol = ncol(diff_table_df))
# Iterate through each position and compare values
for (i in 1:nrow(diff_table_df)) {
  for (j in 1:ncol(diff_table_df)) {
    if (diff_table_df[i, j] > cr_table_df[i, j]) {
      comparison_matrix[i, j] <- "S"
    } else {
      comparison_matrix[i, j] <- "NS"
    }
  }
}
# Create a new data frame with the comparison results
comparison_df <- as.data.frame(comparison_matrix)
row.names(comparison_df) <- data$Condition
colnames(comparison_df) <- data$Condition

#Write output files
write.csv(comparison_df, file.path(local_dir,"Output/Contacts_Marascuilo_results.csv"))