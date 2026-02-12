# Load the necessary packages
library(here)

#Set local directory
local_dir <- here("Fig_S4C/Stat_analysis")

# Read input CSV file located in same folder as script
data <- read.csv(file.path(local_dir, "HeLa-RNAF-intensity.csv"))
data <- na.omit(data)

# Perform the two-factor ANOVA
model <- aov(dependent_variable ~ factor1 * factor2, data = data)

# Check the ANOVA table
anova_table <- anova(model)
write.csv(anova_table, file.path(local_dir,"Output/int-anova.csv"))

# Conduct post hoc tests if significant interaction or main effects are found
if (anova_table$"Pr(>F)"[1] < 0.05 || anova_table$"Pr(>F)"[2] < 0.05 || anova_table$"Pr(>F)"[3] < 0.05) {
  posthoc_factor1 <- TukeyHSD(aov(dependent_variable ~ factor1, data = data))
  posthoc_factor2 <- TukeyHSD(aov(dependent_variable ~ factor2, data = data))
  posthoc_interaction <- TukeyHSD(aov(dependent_variable ~ factor1 * factor2, data = data))

  #Convert to data frames  
  posthoc_factor1_df <- as.data.frame(posthoc_factor1$factor1)
  posthoc_factor2_df <- as.data.frame(posthoc_factor2$factor2)
  posthoc_interaction_df <- as.data.frame(posthoc_interaction$`factor1:factor2`)
  
  #Write output files
  write.csv(posthoc_factor1_df, file.path(local_dir,"Output/int-posthoc_factor1.csv"))
  write.csv(posthoc_factor2_df, file.path(local_dir,"Output/int-posthoc_factor2.csv"))
  write.csv(posthoc_interaction_df, file.path(local_dir,"Output/int-posthoc_interaction.csv"))
}
