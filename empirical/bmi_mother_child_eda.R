# Estimated pre-pregnancy BMI, derived from pre-pregnancy weight (self-reported at Week 11 visit) and height measured at Week 26 visit
# MGD05BMIP


# CWH05WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at birth; age is taken to be 0
# CWH06WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Day 1 visit
# CWH07WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Week 1 visit
# CWH08WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Week 3 visit
# CWH10WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Month 3 visit
# CWH11WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Month 6 visit
# CWH12WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Month 9 visit
# CWH13WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Month 12 visit
# CWH14WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Month 15 visit
# CWH15WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Month 18 visit
# CWH16WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Month 24 visit
# CWH17WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Month 36 visit
# CWH18WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Year 4 visit
# CWH19WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Year 4.5 visit
# CWH20WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Year 5 visit
# CWH21WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Year 5.5 visit
# CWH22WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Year 6 visit
# CWH23WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Year 6.5 visit
# CWH24WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Year 7 visit
# CWH26WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Year 8 visit
# CWH28WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Year 9 visit
# CWH30WBZ	WHO BMI-for-age z-score (adjusted for sex as well) at Year 10 visit
library(readxl)
all_data <- read_excel("empirical/data/FormA564_20221128.xlsx")
head(all_data)

# Extract the relevant columns
bmi_data <- all_data[, c("MGD05BMIP", "CWH05WBZ", "CWH06WBZ", "CWH07WBZ", "CWH08WBZ", 
                         "CWH10WBZ", "CWH11WBZ", "CWH12WBZ", "CWH13WBZ", "CWH14WBZ", 
                         "CWH15WBZ", "CWH16WBZ", "CWH17WBZ", "CWH18WBZ", "CWH19WBZ", 
                         "CWH20WBZ", "CWH21WBZ", "CWH22WBZ", "CWH23WBZ", "CWH24WBZ", 
                         "CWH26WBZ", "CWH28WBZ", "CWH30WBZ")]

# Plot scatterplots of MGD05BMIP with each of the other columns
par(mfrow = c(4, 6))  # Set up the plotting area to have multiple plots
for (col in names(bmi_data)[-1]) {
  plot(bmi_data$MGD05BMIP, bmi_data[[col]], 
       xlab = "MGD05BMIP", ylab = col, 
       main = paste("Scatterplot of MGD05BMIP vs", col))
}

library(ggplot2)

# Reshape the data for ggplot
bmi_data_long <- reshape2::melt(bmi_data, id.vars = "MGD05BMIP")

# Create scatterplots using ggplot
ggplot(bmi_data_long, aes(x = MGD05BMIP, y = value)) +
  geom_point() +
  facet_wrap(~ variable, scales = "free_y") +
  labs(x = "MGD05BMIP", y = "BMI-for-age z-score", 
       title = "Scatterplots of MGD05BMIP vs BMI-for-age z-scores") +
  theme_minimal()

  # Extract data pairs for the 3rd, 5th, 7th, and 9th plot
  data_pairs <- list(
    plot3 = bmi_data[, c("MGD05BMIP", names(bmi_data)[3])],
    plot5 = bmi_data[, c("MGD05BMIP", names(bmi_data)[5])],
    plot7 = bmi_data[, c("MGD05BMIP", names(bmi_data)[7])],
    plot9 = bmi_data[, c("MGD05BMIP", names(bmi_data)[9])]
  )
  
  # Save each data pair as an RDS file
  saveRDS(data_pairs$plot3, file = "empirical/data/plot3_data.rds")
  saveRDS(data_pairs$plot5, file = "empirical/data/plot5_data.rds")
  saveRDS(data_pairs$plot7, file = "empirical/data/plot7_data.rds")
  saveRDS(data_pairs$plot9, file = "empirical/data/plot9_data.rds")














