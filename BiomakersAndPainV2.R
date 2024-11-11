### Step 1: Data Preparation ###

# Load necessary libraries
library(readxl)
library(dplyr)

# Read the data files
biomarkers <- read_excel("C:/Users/reemh/OneDrive/Documents/CAM MsC/Probability and Statistics/Rscript/biomarkers.xlsx")
covariates <- read_excel("C:/Users/reemh/OneDrive/Documents/CAM MsC/Probability and Statistics/Rscript/covariates.xlsx")

# Separate Patient ID and Time Point, filter for week 0, and merge
biomarkers_week0 <- biomarkers %>%
  mutate(PatientID = sub("-.*", "", `Biomarker`),        # Extract Patient ID
         TimePoint = sub(".*-", "", `Biomarker`)) %>%    # Extract Time Point
  filter(TimePoint == "0weeks") %>%                    # Filter for week 0 only
  select(-`Biomarker`)                                   # Remove original combined column

biomarkers_week0$PatientID <- as.numeric(biomarkers_week0$PatientID)


# Merge with covariates data
merged_data <- biomarkers_week0 %>%
  left_join(covariates, by = "PatientID")

# group by high VAS and Low VAS
merged_data$PainGroup <- ifelse(merged_data$`VAS-at-inclusion` >= 5, "High", "Low")



### Step 2: Visually asses the distribution of the random variables for each biomarker ###

### plot histograms for each biomarker data in the High VAS group
library(ggplot2)
library(gridExtra)

# Create an empty list to store plots
plots1 <- list()
# Loop through each biomarker for PainGroup = "High"
for (biomarker in biomarker_columns) {
  
  # Subset the data for the current biomarker and "High" PainGroup
  group_data <- merged_data[[biomarker]][merged_data$PainGroup == "High"]
  
  # Ensure group_data is numeric and remove NA values
  group_data <- as.numeric(na.omit(group_data))
  
  # Plot histogram with normal curve overlay
  p <- ggplot(data.frame(value = group_data), aes(x = value)) +
    geom_histogram(aes(y = ..density..), colour = "black", fill = "lightblue", bins = 15) +
    geom_density(colour = "blue", size = 1) +
    geom_function(fun = dnorm, colour = "red", size = 1,
                  args = list(mean = mean(group_data), sd = sd(group_data))) +
    labs(title = paste("Distribution of", biomarker, "in High Pain Group"),
         x = biomarker, y = "Density") +
    theme_minimal()
  
  plots1[[biomarker]] <- p
}
do.call(grid.arrange, c(plots1, ncol = 3))

### plot histograms for each biomarker data in the Low VAS group
plots2 <- list()
for (biomarker in biomarker_columns) {
  
  # Subset the data for the current biomarker and "Low" PainGroup
  group_data <- merged_data[[biomarker]][merged_data$PainGroup == "Low"]
  
  # Ensure group_data is numeric and remove NA values
  group_data <- as.numeric(na.omit(group_data))
  
  # Plot histogram with normal curve overlay
  p <- ggplot(data.frame(value = group_data), aes(x = value)) +
    geom_histogram(aes(y = ..density..), colour = "black", fill = "lightblue", bins = 15) +
    geom_density(colour = "blue", size = 1) +
    geom_function(fun = dnorm, colour = "red", size = 1,
                  args = list(mean = mean(group_data), sd = sd(group_data))) +
    labs(title = paste("Distribution of", biomarker, "in Low Pain Group"),
         x = biomarker, y = "Density") +
    theme_minimal()
  
  plots2[[biomarker]] <- p
}
do.call(grid.arrange, c(plots2, ncol = 3))

### plot Q-Q plot for each biomarker data in the High VAS group
# Loop through each biomarker for PainGroup = "High"
for (biomarker in biomarker_columns) {
  
  # Subset the data for the current biomarker and "High" PainGroup
  group_data <- merged_data[[biomarker]][merged_data$PainGroup == "High"]
  
  # Ensure group_data is numeric and remove NA values
  group_data <- as.numeric(na.omit(group_data))
  
  # Generate the Q-Q plot
  p <- ggplot(data.frame(sample = group_data), aes(sample = sample)) +
    geom_qq() +
    geom_qq_line() +
    labs(title = paste("Q-Q Plot of", biomarker, "in High Pain Group"),
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  # Display the plot
  print(p)
}

### plot Q-Q plots for each biomarker data in the Low VAS group
# Loop through each biomarker for PainGroup = "Low"
for (biomarker in biomarker_columns) {
  
  # Subset the data for the current biomarker and "Low" PainGroup
  group_data <- merged_data[[biomarker]][merged_data$PainGroup == "Low"]
  
  # Ensure group_data is numeric and remove NA values
  group_data <- as.numeric(na.omit(group_data))
  
  # Generate the Q-Q plot
  p <- ggplot(data.frame(sample = group_data), aes(sample = sample)) +
    geom_qq() +
    geom_qq_line() +
    labs(title = paste("Q-Q Plot of", biomarker, "in Low Pain Group"),
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  # Display the plot
  print(p)
}


library(ggplot2)
library(gridExtra)

# Subset the data for CXCL1 in the "High" PainGroup
cxcl1_data <- merged_data$CXCL1[merged_data$PainGroup == "High"]
cxcl1_data <- na.omit(as.numeric(cxcl1_data))  # Remove NA values and ensure numeric type

# Create a list to store Q-Q plots
qqplots <- list()

# Create the Q-Q plot for the actual CXCL1 data
qqplots[[1]] <- ggplot(data.frame(sample = cxcl1_data), aes(sample = sample)) +
  geom_qq() + geom_qq_line() +
  ggtitle("Actual CXCL1 data in High Pain Group")

# Compute the sample size based on the actual data
n <- length(cxcl1_data)

# Generate 8 simulated datasets from a normal distribution and create Q-Q plots for each
for (i in 2:9) {
  simulated_data <- data.frame(normal_data = rnorm(n, mean(cxcl1_data), sd(cxcl1_data)))
  qqplots[[i]] <- ggplot(simulated_data, aes(sample = normal_data)) +
    geom_qq() + geom_qq_line() +
    ggtitle("Simulated normal data")
}

# Plot the resulting Q-Q plots side-by-side using gridExtra
grid.arrange(grobs = qqplots, ncol = 3)


### Step 3: conduct t.test of the hypothesis ###

# Define biomarker column names
biomarker_columns <- c("IL-8", "VEGF-A", "OPG", "TGF-beta-1", "IL-6", 
                       "CXCL9", "CXCL1", "IL-18", "CSF-1")

library(dplyr)


alpha <- 0.05


# Loop through each biomarker and conduct Welch's t-test
for (biomarker in biomarker_columns) {
  # Perform Welch's t-test
  t_test_result <- t.test(merged_data[[biomarker]] ~ merged_data$PainGroup, var.equal = FALSE)
  
  # Print the biomarker name and p-value only
  cat("Biomarker:", biomarker, "| p-value:", t_test_result$p.value, "\n")
  
  # Interpret the p-value
  if (t_test_result$p.value < alpha) {
    cat("Result: Significant difference between groups\n\n")
  } else {
    cat("Result: No significant difference between groups\n\n")
  }
}

typeI_error<-1-(1-.05)^9
print(typeI_error)

# adjusted alpha for 9 tests
adjusted_alpha <- 0.05 / 9


# Loop through each biomarker and conduct Welch's t-test
for (biomarker in biomarker_columns) {
  # Perform Welch's t-test
  t_test_result <- t.test(merged_data[[biomarker]] ~ merged_data$PainGroup, var.equal = FALSE)
  
  # Print the biomarker name and p-value only
  cat("Biomarker:", biomarker, "| p-value:", t_test_result$p.value, "\n")
  
  # Interpret the p-value
  if (t_test_result$p.value < adjusted_alpha) {
    cat("Result: Significant difference between groups\n\n")
  } else {
    cat("Result: No significant difference between groups\n\n")
  }
}

library(MKinfer)
# Loop through each biomarker and perform permutation t-test
for (biomarker in biomarker_columns) {
  result <- perm.t.test(merged_data[[biomarker]] ~ merged_data$PainGroup)
  print(result)
  #cat("Biomarker:", biomarker, "| Observed Difference:", result$estimate, "| p-value:", result$p.value, "\n\n")
}



### Task 2: Regression Model ####

# Remove rows with NA in Vas-12months for modeling
cleaned_data <- merged_data[!is.na(merged_data$`Vas-12months`), ]
# Convert categorical variables to factors
cleaned_data$`Sex (1=male, 2=female)` <- as.factor(cleaned_data$`Sex (1=male, 2=female)`)
cleaned_data$`Smoker (1=yes, 2=no)` <- as.factor(cleaned_data$`Smoker (1=yes, 2=no)`)

# Split data into training (80%) and test (20%) sets
set.seed(123)  # For reproducibility
sample_indices <- sample(1:nrow(merged_data), size = 0.8 * nrow(cleaned_data))
training_set <- cleaned_data[sample_indices, ]
test_set <- cleaned_data[-sample_indices, ]


# Fit the linear regression model
model <- lm(`Vas-12months` ~ `IL-8` + `VEGF-A` + `OPG` + `TGF-beta-1` + `IL-6` + 
              `CXCL9` + `CXCL1` + `IL-18` + `CSF-1` + `VAS-at-inclusion` + Age + 
              `Smoker (1=yes, 2=no)` + `Sex (1=male, 2=female)`, 
            data = training_set)

# Part (a): Model summary and coefficients
summary(model)
coef(model)

plot(model)


# Predicting on the test set
predictions <- predict(model, newdata = test_set)

mae <- mean(abs(predictions - test_set$`Vas-12months`))
mse <- mean((predictions - test_set$`Vas-12months`)^2)
rmse <- sqrt(mse)

# Print metrics
cat("Mean Absolute Error:", mae, "\n")
cat("Mean Squared Error:", mse, "\n")
cat("Root Mean Squared Error:", rmse, "\n")

# Plot Actual vs. Predicted Values
library(ggplot2)
ggplot(test_set, aes(x = predictions, y = `Vas-12months`)) +
  geom_point() +  # Scatter plot of predicted vs actual values
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # 45-degree line for reference
  labs(title = "Actual vs. Predicted Values of VAS at 12 Months",
       x = "Predicted VAS at 12 Months",
       y = "Actual VAS at 12 Months") +
  theme_minimal()

library(ggfortify)
autoplot(model, which = 1:6, ncol = 2, label.size = 3)

comparison_df <- data.frame(
  Actual_VAS_12months = test_set$`Vas-12months`, # Replace `Vas-12months` with the exact column name if different
  Predicted_VAS_12months = predictions
)

# Display the comparison table
print(comparison_df)












