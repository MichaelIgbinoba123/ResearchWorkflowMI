# Install the googlesheets4 package if not already installed
if (!requireNamespace("googlesheets4", quietly = TRUE)) {
  install.packages("googlesheets4")
}

# Load the required library
library(googlesheets4)

# Read Google Sheets data into R
spreadsheet_url <- "https://docs.google.com/spreadsheets/d/1w01kv29A6fG3sT1Bbo5G6fYOxC525i0z5mR9-iJ5W5c/edit#gid=1532556392"
your_data <- range_read(spreadsheet_url, sheet = "Data Accumulation")

# Remove spaces and special characters from variable names
names(your_data) <- make.names(names(your_data))

# Filter out unwanted variables
selected_variables <- your_data[, !(names(your_data) %in% c("Storm", "Date.Time", "Lon", "ENVP", "Pressure", "Cold.Ring.Eye.Rad"))]

# Check for constant variables
constant_vars <- sapply(selected_variables, function(x) length(unique(x)) == 1)
selected_variables <- selected_variables[, !constant_vars]

# Remove rows with missing values
selected_variables <- na.omit(selected_variables)

# Build the multiple linear regression model
model <- lm(`Pressure.Deficit` ~ ., data = selected_variables)

# Print the summary of the regression model
summary(model)

# Output the correlation matrix for all variables in the model
cor_matrix <- cor(selected_variables)
print(cor_matrix)

