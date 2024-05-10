# load libraries and packages
library(tidyverse)

# Upload all data from csvs,

# Define the folder path
folder_path <- "prescription-analysis/data"

# List all CSV files in the folder
all_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

# Loop through each file
for (file in all_files) {
  # Extract file name without path
  filename <- gsub("\\.csv$", "", basename(file))  # removes extension
  
  # Read, clean, and assign the data directly to environment
  assign(filename, read_csv(file) %>% janitor::clean_names())
}

# List of existing dataframes
filtered_data_list <- c("filtered_2208", 
                        "filtered_2209", 
                        "filtered_2210", 
                        "filtered_2211", 
                        "filtered_2212", 
                        "filtered_2301",
                        "filtered_2302",
                        "filtered_2303",
                        "filtered_2304",
                        "filtered_2305",
                        "filtered_2306",
                        "filtered_2307",
                        "filtered_2308")

# Loop through each data frame
for (file in filtered_data_list) {
  # Extract file name without path
  filename <- basename(file)  # removes path
  
  # Extract year and month (assuming last 4 digits are yymm)
  year_month <- substr(filename, nchar(filename) - 3, nchar(filename))  # last 4 digits
  
  # Create new data frame name with prefix and year_month
  new_df_name <- paste0("YO_", year_month)
  
  # Get the current data frame from the list (replace with your access method)
  current_data <- get(gsub("\\.csv$", "", filename))  # remove extension for access
  
  # Filter data for postcodes starting with "YO"
  yo_filtered_data <- current_data %>% filter(grepl("^YO", postcode))  # filter by postcode
  
  # Assign the filtered data to a new object in the environment
  assign(new_df_name, yo_filtered_data)
}

# List of existing dataframes
YO_data_list <- c("YO_2208", 
                        "YO_2209", 
                        "YO_2210", 
                        "YO_2211", 
                        "YO_2212", 
                        "YO_2301",
                        "YO_2302",
                        "YO_2303",
                        "YO_2304",
                        "YO_2305",
                        "YO_2306",
                        "YO_2307",
                        "YO_2308")

for (data_frame_name in YO_data_list) {
  # Load the current data frame
  df <- get(data_frame_name)
  
  # Filter data for rows with "york" in address_3 or address_4
  filtered_df <- df[tolower(df$address_3) == "york" | tolower(df$address_4) == "york", ]
  
  # Assign new name with "york_" prefix
  new_name <- paste0("york_", data_frame_name)
  
  # Assign the filtered data frame to a new environment variable with the new name
  assign(new_name, filtered_df)
  
  # Optionally, print a message (e.g., to track progress)
  print(paste("Created filtered data frame:", new_name))
}




