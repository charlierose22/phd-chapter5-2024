# load libraries and packages
library(tidyverse)

# Define the folder path (replace with your actual path)
folder_path <- "E:/phd-chapter4-2024/prescription-analysis/data"

# List all CSV files in the folder
all_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

# Loop through each file
for (file in all_files) {
  # Read the current file
  data <- read_csv(file) %>% janitor::clean_names()
  # Extract year and month from filename
  filename <- basename(file)
  year_month <- str_sub(filename, start = 7, end = 10)
  #subset larger csv
  data <- subset(data, bnf_chapter_plus_code == '05: Infections')
  data <- subset(data, icb_name == 'NHS HUMBER AND NORTH YORKSHIRE INTEGRATE')
  # create smaller filtered csv
  write.csv(data, file = paste0("prescription-analysis/filtered_data", "filtered_", year_month, ".csv"),
    row.names = FALSE)
}


data <- read_csv("E:/phd-chapter4-2024/prescription-analysis/data/epd_202306.csv") %>% janitor::clean_names()
# Extract year and month from filename
#subset larger csv
data <- subset(data, bnf_chapter_plus_code == '05: Infections')
data <- subset(data, icb_name == 'NHS HUMBER AND NORTH YORKSHIRE INTEGRATE')
# create smaller filtered csv
write.csv(data, file = "filtered_datafiltered_2306.csv", row.names = FALSE)