# load libraries and packages
library(tidyverse)
library(RColorBrewer)
library(grafify)
library(ggsci)
library(ggstatsplot)
library(ggrepel)
library(ggtext)

# Upload all data from csvs,
BNF_classes <- read_csv("prescription-analysis/API/API_BNF.csv")

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

# Append all dataframes
york_bind <- bind_rows(YO_2208,
                       YO_2209,
                       YO_2210,
                       YO_2211,
                       YO_2212,
                       YO_2301,
                       YO_2302,
                       YO_2303,
                       YO_2304,
                       YO_2305,
                       YO_2306,
                       YO_2307,
                       YO_2308)

# add class information so we can filter out anti-fungals and anti-virals etc
classes_info <- york_bind %>% left_join(BNF_classes, by = "bnf_chemical_substance")
antibiotics_prescriptions <- classes_info %>% drop_na(abx_class)

# set date
antibiotics_prescriptions$year_month <- ym(antibiotics_prescriptions$year_month)

abx <- antibiotics_prescriptions %>%
  group_by(pick(abx_name, year_month, abx_class)) %>%
  summarise(sum = sum(total_quantity))

abx$abx_name <- str_to_title(abx$abx_name)

# split based on target antibiotics for location 
split_prescriptions <- split(abx, abx$abx_class)
aminoglycoside_prescriptions <- split_prescriptions$aminoglycoside
beta_prescriptions <- split_prescriptions$'beta-lactam'
glycopeptide_metronidazole_prescriptions <- split_prescriptions$glycopeptide_metronidazole
macrolide_lincosamide_prescriptions <- split_prescriptions$macrolide_lincosamide
other_prescriptions <- split_prescriptions$other
phenicol_prescriptions <- split_prescriptions$phenicol
quinolone_prescriptions <- split_prescriptions$quinolone
sulfonamide_trimethoprim_prescriptions <- split_prescriptions$sulfonamide_trimethoprim
tetracycline_prescriptions <- split_prescriptions$tetracycline

total_prescriptions <- antibiotics_prescriptions %>% 
  group_by(year_month, abx_class) %>% 
  summarise(
    sum = sum(total_quantity))

darkpalette <- c("firebrick4",
                 "orangered1",
                 "darkgoldenrod2",
                 "yellow",
                 "darkolivegreen4", 
                 "green3",
                 "turquoise3", 
                 "dodgerblue3", 
                 "darkslateblue", 
                 "darkorchid2",
                 "violetred2",
                 "maroon")

# plot
total_prescriptions %>%
  ggplot(aes(x = year_month, y = sum, colour = abx_class)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Total Quantity of Medication Prescribed (units)", colour = "Antibiotic Class") +
  scale_color_d3(palette = "category10", labels = c("Aminoglycoside",
                                                    "Beta-lactam",
                                                    "Glycopeptide and\nMetronidazole",
                                                    "Macrolide and\nLincosamide",
                                                    "Other",
                                                    "Phenicol",
                                                    "Quinolone",
                                                    "Sulfonamide and\nTrimethoprim",
                                                    "Tetracycline")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Aminoglycosides
aminoglycoside_prescriptions %>%
  ggplot(aes(x = year_month, y = sum, colour = abx_name)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Total Quantity of Medication Prescribed (units)", colour = "Antibiotic Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Beta-lactams
beta_prescriptions %>%
  ggplot(aes(x = year_month, y = sum, colour = abx_name)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Total Quantity of Medication Prescribed (units)", colour = "Antibiotic Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Glycopeptide and Metronidazole
glycopeptide_metronidazole_prescriptions %>%
  ggplot(aes(x = year_month, y = sum, colour = abx_name)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Total Quantity of Medication Prescribed (units)", colour = "Antibiotic Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Macrolide and Lincosamide
macrolide_lincosamide_prescriptions %>%
  ggplot(aes(x = year_month, y = sum, colour = abx_name)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Total Quantity of Medication Prescribed (units)", colour = "Antibiotic Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Other
other_prescriptions %>%
  ggplot(aes(x = year_month, y = sum, colour = abx_name)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Total Quantity of Medication Prescribed (units)", colour = "Antibiotic Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Phenicol
phenicol_prescriptions %>%
  ggplot(aes(x = year_month, y = sum, colour = abx_name)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Total Quantity of Medication Prescribed (units)", colour = "Antibiotic Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Quinolone
quinolone_prescriptions %>%
  ggplot(aes(x = year_month, y = sum, colour = abx_name)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Total Quantity of Medication Prescribed (units)", colour = "Antibiotic Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Sulfonamide and Trimethoprim
sulfonamide_trimethoprim_prescriptions %>%
  ggplot(aes(x = year_month, y = sum, colour = abx_name)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Total Quantity of Medication Prescribed (units)", colour = "Antibiotic Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Tetracycline
tetracycline_prescriptions %>%
  ggplot(aes(x = year_month, y = sum, colour = abx_name)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Total Quantity of Medication Prescribed (units)", colour = "Antibiotic Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
