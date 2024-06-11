# load libraries and packages
library(tidyverse)
library(RColorBrewer)
library(ggstatsplot)
library(ggrepel)
library(ggtext)

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
antibiotics <- classes_info %>% drop_na(abx_class)

# set date
antibiotics$year_month <- ym(antibiotics$year_month)

# split based on target antibiotics for location 
split <- split(antibiotics, antibiotics$abx_class)
aminoglycoside <- split$aminoglycoside
beta <- split$'beta-lactam'
glycopeptide_metronidazole <- split$glycopeptide_metronidazole
macrolide_lincosamide <- split$macrolide_lincosamide
other <- split$other
phenicol <- split$phenicol
quinolone <- split$quinolone
sulfonamide_trimethoprim <- split$sulfonamide_trimethoprim
tetracycline <- split$tetracycline

means_months <- antibiotics %>%
  group_by(pick(abx_name, year_month, abx_class)) %>%
  summarise(
    mean = mean(total_quantity),
    std = sd(total_quantity),
    n = length(total_quantity),
    se = std / sqrt(n)
  )

means_total <- antibiotics %>% 
  group_by(year_month, abx_class) %>% 
  summarise(
    mean = mean(total_quantity),
    std = sd(total_quantity),
    n = length(total_quantity),
    se = std / sqrt(n)
  )

darkpalette <- c("firebrick4",
                 "orangered1",
                 "darkgoldenrod2",
                 "yellow3",
                 "darkolivegreen4", 
                 "green3",
                 "turquoise3", 
                 "dodgerblue3", 
                 "darkslateblue", 
                 "darkorchid2",
                 "violetred2",
                 "maroon")



# plot
means_total %>%
  ggplot(aes(x = year_month, y = mean, colour = abx_class)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = year_month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "month", y = "total quantity prescribed") +
  scale_color_manual(values = darkpalette) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  theme_minimal(base_size = 12)
