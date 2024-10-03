# LOAD TIDYVERSE AND CHECK PACKAGE UPDATES
library(tidyverse)
library(fuzzyjoin)
library(stringi)
library(ggforce)
library(RColorBrewer)
library(ggsci)

# IMPORT YOUR CD DATA
seasonal_study <- readr::read_delim(
  "compound-analysis/data/raw-data/Compounds.csv",
  delim = "\t",
  trim_ws = TRUE
) %>%
  janitor::clean_names()

# RENAME DATAFRAME FOR CODE TO WORK WITH MINIMAL CHANGES
basedata <- seasonal_study

# DROP THESE COLUMNS UNLESS THEY HAVE BEEN USED IN YOUR CD WORKFLOW
basedata$tags = NULL
basedata$checked = NULL

# FILTER SAMPLES WITH NO COMPOUND NAMES AND NO MS2 DATA (IF NEEDED)
basedatawithcompoundnames <- with(basedata, basedata[!(name == "" |
                                                         is.na(name)),])
ms2dataonly <- basedatawithcompoundnames[!grepl('No MS2',
                                                basedatawithcompoundnames$ms2), ]

# CREATE A UNIQUE IDENTIFIER FOR EACH FEATURE USING CONCATENATION
# CHANGE BASEDATAWITHCOMPOUNDNAMES TO MS2DATAONLY IF YOU CHOOSE TO RUN THAT CODE LINE
uniqueid <-
  add_column(basedatawithcompoundnames,
             unique_id = NA,
             .after = 0)
uniqueid$unique_id <- str_c(uniqueid$name, "_", uniqueid$rt_min)

# PEAK NUMBERS CAN BE USED AS ANOTHER IDENTIFIER
peaknumber <- add_column(uniqueid, peak_number = NA, .after = 0)
peaknumber$peak_number <- seq.int(nrow(peaknumber))

# REMOVE CD FILE NUMBERS FROM THE END OF SAMPLE NAMES
colnames(peaknumber) <-
  sub("*_raw_f\\d\\d*", "", colnames(peaknumber))

# LENGTHEN THE TABLE TO REMOVE WHITESPACE
longer <- peaknumber %>%
  pivot_longer(
    cols = group_area_b4_a:peak_rating_qc_q,
    names_to = "sample",
    values_to = "result"
  )

# CREATE A SAMPLE NAME COLUMN AND FILL, SO WE CAN GROUP PEAK RATING AND GROUP AREA
samplenames <- add_column(longer, measurement = NA)
samplenames <- mutate(samplenames,
                      measurement = case_when(
                        str_detect(sample,
                                   "group_area") ~
                          "group_area",
                        str_detect(sample,
                                   "peak_rating") ~
                          "peak_rating"
                      ))

# CLEAN SAMPLE NAME COLUMN
samplenames$sample <-
  str_replace_all(samplenames$sample, "group_area_", "")
samplenames$sample <-
  str_replace_all(samplenames$sample, "peak_rating_", "")

# WIDEN TABLE
wider <- samplenames %>%
  pivot_wider(names_from = measurement, values_from = result)

# REMOVE NAS
nona <- drop_na(wider, group_area)

# FILTER DEPENDING ON PEAK RATING NUMBER
peakrating <- subset(nona, peak_rating > 5)

# FILTER DEPENDING ON INTENSITY
grouparea <- subset(peakrating, group_area > 100000)

# FOR TECHNICAL REPLICATES ----
# CREATE A NEW COLUMN
replicates <- add_column(grouparea, replicate = NA)

# MAKE SURE TECHNICAL REPLICATES ARE AT THE END OF THE SAMPLE NAME AND CHANGE A/B/C ACCORDINGLY
replicates <- mutate(
  replicates,
  replicate = case_when(
    str_ends(sample, "a") ~ "a",
    str_ends(sample, "b") ~ "b",
    str_ends(sample, "c") ~ "c",
    str_ends(sample, "d") ~ "d",
    str_ends(sample, "e") ~ "e",
    str_ends(sample, "f") ~ "f",
    str_ends(sample, "g") ~ "g",
    str_ends(sample, "h") ~ "h",
    str_ends(sample, "i") ~ "i",
    str_ends(sample, "j") ~ "j",
    str_ends(sample, "k") ~ "k",
    str_ends(sample, "l") ~ "l",
    str_ends(sample, "m") ~ "m",
    str_ends(sample, "n") ~ "n",
    str_ends(sample, "o") ~ "o",
    str_ends(sample, "p") ~ "p",
    str_ends(sample, "q") ~ "q"
  )
)

# ADD A SAMPLE LOCATION COLUMN AND CLEAN TO REMOVE REPLICATE NAMES SO WE CAN REMOVE SOLOS
replicates$sample_location = replicates$sample
replicates$sample_location <-
  gsub('_.*', '', replicates$sample_location)

# REMOVE PEAKS WITH RESULTS IN ONLY ONE REPLICATE
soloremoved <-
  plyr::ddply(replicates, c("unique_id", "sample_location"),
              function(d) {
                if (nrow(d) > 1)
                  d
                else
                  NULL
              })

# INCLUDE SAMPLE INFO IF NEEDED
sample_info <-
  readr::read_csv("compound-analysis/data/sample-information/naburn-seasonal.csv")
sample_included <- soloremoved %>% left_join(sample_info,
                                             by = "sample_location")
sample_nona <- drop_na(sample_included, day)

#----
# SPLIT RESULTS BASED ON MASS LIST VS MZCLOUD
split <-
  split(sample_nona,
        sample_nona$annot_source_mass_list_search)
mzcloud <- split$"No results"
masslists <- split$"Full match"

# MERGE MASS LISTS INTO ONE COLUMN
masslistmerged <- masslists %>%
  pivot_longer(
    cols = c(starts_with("mass_list_match")) ,
    names_to = "mass_list_name",
    names_prefix = "mass_list_match_",
    values_to = "mass_list_match"
  )

# FILTER FOR NO MATCHES AND INVALID MASS RESULTS
filteredmasslist <- masslistmerged[!grepl('No matches found',
                                          masslistmerged$mass_list_match),]
filteredmzcloud <- mzcloud[!grepl('Invalid Mass',
                                  mzcloud$annot_source_mz_cloud_search),]

# SPLIT THE INDIVIDUAL MASS LISTS
# FIRST FIND MASS LIST NAMES (APPEARING IN CONSOLE)
unique(filteredmasslist$mass_list_name)
# THEN SPLIT BY MASS LIST
splitmasslist <-
  split(filteredmasslist, filteredmasslist$mass_list_name)
antibiotics_compounds <-
  splitmasslist$"antibiotics_itn_msca_answer_160616_w_dtxsi_ds"
metabolites_compounds <- splitmasslist$"itnantibiotic_cyp_metabolites"
psychoactive_compounds <- splitmasslist$"kps_psychoactive_substances_v2"
pharmaceuticals_compounds <- splitmasslist$"kps_pharmaceuticals"

# wide view for samples and fully annotated view
#antibiotics_means <- antibiotics %>%
#  group_by(pick(name, day, month)) %>%
#  summarise(
#    mean = mean(group_area),
#    std = sd(group_area),
#    n = length(group_area),
#    se = std / sqrt(n)
#  )

# clean compound names in antibiotics dataset.
antibiotics_compounds$name <- antibiotics_compounds$name %>%
  fedmatch::clean_strings()

# add antibiotic classes for common antibiotics
class_info <-
  read.csv("compound-analysis/data/antimicrobial_classes.csv")
classes_compounds <- antibiotics_compounds %>%
  fuzzy_left_join(class_info,
                  by = c("name" = "name"),
                  match_fun = str_detect)

# replace NA with 'unknown'
classes_compounds$class <- classes_compounds$class %>%
  replace_na('unknown')

# Remove duplicated rows based on name.x, month, mean
duplicates_compounds <- classes_compounds %>% 
  distinct(name.x, month, group_area, .keep_all = TRUE)

duplicates_compounds$month <- ym(duplicates_compounds$month)

# PRODUCE A CSV OF RESULTS
write.csv(
  antibiotics_compounds,
  "compound-analysis/data/processed-data/itn_antibiotics.csv",
  row.names = FALSE
)
write.csv(
  duplicates_compounds,
  "compound-analysis/data/processed-data/antibiotics_class_organised.csv",
  row.names = FALSE
)
write.csv(
  metabolites_compounds,
  "compound-analysis/data/processed-data/itn_metabolites.csv",
  row.names = FALSE
)
write.csv(
  psychoactive_compounds,
  "compound-analysis/data/processed-data/psychoactive.csv",
  row.names = FALSE
)
write.csv(
  pharmaceuticals_compounds,
  "compound-analysis/data/processed-data/pharmaceuticals.csv",
  row.names = FALSE
)

day29_compounds <- duplicates_compounds[!grepl('1',
                                          duplicates_compounds$day),]
day29_compounds <- day1_compounds[!grepl('unknown',
                                              day1_compounds$class),]

# wide view for samples and fully annotated view
classes_means_compounds <- day29_compounds %>%
  group_by(pick(class, month)) %>%
  summarise(
    mean = mean(group_area),
    std = sd(group_area),
    n = length(group_area),
    se = std / sqrt(n)
  )

means_compounds <- day29_compounds %>%
  group_by(pick(name.x, class, month)) %>%
  summarise(
    mean = mean(group_area),
    std = sd(group_area),
    n = length(group_area),
    se = std / sqrt(n)
  )

means_compounds$name.x <- str_to_title(means_compounds$name.x)

# PLOT THE RESULTS

# split based on target antibiotics for location 
split_compounds <- split(means_compounds, 
                         means_compounds$class)
aminoglycoside_compounds <- split_compounds$aminoglycoside
beta_compounds <- split_compounds$'beta-lactam'
glycopeptide_metronidazole_compounds <- split_compounds$glycopeptide_metronidazole
macrolide_lincosamide_compounds <- split_compounds$macrolide_lincosamide
other_compounds <- split_compounds$other
phenicol_compounds <- split_compounds$phenicol
quinolone_compounds <- split_compounds$quinolone
sulfonamide_trimethoprim_compounds <- split_compounds$sulfonamide_trimethoprim
tetracycline_compounds <- split_compounds$tetracycline

# INDIVIDUAL --------------------------------------------------------------
classes_means_compounds %>%
  ggplot(aes(x = month, y = mean, colour = class)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Class") +
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

aminoglycoside_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

beta_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

macrolide_lincosamide_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

other_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

phenicol_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

quinolone_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

sulfonamide_trimethoprim_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

tetracycline_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

