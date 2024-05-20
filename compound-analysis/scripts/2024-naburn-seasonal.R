# LOAD TIDYVERSE AND CHECK PACKAGE UPDATES
library(tidyverse)
library(fuzzyjoin)
library(stringi)
library(ggforce)
library(RColorBrewer)

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
antibiotics <-
  splitmasslist$"antibiotics_itn_msca_answer_160616_w_dtxsi_ds"
metabolites <- splitmasslist$"itnantibiotic_cyp_metabolites"
psychoactive <- splitmasslist$"kps_psychoactive_substances_v2"
pharmaceuticals <- splitmasslist$"kps_pharmaceuticals"

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
antibiotics$name <- antibiotics$name %>%
  fedmatch::clean_strings()

# add antibiotic classes for common antibiotics
class_info <-
  read.csv("compound-analysis/data/antimicrobial_classes.csv")
classes <- antibiotics %>%
  fuzzy_left_join(class_info,
                  by = c("name" = "name"),
                  match_fun = str_detect)

# replace NA with 'unknown'
classes$class <- classes$class %>%
  replace_na('unknown')

# Remove duplicated rows based on name.x, month, mean
duplicates <- classes %>% 
  distinct(name.x, month, group_area, .keep_all = TRUE)

# PRODUCE A CSV OF RESULTS
write.csv(
  antibiotics,
  "compound-analysis/data/processed-data/itn_antibiotics.csv",
  row.names = FALSE
)
write.csv(
  duplicates,
  "compound-analysis/data/processed-data/antibiotics_class_organised.csv",
  row.names = FALSE
)
write.csv(
  metabolites,
  "compound-analysis/data/processed-data/itn_metabolites.csv",
  row.names = FALSE
)
write.csv(
  psychoactive,
  "compound-analysis/data/processed-data/psychoactive.csv",
  row.names = FALSE
)
write.csv(
  pharmaceuticals,
  "compound-analysis/data/processed-data/pharmaceuticals.csv",
  row.names = FALSE
)


# wide view for samples and fully annotated view
classes_means <- duplicates %>%
  group_by(pick(class, day, month)) %>%
  summarise(
    mean = mean(group_area),
    std = sd(group_area),
    n = length(group_area),
    se = std / sqrt(n)
  )

# PLOT THE RESULTS

# split based on target antibiotics for location 
split <- split(classes_means, classes_means$class)
aminoglycoside <- split$aminoglycoside
beta <- split$'beta-lactam'
glycopeptide_metronidazole <- split$glycopeptide_metronidazole
macrolide_lincosamide <- split$macrolide_lincosamide
other <- split$other
phenicol <- split$phenicol
quinolone <- split$quinolone
sulfonamide_trimethoprim <- split$sulfonamide_trimethoprim
tetracycline <- split$tetracycline

# INDIVIDUAL --------------------------------------------------------------

aminoglycoside %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "peak intensity", fill = "sample day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

beta %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "peak intensity", fill = "sample day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

macrolide_lincosamide %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "peak intensity", fill = "sample day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

other %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "peak intensity", fill = "sample day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

phenicol %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "peak intensity", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~name.x, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

quinolone %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "peak intensity", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~name.x, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

sulfonamide_trimethoprim %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "peak intensity", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~name.x, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

tetracycline %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "peak intensity", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~name.x, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

