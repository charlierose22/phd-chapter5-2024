# LOAD TIDYVERSE AND CHECK PACKAGE UPDATES
library(tidyverse)
library(fuzzyjoin)
library(stringi)
library(ggforce)
library(RColorBrewer)
library(ggsci)
library(grafify)
library(rstatix)

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
day29_compounds <- day29_compounds[!grepl('unknown|antifungal',
                                              day29_compounds$class),]
day29_compounds$name.x <- str_to_title(day29_compounds$name.x)
day29_compounds$class <- str_to_title(day29_compounds$class)

day29_compounds$class <- str_replace_all(day29_compounds$class,'Beta-Lactam','Beta-lactam')
day29_compounds$class <- str_replace_all(day29_compounds$class,'Macrolide_lincosamide','Macrolides and Lincosamides')
day29_compounds$class <- str_replace_all(day29_compounds$class,'Sulfonamide_trimethoprim','Sulfonamides and Trimethoprim')


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

stacked_chem <- day29_compounds
# total
# stacked area plot
stacked_seasonal <- stacked_chem %>% 
  group_by(month, class) %>%
  summarise(n = sum(group_area)) %>%
  mutate(percentage = n / sum(n))
stacked_seasonal$percentage <- stacked_seasonal$percentage * 100

# time
ggplot(stacked_seasonal, aes(x = month,
                             y = percentage,
                             fill = class)) +
  geom_area(alpha=0.6 , linewidth=0.5, colour="black") +
  scale_fill_manual(values = brewer.pal("Spectral", n = 11)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  labs(x = "Month", y = "Proportion of UHPLC Peak Intensity (%)", fill = "Antibiotic Class") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# PLOT THE RESULTS

# split based on target antibiotics for location 
split_compounds <- split(means_compounds, means_compounds$class)
aminoglycoside_compounds <- split_compounds$Aminoglycoside
beta_compounds <- split_compounds$'Beta-lactam'
macrolide_lincosamide_compounds <- split_compounds$'Macrolides and Lincosamides'
other_compounds <- split_compounds$Other
phenicol_compounds <- split_compounds$Phenicol
quinolone_compounds <- split_compounds$Quinolone
sulfonamide_trimethoprim_compounds <- split_compounds$'Sulfonamides and Trimethoprim'
tetracycline_compounds <- split_compounds$Tetracycline

# INDIVIDUAL --------------------------------------------------------------
classes_means_compounds %>%
  ggplot(aes(x = month, y = mean, colour = class)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

aminoglycoside_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

beta_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

macrolide_lincosamide_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

other_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

phenicol_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

quinolone_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

sulfonamide_trimethoprim_compounds <- sulfonamide_trimethoprim_compounds[!grepl('5 4 Tert Butylphenyl Sulfanyl Quinazoline 2 4 Diamine|N 6 Aminohexyl 5 Chlor 1 Naphthalensulfonamid', sulfonamide_trimethoprim_compounds$name.x),]
sulfonamide_trimethoprim_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

tetracycline_compounds %>%
  ggplot(aes(x = month, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  scale_color_grafify(palette = "kelly") +
  labs(x = "Month", y = "Average Compound Intensity (UHPLC)", colour = "Antibiotic Compound Name") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# statistics

# replace days as factor for significance tests
chem_time_factor <- day29_compounds
chem_time_factor$month <- as.factor(chem_time_factor$month)
chem_time_factor$name.x <- as.factor(chem_time_factor$name.x)
chem_time_factor$class <- as.factor(chem_time_factor$class)

chem_split_time_factor <- split(chem_time_factor, chem_time_factor$class)
chem_time_factor_amino <- chem_split_time_factor$Aminoglycoside
chem_time_factor_beta <- chem_split_time_factor$'Beta-lactam'
chem_time_factor_mac <- chem_split_time_factor$'Macrolides and Lincosamides'
chem_time_factor_other <- chem_split_time_factor$Other
chem_time_factor_phen <- chem_split_time_factor$Phenicol
chem_time_factor_quin <- chem_split_time_factor$Quinolone
chem_time_factor_sulf <- chem_split_time_factor$'Sulfonamides and Trimethoprim'
chem_time_factor_tet <- chem_split_time_factor$Tetracycline


chem_time_factor_sulf <- chem_time_factor_sulf[!grepl('5 4 Tert Butylphenyl Sulfanyl Quinazoline 2 4 Diamine|N 6 Aminohexyl 5 Chlor 1 Naphthalensulfonamid', chem_time_factor_sulf$name.x),]

# time study
library(FSA)
# all classes
# one way anova
one_way_all_chem <- aov(group_area ~ class, data = chem_time_factor)
summary(one_way_all_chem)
# two way anova
two_way_all_chem <- aov(group_area ~ month * class, data = chem_time_factor)
summary(two_way_all_chem)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_all_chem)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor, group_area ~ class)
kruskal.test(data = chem_time_factor, group_area ~ month)
# post-hoc test
dunnTest(group_area ~ class,
         data = chem_time_factor,
         method = "bh")
dunnTest(group_area ~ month,
         data = chem_time_factor,
         method = "bh")

ggplot(chem_time_factor, aes(x = class, y = group_area) ) +
  geom_boxplot() +
  scale_x_discrete(name = "Antibiotic Class") +
  scale_y_continuous(name = "UHPLC Peak Intensity") +
  theme_classic()

# beta lactam
# one way anova
chem_one_way_beta <- aov(group_area ~ name.x, data = chem_time_factor_beta)
summary(chem_one_way_beta)
chem_one_way_beta_month <- aov(group_area ~ month, data = chem_time_factor_beta)
summary(chem_one_way_beta_month)
# two way anova
chem_two_way_beta <- aov(group_area ~ month * name.x, data = chem_time_factor_beta)
summary(chem_two_way_beta)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_beta)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_beta, group_area ~ name.x)
kruskal.test(data = chem_time_factor_beta, group_area ~ month)
# post-hoc test
dunnTest(group_area ~ month,
         data = chem_time_factor_beta,
         method = "bh")
dunnTest(group_area ~ name.x,
         data = chem_time_factor_beta,
         method = "bh")

# phenciol
# one way anova
chem_one_way_phen_month <- aov(group_area ~ month, data = chem_time_factor_phen)
summary(chem_one_way_phen_month)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_phen_month)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_phen, group_area ~ month)
# post-hoc test
dunnTest(group_area ~ month,
         data = chem_time_factor_phen,
         method = "bh")

# macrolide
# one way anova
chem_one_way_mac <- aov(group_area ~ name.x, data = chem_time_factor_mac)
summary(chem_one_way_mac)
chem_one_way_mac_month <- aov(group_area ~ month, data = chem_time_factor_mac)
summary(chem_one_way_mac_month)
# two way anova
chem_two_way_mac <- aov(group_area ~ month * name.x, data = chem_time_factor_mac)
summary(chem_two_way_mac)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_mac)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_mac, group_area ~ name.x)
kruskal.test(data = chem_time_factor_mac, group_area ~ month)
# post-hoc test
dunnTest(group_area ~ month,
         data = chem_time_factor_mac,
         method = "bh")
dunnTest(group_area ~ name.x,
         data = chem_time_factor_mac,
         method = "bh")

# quinolone
# one way anova
chem_one_way_quin <- aov(group_area ~ name.x, data = chem_time_factor_quin)
summary(chem_one_way_quin)
chem_one_way_quin_month <- aov(group_area ~ month, data = chem_time_factor_quin)
summary(chem_one_way_quin_month)
# two way anova
chem_two_way_quin <- aov(group_area ~ month * name.x, data = chem_time_factor_quin)
summary(chem_two_way_quin)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_quin)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_quin, group_area ~ name.x)
kruskal.test(data = chem_time_factor_quin, group_area ~ month)
# post-hoc test
dunnTest(group_area ~ month,
         data = chem_time_factor_quin,
         method = "bh")
dunnTest(group_area ~ name.x,
         data = chem_time_factor_quin,
         method = "bh")


# sulfonamide
# one way anova
chem_one_way_sulf <- aov(group_area ~ name.x, data = chem_time_factor_sulf)
summary(chem_one_way_sulf)
chem_one_way_sulf_month <- aov(group_area ~ month, data = chem_time_factor_sulf)
summary(chem_one_way_sulf_month)
# two way anova
chem_two_way_sulf <- aov(group_area ~ month * name.x, data = chem_time_factor_sulf)
summary(chem_two_way_sulf)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_sulf)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_sulf, group_area ~ name.x)
kruskal.test(data = chem_time_factor_sulf, group_area ~ month)
# post-hoc test
dunnTest(group_area ~ month,
         data = chem_time_factor_sulf,
         method = "bh")
dunnTest(group_area ~ name.x,
         data = chem_time_factor_sulf,
         method = "bh")
ggplot(chem_time_factor_sulf, aes(x = name.x, y = group_area) ) +
  geom_boxplot() +
  scale_x_discrete(name = "Antibiotic Name") +
  scale_y_continuous(name = "UHPLC Peak Intensity") +
  theme_classic()

# tetracycline
# one way anova
chem_one_way_tet_month <- aov(group_area ~ month, data = chem_time_factor_tet)
summary(chem_one_way_tet_month)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_tet_month)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_tet, group_area ~ month)
# post-hoc test
dunnTest(group_area ~ month,
         data = chem_time_factor_tet,
         method = "bh")

# other
# one way anova
chem_one_way_other <- aov(group_area ~ name.x, data = chem_time_factor_other)
summary(chem_one_way_other)
chem_one_way_other_month <- aov(group_area ~ month, data = chem_time_factor_other)
summary(chem_one_way_other_month)
# two way anova
chem_two_way_other <- aov(group_area ~ month * name.x, data = chem_time_factor_other)
summary(chem_two_way_other)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_other)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_other, group_area ~ name.x)
kruskal.test(data = chem_time_factor_other, group_area ~ month)
# post-hoc test
dunnTest(group_area ~ month,
         data = chem_time_factor_other,
         method = "bh")
dunnTest(group_area ~ name.x,
         data = chem_time_factor_other,
         method = "bh")
