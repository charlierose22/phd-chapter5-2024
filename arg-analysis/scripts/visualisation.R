# GRAPH ----------------------------------------------------------------
library(RColorBrewer)
library(ggsci)
library(grafify)


means$month <- ym(means$month)
means_class_genes$month <- ym(means_class_genes$month)
# Create a list of target antibiotics.
class <- unique(means$class)



# split based on target antibiotics for location 
split_seasonal <- split(means, means$class)
aminoglycoside_seasonal <- split_seasonal$aminoglycoside
beta_seasonal <- split_seasonal$'beta-lactam'
glycopeptide_metronidazole_seasonal <- split_seasonal$glycopeptide_metronidazole
macrolide_lincosamide_seasonal <- split_seasonal$macrolide_lincosamide
other_seasonal <- split_seasonal$other
phenicol_seasonal <- split_seasonal$phenicol
quinolone_seasonal <- split_seasonal$quinolone
sulfonamide_trimethoprim_seasonal <- split_seasonal$sulfonamide_trimethoprim
tetracycline_seasonal <- split_seasonal$tetracycline
mdr_seasonal <- split_seasonal$mdr
mge_integrons_seasonal <- split_seasonal$mge_integrons

# total
# stacked area plot
stacked_seasonal <- day29_genes %>% 
  group_by(month, class) %>%
  summarise(n = sum(delta_ct)) %>%
  mutate(percentage = n / sum(n))

# time
ggplot(stacked_seasonal, aes(x = month,
                            y = percentage,
                            fill = class)) + 
  geom_area(alpha = 0.6 , size = 0.5, colour = "black") +
  scale_fill_manual(values = brewer.pal("Spectral", n = 11)) +
  labs(x = "Month", y = "Proportion of gene abundance", fill = "Target antibiotic class") +
  theme_bw(base_size = 12)


# INDIVIDUAL --------------------------------------------------------------
# total gene
means_class_genes %>%
  ggplot(aes(x = month, y = mean, colour = class)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Aminoglycosides
aminoglycoside_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Gene") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Beta-lactam
beta_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Gene") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Glycopeptide and Metronidazole
glycopeptide_metronidazole_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Gene") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Macrolide and Lincosamide
macrolide_lincosamide_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Gene") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Other
other_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Gene") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Phenicol
phenicol_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Gene") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Quinolone
quinolone_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Gene") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Sulfonamide and Trimethoprim
sulfonamide_trimethoprim_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Gene") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Tetracycline
tetracycline_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Gene") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# MDR
mdr_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Gene") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# MGE and Integrons
mge_integrons_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Gene") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

