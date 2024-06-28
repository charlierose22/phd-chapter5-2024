
# GRAPH ----------------------------------------------------------------
library(RColorBrewer)


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


# INDIVIDUAL --------------------------------------------------------------
means_class_genes %>%
  ggplot(aes(x = month, y = mean, colour = class)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Month", y = "Average Gene Abundance Relative to 16S", colour = "Antibiotic Class") +
  scale_color_d3(palette = "category10", labels = c("Aminoglycoside",
                                                    "Beta-lactam",
                                                    "Glycopeptide and\nMetronidazole",
                                                    "Macrolide and\nLincosamide",
                                                    "Other",
                                                    "Phenicol",
                                                    "Quinolone",
                                                    "Sulfonamide and\nTrimethoprim",
                                                    "Tetracycline")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")


aminoglycoside_seasonal %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

beta_seasonal %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

glycopeptide_metronidazole_seasonal %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

macrolide_lincosamide_seasonal %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

other_seasonal %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

phenicol_seasonal %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

quinolone_seasonal %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

sulfonamide_trimethoprim_seasonal %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

tetracycline_seasonal %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

mdr_seasonal %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

mge_integrons_seasonal %>%
  ggplot(aes(x = month, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", fill = "sample_day") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

