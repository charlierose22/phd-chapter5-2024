
# GRAPH ----------------------------------------------------------------
library(RColorBrewer)

# Create a list of target antibiotics.
class <- unique(means$class)


# split based on target antibiotics for location 
split <- split(means, means$class)
aminoglycoside <- split$aminoglycoside
beta <- split$'beta-lactam'
glycopeptide_metronidazole <- split$glycopeptide_metronidazole
macrolide_lincosamide <- split$macrolide_lincosamide
other <- split$other
phenicol <- split$phenicol
quinolone <- split$quinolone
sulfonamide_trimethoprim <- split$sulfonamide_trimethoprim
tetracycline <- split$tetracycline
mdr <- split$mdr
mge_integrons <- split$mge_integrons


# INDIVIDUAL --------------------------------------------------------------

aminoglycoside %>%
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

beta %>%
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

glycopeptide_metronidazole %>%
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

macrolide_lincosamide %>%
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

other %>%
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

phenicol %>%
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

quinolone %>%
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

sulfonamide_trimethoprim %>%
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

tetracycline %>%
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

mdr %>%
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

mge_integrons %>%
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

