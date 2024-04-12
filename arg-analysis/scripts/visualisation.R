
# LOCATION ----------------------------------------------------------------

library(RColorBrewer)
darkpalette <- c("mediumblue", 
                 "darkred", 
                 "darkorange4", 
                 "orange", 
                 "darkgreen", 
                 "green", 
                 "darkmagenta", 
                 "hotpink", 
                 "darkturquoise", 
                 "tomato")

# Create a list of target antibiotics.
target_antibiotics <- unique(assay_samples$target_antibiotics_major)

# split based on target antibiotics for location 
split <- split(means, means$target_antibiotics_major)
amino <- split$Aminoglycoside
beta <- split$`Beta Lactam`
int <- split$Integron
mdr <- split$MDR
mlsb <- split$MLSB
mge <- split$MGE
other <- split$Other
phen <- split$Phenicol
sulf <- split$Sulfonamide
tet <- split$Tetracycline
quin <- split$Quinolone
vanc <- split$Vancomycin
trim <- split$Trimethoprim


# INDIVIDUAL --------------------------------------------------------------

amino %>%
  ggplot(aes(x = month, y = mean, fill = age)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", color = "age") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

beta %>% 
  ggplot(aes(x = month, y = mean, fill = age)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", color = "age") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

mdr %>% 
  ggplot(aes(x = month, y = mean, fill = age)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", color = "age") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

mge %>% 
  ggplot(aes(x = month, y = mean, fill = age)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", color = "age") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

mlsb %>%
  ggplot(aes(x = month, y = mean, fill = age)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", color = "age") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

other %>%
  ggplot(aes(x = month, y = mean, fill = age)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", color = "age") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

phen %>%
  ggplot(aes(x = month, y = mean, fill = age)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", color = "age") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

quin %>% 
  ggplot(aes(x = month, y = mean, fill = age)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", color = "age") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

sulf %>% 
  ggplot(aes(x = month, y = mean, fill = age)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", color = "age") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

tet %>% 
  ggplot(aes(x = month, y = mean, fill = age)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", color = "age") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

trim %>% 
  ggplot(aes(x = month, y = mean, fill = age)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", color = "age") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

vanc %>% 
  ggplot(aes(x = month, y = mean, fill = age)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .2,
                position = position_dodge(width = 0.6)) +
  labs(x = "month", y = "relative abundance", color = "age") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  facet_wrap(~gene, scales = "free", ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
