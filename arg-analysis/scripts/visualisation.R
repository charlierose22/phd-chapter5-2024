# GRAPH ----------------------------------------------------------------
library(RColorBrewer)
library(ggsci)
library(grafify)


means$month <- ym(means$month)
means_class_genes$month <- ym(means_class_genes$month)
means_mechanism$month <- ym(means_mechanism$month)

# Create a list of target antibiotics.
class <- unique(means$class)

# split based on target antibiotics for location 
split_seasonal <- split(means, means$class)
aminoglycoside_seasonal <- split_seasonal$Aminoglycoside
beta_seasonal <- split_seasonal$'Beta-lactam'
glycopeptide_metronidazole_seasonal <- split_seasonal$'Glycopeptides and Metronidazole'
macrolide_lincosamide_seasonal <- split_seasonal$'MLSB'
other_seasonal <- split_seasonal$Other
phenicol_seasonal <- split_seasonal$Phenicol
quinolone_seasonal <- split_seasonal$Quinolone
sulfonamide_trimethoprim_seasonal <- split_seasonal$'Sulfonamides and Trimethoprim'
tetracycline_seasonal <- split_seasonal$Tetracycline
mdr_seasonal <- split_seasonal$'Multi-Drug Resistance'
mge_integrons_seasonal <- split_seasonal$'MGEs and Integrons'

# replace days as factor for significance tests
month_factor <- day29_genes
month_factor$month <- as.factor(month_factor$month)

stacked_genes <- day29_genes
stacked_genes$month <- ym(stacked_genes$month)
# total
# stacked area plot
stacked_seasonal <- stacked_genes %>% 
  group_by(month, class) %>%
  summarise(n = sum(delta_ct)) %>%
  mutate(percentage = n / sum(n))
stacked_seasonal$percentage <- stacked_seasonal$percentage * 100

stacked_mechanism <- stacked_genes %>% 
  group_by(month, mechanism) %>%
  summarise(n = sum(delta_ct)) %>%
  mutate(percentage = n / sum(n))
stacked_mechanism$percentage <- stacked_mechanism$percentage * 100

# time
ggplot(stacked_seasonal, aes(x = month,
                            y = percentage,
                            fill = class)) +
  geom_area(alpha=0.6 , linewidth=0.5, colour="black") +
  scale_fill_manual(values = brewer.pal("Spectral", n = 11)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  labs(x = "Month", y = "Proportion of Relative Gene Abundance (%)", fill = "Target Antibiotic Class") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggplot(stacked_mechanism, aes(x = month,
                             y = percentage,
                             fill = mechanism)) +
  geom_area(alpha=0.6 , linewidth=0.5, colour="black") +
  scale_fill_manual(values = brewer.pal("Spectral", n = 11)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  labs(x = "Month", y = "Proportion of Relative Gene Abundance (%)", fill = "Resistance Mechanism Type") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")


# INDIVIDUAL --------------------------------------------------------------
# total gene
# total gene over time
means_class_genes %>%
  ggplot(aes(x = month, y = mean, colour = class)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance (log10)", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  scale_y_continuous(trans='log10') +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# mechanisms
means_mechanism %>%
  ggplot(aes(x = month, y = mean, colour = mechanism)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance (log10)", color = "Resistance Mechanism Type") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  scale_y_continuous(trans='log10') +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Aminoglycosides
aminoglycoside_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Beta-lactam
beta_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Glycopeptide and Metronidazole
glycopeptide_metronidazole_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Macrolide and Lincosamide
macrolide_lincosamide_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Other
other_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Phenicol
phenicol_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Quinolone
quinolone_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Sulfonamide and Trimethoprim
sulfonamide_trimethoprim_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Tetracycline
tetracycline_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# MDR
mdr_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# MGE and Integrons
mge_integrons_seasonal %>%
  ggplot(aes(x = month, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = month,
                    ymin = mean - se,
                    ymax = mean + se), width = 5) +
  geom_line() +
  labs(x = "Month", y = "Average Relative Gene Abundance", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly") +
  scale_x_date(date_breaks = "1 month", date_labels = "%y-%b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# statistics
# replace days as factor for significance tests
month_factor <- day29_genes
month_factor$month <- as.factor(month_factor$month)
month_factor$name.x <- as.factor(month_factor$gene)
month_factor$class <- as.factor(month_factor$class)

split_month_time_factor <- split(month_factor, month_factor$class)
month_factor_amino <- split_month_time_factor$Aminoglycoside
month_factor_beta <- split_month_time_factor$'Beta-lactam'
month_factor_glyc <- split_month_time_factor$'Glycopeptides and Metronidazole'
month_factor_mac <- split_month_time_factor$'MLSB'
month_factor_other <- split_month_time_factor$Other
month_factor_phen <- split_month_time_factor$Phenicol
month_factor_quin <- split_month_time_factor$Quinolone
month_factor_sulf <- split_month_time_factor$'Sulfonamides and Trimethoprim'
month_factor_tet <- split_month_time_factor$Tetracycline
month_factor_mdr <- split_month_time_factor$'Multi-Drug Resistance'
month_factor_mge <- split_month_time_factor$'MGEs and Integrons'

# time study
library(FSA)
# all classes
# one way anova
one_way_all <- aov(delta_ct ~ class, data = month_factor)
summary(one_way_all)
# two way anova
two_way_all <- aov(delta_ct ~ month * class, data = month_factor)
summary(two_way_all)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_all)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor, delta_ct ~ class)
kruskal.test(data = month_factor, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ class,
         data = month_factor,
         method = "bh")
dunnTest(delta_ct ~ month,
         data = month_factor,
         method = "bh")

ggplot(month_factor, aes(x = class, y = delta_ct) ) +
  geom_boxplot() +
  theme_classic()

# one way anova
one_way_mech <- aov(delta_ct ~ mechanism, data = month_factor)
summary(one_way_mech)
# two way anova
two_way_mech <- aov(delta_ct ~ month * class, data = month_factor)
summary(two_way_mech)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_mech)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor, delta_ct ~ mechanism)
kruskal.test(data = month_factor, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ mechanism,
         data = month_factor,
         method = "bh")
dunnTest(delta_ct ~ month,
         data = month_factor,
         method = "bh")

ggplot(month_factor, aes(x = class, y = delta_ct) ) +
  geom_boxplot() +
  theme_classic()

# aminoglycoside
# one way anova
one_way_amino <- aov(delta_ct ~ gene, data = month_factor_amino)
summary(one_way_amino)
one_way_amino_month <- aov(delta_ct ~ month, data = month_factor_amino)
summary(one_way_amino_month)
# two way anova
two_way_amino <- aov(delta_ct ~ month * gene, data = month_factor_amino)
summary(two_way_amino)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_amino)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor_amino, delta_ct ~ gene)
kruskal.test(data = month_factor_amino, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ month,
         data = month_factor_amino,
         method = "bh")
dunnTest(delta_ct ~ gene,
         data = month_factor_amino,
         method = "bh")

ggplot(month_factor_amino, aes(x = gene, y = delta_ct) ) +
  geom_boxplot() +
  theme_classic()

# beta lactam
# one way anova
one_way_beta <- aov(delta_ct ~ gene, data = month_factor_beta)
summary(one_way_beta)
one_way_beta_month <- aov(delta_ct ~ month, data = month_factor_beta)
summary(one_way_beta_month)
# two way anova
two_way_beta <- aov(delta_ct ~ month * gene, data = month_factor_beta)
summary(two_way_beta)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_beta)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor_beta, delta_ct ~ gene)
kruskal.test(data = month_factor_beta, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ month,
         data = month_factor_beta,
         method = "bh")
dunnTest(delta_ct ~ gene,
         data = month_factor_beta,
         method = "bh")
ggplot(month_factor_beta, aes(x = gene, y = delta_ct) ) +
  geom_boxplot() +
  theme_classic()

# glycopeptide and metronidazole
# one way anova
one_way_glyc <- aov(delta_ct ~ gene, data = month_factor_glyc)
summary(one_way_glyc)
one_way_glyc_month <- aov(delta_ct ~ month, data = month_factor_glyc)
summary(one_way_glyc_month)
# two way anova
two_way_glyc <- aov(delta_ct ~ month * gene, data = month_factor_glyc)
summary(two_way_glyc)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_glyc)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor_glyc, delta_ct ~ gene)
kruskal.test(data = month_factor_glyc, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ month,
         data = month_factor_glyc,
         method = "bh")
dunnTest(delta_ct ~ gene,
         data = month_factor_glyc,
         method = "bh")

# macrolide and lincosamide
# one way anova
one_way_mac <- aov(delta_ct ~ gene, data = month_factor_mac)
summary(one_way_mac)
one_way_mac_month <- aov(delta_ct ~ month, data = month_factor_mac)
summary(one_way_mac_month)
# two way anova
two_way_mac <- aov(delta_ct ~ month * gene, data = month_factor_mac)
summary(two_way_mac)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_mac)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor_mac, delta_ct ~ gene)
kruskal.test(data = month_factor_mac, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ month,
         data = month_factor_mac,
         method = "bh")
dunnTest(delta_ct ~ gene,
         data = month_factor_mac,
         method = "bh")

# other
# one way anova
one_way_other <- aov(delta_ct ~ gene, data = month_factor_other)
summary(one_way_other)
one_way_other_month <- aov(delta_ct ~ month, data = month_factor_other)
summary(one_way_other_month)
# two way anova
two_way_other <- aov(delta_ct ~ month * gene, data = month_factor_other)
summary(two_way_other)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_other)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor_other, delta_ct ~ gene)
kruskal.test(data = month_factor_other, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ month,
         data = month_factor_other,
         method = "bh")
dunnTest(delta_ct ~ gene,
         data = month_factor_other,
         method = "bh")

# phenicol
# one way anova
one_way_phen <- aov(delta_ct ~ gene, data = month_factor_phen)
summary(one_way_phen)
one_way_phen_month <- aov(delta_ct ~ month, data = month_factor_phen)
summary(one_way_phen_month)
# two way anova
two_way_phen <- aov(delta_ct ~ month * gene, data = month_factor_phen)
summary(two_way_phen)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_phen)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor_phen, delta_ct ~ gene)
kruskal.test(data = month_factor_phen, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ month,
         data = month_factor_phen,
         method = "bh")
dunnTest(delta_ct ~ gene,
         data = month_factor_phen,
         method = "bh")

# quinolone
# one way anova
one_way_quin <- aov(delta_ct ~ gene, data = month_factor_quin)
summary(one_way_quin)
one_way_quin_month <- aov(delta_ct ~ month, data = month_factor_quin)
summary(one_way_quin_month)
# two way anova
two_way_quin <- aov(delta_ct ~ month * gene, data = month_factor_quin)
summary(two_way_quin)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_quin)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor_quin, delta_ct ~ gene)
kruskal.test(data = month_factor_quin, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ month,
         data = month_factor_quin,
         method = "bh")
dunnTest(delta_ct ~ gene,
         data = month_factor_quin,
         method = "bh")
ggplot(month_factor_quin, aes(x = gene, y = delta_ct) ) +
  geom_boxplot() +
  theme_classic()

# sulfonamide and trimethoprim
# one way anova
one_way_sulf <- aov(delta_ct ~ gene, data = month_factor_sulf)
summary(one_way_sulf)
one_way_sulf_month <- aov(delta_ct ~ month, data = month_factor_sulf)
summary(one_way_sulf_month)
# two way anova
two_way_sulf <- aov(delta_ct ~ month * gene, data = month_factor_sulf)
summary(two_way_sulf)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_sulf)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor_sulf, delta_ct ~ gene)
kruskal.test(data = month_factor_sulf, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ month,
         data = month_factor_sulf,
         method = "bh")
dunnTest(delta_ct ~ gene,
         data = month_factor_sulf,
         method = "bh")
ggplot(month_factor_sulf, aes(x = gene, y = delta_ct) ) +
  geom_boxplot() +
  theme_classic()

# tetracycline
# one way anova
one_way_tet <- aov(delta_ct ~ gene, data = month_factor_tet)
summary(one_way_tet)
one_way_tet_month <- aov(delta_ct ~ month, data = month_factor_tet)
summary(one_way_tet_month)
# two way anova
two_way_tet <- aov(delta_ct ~ month * gene, data = month_factor_tet)
summary(two_way_tet)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_tet)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor_tet, delta_ct ~ gene)
kruskal.test(data = month_factor_tet, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ month,
         data = month_factor_tet,
         method = "bh")
dunnTest(delta_ct ~ gene,
         data = month_factor_tet,
         method = "bh")

# mdr
# one way anova
one_way_mdr <- aov(delta_ct ~ gene, data = month_factor_mdr)
summary(one_way_mdr)
one_way_mdr_month <- aov(delta_ct ~ month, data = month_factor_mdr)
summary(one_way_mdr_month)
# two way anova
two_way_mdr <- aov(delta_ct ~ month * gene, data = month_factor_mdr)
summary(two_way_mdr)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_mdr)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor_mdr, delta_ct ~ gene)
kruskal.test(data = month_factor_mdr, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ month,
         data = month_factor_mdr,
         method = "bh")
dunnTest(delta_ct ~ gene,
         data = month_factor_mdr,
         method = "bh")

# mge and integrons
# one way
one_way_mge_month <- aov(delta_ct ~ month, data = month_factor_mge)
summary(one_way_mge_month)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_mge_month)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = month_factor_mge, delta_ct ~ month)
# post-hoc test
dunnTest(delta_ct ~ month,
         data = month_factor_mge,
         method = "bh")
ggplot(month_factor_mge, aes(x = month, y = delta_ct) ) +
  geom_boxplot() +
  labs(x = "Month", y = "Relative Gene Abundance") +
  theme_bw(base_size = 12)
