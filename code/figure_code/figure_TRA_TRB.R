library(stringr)
library(purrr)
library(dplyr)
library(viridis)
library(ggsci)
library(ggplot2)
library(ggpattern)
library(ggpubr)
library(parallel)
library(ggpval)

folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"
file_list = list.files(folder_path)
file_GSE = file_list[str_detect(file_list, "GSE|EMTAB9357_pbmc|PRJCA002413_pbmc|10X|Zhang_2020")]

sample_information_file = paste0(folder_path, "/", file_GSE,"/sample_information_", file_GSE, ".csv")
sample_information = sample_information_file %>% map_dfr(read.csv)

df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")
df_immu_all = df_immu_file %>% map_dfr(readRDS)

# between clones
median_distance_two_clones_CD4_file = paste0(folder_path, "/", file_GSE, "/output_CD4/median_distance_two_clones_CD4_TRA_TRB.rds")
median_distance_two_clones_CD4 = median_distance_two_clones_CD4_file %>% map_dfr(readRDS)
median_distance_two_clones_CD8_file = paste0(folder_path, "/", file_GSE, "/output_CD8/median_distance_two_clones_CD8_TRA_TRB.rds")
median_distance_two_clones_CD8 = median_distance_two_clones_CD8_file %>% map_dfr(readRDS)

sample_median_CD4_file = paste0(folder_path, "/", file_GSE, "/output_CD4/sample_median_TRA_TRB.rds")
sample_median_CD8_file = paste0(folder_path, "/", file_GSE, "/output_CD8/sample_median_TRA_TRB.rds")
sample_median_CD4 = 
  sample_median_CD4_file %>% map_dfr(readRDS) %>% 
  merge(median_distance_two_clones_CD4) %>%
  mutate(median_norm = median/median_distance_two_clones) %>%
  mutate(clone_type = "CD4+CD8-")
sample_median_CD8 = 
  sample_median_CD8_file %>% map_dfr(readRDS) %>% 
  merge(median_distance_two_clones_CD8) %>%
  mutate(median_norm = median/median_distance_two_clones) %>%
  mutate(clone_type = "CD4-CD8+")

sample_median = rbind(sample_median_CD4, sample_median_CD8) %>% merge(sample_information)

Ran_sample_median_CD4_file = paste0(folder_path, "/", file_GSE, "/output_CD4/Ran_sample_median_TRA_TRB.rds")
Ran_sample_median_CD8_file = paste0(folder_path, "/", file_GSE, "/output_CD8/Ran_sample_median_TRA_TRB.rds")
Ran_sample_median_CD4 = 
  Ran_sample_median_CD4_file %>% map_dfr(readRDS) %>% 
  merge(median_distance_two_clones_CD4) %>%
  mutate(median_norm = median/median_distance_two_clones) %>%
  mutate(clone_type = "CD4+CD8-")
Ran_sample_median_CD8 = 
  Ran_sample_median_CD8_file %>% map_dfr(readRDS) %>% 
  merge(median_distance_two_clones_CD8) %>%
  mutate(median_norm = median/median_distance_two_clones) %>%
  mutate(clone_type = "CD4-CD8+")
Ran_sample_median = rbind(Ran_sample_median_CD4, Ran_sample_median_CD8) %>% merge(sample_information)
sample_median = rbind(sample_median %>% mutate(type = "Estimated"),
                      Ran_sample_median %>% mutate(type ="Random"))

# supp2
pdf(paste0(folder_path, "/plots_paper/figure1_TRA_TRB/trend.pdf"), width = 20, height = 12)
ggplot(data = sample_median %>% filter(distance_TRA_TRB <= 20)) +
  geom_boxplot(aes(x=factor(distance_TRA_TRB, levels = sort(unique(distance_TRA_TRB))), 
                   y=median, fill=type), outlier.shape = NA) +
  facet_wrap(~clone_type+disease, ncol=5) +
  theme_bw() +
  labs(x = "Levenshtein distance (TRA+TRB)", y = "GEX dissimilarity", 
       color="") +
  scale_y_continuous(limits = c(0, 25)) +
  scale_fill_manual(values = alpha(c('pink', 'skyblue'), alpha = 0.5), 
                    name = "",
                    labels = c('real data', 'permuted data'), 
                    breaks = c('Estimated', 'Random')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        panel.border = element_rect(colour = "black", linewidth = 0.8))
dev.off()

trend_df = 
  sample_median %>%
  group_by(disease, distance_TRA_TRB, clone_type, type) %>%
  summarise(median = median(median_norm),
            .groups = "drop")

test_results <- trend_df %>%
  filter(distance_TRA_TRB <= 20) %>%
  arrange(desc(disease)) %>%
  group_by(clone_type, distance_TRA_TRB) %>%
  summarize(p.value = wilcox.test(median[type == 'Estimated'], 
                                  median[type == 'Random'],
                                  alternative = "less")$p.value) %>%
  ungroup() %>%
  mutate(signif = case_when(
    p.value < 0.001 ~ '***',
    p.value < 0.01 ~ '**',
    p.value < 0.05 ~ '*',
    TRUE ~ 'ns'
  ))

# Benjamini-Hochberg (BH) multiple testing (try, not used finally)
#bh = BH(test_results %>% filter(clone_type == "CD4+CD8-") %>% pull(p.value))
#bh = BH(test_results$p.value)
# Convert distances to factors with the same levels as in the plot
test_results$distance_TRA_TRB <- factor(test_results$distance_TRA_TRB, levels = sort(unique(trend_df$distance_TRA_TRB)))
all_diseases <- unique(trend_df$disease) %>% sort()

# Create a data frame for annotations
annotations_df <- test_results %>%
  mutate(y = 1.6,  # Adjusted y_position for text
         yend = 1.65,  # Adjusted y_position for line end
         x = as.numeric(distance_TRA_TRB),  # Convert factor to numeric for positioning
         xmin = x - 0.2,
         xmax = x + 0.2)


library(ggsignif)
# new
pdf(paste0(folder_path, "/plots_paper/figure1_TRA_TRB/figure1d.pdf"), width = 8, height =4.7)
ggplot(data = trend_df %>% filter(distance_TRA_TRB <= 20)) +
  geom_boxplot(aes(x = factor(distance_TRA_TRB), 
                   y = median,
                   fill=type), alpha = 0.8,
               outlier.shape = NA,
               position = position_dodge2(width = 0.75)) +
  geom_point(data = trend_df %>% filter(type == "Estimated" &
                                          distance_TRA_TRB <= 20),
             aes(x = factor(distance_TRA_TRB),
                 y = median, 
                 col = disease,
                 group = type),
             alpha = 0.5,
             position = position_nudge(x = -0.2)) +
  facet_grid(clone_type~.) +
  geom_text(data = annotations_df, aes(x = x, y = y, label = signif), vjust = -1) +
  geom_segment(data = annotations_df, aes(x = xmin, xend = xmax, y = yend, yend = yend)) +
  theme_bw() +
  labs(x = "Levenshtein distance (TRA+TRB)", y = "GEX dissimilarity", 
       fill="") +
  scale_color_d3(palette = "category20", name = "Disease", limits = all_diseases) +
  scale_fill_manual(values = alpha(c('pink', 'skyblue'), alpha = 0.5), 
                    name = "",
                    labels = c('real data', 'permuted data'), 
                    breaks = c('Estimated', 'Random')) +
  scale_y_continuous(limits = c(0.5, 1.9)) +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.spacing = unit(0, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        strip.background =element_rect(fill="white"),
        panel.border = element_rect(colour = "black", linewidth = 0.8))
dev.off()


