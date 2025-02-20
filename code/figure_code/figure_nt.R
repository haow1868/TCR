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

# criterion: subtype percentage

total_num = df_immu_all %>% 
  dplyr::count(sample, ID_nt, sample_ID_nt, name="total_num")

df_immu_all = df_immu_all %>% merge(total_num)

df_immu_all = df_immu_all %>%
  mutate(category = case_when(
    total_num == 1 ~ "Single",
    total_num >= 2  & total_num < 11  ~ "Small",
    total_num >= 11 & total_num < 100 ~ "Medium",
    total_num >= 100               ~ "Large"
  ))
df_immu_all$category = factor(df_immu_all$category, 
                              levels = c("Large", "Medium", "Small", "Single"))
count_pattern = ggplot(data = df_immu_all, 
                       aes(x=disease, fill=subtype, pattern=category)) +
  geom_bar_pattern(
    color = "black", 
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.025,
    pattern_key_scale_factor = 0.6
  ) +
  coord_flip() +
  theme_bw() +
  scale_fill_npg() +
  labs(x = "", y = "", pattern = "Category", fill = "Subtype") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(axis.title = element_text(size = 12),
        #  axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.8)) +
  scale_pattern_manual(values = c("stripe", "crosshatch", "circle", "none"))

# figure 1a
pdf(paste0(folder_path, "/plots_paper/figure1_nt/count_pattern.pdf"), width = 7, height = 5)
count_pattern
dev.off()


df_immu_all_per = readRDS(paste0(folder_path, "/data_allStudy/df_immu_all_per_nt.rds"))

clone_prop = function(df){
  sub_num = df %>% 
    dplyr::count(sample, ID_nt, sample_ID_nt, subtype) %>%
    merge(total_num) %>%
    filter(total_num > 1) %>%
    mutate(per = n/total_num) %>%
    mutate(diff = abs(per-0.5)) 
  
  sub_num = sub_num %>%
    dplyr::mutate(category = case_when(
      total_num == 1 ~ "Single",
      total_num >= 2  & total_num < 11  ~ "Small",
      total_num >= 11 & total_num < 100 ~ "Medium",
      total_num >= 100               ~ "Large"
    ))
  
  sub_num = sub_num %>% 
    dplyr::mutate(clone_subtype = ifelse(per == 1, "Same celltype", "Mixed celltype")) %>%
    dplyr::mutate(clone_subtype = case_when(
      clone_subtype == "Same celltype" ~ subtype,
      clone_subtype == "Mixed celltype" ~ "Mixture",
    )) 
  
  sub_num = sub_num %>%
    group_by(sample_ID_nt) %>%
    mutate(max_prop = max(per)) %>%
    distinct(sample, ID_nt, sample_ID_nt, category, clone_subtype, total_num, max_prop) %>%
    ungroup()
  return(sub_num)
}

sub_num_perm_list = mclapply(df_immu_all_per, clone_prop, mc.cores = 23)
sub_num_perm = sub_num_perm_list %>% bind_rows()
sub_num_perm = sub_num_perm %>%
  group_by(sample, ID_nt, sample_ID_nt, category, total_num) %>%
  summarise(max_prop = mean(max_prop)) %>%
  mutate(type = "Permuted")
sub_num_est = clone_prop(df_immu_all) %>% mutate(type = "Estimated")
sub_num = rbind(sub_num_perm, sub_num_est)
saveRDS(sub_num, paste0(folder_path, "/data_allStudy/sub_num_nt.rds"))

sub_num <- readRDS(paste0(folder_path,"/data_allStudy/sub_num_nt.rds"))
sub_num_est = sub_num %>% filter(type == "Estimated")
sum(sub_num_est$max_prop > 0.9)/nrow(sub_num_est)
sum(sub_num_est$max_prop > 0.95)/nrow(sub_num_est)
sum(sub_num_est$max_prop > 0.99)/nrow(sub_num_est)
sum(sub_num_est$max_prop == 1)/nrow(sub_num_est)


# figure1b
pdf(paste0(folder_path, "/plots_paper/figure1_nt/clone_perc_permuted.pdf"), width = 5, height = 4)
ggplot(data = sub_num) +
  geom_boxplot(aes(x=factor(category, levels = c("Small", "Medium", "Large")), 
                   y=max_prop, fill=type), outlier.shape = NA) +
  theme_bw() +
  labs(x="clone size", y="purity") +
  scale_fill_manual(#values = c('chocolate3','goldenrod'),
                    values = alpha(c('pink', 'skyblue'), alpha = 0.5),
                    name = "",
                    labels = c('real data', 'permuted data'), 
                    breaks = c('Estimated', 'Permuted')) +
  theme(panel.border = element_rect(colour = "black", linewidth = 0.8))
dev.off()


sub_num_mix %>% ungroup() %>% count(total_num < 11) %>% mutate(n/nrow(sub_num_mix))
sub_num_mix %>% filter(total_num >= 20) %>% count(total_num >= 20 & max_prop > 0.7)

sub_num_mix %>% count(max_prop <= 0.9) %>% mutate(n/sum(n))


# within same clones
distance_full_CD4_file = paste0(folder_path, "/", file_GSE, "/DNA/output_CD4/distance_full.rds")
distance_full_CD8_file = paste0(folder_path, "/", file_GSE, "/DNA/output_CD8/distance_full.rds")
distance_full_CD4 = distance_full_CD4_file %>% map_dfr(readRDS) %>% mutate(clone_type = "CD4+CD8-")
distance_full_CD8 = distance_full_CD8_file %>% map_dfr(readRDS) %>% mutate(clone_type = "CD4-CD8+")

distance_full = rbind(distance_full_CD4, distance_full_CD8) %>% merge(sample_information)

distance_full = distance_full %>%
  mutate(category = case_when(
    total_num >= 2  & total_num < 11  ~ "Small",
    total_num >= 11 & total_num < 100 ~ "Medium",
    total_num >= 100               ~ "Large"
  ))

distance_full = distance_full %>%
  filter(type %in% c("estimated", "permuted_within"))
distance_full$type = if_else(distance_full$type == "estimated", "Estimated", "Permuted")
distance_full$type = factor(distance_full$type)
distance_full$category = factor(distance_full$category, levels = c("Small", "Medium", "Large"))

permuted_df = 
  distance_full %>%
  filter(type %in% c("Estimated", "Permuted")) %>%
  group_by(sample, ID_nt, sample_ID_nt, type) %>%
  mutate(PC_dis = mean(PC_dis)) %>% #average for 100 times permuted distance
  distinct() %>%
  ungroup() %>%
  group_by(disease, type, clone_type, category) %>%
  summarise(median = median(PC_dis),
            .groups = "drop")
all_diseases <- unique(permuted_df$disease) %>% sort()

permuted_df_sam = 
  distance_full %>%
  filter(type %in% c("Estimated", "Permuted")) %>%
  group_by(sample, ID_nt, sample_ID_nt, type) %>%
  mutate(PC_dis = mean(PC_dis)) %>% #average for 100 times permuted distance
  distinct() %>%
  ungroup() %>%
  group_by(sample, type, clone_type, disease, category) %>%
  dplyr::summarise(median = median(PC_dis),
                   .groups = "drop")

permuted_df_disease = 
  permuted_df_sam %>%
  group_by(type, clone_type, disease, category) %>%
  dplyr::summarise(median = median(median),
                   .groups = "drop")

# figure1c
figure1c = ggplot(data = permuted_df_disease) +
  geom_boxplot(aes(x = factor(category),
                   y = median, 
                   fill = type),
               outlier.shape = NA, position = position_dodge2(width = 0.75)) +
  geom_jitter(aes(x = factor(category),
                  y = median, 
                  col = disease,
                  group = type),
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75),
              alpha=0.5) +
  theme_bw() +
  facet_grid(clone_type ~.) +
  labs(x = "# of cells", y = "GEX dissimilarity") +
  scale_color_d3(palette = "category20", name = "Disease", limits = all_diseases) +
  scale_fill_manual(#values = c('chocolate3','goldenrod'),
                    values = alpha(c('pink', 'skyblue'), alpha = 0.5),
                    name = "",
                    labels = c('real data', 'permuted data'), 
                    breaks = c('Estimated', 'Permuted')) +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.spacing = unit(0, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        strip.background =element_rect(fill="white"),
        panel.border = element_rect(colour = "black", linewidth = 0.8))


pdf(paste0(folder_path, "/plots_paper/figure1_nt/figure1c.pdf"), width = 6, height = 4.6)
figure1c
dev.off()


permuted_df_facet = 
  distance_full %>%
  filter(type %in% c("Estimated", "Permuted")) %>%
  group_by(sample, ID_nt, sample_ID_nt, type) %>%
  mutate(PC_dis = mean(PC_dis)) %>% #average for 100 times permuted distance
  distinct() %>%
  ungroup() %>%
  group_by(sample, type, clone_type, disease) %>%
  dplyr::summarise(median = median(PC_dis),
                   .groups = "drop")
permuted_df_facet$clone_type = factor(permuted_df_facet$clone_type, levels = c("CD4+CD8-", "CD4-CD8+"))

remove_outliers <- function(df, y) {
  Q1 <- quantile(df[[y]], 0.25)
  Q3 <- quantile(df[[y]], 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  df %>% 
    filter(df[[y]] >= lower_bound & df[[y]] <= upper_bound)
}

filtered_ids <- permuted_df_facet %>%
  filter(type == "Estimated") %>%
  group_by(clone_type, disease) %>%
  do(remove_outliers(., "median")) %>%
  pull(sample)

# Filter the original dataset to include only the pairs that are left
final_df <- permuted_df_facet %>%
  filter(sample %in% filtered_ids)

# supp1
permuted_df_facet = permuted_df_facet %>% mutate(type=ifelse(type=="Permuted", "Permuted", "Real"))
pdf(paste0(folder_path, "/plots_paper/figure1_nt/clones_permuted_facet.pdf"), width = 15, height = 18)
ggpaired(permuted_df_facet, x = "type", y = "median", id = "sample",
         line.color = "grey", point.size = 0.1, line.size = 0.001, 
         xlab = "Type", ylab = "GEX dissimilarity", repel = T) +
  stat_compare_means(paired = TRUE, method = "t.test") +
  facet_wrap(c("clone_type", "disease"), nrow = 5, ncol=5) +
#  scale_y_continuous(limits = c(0,35)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        panel.border = element_rect(colour = "black", linewidth = 0.8))
dev.off()


# between clones
median_distance_two_clones_CD4_file = paste0(folder_path, "/", file_GSE, "/DNA/output_CD4/median_distance_two_clones_CD4.rds")
median_distance_two_clones_CD4 = median_distance_two_clones_CD4_file %>% map_dfr(readRDS)
median_distance_two_clones_CD8_file = paste0(folder_path, "/", file_GSE, "/DNA/output_CD8/median_distance_two_clones_CD8.rds")
median_distance_two_clones_CD8 = median_distance_two_clones_CD8_file %>% map_dfr(readRDS)

sample_median_CD4_file = paste0(folder_path, "/", file_GSE, "/DNA/output_CD4/sample_median.rds")
sample_median_CD8_file = paste0(folder_path, "/", file_GSE, "/DNA/output_CD8/sample_median.rds")
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

Ran_sample_median_CD4_file = paste0(folder_path, "/", file_GSE, "/DNA/output_CD4/Ran_sample_median.rds")
Ran_sample_median_CD8_file = paste0(folder_path, "/", file_GSE, "/DNA/output_CD8/Ran_sample_median.rds")
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
pdf(paste0(folder_path, "/plots_paper/figure1_nt/trend.pdf"), width = 14, height = 12)
ggplot(data = sample_median %>% filter(distance <= 10)) +
  geom_boxplot(aes(x=factor(distance, levels = sort(unique(distance))), 
                   y=median, fill=type), outlier.shape = NA) +
  facet_wrap(~clone_type+disease, ncol=5) +
  theme_bw() +
  labs(x = "Levenshtein distance", y = "GEX dissimilarity", 
       color="") +
  scale_y_continuous(limits = c(0, 25)) +
  scale_fill_manual(#values = c('chocolate3','goldenrod'),
                    values = alpha(c('pink', 'skyblue'), alpha = 0.5),
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
  group_by(disease, distance, clone_type, type) %>%
  summarise(median = median(median_norm),
            .groups = "drop")

test_results <- trend_df %>%
  filter(distance <= 10) %>%
  arrange(desc(disease)) %>%
  group_by(clone_type, distance) %>%
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

# Benjamini-Hochberg (BH) multiple testing
bh = BH(test_results %>% filter(clone_type == "CD4+CD8-") %>% pull(p.value))
bh = BH(test_results$p.value)
# Convert distances to factors with the same levels as in the plot
test_results$distance <- factor(test_results$distance, levels = sort(unique(trend_df$distance)))

# Create a data frame for annotations
annotations_df <- test_results %>%
  mutate(y = 1.2,  # Adjusted y_position for text
         yend = 1.25,  # Adjusted y_position for line end
         x = as.numeric(distance),  # Convert factor to numeric for positioning
         xmin = x - 0.2,
         xmax = x + 0.2)


library(ggsignif)
# new
pdf(paste0(folder_path, "/plots_paper/figure1_nt/figure1d.pdf"), width = 6.5, height =4.7)
ggplot(data = trend_df %>% filter(distance <= 10)) +
  geom_boxplot(aes(x = factor(distance), 
                   y = median,
                   fill=type), alpha = 0.8,
               outlier.shape = NA,
               position = position_dodge2(width = 0.75)) +
  geom_point(data = trend_df %>% filter(type == "Estimated" &
                                          distance <= 10),
             aes(x = factor(distance),
                 y = median, 
                 col = disease,
                 group = type),
             alpha = 0.5,
             position = position_nudge(x = -0.2)) +
  facet_grid(clone_type~.) +
  geom_text(data = annotations_df, aes(x = x, y = y, label = signif), vjust = -1) +
  geom_segment(data = annotations_df, aes(x = xmin, xend = xmax, y = yend, yend = yend)) +
  theme_bw() +
  labs(x = "Levenshtein distance", y = "GEX dissimilarity", 
       fill="") +
  scale_color_d3(palette = "category20", name = "Disease", limits = all_diseases) +
  scale_fill_manual(#values = c('chocolate3','goldenrod'),
                    values = alpha(c('pink', 'skyblue'), alpha = 0.5),
                    name = "",
                    labels = c('real data', 'permuted data'), 
                    breaks = c('Estimated', 'Random')) +
  scale_y_continuous(limits = c(0.5,1.4)) +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.spacing = unit(0, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        strip.background =element_rect(fill="white"),
        panel.border = element_rect(colour = "black", linewidth = 0.8))
dev.off()


