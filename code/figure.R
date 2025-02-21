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
library(scales)

folder_path = "/home/h/hw309/projects/Research/TCR and Gene Expression/PC distance analysis"
file_list = list.files(folder_path)
file_GSE = file_list[str_detect(file_list, "GSE|EMTAB9357_pbmc|PRJCA002413_pbmc|10X|Zhang_2020")]

sample_information_file = paste0(folder_path, "/", file_GSE,"/sample_information_", file_GSE, ".csv")
sample_information = sample_information_file %>% map_dfr(read.csv)

df_immu_file = paste0(folder_path, "/", file_GSE, "/df_immu_", file_GSE, ".rds")
df_immu_all = df_immu_file %>% map_dfr(readRDS)

# criterion: subtype percentage
cell_count = function(ID, sample_ID){
  total_num = df_immu_all %>% 
    dplyr::count(sample, !!sym(ID), !!sym(sample_ID), name="total_num")
  
  df_immu_all = df_immu_all %>% merge(total_num)
  df_immu_all = df_immu_all %>%
    mutate(category = case_when(
      total_num == 1 ~ "Single",
      total_num >= 2  & total_num < 11  ~ "Small",
      total_num >= 11 & total_num < 100 ~ "Medium",
      total_num >= 100               ~ "Large"
    ))
  df_immu_all$category = factor(df_immu_all$category, 
                                levels = c("Single", "Small", "Medium", "Large"))
  count_df = df_immu_all %>%
    dplyr::count(disease, subtype, category)
  
  return(count_df)
}

count_df_TRA = cell_count("ID_TRA", "sample_ID_TRA") %>% mutate(clone_type = "TRA")
count_df_TRB = cell_count("ID_TRB", "sample_ID_TRB") %>% mutate(clone_type = "TRB")
count_df = rbind(count_df_TRA, count_df_TRB)

all_diseases = unique(sample_information$disease) %>% sort()

count_pattern = ggplot(data = count_df) +
  geom_bar(aes(x=subtype, y=n, fill=disease), stat = "identity", position = "dodge", alpha=0.5) +
  facet_grid(clone_type ~ category, scales = "free_x") +
  scale_fill_d3(palette = "category20", name = "Disease", limits = all_diseases) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.spacing = unit(0.05, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        strip.background =element_rect(fill="white", colour = "black", linewidth = rel(2)),
        panel.border = element_rect(colour = "black")) +
  labs(x="Subtype", y="Count")

# figure 1a
pdf(paste0(folder_path, "/plots_paper/figure1/count_pattern.pdf"), width = 13, height = 5)
count_pattern
dev.off()


# figure 1b

purity_perm_TRA = readRDS(paste0(folder_path, "/data_allStudy/aa/purity_perm_TRA.rds")) %>%
  mutate(clone_type = "TRA") %>% select(sample, category, max_prop, type, clone_type)
purity_perm_TRB = readRDS(paste0(folder_path, "/data_allStudy/aa/purity_perm_TRB.rds")) %>%
  mutate(clone_type = "TRB") %>% select(sample, category, max_prop, type, clone_type)

purity_est_TRA = readRDS(paste0(folder_path, "/data_allStudy/aa/purity_TRA.rds")) %>%
  mutate(clone_type = "TRA") %>% select(sample, category, max_prop, type, clone_type)
purity_est_TRB = readRDS(paste0(folder_path, "/data_allStudy/aa/purity_TRB.rds")) %>%
  mutate(clone_type = "TRB") %>% select(sample, category, max_prop, type, clone_type)

purity = rbind(purity_perm_TRA, purity_est_TRA,
                purity_perm_TRB, purity_est_TRB)

purity_est <- readRDS(paste0(folder_path, "/data_allStudy/aa/purity.rds"))
purity_est %>% filter(clone_subtype != "Mixture") %>% 
  count(total_num < 11) %>% mutate(n/nrow(purity))
sum(purity_est$max_prop > 0.9)/nrow(purity_est)
sum(purity_est$max_prop > 0.95)/nrow(purity_est)
sum(purity_est$max_prop > 0.99)/nrow(purity_est)
sum(purity_est$max_prop == 1)/nrow(purity_est)

figure1b = ggplot(data = purity) +
  geom_boxplot(aes(x=factor(category, levels = c("Small", "Medium", "Large")), 
                   y=max_prop, fill=type), outlier.shape = NA) +
  facet_grid(clone_type~.) +
  theme_bw() +
  labs(x="Clone size", y="Purity") +
  scale_fill_manual(#values = c('chocolate3','goldenrod'),
    values = alpha(c('pink', 'skyblue'), alpha = 0.5),
    name = "",
    labels = c('real data', 'permuted data'), 
    breaks = c('Estimated', 'Permuted')) +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.spacing = unit(0.05, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        strip.background =element_rect(fill="white", colour = "black", linewidth = rel(2)),
        panel.border = element_rect(colour = "black"))

# figure1b
pdf(paste0(folder_path, "/plots_paper/figure1/clone_perc_permuted.pdf"), width = 7, height = 3)
figure1b
dev.off()


# within same clones

figure1c_df = function(CD4_file, CD8_file, clone_type){
  distance_full_CD4 = CD4_file %>% map_dfr(readRDS) %>% mutate(clone_subtype = "CD4+CD8-")
  distance_full_CD8 = CD8_file %>% map_dfr(readRDS) %>% mutate(clone_subtype = "CD4-CD8+")
  
  distance_full = rbind(distance_full_CD4, distance_full_CD8) %>% merge(sample_information)
  
  distance_full = distance_full %>%
    mutate(category = case_when(
      total_num >= 2  & total_num < 11  ~ "Small",
      total_num >= 11 & total_num < 100 ~ "Medium",
      total_num >= 100               ~ "Large"
    ))
  
  distance_full$type = factor(distance_full$type)
  distance_full$category = factor(distance_full$category, levels = c("Small", "Medium", "Large"))

  # Median for each sample
  df_sam = 
    distance_full %>%
    group_by(sample, type, clone_subtype, disease, category) %>%
    dplyr::summarise(median = median(PC_dis),
                     .groups = "drop")
  # Median of median 
  df_disease = 
    df_sam %>%
    group_by(type, clone_subtype, disease, category) %>%
    dplyr::summarise(median = median(median),
                     .groups = "drop") %>%
    mutate(clone_type = clone_type)
  return(df_disease)
}

distance_full_CD4_file = paste0(folder_path, "/", file_GSE, "/aa/output_CD4/TRA/distance_full.rds")
distance_full_CD8_file = paste0(folder_path, "/", file_GSE, "/aa/output_CD8/TRA/distance_full.rds")
permuted_df_disease_TRA = figure1c_df(distance_full_CD4_file, distance_full_CD8_file, "TRA")

distance_full_CD4_file = paste0(folder_path, "/", file_GSE, "/aa/output_CD4/TRB/distance_full.rds")
distance_full_CD8_file = paste0(folder_path, "/", file_GSE, "/aa/output_CD8/TRB/distance_full.rds")
permuted_df_disease_TRB = figure1c_df(distance_full_CD4_file, distance_full_CD8_file, "TRB")

distance_full_CD4_file = paste0(folder_path, "/", file_GSE, "/aa/output_CD4/TRA+TRB/distance_full.rds")
distance_full_CD8_file = paste0(folder_path, "/", file_GSE, "/aa/output_CD8/TRA+TRB/distance_full.rds")
permuted_df_disease = figure1c_df(distance_full_CD4_file, distance_full_CD8_file, "TRA+TRB")


permuted_df_disease = rbind(permuted_df_disease_TRA, permuted_df_disease_TRB)
permuted_df_disease$clone_subtype = factor(permuted_df_disease$clone_subtype,
                                           levels = c("CD4+CD8-", "CD4-CD8+"))
all_diseases <- unique(sample_information$disease) %>% sort()

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
  facet_grid(clone_type~clone_subtype) +
  labs(x = "# of cells", y = "GEX dissimilarity") +
  scale_color_d3(palette = "category20", name = "Disease", limits = all_diseases) +
  scale_fill_manual(#values = c('chocolate3','goldenrod'), 
                    values = alpha(c('pink', 'skyblue'), alpha = 0.5),
                    name = "",
                    labels = c('real data', 'permuted data'), 
                    breaks = c('estimated', 'permuted')) +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.spacing = unit(0.05, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        strip.background =element_rect(fill="white", colour = "black", linewidth = rel(2)),
        panel.border = element_rect(colour = "black"))

pdf(paste0(folder_path, "/plots_paper/figure1/figure1c.pdf"), width = 10, height = 6)
figure1c
dev.off()


permuted_df_facet = 
  distance_full %>%
  group_by(sample, type, clone_type, disease) %>%
  dplyr::summarise(median = median(PC_dis),
                   .groups = "drop")
permuted_df_facet$clone_type = factor(permuted_df_facet$clone_type, levels = c("CD4+CD8-", "CD4-CD8+"))

# supp1
permuted_df_facet = permuted_df_facet %>% mutate(type=ifelse(type=="Permuted", "Permuted", "Real"))
pdf(paste0(folder_path, "/plots_paper/figure1/clones_permuted_facet.pdf"), width = 15, height = 18)
ggpaired(permuted_df_facet, x = "type", y = "median", id = "sample",
         line.color = "grey", point.size = 0.1, line.size = 0.001, 
         xlab = "Type", ylab = "GEX dissimilarity", repel = T) +
  stat_compare_means(paired = TRUE, method = "t.test") +
  facet_wrap(c("clone_type", "disease"), nrow = 5, ncol=5) +
  scale_y_continuous(limits = c(0,30)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        panel.border = element_rect(colour = "black", linewidth = 0.8))
dev.off()


# between clones

bet_clones_df = function(clone_type){
  med_dist_bet_clones_CD4 = paste0(folder_path, "/", file_GSE, "/aa/output_CD4/", clone_type, "/med_dist_bet_clones.rds") %>% map_dfr(readRDS)
  med_dist_bet_clones_CD8 = paste0(folder_path, "/", file_GSE, "/aa/output_CD8/", clone_type, "/med_dist_bet_clones.rds") %>% map_dfr(readRDS)
  
  sample_median_CD4_file = paste0(folder_path, "/", file_GSE, "/aa/output_CD4/", clone_type, "/sample_median.rds")
  sample_median_CD8_file = paste0(folder_path, "/", file_GSE, "/aa/output_CD8/", clone_type, "/sample_median.rds")
  sample_median_CD4 = 
    sample_median_CD4_file %>% map_dfr(readRDS) %>% merge(med_dist_bet_clones_CD4) %>%
    mutate(median_norm = median/median_distance_two_clones) %>% mutate(clone_subtype = "CD4+CD8-")
  sample_median_CD8 = 
    sample_median_CD8_file %>% map_dfr(readRDS) %>% merge(med_dist_bet_clones_CD8) %>%
    mutate(median_norm = median/median_distance_two_clones) %>% mutate(clone_subtype = "CD4-CD8+")
  sample_median = rbind(sample_median_CD4, sample_median_CD8) %>% merge(sample_information)
  
  Ran_sample_median_CD4_file = paste0(folder_path, "/", file_GSE, "/aa/output_CD4/", clone_type, "/Ran_sample_median.rds")
  Ran_sample_median_CD8_file = paste0(folder_path, "/", file_GSE, "/aa/output_CD8/", clone_type, "/Ran_sample_median.rds")
  Ran_sample_median_CD4 = 
    Ran_sample_median_CD4_file %>% map_dfr(readRDS) %>% merge(med_dist_bet_clones_CD4) %>%
    mutate(median_norm = median/median_distance_two_clones) %>% mutate(clone_subtype = "CD4+CD8-")
  Ran_sample_median_CD8 = 
    Ran_sample_median_CD8_file %>% map_dfr(readRDS) %>% merge(med_dist_bet_clones_CD8) %>%
    mutate(median_norm = median/median_distance_two_clones) %>% mutate(clone_subtype = "CD4-CD8+")
  Ran_sample_median = rbind(Ran_sample_median_CD4, Ran_sample_median_CD8) %>% merge(sample_information)
  sample_median = rbind(sample_median %>% mutate(type = "Estimated"),
                        Ran_sample_median %>% mutate(type ="Random"))
  return(sample_median)
}

sample_median = rbind(bet_clones_df(clone_type = "TRA") %>% mutate(clone_type = "TRA"),
                      bet_clones_df(clone_type = "TRB") %>% mutate(clone_type = "TRB"))


# supp2
pdf(paste0(folder_path, "/plots_paper/figure1/trend.pdf"), width = 14, height = 12)
  ggplot(data = sample_median %>% filter(distance <= 10)) +
  geom_boxplot(aes(x=factor(distance, levels = sort(unique(distance))), 
                   y=median, fill=type), outlier.shape = NA) +
  facet_wrap(~clone_subtype+disease, ncol=5) +
  theme_bw() +
  labs(x = "Levenshtein distance (TRB)", y = "GEX dissimilarity", 
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
  group_by(disease, distance, clone_subtype, clone_type, type) %>%
  summarise(median = median(median_norm),
            .groups = "drop")

test_results <- trend_df %>%
  filter(distance <= 10) %>%
  arrange(desc(disease)) %>%
  group_by(clone_subtype, clone_type, distance) %>%
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


# Convert distances to factors with the same levels as in the plot
test_results$distance <- factor(test_results$distance, levels = sort(unique(trend_df$distance)))
all_diseases <- unique(trend_df$disease) %>% sort()

# Create a data frame for annotations
annotations_df <- test_results %>%
  mutate(y = 1.25,  # Adjusted y_position for text
         yend = 1.3,  # Adjusted y_position for line end
         x = as.numeric(distance),  # Convert factor to numeric for positioning
         xmin = x - 0.2,
         xmax = x + 0.2)

library(ggsignif)
figure1e = ggplot(data = trend_df %>% filter(distance <= 10)) +
  geom_boxplot(aes(x = factor(distance), 
                   y = median,
                   fill=type),
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
  facet_grid(clone_type ~ factor(clone_subtype, levels = c("CD4+CD8-", "CD4-CD8+"))) +
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
  scale_y_continuous(limits = c(0.5,1.45)) +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.box.spacing = unit(0.05, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        strip.background =element_rect(fill="white", colour = "black", linewidth = rel(2)),
        panel.border = element_rect(colour = "black"))
# new
pdf(paste0(folder_path, "/plots_paper/figure1/figure1e.pdf"), width = 12, height =6)
figure1e
dev.off()


# supptable 1
Supp_table1 = df_immu_all %>%
  select(sample, GEO, gender, tissue, disease) %>%
  distinct()
write.csv(Supp_table1, 
          "~/projects/Research/TCR and Gene Expression/PC distance analysis/Supplementary Table S1.csv")





