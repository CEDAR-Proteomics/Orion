## MaxIDs

## Preliminary Analysis on Orion Plates 1-30
# Note: no sample 38 so 164 total samples

# Load packages
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(stringr)
library(pcaMethods)
library(ggrepel)
library(readxl)
library(zoo)
library(scico)


## Import data ----------------------------------------------------------------------------------------------------------
peptide_panel <- read.delim("./Peptide_Panel.tsv", sep = "\t", header = TRUE)
sample_annotation <- read_xlsx("./Prep Info and Observations Log Final.xlsx")

# create a new column that combines the parent protein name for each peptide and the associated peptide sequence
peptide_panel$peptide_ID <- paste(peptide_panel$Protein.Group, peptide_panel$Peptide, sep = "_")

# extract the sample number into a new column from the Sample.Name column
peptide_panel$SampleNum <- as.numeric(str_sub(peptide_panel$Sample.Name, 7))

blank_panel <- peptide_panel[0, ]
blank_panel[nrow(blank_panel) + 48, ] <- NA

full_sequence <- 1:450
x_num <- unique(peptide_panel$SampleNum)
missing_samples <- full_sequence[!full_sequence %in% x_num]
blank_panel$SampleNum <- missing_samples

peptide_panel <- rbind(peptide_panel, blank_panel)


# some sample formatting of the annotation file
sample_annotation <- sample_annotation[-c(1,2), ]
sample_annotation$Plate <- na.locf(sample_annotation$`Plate Name`)

sample_annotation_selected <- sample_annotation %>% 
  select(`OHSU ID`, `UT ID`, `Qualitative Observations (supplementary info)`, Color, `Solids observed?`, Plate) %>% 
  filter(substr(`OHSU ID`, 1, 1) == "S")

# extract sample number from the annotation file
sample_annotation_selected$SampleNum <- as.numeric(str_sub(sample_annotation_selected$`OHSU ID`, 2))

# combine the annotation describing the serum color and presence of solids into one column
sample_annotation_selected$Color_Plus_Precipitates_Observed <- paste(sample_annotation_selected$Color, 
                                                                     sample_annotation_selected$`Solids observed?`,
                                                                     sep = "|")

# qualifying outlier status for color and precipitate
sample_annotation_selected <- sample_annotation_selected %>% 
  mutate(Outlier_Color = ifelse(Color %in% c("slightly peach", "peach", "darker peach", "orangeish yellow", "reddish orange", "red"), "outlier", "normal"),
         Outlier_Solids = ifelse(`Solids observed?` == "abundant", "outlier", "normal")
  )

sample_annotation_selected <- sample_annotation_selected %>% 
  mutate(Outlier_Status = ifelse(Outlier_Color == "outlier" | Outlier_Solids == "outlier", "outlier", "normal"),
         Outlier_Status = factor(Outlier_Status)
  )

outlier_list_full <- sample_annotation_selected %>% 
  filter(Outlier_Status == "outlier") %>% 
  select(SampleNum)


## PCA -----------------------------------------------------------------------------------------------------------
peptides_wide <- peptide_panel %>% 
  drop_na() %>% 
  select(Sample.Name, Intensity..Log10., peptide_ID) %>% 
  pivot_wider(names_from = peptide_ID,
              values_from = Intensity..Log10.)

peptides_wide <- peptides_wide %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Sample.Name")

peptides_ctr <- prep(peptides_wide, scale = "none", center = TRUE)

peptides_pca <- pca(peptides_ctr, method = "nipals", center = FALSE)

pca_scores <- scores(peptides_pca)
pca_scores <- cbind(rownames(pca_scores), data.frame(pca_scores, row.names = NULL))

# extract the sample number into a new column from the Sample.Name column
pca_scores$SampleNum <- as.numeric(str_sub(pca_scores$`rownames(pca_scores)`, 7))

pca_scores <- merge(pca_scores, sample_annotation_selected, by = "SampleNum", all.x = TRUE)
pca_scores <- merge(pca_scores, peptide_ids_annotation, by = "SampleNum", all.x = TRUE)

pca_scores$PlateNum_old <- str_sub(pca_scores$Plate, 14, 15)
pca_scores$PlateNum_old <- str_remove(pca_scores$PlateNum_old, "_")
pca_scores$PlateNum_old[pca_scores$PlateNum_old == "B"] <- 0
pca_scores$PlateNum_old[pca_scores$PlateNum_old == "Y"] <- 0

pca_scores$PlateNum_new <- str_sub(pca_scores$Plate, 12, 13)
pca_scores$PlateNum_new <- str_remove(pca_scores$PlateNum_new, "_")
pca_scores$PlateNum_new[pca_scores$PlateNum_new == "P"] <- 0

pca_scores$PlateNum_old <- as.numeric(pca_scores$PlateNum_old)
pca_scores$PlateNum_new <- as.numeric(pca_scores$PlateNum_new)

pca_scores$PlateNum <- pca_scores$PlateNum_new + pca_scores$PlateNum_old

# pca plot colored by plate number
pca_plot <- pca_scores %>% 
  ggplot(aes(x = PC1, y = PC2, color = as.factor(PlateNum), text = paste(`OHSU ID`, Num_Peptides, Color_Plus_Precipitates_Observed, sep = ", "))) +
  #geom_point(size = 2, color = '#22c6c9') +
  geom_point(size = 2, aes(color = as.factor(PlateNum))) +
  theme_bw() +
  scale_color_scico_d(palette = 'batlow', name = "Plate Number") +
  #theme(legend.position = "none") +
  labs(title = "PCA on Peptides")
pca_plot


# pca plot colored by description of samples
pca_plot_color_ppt <- pca_scores %>% 
  ggplot(aes(x = PC1, y = PC2, color = as.factor(Color_Plus_Precipitates_Observed), text = paste(`OHSU ID`, Num_Peptides, Color_Plus_Precipitates_Observed, sep = ", "))) +
  #geom_point(size = 2, color = '#22c6c9') +
  geom_point(size = 2, aes(color = Color_Plus_Precipitates_Observed)) +
  theme_bw() +
  #theme(legend.position = "none") +
  labs(title = "PCA on Peptides")
pca_plot_color_ppt


# pca plot colored by outlier status
pca_plot_outlier <- pca_scores %>% 
  ggplot(aes(x = PC1, y = PC2, color = Outlier_Status, text = paste(`OHSU ID`, Num_Peptides, Color_Plus_Precipitates_Observed, sep = ", "))) +
  #geom_point(size = 2, color = '#22c6c9') +
  geom_point(size = 2, alpha = 0.7, aes(color = Outlier_Status)) +
  theme_bw() +
  #theme(legend.position = "none") +
  labs(title = "PCA on Peptides") +
  scale_color_manual(values = c('#22c6c9', 'orange'), name = "Outlier Status")
pca_plot_outlier


## Sample selection using PCA ----------------------------------------------------------------------------------------
# filter out outliers from sample descriptions in the PCA scores df
pca_scores_normals <- pca_scores %>% 
  filter(Outlier_Status == "normal")

# order by PC1 value for each sample and designate the ranking for each
pca_scores_normals$Order_PC1 <- rank(-pca_scores_normals$PC1, ties.method = "first")

# generate a sequence of numbers that are evenly spaced between the number of samples in the set
filtering_seq <- round(seq(1, length(pca_scores_normals$Order_PC1), length = 50))

# filter the dataframe of PCA scores to only include the 50 samples from the sequence of ranks of PCA scores
sample_list_PC1 <- pca_scores_normals %>% 
  filter(pca_scores_normals$Order_PC1 %in% filtering_seq)

# plot PC1 to see arrangement of samples
PC1_sample_list_plot <- sample_list_PC1 %>% 
  ggplot(aes(x = PC1, y = 0)) +
  geom_point(size = 2, alpha = 0.7, color = '#22c6c9') +
  theme_bw() +
  labs(title = "Distribution of PC1 for Selected Samples",
       y = "")
PC1_sample_list_plot

# PUll out the corresponding sample numbers that map back to the chosen PCA scores 
sample_list_seq <- sample_list_PC1$SampleNum

# Generate peptide list from selected samples
peptide_panel_filtered_normals <- peptide_panel_filtered %>% 
  filter(SampleNum %in% sample_list_seq)

# generate the list of peptides that are identified by samples selected from PCA
peptide_list_filtered_normals <- peptide_panel_filtered_normals$peptide_ID %>% 
  unique()

# Adding annotation, specifically outlier status, to the peptide panel dataframe
peptide_panel_normal_annotated <- merge(peptide_panel_filtered, sample_annotation_selected, by = "SampleNum")

# This dataframe has the peptides for each 'normal' sample
peptide_panel_normal_annotated <- peptide_panel_normal_annotated %>% 
  filter(Outlier_Status == "normal")

# generate the list of peptides that are identified by all normal samples
peptide_list_normals <- peptide_panel_normal_annotated$peptide_ID %>% 
  unique()

list_compare_50 <- list('Selected Samples (50)' = peptide_list_filtered_normals, 'Full Set' = peptide_list_normals)
ggvenn(list_compare_50)


## Random sample selection ------------------------------------------------------------------------------------------
set.seed(123)

# create a vector of the sample numbers designated as 'normal'
list_sample_num_normals <- peptide_panel_normal_annotated$SampleNum %>% 
  unique()

# randomly sample 50 times from the list of possible 'normal' sample numbers
rand01 <- sample(list_sample_num_normals, 50)
rand02 <- sample(list_sample_num_normals, 50)
rand03 <- sample(list_sample_num_normals, 50)
rand04 <- sample(list_sample_num_normals, 50)
rand05 <- sample(list_sample_num_normals, 50)
rand06 <- sample(list_sample_num_normals, 50)
rand07 <- sample(list_sample_num_normals, 50)
rand08 <- sample(list_sample_num_normals, 50)
rand09 <- sample(list_sample_num_normals, 50)
rand10 <- sample(list_sample_num_normals, 50)

peptide_panel_rand01 <- peptide_panel_filtered %>% 
  filter(SampleNum %in% rand01)

peptide_panel_rand02 <- peptide_panel_filtered %>% 
  filter(SampleNum %in% rand02)

peptide_panel_rand03 <- peptide_panel_filtered %>% 
  filter(SampleNum %in% rand03)

peptide_panel_rand04 <- peptide_panel_filtered %>% 
  filter(SampleNum %in% rand04)

peptide_panel_rand05 <- peptide_panel_filtered %>% 
  filter(SampleNum %in% rand05)

peptide_panel_rand06 <- peptide_panel_filtered %>% 
  filter(SampleNum %in% rand06)

peptide_panel_rand07 <- peptide_panel_filtered %>% 
  filter(SampleNum %in% rand07)

peptide_panel_rand08 <- peptide_panel_filtered %>% 
  filter(SampleNum %in% rand08)

peptide_panel_rand09 <- peptide_panel_filtered %>% 
  filter(SampleNum %in% rand09)

peptide_panel_rand10 <- peptide_panel_filtered %>% 
  filter(SampleNum %in% rand10)

peptide_panel_rand01$Set_ID <- 1
peptide_panel_rand02$Set_ID <- 2
peptide_panel_rand03$Set_ID <- 3
peptide_panel_rand04$Set_ID <- 4
peptide_panel_rand05$Set_ID <- 5
peptide_panel_rand06$Set_ID <- 6
peptide_panel_rand07$Set_ID <- 7
peptide_panel_rand08$Set_ID <- 8
peptide_panel_rand09$Set_ID <- 9
peptide_panel_rand10$Set_ID <- 10
peptide_panel_filtered_normals$Set_ID <- "PC1"

peptide_lists_combined <- rbind(peptide_panel_rand01, peptide_panel_rand02)
peptide_lists_combined <- rbind(peptide_lists_combined, peptide_panel_rand03)
peptide_lists_combined <- rbind(peptide_lists_combined, peptide_panel_rand04)
peptide_lists_combined <- rbind(peptide_lists_combined, peptide_panel_rand05)
peptide_lists_combined <- rbind(peptide_lists_combined, peptide_panel_rand06)
peptide_lists_combined <- rbind(peptide_lists_combined, peptide_panel_rand07)
peptide_lists_combined <- rbind(peptide_lists_combined, peptide_panel_rand08)
peptide_lists_combined <- rbind(peptide_lists_combined, peptide_panel_rand09)
peptide_lists_combined <- rbind(peptide_lists_combined, peptide_panel_rand10)
peptide_lists_combined <- rbind(peptide_lists_combined, peptide_panel_filtered_normals)

peptide_lists_combined_summary <- peptide_lists_combined %>% 
  group_by(Set_ID) %>% 
  summarise(Num_Peptides = n_distinct(peptide_ID, na.rm = TRUE))

peptide_lists_combined_summary <- peptide_lists_combined_summary %>% 
  mutate(Percent_Coverage = Num_Peptides/length(peptide_list_normals))

peptide_lists_combined_summary$Set_ID <- factor(peptide_lists_combined_summary$Set_ID, 
                                                levels = c("PC1", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

peptide_list_plot <- peptide_lists_combined_summary %>% 
  ggplot(aes(x = Set_ID, y = Num_Peptides, text = round(Percent_Coverage, digits = 3))) +
  geom_bar(stat = "identity", fill = '#22c6c9') +
  geom_text(aes(label = round(Percent_Coverage, digits = 3)), vjust = 1.5, hjust = 0.5, color = "white") +
  geom_text(aes(label = Num_Peptides), vjust = -0.5, hjust = 0.5) +
  theme_bw() +
  labs(title = "Number of Peptides for Samples Randomly Selected vs PCA",
       x = "Set List",
       y = "Peptides")
peptide_list_plot


## Sample Selection with Rolling Highest IDs ------------------------------------------------------------------------
# set size of sample set
set_size <- 57

# create empty list with length equal to the number of samples that will comprise the library set
sample_list_max <- vector(mode = "list", length = 0)

# list of samples that cannot be used for spectral library generation because of insufficient volume or remaining peptide
sample_no_go <- c(6, 61, 74, 114, 362, 429, 141, 151)

# create df of peptides and associated samples that will be sequentially filtered down as samples are removed
peptide_panel_loop_filtered <- peptide_panel_normal_annotated %>% 
  filter(!SampleNum %in% sample_no_go)

# loop that sequentially filters for successive samples starting at 1 and determines the number of unique peptides for each group
for(i in 1:set_size) { 
  # calculate the number of peptide ids for each sample that is remaining in the dataframe of peptides
  peptide_ids_loop_filtered <- peptide_panel_loop_filtered %>% 
    group_by(SampleNum) %>% 
    summarize(Num_Peptides = n_distinct(peptide_ID, na.rm = TRUE))
  
  # append the sample number with the maximum number of peptide ids in the remaining list
  sample_list_max <- append(sample_list_max, as.integer(peptide_ids_loop_filtered[which.max(peptide_ids_loop_filtered$Num_Peptides), 1]))  
  
  list_peptides_to_filter <- peptide_panel_loop_filtered %>% 
    filter(SampleNum == as.integer(peptide_ids_loop_filtered[which.max(peptide_ids_loop_filtered$Num_Peptides), 1])) %>% 
    pull(peptide_ID)
  
  # check that the list of peptides pulled from the list matches the number of peptide IDs calculated in the first step of the loop
  if (length(list_peptides_to_filter) != as.integer(peptide_ids_loop_filtered[which.max(peptide_ids_loop_filtered$Num_Peptides), 2])) {
    print(paste("Number of peptide IDs in filtered list did not match calculated IDs at i =", i, sep = " "))
    print(paste("The number of peptides in the list of peptides to filter is", length(list_peptides_to_filter), sep = " "))
    print(paste("The number of peptides determined from the panel is", as.integer(peptide_ids_loop_filtered[which.max(peptide_ids_loop_filtered$Num_Peptides), 1]), sep = " "))
  }
  
  # filter peptide_panel_loop_filtered for all peptides that are covered by the sample number with the highest number of IDs
  peptide_panel_loop_filtered <- peptide_panel_loop_filtered %>% 
    filter(!peptide_ID %in% list_peptides_to_filter) # filters out the peptides covered by the sample with the highest number of IDs
  
}

sample_list_max57_df <- as.data.frame(sample_list_max)
sample_list_max57_df <- t(sample_list_max57_df)

sample_list_max75_df <- as.data.frame(sample_list_max)
sample_list_max75_df <- t(sample_list_max75_df)

outlier_list <- pca_scores %>% 
  filter(Outlier_Status == "outlier") %>% 
  select(SampleNum)

# Use the list of sample numbers from the loop to filter the panel of peptides and samples
peptide_panel_from_loop <- peptide_panel_filtered %>% 
  filter(SampleNum %in% sample_list_max)

# Annotate these entries as from the loop that uses the highest IDs
peptide_panel_from_loop$Set_ID <- "Max_IDs"

# Combined with the dataframe with the prior lists of peptides from each sample selection methodology
peptide_lists_combined_max <- rbind(peptide_lists_combined, peptide_panel_from_loop)

# Calculate the number of peptide IDs from each method (Set ID)
peptide_lists_combined_max_summary <- peptide_lists_combined_max %>% 
  group_by(Set_ID) %>% 
  summarise(Num_Peptides = n_distinct(peptide_ID, na.rm = TRUE))

# Calculate the percent coverage for each method set
peptide_lists_combined_max_summary <- peptide_lists_combined_max_summary %>% 
  mutate(Percent_Coverage = Num_Peptides/length(peptide_list_normals))

# Convert Set_ID to factor and reorder factor level
peptide_lists_combined_max_summary$Set_ID <- factor(peptide_lists_combined_max_summary$Set_ID, 
                                                    levels = c("Max_IDs", "PC1", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

peptide_list_max_figure <- peptide_lists_combined_max_summary %>% 
  ggplot(aes(x = Set_ID, y = Num_Peptides, text = round(Percent_Coverage, digits = 3))) +
  geom_bar(stat = "identity", fill = '#22c6c9') +
  geom_text(aes(label = round(Percent_Coverage, digits = 3)), vjust = 1.5, hjust = 0.5, color = "white", fontface = "bold", size = 6) +
  geom_text(aes(label = Num_Peptides), vjust = -0.5, hjust = 0.5, fontface = "bold", size = 6, color = '#22c6c9') +
  theme_bw() +
  labs(title = "Comparison of Peptide Coverage Between Selection Strategies",
       x = "Set List",
       y = "Peptides") +
  theme(plot.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 22, face = "bold")
  ) +
  ylim(0, 20000)
peptide_list_max_figure


## Modeling Sample Selection with Rolling Highest IDs ------------------------------------------------------------------------
# set size of iterations to model
model_size <- 100

# create empty list with length equal to the number of samples that will comprise the library set
sample_list_max <- vector(mode = "list", length = 0)

# create df of peptides and associated samples that will be sequentially filtered down as samples are removed
peptide_panel_loop_filtered <- peptide_panel_normal_annotated

# initialize new df that stores ID numbers for each set size
model_ids <- data.frame(matrix(0, ncol = model_size))

# loop that runs through iterations of the model/library size
for (j in 1:model_size) {
  # reset the temporary peptide list to its original state that will be filtered in the inner loop
  peptide_panel_loop_filtered_temp <- peptide_panel_loop_filtered
  
  # loop that sequentially filters for successive samples starting at 1 and determines the number of unique peptides for each group
  for(i in 1:j) {
    # calculate the number of peptide ids for each sample that is remaining in the dataframe of peptides
    peptide_ids_loop_filtered <- peptide_panel_loop_filtered_temp %>% 
      group_by(SampleNum) %>% 
      summarize(Num_Peptides = n_distinct(peptide_ID, na.rm = TRUE))
    
    # append the sample number with the maximum number of peptide ids in the remaining list
    sample_list_max <- append(sample_list_max, as.integer(peptide_ids_loop_filtered[which.max(peptide_ids_loop_filtered$Num_Peptides), 1]))  
    
    # create a list of peptides that are id'ed by the sample number with the highest number of IDs
    list_peptides_to_filter <- peptide_panel_loop_filtered_temp %>% 
      filter(SampleNum == as.integer(peptide_ids_loop_filtered[which.max(peptide_ids_loop_filtered$Num_Peptides), 1])) %>% 
      pull(peptide_ID)
    
    # check that the list of peptides pulled from the list matches the number of peptide IDs calculated in the first step of the loop
    if (length(list_peptides_to_filter) != as.integer(peptide_ids_loop_filtered[which.max(peptide_ids_loop_filtered$Num_Peptides), 2])) {
      print(paste("Number of peptide IDs in filtered list did not match calculated IDs at i =", i, sep = " "))
      print(paste("The number of peptides in the list of peptides to filter is", length(list_peptides_to_filter), sep = " "))
      print(paste("The number of peptides determined from the panel is", as.integer(peptide_ids_loop_filtered[which.max(peptide_ids_loop_filtered$Num_Peptides), 1]), sep = " "))
    }
    
    # filter peptide_panel_loop_filtered for all peptides that are covered by the sample number with the highest number of IDs
    peptide_panel_loop_filtered_temp <- peptide_panel_loop_filtered_temp %>% 
      filter(!peptide_ID %in% list_peptides_to_filter) # filters out the peptides covered by the sample with the highest number of IDs
    
  }
  
  peptide_panel_from_each_iteration <- peptide_panel_filtered %>% 
    filter(SampleNum %in% sample_list_max) # filter the samples chosen from the selection loop from the initial peptide dataset
  model_ids[, j] <- length(unique(peptide_panel_from_each_iteration$peptide_ID)) # adds the number of peptide ids from each iteration of the outer loop
  sample_list_max <- vector(mode = "list", length = 0) # reset the sample_list_max list
}

model_ids_summary <- model_ids %>% 
  pivot_longer(cols = starts_with("X") ,names_to = "Code", values_to = "Peptides")

model_ids_summary$Num_Samples <- as.numeric(str_sub(model_ids_summary$Code, 2))

model_ids_summary <- model_ids_summary %>% 
  mutate(Percent_Coverage = Peptides/length(peptide_list_normals))

# figure
plot_model_figure <- model_ids_summary %>% 
  ggplot(aes(x = Num_Samples, y = Peptides, text = Percent_Coverage)) +
  geom_line(linetype = "dashed") +
  geom_point(size = 4, alpha = 0.7, color = '#22c6c9') +
  theme_bw() +
  labs(title = "Peptides Per Sample Set Size",
       x = "Set Size") +
  theme(plot.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 22, face = "bold"))

plot_model_figure
