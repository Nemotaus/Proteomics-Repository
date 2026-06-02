#Remove all data from environment
rm(list = ls())

#unload all packages
lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))

#attach dependencies
library(tidyverse)
library(reshape2)
library(svglite)
library(car)
library(multcompView)
library(broom)
library(VIM)

#load a single csv file
daf1 <- read.csv(file.choose())

#This line of code is for loading in multiple CSV files from a single
#master directory. Use this code for combining multiple MS Runs.

# Function to select and read multiple files
read_multiple_files <- function() {
  file_paths <- character(0)  # Initialize an empty character vector
  while(TRUE) {
    file_path <- tryCatch(file.choose(), error = function(e) return(NULL))
    if (is.null(file_path) || identical(file_path, "")) break  # Check if selection was cancelled
    file_paths <- c(file_paths, file_path)
  }
  return(file_paths)
}

#These arguments will create a separate data frame for each CSV
#read multiple files will collect and store the file paths
file_list <- read_multiple_files()

#This code will make each CSV file its own data frame
df_list <- lapply(file_list, read.csv)

# Merge all data frames on specific columns, handle single file case
if (length(df_list) == 1) {
  daf1 <- df_list[[1]]
} else {
  daf1 <- df_list %>%
    reduce(function(x, y) merge(x,
                                y,
                                by = c("Protein.Group", "Genes", "Protein.Names"),
                                all = TRUE))
}

#Median Normalization
#Remove HeLa Column
daf1 <- daf1 %>%
  select(-HeLa)

#Remove Protein columns. This does not remove the protein columns? But it does rearrange them
daf1 <- daf1 %>%
  select(-Protein...x, -Protein...y)

#Remove rows with spontaneous hits from the HeLa column
daf1 <- daf1 %>%
  rowwise() %>%
  filter(!if_all(4:15, is.na)) %>%
  ungroup()

# Filter the data frame as a coincidence finder where original mice have NAs while other mice have protein abundnace
filtered_dafOri <- daf1 %>%
  rowwise() %>%
  filter(
    # Check if columns 4, 5, 10, and 11 have any NA
    if_all(c(4, 5, 10, 11), is.na) &
      # Check if all other columns (6-9 and 12-15) are not NA
      if_all(c(6:9, 12:15), ~ !is.na(.))
  ) %>%
  ungroup()

# Filter the data frame as a coincidence finder where original mice have NAs while other mice have protein abundnace
filtered_dafCor <- daf1 %>%
  rowwise() %>%
  filter(
    # Check if columns 4, 5, 10, and 11 have any NA
    if_all(c(6, 7, 12, 13), is.na) &
      # Check if all other columns (6-9 and 12-15) are not NA
      if_all(c(4:5, 8:11, 14:15), ~ !is.na(.))
  ) %>%
  ungroup()

# Filter the data frame as a coincidence finder where original mice have NAs while other mice have protein abundnace
filtered_dafWT <- daf1 %>%
  rowwise() %>%
  filter(
    # Check if columns 4, 5, 10, and 11 have any NA
    if_all(c(8, 9, 14, 15), is.na) &
      # Check if all other columns (6-9 and 12-15) are not NA
      if_all(c(4:7, 10:13), ~ !is.na(.))
  ) %>%
  ungroup()

# Export the dataframe 'daf1' as a CSV file
write.csv(daf1, "daf1_exported.csv", row.names = FALSE)

# Calculate the median of columns 4-15, ignoring NAs
medians <- daf1 %>%
  select(4:15) %>%
  map_dbl(~ median(., na.rm = TRUE))

# Print the medians
print(medians)

# Reshape the data into long format
long_daf1 <- daf1 %>%
  select(4:15) %>%
  pivot_longer(cols = everything(), names_to = "variable", 
               values_to = "Log2 Abundance (AU)")

long_daf1 <- long_daf1 %>%
  mutate(variable = factor(variable,
                           levels = c('MA1', 'MA1.1','MA3','MA4',
                                      'X3280...', 'X3289...','X3224...','X3291...', 
                                      'X3238...', 'X3277...', 'X2093...','X3225...')))
                                      
                                    
# Create the box plot of un normalized data
premedian <- ggplot(long_daf1, aes(x = variable, y = `Log2 Abundance (AU)`)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 35)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),   # removes the "Genotype" label
    axis.ticks.x = element_blank(),    # removes the tick marks themselves
    axis.text.y = element_text(size = 10),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(premedian)

# Save as SVG
ggsave("Medians Prior to Norm.svg", plot = premedian, width = 13.333, height = 7.5, units = "in", dpi = 600)

# Select columns 4-15 for calculation
cols_to_normalize <- daf1[, 4:15]
cols_to_normalize <- apply(cols_to_normalize, 2, as.numeric)

# Step 1: Calculate global median of columns 4-15
global_median <- median(as.numeric(unlist(cols_to_normalize)), na.rm = TRUE)

# Step 2: Calculate median for each column individually
column_medians <- apply(cols_to_normalize, 2, median, na.rm = TRUE)

# Step 3: Initialize a vector to store column median differences
col_median_diffs <- numeric(ncol(cols_to_normalize))

# Step 4: Normalize the data
normalized_data <- cols_to_normalize
for (col in 1:ncol(cols_to_normalize)) {
  col_median_diffs[col] <- global_median - column_medians[col]
  normalized_data[, col] <- normalized_data[, col] + col_median_diffs[col]
}

# Combine normalized data with the original non-normalized columns (assuming you want to keep other columns)
normalized_daf1 <- cbind(daf1[, -c(4:15)], normalized_data)

# Calculate the median of columns 2-7, ignoring NAs
normmedians <- normalized_daf1 %>%
  select(4:15) %>%
  map_dbl(~ median(., na.rm = TRUE))

# Print the medians
print(normmedians)

# Export the dataframe 'daf1' as a CSV file
write.csv(normalized_daf1, "normalized_daf1_exported.csv", row.names = FALSE)

# Get the column names of the normalized data
normalized_cols <- colnames(normalized_data)

# Reshape the data into long format
long_normalized_daf1 <- normalized_daf1 %>%
  pivot_longer(cols = all_of(normalized_cols), names_to = "variable", values_to = "value")

long_normalized_daf1 <- long_normalized_daf1 %>%
  mutate(variable = factor(variable,
                           levels = c('MA1', 'MA1.1','MA3','MA4',
                                      'X3280...', 'X3289...','X3224...','X3291...', 
                                      'X3238...', 'X3277...', 'X2093...','X3225...')))

# Create the box plot of un normalized data
postmedian <- ggplot(long_normalized_daf1, aes(x = variable, y = value)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 35)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),   # removes the "Genotype" label
    axis.ticks.x = element_blank(),    # removes the tick marks themselves
    axis.text.y = element_text(size = 10),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(postmedian)

# Save as SVG
ggsave("Medians After Norm.svg", plot = postmedian, width = 13.333, height = 7.5, units = "in", dpi = 600)

#Keep rows with murine p4 Atpases which interact with CDC50a
df2 <- normalized_daf1 %>% 
  filter(Genes %in% c('Atp8a1', 'Atp8a2', 'Atp8b1', 'Atp8b2', 'Atp8b3',
                      'Atp8b4', 'Atp10a', 'Atp10b', 'Atp10d', 'Atp11a',
                      'Atp11b', 'Atp11c', 'Tmem30a'))

#Alter the naming of the mice in the columns
#Edit the names for each set of mice that pass through this code
#Define new names
newN <-c('Protein','3280', '3289', '3238',
        '3277', 'MA1', 'MA2','3224','3291',
        '2093','3225','MA3','MA4')

#Define old names
oldN <- c('Genes', 'X3280...', 'X3289...', 'X3238...',
         'X3277...', 'MA1', 'MA1.1','X3224...','X3291...',
         'X2093...','X3225...','MA3','MA4')

#Rename the columns of df2
df2 <- df2 %>%
  rename_with(~newN, all_of(oldN))

#remove the Protein group and Protein names columns
df2 <- df2 %>%
  select(-Protein.Group, -Protein.Names)

#Convert the dataframe to long type
long_df2 <- melt(df2, id='Protein')

#add genotypes
#check this for each mass spec run
#I load samples as -/-, +/+, WT/WT to prevent cross contamination
long_df3 <- long_df2 %>%
  mutate(Genotype = case_when(
    variable %in% c('3280','3289','3224','3291') ~ "Original",
    variable %in% c('3238', '3277','2093','3225') ~ "Corrected",
    variable %in% c('MA1', 'MA2','MA3','MA4') ~ "Wildtype",
    TRUE ~ NA_character_
  ))

#rename columns
long_df3 <- long_df3 %>%
  rename("Mouse_ID" = variable,
         "Log2_Abundance" = value)

#Change the NA values to 0's
#This will allow us to divide the log2 abundance from originals by wildtypes and get 0 or none
long_df3[is.na(long_df3)] <- 0

#Reframe the data to wide format lol
wide_df4 <- long_df3 %>%
  spread(key = "Protein", value = "Log2_Abundance")

#The code breaks the numeric variables so they have to be reset
#edit this to incorporate all flippases encountered
wide_df4 <- wide_df4 %>%
  mutate(across(starts_with("Atp") | starts_with("Tmem"), as.numeric))
wide_df4$Mouse_ID <- as.character(wide_df4$Mouse_ID)

is_numeric <- sapply(wide_df4, is.numeric)
print(is_numeric)

# Calculate overall median
median_Tmem <- median(wide_df4$Tmem30a, na.rm = TRUE)
median_Atp10d <- median(wide_df4$Atp10d, na.rm = TRUE)

#Homoscedasticity of TMEM30a values to support its use in normalization
bartlett.test(data = wide_df4, Tmem30a ~ Genotype)

fligner.test(Tmem30a ~ Genotype, data = wide_df4)

#variance and standard deviation
variance_Tmem30a <- var(wide_df4$Tmem30a)
variance_Tmem30a

sd_Tmem30a <- sd(wide_df4$Tmem30a)
sd_Tmem30a


# Reorder Genotype factor
long_df3 <- long_df3 %>%
  mutate(Genotype = factor(Genotype, levels = c("Wildtype", 
                                                "Original", 
                                                "Corrected")))

# Protein list
proteins_of_interest <- c('Atp8a1', 'Atp8a2', 'Atp8b1', 'Atp8b2', 'Atp8b3',
                          'Atp8b4', 'Atp10a', 'Atp10b', 'Atp10d', 'Atp11a',
                          'Atp11b', 'Atp11c', 'Tmem30a')

# Custom genotype colors
genotype_colors <- c("Original" = "#1b1b1b",
                     "Corrected" = "#4c90ff",
                     "Wildtype" = "#ff9f1c")

# Compute mean and SEM per genotype per protein
raw_summary_stats <- long_df3 %>%
  group_by(Protein, Genotype) %>%
  summarise(
    mean_val = mean(Log2_Abundance, na.rm = TRUE),
    sem_val  = sd(Log2_Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Loop over proteins
for(prot in proteins_of_interest){
  
  # Filter data for this protein
  df_prot <- long_df3 %>% filter(Protein == prot)
  
  # Calculate median
  median_val <- median(df_prot$Log2_Abundance, na.rm = TRUE)
  
  # Calculate max with 10% extra headroom
  max_val <- max(df_prot$Log2_Abundance, na.rm = TRUE) * 1.1
  
  # Get mean ± SEM for this protein
  stat_prot <- raw_summary_stats %>% filter(Protein == prot)

  # Create ggplot
  p <- ggplot(df_prot, aes(x = Genotype, y = Log2_Abundance, color = Genotype)) +
    geom_jitter(width = 0.15, size = 4) +
    geom_errorbar(
      data = stat_prot,
      aes(x = Genotype, ymin = mean_val - sem_val, ymax = mean_val + sem_val),
      width = 0.2, color = "#4d4d4d", size = 0.7,
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(x = 3.0, y = median_val - 2, 
          label = paste0("Median: ", round(median_val,2))),
      color = "darkred",
      hjust = 0,
      inherit.aes = FALSE,
      size = 3.5
    ) +
    geom_hline(yintercept = median_val, color = "darkred", 
               linetype = "dashed", 
               linewidth = 0.7) +
    guides(color = guide_legend(title = "Genotype"))+  # ensures only one legend
    scale_color_manual(values = genotype_colors, name = "Genotype") +  # simple legend title
    scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, max_val)) +
    labs(title = prot, x = "Genotype", y = "Log2 Abundance") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),   # removes the "Genotype" label
      axis.text.x = element_blank(),    # removes the tick labels (e.g., Wildtype, Original, Corrected)
      axis.ticks.x = element_blank(),    # removes the tick marks themselves
      axis.text.y = element_text(size = 10),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "top",                # move legend to the top
      legend.title = element_text(size = 12), # legend title size
      legend.text = element_text(size = 10),  # legend text size
      legend.key.size = unit(0.5, "cm"),      # smaller legend key
      legend.background = element_blank(),    # remove grey box
      legend.key = element_blank()            # remove grey behind symbols
    )
  
  print(p)
  
  # Optional: save each as SVG
  ggsave(filename = paste0(prot, "_dotplot.svg"), plot = p,
         width = 5, height = 5, units = "in", dpi = 600)
}

#check this code as it may depend on the flippases captured from the mass spec
normalized_df5 <- wide_df4 %>%
  mutate(across(Atp10d:Tmem30a, ~ . / Tmem30a))

#statistics on the normalized flippase abundances

#Assuming df is already normalized and in long format, otherwise pivot_longer if needed
df5_long <- normalized_df5 %>%
  pivot_longer(cols = starts_with("Atp"), 
               names_to = "P4_type_ATPase", 
               values_to = "Normalized_Value")

# Separate data frames for Original/Wild Type and Corrected/Wild Type comparisons
df5_original <- df5_long %>%
  filter(Genotype %in% c("Original", "Wildtype")) %>%
  filter(!(P4_type_ATPase == "Atp10d"))

df5_corrected <- df5_long %>%
  filter(Genotype %in% c("Corrected", "Wildtype"))

#convert normalized data to long format
#Prepare data for long format
normalized_df6 <- normalized_df5 %>%
select(-Mouse_ID, -Tmem30a)

#convert to long format
Normalized_long <- melt(normalized_df6, 
                      id.vars = "Genotype", 
                      variable.name = "Flippase", 
                      value.name = "Fold_Change")

# Adjust the order of the Genotype factor levels
Normalized_long <- Normalized_long %>%
  mutate(Genotype = factor(Genotype, levels = c("Original", 
                                                "Corrected", 
                                                "Wildtype")))

#Adjust Flippase Order
Normalized_long <- Normalized_long %>%
  mutate(Flippase = factor(Flippase, levels = c("Atp8a1", 
                                                "Atp8b1", 
                                                "Atp10d",
                                                "Atp11a",
                                                "Atp11b",
                                                "Atp11c")
                           )
         )

#ggplot of individual data points
normalized_points <- ggplot(Normalized_long, aes(x = Flippase, 
                                                 y = Fold_Change, 
                                                 fill = Genotype)) +
  geom_boxplot(width = 0.5,
               outlier.shape = NA) +
  geom_dotplot(aes(group = interaction(Flippase, Genotype)), 
               binaxis = "y", 
               stackdir = "center", 
               dotsize = 0.5,
               position = position_dodge(width = 0.5)) + 
  labs(x = "Flippase", y = "Normalized Average", title = "Fold Change of P4-Type ATPases") +
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c("Original" = "grey", 
                               "Corrected" = "blue", 
                               "Wildtype" = "orange"))+
  theme_bw() + # Select theme with a white background  
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = 'black')
  )
print(normalized_points)

#normalized data
# Extract the mean wildtype values for normalization
wildtype_average <- Normalized_long %>%
  filter(Genotype == "Wildtype") %>%
  group_by(Flippase) %>%
  summarize(Wildtype_Average = mean(Fold_Change, na.rm = TRUE)) %>%
  ungroup()

# Normalize the data to the wildtype values for each flippase, including wildtype
df_normalized <- Normalized_long %>%
  left_join(wildtype_average, by = "Flippase") %>%
  mutate(Normalized_Value = Fold_Change / Wildtype_Average)

# Perform ANOVA followed by Tukey's post-hoc test
anova_results <- df_normalized %>%
  group_by(Flippase) %>%
  do({
    # Check if all values in the "Original" group are zero for Atp10D
    if (all(.$Normalized_Value[.$Genotype == "Original"] == 0) & unique(.$Flippase) == "Atp10D") {
      data.frame(comparison = NA, p.adj = NA, Flippase = unique(.$Flippase))
    } else {
      # Perform ANOVA
      model <- aov(Normalized_Value ~ Genotype, data = .)
      
      # Perform Tukey's HSD test
      tukey_results <- TukeyHSD(model)
      
      # Extract and format p-values from Tukey's results
      tukey_p_values <- as.data.frame(tukey_results$Genotype) %>%
        rownames_to_column(var = "comparison") %>%
        mutate(Flippase = unique(.$Flippase))
      
      # Return results
      tukey_p_values
    }
  }) %>%
  ungroup()

# View the ANOVA and Tukey's HSD results
print(anova_results)

# Adjust the order of the Genotype factor levels
df_normalized <- df_normalized %>%
  mutate(Genotype = factor(Genotype, levels = c("Original", 
                                                "Corrected", 
                                                "Wildtype")))

#Adjust Flippase Order
df_normalized <- df_normalized %>%
  mutate(Flippase = factor(Flippase, levels = c("Atp8a1", 
                                                "Atp8b1", 
                                                "Atp10d",
                                                "Atp11a",
                                                "Atp11b",
                                                "Atp11c")
  )
  )


# Create the plot
Normalized_to_WT <- ggplot(df_normalized, aes(x = Flippase, y = Normalized_Value, fill = Genotype)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_dotplot(aes(group = interaction(Flippase, Genotype)), 
               binaxis = "y", 
               stackdir = "center", 
               dotsize = 0.5,
               position = position_dodge(width = 0.5)) + 
  labs(title = "Normalized ATPase Expression",
       x = "Flippase",
       y = "Fold Change ATPase Expression") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("Original" = "grey", 
                               "Corrected" = "blue", 
                               "Wildtype" = "orange")) +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = 'black')
  )

#print
print(Normalized_to_WT)

#
# Export the dataframe 'wt_normalized' as a CSV file
write.csv(df_normalized, "df_normalized_WT.csv", row.names = FALSE)

# Save as SVG
ggsave("Normalized_to_Wildtype.svg", plot = Normalized_to_WT, width = 13.333, height = 7.5, units = "in", dpi = 600)

#volcano plot based off normalized daf1
# Identify and modify the row with Atp10d while dropping all NA value rows
normalized_daf2 <- normalized_daf1 %>%
  rowwise() %>%
  mutate(
    across(c(4, 5, 10, 11), ~ if_else(Genes == "Atp10d" & is.na(.), 0, .))
  ) %>%
  ungroup() %>%
  # Remove all rows containing at least 1 NA value between columns 4-15
  filter(if_all(c(4:15), ~ !is.na(.))) %>%
  # Drop columns 8, 9, 14, and 15
  select(-c(1,3, 8, 9, 14, 15))
  
#Alter the naming of the mice in the columns
#Edit the names for each set of mice that pass through this code
#Define new names
newN2 <-c('Protein','3280', '3289', '3238',
         '3277', '3224','3291',
         '2093','3225')

#Define old names
oldN2 <- c('Genes', 'X3280...', 'X3289...', 'X3238...',
          'X3277...', 'X3224...','X3291...',
          'X2093...','X3225...')

#Rename the columns of df2
normalized_daf2 <- normalized_daf2 %>%
  rename_with(~newN2, all_of(oldN2))

#Convert the dataframe to long type
long_normalized_daf2 <- melt(normalized_daf2, id='Protein')

#add genotypes
#check this for each mass spec run
long_normalized_daf3 <- long_normalized_daf2 %>%
  mutate(Genotype = case_when(
    variable %in% c('3280','3289','3224','3291') ~ "Original",
    variable %in% c('3238', '3277','2093','3225') ~ "Corrected",
    TRUE ~ NA_character_
  )
  )

# Remove rows where Protein is blank or NA
long_normalized_daf3 <- long_normalized_daf3 %>%
  filter(!is.na(Protein) & Protein != "")

#rename columns
long_normalized_daf3 <- long_normalized_daf3 %>%
  rename("Mouse_ID" = variable,
         "Log2_Abundance" = value)

#Summary statistics
results_volcano <- long_normalized_daf3 %>%
  group_by(Protein) %>%
  summarise(
    t_test = list(t.test(Log2_Abundance ~ Genotype)),
    .groups = 'drop'
  ) %>%
  mutate(
    # Extract p-value from the T-test result
    p_value = map_dbl(t_test, ~ tidy(.)$p.value[1])
  ) %>%
  select(Protein, p_value)

# Calculate fold changes
fold_changes <- long_normalized_daf3 %>%
  group_by(Protein, Genotype) %>%
  summarise(
    avg_abundance = mean(Log2_Abundance, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_wider(names_from = Genotype, values_from = avg_abundance) %>%
  # Calculate fold change: Corrected / Original
  mutate(
    fold_change = Corrected / Original
  )

# Combine T-test results with fold changes
volcano_daf <- results_volcano %>%
  left_join(fold_changes, by = "Protein")

#Calculate TwoFC and Log10p
volcano_daf2 <- volcano_daf %>%
  mutate(TwoFC = log2(fold_change)) %>%
  mutate(Log10p = -log10(p_value))

# Export the dataframe 'volcano_daf2' as a CSV file
write.csv(volcano_daf2, "volcano_daf2.csv", row.names = FALSE)

# Replace Inf with 8 for plotting
volcano_daf2 <- volcano_daf2 %>%
  mutate(
    TwoFC = case_when(
      is.infinite(TwoFC) & TwoFC > 0 ~ 8,   # Replace positive Inf with 8
      is.infinite(TwoFC) & TwoFC < 0 ~ -8,  # Replace negative Inf with -8
      TRUE ~ TwoFC  # Leave other values unchanged
    )
  )

#Coloring points based off threshold
volcano_daf2 <- volcano_daf2 %>%
  mutate(Change = case_when(TwoFC >= 1 & Log10p >= 1.30103 ~ "up",
                                TwoFC <= -1 & Log10p >= 1.30103 ~ "down",
                                TRUE ~ "ns")) 
#establishes classifications
volcano_daf2$Class <- "Other"  # Default class
volcano_daf2$Class[volcano_daf2$Protein == "Atp8a1"] <- "Class 1a"
volcano_daf2$Class[volcano_daf2$Protein == "Atp8b1"] <- "Class 1b"
volcano_daf2$Class[volcano_daf2$Protein == "Atp10d"] <- "Class 5"
volcano_daf2$Class[volcano_daf2$Protein %in% c("Atp11a", "Atp11b", "Atp11c")] <- "Class 6"
volcano_daf2$Class[volcano_daf2$Protein == "Tmem30a"] <- "Cdc"

# Define the color mapping 
color_mapping <- c("Class 1a" = "#004949",
                   "Class 1b" = "#db6d00",
                   "Class 5" = "#6db6ff",
                   "Class 6" = "#490092",
                   "Cdc" = "#920000",
                   "Other" = "grey")

# Define the size mapping
sizes_v <- c("Class 1a" = 3, 
             "Class 1b" = 3, 
             "Class 5" = 3, 
             "Class 6" = 3, 
             "Cdc" = 3,
             "Other" = 1)

# Define the alpha mapping
alphas_v <- c("Class 1a" = 1, 
              "Class 1b" = 1, 
              "Class 5" = 1, 
              "Class 6" = 1, 
              "Cdc" = 1,
              "Other" = 0.5)

# Separate flippases of interest and other proteins
flippases_of_interest <- volcano_daf2[volcano_daf2$Class %in% c("Class 1a", "Class 1b", "Class 5", "Class 6", "Cdc"), ]
other_proteins <- volcano_daf2[volcano_daf2$Class == "Other", ]

# Generate Plot
volcano <- ggplot() +
  # Plot other proteins first
  geom_point(data = other_proteins, 
             aes(TwoFC, Log10p, fill = Class, size = Class, alpha = Class), 
             shape = 21, colour = "black") +
  # Plot flippases of interest on top
  geom_point(data = flippases_of_interest, 
             aes(TwoFC, Log10p, fill = Class, size = Class, alpha = Class), 
             shape = 21, colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth=0.5) + 
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed", linewidth=0.5) +
  scale_fill_manual(values = color_mapping) + # Modify point colour with class colors
  scale_size_manual(values = sizes_v) + # Modify point size with class sizes
  scale_alpha_manual(values = alphas_v) + # Modify point transparency with class alphas
  scale_x_continuous(breaks = c(seq(-8, 8, 2)), # Modify x-axis tick intervals    
                     limits = c(-8, 8)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) +
  labs(x = "log2(Fold Change)",
       y = "-log10(P-Value)",
       fill = "Flippase Class",  # Update legend title
       colour = "Expression \nchange") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        legend.position="none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_line(linewidth = 0.5),
        axis.ticks.length = unit(.25, "cm"))

# Display the plot
print(volcano)

# Save as SVG
ggsave("Volcano.svg", plot = volcano, width = 13.333, height = 7.5, units = "in", dpi = 600)

##############################################
###############################################
################################################
#############################################
################################################
#Secret code

#Remove rows with spontaneous hits from the HeLa column
df3 <- df2 %>%
  rowwise() %>%
  filter(!if_any(2:13, is.na)) %>%
  ungroup()

#Export the dataframe 'daf1' as a CSV file
write.csv(df3, "df3 all common.csv", row.names = FALSE)
