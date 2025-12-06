###### **********************
## Assignment 4
##
## Gia Ly
##
## 2025-12-05
##
##***************************

## _ Packages used -------
## ONLY INSTALL IF NEEDED ##
# install.packages("tidyverse")
# install.packages("viridis")
# install.packages("ggplot2")
# install.packages("rentrez")
# insatll.packages("Biosintings")
# install.packages("stringr")
# install.packages("randomForest")
# install.packages("extrafont")
# install.packages("pROC")

library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("as.factor", "base")
library(viridis)
theme_set(theme_light())
library(ggplot2)
library(rentrez)
library(Biostrings)
library(stringr)
library(randomForest)
library(e1071)
library(pROC)
#library(extrafont)
#extrafont::font_import(prompt = FALSE) #sorry... I really dislike Arial :,)
#extrafont::loadfonts(device = "win")

## _ Creating Functions To Import Data From NCBI -------

# creating function to fetch FASTA files from NCBI nuccore based on a provided search term
FetchFastaFiles <- function(searchTerm, seqsPerFile = 100, fastaFileName) {
  search1 <- entrez_search(db = "nuccore", term = searchTerm)
  search2 <- entrez_search(db = "nuccore", term = searchTerm, retmax = search1$count, use_history = T)
  for (start_rec in seq(0, search2$retmax, seqsPerFile)) {
    fname <- paste(fastaFileName, start_rec, ".fasta", sep = "")
    recs <- entrez_fetch(
      db = "nuccore", web_history = search2$web_history,
      Sys.sleep(0.34),
      rettype = "fasta", retstart = start_rec, retmax = seqsPerFile
    )
    write(recs, fname)
    print(paste("Wrote records to ", fname, sep = ""))
  }

  return(search2)
}

# creating function to merge multiple FASTA files into one dataframe
MergeFastaFiles <- function(filePattern) {
  fastaFiles <- list.files(pattern = filePattern)
  l_fastaFiles <- lapply(fastaFiles, readDNAStringSet)
  l_dfFastaFiles <- lapply(l_fastaFiles, function(x) data.frame(Title = names(x), Sequence = paste(x)))
  dfSeqs <- do.call("rbind", l_dfFastaFiles)
  return(dfSeqs)
}

## __ Importing + Merging FAMILY -------

# solanaceae
FetchFastaFiles(
  searchTerm = "Solanaceae[Organism] AND rbcL[Gene] AND 400:2000[Sequence Length]",
  fastaFileName = "solanaceae_rbcL_"
)

solanaceae_df <- MergeFastaFiles("solanaceae_rbcL_")

# convolvulaceae
FetchFastaFiles(
  searchTerm = "Convolvulaceae[Organism] AND rbcL[Gene] AND 400:2000[Sequence Length]",
  fastaFileName = "convolvulaceae_rbcL_"
)

convolvulaceae_df <- MergeFastaFiles("convolvulaceae_rbcL_")

## __ Importing GENUS -------

# solanum
FetchFastaFiles(
  searchTerm = "Solanum[Organism] AND rbcL[Gene] AND 400:2000[Sequence Length]",
  fastaFileName = "solanum_rbcL_"
)

solanum_df <- MergeFastaFiles("solanum_rbcL_")

# capsicum
FetchFastaFiles(
  searchTerm = "Capsicum[Organism] AND rbcL[Gene] AND 400:2000[Sequence Length]",
  fastaFileName = "capsicum_rbcL_"
)

capsicum_df <- MergeFastaFiles("capsicum_rbcL_")

## __ Importing SPECIES -------

# solanum lycopersicum
FetchFastaFiles(
  searchTerm = "\"Solanum lycopersicum\"[Organism] AND rbcL[Gene] AND 400:2000[Sequence Length]",
  fastaFileName = "solanum_lycopersicum_rbcL_"
)

lycopersicum_df <- MergeFastaFiles("solanum_lycopersicum_rbcL_")

# solanum tuberosum
FetchFastaFiles(
  searchTerm = "\"Solanum tuberosum\"[Organism] AND rbcL[Gene] AND 400:2000[Sequence Length]",
  fastaFileName = "solanum_tuberosum_rbcL_"
)

tuberosum_df <- MergeFastaFiles("solanum_tuberosum_rbcL_")

# solanum melongena
FetchFastaFiles(
  searchTerm = "\"Solanum melongena\"[Organism] AND rbcL[Gene] AND 400:2000[Sequence Length]",
  fastaFileName = "solanum_melongena_rbcL_"
)

melongena_df <- MergeFastaFiles("solanum_melongena_rbcL_")

## _ Cleaning and Sorting Data -------

# putting all of the df in one list to streamline cleaning process
raw_dfs <- list(
  solanaceae = solanaceae_df,
  convolvulaceae = convolvulaceae_df,
  solanum = solanum_df,
  capsicum = capsicum_df,
  lycopersicum = lycopersicum_df,
  tuberosum = tuberosum_df,
  melongena = melongena_df
)

# below are concepts borrowed from Karl's Random Forest script (oct 21st)
clean_and_feature <- function(df) {
  # adding labels from title
  df <- df %>%
    mutate(
      unique_id    = word(Title, 1),
      genus_name   = word(Title, 2),
      species_name = paste(word(Title, 2), word(Title, 3)),
      nucleotides  = Sequence
    ) %>%
    filter(!is.na(nucleotides))

  # clean sequences
  df <- df %>%
    mutate(nucleotides2 = str_remove(nucleotides, "^[-N]+")) %>%
    mutate(nucleotides2 = str_remove(nucleotides2, "[-N]+$")) %>%
    mutate(nucleotides2 = str_remove_all(nucleotides2, "-+")) %>%
    filter(str_count(nucleotides2, "N") <= 0.05 * nchar(nucleotides2))

  # quartile-based length filter
  q1 <- quantile(nchar(df$nucleotides2), 0.25)
  q3 <- quantile(nchar(df$nucleotides2), 0.75)

  df <- df %>%
    filter(nchar(nucleotides2) >= q1 & nchar(nucleotides2) <= q3)

  # convert to DNAStringSet
  df$nucleotides2 <- DNAStringSet(df$nucleotides2)

  # add nucleotide frequencies
  df <- cbind(df, as.data.frame(letterFrequency(df$nucleotides2,
    letters = c("A", "C", "G", "T")
  )))

  # add proportions
  df <- df %>%
    mutate(
      Aprop = A / (A + C + G + T),
      Tprop = T / (A + C + G + T),
      Gprop = G / (A + C + G + T)
    )

  # add di- and trinucleotide frequencies
  df <- cbind(df, as.data.frame(dinucleotideFrequency(df$nucleotides2, as.prob = TRUE)))
  df <- cbind(df, as.data.frame(trinucleotideFrequency(df$nucleotides2, as.prob = TRUE)))

  # convert back to characters for tidyverse use
  df$nucleotides2 <- as.character(df$nucleotides2)

  return(df)
}

# applying the function to the whole list
cleaned_dfs <- lapply(raw_dfs, clean_and_feature)

# sorting taxons
### FAMILY-LEVEL CLASSIFICATION###
family_df <- bind_rows(
  cleaned_dfs$solanaceae %>% mutate(group = "Solanaceae"),
  cleaned_dfs$convolvulaceae %>% mutate(group = "Convolvulaceae")
)

### GENUS-LEVEL CLASSIFICATION###
genus_df <- bind_rows(
  cleaned_dfs$solanum %>% mutate(group = "Solanum"),
  cleaned_dfs$capsicum %>% mutate(group = "Capsicum")
)

### SPECIES-LEVEL CLASSIFICATION###
species_df <- bind_rows(
  cleaned_dfs$lycopersicum %>% mutate(group = "lycopersicum"),
  cleaned_dfs$tuberosum %>% mutate(group = "tuberosum"),
  cleaned_dfs$melongena %>% mutate(group = "melongena")
)
## _ Sequence Length (Figure 1) ----
family_df %>%
  mutate(seq_length = nchar(.data$Sequence)) %>%
  ggplot(aes(x = seq_length, fill = group)) +
  geom_histogram(alpha = 0.7, bins = 30, color = "white") +
  theme_minimal(base_size = 16) +
  theme(text = element_text(family = "Calibri", size = 16)) +
  labs(
    title = "Sequence Length Distribution — Family Level",
    x = "Sequence Length (bp)",
    y = "Count",
    fill = "Group"
  ) +
  scale_fill_manual(values = c(
    "Solanaceae" = "#C8A2C8",
    "Convolvulaceae" = "#8A9A5B"
  ))

## _ Creating Training Functions -------

split_train_validation <- function(df, label_col, prop = 0.2) {
  min_n <- min(table(df[[label_col]]))

  set.seed(123)
  df_validation <- df %>%
    group_by(.data[[label_col]]) %>%
    sample_n(floor(prop * min_n))

  df_training <- df %>%
    filter(!unique_id %in% df_validation$unique_id) %>%
    group_by(.data[[label_col]]) %>%
    sample_n(ceiling((1 - prop) * min_n))

  return(list(train = df_training, validation = df_validation))
}

# random forest
train_rf <- function(train_df, label_col, feature_cols) {
  rf <- randomForest(
    x = train_df[, feature_cols],
    y = base::as.factor(train_df[[label_col]]),
    ntree = 300,
    importance = TRUE
  )

  return(rf)
}

# evaluating random forest
evaluate_rf <- function(model, validation_df, label_col, feature_cols) {
  preds <- predict(model, validation_df[, feature_cols])
  cm <- table(
    observed = validation_df[[label_col]],
    predicted = preds
  )
  return(cm)
}

# SVM classifier (RBF kernel)
train_svm <- function(train_df, label_col, feature_cols) {
  svm(
    x = train_df[, feature_cols],
    y = as.factor(train_df[[label_col]]),
    kernel = "radial",
    probability = TRUE
  )
}

# evaluating SVM
evaluate_svm <- function(model, validation_df, label_col, feature_cols) {
  preds <- predict(model, validation_df[, feature_cols])
  cm <- table(
    observed = validation_df[[label_col]],
    predicted = preds
  )
  return(cm)
}


## _ Running Random Forest Classifier (1) -------

# prepare feature columns (numeric only)

get_feature_cols <- function(df) {
  remove <- c(
    "Title", "Sequence", "unique_id", "species_name",
    "genus_name", "group", "nucleotides", "nucleotides2"
  )
  
  # Keep only numeric columns
  numeric_df <- df %>% select(-all_of(remove))
  numeric_cols <- names(numeric_df)[sapply(numeric_df, is.numeric)]
  
  return(numeric_cols)
}

## Make sure all datasets have NO NAs
family_df  <- family_df %>% drop_na()
genus_df   <- genus_df %>% drop_na()
species_df <- species_df %>% drop_na()

## Determine feature columns ONCE per dataframe
family_features  <- get_feature_cols(family_df)
genus_features   <- get_feature_cols(genus_df)
species_features <- get_feature_cols(species_df)

# === RANDOM FOREST SECTION ===

## FAMILY (RF)
fam_split <- split_train_validation(family_df, "group")
rf_family <- train_rf(fam_split$train, "group", family_features)
evaluate_rf(rf_family, fam_split$validation, "group", family_features)

## GENUS (RF)
gen_split <- split_train_validation(genus_df, "group")
rf_genus <- train_rf(gen_split$train, "group", genus_features)
evaluate_rf(rf_genus, gen_split$validation, "group", genus_features)

## SPECIES (RF)
sp_split <- split_train_validation(species_df, "group")
rf_species <- train_rf(sp_split$train, "group", species_features)
evaluate_rf(rf_species, sp_split$validation, "group", species_features)

# === SVM SECTION  ===

## FAMILY (SVM)
svm_family <- train_svm(fam_split$train, "group", family_features)
evaluate_svm(svm_family, fam_split$validation, "group", family_features)

## GENUS (SVM)
svm_genus <- train_svm(gen_split$train, "group", genus_features)
evaluate_svm(svm_genus, gen_split$validation, "group", genus_features)

## SPECIES (SVM)
svm_species <- train_svm(sp_split$train, "group", species_features)
evaluate_svm(svm_species, sp_split$validation, "group", species_features)

## _ HEATMAP (Figure 2) ----
# Generate predictions for family model
# Generate predictions for family model
fam_preds <- predict(
  rf_family,
  fam_split$validation[, family_features]
)

# Build confusion matrix
cm <- table(
  observed = fam_split$validation$group,
  predicted = fam_preds
)

# Convert to tidy dataframe
cm_df <- as.data.frame(cm) %>%
  group_by(observed) %>%
  mutate(
    row_total = sum(Freq),
    prop = Freq / row_total
  )

# Plot
ggplot(cm_df, aes(x = predicted, y = observed, fill = prop)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(
    aes(label = paste0(Freq, "\n(", scales::percent(prop, accuracy = 1), ")")),
    size = 6, color = "black"
  ) +
  scale_fill_gradient(
    low = "#E0E5D2",
    high = "#8A9A5B",
    name = "Row % Correct"
  ) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Confusion Matrix — Family Level",
    subtitle = paste(
      "Random Forest Accuracy:",
      scales::percent(mean(fam_preds == fam_split$validation$group), accuracy = 1)
    ),
    x = "Predicted",
    y = "Observed"
  )
## _ TOP 20 (Figure 3) ----

# identify k-mer feature columns safely
feature_cols <- base::setdiff(
  colnames(family_df),
  c("title", "sequence", "unique_id", "species_name", "group")
)

# extract importance values for only k-mer features
imp <- importance(rf_family, type = 1)

imp_df <- data.frame(
  feature = rownames(imp),
  importance = imp[, 1]
)

# remove any accidental non-kmer features (safety check)
imp_df <- imp_df %>%
  filter(feature %in% feature_cols)

# Select the top 20 features (use dplyr::slice explicitly)
top20 <- imp_df %>%
  arrange(desc(importance)) %>%
  dplyr::slice(1:20)

# plot the results
ggplot(top20, aes(x = reorder(feature, importance), y = importance)) +
  geom_col(fill = "#8A9A5B") +
  coord_flip() +
  theme_minimal() +
  theme_set(theme_minimal(base_family = "Calibri", base_size = 16)) +
  labs(
    title = "Top 20 Most Informative k-mers — Family Classification",
    x = "k-mer Feature",
    y = "Importance (Mean Decrease Accuracy)"
  )

## _ PCA (Figure 4)----
# removing non-numeric columns and prepare PCA dataframe
remove_cols <- c("Title", "title", "Sequence", "sequence", "unique_id", "species_name", "genus_name", "group", "nucleotides", "nucleotides2")

feature_cols <- base::setdiff(colnames(family_df), remove_cols)

pca_df <- family_df %>%
  select(all_of(feature_cols), group) %>%
  mutate(across(all_of(feature_cols), as.numeric)) %>%
  drop_na()

# running PCA with prcomp
pca <- prcomp(pca_df[, feature_cols], center = TRUE, scale. = TRUE)

scores <- as.data.frame(pca$x)
scores$group <- pca_df$group

# plotting PCA with ggplot
ggplot(scores, aes(PC1, PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(linewidth = 1.2) +
  theme_minimal(base_size = 16) +
  labs(
    title = "PCA of rbcL k-mer Features — Family Level",
    x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"),
    color = "Family"
  ) +
  scale_color_manual(values = c(
    "Solanaceae" = "#C8A2C8",
    "Convolvulaceae" = "#8A9A5B"
  )) +
  theme(
    text = element_text(family = "Calibri", size = 16)
  )

## _ GC Content (Figure 5) ----

family_df %>%
  mutate(
    G_count = stringr::str_count(.data$Sequence, "G"),
    C_count = stringr::str_count(.data$Sequence, "C"),
    GC = (G_count + C_count) / nchar(.data$Sequence)
  ) %>%
  ggplot(aes(x = GC, fill = group)) +
  geom_density(alpha = 0.6) +
  theme_minimal(base_size = 16) +
  theme(text = element_text(family = "Calibri", size = 16)) +
  labs(
    title = "GC Content Distribution — Family Level",
    x = "GC Proportion",
    y = "Density",
    fill = "Group"
  ) +
  scale_fill_manual(values = c(
    "Solanaceae" = "#C8A2C8",
    "Convolvulaceae" = "#8A9A5B"
  ))


## _ SVM Classifier (2) ----
# since my RF training turned out... too perfect, it is recommended in the instructions to add another classifier:
train_svm <- function(train_df, label_col, feature_cols) {
  svm_model <- svm(
    x = train_df[, feature_cols],
    y = as.factor(train_df[[label_col]]),
    kernel = "radial",
    probability = TRUE
  )
  return(svm_model)
}

evaluate_svm <- function(model, validation_df, label_col, feature_cols) {
  preds <- predict(model, validation_df[, feature_cols])
  cm <- table(
    observed = validation_df[[label_col]],
    predicted = preds
  )
  return(cm)
}

svm_family <- train_svm(fam_split$train, "group", feature_cols)
evaluate_svm(svm_family, fam_split$validation, "group", feature_cols)

## _ ROC Curve (Figure 6) ----

# probabilities for RF
rf_prob <- predict(
  rf_family,
  fam_split$validation[, feature_cols],
  type = "prob"
)[, "Solanaceae"] # choose one class for ROC

# probabilities for SVM
svm_raw <- predict(
  svm_family,
  fam_split$validation[, feature_cols],
  probability = TRUE
)
svm_prob <- attr(svm_raw, "probabilities")[, "Solanaceae"]

# true labels: 1 = Solanaceae, 0 = Convolvulaceae
true_labels <- ifelse(
  fam_split$validation$group == "Solanaceae",
  1,
  0
)

# build ROCR objects
rf_pred <- ROCR::prediction(rf_prob, true_labels)
svm_pred <- ROCR::prediction(svm_prob, true_labels)

rf_perf <- ROCR::performance(rf_pred, "tpr", "fpr")
svm_perf <- ROCR::performance(svm_pred, "tpr", "fpr")

# plot ROC
plot(
  rf_perf,
  col = "#8A9A5B",
  lwd = 3,
  main = "ROC Curve — Family Level Classification"
)
plot(
  svm_perf,
  col = "#C8A2C8",
  lwd = 3,
  add = TRUE
)

# compute AUCs
rf_auc <- ROCR::performance(rf_pred, "auc")@y.values[[1]]
svm_auc <- ROCR::performance(svm_pred, "auc")@y.values[[1]]

legend(
  "bottomright",
  legend = c(
    paste0("Random Forest (AUC = ", round(rf_auc, 3), ")"),
    paste0("SVM (AUC = ", round(svm_auc, 3), ")")
  ),
  col = c("#8A9A5B", "#C8A2C8"),
  lwd = 3
)

