# PCA Functions for Compositional Data Analysis
# These functions are used for Principal Component Analysis on compositional microbiome data
# Functions handle zero replacement, CLR transformation, and visualization

# Function to remove taxa/samples with low prevalence
# Input: A wide df with samples as rows and taxa as columns
# Result: Remove taxa that are in too few samples, or samples that have too few taxa
remove_samples_or_taxa <- function(df, min_samples_per_taxon, min_taxa_per_sample){
  # Store & print original df dimensions
  original_rows <- nrow(df)
  original_cols <- ncol(df)
  cat("Initial df: Rows (samples):", original_rows, ", Columns (taxa):", original_cols, "\n")

  # Remove taxa that are in less than x samples
  df <- df %>% select(where(~ sum(. != 0) > min_samples_per_taxon))
  
  # Remove samples that contain less than x taxa
  df <- df %>% filter(rowSums(. != 0) > min_taxa_per_sample)
  
  df <- df %>% select(where(~ any(. != 0)) | where(is.character) | where(is.factor))  

  removed_rows <- original_rows - nrow(df)
  removed_cols <- original_cols - ncol(df)
  cat("Removed: Rows (samples):", removed_rows, ", Columns (taxa):", removed_cols, "\n")

  cat("Resulting df: Rows (samples):", nrow(df), ", Columns (taxa):", ncol(df), "\n")
  
  return(df)
}

# CLR (Centered Log-Ratio) transformation function
clr_transform <- function(x) {
  log(x) - mean(log(x), na.rm = TRUE)
}

# Perform PCA on compositional data
# Handles zero replacement, centering, scaling, and CLR transformation
perform_pca <- function(df, zero_method = "GBM", z_delete = TRUE) {
  # Store original dimensions
  original_rows <- nrow(df)
  original_cols <- ncol(df)

  # 1. Zero replacement
  if (any(df == 0)) { # cmultRepl already checks for this
    print("Zeros found")
    df <- cmultRepl(df,
      method = zero_method, output = "prop",
      z.warning = 0.8, z.delete = z_delete
    )
    df <- df * 100
  }

  # Print removed rows and columns
  removed_rows <- original_rows - nrow(df)
  removed_cols <- original_cols - ncol(df)
  cat("Rows (samples) removed after zero replacement:", removed_rows, "\n")
  cat("Columns (taxa) removed after zero replacement:", removed_cols, "\n")

  # Geometric mean function
  geometric_mean <- function(x) {
    # Use log to avoid underflow
    exp(mean(log(x), na.rm = TRUE))
  }

  # 2. Calculate geometric mean of the parts (taxa) of the data set.
  taxa_geometric_means <- apply(df, 2, geometric_mean)

  # 3. Center data
  df_centered <- sweep(df, 2, taxa_geometric_means, FUN = "/")

  df_centered <- as.matrix(df_centered)

  # Compute the Variation Matrix
  variation_matrix <- outer(
    1:ncol(df_centered), 1:ncol(df_centered),
    Vectorize(function(i, j) var(log(df_centered[, i] / df_centered[, j]), na.rm = TRUE))
  )

  # Calculate Total Variance
  D <- ncol(df_centered) # Number of taxa (columns)
  totvar <- (1 / (2 * D)) * sum(variation_matrix, na.rm = TRUE)

  # 4. Scale data
  power_exponent <- 1 / sqrt(totvar)
  df_scaled <- df_centered^power_exponent

  # CLR transform data
  df_clr <- as.data.frame(t(apply(df_scaled, 1, clr_transform)))

  df_clr_dist <- as.data.frame(t(apply(df, 1, clr_transform)))

  # Perform PCA on zero replaced, centered, scaled, and CLR transformed df
  pca_result <- prcomp(df_clr, center = FALSE, scale. = FALSE)

  return(list(
    df_clr = df_clr,
    df_clr_dist = df_clr_dist,
    pca_result = pca_result
  ))
}

# Plot PCA ordination
# Use the 'pca_result' df produced from 'perform_pca' function to make the PCA plot
plot_pca <- function(df, 
                     samples_color_metadata, samples_shape_metadata, 
                     samples_color_value, loadings_color_metadata, 
                     loadings_color_value, loadings_taxon_level,
                     sample_metadata, genome_metadata, order_colors, 
                     custom_ggplot_theme, scaling_factor_value = 1.5, 
                     loadings_number = 10, show_labels = FALSE, add_arrows = FALSE, add_centroids = FALSE) {
  
  # Extract scores from PCA results
  scores <- rownames_to_column(as.data.frame(df$x), var = "microsample")
  scores <- left_join(scores, sample_metadata, by = join_by(microsample == microsample))
  
  # Calculate limits for x and y axes
  x_limit <- max(abs(scores$PC1))
  y_limit <- max(abs(scores$PC2))
  
  # Calculate variance explained by each PC (principal component) & create labels for plot
  variance_explained <- (df$sdev^2) / sum(df$sdev^2) * 100
  pc1_label <- paste0("PC1: ", round(variance_explained[1], 2), "% variance explained")
  pc2_label <- paste0("PC2: ", round(variance_explained[2], 2), "% variance explained")
  
  # Set a scaling factor for loadings (arrows)
  scaling_factor <- scaling_factor_value
  
  # Extract and scale the loadings
  loadings <- df$rotation[, 1:2] %>%
    as.data.frame() %>%
    mutate(genome = rownames(.)) %>%  
    mutate(PC1 = PC1 * scaling_factor,
           PC2 = PC2 * scaling_factor) %>%
    left_join(genome_metadata, by = join_by(genome == genome)) %>%  
    mutate(abs_loading = sqrt(PC1^2 + PC2^2)) %>%  
    arrange(desc(abs_loading)) %>%
    slice_max(order_by = abs_loading, n = loadings_number) %>%
    mutate(order_color = order_colors[order])
  
  # Create ggplot
  p <- ggplot() +
    # Plot the samples (points)
    geom_point(data = scores, 
               aes(x = PC1, y = PC2, 
                   fill = .data[[samples_color_metadata]],
                   shape = .data[[samples_shape_metadata]]), #.data[[]] tells ggplot2 to look up the column dynamically.
               size = 2.5, alpha = 0.8,
               color = "black", stroke = 0.3) +
    
    scale_fill_manual(values = samples_color_value) + 
    scale_shape_manual(values = c(21, 24, 23, 22, 25)) +
    new_scale_color() + #  new_scale_fill() ?
    # Plot the loadings (taxa, i.e. arrows)
    geom_segment(data = loadings, 
                 aes(x = 0, y = 0, xend = PC1, yend = PC2, color = .data[[loadings_color_metadata]]),
                 arrow = arrow(length = unit(0.2, "cm")), 
                 size = 0.7, alpha = 0.9) +
    scale_color_manual(name = "Classification", values = loadings_color_value) +
    geom_text_repel(data = loadings, 
                    aes(x = PC1, y = PC2, label = .data[[loadings_taxon_level]]),
                    color = "black", size = 3, vjust = -0.5, alpha=0.7, max.overlaps = 20) +
    labs(title = "PCA Ordination Plot",
         x = pc1_label,
         y = pc2_label) +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    scale_y_continuous(limits = c(-y_limit, y_limit)) +
    geom_hline(yintercept = 0, color = "darkgrey") +
    geom_vline(xintercept = 0, color = "darkgrey") +
    theme_minimal() +
    custom_ggplot_theme +
    # coord_fixed(ratio = 1) +
    guides(fill = guide_legend(override.aes = list(shape = 21, color = "black"))) # to print the legend for colour correctly!
  
  
      # ADD text to points
     if (show_labels) {
        p <- p + 
          geom_text(
            data = scores, 
            aes(x = PC1, y = PC2, label = animal),
            position = position_nudge(y = -0.02),
            size = 3.5, alpha = 0.7
          )
     }
  
  if (add_arrows) {
  # Prepare arrow data
  arrows_df <- scores %>%
    arrange(treatment, age) %>%
    group_by(treatment) %>%
    mutate(
      PC1_next = lead(PC1),
      PC2_next = lead(PC2)
    ) %>%
    filter(!is.na(PC1_next))
  
  # Add arrows to the plot
  p <- p + 
    geom_segment(
      data = arrows_df,
      aes(
        x = PC1, y = PC2,
        xend = PC1_next, yend = PC2_next,
        color = treatment
      ),
      arrow = arrow(length = unit(0.15, "cm")),
      linewidth = 0.6,
      alpha = 0.8,
      inherit.aes = FALSE
    )
  }
  
if (add_centroids) {

  # Compute centroids per treatment × age
  centroids <- scores %>%
    group_by(.data[[samples_color_metadata]], .data[[samples_shape_metadata]]) %>%
    summarise(
      PC1 = mean(PC1),
      PC2 = mean(PC2),
      .groups = "drop"
    ) %>%
    rename(
      color_group = 1,
      shape_group = 2
    ) %>%
    arrange(color_group, shape_group)

  centroid_arrows <- centroids %>%
    group_by(color_group) %>%
    mutate(
      PC1_next = lead(PC1),
      PC2_next = lead(PC2)
    ) %>%
    filter(!is.na(PC1_next))

  # Add new color scale so centroids use treatment palette again
  p <- p + 
    new_scale_color() +
    new_scale_fill() +

    geom_segment(
      data = centroid_arrows,
      aes(
        x = PC1, y = PC2,
        xend = PC1_next, yend = PC2_next,
        color = color_group
      ),
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      linewidth = 1.2,
      alpha = 0.7,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = centroids,
      aes(
        x = PC1, y = PC2,
        shape = factor(shape_group),
        color = color_group,
        fill  = color_group
      ),
      size = 4.5,
      stroke = 1,
      alpha = 0.5,
      fill = "white",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = centroids,
      aes(
        x = PC1, y = PC2,
        label = shape_group,
        color = color_group
      ),
      size = 3.3,
      vjust = -1.2,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = samples_color_value) +
    scale_fill_manual(values = samples_color_value)
}
  
  return(p)
}

# Plot scree plot for PCA
plot_scree <- function(pca_result_list) {
  # Extract variance explained by each principal component
  variance_explained <- (pca_result_list$pca_result$sdev^2) / sum(pca_result_list$pca_result$sdev^2) * 100
  
  # Create data frame for plotting
  scree_data <- data.frame(
    PC = paste0("PC", 1:length(variance_explained)),
    Variance = variance_explained,
    Cumulative = cumsum(variance_explained)
  )
  
  # Create scree plot
  p <- ggplot(scree_data %>% slice_head(n = 10), aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = Cumulative, group = 1), color = "red", linewidth = 1) +
    geom_point(aes(y = Cumulative), color = "red", size = 2) +
    labs(
      x = "Principal Component",
      y = "Variance Explained (%)",
      title = "PCA Scree Plot"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

