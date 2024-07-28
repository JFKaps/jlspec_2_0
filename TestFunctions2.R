# Load necessary packages ---------------------------
library(ggplot2)
library(dbscan)
library(fpc)
library(plotly)
library(factoextra)
library(caret)
library(ggdendro)
library(reshape2)
library(car)
library(heatmaply)
library(party)
library(rpart)
library(rpart.plot)
library(ROCR)
library(randomForest)

# Define the extract_pattern and format_strings functions
extract_pattern <- function(name) {
  pattern <- "([A-Za-z]+\\s[0-9]+:[0-9]+)|([A-Za-z]+\\s[[:alpha:]]?-?[0-9]+:[0-9]+)"
  matches <- regmatches(name, gregexpr(pattern, name))
  
  if (length(matches[[1]]) > 0) {
    return(matches[[1]][1])
  } else {
    return(name)
  }
}

format_strings <- function(input_strings) {
  formatted_strings <- gsub("\\s+", "", input_strings)
  formatted_strings <- gsub("([A-Za-z]*)(\\d+):(\\d+)", "\\1(\\2:\\3)", formatted_strings)
  return(formatted_strings)
}

# Principal Component Analysis (PCA) ---------------------------
perform_pca_analysis <- function(data, sequence) {
  PCA <- prcomp(t(data), center = TRUE)
  var <- PCA$sdev^2
  var_percentage <- var / sum(var) * 100
  cumulative_var_percentage <- cumsum(var_percentage)
  var_explained_first_two <- sum(var[1:2]) / sum(var) * 100
  cat("Percentage of variance explained by first two components:", var_explained_first_two, "%\n")
  
  PCA.df <- data.frame(
    Sample = row.names(t(data)),
    PC1 = PCA$x[, 1],
    PC2 = PCA$x[, 2],
    class = as.factor(sequence$class),
    PC = 1:length(var_percentage),
    VarianceExplained = var_percentage,
    CumulativeVarianceExplained = cumulative_var_percentage
  )
  
  Scree_plot <- ggplot(PCA.df, aes(x = PC)) +
    geom_point(aes(y = VarianceExplained, color = "Variance Explained"), size = 2) +
    geom_line(aes(y = VarianceExplained, color = "Variance Explained", linetype = "Variance Explained"), linewidth = 1) +
    geom_point(aes(y = CumulativeVarianceExplained, color = "Cumulative Variance Explained"), size = 2) +
    geom_line(aes(y = CumulativeVarianceExplained, color = "Cumulative Variance Explained", linetype = "Cumulative Variance Explained"), linewidth = 1) +
    labs(x = "Principal Component", y = "Variance Explained (%)", title = "Variance Explained by Principal Components") +
    scale_color_manual(values = c("Variance Explained" = "black", "Cumulative Variance Explained" = "red"), labels = c("Variance Explained", "Cumulative Variance Explained")) +
    scale_linetype_manual(values = c("Variance Explained" = "solid", "Cumulative Variance Explained" = "dashed"), labels = c("Variance Explained", "Cumulative Variance Explained")) +
    theme_bw() +
    theme(legend.position = "top") +
    guides(color = guide_legend(title = "Explained Type", override.aes = list(linetype = c("solid", "dashed")))) +
    scale_y_continuous(sec.axis = sec_axis(~., name = "Cumulative Variance Explained (%)"))
  
  list(Scree_plotly = ggplotly(Scree_plot), PCA.df = PCA.df)
}

create_pca_plot <- function(PCA.df, custom_colors) {
  plot_PCA <- ggplot(PCA.df, aes(x = PC1, y = PC2, sample = Sample, color = class)) +
    geom_point(size = 2) +
    labs(title = "PCA", color = "Condition", x = "PC1", y = "PC2") +
    scale_color_manual(values = custom_colors) +
    theme_bw()
  
  ggplotly(plot_PCA)
}

# Local Outlier Factor (LOF) ---------------------------
calculate_and_plot_lof <- function(PCA.df, threshold = 1) {
  PCA.df$LOF <- lof(PCA.df[, c("PC1", "PC2")], minPts = 4)
  
  LOF_plot <- ggplot(PCA.df, aes(x = reorder(Sample, LOF), y = LOF)) +
    geom_line() +
    geom_point(shape = 19) +
    labs(x = "Sample", y = "LOF", title = "Local Outlier Factor (LOF) Scores") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  LOF_plotly <- ggplotly(LOF_plot)
  
  PCA.df$Outlier_Detection <- ifelse(PCA.df$LOF > threshold, "Outlier", "Inlier")
  
  LOF_OD_plot <- ggplot(PCA.df, aes(x = PC1, y = PC2, sample = Sample, color = Outlier_Detection)) +
    geom_point() +
    scale_color_manual(values = c("Inlier" = "blue", "Outlier" = "red")) +
    labs(x = "PC1", y = "PC2", title = "Outlier Detection with LOF", color = "Outlier Detection") +
    theme_bw()
  
  ggplotly(LOF_OD_plot)
}

# DBSCAN ---------------------------
perform_dbscan <- function(PCA.df) {
  kNN_plot <- kNNdistplot(PCA.df[, 2:3], k = 4)
  abline(h = 5, col = "red", lty = 2)
  
  res <- dbscan::dbscan(as.matrix(PCA.df[, 2:3]), eps = 5, minPts = 4)
  
  plot_data <- data.frame(Sample = rownames(PCA.df), 
                          PC1 = PCA.df$PC1,
                          PC2 = PCA.df$PC2,
                          Cluster = factor(res$cluster),
                          Outlier = ifelse(res$cluster == 0, "Outlier", "Inlier"))
  
  DBSCAN_plot <- ggplot(plot_data, aes(text1 = Sample,
                                       x = PC1,
                                       y = PC2,
                                       color = Cluster,
                                       text2 = Outlier)) +
    geom_point() +
    labs(title = "DBSCAN Clustering - Eps 5 - MinPts 4",
         x = "PC1",
         y = "PC2",
         color = "Cluster") +
    theme_bw()
  
  ggplotly(DBSCAN_plot)
}

# Hierarchical DBSCAN ---------------------------
perform_hdbscan <- function(PCA.df) {
  res_hdbscan <- hdbscan(dist(PCA.df[, 2:3], method = "euclidean"), minPts = 5)
  plot(res_hdbscan)
  plot(res_hdbscan, show_flat = TRUE)
  
  plot_data <- data.frame(Sample = rownames(PCA.df), PC1 = PCA.df[,2], PC2 = PCA.df[,3], Cluster = factor(res_hdbscan$cluster), OutlierScore = res_hdbscan$outlier_scores)
  
  HDBSCAN_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Cluster, text = Sample, outlierscore = OutlierScore)) +
    geom_point(pch = 19) +
    labs(title = "Hierarchical DBSCAN - MinPts 4", x = "PC1", y = "PC2", color = "Cluster", size = "Outlier Score") +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size = 4)), size = guide_legend(override.aes = list(alpha = 0.5)))
  
  ggplotly(HDBSCAN_plot)
}

# OPTICS ---------------------------
perform_optics <- function(PCA.df) {
  kNN <- kNNdist(PCA.df[,2:3], k = 3, all = FALSE)
  kNNdistplot(PCA.df[,2:3], k = 3)
  abline(h = 5, col = "red", lty = 2)
  text(x = 9, y = 6, "Knee point")
  
  opt <- dbscan::optics(PCA.df[, 2:3], minPts = 4)
  
  plot(opt,
       main = "Reachability Plot",
       ylim = c(0.55, max(opt$reachdist[!is.infinite(opt$reachdist)], na.rm = TRUE) + 1),
       ylab = "Reachability Distance",
       xlab = "Order of Points")
  
  opt_threshold <- extractDBSCAN(opt, eps_cl = 6)
  plot(opt_threshold,
       main = "Reachability Plot kt = 6",
       ylim = c(0.55, max(opt_threshold$reachdist[!is.infinite(opt_threshold$reachdist)], na.rm = TRUE) + 1))
  
  sorted_clusters <- sort(unique(opt_threshold$cluster), decreasing = TRUE)
  hullplot(PCA.df[, 2:3], opt_threshold, ylab = "PC2", xlab = "PC1", pch = 19)
  legend("topright", legend = paste("Cluster", sorted_clusters), col = sorted_clusters + 1, pch = 19)
  
  opt <- optics(PCA.df[,2:3], eps = 5, minPts = 4)
  plot(opt, main = "Reachability Plot",
       ylim = c(0.55, max(opt$reachdist[!is.infinite(opt$reachdist)], na.rm = TRUE) + 1))
  
  plot(PCA.df[opt$order, 2:3], col = "red", main = "OPTICS Traversal", ylab = "PC2", xlab = "PC1")
  polygon(PCA.df[,2:3][opt$order, ])
  
  opt1 <- extractDBSCAN(opt, eps_cl = 4)
  plot(opt1,
       main = "Reachability Plot kt = 4",
       ylim = c(0.55, max(opt$reachdist[!is.infinite(opt$reachdist)], na.rm = TRUE) + 1))
  hullplot(PCA.df[,2:3], opt1, ylab = "PC2", xlab = "PC1")
  legend("topright", legend = paste("Cluster", unique(opt1$cluster)), col = opt1$cluster + 1, pch = 19)
}

# K-means Clustering ---------------------------
perform_kmeans <- function(PCA.df, custom_colors) {
  wss <- fviz_nbclust(PCA.df[,2:3], kmeans, method = "wss", k.max = nrow(PCA.df) - 1)
  ggplotly(wss)
  sil <- fviz_nbclust(PCA.df[,2:3], kmeans, method = "silhouette", k.max = nrow(PCA.df) - 1)
  ggplotly(sil)
  gap <- fviz_nbclust(PCA.df[,2:3], kmeans, method = "gap_stat", k.max = nrow(PCA.df) - 1)
  ggplotly(gap)
  
  layers <- sil$layers
  for (layer in layers) {
    if (!is.null(layer$data) && "xintercept" %in% names(layer$data)) {
      x_intercepts <- layer$data$xintercept
      print(x_intercepts)
    }
  }
  
  k <- x_intercepts
  kmeans_res <- kmeans(PCA.df[,2:3], centers = k)
  cluster_labels <- kmeans_res$cluster[rownames(PCA.df[,2:3])]
  
  kmeans.df <- data.frame(x = PCA.df[,2], y = PCA.df[,3], cluster = factor(cluster_labels), Sample = PCA.df$Sample)
  plot_kmeans <- ggplot(kmeans.df, aes(x = x, y = y, color = cluster, sample = Sample)) +
    ggtitle("k-means Clustering") +
    geom_point(size = 2) +
    labs(color = "Clusters") +
    xlab("PC1") +
    ylab("PC2") +
    scale_color_manual(values = custom_colors) +
    theme_bw()
  
  ggplotly(plot_kmeans)
}

# Hierarchical Clustering ---------------------------
perform_hierarchical_clustering <- function(PCA.df, sequence, custom_colors) {
  TrueLabel_plot <- ggplot(PCA.df, aes(x = PC1, y = PC2, color = class, text = Sample)) +
    geom_point(shape = 16) +
    labs(title = "Dimensionally reduced data - True labels", x = "PC1", y = "PC2") +
    scale_color_discrete(name = "Class") +
    guides(color = guide_legend(title = "Class", override.aes = list(shape = 16))) +
    theme_bw()
  
  ggplotly(TrueLabel_plot)
  
  data_dist <- dist(PCA.df[, 2:3])
  data_complete <- hclust(data_dist, method = "complete")
  completetree <- cutree(data_complete, k = 3)
  completetree_factor <- factor(completetree, levels = levels(as.factor(sequence$class)))
  con_mat <- confusionMatrix(completetree_factor, as.factor(sequence$class))
  print(con_mat)
  
  plt <- as.data.frame(con_mat$table)
  plt$Prediction <- factor(plt$Prediction, levels = rev(levels(plt$Prediction)))
  
  con_mat_comp_plot <- ggplot(plt, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) +
    scale_fill_gradient(low = "white", high = "firebrick") +
    labs(x = "Reference", y = "Prediction") +
    scale_x_discrete(labels = levels(completetree_factor)) +
    scale_y_discrete(labels = rev(levels(as.factor(sequence$class)))) +
    theme_bw() +
    ggtitle("Confusion Matrix - Complete")
  
  ggplotly(con_mat_comp_plot)
  
  dend_data_comp <- as.dendrogram(data_complete)
  dend_data <- dendro_data(dend_data_comp, type = "rectangle")
  
  dend_plot_comp <- ggplot(dend_data$segments) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = dend_data$labels, aes(x, y, label = label), hjust = 1, angle = 90, size = 3) +
    ylim(-4, max(dend_data$segments[, 4] + 1)) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title = "Dendrogram - Complete")
  
  dend_plot_comp
  
  DimR_ComL_plot <- ggplot(PCA.df, aes(x = PC1, y = PC2, sample = Sample, color = completetree_factor)) +
    geom_point(shape = 16) +
    labs(title = "Dimensionally reduced data - Complete linkage", x = "PC1", y = "PC2") +
    scale_color_discrete(name = "Cluster") +
    theme_bw()
  
  ggplotly(DimR_ComL_plot)
}
