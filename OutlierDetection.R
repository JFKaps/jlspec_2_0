# Load neccesary packages  ----
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


# Load metadata set and merged sequence file ----
# Eva data 
setwd("~/Desktop/SDU/Cand/2.Sem/ISA/Jakob/Bsc-R-scripts--main/Data working")
sequence <- read.csv("Eva sequence_lipidomics_pos.csv")
data <- read.csv("Eva-pos-export-from-profinder.csv")

head(sequence)
head(data)

data <- read.csv("OG lab_data.csv")
sequence <- read.csv("OG seq_OBS_manipuleret.csv")

head(sequence)
head(data)

data <- read.csv("Signe lipid mLiver jackob_new.csv")
sequence <- read.csv("Signe Mliver liopid seq jacob1.csv")

head(sequence)
head(data)

data <- read.csv("helle fix.csv")
sequence <- read.csv("helle seq fix.csv")

head(sequence)
head(data)

#### Data Cleaning ---- 
# Subset data and sequence to keep only 'Name' and 'Sample' columns
data <- data[, sequence[ , 'labels'] %in% c("Name","Sample")]
sequence <- sequence[sequence[ , 'labels'] %in% c("Name","Sample"), ]


# Check if sequence is loaded successfully
if (is.null(sequence)) {
  return(NULL)  # Stop the observeEvent if no file is uploaded
}

# Capture number of rows before filtering
number_of_rows_before <- nrow(data)

# Clean compound names using custom functions
data[, 1] <- sapply(data[, 1], extract_pattern)
data[, 1] <- sapply(data[, 1], format_strings)

# Filter rows based on specific pattern
# data <- filter_data_by_pattern(data)

# Make unique row names and perform additional data transformations
unique_row_names <- make.unique(as.character(data[,"Compound.Name"]))
rownames(data) <- unique_row_names

# Remove the Compound.Name column as it is now used for row names
data <- data[,-1]
# Perform logarithmic transformation (e.g., log2) on the  data
data <- log2(data)
# Remove rows with NA (missing) values
data <- na.omit(data)
# Scale the data (standardization) for further analysis
data_scaled <- as.data.frame(scale(data))
# Transpose the data matrix for easier manipulation (columns as variables)
data_T <- as.data.frame(t(data))
# Transpose the scaled data matrix for easier manipulation (columns as variables)
data_scaled_T <- as.data.frame(t(data_scaled))
# Remove the first row of 'sequence'
sequence <- sequence[-1,]

#### Visualization of data ----
# Define custom colors
custom_colors <- c(
  "1" = '#FF5733',  # Brighter Red
  "2" = '#33B5E5',  # Brighter Blue
  "3" = '#66BB6A',  # Brighter Green
  "4" = '#BA68C8',  # Brighter Purple
  "5" = '#FFA500',  # Brighter Orange
  "6" = '#FFFF66'   # Brighter Yellow
)
# Define functions for data visualization
Visualize_histograms <- function(data, sequence) {
  par(mfrow = c(2, ceiling(ncol(data) / 4)))
  for (i in 1:ncol(data)) {
    hist(data[, i],
         main = colnames(data)[i],
         xlab = "",
         ylab = "Expression")
    abline(v = mean(data[,i]),
           col = "red")
    abline(v = median(data[,i]),
           col = "blue")
    abline(v = mean(data[,i]) + sd(data[,i]),
           col = "green")
    abline(v = mean(data[,i]) - sd(data[,i]),
           col = "green")
  }
  par(mfrow = c(1, 1))
}

Visualize_boxplot <- function(data, sequence) {
  # Get unique class labels from the sequence$class column
  class_labels <- unique(na.omit(sequence$class))
  
  # Assign colors based on class labels using custom_colors mapping
  colors <- custom_colors[class_labels]
  
  # Plot the boxplot with custom colors based on class
  boxplot(data,
          las = 2,
          main = "Boxplot",
          col = colors[sequence$class],  # Map class labels to custom colors
          ylim = c(min(data, na.rm = TRUE), max(data, na.rm = TRUE)))
}

# Perform data visualization
Visualize_histograms(data, sequence)
Visualize_boxplot(data,sequence)

#### Principal Component Analysis (PCA) ---- 
perform_pca_analysis <- function(data, sequence) {
  # Perform PCA
  PCA <- prcomp(t(data), center = TRUE)
  
  # Calculate the variance explained by each component
  var <- PCA$sdev^2
  
  # Calculate the percentage of variance explained by each component
  var_percentage <- var / sum(var) * 100
  
  # Calculate cumulative variance explained
  cumulative_var_percentage <- cumsum(var_percentage)
  
  # Calculate the percentage of variance explained by the first two components
  var_explained_first_two <- (sum(var[1:2]) / sum(var)) * 100
  cat("Percentage of variance explained by first two components:", var_explained_first_two, "%\n")
  
  # Create a data frame for plotting PCA
  PCA.df <- data.frame(
    Sample = row.names(t(data)),
    PC1 = PCA$x[, 1],
    PC2 = PCA$x[, 2],
    class = as.factor(sequence$class),
    PC = 1:length(var_percentage), 
    VarianceExplained = var_percentage,
    CumulativeVarianceExplained = cumulative_var_percentage
  )
  
  # Get the PCs
  PCs <- PCA$x[, 1:2]
  
  # Create the PCA variance explained plot
  Scree_plot <- ggplot(PCA.df, aes(x = PC)) +
    # Plot individual variance explained
    geom_point(aes(y = VarianceExplained,
                   color = "Variance Explained"), size = 2) +
    geom_line(aes(y = VarianceExplained,
                  color = "Variance Explained",
                  linetype = "Variance Explained"), linewidth = 1) +
    # Plot cumulative variance explained
    geom_point(aes(y = CumulativeVarianceExplained,
                   color = "Cumulative Variance Explained"), size = 2) +
    geom_line(aes(y = CumulativeVarianceExplained,
                  color = "Cumulative Variance Explained",
                  linetype = "Cumulative Variance Explained"), linewidth = 1) +
    # Customize labels and title
    labs(
      x = "Principal Component",
      y = "Variance Explained (%)",
      title = "Variance Explained by Principal Components"
    ) +
    # Customize colors and linetypes
    scale_color_manual(
      values = c("Variance Explained" = "black",
                 "Cumulative Variance Explained" = "red"),
      labels = c("Variance Explained",
                 "Cumulative Variance Explained")
    ) +
    scale_linetype_manual(
      values = c("Variance Explained" = "solid",
                 "Cumulative Variance Explained" = "dashed"),
      labels = c("Variance Explained",
                 "Cumulative Variance Explained")
    ) +
    # Apply bw theme and adjust legend position
    theme_bw() +
    theme(
      legend.position = "top"
    ) +
    # Customize legend labels and linetypes
    guides(
      color = guide_legend(
        title = "Explained Type",
        override.aes = list(
          linetype = c("solid", "dashed")
        )
      )
    ) +
    # Add secondary y-axis for cumulative variance explained
    scale_y_continuous(
      sec.axis = sec_axis(
        ~., 
        name = "Cumulative Variance Explained (%)"
      )
    )
  
  # Convert ggplot to plotly
  Scree_plotly <- ggplotly(Scree_plot)
  
  return(list(Scree_plotly = Scree_plotly,
              PCA.df = PCA.df))
}

# Function to create the PCA plot
create_pca_plot <- function(PCA.df, custom_colors) {
  plot_PCA <- ggplot(PCA.df, aes(x = PC1,
                                 y = PC2,
                                 sample = Sample,
                                 color = class)) +
    # Customize plot aesthetics
    geom_point(size = 2) +  # Scatter plot of PCA points
    # Customize plot labels and title
    labs(
      title = "PCA",
      color = "Condition",  # Legend label for color
      x = "PC1",             # X-axis label
      y = "PC2"              # Y-axis label
    ) +
    # Set manual color scale based on custom_colors
    scale_color_manual(values = custom_colors) +
    # Use a black and white theme
    theme_bw()
  
  # Convert ggplot to plotly
  plot_PCAly <- ggplotly(plot_PCA)
  
  return(plot_PCAly)
}
# Example usage:
result <- perform_pca_analysis(data, sequence)
scree_plotly <- result$Scree_plotly  # To display the scree plot
scree_plotly
pca_plotly <- create_pca_plot(result$PCA.df, custom_colors)  # To display the PCA plot
pca_plotly

result$PCA.df
#### Clustering-Based Methods ----
#### K-means clustering ---- 
# k-means clustering
# for k-means we need to set a number of clusters k. To find a number we try three different ways of doing so (notice we now use the objects PCs we created based on our PCA) 
wss <- fviz_nbclust(result$PCA.df[,2:3],
                    kmeans,
                    method = "wss",
                    k.max = nrow(result$PCA.df)-1)# 1st one is within-cluster sum of squares, the optimal one is where the elbow lies
ggplotly(wss) # show the plot
sil <- fviz_nbclust(result$PCA.df[,2:3],
                    kmeans,
                    method = "silhouette",
                    k.max = nrow(result$PCA.df)-1) # 2nd is silhouette score. This actually shows you the optimal cluster number
ggplotly(sil) # show the plot
gap <- fviz_nbclust(result$PCA.df[,2:3],
                    kmeans,
                    method = "gap_stat",
                    k.max = nrow(result$PCA.df)-1) # 3rd is gap statistiscs. This also shows you the optimal cluster number
ggplotly(gap) # show the plot

# Extract layers from the ggplot object
layers <- sil$layers

# Iterate over layers to find xintercept values
for (layer in layers) {
  if (!is.null(layer$data) && "xintercept" %in% names(layer$data)) {
    x_intercepts <- layer$data$xintercept
    # Now 'x_intercepts' contains the numeric values associated with xintercept
    print(x_intercepts)
  }
}


k = x_intercepts# set k = optimal value from plots
kmeans <- kmeans(result$PCA.df[,2:3], centers = k) # this line performs the k-means clustering
cluster_labels <- kmeans$cluster[rownames(result$PCA.df[,2:3])] # this matches our clustering with our PCA

kmeans.df <- data.frame(x = result$PCA.df[,2],
                        y = result$PCA.df[,3],
                        cluster = factor(cluster_labels),
                        Sample = plot_data$Sample) # as before we create a dataframe for plotting our clustering
plot_kmeans = ggplot(kmeans.df,
                     aes(x = x,
                         y = y,
                         color = cluster,
                         sample = Sample)) + # again we use ggplot to customize our plot
  ggtitle("k-means Clustering") +
  geom_point(size = 2) +
  labs(color = "Clusters") +
  xlab("PC1") +
  ylab("PC2") +
  scale_color_manual(values = custom_colors) +
  theme_bw()
ggplotly(plot_kmeans)

#### Hierarchical clustering ---- 
#### Density-Based Methods ----
#### Local Outlier Factor (LOF) ---- 
calculate_and_plot_lof <- function(PCA.df) {
  # Calculate LOF scores
  result$LOF <- lof(result$PCA.df[, c("PC1", "PC2")], minPts = 4)
  
  # Plot LOF scores using ggplot
  LOF_plot <- ggplot(PCA.df, aes(x = reorder(Sample, LOF), y = LOF)) +
    geom_line() +
    geom_point(shape = 19) +
    labs(
      x = "Sample",
      y = "LOF",
      title = "Local Outlier Factor (LOF) Scores"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels vertically
  
  # Convert ggplot to plotly
  LOF_plotly <- ggplotly(LOF_plot)
  
  return(list(LOF_plotly = LOF_plotly,
              result = result))
}

calculate_and_plot_lof(result$PCA.df)

threshold <- 1.5 

# Calculate outlier detection based on LOF scores and threshold
PCA.df$`Outlier Detection` <- ifelse(PCA.df$LOF > threshold, "Outlier", "Inlier")

# Plot using ggplot to visualize Outlier Detection with LOF
LOF_OD_plot <- ggplot(PCA.df, aes(x = x,
                                  y = y,
                                  sample = Sample,
                                  lof = LOF,
                                  color = `Outlier Detection`)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red"), labels = c("Outlier", "Inlier")) +
  theme_bw() +
  labs(
    x = "PC1",  # Set x-axis label to "PC1"
    y = "PC2",  # Set y-axis label to "PC2"
    title = "Outlier Detection with LOF",
    color = "Outlier Detection"
  )

# Plot LOF outliers using plotly
ggplotly(LOF_OD_plot)


#### Density-Based Spatial Clustering of Applications with Noise (DBSCAN) ---- 
kNN_plot <- kNNdistplot(PCs, minPts = 4)
abline(h = 5, col = "red", lty = 2)

res <- dbscan(as.matrix(result$PCA.df[,2:3]), eps = 5, MinPts = 4)

# Combine PCA results and cluster assignments into a data frame
plot_data <- data.frame(Sample = rownames(result$PCA.df), 
                        PC1 = result$PCA.df[,2],
                        PC2 = result$PCA.df[,3],
                        Cluster = factor(res$cluster),
                        Corepoint = res$isseed)

# Create a ggplot scatter plot
DBSCAN_plot <- ggplot(plot_data, aes(x = PC1,
                                     y = PC2,
                                     color = Cluster,
                                     sample = Sample,
                                     core = Corepoint)) +
  geom_point(shape = 19) +
  labs(title = "DBSCAN Clustering - Eps 3.5 - MinPts 4",
       x = "PC1",
       y = "PC2",
       color = "Cluster") +
  theme_bw()
ggplotly(DBSCAN_plot)

#### Hierarchical DBSCAN ----
res <- hdbscan(dist(PCs, method = "euclidean"), minPts = 4)
print(res)

plot(res)
plot(res, show_flat = TRUE)

# Extract the cluster membership and outlier scores
clusters <- res$cluster
outlier_scores <- res$outlier_scores

# Combine PCA results and cluster assignments into a data frame
plot_data <- data.frame(Sample = rownames(PCA.df), 
                        PC1 = PCA.df[,2],
                        PC2 = PCA.df[,3],
                        Cluster = factor(res$cluster),
                        OutlierScore = outlier_scores)

# Create a ggplot scatter plot for Euclidean distance clustering
HDBSCAN_plot <- ggplot(plot_data, aes(x = PC1,
                                      y = PC2,
                                      color = Cluster,
                                      text = Sample,
                                      outlierscore = OutlierScore)) +  # use 'text' instead of 'sample' for tooltip
  geom_point(pch = 19) +
  labs(title = "Hierarchical DBSCAN - MinPts 4",
       x = "PC1",
       y = "PC2",
       color = "Cluster",
       size = "Outlier Score") +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 4)),
         size = guide_legend(override.aes = list(alpha = 0.5))) # fix legend point size

# Convert to ggplotly
ggplotly(HDBSCAN_plot)


#### OPTICS ---- 
# kNN distances 
kNN <- kNNdist(result$PCA.df[,2:3], k = 3, all = FALSE)
kNN

kNNdistplot(result$PCA.df[,2:3], k = 3)
abline(h = 5, col = "red", lty = 2)
text(x = 9, y = 6, "Knee point")

# OPTICS 
opt <- optics(result$PCA.df[,2:3], minPts = 4)
plot(opt,
     main = "Reachability Plot",
     ylim = c(min(opt$reachdist[!is.infinite(opt$reachdist)]),
              max(opt$reachdist[!is.infinite(opt$reachdist)]+1)))
# OPTICS with cutoff 
opt1 <- extractDBSCAN(opt, eps_cl = 6)
opt1
plot(opt1,
     main = "Reachability Plot kt = 6",
     ylim = c(min(opt$reachdist[!is.infinite(opt$reachdist)]),
              max(opt$reachdist[!is.infinite(opt$reachdist)]+1)))
hullplot(result$PCA.df[,2:3], opt1,
         ylab = "PC2",
         xlab = "PC1",
         pch = 19)
legend("topright",
       legend = paste("Cluster", unique(opt1$cluster)),
       col = opt1$cluster+1,
       pch = 19)


# OPTICS with specified eps
opt <- optics(result$PCA.df[,2:3], eps = 5, minPts = 4)
plot(opt,
     main = "Reachability Plot",
     ylim = c(min(opt$reachdist[!is.infinite(opt$reachdist)]),
              max(opt$reachdist[!is.infinite(opt$reachdist)]+1)))
plot(result$PCA.df[opt$order, 2:3], 
     col = "red",
     main = "OPTICS Traversal",
     ylab = "PC2",
     xlab = "PC1")
polygon(result$PCA.df[,2:3][opt$order,])

# OPTICS with cutoff 
opt1 <- extractDBSCAN(opt, eps_cl = 4)
opt1
plot(opt1,
     main = "Reachability Plot kt = 4",
     ylim = c(min(opt$reachdist[!is.infinite(opt$reachdist)]),
              max(opt$reachdist[!is.infinite(opt$reachdist)]+1)))
hullplot(result$PCA.df[,2:3], opt1,
         ylab = "PC2",
         xlab = "PC1")
legend("topright",
       legend = paste("Cluster", unique(opt1$cluster)),
       col = opt1$cluster+1,
       pch = 19)

#### Complete ---- 
# Create a ggplot object to visualize dimensionally reduced data with true class labels
TrueLabel_plot <- ggplot(data = PCA.df,
                         aes(x = x,
                             y = y,
                             sample = Sample,
                             color = class)) +
  geom_point(shape = 16) +  # Add points with shape 16 (solid circle)
  labs(
    title = "Dimensionally reduced data - True labels",
    x = "PC1",
    y = "PC2"
  ) +  # Customize plot labels and title
  scale_color_discrete(name = "Class") +  # Customize point color legend
  guides(color = guide_legend(
    title = "Class",
    override.aes = list(shape = 16)
  )) +  # Add a manual legend with points shaped as solid circles
  theme_bw()  # Apply a black and white theme
# Convert ggplot object to Plotly for interactive visualization
ggplotly(TrueLabel_plot)

# Calculate Euclidean distance between data points based on first two principal components (PC1 and PC2)
data_dist <- dist(PCA.df[, 2:3])
# Perform hierarchical clustering using complete linkage method on the distance matrix
data_complete <- hclust(data_dist, method = "complete")
# Cut the hierarchical clustering tree into 3 clusters
completetree <- cutree(data_complete, k = 3)
# Convert the cluster labels to a factor with the same levels as the original class labels
completetree_factor <- factor(completetree,
                              levels = levels(as.factor(sequence$class)))
# Compute confusion matrix between predicted clusters and true class labels
con_mat <- confusionMatrix(completetree_factor,
                           as.factor(sequence$class))
# Print the confusion matrix
print(con_mat)

# Create a data frame from confusion matrix
plt <- as.data.frame(con_mat$table)
plt$Prediction <- factor(plt$Prediction, levels = rev(levels(plt$Prediction)))

# Plot confusion matrix as a heatmap
con_mat_comp_plot <- ggplot(plt,
                            aes(x = Reference,
                                y = Prediction,
                                fill = Freq)) +
  geom_tile() +  # Add heatmap tiles
  geom_text(aes(label = Freq)) +  # Add text labels
  scale_fill_gradient(low = "white", high = "firebrick") +  # Gradient fill colors
  labs(x = "Reference", y = "Prediction") +  # Customize axis labels
  scale_x_discrete(labels = c(levels(completetree_factor))) +  # Set x-axis labels based on factor levels
  scale_y_discrete(labels = c(rev(levels(as.factor(sequence$class))))) +  # Set y-axis labels in reverse order
  theme_bw() +  # Apply a black and white theme
  ggtitle("Confusion Matrix - Complete")  # Add a title to the plot
# Convert ggplot object to Plotly for interactive visualization
ggplotly(con_mat_comp_plot)

# Convert data to dendrogram format
dend_data_comp <- as.dendrogram(data_complete)
# Extract dendrogram data for plotting
dend_data <- dendro_data(dend_data_comp, type = "rectangle")
# Plotting the dendrogram without x-axis
dend_plot_comp <- ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +  # Add segments for dendrogram branches
  geom_text(data = dend_data$labels, aes(x, y, label = label),  # Add labels to dendrogram
            hjust = 1, angle = 90, size = 3) +  # Align labels to the right with 90-degree rotation (horizontal)
  ylim(-3, max(dend_data$segments[,4] + 1)) +  # Set y-axis limits
  theme_bw() +  # Apply a black and white theme
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank()) +  # Remove x-axis ticks 
  #scale_color_manual(values = c(dend_data$labels)) +  # Manually specify colors based on factor levels
  labs(title = "Dendrogram - Complete")  # Add a title to the plot
# Display the dendrogram plot
dend_plot_comp


# Create a ggplot object to visualize dimensionally reduced data with complete linkage labels
DimR_ComL_plot <- ggplot(data = PCA.df,
                         aes(x = x,
                             y = y,
                             sample = Sample,
                             color = completetree_factor)) +
  geom_point(shape = 16) +  # Add points with shape 16 (solid circle)
  labs(
    title = "Dimensionally reduced data - Complete linkage",
    x = "PC1", y = "PC2"
  ) +  # Customize plot labels and title
  scale_color_discrete(name = "Cluster") +  # Customize point color legend
  theme_bw()  # Apply a black and white theme
# Convert ggplot object to Plotly for interactive visualization
ggplotly(DimR_ComL_plot)


#### single ---- 
data_single <- hclust(data_dist, method = "single")
singletree <- cutree(data_single, k = 3)

singletree_factor <- factor(singletree, levels = levels(as.factor(sequence$class)))

con_mat <- confusionMatrix(singletree_factor, as.factor(sequence$class))
print(con_mat)

plt <- as.data.frame(con_mat$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

con_mat_sing_plot <- ggplot(plt, aes(Reference,Prediction,  fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="firebrick") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels=c(levels(singletree_factor))) +
  scale_y_discrete(labels=c(rev(levels(as.factor(sequence$class))))) +
  theme_bw() +
  ggtitle("Confusion Matrix - Single")
ggplotly(con_mat_sing_plot)

# Convert data to dendrogram format
dend_data_sing <- as.dendrogram(data_single)

# Extract dendrogram data for plotting
dend_data_single <- dendro_data(dend_data_sing, type = "rectangle")

# Plotting the dendrogram without x-axis
dend_plot_single <- ggplot(dend_data_single$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +  # Add segments for dendrogram branches
  geom_text(data = dend_data_single$labels, aes(x, y, label = label),  # Add labels to dendrogram
            hjust = 1, angle = 90, size = 3) +  # Align labels to the right with 90-degree rotation (horizontal)
  ylim(-2, max(dend_data_single$segments[,4] + 1)) +  # Set y-axis limits
  theme_bw() +  # Apply a black and white theme
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank()) +  # Remove x-axis ticks
  labs(title = "Dendrogram - Single")  # Add a title to the plot

# Display the dendrogram plot
dend_plot_single
# ggplotly(dend_plot_single)

DimR_SinL_plot <- ggplot(data = PCA.df, aes(x = x, y = y, color = singletree_factor)) +
  geom_point(shape = 16) +
  labs(title = "Dimensionally reduced data - Single linkage",
       x = "PC1", y = "PC2") +
  scale_color_discrete(name = "Cluster") +
  theme_bw()
ggplotly(DimR_SinL_plot)

#### Average ---- 
data_average <- hclust(data_dist, method = "average")
averagetree <- cutree(data_average, k = 3)

averagetree_factor <- factor(averagetree, levels = levels(as.factor(sequence$class)))

con_mat <- confusionMatrix(averagetree_factor, as.factor(sequence$class))
print(con_mat)

plt <- as.data.frame(con_mat$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

con_mat_aver_plot <- ggplot(plt, aes(Reference,Prediction,  fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="firebrick") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels=c(levels(averagetree_factor))) +
  scale_y_discrete(labels=c(rev(levels(as.factor(sequence$class))))) + 
  theme_bw() + 
  ggtitle("Confusion Matrix - Average")
ggplotly(con_mat_aver_plot)

# Convert data to dendrogram format
dend_data_aver <- as.dendrogram(data_average)

# Extract dendrogram data for plotting
dend_data_aver <- dendro_data(dend_data_aver, type = "rectangle")

# Plotting the dendrogram without x-axis
dend_plot_aver <- ggplot(dend_data_aver$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +  # Add segments for dendrogram branches
  geom_text(data = dend_data_aver$labels, aes(x, y, label = label),  # Add labels to dendrogram
            hjust = 1, angle = 90, size = 3) +  # Align labels to the right with 90-degree rotation (horizontal)
  ylim(-2, max(dend_data_aver$segments[,4] + 1)) +  # Set y-axis limits
  theme_bw() +  # Apply a black and white theme
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank()) +  # Remove x-axis ticks
  labs(title = "Dendrogram - Average")  # Add a title to the plot

# Display the dendrogram plot
dend_plot_aver
# ggplotly(dend_plot_aver)

DimR_AveL_plot <- ggplot(data = PCA.df, aes(x = x, y = y, color = averagetree_factor)) +
  geom_point(shape = 16) +
  labs(title = "Dimensionally reduced data - Average linkage",
       x = "PC1", y = "PC2") +
  scale_color_discrete(name = "Cluster") +
  theme_bw()
ggplotly(DimR_AveL_plot)

#### Ward.D2 ----
data_wardD2 <- hclust(data_dist, method = "ward.D2")
wardD2tree <- cutree(data_wardD2, k = 3)

wardD2tree_factor <- factor(wardD2tree, levels = levels(as.factor(sequence$class)))

con_mat <- confusionMatrix(wardD2tree_factor, as.factor(sequence$class))
print(con_mat)

plt <- as.data.frame(con_mat$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

con_mat_aver_plot <- ggplot(plt, aes(Reference,Prediction,  fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="firebrick") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels=c(levels(wardD2tree_factor))) +
  scale_y_discrete(labels=c(rev(levels(as.factor(sequence$class))))) + 
  theme_bw() + 
  ggtitle("Confusion Matrix - Ward D2")
ggplotly(con_mat_aver_plot)

# Convert data to dendrogram format
dend_data_wardD2 <- as.dendrogram(data_wardD2)

# Extract dendrogram data for plotting
dend_data_wardD2 <- dendro_data(dend_data_wardD2, type = "rectangle")

# Plotting the dendrogram without x-axis
dend_plot_wardD2 <- ggplot(dend_data_wardD2$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +  # Add segments for dendrogram branches
  geom_text(data = dend_data_wardD2$labels, aes(x, y, label = label),  # Add labels to dendrogram
            hjust = 1, angle = 90, size = 3) +  # Align labels to the right with 90-degree rotation (horizontal)
  ylim(-2, max(dend_data_wardD2$segments[,4] + 1)) +  # Set y-axis limits
  theme_bw() +  # Apply a black and white theme
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank()) +  # Remove x-axis ticks
  labs(title = "Dendrogram - Ward D2")  # Add a title to the plot

# Display the dendrogram plot
dend_plot_wardD2
# ggplotly(dend_plot_wardD2)

DimR_wardD2_plot <- ggplot(data = PCA.df, aes(x = x, y = y, color = wardD2tree_factor)) +
  geom_point(shape = 16) +
  labs(title = "Dimensionally reduced data - Ward D2 linkage",
       x = "PC1", y = "PC2") +
  scale_color_discrete(name = "Cluster") +
  theme_bw()
ggplotly(DimR_wardD2_plot)




#### Visualization techniques ---- 

#### Heatmap ---- 
gradient_col <- ggplot2::scale_fill_gradient2(
  low = "blue", high = "red", 
  midpoint = 0.5, limits = c(-3, 3)
)

heatmaply(
  data_scaled_T,
  k_col = 10,  # Number of features to display
  k_row = 4,   # Number of clusters for rows
  scale_fill_gradient_fun = gradient_col,  # Apply defined color gradient
  row_side_colors = PCA.df$class,  # Side bar indicating different classes
  clustering_distance_rows = "euclidean",  # Distance metric for row clustering
  clustering_distance_cols = "euclidean",  # Distance metric for column clustering
  dendrogram = "both",  # Display both row and column dendrograms
  margins = c(NA,NA,NA,NA),  # Margin settings
  main = "Interactive Heatmap",  # Main title
  scale = "none"  # No scaling applied
)

# Calculate variability (e.g., standard deviation) for each feature
feature_variability <- apply(data_scaled_T, 2, sd)

# Identify top 20 most and least changing features based on variability
top_n <- 20  # Number of top features to select

# Sort features by variability (standard deviation)
most_variable_features <- names(sort(feature_variability, decreasing = TRUE)[1:top_n])
least_variable_features <- names(sort(feature_variability)[1:top_n])

# Combine most and least variable feature names
selected_features <- c(most_variable_features, least_variable_features)

# Subset data to include only selected features
subset_data <- data_scaled_T[, selected_features]

# Create heatmap with selected features
heatmaply(
  subset_data,
  k_col = length(selected_features),  # Number of selected features
  k_row = 4,  # Number of clusters for rows
  scale_fill_gradient_fun = gradient_col,  # Apply defined color gradient
  #  row_side_colors = PCA.df$class,  # Side bar indicating different classes
  dendrogram = "both",  # Display both row and column dendrograms
  margins = c(NA, NA, NA, NA),  # Margin settings
  main = "Interactive Heatmap with Top Variable Features",  # Main title
  scale = "none"  # No scaling applied
)

#### Volcanoplot ---- 
#### Perform statistical testing
# Loop over unique conditions
for (condition in unique(sequence$class)) {
  # Extract the matching runs for the current condition
  condition_runs <- sequence$sample[sequence$class == condition]
  # Extract columns from Data corresponding to the current condition
  condition_data <- data[, condition_runs, drop = FALSE]
  # Create a dataframe for the current condition
  dataframe_name <- paste("df_", gsub(" ", "_", condition), sep = "")
  assign(dataframe_name, condition_data)
}

ANOVA_data <- data # first we copy our data into a new object
ANOVA_data$Metabolite <- rownames(ANOVA_data) # then we add the proteins as a new column
ANOVA_data_long <- melt(ANOVA_data,
                        id.vars = "Metabolite",
                        variable.name = "sample",
                        value.name = "Expression")  # now we reshape the data table into a long format
ANOVA_data_long <- merge(ANOVA_data_long,
                         sequence[, c("sample", "class")],
                         by.x = "sample",
                         by.y = "sample",
                         all.x = TRUE) # finally we add the conditions as the last column

# to check distribution 
qq_plot <- ggplot(data = ANOVA_data_long, aes(sample = Expression)) +
  stat_qq(distribution = qnorm) +  # Use qnorm for standard normal distribution
  geom_abline(slope = 1, intercept = median(ANOVA_data_long$Expression), color = "red", linetype = "dashed") +  # Add line y = x
  labs(title = "Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_bw()
ggplotly(qq_plot)

# Create an empty dataframe to store Levene's test results
levene_results <- data.frame(Metabolite = character(),
                             p_value = numeric(),
                             stringsAsFactors = FALSE)

# Loop over unique metabolites
for (metabolite in unique(ANOVA_data_long$Metabolite)) {
  # Subset data for the current metabolite
  metabolite_data <- subset(ANOVA_data_long, Metabolite == metabolite)
  # Check if 'class' is a factor; if not, convert it to a factor
  if (!is.factor(metabolite_data$class)) {
    metabolite_data$class <- factor(metabolite_data$class)
  }
  # Perform Levene's test
  levene_result <- leveneTest(Expression ~ class, data = metabolite_data)
  # Extract p-value from Levene's test result
  p_value_levene <- levene_result$Pr[1]
  # Store Levene's test results in the dataframe
  levene_results <- rbind(levene_results,
                          data.frame(Metabolite = rep(metabolite, length(p_value_levene)),
                                     p_value = p_value_levene,
                                     stringsAsFactors = FALSE))
}
print(head(levene_results)) # the Levene's test tells you if you violate the assumption of equal variances among groups 

# Create an empty dataframe to store ANOVA test results
results_df_anova <- data.frame(Metabolite = character(),
                               p_value = numeric(),
                               stringsAsFactors = FALSE)
# Loop over unique proteins - SECOND FOR LOOP
for (metabolite in unique(ANOVA_data_long$Metabolite)) {
  # Subset data for the current protein
  metabolite_data <- subset(ANOVA_data_long, Metabolite == metabolite)
  # Fit the ANOVA (this anova assumes var is equal, if that is not the case from the Levene's test you set var.equal = F )
  model <-  oneway.test(Expression ~ class,
                        data = metabolite_data,
                        var.equal = TRUE) # set var.equal = FALSE here if you violate the homogeneity of variance
  # Extract p-value from the ANOVA model
  p_value <- model$p.value
  # Store results in the dataframe
  results_df_anova <- rbind(results_df_anova,
                            data.frame(Metabolite = rep(metabolite, length(p_value)),
                                       p_value = p_value))
}

non_significant_proteins <- results_df_anova$Metabolite[results_df_anova$p_value >= 0.05] # Let us create a value with the names of the non-significant proteins
print(head(results_df_anova)) # the ANOVA tells you if any of the group means differ but NOT which ones differ from one another for the given protein

# Create an empty dataframe to store Tukey's test results
results_df_tukey <- data.frame(
  Metabolite = character(),
  Comparison = character(),
  diff = numeric(),
  lwr = numeric(),
  upr = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop over unique proteins
for (metabolite in unique(ANOVA_data_long$Metabolite)) {
  # Subset data for the current protein
  metabolite_data <- subset(ANOVA_data_long, Metabolite == metabolite)
  # Check if 'class' is a factor; if not, convert it to a factor
  if (!is.factor(metabolite_data$class)) {
    metabolite_data$class <- factor(metabolite_data$class)
  }
  # Check if the protein is significant in the ANOVA
  is_significant <- results_df_anova$p_value[results_df_anova$Metabolite == metabolite] < 0.05
  if (is_significant) {
    # Perform Tukey's post hoc test
    tukey_result <- TukeyHSD(aov(Expression ~ class, data = metabolite_data))
    # Extract relevant information from the Tukey's test result
    comparisons <- as.data.frame(tukey_result$class[, c("diff", "lwr", "upr", "p adj")])
    # Add protein name and group information to the output
    comparisons$Metabolite <- rep(metabolite, nrow(comparisons))
    comparisons$Comparison <- rownames(comparisons)
    # Store results in the dataframe
    results_df_tukey <- rbind(results_df_tukey, comparisons)
  }
}

# Print the Tukey's test results
print(head(results_df_tukey))

# Now we did our statistical test correctly but we want to visualize the result.For that we create volcano plots for each comparison
# Split the dataframe by the 'Comparison' column
list_of_comparison_dataframes <- split(results_df_tukey, results_df_tukey$Comparison)
# Loop over the list of dataframes
for (i in seq_along(list_of_comparison_dataframes)) {
  # Extract the current dataframe
  current_df <- list_of_comparison_dataframes[[i]]
  # Extract comparison name
  comparison_name <- unique(current_df$Comparison)
  # Create a new dataframe with the same content
  new_df_name <- paste("df_", gsub("-", "_", comparison_name), sep = "")
  assign(new_df_name, current_df)
}

# we now use the Condition specific dataframes and first filter away the non-significant proteins
df_1 <- df_1[!rownames(df_1) %in% non_significant_proteins, ]
df_2 <- df_2[!rownames(df_2) %in% non_significant_proteins, ]
df_3 <- df_3[!rownames(df_3) %in% non_significant_proteins, ]
df_4 <- df_4[!rownames(df_4) %in% non_significant_proteins, ]

# AND THEN we calculate the Log2 Fold change for each significant protein for our comparisons
df_2_1$Log2FC<- rowMeans(df_2) - rowMeans(df_1)
rownames(df_2_1) <- df_2_1$Metabolite
df_3_1$Log2FC<- rowMeans(df_3) - rowMeans(df_1)
rownames(df_3_1) <- df_3_1$Metabolite
df_4_1$Log2FC<- rowMeans(df_4) - rowMeans(df_1)
rownames(df_4_1) <- df_4_1$Metabolite

df_3_2$Log2FC<- rowMeans(df_3) - rowMeans(df_2)
rownames(df_3_2) <- df_3_2$Metabolite
df_4_2$Log2FC<- rowMeans(df_4) - rowMeans(df_2)
rownames(df_4_2) <- df_4_2$Metabolite

df_4_3$Log2FC<- rowMeans(df_4) - rowMeans(df_3)
rownames(df_4_3) <- df_4_3$Metabolite


# Determine the top 10 upregulated and downregulated proteins
top_upregulated_2_1 <- head(df_2_1[order(-df_2_1$Log2FC), ], 10)
top_downregulated_2_1 <- head(df_2_1[order(df_2_1$Log2FC), ], 10)

# Create the volcano plot
volcano_plot_2_1 <- ggplot(df_2_1, aes(x = Log2FC, y = -log10(`p adj`))) +
  geom_point(aes(color = ifelse(`p adj` < 0.05 & Log2FC > 0.3785, "Upregulated", 
                                ifelse(`p adj` < 0.05 & Log2FC < -0.3785, "Downregulated", "Below Threshold"))), size = 1.5) +
  xlim(c(min(df_2_1$Log2FC) - 1,
         max(df_2_1$Log2FC) + 1)) +
  ylim(c(0, max(-log10(df_2_1$`p adj`)) + 1)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot Group2 vs Group1",
       x = "Log2 Fold Change",
       y = "-log10(p-value)") +
  theme_bw() +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Below Threshold" = "grey")) +
  guides(color = guide_legend(title = NULL)) +
  geom_text(data = rbind(top_upregulated_2_1, top_downregulated_2_1), 
            aes(label = rownames(rbind(top_upregulated_2_1, top_downregulated_2_1))), 
            size = 3, nudge_y = 0.5)

# Convert to plotly
ggplotly(volcano_plot_2_1)








