library(shiny)
library(ggplot2)
library(plotly)
library(dbscan)
library(factoextra)
library(caret)
library(ggdendro)
library(dplyr)

# Source the external functions file
source("TestFunctions.R")

# Define UI for application
ui <- fluidPage(
  titlePanel("Clustering and Visualization with Shiny App"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("sequenceFile", "Choose Sequence CSV File", accept = ".csv"),
      fileInput("dataFile", "Choose Data CSV File", accept = ".csv"),
      actionButton("loadData", "Load Data"),
      actionButton("runDataCleaning", "Run Data Cleaning"),
      hr(),
      uiOutput("tabSelection")
    ),
    
    mainPanel(
      tabsetPanel(id = "tabs",
                  tabPanel("Data",
                           tabsetPanel(
                             tabPanel("Sequence File",
                                      actionButton("load_data", "Load Data"),                                   # Button to load data
                                      DTOutput("sequence_table")),        # Tab for sequence file data
                             tabPanel("Data File", DTOutput("data_table")),                # Tab for data file
                             tabPanel("Cleaned Sequence File",
                                      actionButton("run_data_cleaning", "Run Data Cleaning"),                   # Button to run data cleaning
                                      DTOutput("cleaned_sequence_table")),  # Tab for cleaned sequence file
                             tabPanel("Cleaned Data File", DTOutput("cleaned_data_table"))           # Tab for cleaned data file
                           )),
                  tabPanel("Dimensional Reduction",
                           actionButton("run_pca", "Run PCA"),                # Button to run PCA
                           plotlyOutput("pca_plot"),                          # Output for PCA plot
                           plotlyOutput("scree_plot")                         # Output for Scree plot
                  ),
                  tabPanel("Outlier Detection",
                           tabsetPanel(
                             tabPanel("K-means",
                                      selectInput("kmeans_eval_method", "Choose Evaluation Method:",  # Dropdown to select evaluation method for K-means
                                                  choices = c("Within Sum of Square (WSS)" = "wss", "Silhouette" = "silhouette", "Gap Statistic" = "gap_stat")),
                                      actionButton("compute_kmeans_eval", "Compute Evaluation"),      # Button to compute K-means evaluation
                                      plotlyOutput("kmeans_eval_plot"),                                # Output for K-means evaluation plot
                                      numericInput("num_clusters", "Number of Clusters (k):", value = 3, min = 1),  # Input for number of clusters
                                      numericInput("percentile_threshold", "Percentile Threshold:", value = 95, min = 0, max = 100, step = 1), # Input for percentile threshold
                                      actionButton("run_kmeans", "Run K-means"),                      # Button to run K-means
                                      plotlyOutput("kmeans_plot"),                                    # Output for K-means plot
                                      DTOutput("kmeans_outliers"),                                    # Output for K-means outliers table
                                      actionButton("remove_kmeans_outliers", "Remove Selected Outliers"),  # Button to remove selected K-means outliers
                                      actionButton("save_cleaned_kmeans", "Save Cleaned Data")        # Button to save cleaned K-means data
                             ),
                             tabPanel("Hierarchical",
                                      selectInput("clustering_method", "Select Clustering Method:",  # Dropdown to select clustering method
                                                  choices = c("Single" = "single",
                                                              "Complete" = "complete",
                                                              "Average" = "average",
                                                              "Ward's D2" = "ward.D2")),
                                      numericInput("num_clusters_hierarchical", "Number of Clusters (k):", value = 3, min = 1),  # Input for number of clusters
                                      numericInput("threshold", "Dendrogram Threshold (Distance):", value = 5, min = 0),  # Input for dendrogram threshold
                                      actionButton("run_hierarchical", "Run Hierarchical Clustering"),  # Button to run hierarchical clustering
                                      plotlyOutput("hclust_plot"),                                      # Output for hierarchical clustering plot
                                      plotlyOutput("conf_matrix_plot"),                                 # Output for confusion matrix plot
                                      plotlyOutput("dendrogram_plot"),                                  # Output for dendrogram plot
                                      DTOutput("hierarchical_outliers"),                                # Output for hierarchical outliers table
                                      actionButton("remove_hierarchical_outliers", "Remove Selected Outliers"),  # Button to remove selected hierarchical outliers
                                      actionButton("save_cleaned_hierarchical", "Save Cleaned Data")    # Button to save cleaned hierarchical data
                             ),
                             tabPanel("DBSCAN",
                                      numericInput("knn", "Choose k for kNN Distance Plot:", value = 5, min = 1),  # Input for k in kNN
                                      actionButton("compute_knn", "Compute kNN Distance Plot"),                   # Button to compute kNN distance plot
                                      plotlyOutput("knn_plot"),                                                    # Output for kNN plot
                                      numericInput("eps", "Choose epsilon for DBSCAN:", value = 0.5, min = 0.01, step = 0.1),  # Input for epsilon in DBSCAN
                                      numericInput("min_pts_dbscan", "Choose minPts for DBSCAN:", value = 5, min = 1),  # Input for minPts in DBSCAN
                                      actionButton("run_dbscan", "Run DBSCAN"),                                   # Button to run DBSCAN
                                      plotlyOutput("dbscan_plot"),                                                # Output for DBSCAN plot
                                      DTOutput("dbscan_outliers"),                                                # Output for DBSCAN outliers table
                                      actionButton("remove_dbscan_outliers", "Remove Selected Outliers"),         # Button to remove selected DBSCAN outliers
                                      actionButton("save_cleaned_dbscan", "Save Cleaned Data")                   # Button to save cleaned DBSCAN data
                             ),
                             tabPanel("HDBSCAN",
                                      numericInput("min_pts_hdbscan", "Choose minPts for HDBSCAN:", value = 5, min = 1),  # Input for minPts in HDBSCAN
                                      numericInput("threshold_hdbscan", "Outlier Threshold for HDBSCAN:", value = 0.85, min = 0.01, max = 1),  # Input for outlier threshold in HDBSCAN
                                      actionButton("run_hdbscan", "Run HDBSCAN"),                                          # Button to run HDBSCAN
                                      plotlyOutput("hdbscan_plot"),                                                        # Output for HDBSCAN plot
                                      DTOutput("hdbscan_outliers"),                                                        # Output for HDBSCAN outliers table
                                      actionButton("remove_hdbscan_outliers", "Remove Selected Outliers"),                 # Button to remove selected HDBSCAN outliers
                                      actionButton("save_cleaned_hdbscan", "Save Cleaned Data")                            # Button to save cleaned HDBSCAN data
                             ),
                             tabPanel("OPTICS",
                                      numericInput("min_pts_optics", "Choose minPts for OPTICS:", value = 5, min = 1),  # Input for minPts in OPTICS
                                      numericInput("eps_optics", "Choose eps for OPTICS (optional):", value = NA, min = 0.1, step = 0.1),  # Input for eps in OPTICS
                                      numericInput("eps_cl_optics", "Choose cutoff (eps_cl) for OPTICS:", value = 0.5, min = 0.1, step = 0.1),  # Input for eps_cl in OPTICS
                                      actionButton("run_optics", "Run OPTICS"),                                     # Button to run OPTICS
                                      plotOutput("optics_reachability_plot"),                                      # Output for OPTICS reachability plot
                                      plotOutput("reachability_plot_threshold"),                                   # Output for reachability plot threshold
                                      plotlyOutput("cluster_plot"),                                                # Output for cluster plot
                                      DTOutput("optics_outliers"),                                                 # Output for OPTICS outliers table
                                      actionButton("remove_optics_outliers", "Remove Selected Outliers"),          # Button to remove selected OPTICS outliers
                                      actionButton("save_cleaned_optics", "Save Cleaned Data")                     # Button to save cleaned OPTICS data
                             ),
                             tabPanel("LOF",
                                      numericInput("lof_threshold", "Threshold for LOF:", value = 1.5, min = 0, step = 0.1),  # Input for threshold in LOF
                                      numericInput("lof_k", "k for LOF:", value = 4, min = 1),                                 # Input for k in LOF
                                      actionButton("run_lof", "Run LOF"),                                                    # Button to run LOF
                                      plotlyOutput("lof_plot"),                                                             # Output for LOF plot
                                      plotlyOutput("lof_od_plot"),                                                          # Output for LOF outlier detection plot
                                      DTOutput("lof_outliers"),                                                             # Output for LOF outliers table
                                      actionButton("remove_lof_outliers", "Remove Selected Outliers"),                      # Button to remove selected LOF outliers
                                      actionButton("save_cleaned_lof", "Save Cleaned Data")                                 # Button to save cleaned LOF data
                             )
                           )),
                  tabPanel("Visualization",
                           tabsetPanel(
                             tabPanel("Heatmap",
                                      numericInput("num_top_features", "Number of Top Features:", value = 20, min = 1),  # Input for number of top features in heatmap
                                      selectInput("dendrogram_option", "Dendrogram Option:",                          # Dropdown to select dendrogram option
                                                  choices = c("both", "row", "column", "none"), selected = "both"),
                                      actionButton("run_heatmap", "Generate Heatmap"),                                # Button to generate heatmap
                                      plotlyOutput("heatmap_plot"),                                                   # Output for heatmap plot
                                      plotlyOutput("selected_features_heatmap_plot")                                  # Output for selected features heatmap plot
                             ),
                             tabPanel("Volcano Plot",
                                      selectInput("group1", "Select Group 1", choices = NULL),                        # Dropdown to select group 1 for volcano plot
                                      selectInput("group2", "Select Group 2", choices = NULL),                        # Dropdown to select group 2 for volcano plot
                                      numericInput("log2fc_threshold", "Log2 Fold Change Threshold:", value = 2),     # Input for log2 fold change threshold
                                      numericInput("pval_threshold", "P-value Threshold:", value = 0.05),             # Input for p-value threshold
                                      actionButton("run_volcano_plot", "Generate Volcano Plot"),                      # Button to generate volcano plot
                                      plotlyOutput("volcano_plot"),                                                   # Output for volcano plot
                                      DTOutput("volcano_table")                                                       # Output for volcano plot table
                             )
                           )
                  )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  data <- reactiveVal()          # Holds the raw data file
  sequence <- reactiveVal()      # Holds the sequence file
  cleaned_data <- reactiveVal()  # Holds the cleaned data
  cleaned_sequence <- reactiveVal() # Holds the cleaned sequence
  pca_result <- reactiveVal()    # Holds the PCA result
  
  # Load Data from CSV Files
  observeEvent(input$load_data, {
    req(input$sequence_file, input$data_file)
    
    sequence_file <- input$sequence_file$datapath
    data_file <- input$data_file$datapath
    
    seq_data <- read.csv(sequence_file)
    raw_data <- read.csv(data_file)
    
    sequence(seq_data)
    data(raw_data)
    
    print("Data loaded successfully.")
    print(paste("Sequence dimensions:", dim(seq_data)))
    print(paste("Data dimensions:", dim(raw_data)))
  })
  
  # Run Data Cleaning
  observeEvent(input$run_data_cleaning, {
    req(data(), sequence())
    
    # Filter data and sequence based on labels
    filtered_data <- data()[, sequence()[, 'labels'] %in% c("Name", "Sample")]
    filtered_sequence <- sequence()[sequence()[, 'labels'] %in% c("Name", "Sample"), ]
    
    number_of_rows_before <- nrow(filtered_data)
    
    # Clean compound names using custom functions
    filtered_data[, 1] <- sapply(filtered_data[, 1], extract_pattern)
    filtered_data[, 1] <- sapply(filtered_data[, 1], format_strings)
    
    # Make unique row names and perform additional data transformations
    unique_row_names <- make.unique(as.character(filtered_data[, "Compound.Name"]))
    rownames(filtered_data) <- unique_row_names
    filtered_data <- filtered_data[, -1]
    
    # Perform logarithmic transformation and clean data
    filtered_data <- log2(filtered_data)
    filtered_data <- na.omit(filtered_data)
    data_scaled <- as.data.frame(scale(filtered_data))
    data_scaled_T <- as.data.frame(t(data_scaled))
    filtered_sequence <- filtered_sequence[-1, ]
    
    # Debugging information
    print(paste("Number of rows before cleaning:", number_of_rows_before))
    print(paste("Number of rows after cleaning:", nrow(filtered_data)))
    print(paste("Filtered data dimensions (scaled and transposed):", dim(data_scaled_T)))
    print(paste("Filtered sequence dimensions:", dim(filtered_sequence)))
    
    # Check for empty data frames
    if (nrow(data_scaled_T) == 0 || ncol(data_scaled_T) == 0) {
      showNotification("Error: Cleaned data has zero dimensions after scaling and transposing.", type = "error")
      return(NULL)
    }
    
    # Update the reactive values
    cleaned_data(data_scaled_T)
    cleaned_sequence(filtered_sequence)
  })
  
  # Display Raw Data and Cleaned Data
  output$sequence_table <- renderTable({
    req(sequence())
    head(sequence(), 20)
  })
  
  output$data_table <- renderTable({
    req(data())
    head(data(), 20)
  })
  
  output$cleaned_sequence_table <- renderTable({
    req(cleaned_sequence())
    head(cleaned_sequence(), 20)
  })
  
  output$cleaned_data_table <- renderTable({
    req(cleaned_data())
    head(cleaned_data(), 20)
  }, rownames = TRUE)
  
  # Perform PCA analysis
  observeEvent(input$run_pca, {
    req(cleaned_data(), cleaned_sequence())
    
    cleaned_data_df <- cleaned_data()
    cleaned_sequence_df <- cleaned_sequence()
    
    print(paste("Data dimensions before PCA:", dim(cleaned_data_df)))
    print(paste("Sequence dimensions before PCA:", dim(cleaned_sequence_df)))
    
    # Check if the data has non-zero dimensions
    if (nrow(cleaned_data_df) == 0 || ncol(cleaned_data_df) == 0) {
      showNotification("Error: Cleaned data has zero dimensions and cannot be used for PCA.", type = "error")
      return(NULL)
    }
    
    # Ensure no missing values
    if (anyNA(cleaned_data_df)) {
      showNotification("Error: Cleaned data contains NA values and cannot be used for PCA.", type = "error")
      return(NULL)
    }
    
    pca_result(perform_pca_analysis(cleaned_data_df, cleaned_sequence_df))
  })
  
  # PCA Plot
  output$pca_plot <- renderPlotly({
    req(pca_result())
    PCA.df <- pca_result()$PCA.df
    create_pca_plot(PCA.df, custom_colors)
  })
  
  # Scree Plot
  output$scree_plot <- renderPlotly({
    req(pca_result())
    pca_result()$Scree_plotly
  })
  
  # Compute and display WSS plot
  observeEvent(input$compute_wss, {
    req(pca_result())
    PCA.df <- pca_result()$PCA.df
    wss <- fviz_nbclust(PCA.df[,2:3], FUNcluster = stats::kmeans, method = "wss", k.max = nrow(PCA.df)-1)
    output$wss_plot <- renderPlotly(ggplotly(wss))
  })
  
  # Compute and display Silhouette plot
  observeEvent(input$compute_silhouette, {
    req(pca_result())
    PCA.df <- pca_result()$PCA.df
    sil <- fviz_nbclust(PCA.df[,2:3], FUNcluster = stats::kmeans, method = "silhouette", k.max = nrow(PCA.df)-1)
    output$sil_plot <- renderPlotly(ggplotly(sil))
  })
  
  # Compute and display Gap Statistics plot
  observeEvent(input$compute_gap, {
    req(pca_result())
    PCA.df <- pca_result()$PCA.df
    gap <- fviz_nbclust(PCA.df[,2:3], FUNcluster = stats::kmeans, method = "gap_stat", k.max = nrow(PCA.df)-1)
    output$gap_plot <- renderPlotly(ggplotly(gap))
  })
  
  # Run K-means Clustering
  observeEvent(input$run_kmeans, {
    req(pca_result())
    PCA.df <- pca_result()$PCA.df
    k <- input$num_clusters
    
    kmeans_res <- stats::kmeans(PCA.df[,2:3], centers = k)
    cluster_labels <- kmeans_res$cluster
    
    kmeans.df <- data.frame(
      x = PCA.df[,2],
      y = PCA.df[,3],
      cluster = factor(cluster_labels),
      Sample = PCA.df$Sample
    )
    
    plot_kmeans <- ggplot(kmeans.df, aes(x = x, y = y, color = cluster, text = Sample)) +
      ggtitle("K-means Clustering") +
      geom_point(size = 2) +
      labs(color = "Clusters") +
      xlab("PC1") +
      ylab("PC2") +
      scale_color_manual(values = custom_colors) +
      theme_bw()
    
    output$kmeans_plot <- renderPlotly(ggplotly(plot_kmeans))
  })
  
  # Hierarchical Clustering Plot
  observeEvent(input$run_hierarchical, {
    req(pca_result())
    PCA.df <- pca_result()$PCA.df
    k <- input$num_clusters_hierarchical
    method <- input$clustering_method
    
    output$hclust_plot <- renderPlotly({
      perform_hierarchical_clustering(PCA.df, cleaned_sequence(), custom_colors, method, k)
    })
    output$conf_matrix_plot <- renderPlotly({
      perform_confusion_matrix(PCA.df, cleaned_sequence(), custom_colors, method, k)
    })
    output$dendrogram_plot <- renderPlotly({
      perform_dendrogram(PCA.df, cleaned_sequence(), custom_colors, method, k)
    })
  })
  
  # Compute and display kNN distance plot
  observeEvent(input$compute_knn, {
    req(pca_result())
    PCA.df <- pca_result()$PCA.df
    k <- input$knn
    
    # Debug print
    print(paste("Computing kNN distance plot with k =", k))
    
    output$knn_plot <- renderPlotly({
      perform_kNN_dist_plot(PCA.df, k)
    })
  })
  
  # Run DBSCAN
  observeEvent(input$run_dbscan, {
    req(pca_result())
    PCA.df <- pca_result()$PCA.df
    eps <- input$eps
    min_pts <- input$min_pts_dbscan
    
    # Debug print
    print(paste("Running DBSCAN with eps =", eps, "and minPts =", min_pts))
    
    output$dbscan_plot <- renderPlotly({
      perform_dbscan_clustering(PCA.df, eps, min_pts)
    })
  })
  
  # Run HDBSCAN
  observeEvent(input$run_hdbscan, {
    req(pca_result())
    PCA.df <- pca_result()$PCA.df
    min_pts <- input$min_pts_hdbscan
    
    # Debug print
    print(paste("Running HDBSCAN with minPts =", min_pts))
    
    output$hdbscan_plot <- renderPlotly({
      perform_hdbscan_clustering(PCA.df, min_pts)
    })
  })
  
  # Run OPTICS
  observeEvent(input$run_optics, {
    req(pca_result())
    PCA.df <- pca_result()$PCA.df
    min_pts <- input$min_pts_optics
    eps <- if (is.na(input$eps_optics)) NULL else input$eps_optics
    eps_cl <- input$eps_cl_optics
    
    # Debug print
    showNotification(paste("Running OPTICS with minPts =", min_pts, ", eps =", eps, ", and eps_cl =", eps_cl))
    
    tryCatch({
      optics_results <- perform_optics_analysis(PCA.df, k = input$knn, eps, min_pts, eps_cl)
      
      output$optics_plot <- renderPlotly({
        optics_results$kNN_plot
      })
      
      output$optics_reachability_plot <- renderPlotly({
        optics_results$reachability_plot
      })
      
      output$reachability_plot_threshold <- renderPlotly({
        optics_results$reachability_plot_threshold
      })
      
      output$cluster_plot <- renderPlotly({
        optics_results$cluster_plot
      })
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # LOF Plot
  observeEvent(input$run_lof, {
    req(pca_result())
    lof_result <- calculate_and_plot_lof(pca_result()$PCA.df, threshold = 1.5)
    output$lof_plot <- renderPlotly(lof_result$lof_plotly)
    output$lof_od_plot <- renderPlotly(lof_result$lof_od_plotly)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)