# Load required libraries
library(shiny)
library(shinydashboard)
library(bslib)
library(DESeq2)
library(tximport)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(matrixStats)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(dplyr)
library(rlang)
library(rhdf5)

# Load utility functions from an external file
source("utils.R")

# Define the User Interface (UI)
ui <- fluidPage(
  # Set the theme for the app using bslib
  theme = bs_theme(bootswatch = "darkly",
                   secondary = "#BA0C2F",
                   "table-bg" = "primary"),
  titlePanel("DESeq2 Analysis and Visualizations"),
  
  # Define a sidebar layout with input controls and a main panel for plots
  sidebarLayout(
    # Sidebar Panel: Contains data input, analysis parameters, and multiple tabs
    sidebarPanel(
      # Display a logo image
      HTML('<img src="logo.png" width="100%" height="auto">'),
      br(), br(),
      tabsetPanel(
        # Tab for Data Input
        tabPanel("Data Input",
                 # Select the type of count data to use
                 radioButtons("data_type", "Select Data Type:",
                              choices = c("Kallisto (tximport)" = "kallisto",
                                          "CSV Count Matrix" = "csv")),
                 # Conditional panel: file input for Kallisto files
                 conditionalPanel(
                   condition = "input.data_type == 'kallisto'",
                   fileInput("kallisto_files", "Upload Kallisto Abundance Files (.h5)", 
                             multiple = TRUE, accept = ".h5")
                 ),
                 # Conditional panel: file input for CSV count matrix
                 conditionalPanel(
                   condition = "input.data_type == 'csv'",
                   fileInput("csv_counts", "Upload CSV Count Matrix", accept = ".csv")
                 ),
                 # Button and input for metadata
                 downloadButton("download_metadata_template", "Download Metadata Template"),
                 fileInput("metadata_file", "Upload Metadata CSV", accept = ".csv"),
                 br(),
                 # Design variables for the DESeq2 model
                 selectizeInput("design_vars", "Select design variables:",
                                choices = NULL, multiple = TRUE),
                 # Show a text preview of the design formula
                 textOutput("design_preview"),
                 br(),
                 # Button to run DESeq2 analysis
                 actionButton("run_deseq", "Run DESeq2 Analysis")
        ),
        
        # Tab for Pairwise Comparison
        tabPanel("Pairwise Comparison",
                 # Inputs for selecting the contrast variable and groups
                 selectInput("contrast_var", "Select variable for comparison", choices = NULL),
                 selectInput("group1", "Test Group", choices = NULL),
                 selectInput("group2", "Control Group", choices = NULL),
                 # Options for labeling genes in the volcano plot
                 radioButtons("label_option", "Label genes in Volcano Plot", 
                              choices = c("Top N genes" = "top", "Specific gene" = "specific"),
                              selected = "top", inline = TRUE),
                 # Conditional input for top N genes
                 conditionalPanel(
                   condition = "input.label_option == 'top'",
                   numericInput("top_n", "Number of top genes to label", value = 10, min = 1)
                 ),
                 # Conditional input for a specific gene label
                 conditionalPanel(
                   condition = "input.label_option == 'specific'",
                   selectizeInput("specific_gene", "Gene to label", choices = NULL, options = list(placeholder = "Select gene"))
                 ),
                 # Cutoff inputs for adjusted p-value and log2 fold change
                 numericInput("padj_cutoff", "Adjusted P-value cutoff", value = 0.05, 
                              min = 0, max = 1, step = 0.01),
                 numericInput("log2fc_cutoff", "Log2 Fold Change cutoff", value = 1,
                              min = 0, step = 0.1),
                 # Button to run the pairwise contrast
                 actionButton("run_contrast", "Run Pairwise Comparison")
        ),
        
        # Tab for Gene-level Analysis
        tabPanel("Gene-level Analysis",
                 # Input for selecting a gene for detailed expression analysis
                 selectizeInput("gene_select", "Select Gene", choices = NULL, options = list(placeholder = "Search gene...")),
                 # Button to plot gene expression
                 actionButton("plot_gene", "Plot Gene Expression")
        )
      ),
      br(),
      # Footer information
      tags$p("Created by Andy Ring"),
      tags$p("Version 1.0.2 | May 16th, 2025")
    ),
    
    # Main Panel: Contains the output plots and download buttons
    mainPanel(
      tabsetPanel(
        # Tab for Exploratory Analysis (PCA and Heatmap)
        tabPanel("Exploratory Analysis", 
                 # PCA plot card
                 card(
                   card_header("PCA"),
                   checkboxInput("show_pca_labels", "Show PCA Sample Labels", value = FALSE),
                   plotOutput("pca_plot"),
                   fluidRow(
                     column(3, downloadButton("download_pca", "Download PCA Plot"))
                   )
                 ),
                 # Heatmap card with an added option to toggle column clustering
                 card(
                   card_header("Heatmap"),
                   checkboxInput("show_gene_names", "Show Gene Names", value = FALSE),
                   checkboxInput("cluster_cols", "Cluster Columns", value = TRUE),  # New checkbox for clustering columns
                   plotOutput("heatmap_plot"),
                   fluidRow(
                     column(3, downloadButton("download_heatmap", "Download Heatmap"))
                   )
                 )
        ),
        
        # Tab for Pairwise Comparison (Volcano Plot)
        tabPanel("Pairwise Comparison",
                 card(
                   card_header("Volcano Plot"),
                   plotOutput("volcano_plot"),
                   fluidRow(
                     column(3, downloadButton("download_volcano", "Download Volcano Plot")),
                     column(3, downloadButton("download_results", "Download CSV Results"))
                   )
                 ),
                 br(), br(),
                 fluidRow(
                   column(6, value_box(title = "Significant Genes", textOutput("significant_genes_count_box"), theme = "success"))
                 )
        ),
        
        # Tab for Gene-level Analysis (Gene Expression)
        tabPanel("Gene-level Analysis",
                 card(
                   card_header("Gene Expression"),
                   plotOutput("gene_barplot"),
                   fluidRow(
                     column(6, downloadButton("download_gene", "Download Gene Expression Plot"))
                   )
                 ),
                 br(), br(),
                 fluidRow(
                   column(6, value_box(title = "Adjusted P-value", textOutput("gene_adjP_box"), theme = "success"))
                 )
        )
      )
    )
  )
)

# Define the Server logic
server <- function(input, output, session) {
  
  # Reactive expression to load and process count data based on selected data type
  data_counts <- reactive({
    req(input$data_type)
    if (input$data_type == "kallisto") {
      req(input$kallisto_files)
      files <- input$kallisto_files$datapath
      names(files) <- sapply(input$kallisto_files$name, function(x) tools::file_path_sans_ext(x))
      txi <- tximport(files, type = "kallisto", txOut = TRUE, countsFromAbundance = "scaledTPM")
      countData <- as.data.frame(txi$counts)
      countData <- convertIDs(countData, session)
      return(countData)
    } else {
      req(input$csv_counts)
      countData <- read.csv(input$csv_counts$datapath, row.names = 1, check.names = FALSE)
      countData <- convertIDs(countData, session)
      return(countData)
    }
  })
  
  # Reactive expression to load and process metadata
  metadata <- reactive({
    req(input$metadata_file)
    meta <- read.csv(input$metadata_file$datapath, row.names = 1, check.names = FALSE)
    meta[] <- lapply(meta, function(x) if (is.character(x)) as.factor(x) else x)
    return(meta)
  })
  
  # Render a text preview of the design formula for DESeq2 analysis
  output$design_preview <- renderText({
    req(input$design_vars)
    paste("~", paste(input$design_vars, collapse = " + "))
  })
  
  # Update design variable choices based on metadata columns
  observeEvent(metadata(), {
    meta <- metadata()
    updateSelectizeInput(session, "design_vars", choices = names(meta))
  })
  
  # Create the DESeq2 dataset and run the differential expression analysis
  dds <- eventReactive(input$run_deseq, {
    req(data_counts(), metadata(), input$design_vars)
    countData <- data_counts()
    colData <- metadata()
    if (input$data_type == "kallisto") { 
      countData <- round(countData) 
    }
    common <- intersect(rownames(colData), colnames(countData))
    if (length(common) == 0) {
      stop("No matching samples between count data and metadata. Please check your sample names.")
    }
    countData <- countData[, common, drop = FALSE]
    colData <- colData[common, , drop = FALSE]
    design_formula <- as.formula(paste("~", paste(input$design_vars, collapse = " + ")))
    dds_obj <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = design_formula)
    dds_obj <- dds_obj[rowSums(counts(dds_obj)) > 1, ]
    withProgress(message = "Running DESeq2", value = 0, {
      incProgress(0.5)
      dds_obj <- DESeq(dds_obj)
      incProgress(0.5)
    })
    return(dds_obj)
  })
  
  # Update gene selection inputs once DESeq2 analysis is complete
  observeEvent(dds(), {
    updateSelectizeInput(session, "specific_gene", choices = rownames(dds()))
    updateSelectizeInput(session, "gene_select", choices = rownames(dds()))
  })
  
  # Update contrast variable choices based on the design variables that are factors or characters
  observeEvent(dds(), {
    req(input$design_vars)
    colData_df <- as.data.frame(colData(dds()))
    factor_vars <- intersect(input$design_vars, names(colData_df)[sapply(colData_df, function(x) is.factor(x) || is.character(x))])
    updateSelectInput(session, "contrast_var", choices = factor_vars)
  })
  
  # Update the group selection inputs based on the selected contrast variable
  observeEvent(input$contrast_var, {
    req(dds(), input$contrast_var)
    colData_df <- as.data.frame(colData(dds()))
    groups <- unique(as.character(colData_df[[input$contrast_var]]))
    updateSelectInput(session, "group1", choices = groups)
    updateSelectInput(session, "group2", choices = groups, selected = if(length(groups) > 1) groups[2] else groups[1])
  })
  
  # Reactive expression to perform a variance stabilizing transformation (vst) for visualization
  vst_data <- reactive({
    req(dds())
    tryCatch({
      vst(dds(), blind = TRUE)
    }, error = function(e) {
      showNotification("vst transformation failed, switching to rlog transformation", type = "warning", duration = 5)
      rlog(dds(), blind = TRUE)
    })
  })
  
  # Render the PCA plot using transformed data
  output$pca_plot <- renderPlot({
    req(vst_data(), input$design_vars)
    pcaData <- plotPCA(vst_data(), intgroup = input$design_vars, returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    colorVar <- input$design_vars[1]
    p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = !!sym(colorVar))) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      theme_light()
    if (input$show_pca_labels) {
      p <- p + geom_text(aes(label = rownames(pcaData)), vjust = -1, size = 3)
    }
    p
  })
  
  # Render the heatmap with the option to toggle column clustering
  output$heatmap_plot <- renderPlot({
    req(vst_data())
    mat <- assay(vst_data())
    mat <- mat - rowMeans(mat)
    rv <- rowVars(mat)
    select <- order(rv, decreasing = TRUE)[1:min(50, length(rv))]
    annotation <- as.data.frame(colData(dds()))
    if (!is.null(input$design_vars) && length(input$design_vars) > 0) {
      annotation <- annotation[, input$design_vars, drop = FALSE]
    }
    pheatmap(mat[select, ],
             annotation_col = annotation,
             show_rownames = input$show_gene_names,
             cluster_cols = input$cluster_cols)  # Use the clustering option from UI
  })
  
  # Download handler for saving the PCA plot as a high-resolution image
  output$download_pca <- downloadHandler(
    filename = function() { "PCA_highres.png" },
    content = function(file) {
      req(vst_data(), input$design_vars)
      pcaData <- plotPCA(vst_data(), intgroup = input$design_vars, returnData = TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      colorVar <- input$design_vars[1]
      p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = !!sym(colorVar))) +
        geom_point(size = 3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        theme_light()
      if (input$show_pca_labels) {
        p <- p + geom_text(aes(label = rownames(pcaData)), vjust = -1, size = 3)
      }
      ggsave(filename = file, plot = p, width = 18, height = 12, dpi = 800)
    }
  )
  
  # Download handler for saving the heatmap as a high-resolution image,
  # ensuring the clustering option is applied
  output$download_heatmap <- downloadHandler(
    filename = function() { "heatmap_highres.png" },
    content = function(file) {
      req(vst_data())
      mat <- assay(vst_data())
      mat <- mat - rowMeans(mat)
      rv <- rowVars(mat)
      select <- order(rv, decreasing = TRUE)[1:min(50, length(rv))]
      annotation <- as.data.frame(colData(dds()))
      if (!is.null(input$design_vars) && length(input$design_vars) > 0) {
        annotation <- annotation[, input$design_vars, drop = FALSE]
      }
      p <- pheatmap(mat[select, ],
                    annotation_col = annotation,
                    show_rownames = input$show_gene_names,
                    cluster_cols = input$cluster_cols,  # Use clustering option from UI
                    silent = TRUE)
      ggsave(filename = file, plot = p$gtable, width = 12, height = 18, dpi = 800)
    }
  )
  
  # Download handler for the volcano plot image
  output$download_volcano <- downloadHandler(
    filename = function() { "volcano_plot_highres.png" },
    content = function(file) {
      req(contrast_res())
      results <- contrast_res()$res
      x_min <- min(results$log2FoldChange, na.rm = TRUE) - 1
      x_max <- max(results$log2FoldChange, na.rm = TRUE) + 1
      y_max <- max(-log10(results$padj), na.rm = TRUE) + 1
      gene_labs <- NULL
      sig_genes <- contrast_res()$sig
      if (nrow(sig_genes) > 0) {
        sig_genes <- sig_genes[order(sig_genes$padj), ]
        if (input$label_option == "top") {
          gene_labs <- rownames(sig_genes)[1:min(input$top_n, nrow(sig_genes))]
        } else if (input$label_option == "specific" && nzchar(input$specific_gene)) {
          gene_labs <- input$specific_gene
        }
      }
      p <- EnhancedVolcano(results,
                           x = "log2FoldChange",
                           y = "padj",
                           lab = rownames(results),
                           selectLab = gene_labs,
                           FCcutoff = input$log2fc_cutoff,
                           pCutoff = input$padj_cutoff,
                           legendPosition = "none",
                           title = "",
                           drawConnectors = TRUE,
                           xlim = c(x_min, x_max),
                           ylim = c(0, y_max),
                           labSize = 5.0,
                           pointSize = 5.0,
                           boxedLabels = TRUE)
      ggsave(filename = file, plot = p, width = 18, height = 12, dpi = 800)
    }
  )
  
  # Download handler for the gene expression plot
  output$download_gene <- downloadHandler(
    filename = function() { "gene_expression_plot.png" },
    content = function(file) {
      req(gene_expr())
      df <- gene_expr()$df
      df_summary <- df %>% group_by(group) %>% 
        summarise(mean_count = mean(count), se = sd(count)/sqrt(n()))
      p <- ggplot() +
        geom_jitter(data = df, aes(x = group, y = count), width = 0.1, color = "black", alpha = 0.5) +
        geom_col(data = df_summary, aes(x = group, y = mean_count, fill = group), 
                 position = position_dodge(width = 0.9), alpha = 0.7) +
        geom_errorbar(data = df_summary, aes(x = group, ymin = mean_count - se, ymax = mean_count + se), 
                      width = 0.2, position = position_dodge(width = 0.9)) +
        theme_light() +
        labs(title = paste("Gene Expression for", input$gene_select),
             y = "Normalized Count", x = "Group")
      ggsave(filename = file, plot = p, device = "png", width = 12, height = 18, dpi = 800)
    }
  )
  
  # Download handler for the CSV results from the pairwise contrast analysis
  output$download_results <- downloadHandler(
    filename = function() { "significant_genes.csv" },
    content = function(file) {
      req(contrast_res())
      write.csv(contrast_res()$sig, file)
    }
  )
  
  # Reactive expression to compute differential expression results for a given contrast
  contrast_res <- eventReactive(input$run_contrast, {
    req(dds(), input$contrast_var, input$group1, input$group2)
    res <- results(dds(), contrast = c(input$contrast_var, input$group1, input$group2))
    sig <- res[which(res$padj < input$padj_cutoff & abs(res$log2FoldChange) > input$log2fc_cutoff), ]
    list(res = as.data.frame(res), sig = as.data.frame(sig))
  })
  
  # Render the number of significant genes as text output
  output$significant_genes_count_box <- renderText({
    req(contrast_res())
    num_sig <- nrow(contrast_res()$sig)
  })
  
  # Render the volcano plot for the pairwise contrast
  output$volcano_plot <- renderPlot({
    req(contrast_res())
    results <- contrast_res()$res
    x_min <- min(results$log2FoldChange, na.rm = TRUE) - 1
    x_max <- max(results$log2FoldChange, na.rm = TRUE) + 1
    y_max <- max(-log10(results$padj), na.rm = TRUE) + 1
    gene_labs <- NULL
    sig_genes <- contrast_res()$sig
    if (nrow(sig_genes) > 0) {
      sig_genes <- sig_genes[order(sig_genes$padj), ]
      if (input$label_option == "top") {
        gene_labs <- rownames(sig_genes)[1:min(input$top_n, nrow(sig_genes))]
      } else if (input$label_option == "specific" && nzchar(input$specific_gene)) {
        gene_labs <- input$specific_gene
      }
    }
    p <- EnhancedVolcano(results,
                         x = "log2FoldChange",
                         y = "padj",
                         lab = rownames(results),
                         selectLab = gene_labs,
                         FCcutoff = input$log2fc_cutoff,
                         pCutoff = input$padj_cutoff,
                         legendPosition = "none",
                         title = "",
                         drawConnectors = TRUE,
                         xlim = c(x_min, x_max),
                         ylim = c(0, y_max),
                         labSize = 5.0,
                         pointSize = 5.0,
                         boxedLabels = TRUE)
    p
  })
  
  # Update gene selection input for gene-level analysis when DESeq2 results are available
  observeEvent(dds(), {
    updateSelectizeInput(session, "gene_select", choices = rownames(dds()))
  })
  
  # Reactive expression to get gene expression data for the selected gene and group information
  gene_expr <- eventReactive(input$plot_gene, {
    req(dds(), input$gene_select)
    gene <- input$gene_select
    norm_counts <- counts(dds(), normalized = TRUE)[gene, ]
    df <- data.frame(
      sample = names(norm_counts),
      count = as.numeric(norm_counts),
      stringsAsFactors = FALSE
    )
    if (!is.null(input$contrast_var)) {
      df$group <- as.factor(colData(dds())[[input$contrast_var]])
    } else {
      df$group <- "All"
    }
    adjP <- NA
    if (!is.null(contrast_res())) {
      res <- contrast_res()$res
      if (gene %in% rownames(res)) {
        adjP <- res[gene, "padj"]
      }
    }
    list(df = df, adjP = adjP)
  })
  
  # Render the gene expression barplot
  output$gene_barplot <- renderPlot({
    req(gene_expr())
    df <- gene_expr()$df
    df_summary <- df %>% group_by(group) %>% 
      summarise(mean_count = mean(count), se = sd(count)/sqrt(n()))
    ggplot() +
      geom_jitter(data = df, aes(x = group, y = count), width = 0.1, color = "black", alpha = 0.5) +
      geom_col(data = df_summary, aes(x = group, y = mean_count, fill = group), 
               position = position_dodge(width = 0.9), alpha = 0.7) +
      geom_errorbar(data = df_summary, aes(x = group, ymin = mean_count - se, ymax = mean_count + se),
                    width = 0.2, position = position_dodge(width = 0.9)) +
      theme_light() +
      labs(title = paste("Gene Expression for", input$gene_select),
           y = "Normalized Count", x = "Group")
  })
  
  # Render the adjusted p-value for the selected gene
  output$gene_adjP_box <- renderText({
    req(gene_expr())
    adjP <- gene_expr()$adjP
  })
  
  # Download handler for the metadata template CSV file
  output$download_metadata_template <- downloadHandler(
    filename = function() { "metadata_template.csv" },
    content = function(file) {
      template <- data.frame(
        sample = c("sample1", "sample2", "sample3"),
        Treatment = c("Control", "Treatment", "Control"),
        Condition = c("A", "A", "B"),
        stringsAsFactors = FALSE
      )
      write.csv(template, file, row.names = FALSE)
    }
  )
}

# Run the Shiny App
shinyApp(ui, server)



















