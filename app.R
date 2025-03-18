library(dplyr)
library(ggplot2)
library(tidyverse)

# Set the app directory
if (requireNamespace("rprojroot", quietly = TRUE)) {
  app_dir <- rprojroot::find_rstudio_root_file()
  print(app_dir)
}

# Define paths relative to the app directory
data_dir <- file.path(app_dir, "gtex_v10_shiny/data")  # Go one level up to access gtex_v10_shiny
raw_data_dir <- file.path(data_dir, "raw_data")

# Load tissue names
tissue_file <- file.path(data_dir, "tissue_names.txt")
if (file.exists(tissue_file)) {
  tissues <- readLines(tissue_file)
  tissues <- gsub("_", " ", tissues)
} else {
  stop("Tissue names file not found!")
}

# Load gene names
gene_file <- file.path(data_dir, "gene_names.txt")
if (file.exists(gene_file)) {
  genes <- readLines(gene_file)
} else {
  stop("Gene names file not found!")
}

# Define UI
ui <- fluidPage(
  titlePanel("GTEx Gene Expression Analysis"),
  sidebarLayout(
    sidebarPanel(
      selectInput("gene", "Select gene:", choices = genes),
      selectInput("tissue", "Select tissue:", choices = tissues),
      actionButton("plot", "Generate Plot"),
      checkboxInput("plot_type", "Show Violin & Box Plot (Uncheck for Regression Plot)",value = TRUE)  # New checkbox for plot type
      
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("tpmPlot")),
        tabPanel("P-Value Table", 
                 DT::dataTableOutput("pvalueTable"),
                 downloadButton("downloadTable", "Download P-Value Table"))
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  # Function to read and preprocess data
  read_and_preprocess_data <- function(gene, tissue) {
    req(gene, tissue)  # Ensure both inputs are available
    
    # # Define the URL for the gene expression data
    # exp_url <- sprintf("https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_%s.gct.gz", gsub(" ", "_", tissue))
    # 
    # # Check if the URL is accessible
    # if (http_error(exp_url)) {
    #   stop("Error: Unable to access data at ", exp_url)
    # }
    # 
    # # Open the connection to the gzipped file from the URL
    # con <- gzcon(url(exp_url, "rb"))
    # 
    # # Read the GCT file
    # exp <- read.table(con, sep = "\t", skip = 2, header = TRUE)
    # close(con)  # Close the connection
    exp.path <- file.path(raw_data_dir, sprintf("gene_tpm_v10_%s.gct.gz", gsub(" ", "_", tissue)))
    if (!file.exists(exp.path)) return(paste("Error: Expression data not found at", exp.path))
    exp <- read.table(gzfile(exp.path), sep = "\t", skip = 2, header = TRUE)
    
    # Load metadata file
    metadata_path <- file.path(raw_data_dir, "Updated_GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")
    if (!file.exists(metadata_path)) {
      stop("Error: Metadata file missing at ", metadata_path)
    }
    metadata <- read.table(metadata_path, sep = "\t", header = TRUE)
    
    # Extract data for the selected gene and tissue
    X <- exp %>% filter(Description == gene)
    if (nrow(X) == 0) return(NULL)
    
    # Prepare data
    X <- X %>%
      select(-Name, -Description) %>%
      pivot_longer(cols = everything(), names_to = "sample", values_to = "TPM") %>%
      mutate(donor = sub("^([^.]+)\\.([^.]+).*", "\\1-\\2", sample)) %>% 
      mutate(logTPM = log(TPM + 1)) %>% 
      select(donor, logTPM)
    
    # Merge with metadata
    mergedData <- left_join(X, metadata, by = "donor")
    if (nrow(mergedData) == 0) return(NULL)
    
    return(mergedData)
  }
  
  # Reactive expression to generate plot data
  plotData <- eventReactive(input$plot, {
    tryCatch({
      df <- read_and_preprocess_data(input$gene, input$tissue)
      validate(need(!is.null(df), "Error: No matching data for this gene-tissue pair."))
      df
    }, error = function(e) {
      return(data.frame(error_msg = e$message))
    })
  })
  
  # Render the plot
  output$tpmPlot <- renderPlot({
    df <- plotData()
    
    # Check for errors in the reactive data
    validate(need(!"error_msg" %in% colnames(df) || is.null(df$error_msg), df$error_msg))
    
    # Plot using ggplot2
    if (input$plot_type) {
      # Violin plot
      ggplot(df, aes(x = age_plot, y = logTPM, fill = sex_plot)) +
        geom_violin(trim = FALSE, alpha = 0.3) +  # Violin plot with transparency
        geom_boxplot(width = 0.2, alpha = 0.6) +  # Box plot on top with transparency
        scale_fill_manual(name = "Sex", values = c("Male" = "steelblue", "Female" = "red")) +
        ggtitle(sprintf("Expression of %s in %s (Box & Violin Plot)", input$gene, input$tissue)) +
        xlab("Age") + ylab("logTPM") + 
        theme_minimal()
    } else {
      ggplot(df, aes(x = age_plot, y = logTPM, color = sex_plot)) +
        geom_smooth(method = "lm", formula = y ~ x, fill = "lightgray", alpha = 0.3) +
        geom_point(alpha = 0.7, size = 2) +
        scale_color_manual(name = "Sex", values = c("Male" = "steelblue", "Female" = "red")) +
        ggtitle(sprintf("Expression of %s in %s", input$gene, input$tissue)) +
        xlab("Age") + ylab("logTPM") + theme_minimal()
    }
    
  })
  
  # Function to load the p-value table
  load_pvalue_table <- function(tissue) {
    output_path <- file.path(here::here("gtex_v10_shiny/data"), paste0(gsub(" ", "_", tissue), "_pvalue_results.csv"))
    if (!file.exists(output_path)) {
      return(NULL)
    }
    read.csv(output_path)
  }
  
  # Render the p-value table
  output$pvalueTable <- DT::renderDataTable({
    req(input$tissue)
    pvalue_data <- load_pvalue_table(input$tissue)
    
    validate(need(!is.null(pvalue_data), "P-Value table not found for the selected tissue."))
    
    DT::datatable(pvalue_data, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  # Download handler
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste0(input$tissue, "_pvalue_results.csv")
    },
    content = function(file) {
      pvalue_data <- load_pvalue_table(input$tissue)
      if (is.null(pvalue_data)) {
        stop("P-Value table not found.")
      }
      write.csv(pvalue_data, file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)