# =============================================================================
# mod_qc_imputation.R
# Module: "QC & Imputation" tab
#
# Converted from BiomiX_preview_Imputation.r (standalone Shiny app) into a
# Shiny module so it can live as a tab inside the main BiomiX interface.
#
# All original logic is preserved exactly as-is. Changes made:
#   - Wrapped UI and server into moduleUI / moduleServer functions
#   - Removed: renv::load(), get_default_browser(), shinyApp() launcher,
#     stopApp() / closeApp button (not needed inside a tab)
#   - Namespaced all input/output IDs with ns()
#
# Dependencies (add to Dockerfile / install.packages if missing):
#   DT, ggplot2, mice, mixOmics, data.table, readxl, tidyr, dplyr
#   Optional: randomForest (for Random Forest imputation)
#             glmnet       (for Lasso imputation)
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(mice)
library(mixOmics)
library(data.table)
library(readxl)
library(DT)

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------

qcImputationUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    
    # Page title (styled to match BiomiX theme)
    div(
      class = "section-label",
      style = "font-size: 15px; margin: 8px 0 16px 0;",
      icon("flask"), " QC & Imputation"
    ),
    
    sidebarLayout(
      
      # -----------------------------------------------------------------------
      # Sidebar: controls
      # -----------------------------------------------------------------------
      sidebarPanel(
        
        # Imputation method
        selectInput(
          ns("imputation_method"),
          "Choose Imputation Method:",
          choices = c(
            "Replace with 0",
            "Mean",
            "Median",
            "Lasso",
            "Random Forest",
            "NIPALS",
            "Replace NA with 0"
          )
        ),
        
        # Missingness thresholds
        numericInput(ns("sample_threshold"),
                     "Variable Missingness Threshold (%)",
                     value = 50, min = 0, max = 100),
        
        numericInput(ns("var_threshold"),
                     "Sample Missingness Threshold (%)",
                     value = 50, min = 0, max = 100),
        
        br(),
        
        # Action buttons
        actionButton(ns("Impute"), "Impute the matrix",
                     icon = icon("wand-magic-sparkles"), class = "btn-primary",
                     style = "width:100%; margin-bottom:8px;"),
        
        actionButton(ns("Filter"), "Filter Variables and Samples",
                     icon = icon("filter"),
                     style = "width:100%; margin-bottom:8px;"),
        
        br(),
        
        actionButton(ns("resetBtn"), "Reset to Original Data",
                     icon = icon("rotate-left"),
                     style = "width:100%;")
      ),
      
      # -----------------------------------------------------------------------
      # Main panel: tabset (identical structure to original)
      # -----------------------------------------------------------------------
      mainPanel(
        tabsetPanel(
          
          # Tab 1: Upload
          tabPanel(
            "Upload Data",
            br(),
            fileInput(ns("file"), "Upload a File",
                      accept = c(".csv", ".tsv", ".xls", ".xlsx")),
            checkboxInput(ns("transpose"),         "Transpose Data",        value = FALSE),
            checkboxInput(ns("remove_first_row"),  "Remove First Row",      value = FALSE),
            checkboxInput(ns("remove_first_col"),  "Remove First Column",   value = FALSE),
            actionButton(ns("load_data"), "Load Data",
                         icon = icon("upload"), class = "btn-primary"),
            br(), br(),
            uiOutput(ns("warningMessage"))
          ),
          
          # Tab 2: Boxplot by variables
          tabPanel("Boxplot by Variables",
                   br(), plotOutput(ns("barplot"), height = "500px")),
          
          # Tab 3: Barplot by samples
          tabPanel("Barplot by Samples",
                   br(), plotOutput(ns("boxplot"), height = "500px")),
          
          # Tab 4: Summary table
          tabPanel("Summary Table",
                   br(), DT::dataTableOutput(ns("summary_table"))),
          
          # Tab 5: Full data table
          tabPanel("Data Table",
                   br(), DT::dataTableOutput(ns("dataTable"))),
          
          # Tab 6: Download
          tabPanel(
            "Download",
            br(),
            div(
              style = "padding: 16px;",
              p("Download the filtered / imputed matrix as a TSV file."),
              downloadButton(ns("downloadData"), "Download processed data",
                             icon = icon("download"))
            )
          )
        )
      )
    )
  )
}


# -----------------------------------------------------------------------------
# IMPUTATION FUNCTION (unchanged from original)
# -----------------------------------------------------------------------------

.impute_data <- function(data, method) {
  df <- data
  df <- t(df)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  if (method == "Replace with 0") {
    df[is.na(df)] <- 0
    
  } else if (method == "Mean") {
    df <- df %>%
      mutate_all(~ifelse(is.na(.), mean(., na.rm = TRUE), .))
    
  } else if (method == "Median") {
    df <- df %>%
      mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))
    
  } else if (method == "Lasso") {
    pred <- quickpred(df, mincor = 0.4, minpuc = 0.5)
    imputed_data <- mice(df, method = "lasso.norm", m = 1,
                         predictorMatrix = pred, maxit = 5, printFlag = TRUE)
    df <- complete(imputed_data)
    
  } else if (method == "Random Forest") {
    df <- missForest(df)$ximp
    
  } else if (method == "NIPALS") {
    if (any(is.na(df))) {
      X.na <- as.matrix(df)
      df   <- impute.nipals(X = X.na, ncomp = 10)
    } else {
      message("No missing values detected, skipping NIPALS imputation.")
    }
    
  } else if (method == "Replace NA with 0") {
    df[is.na(df)] <- 0
  }
  
  df <- t(df)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  return(df)
}


# -----------------------------------------------------------------------------
# SERVER
# -----------------------------------------------------------------------------

qcImputationServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive state: current working data + original snapshot for reset
    combined_data <- reactiveValues(numeric = NULL)
    original_data <- reactiveValues(numeric = NULL)
    
    # -------------------------------------------------------------------------
    # Load data (unchanged logic from original)
    # -------------------------------------------------------------------------
    observeEvent(input$load_data, {
      req(input$file)
      
      file_ext <- tools::file_ext(input$file$name)
      
      df <- tryCatch({
        if (file_ext == "csv") {
          data.table::fread(input$file$datapath, header = TRUE)
        } else if (file_ext == "tsv") {
          data.table::fread(input$file$datapath, header = TRUE, sep = "\t")
        } else if (file_ext %in% c("xlsx", "xls")) {
          readxl::read_excel(input$file$datapath, col_names = TRUE)
        } else {
          showNotification(
            "Unsupported file format. Please upload CSV, TSV, or Excel files.",
            type = "error"
          )
          return()
        }
      }, error = function(e) {
        showNotification(paste("Error reading file:", e$message), type = "error")
        NULL
      })
      
      if (is.null(df)) return()
      
      # Fix NA column names
      rows_with_na <- which(is.na(colnames(df)))
      if (length(rows_with_na) > 1)
        colnames(df)[rows_with_na] <- paste("missing", seq_along(rows_with_na), sep = "_")
      colnames(df) <- make.unique(colnames(df), sep = "_")
      
      # Transpose if requested
      if (input$transpose) {
        df <- t(df)
        colnames(df) <- as.character(unlist(df[1, ]))
        df <- df[-1, ]
        df <- as.data.frame(df, stringsAsFactors = FALSE)
      }
      
      # Remove first row if requested
      if (input$remove_first_row)
        df <- df[-1, ]
      
      # Remove first column if requested (use it as rownames)
      if (input$remove_first_col) {
        sav <- df[, 1]
        df  <- df[, -1]
        df  <- as.data.frame(df)
        sav <- as.vector(sav[[1]])
        rows_with_na <- which(is.na(sav))
        if (length(rows_with_na) > 1)
          sav[rows_with_na] <- paste("missing", seq_along(rows_with_na), sep = "_")
        sav <- make.unique(sav, sep = "_")
        row.names(df) <- as.vector(sav)
      }
      
      combined_data$numeric <- df
      original_data$numeric <- df
      
      showNotification(
        paste0("Data loaded: ", nrow(df), " variables × ", ncol(df), " samples"),
        type = "message", duration = 4
      )
    })
    
    # -------------------------------------------------------------------------
    # Warning: non-numeric columns (unchanged)
    # -------------------------------------------------------------------------
    output$warningMessage <- renderUI({
      req(combined_data$numeric)
      plot_data  <- combined_data$numeric
      is_numeric <- sapply(plot_data, is.numeric)
      if (any(!is_numeric)) {
        warning_cols <- paste(names(plot_data)[!is_numeric], collapse = ", ")
        div(
          style = "color: red; font-weight: bold;",
          paste("Warning: Non-numeric values found in columns:", warning_cols)
        )
      } else {
        NULL
      }
    })
    
    # -------------------------------------------------------------------------
    # Barplot: missing/zero by samples (unchanged)
    # -------------------------------------------------------------------------
    output$boxplot <- renderPlot({
      req(combined_data$numeric)
      df            <- combined_data$numeric
      total_samples <- nrow(df)
      
      stats <- df %>%
        summarise_all(~sum(is.na(.) | . == 0)) %>%
        tidyr::pivot_longer(everything(),
                            names_to  = "samples",
                            values_to = "count_missing") %>%
        mutate(percent_missing = (count_missing / total_samples) * 100)
      
      ggplot(stats, aes(x = samples, y = percent_missing)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        labs(
          title = "Percentage of Missing/Zero Values by Samples",
          y     = "% Missing or Zero",
          x     = "Samples"
        ) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    })
    
    # -------------------------------------------------------------------------
    # Boxplot: missing/zero by variables (unchanged)
    # -------------------------------------------------------------------------
    output$barplot <- renderPlot({
      req(combined_data$numeric)
      df        <- combined_data$numeric
      total_vars <- ncol(df)
      variable  <- rownames(df)
      
      stats <- df %>%
        rowwise() %>%
        mutate(count_missing = sum(is.na(c_across()) | c_across() == 0)) %>%
        ungroup() %>%
        mutate(percent_missing = (count_missing / total_vars) * 100)
      
      stats <- stats[, c(ncol(stats) - 1, ncol(stats))]
      rownames(stats) <- variable
      
      ggplot(stats, aes(x = factor(rownames(stats)), y = percent_missing)) +
        geom_bar(stat = "identity", fill = "darkorange") +
        labs(
          title = "Percentage of Missing/Zero Values by Variables",
          y     = "% Missing or Zero",
          x     = "Variables"
        ) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    })
    
    # -------------------------------------------------------------------------
    # Summary table — max 50 rows/cols (unchanged)
    # -------------------------------------------------------------------------
    output$summary_table <- DT::renderDataTable({
      req(combined_data$numeric)
      plot_data <- combined_data$numeric
      
      if (ncol(plot_data) > 50) {
        showNotification(
          "More than 50 variables detected. Displaying only the first 50.",
          type = "warning", duration = 5
        )
        plot_data <- plot_data[, 1:50]
      }
      if (nrow(plot_data) > 50) {
        showNotification(
          "More than 50 samples detected. Displaying only the first 50.",
          type = "warning", duration = 5
        )
        plot_data <- plot_data[1:50, ]
      }
      plot_data
    })
    
    # -------------------------------------------------------------------------
    # Full data table (unchanged)
    # -------------------------------------------------------------------------
    output$dataTable <- DT::renderDataTable({
      req(combined_data$numeric)
      combined_data$numeric
    })
    
    # -------------------------------------------------------------------------
    # Reset to original data (unchanged)
    # -------------------------------------------------------------------------
    observeEvent(input$resetBtn, {
      req(original_data$numeric)
      combined_data$numeric <- original_data$numeric
      showNotification("Data has been reset to its original state.",
                       type = "message", duration = 3)
    })
    
    # -------------------------------------------------------------------------
    # Filter variables and samples by missingness threshold (unchanged)
    # -------------------------------------------------------------------------
    observeEvent(input$Filter, {
      req(combined_data$numeric)
      df            <- combined_data$numeric
      total_samples <- nrow(df)
      
      # Sample-level filtering
      stats <- df %>%
        summarise_all(~sum(is.na(.) | . == 0)) %>%
        tidyr::pivot_longer(everything(),
                            names_to  = "samples",
                            values_to = "count_missing") %>%
        mutate(percent_missing = (count_missing / total_samples) * 100)
      
      sam_to_keep <- stats %>%
        filter(percent_missing <= input$sample_threshold) %>%
        pull(samples)
      
      df_filtered <- df[, sam_to_keep, drop = FALSE]
      
      # Variable-level filtering
      total_vars <- ncol(df)
      variable   <- rownames(df)
      
      var_stats <- df %>%
        rowwise() %>%
        mutate(count_missing = sum(is.na(c_across()) | c_across() == 0)) %>%
        ungroup() %>%
        mutate(percent_missing = (count_missing / total_vars) * 100)
      
      var_stats <- var_stats[, c(ncol(var_stats) - 1, ncol(var_stats))]
      rownames(var_stats) <- variable
      
      vars_to_keep  <- which(var_stats$percent_missing <= input$var_threshold)
      df_filtered   <- df_filtered[vars_to_keep, , drop = FALSE]
      
      combined_data$numeric <- df_filtered
      
      showNotification(
        paste0("Filtered: ", nrow(df_filtered), " variables × ",
               ncol(df_filtered), " samples remaining."),
        type = "message", duration = 4
      )
    })
    
    # -------------------------------------------------------------------------
    # Imputation (unchanged)
    # -------------------------------------------------------------------------
    observeEvent(input$Impute, {
      req(combined_data$numeric)
      method <- input$imputation_method
      
      withProgress(message = paste("Imputing with", method, "..."), value = 0.5, {
        combined_data$numeric <- .impute_data(combined_data$numeric, method)
      })
      
      showNotification(
        paste0("Imputation complete (", method, ")."),
        type = "message", duration = 4
      )
    })
    
    # -------------------------------------------------------------------------
    # Download processed data as TSV (unchanged)
    # -------------------------------------------------------------------------
    output$downloadData <- downloadHandler(
      filename    = function() paste0("imputed_data_", Sys.Date(), ".tsv"),
      content     = function(file) {
        saved <- cbind(rownames(combined_data$numeric), combined_data$numeric)
        colnames(saved)[1] <- "ID"
        write.table(saved, file, row.names = FALSE, quote = FALSE, sep = "\t")
      },
      contentType = "text/tab-separated-values"
    )
    
  })
}