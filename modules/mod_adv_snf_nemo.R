# =============================================================================
# mod_adv_snf_nemo.R
# Module: "SNF / NEMO" tab of the Advanced Options panel
#
# Sections:
#   1. SNF parameters:  N neighbors, variance affinity matrix, N iterations
#   2. NEMO parameters: N neighbors, variance affinity matrix
#   3. Shared metadata variables: numerical columns, categorical columns,
#      condition column (all populated from the loaded metadata)
#   4. Numeric options: cluster range, N variables heatmap, N max features
#
# Returned value: adv_snf_nemo_data() -> list with all parameters
# =============================================================================

advSnfNemoUI <- function(id) {
  ns <- NS(id)

  tagList(

    # ---- 1. SNF parameters -------------------------------------------------
    .adv_section_header("SNF parameters"),

    fluidRow(
      column(3,
        .adv_param_block("N° neighbors",
          numericInput(ns("snf_neighbors"), NULL, value = 5, min = 2, step = 1, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("Variance affinity matrix",
          numericInput(ns("snf_variance"), NULL, value = 0.5, min = 0.01, max = 1, step = 0.05, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("N° iterations",
          numericInput(ns("snf_iterations"), NULL, value = 10, min = 1, step = 1, width = "100%")
        )
      )
    ),

    # ---- 2. NEMO parameters ------------------------------------------------
    .adv_section_header("NEMO parameters"),

    fluidRow(
      column(3,
        .adv_param_block("N° neighbors",
          numericInput(ns("nemo_neighbors"), NULL, value = 10, min = 2, step = 1, width = "100%")
        )
      ),
      column(3,
             .adv_param_block("Only common samples",
                              selectInput(ns("only_common_samples"), NULL,
                                           choices  = c("FALSE","TRUE"),
                                           selected = "TRUE", width = "100%")
             )
      )
    ),
    

    # ---- 3. Metadata variables (populated from loaded metadata) ------------
    .adv_section_header("Metadata variables for SNF / NEMO"),

    fluidRow(
      column(4,
        .adv_param_block("Columns for Variable Enrichment test",
          selectInput(ns("snf_enrich_cols"), NULL,
                      choices = NULL, multiple = TRUE, width = "100%")
        )
      ),
      column(4,
        .adv_param_block("Columns for Survival variable test columns",
          selectInput(ns("snf_surv_cols"), NULL,
                      choices = NULL, multiple = TRUE, width = "100%")
        )
      ),
      column(4,
        .adv_param_block("Condition column",
          selectInput(ns("snf_condition_col"), NULL,
                      choices = NULL, selected = NULL, width = "100%")
        )
      )
    ),

    # ---- 4. Clustering and feature options ---------------------------------
    .adv_section_header("Clustering and feature options"),

    fluidRow(
      column(3,
        .adv_param_block("Max cluster range",
          numericInput(ns("snf_cluster_range"), NULL, value = 4, min = 2, step = 1, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("N° variables in heatmap",
          numericInput(ns("snf_var_heatmap"), NULL, value = 20, min = 1, step = 1, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("Max N° features SNF / NEMO",
          numericInput(ns("snf_max_features"), NULL, value = 5000, min = 100, step = 100, width = "100%")
        )
      ),
      column(3,
             .adv_param_block("Scaling the SNF/NEMO data",
                              selectInput(ns("apply_scaling_SNF"), NULL,
                                           choices  = c("FALSE","TRUE"),
                                           selected = "TRUE", width = "100%")
             )
      ),
    )
  )
}


advSnfNemoServer <- function(id, metadata_df) {
  moduleServer(id, function(input, output, session) {

    # Populate column selectors when metadata is loaded
    observe({
      df <- metadata_df()
      if (!is.null(df)) {
        cols     <- colnames(df)
        num_cols <- cols[sapply(df, is.numeric)]
        cat_cols <- cols[!sapply(df, is.numeric)]
        updateSelectInput(session, "snf_enrich_cols",      choices = num_cols, selected = NULL)
        updateSelectInput(session, "snf_surv_cols",      choices = cat_cols, selected = NULL)
        updateSelectInput(session, "snf_condition_col", choices = cols,     selected = "CONDITION")
      }
    })

    reactive({
      # Selected columns -> "/" separated string (matches JSON format)
      num_cols <- paste(input$snf_enrich_cols %||% character(0), collapse = "/")
      cat_cols <- paste(input$snf_surv_cols %||% character(0), collapse = "/")

      list(
        snf_neighbors     = input$snf_neighbors     %||% 10L,
        snf_variance      = input$snf_variance      %||% 0.5,
        snf_iterations    = input$snf_iterations    %||% 10L,
        nemo_neighbors    = input$nemo_neighbors    %||% 10L,
        apply_scaling_SNF     = input$apply_scaling_SNF     %||% "TRUE",
        only_common_samples = input$only_common_samples %||% "TRUE",
        snf_max_features  = input$snf_max_features  %||% 5000L,
        snf_enrich_cols      = if (nchar(num_cols) > 0) num_cols else EMPTY_VALUE,
        snf_surv_cols      = if (nchar(cat_cols) > 0) cat_cols else EMPTY_VALUE,
        snf_condition_col = input$snf_condition_col %||% "CONDITION",
        snf_cluster_range = input$snf_cluster_range %||% 8L,
        snf_var_heatmap   = input$snf_var_heatmap   %||% 20L,
        snf_max_features  = input$snf_max_features  %||% 5000L
      )
    })

  })
}
