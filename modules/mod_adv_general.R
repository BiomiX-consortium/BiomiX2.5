# =============================================================================
# mod_adv_general.R
# Module: "General" tab of the Advanced Options panel (SecondWindow -> tab_13)
#
# Sections:
#   1. Transcriptomics: Log2FC threshold, padj threshold, optional gene panel file
#   2. Metabolomics:    Log2FC threshold, padj threshold
#   3. Methylomics:     padj threshold, Log2FC threshold, array type
#   4. Subgrouping:     N genes Z-score>2, N genes Z-score>1, N top DE genes heatmap
#   5. Clustering:      distance method, linkage method
#   6. Computational:   CPU threads, N° input features for integration
#
# Returned value: adv_general_data() -> list with all parameters
#
# NOTE: this file also defines the shared UI helpers (.adv_section_header,
#       .adv_param_block) used by all other Advanced modules.
# =============================================================================

advGeneralUI <- function(id) {
  ns <- NS(id)

  tagList(

    # ---- 1. Transcriptomics ------------------------------------------------
    .adv_section_header("Transcriptomics Options"),

    fluidRow(
      column(3,
        .adv_param_block("Log2FC threshold",
          numericInput(ns("tx_fc"), NULL, value = 0.5, min = 0, step = 0.25, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("P.adj threshold",
          numericInput(ns("tx_padj"), NULL, value = 0.05, min = 0, max = 1, step = 0.01, width = "100%")
        )
      ),
      column(6,
        .adv_param_block("Gene panel file (optional)",
          div(style = "display:flex; gap:10px; align-items:center;",
            shinyFilesButton(ns("tx_panel_btn"), "Browse...", "Select gene panel file",
                             multiple = FALSE, class = "btn-browse-sm"),
            uiOutput(ns("tx_panel_label"))
          )
        )
      )
    ),

    # ---- 2. Metabolomics ---------------------------------------------------
    .adv_section_header("Metabolomics Options"),

    fluidRow(
      column(3,
        .adv_param_block("Log2FC threshold",
          numericInput(ns("met_fc"), NULL, value = 0.5, min = 0, step = 0.25, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("P.adj threshold",
          numericInput(ns("met_padj"), NULL, value = 0.05, min = 0, max = 1, step = 0.01, width = "100%")
        )
      )
    ),

    # ---- 3. Methylomics ----------------------------------------------------
    .adv_section_header("Methylomics Options"),

    fluidRow(
      column(3,
        .adv_param_block("P.adj threshold",
          numericInput(ns("meth_padj"), NULL, value = 0.05, min = 0, max = 1, step = 0.01, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("Log2FC threshold",
          numericInput(ns("meth_fc"), NULL, value = 0.15, min = 0, step = 0.05, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("Array type",
          selectInput(ns("meth_array"), NULL,
                      choices = c("450K", "EPIC"), selected = "450K", width = "100%")
        )
      )
    ),

    # ---- 4. Subgrouping / Gene panel positivity criteria -------------------
    .adv_section_header("Gene panel and positivity criteria"),

    fluidRow(
      column(3,
        .adv_param_block("N° genes Z-score > 2",
          numericInput(ns("sub_zscore2"), NULL, value = 10, min = 1, step = 1, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("N° genes Z-score > 1",
          numericInput(ns("sub_zscore1"), NULL, value = 10, min = 1, step = 1, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("N° top DE genes in heatmap",
          numericInput(ns("sub_heatmap"), NULL, value = 20, min = 1, step = 1, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("Subgrouping",
          checkboxInput(ns("sub_enable"), "Enable subgrouping", value = FALSE)
        )
      )
    ),

    # ---- 5. Clustering -----------------------------------------------------
    .adv_section_header("Clustering Options"),

    fluidRow(
      column(4,
        .adv_param_block(
          tagList("Clustering distance",
                  tags$a(href = "https://ixi-97.github.io/Parameters_guidelines.html#Clustering_Distance_Recommendation",
                         target = "_blank", icon("circle-question"),
                         style = "margin-left:6px; font-size:11px;")),
          selectInput(ns("clust_distance"), NULL,
                      choices  = c("euclidean","maximum","manhattan","canberra",
                                   "binary","minkowski","pearson","spearman"),
                      selected = "euclidean", width = "100%")
        )
      ),
      column(4,
        .adv_param_block(
          tagList("Clustering method",
                  tags$a(href = "https://ixi-97.github.io/Parameters_guidelines.html#Clustering_Methods_Recommendation",
                         target = "_blank", icon("circle-question"),
                         style = "margin-left:6px; font-size:11px;")),
          selectInput(ns("clust_method"), NULL,
                      choices  = c("ward.D2","ward.D","single","complete",
                                   "average","mcquitty","median","centroid"),
                      selected = "ward.D2", width = "100%")
        )
      )
    ),

    # ---- 6. Computational resources ----------------------------------------
    .adv_section_header("Computational resources / General integration options"),

    fluidRow(
      column(3,
        .adv_param_block("CPU threads",
          numericInput(ns("cpu_threads"), NULL, value = 3, min = 1, max = 128, step = 1, width = "100%")
        )
      ),
      column(4,
        .adv_param_block("N° input features for integration",
          selectInput(ns("mofa_features"), NULL,
                      choices  = c("500","1000","2000","3000","5000","10000"),
                      selected = "5000", width = "100%")
        )
      )
    )
  )
}


advGeneralServer <- function(id, roots) {
  moduleServer(id, function(input, output, session) {

    panel_path <- reactiveVal(EMPTY_VALUE)

    shinyFileChoose(input, "tx_panel_btn", roots = roots, session = session,
                    filetypes = c("tsv","csv","txt","xlsx","xls","gmt"))

    observeEvent(input$tx_panel_btn, {
      info <- parseFilePaths(roots, input$tx_panel_btn)
      if (nrow(info) > 0) panel_path(as.character(info$datapath[1]))
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    output$tx_panel_label <- renderUI({
      p <- panel_path()
      if (p == EMPTY_VALUE)
        span(class = "file-label-empty", "No file (optional)")
      else
        span(class = "file-label-ok", icon("circle-check"), " ", basename(p))
    })

    reactive({
      list(
        tx_fc         = input$tx_fc         %||% 0.5,
        tx_padj       = input$tx_padj       %||% 0.05,
        tx_panel      = panel_path(),
        met_fc        = input$met_fc        %||% 0.5,
        met_padj      = input$met_padj      %||% 0.05,
        meth_fc       = input$meth_fc       %||% 0.15,
        meth_padj     = input$meth_padj     %||% 0.05,
        meth_array    = input$meth_array    %||% "450K",
        sub_zscore2   = input$sub_zscore2   %||% 10L,
        sub_zscore1   = input$sub_zscore1   %||% 10L,
        sub_heatmap   = input$sub_heatmap   %||% 20L,
        sub_enable    = if (isTRUE(input$sub_enable)) "YES" else "NO",
        clust_distance = input$clust_distance %||% "euclidean",
        clust_method   = input$clust_method   %||% "ward.D2",
        cpu_threads    = input$cpu_threads    %||% 3L,
        mofa_features  = input$mofa_features  %||% "5000"
      )
    })

  })
}


# =============================================================================
# Shared UI helpers — available to all Advanced modules after source()
# =============================================================================

# Section header with lateral divider lines (matches PyQt style)
.adv_section_header <- function(title) {
  div(
    class = "adv-section-header",
    tags$hr(class = "adv-hr-left"),
    span(class = "adv-section-title", title),
    tags$hr(class = "adv-hr-right")
  )
}

# Label + widget block with consistent styling
.adv_param_block <- function(label, widget) {
  div(
    class = "adv-param-block",
    div(class = "adv-param-label", label),
    widget
  )
}
