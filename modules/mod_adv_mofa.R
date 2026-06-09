# =============================================================================
# mod_adv_mofa.R
# Modulo: tab "MOFA & Diablo" delle opzioni avanzate
#
# Contiene in un unico tab:
#   1. MOFA parameters:   max iterations, convergence mode, contribution threshold
#   2. Diablo parameters: feature selection, N iterations, design matrix file
#   3. Factor interpretation (condivisa MOFA + Diablo):
#        - Bibliography: search type, N articles, N features, N contributors
#        - Pathway:      padj threshold, N pathways
#        - Clinical:     numerical/binary + file upload
#
# Il metadata filtering è nel tab separato "Metadata Filtering".
#
# Output esposti: adv_mofa_data() → lista con tutti i parametri
# =============================================================================

advMofaUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    
    # =========================================================================
    # SEZIONE MOFA
    # =========================================================================
    .adv_section_header("MOFA+ algorithm parameters"),
    
    fluidRow(
      column(3,
             .adv_param_block("Max iterations",
                              selectInput(ns("mofa_max_iter"), NULL,
                                          choices  = c("1000","2000","5000","10000"),
                                          selected = "5000", width = "100%")
             )
      ),
      column(3,
             .adv_param_block("Convergence mode",
                              selectInput(ns("mofa_convergence"), NULL,
                                          choices  = c("fast","medium","slow"),
                                          selected = "fast", width = "100%")
             )
      ),
      column(3,
             .adv_param_block("Contribution threshold",
                              numericInput(ns("mofa_threshold"), NULL,
                                           value = 0.5, min = 0, max = 1, step = 0.05, width = "100%")
             )
      )
    ),
    
    # =========================================================================
    # SEZIONE DIABLO
    # =========================================================================
    .adv_section_header("Diablo parameters"),
    
    fluidRow(
      column(3,
             .adv_param_block("Feature selection",
                              selectInput(ns("diablo_feature_sel"), NULL,
                                          choices  = c("NO","YES"),
                                          selected = "NO", width = "100%")
             )
      ),
      column(3,
             .adv_param_block("N° iterations",
                              numericInput(ns("diablo_n_iter"), NULL,
                                           value = 99, min = 1, step = 1, width = "100%")
             )
      ),
      column(6,
             .adv_param_block("Design matrix file (optional)",
                              div(style = "display:flex; gap:10px; align-items:center;",
                                  shinyFilesButton(ns("diablo_design_btn"), "Browse...",
                                                   "Select Diablo design matrix file",
                                                   multiple = FALSE, class = "btn-browse-sm"),
                                  uiOutput(ns("diablo_design_label"))
                              )
             )
      )
    ),
    
    # =========================================================================
    # SEZIONE FACTOR INTERPRETATION — condivisa MOFA + Diablo
    # =========================================================================
    div(
      style = "margin-top: 24px; padding: 16px; background: var(--biomix-blue-light);
               border-radius: 8px; border-left: 4px solid var(--biomix-blue);",
      div(
        style = "font-family: 'Space Mono', monospace; font-size: 11px; font-weight: 700;
                 text-transform: uppercase; letter-spacing: 0.08em;
                 color: var(--biomix-blue); margin-bottom: 12px;",
        icon("magnifying-glass-chart"),
        " Factor Interpretation — used by both MOFA and Diablo"
      ),
      
      # ---- Bibliography ---------------------------------------------------
      .adv_section_header("Bibliography"),
      
      fluidRow(
        column(3,
               .adv_param_block("Search type",
                                selectInput(ns("biblio_type"), NULL,
                                            choices  = c("Title/Abstract","Title","Abstract"),
                                            selected = "Title/Abstract", width = "100%")
               )
        ),
        column(3,
               .adv_param_block("N° articles retrieved",
                                selectInput(ns("biblio_n_articles"), NULL,
                                            choices  = c("10","50","100","200","300"),
                                            selected = "10", width = "100%")
               )
        ),
        column(3,
               .adv_param_block("N° top Keywords extracted",
                                selectInput(ns("biblio_n_keywords"), NULL,
                                            choices  = c("100","500","1000","2000","5000"),
                                            selected = "10", width = "100%")
               )
        ),
        column(3,
               .adv_param_block("N° top contributors",
                                numericInput(ns("biblio_top_contrib"), NULL,
                                             value = 10, min = 1, step = 1, width = "100%")
               )
        )
      ),
      
      # ---- Pathway --------------------------------------------------------
      .adv_section_header("Pathway analysis"),
      
      fluidRow(
        column(3,
               .adv_param_block("P.adj threshold",
                                numericInput(ns("pathway_padj"), NULL,
                                             value = 0.05, min = 0, max = 1, step = 0.01, width = "100%")
               )
        ),
        column(3,
               .adv_param_block("N° pathways shown",
                                numericInput(ns("pathway_n"), NULL,
                                             value = 10, min = 1, step = 1, width = "100%")
               )
        )
      ),
      
      # ---- Clinical -------------------------------------------------------
      .adv_section_header("Clinical data integration"),
      
      fluidRow(
        column(3,
               .adv_param_block("Data types to correlate",
                                tagList(
                                  checkboxInput(ns("clinical_numerical"), "Numerical metadata", value = FALSE),
                                  checkboxInput(ns("clinical_binary"),    "Binary metadata",    value = FALSE)
                                )
               )
        ),
        column(4,
               .adv_param_block("Numerical clinical data file (optional)",
                                div(style = "display:flex; gap:10px; align-items:center;",
                                    shinyFilesButton(ns("clinic_numerical_btn"), "Browse...",
                                                     "Select numerical clinical data file",
                                                     multiple = FALSE, class = "btn-browse-sm"),
                                    uiOutput(ns("clinic_numerical_label"))
                                )
               )
        ),
        column(4,
               .adv_param_block("Binary clinical data file (optional)",
                                div(style = "display:flex; gap:10px; align-items:center;",
                                    shinyFilesButton(ns("clinic_binary_btn"), "Browse...",
                                                     "Select binary clinical data file",
                                                     multiple = FALSE, class = "btn-browse-sm"),
                                    uiOutput(ns("clinic_binary_label"))
                                )
               )
        )
      )
    ) # fine box Factor Interpretation
  )
}


advMofaServer <- function(id, roots) {
  moduleServer(id, function(input, output, session) {
    
    clinic_file_numerical_path <- reactiveVal(EMPTY_VALUE)
    clinic_file_binary_path    <- reactiveVal(EMPTY_VALUE)
    diablo_design_path         <- reactiveVal(EMPTY_VALUE)
    
    # Numerical clinical file
    shinyFileChoose(input, "clinic_numerical_btn", roots = roots, session = session,
                    filetypes = c("tsv","csv","txt","xlsx","xls"))
    observeEvent(input$clinic_numerical_btn, {
      info <- parseFilePaths(roots, input$clinic_numerical_btn)
      if (nrow(info) > 0) clinic_file_numerical_path(as.character(info$datapath[1]))
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    output$clinic_numerical_label <- renderUI({
      p <- clinic_file_numerical_path()
      if (p == EMPTY_VALUE) span(class = "file-label-empty", "No file (optional)")
      else span(class = "file-label-ok", icon("circle-check"), " ", basename(p))
    })
    
    # Binary clinical file
    shinyFileChoose(input, "clinic_binary_btn", roots = roots, session = session,
                    filetypes = c("tsv","csv","txt","xlsx","xls"))
    observeEvent(input$clinic_binary_btn, {
      info <- parseFilePaths(roots, input$clinic_binary_btn)
      if (nrow(info) > 0) clinic_file_binary_path(as.character(info$datapath[1]))
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    output$clinic_binary_label <- renderUI({
      p <- clinic_file_binary_path()
      if (p == EMPTY_VALUE) span(class = "file-label-empty", "No file (optional)")
      else span(class = "file-label-ok", icon("circle-check"), " ", basename(p))
    })
    
    # Diablo design matrix file
    shinyFileChoose(input, "diablo_design_btn", roots = roots, session = session,
                    filetypes = c("tsv","csv","txt","xlsx","xls"))
    observeEvent(input$diablo_design_btn, {
      info <- parseFilePaths(roots, input$diablo_design_btn)
      if (nrow(info) > 0) diablo_design_path(as.character(info$datapath[1]))
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    output$diablo_design_label <- renderUI({
      p <- diablo_design_path()
      if (p == EMPTY_VALUE) span(class = "file-label-empty", "Nessun file (opzionale)")
      else span(class = "file-label-ok", icon("circle-check"), " ", basename(p))
    })
    
    reactive({
      list(
        # MOFA
        mofa_max_iter    = input$mofa_max_iter    %||% "5000",
        mofa_convergence = input$mofa_convergence %||% "fast",
        mofa_threshold   = input$mofa_threshold   %||% 0.5,
        
        # Diablo
        diablo_feature_sel = input$diablo_feature_sel %||% "NO",
        diablo_n_iter      = input$diablo_n_iter      %||% 99L,
        diablo_design_file = diablo_design_path(),
        
        # Bibliography
        biblio_type        = input$biblio_type        %||% "Title/Abstract",
        biblio_n_articles  = input$biblio_n_articles  %||% "10",
        biblio_n_keywords  = input$biblio_n_keywords  %||% "10",
        biblio_top_contrib = input$biblio_top_contrib %||% 10L,
        
        # Pathway
        pathway_padj = input$pathway_padj %||% 0.05,
        pathway_n    = input$pathway_n    %||% 10L,
        
        # Clinical
        clinical_numerical      = if (isTRUE(input$clinical_numerical)) "YES" else "NO",
        clinical_binary         = if (isTRUE(input$clinical_binary))    "YES" else "NO",
        clinic_file_numerical   = clinic_file_numerical_path(),
        clinic_file_binary      = clinic_file_binary_path()
      )
    })
    
  })
}