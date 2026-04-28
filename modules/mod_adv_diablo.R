# =============================================================================
# mod_adv_diablo.R
# Modulo: tab "Diablo" delle opzioni avanzate
#
# Sezioni:
#   1. Feature selection method
#   2. N° iterations
#   3. Design matrix file (opzionale)
#
# Output esposti: adv_diablo_data() → lista con tutti i parametri
# =============================================================================

advDiabloUI <- function(id) {
  ns <- NS(id)

  tagList(

    .adv_section_header("Diablo parameters"),

    fluidRow(
      column(3,
        .adv_param_block("Feature selection",
          selectInput(ns("feature_selection"), NULL,
                      choices  = c("NO", "YES"),
                      selected = "NO", width = "100%")
        )
      ),
      column(3,
        .adv_param_block("N° iterations",
          numericInput(ns("n_iter"), NULL, value = 99, min = 1, step = 1, width = "100%")
        )
      )
    ),

    .adv_section_header("Design matrix file (optional)"),

    fluidRow(
      column(8,
        .adv_param_block("Design matrix file",
          div(style = "display:flex; gap:10px; align-items:center;",
            shinyFilesButton(ns("design_file_btn"), "Browse...",
                             "Select design matrix file", multiple = FALSE,
                             class = "btn-browse-sm"),
            uiOutput(ns("design_file_label"))
          )
        )
      )
    )
  )
}


advDiabloServer <- function(id, roots) {
  moduleServer(id, function(input, output, session) {

    design_file_path <- reactiveVal(EMPTY_VALUE)

    shinyFileChoose(input, "design_file_btn", roots = roots, session = session,
                    filetypes = c("tsv", "csv", "txt", "xlsx", "xls"))

    observeEvent(input$design_file_btn, {
      info <- parseFilePaths(roots, input$design_file_btn)
      if (nrow(info) > 0) design_file_path(as.character(info$datapath[1]))
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    output$design_file_label <- renderUI({
      p <- design_file_path()
      if (p == EMPTY_VALUE) span(class = "file-label-empty", "Nessun file (opzionale)")
      else span(class = "file-label-ok", icon("circle-check"), " ", basename(p))
    })

    reactive({
      list(
        feature_selection = input$feature_selection %||% "NO",
        n_iter            = input$n_iter            %||% 99L,
        design_file       = design_file_path()
      )
    })

  })
}
