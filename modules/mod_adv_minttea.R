# =============================================================================
# mod_adv_minttea.R
# Module: "MintTea" tab of the Advanced Options panel
#
# Sections:
#   1. Cross-validation: n_repeats, n_folds
#   2. Network:          edge_threshold
#   3. Feature selection: keepX (features per component)
#   4. Design:           design_val (off-diagonal value in the design matrix)
#
# Returned value: adv_minttea_data() -> list with all parameters
# =============================================================================

advMintteaUI <- function(id) {
  ns <- NS(id)

  tagList(

    .adv_section_header("MintTea — Cross-validation"),

    fluidRow(
      column(3,
        .adv_param_block("N° repeats",
          numericInput(ns("n_repeats"), NULL, value = 10, min = 1, step = 1, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("N° folds",
          numericInput(ns("n_folds"), NULL, value = 10, min = 2, step = 1, width = "100%")
        )
      )
    ),

    .adv_section_header("MintTea — Network"),

    fluidRow(
      column(3,
        .adv_param_block("Edge threshold",
          numericInput(ns("edge_threshold"), NULL,
                       value = 0.5, min = 0, max = 1, step = 0.05, width = "100%")
        )
      )
    ),

    .adv_section_header("MintTea — Feature selection"),

    fluidRow(
      column(3,
        .adv_param_block("keepX (features per component)",
          numericInput(ns("keep_x"), NULL, value = 20, min = 1, step = 1, width = "100%")
        )
      ),
      column(3,
        .adv_param_block("Design value (off-diagonal)",
          numericInput(ns("design_val"), NULL,
                       value = 0.5, min = 0, max = 1, step = 0.05, width = "100%")
        )
      )
    )
  )
}


advMintteaServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    reactive({
      list(
        n_repeats      = input$n_repeats      %||% 10L,
        n_folds        = input$n_folds        %||% 10L,
        edge_threshold = input$edge_threshold %||% 0.5,
        keep_x         = input$keep_x         %||% 20L,
        design_val     = input$design_val     %||% 0.5
      )
    })

  })
}
