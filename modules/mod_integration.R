# =============================================================================
# mod_integration.R
# Module: "Integration Options" section of the main screen.
#
# Responsibilities:
#   - Integration toggle: Yes / No radio buttons
#   - If Yes: method selection (MOFA / Diablo / SNF / NEMO / MintTea)
#   - Numeric parameters shown conditionally by method:
#       MOFA, Diablo -> N° factors + Factor to explore + Omics overlay
#       SNF, NEMO, MintTea -> no numeric params here
#                             (their specific params are in Advanced Options)
#
# Returned value: integration_data() -> list with:
#                    $data_integration  ("YES"/"NO")
#                    $integration_type  (e.g. "MOFA" or "NO_INTEGRATION")
#                    $n_factors         (integer, 0 = auto)
#                    $factor_explore    (integer)
#                    $omics_overlay     (integer)
# =============================================================================

# Methods that expose numeric parameters on the main screen
METHODS_WITH_FACTORS <- c("MOFA", "Diablo")


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------

integrationUI <- function(id) {
  ns <- NS(id)

  div(
    class = "integration-panel",

    div(class = "panel-section-title", icon("circle-nodes"), " Integration Options"),

    fluidRow(

      # Column 1: Integration Yes/No
      column(3,
        div(class = "integration-block",
          div(class = "integration-block-label", "Integration"),
          radioButtons(
            ns("data_integration"),
            label    = NULL,
            choices  = c("Yes" = "YES", "No" = "NO"),
            selected = "NO",
            inline   = TRUE
          )
        )
      ),

      # Column 2: Method (visible only when Integration = YES)
      column(3,
        conditionalPanel(
          condition = paste0("input['", ns("data_integration"), "'] == 'YES'"),
          div(class = "integration-block",
            div(class = "integration-block-label", "Method"),
            radioButtons(
              ns("method"),
              label    = NULL,
              choices  = c("MOFA", "Diablo", "SNF", "NEMO", "MintTea"),
              selected = "MOFA",
              inline   = FALSE
            )
          )
        )
      ),

      # Column 3: Numeric parameters (MOFA and Diablo only)
      column(6,
        # Show parameters for MOFA / Diablo
        conditionalPanel(
          condition = paste0(
            "input['", ns("data_integration"), "'] == 'YES' && ",
            "(['MOFA','Diablo'].indexOf(input['", ns("method"), "']) !== -1)"
          ),
          div(class = "integration-block",
            div(class = "integration-block-label", "Parameters"),
            fluidRow(
              column(4,
                div(class = "param-label", "N° factors",
                    br(), span(style = "font-size:10px; color:#aaa;", "(0 = auto)")),
                numericInput(ns("n_factors"), label = NULL,
                             value = 0, min = 0, max = 100, step = 1, width = "100%")
              ),
              column(4,
                div(class = "param-label", "Factor to explore"),
                numericInput(ns("factor_explore"), label = NULL,
                             value = 1, min = 1, max = 100, step = 1, width = "100%")
              ),
              column(4,
                div(class = "param-label", "Omics overlay"),
                numericInput(ns("omics_overlay"), label = NULL,
                             value = 0, min = 0, max = 10, step = 1, width = "100%")
              )
            )
          )
        ),

        # Informational message for SNF / NEMO / MintTea
        conditionalPanel(
          condition = paste0(
            "input['", ns("data_integration"), "'] == 'YES' && ",
            "(['SNF','NEMO','MintTea'].indexOf(input['", ns("method"), "']) !== -1)"
          ),
          div(
            class = "integration-block",
            div(class = "integration-block-label", "Parameters"),
            div(
              class = "status-warn",
              style = "margin-top: 8px; font-size: 12px;",
              icon("circle-info"),
              " Method-specific parameters for this integration method are in ",
              strong("Advanced Options"), "."
            )
          )
        )
      )

    ) # end fluidRow
  )
}


# -----------------------------------------------------------------------------
# SERVER
# -----------------------------------------------------------------------------

integrationServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    reactive({
      is_integration <- input$data_integration == "YES"

      list(
        data_integration = input$data_integration %||% "NO",

        # If Integration = NO -> integration_type = "NO_INTEGRATION"
        integration_type = if (is_integration) {
          input$method %||% "MOFA"
        } else {
          "NO_INTEGRATION"
        },

        n_factors      = if (is_integration) as.integer(input$n_factors      %||% 0) else 0L,
        factor_explore = if (is_integration) as.integer(input$factor_explore %||% 1) else 1L,
        omics_overlay  = if (is_integration) as.integer(input$omics_overlay  %||% 0) else 0L
      )
    })

  })
}
