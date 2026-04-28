# =============================================================================
# mod_inputs.R
# Module: six omics input rows (main screen).
#
# Responsibilities:
#   - Renders one row per input (1–6) with columns:
#       Preview-QC | Single Omics Analysis | Data Type | Integration |
#       Label | Filter by Label | Data Upload
#   - Handles interactive file selection via shinyFiles for each input
#   - Turns the Upload button green when a file is selected
#
# Additional dependency: shinyFiles
#
# Input arguments: roots = named vector of accessible root directories (from app.R)
# Returned value:  inputs_data() -> list of 6 elements, each containing:
#                    $preview    ("YES"/"NO")
#                    $analysis   ("YES"/"NO")
#                    $data_type  (string, e.g. "Metabolomics")
#                    $integration("YES"/"NO")
#                    $label      (free text string)
#                    $selection  ("YES"/"NO")
#                    $path       (absolute file path or "X" if not set)
# =============================================================================

library(shinyFiles)


# -----------------------------------------------------------------------------
# Helper: generate a single input row for index i (1..6)
# -----------------------------------------------------------------------------

.input_row_ui <- function(ns, i) {
  fluidRow(
    class = paste0("input-row", if (i %% 2 == 0) " input-row-alt" else ""),

    # Row label
    column(1, div(class = "input-row-label", paste0("Input ", i))),

    # Preview-QC checkbox
    column(1, div(class = "checkbox-center",
      checkboxInput(ns(paste0("preview_", i)), label = NULL, value = FALSE)
    )),

    # Single Omics Analysis checkbox
    column(1, div(class = "checkbox-center",
      checkboxInput(ns(paste0("analysis_", i)), label = NULL, value = FALSE)
    )),

    # Data Type dropdown
    column(2,
      selectInput(
        ns(paste0("data_type_", i)),
        label   = NULL,
        choices = OMICS_TYPES,
        width   = "100%"
      )
    ),

    # Integration checkbox
    column(1, div(class = "checkbox-center",
      checkboxInput(ns(paste0("integration_", i)), label = NULL, value = FALSE)
    )),

    # Label text input
    column(2,
      textInput(
        ns(paste0("label_", i)),
        label       = NULL,
        placeholder = paste0("Label input ", i),
        width       = "100%"
      )
    ),

    # Filter by Label (Selection) checkbox
    column(1, div(class = "checkbox-center",
      checkboxInput(ns(paste0("selection_", i)), label = NULL, value = FALSE)
    )),

    # File upload — shinyFilesButton opens a native filesystem browser
    column(3,
      div(class = "upload-btn-wrapper",
        shinyFilesButton(
          id       = ns(paste0("upload_btn_", i)),
          label    = "Browse file...",
          title    = paste0("Select file for Input ", i),
          multiple = FALSE,
          icon     = icon("folder-open"),
          class    = "btn-upload"
        ),
        # Shows selected path or empty state
        uiOutput(ns(paste0("file_label_", i)))
      )
    )
  )
}


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------

inputsUI <- function(id) {
  ns <- NS(id)

  tagList(

    # Column headers
    fluidRow(
      class = "inputs-header",
      column(1, "Input"),
      column(1, div(class = "col-header", "Preview", br(), "QC")),
      column(1, div(class = "col-header", "Single Omics", br(), "Analysis")),
      column(2, div(class = "col-header", "Data Type")),
      column(1, div(class = "col-header", "Integration")),
      column(2, div(class = "col-header", "Label")),
      column(1, div(class = "col-header", "Filter by", br(), "Label")),
      column(3, div(class = "col-header", "Data Upload"))
    ),

    # 6 input rows
    lapply(1:6, function(i) .input_row_ui(ns, i))
  )
}


# -----------------------------------------------------------------------------
# SERVER
# -----------------------------------------------------------------------------

inputsServer <- function(id, roots) {
  moduleServer(id, function(input, output, session) {

    # Reactive state: absolute file paths for each input
    # Initialised to EMPTY_VALUE = "X"
    file_paths <- reactiveValues(
      p1 = EMPTY_VALUE, p2 = EMPTY_VALUE, p3 = EMPTY_VALUE,
      p4 = EMPTY_VALUE, p5 = EMPTY_VALUE, p6 = EMPTY_VALUE
    )

    # Factory: register shinyFiles + observer + label for input i
    make_upload_handler <- function(i) {
      field  <- paste0("p", i)
      btn_id <- paste0("upload_btn_", i)

      # Bind the button to the filesystem using the configured roots
      shinyFileChoose(
        input     = input,
        id        = btn_id,
        roots     = roots,
        session   = session,
        filetypes = c("tsv", "csv", "txt", "xlsx", "xls")
      )

      # When the user confirms a file selection
      observeEvent(input[[btn_id]], {
        # parseFilePaths converts the shinyFiles selection to a dataframe with $datapath
        info <- parseFilePaths(roots, input[[btn_id]])
        if (nrow(info) > 0)
          file_paths[[field]] <- as.character(info$datapath[1])
      }, ignoreNULL = TRUE, ignoreInit = TRUE)

      # Label below the button: show filename and full path
      output[[paste0("file_label_", i)]] <- renderUI({
        p <- file_paths[[field]]
        if (p == EMPTY_VALUE) {
          span(class = "file-label-empty", icon("file"), " No file selected")
        } else {
          tagList(
            span(class = "file-label-ok", icon("circle-check"), " ", basename(p)),
            br(),
            span(class = "file-label-path", p)
          )
        }
      })

      # Update button CSS class (blue -> green) via shinyjs when file is set
      observe({
        p <- file_paths[[field]]
        if (p != EMPTY_VALUE) {
          shinyjs::runjs(sprintf(
            "var el = document.getElementById('%s');
             if(el){ el.classList.add('btn-upload-done'); el.classList.remove('btn-upload'); }",
            session$ns(btn_id)
          ))
        }
      })
    }

    # Apply the factory to all 6 inputs
    lapply(1:6, make_upload_handler)

    # Value returned to the parent (app.R)
    reactive({
      lapply(1:6, function(i) {
        field <- paste0("p", i)
        list(
          preview     = if (isTRUE(input[[paste0("preview_",     i)]])) "YES" else "NO",
          analysis    = if (isTRUE(input[[paste0("analysis_",    i)]])) "YES" else "NO",
          data_type   = input[[paste0("data_type_",   i)]] %||% "Undefined",
          integration = if (isTRUE(input[[paste0("integration_", i)]])) "YES" else "NO",
          label       = input[[paste0("label_",       i)]] %||% "",
          selection   = if (isTRUE(input[[paste0("selection_",   i)]])) "YES" else "NO",
          path        = file_paths[[field]]
        )
      })
    })

  })
}
