# =============================================================================
# mod_startup.R
# Module: BiomiX startup screen.
#
# Responsibilities:
#   - Interactive metadata file selection via shinyFiles
#   - Read the CONDITION column and populate GROUP_1 / GROUP_2 dropdowns
#   - Interactive output folder selection via shinyDirButton
#
# Input arguments: roots = named vector of accessible root directories (from app.R)
# Returned value:  startup_data() -> reactive list with:
#                    $metadata_path  (string: absolute path)
#                    $metadata_df    (dataframe)
#                    $group1         (string)
#                    $group2         (string)
#                    $output_dir     (string: absolute path)
#                    $ready          (logical: TRUE when all fields are filled)
# =============================================================================

library(shinyFiles)


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------

startupUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    
    # --- Top banner ---------------------------------------------------------
    div(
      class = "startup-header",
      style = "display: flex; align-items: center; justify-content: space-between;",
      
      # Left side: title + subtitle
      div(
        div(class = "startup-title", "BiomiX",
            span(class = "startup-version", paste0("v", BIOMIX_VERSION))),
        div(class = "startup-subtitle", "Multiomics analysis pipeline")
      ),
      
      # Right side: logo (embedded as base64 via global.R — no path issues)
      div(style = "display: flex; align-items: center;", BIOMIX_LOGO_TAG)
    ),
    
    # --- Main card ----------------------------------------------------------
    div(
      class = "startup-card",
      
      # Section 1: Metadata file selection
      div(
        class = "startup-section",
        div(class = "section-label", "1. Metadata file"),
        div(class = "section-label-sub",
            "File must contain a 'CONDITION' column (.xlsx, .tsv, .csv)"),
        div(
          style = "display: flex; align-items: center; gap: 12px;",
          shinyFilesButton(
            id       = ns("metadata_file"),
            label    = "Select metadata file...",
            title    = "Select the metadata file",
            multiple = FALSE,
            icon     = icon("file-medical"),
            class    = "btn-browse"
          ),
          uiOutput(ns("metadata_path_display"))
        ),
        uiOutput(ns("metadata_status"))
      ),
      
      hr(class = "startup-divider"),
      
      # Section 2: Group selection (populated after metadata is loaded)
      div(
        class = "startup-section",
        div(class = "section-label", "2. Comparison groups"),
        fluidRow(
          column(6,
                 div(class = "group-label condition-label", "Condition (GROUP_1)"),
                 selectInput(ns("group1"), label = NULL, choices = NULL, width = "100%")
          ),
          column(6,
                 div(class = "group-label control-label", "Control (GROUP_2)"),
                 selectInput(ns("group2"), label = NULL, choices = NULL, width = "100%")
          )
        ),
        uiOutput(ns("groups_warning"))
      ),
      
      hr(class = "startup-divider"),
      
      # Section 3: Output folder selection
      div(
        class = "startup-section",
        div(class = "section-label", "3. Output folder"),
        div(class = "section-label-sub",
            "Folder where BiomiX will save results and the JSON configuration file"),
        div(
          style = "display: flex; align-items: center; gap: 12px;",
          shinyDirButton(
            id    = ns("output_dir"),
            label = "Select folder...",
            title = "Select the output folder",
            icon  = icon("folder-open"),
            class = "btn-browse"
          ),
          uiOutput(ns("output_status"))
        )
      )
      
    ) # end startup-card
  )
}


# -----------------------------------------------------------------------------
# SERVER
# -----------------------------------------------------------------------------

startupServer <- function(id, roots) {
  moduleServer(id, function(input, output, session) {
    
    # Internal reactive state
    metadata_df   <- reactiveVal(NULL)
    metadata_path <- reactiveVal("")
    output_dir    <- reactiveVal("")
    
    # Connect metadata file button to the filesystem
    shinyFileChoose(
      input     = input,
      id        = "metadata_file",
      roots     = roots,
      session   = session,
      filetypes = c("xlsx", "xls", "tsv", "csv", "txt")
    )
    
    # Connect output folder button to the filesystem
    shinyDirChoose(
      input   = input,
      id      = "output_dir",
      roots   = roots,
      session = session
    )
    
    # React to metadata file selection
    observeEvent(input$metadata_file, {
      info <- parseFilePaths(roots, input$metadata_file)
      if (nrow(info) == 0) return()
      
      path <- as.character(info$datapath[1])
      metadata_path(path)
      
      df <- tryCatch(
        read_metadata(path),
        error = function(e) {
          showNotification(e$message, type = "error", duration = 8)
          NULL
        }
      )
      
      if (!is.null(df)) {
        metadata_df(df)
        conditions <- get_conditions(df)
        updateSelectInput(session, "group1", choices = conditions, selected = conditions[1])
        updateSelectInput(session, "group2", choices = conditions,
                          selected = if (length(conditions) >= 2) conditions[2] else conditions[1])
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    # React to output folder selection
    observeEvent(input$output_dir, {
      path <- parseDirPath(roots, input$output_dir)
      if (length(path) > 0 && nchar(path) > 0)
        output_dir(as.character(path))
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    # Display selected metadata path
    output$metadata_path_display <- renderUI({
      p <- metadata_path()
      if (nchar(p) == 0)
        span(class = "file-label-empty", "No file selected")
      else
        span(class = "file-label-ok", icon("circle-check"), " ", basename(p))
    })
    
    # Feedback: metadata info after successful read
    output$metadata_status <- renderUI({
      df <- metadata_df()
      if (is.null(df)) return(NULL)
      n_samples    <- nrow(df)
      n_conditions <- length(unique(df$CONDITION))
      conditions   <- paste(unique(df$CONDITION), collapse = ", ")
      div(class = "status-ok",
          icon("circle-check"),
          sprintf(" %d samples · %d conditions: %s", n_samples, n_conditions, conditions))
    })
    
    # Warning if both groups are identical
    output$groups_warning <- renderUI({
      req(input$group1, input$group2)
      if (input$group1 == input$group2)
        div(class = "status-warn", icon("triangle-exclamation"),
            " GROUP_1 and GROUP_2 are identical — please select two different conditions.")
    })
    
    # Feedback: selected output folder
    output$output_status <- renderUI({
      d <- output_dir()
      if (nchar(d) == 0)
        span(class = "file-label-empty", "No folder selected")
      else
        span(class = "file-label-ok", icon("folder-open"), " ", d)
    })
    
    # Value returned to the parent (app.R)
    reactive({
      list(
        metadata_path = metadata_path(),
        metadata_df   = metadata_df(),
        group1        = input$group1 %||% "",
        group2        = input$group2 %||% "",
        output_dir    = output_dir(),
        ready         = !is.null(metadata_df()) &&
          nchar(output_dir()) > 0 &&
          (input$group1 %||% "") != (input$group2 %||% "")
      )
    })
    
  })
}