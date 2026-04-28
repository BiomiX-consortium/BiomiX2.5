# =============================================================================
# mod_adv_metadata_filter.R
# Module: "Metadata Filtering" tab of the Advanced Options panel
#
# Manages up to 4 sample filters based on metadata columns.
# Column choices are populated automatically when the metadata file is loaded
# in the Setup screen.
#
# Each filter row contains:
#   - Column:    metadata column to filter on
#   - Data type: ID / Numerical / Categorical
#   - Threshold: numeric cutoff or exact category value
#
# A green indicator and an active-filter summary are shown to help the user
# track which filters are currently configured.
#
# Returned value: adv_metadata_filter_data() -> list of 4 filter entries
# =============================================================================

advMetadataFilterUI <- function(id) {
  ns <- NS(id)

  tagList(

    div(
      style = "margin-bottom: 16px; font-size: 12px; color: var(--biomix-text-muted);",
      icon("circle-info"),
      " Column choices are populated automatically from the metadata file loaded in the
        Setup screen. Load the metadata file first if the dropdowns are empty."
    ),

    .adv_section_header("Sample filtering rules"),

    # Column header row
    fluidRow(
      style = "font-family: 'Space Mono', monospace; font-size: 11px; font-weight: 700;
               text-transform: uppercase; color: var(--biomix-text-muted);
               padding: 6px 15px; background: var(--biomix-blue-light);
               border-radius: 6px; margin-bottom: 4px;",
      column(1, ""),
      column(3, "Column"),
      column(3, "Data type"),
      column(4, "Threshold / value"),
      column(1, "")
    ),

    # 4 filter rows
    lapply(1:4, function(i) {
      div(
        style = paste0(
          "background:", if (i %% 2 == 0) "var(--biomix-blue-light)" else "white", ";",
          "border: 1px solid var(--biomix-border); border-radius: 6px;",
          "padding: 8px 12px; margin-bottom: 6px;"
        ),
        fluidRow(
          column(1,
            div(
              style = "padding-top: 6px; font-family: 'Space Mono', monospace;
                       font-size: 12px; font-weight: 700; color: var(--biomix-blue);",
              paste0("F", i)
            )
          ),
          column(3,
            selectInput(ns(paste0("filt_col_", i)), NULL,
                        choices = c("ID"), selected = "ID", width = "100%")
          ),
          column(3,
            selectInput(ns(paste0("filt_type_", i)), NULL,
                        choices  = c("ID","Numerical","Categorical"),
                        selected = "ID", width = "100%")
          ),
          column(4,
            textInput(ns(paste0("filt_threshold_", i)), NULL,
                      placeholder = "e.g. 50  or  GroupA", width = "100%")
          ),
          # Visual indicator: green when filter is configured
          column(1,
            div(style = "padding-top: 6px;",
              uiOutput(ns(paste0("filt_status_", i)))
            )
          )
        )
      )
    }),

    # Active filter summary
    div(style = "margin-top: 16px;", uiOutput(ns("filters_summary")))
  )
}


advMetadataFilterServer <- function(id, metadata_df) {
  moduleServer(id, function(input, output, session) {

    # Populate column dropdowns when metadata is loaded
    observe({
      df <- metadata_df()
      if (!is.null(df)) {
        cols <- c("ID", colnames(df))
        for (i in 1:4)
          updateSelectInput(session, paste0("filt_col_", i),
                            choices = cols, selected = "ID")
      }
    })

    # Status indicator per row (green = configured, grey = empty)
    lapply(1:4, function(i) {
      output[[paste0("filt_status_", i)]] <- renderUI({
        col    <- input[[paste0("filt_col_",       i)]] %||% "ID"
        thr    <- input[[paste0("filt_threshold_", i)]] %||% ""
        active <- col != "ID" && nchar(trimws(thr)) > 0
        if (active)
          span(style = "color: var(--biomix-green); font-size: 16px;", icon("circle-check"))
        else
          span(style = "color: var(--biomix-border); font-size: 16px;", icon("circle"))
      })
    })

    # Summary of active filters
    output$filters_summary <- renderUI({
      active <- Filter(Negate(is.null), lapply(1:4, function(i) {
        col  <- input[[paste0("filt_col_",       i)]] %||% "ID"
        type <- input[[paste0("filt_type_",      i)]] %||% "ID"
        thr  <- input[[paste0("filt_threshold_", i)]] %||% ""
        if (col != "ID" && nchar(trimws(thr)) > 0)
          sprintf("Filter %d: %s [%s] %s", i, col, type, thr)
        else NULL
      }))

      if (length(active) == 0) {
        div(class = "file-label-empty", icon("filter"), " No active filters")
      } else {
        div(
          class = "status-ok",
          icon("filter"), sprintf(" %d active filter(s):", length(active)),
          tags$ul(style = "margin: 4px 0 0 16px;",
            lapply(active, function(x) tags$li(style = "font-size:12px;", x))
          )
        )
      }
    })

    # Value returned to the parent (app.R)
    reactive({
      lapply(1:4, function(i) {
        list(
          column    = input[[paste0("filt_col_",       i)]] %||% "ID",
          data_type = input[[paste0("filt_type_",      i)]] %||% "ID",
          threshold = input[[paste0("filt_threshold_", i)]] %||% ""
        )
      })
    })

  })
}
