# =============================================================================
# mod_adv_metabolomics.R
# Module: "Metabolomics" tab of the Advanced Options panel
#
# Replaces the previous mod_adv_metabolomics_ms1.R and mod_adv_metabolomics_ms2.R
# which were split incorrectly. All annotation modes now live in one place.
#
# Three annotation modes (radio buttons at the top):
#
#   ANNOTATED   — user provides a pre-annotated file with compound_name/HMDB/KEGG
#                 columns. No database matching needed.
#
#   MS1         — non-targeted m/z matching against databases.
#                 Requires: ion mode, adducts, ppm tolerance, databases,
#                 chromatography column, optional annotation file.
#
#   MS1+MS2     — m/z matching + spectral matching against MS2 databases.
#                 Everything from MS1 plus: .mzML directory, MS2-specific
#                 adducts, MS2 databases with priority ranking, MS2 ppm tolerance,
#                 optional MS2 annotation file.
#
# Sections appear/disappear via conditionalPanel based on the selected mode.
# MS1 adducts and MS2 adducts are fully separate widgets.
#
# Returned value: adv_metabolomics_data() -> list with all parameters
# =============================================================================

library(shinyFiles)

advMetabolomicsUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    
    # =========================================================================
    # MODE SELECTOR — drives all conditional panels below
    # =========================================================================
    div(
      style = "background: var(--biomix-blue-light); border-radius: 8px;
               padding: 16px 20px; margin-bottom: 20px;
               border-left: 4px solid var(--biomix-blue);",
      div(class = "adv-param-label", style = "margin-bottom: 10px;",
          "Annotation mode"),
      radioButtons(
        ns("annot_mode"), NULL,
        choices = c(
          "Annotated (pre-annotated file with compound_name / HMDB / KEGG)" = "Annotated",
          "MS1 only (non-targeted m/z matching against databases)"           = "MS1",
          "MS1 + MS2 (m/z matching + spectral matching)"                    = "MS1+MS2"
        ),
        selected = "MS1",
        inline   = FALSE
      )
    ),
    
    # =========================================================================
    # ANNOTATED MODE
    # =========================================================================
    conditionalPanel(
      condition = paste0("input['", ns("annot_mode"), "'] == 'Annotated'"),
      
      .adv_section_header("Metabolite annotation"),
      
      fluidRow(
        column(4,
               .adv_param_block("Metabolite annotation column",
                                selectInput(ns("annot_compound_col"), NULL,
                                            choices  = c("compound_name", "HMDB", "KEGG"),
                                            selected = "HMDB", width = "100%")
               )
        )
      )
    ),
    
    # =========================================================================
    # MS1 MODE (visible ONLY when MS1 is selected, not MS1+MS2)
    # =========================================================================
    conditionalPanel(
      condition = paste0("input['", ns("annot_mode"), "'] == 'MS1'"),
      
      .adv_section_header("MS1 matching parameters"),
      
      fluidRow(
        column(3,
               .adv_param_block("Ion mode",
                                selectInput(ns("ms1_ion_mode"), NULL,
                                            choices  = c("positive", "negative", "neutral"),
                                            selected = "positive", width = "100%")
               )
        ),
        column(3,
               .adv_param_block("M/Z tolerance (ppm)",
                                numericInput(ns("ms1_ppm"), NULL,
                                             value = 15, min = 1, max = 100, step = 1, width = "100%")
               )
        )
      ),
      
      # MS1 adducts positive
      .adv_section_header("MS1 adducts — positive mode"),
      div(class = "adv-param-block",
          fluidRow(
            column(2, checkboxInput(ns("ms1_pos_MH"),    "M+H",     value = TRUE)),
            column(2, checkboxInput(ns("ms1_pos_M2H"),   "M+2H",    value = TRUE)),
            column(2, checkboxInput(ns("ms1_pos_MNa"),   "M+Na",    value = TRUE)),
            column(2, checkboxInput(ns("ms1_pos_MK"),    "M+K",     value = TRUE)),
            column(2, checkboxInput(ns("ms1_pos_MNH4"),  "M+NH4",   value = TRUE)),
            column(2, checkboxInput(ns("ms1_pos_MHH2O"), "M+H-H2O", value = TRUE))
          )
      ),
      
      # MS1 adducts negative
      .adv_section_header("MS1 adducts — negative mode"),
      div(class = "adv-param-block",
          fluidRow(
            column(2, checkboxInput(ns("ms1_neg_MH"),    "M-H",     value = FALSE)),
            column(2, checkboxInput(ns("ms1_neg_MCl"),   "M-Cl",    value = FALSE)),
            column(2, checkboxInput(ns("ms1_neg_MFAH"),  "M+FA-H",  value = FALSE)),
            column(2, checkboxInput(ns("ms1_neg_MHH2O"), "M-H-H2O", value = FALSE))
          )
      ),
      
      # MS1 databases
      .adv_section_header("Databases MS1"),
      div(class = "adv-param-block",
          fluidRow(
            column(2, checkboxInput(ns("ms1_db_hmdb"),      "HMDB",      value = TRUE)),
            column(2, checkboxInput(ns("ms1_db_lipidmaps"), "LipidMaps", value = TRUE)),
            column(2, checkboxInput(ns("ms1_db_metlin"),    "Metlin",    value = TRUE)),
            column(2, checkboxInput(ns("ms1_db_kegg"),      "KEGG",      value = TRUE))
          )
      ),
      
      # MS1 optional annotation file
      .adv_section_header("MS1 annotation file (optional)"),
      fluidRow(
        column(6,
               .adv_param_block("Annotation file",
                                div(style = "display:flex; gap:10px; align-items:center;",
                                    shinyFilesButton(ns("ms1_annot_file_btn"), "Browse...",
                                                     "Select MS1 annotation file", multiple = FALSE,
                                                     class = "btn-browse-sm"),
                                    uiOutput(ns("ms1_annot_file_label"))
                                )
               )
        ),
        column(3,
               .adv_param_block("Input index",
                                selectInput(ns("ms1_input_idx"), NULL,
                                            choices  = c("No", paste0(1:6, "\u00b0 Input")),
                                            selected = "No", width = "100%")
               )
        )
      )
    ),
    
    # =========================================================================
    # MS2 EXTRA SECTION (only visible in MS1+MS2 mode)
    # =========================================================================
    conditionalPanel(
      condition = paste0("input['", ns("annot_mode"), "'] == 'MS1+MS2'"),
      
      div(
        style = "margin-top: 24px; padding: 16px;
                 background: var(--biomix-blue-light);
                 border-radius: 8px;
                 border-left: 4px solid var(--biomix-blue);",
        
        div(
          style = "font-family: 'Space Mono', monospace; font-size: 11px;
                   font-weight: 700; text-transform: uppercase;
                   letter-spacing: 0.08em; color: var(--biomix-blue);
                   margin-bottom: 12px;",
          icon("circle-plus"), " MS2 — additional parameters"
        ),
        
        # MS2 directory + ion mode
        .adv_section_header("MS2 raw files directory (.mzML)"),
        fluidRow(
          column(3,
                 .adv_param_block("M/Z tolerance MS2 (ppm)",
                                  numericInput(ns("ms2_ppm"), NULL,
                                               value = 25, min = 1, max = 100, step = 1, width = "100%")
                 )
          ),
          column(3,
                 .adv_param_block("RT match MS1-MS2 (sec)",
                                  numericInput(ns("ms2_rt"), NULL,
                                               value = 10, min = 1, max = 100, step = 1, width = "100%")
                 )
          ),
          column(3,
                 .adv_param_block("Chromatography column",
                                  selectInput(ns("ms2_col_type"), NULL,
                                              choices  = c("rp", "hilic"),
                                              selected = "rp", width = "100%")
                 )
          ),
          column(3,
                 .adv_param_block("Ion mode MS2",
                                  selectInput(ns("ms2_ion_mode"), NULL,
                                              choices  = c("positive", "negative"),
                                              selected = "positive", width = "100%")
                 )
          )
        ),
        
        # MS2 adducts positive
        .adv_section_header("MS2 adducts — positive mode"),
        div(class = "adv-param-block",
            fluidRow(
              column(2, checkboxInput(ns("ms2_pos_MH"),    "M+H",     value = TRUE)),
              column(2, checkboxInput(ns("ms2_pos_M2H"),   "M+2H",    value = TRUE)),
              column(2, checkboxInput(ns("ms2_pos_MNa"),   "M+Na",    value = TRUE)),
              column(2, checkboxInput(ns("ms2_pos_MK"),    "M+K",     value = TRUE)),
              column(2, checkboxInput(ns("ms2_pos_MNH4"),  "M+NH4",   value = TRUE)),
              column(2, checkboxInput(ns("ms2_pos_MHH2O"), "M+H-H2O", value = TRUE))
            )
        ),
        
        # MS2 adducts negative
        .adv_section_header("MS2 adducts — negative mode"),
        div(class = "adv-param-block",
            fluidRow(
              column(2, checkboxInput(ns("ms2_neg_MH"),    "M-H",     value = FALSE)),
              column(2, checkboxInput(ns("ms2_neg_MCl"),   "M-Cl",    value = FALSE)),
              column(2, checkboxInput(ns("ms2_neg_MFAH"),  "M+FA-H",  value = FALSE)),
              column(2, checkboxInput(ns("ms2_neg_MHH2O"), "M-H-H2O", value = FALSE))
            )
        ),
        
        # MS2 databases with priority
        .adv_section_header("Databases MS2 (assign search priority)"),
        fluidRow(
          column(3,
                 .adv_param_block("HMDB",
                                  selectInput(ns("ms2_db_hmdb"), NULL,
                                              choices  = c("No","1 (priority)","2 (priority)","3 (priority)"),
                                              selected = "1 (priority)", width = "100%")
                 )
          ),
          column(3,
                 .adv_param_block("MONA",
                                  selectInput(ns("ms2_db_mona"), NULL,
                                              choices  = c("No","1 (priority)","2 (priority)","3 (priority)"),
                                              selected = "2 (priority)", width = "100%")
                 )
          ),
          column(3,
                 .adv_param_block("MassBank",
                                  selectInput(ns("ms2_db_massbank"), NULL,
                                              choices  = c("No","1 (priority)","2 (priority)","3 (priority)"),
                                              selected = "3 (priority)", width = "100%")
                 )
          )
        ),
        
        # MS1 databases used in MS2 workflow
        .adv_section_header("MS1 databases used in MS2 workflow"),
        div(class = "adv-param-block",
            fluidRow(
              column(2, checkboxInput(ns("ms2_ms1_hmdb"),      "HMDB",      value = TRUE)),
              column(2, checkboxInput(ns("ms2_ms1_lipidmaps"), "LipidMaps", value = TRUE)),
              column(2, checkboxInput(ns("ms2_ms1_metlin"),    "Metlin",    value = TRUE)),
              column(2, checkboxInput(ns("ms2_ms1_kegg"),      "KEGG",      value = TRUE))
            )
        ),
        
        # MS1 tolerance for MS2 + annotation file
        .adv_section_header("MS1 matching in MS2 workflow"),
        fluidRow(
          column(3,
                 .adv_param_block("MS1 tolerance (ppm) for MS2",
                                  numericInput(ns("ms2_ms1_ppm"), NULL,
                                               value = 15, min = 1, max = 100, step = 1, width = "100%")
                 )
          ),
          column(3,
                 .adv_param_block("MS2 annotation file index",
                                  selectInput(ns("ms2_annot_idx"), NULL,
                                              choices  = c("No", paste0(1:6, "\u00b0 Input")),
                                              selected = "No", width = "100%")
                 )
          ),
          column(3,
                 .adv_param_block("MS2 annotation file (optional)",
                                  div(style = "display:flex; gap:10px; align-items:center;",
                                      shinyFilesButton(ns("ms2_annot_file_btn"), "Browse...",
                                                       "Select MS2 annotation file", multiple = FALSE,
                                                       class = "btn-browse-sm"),
                                      uiOutput(ns("ms2_annot_file_label"))
                                  )
                 )
          ),
          column(3,
                 .adv_param_block("Folder containing .mzML files",
                                  div(style = "display:flex; gap:10px; align-items:center;",
                                      shinyDirButton(ns("ms2_dir_btn"), "Browse folder...",
                                                     "Select .mzML folder", class = "btn-browse-sm"),
                                      uiOutput(ns("ms2_dir_label"))
                                  )
                 )
          )
        )
      ) # end MS2 box
    )   # end MS2 conditionalPanel
  )
}


# =============================================================================
# SERVER
# =============================================================================

advMetabolomicsServer <- function(id, roots) {
  moduleServer(id, function(input, output, session) {
    
    # File path reactive values
    ms1_annot_file     <- reactiveVal(EMPTY_VALUE)
    ms2_dir            <- reactiveVal(EMPTY_VALUE)
    ms2_annot_file     <- reactiveVal(EMPTY_VALUE)
    
    # --- File browsers -------------------------------------------------------
    
    shinyFileChoose(input, "ms1_annot_file_btn", roots = roots, session = session,
                    filetypes = c("tsv","csv","txt","xlsx","xls"))
    observeEvent(input$ms1_annot_file_btn, {
      info <- parseFilePaths(roots, input$ms1_annot_file_btn)
      if (nrow(info) > 0) ms1_annot_file(as.character(info$datapath[1]))
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    output$ms1_annot_file_label <- renderUI({
      p <- ms1_annot_file()
      if (p == EMPTY_VALUE) span(class="file-label-empty","No file (optional)")
      else span(class="file-label-ok", icon("circle-check"), " ", basename(p))
    })
    
    shinyDirChoose(input, "ms2_dir_btn", roots = roots, session = session)
    observeEvent(input$ms2_dir_btn, {
      path <- parseDirPath(roots, input$ms2_dir_btn)
      if (length(path) > 0 && nchar(path) > 0) ms2_dir(as.character(path))
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    output$ms2_dir_label <- renderUI({
      p <- ms2_dir()
      if (p == EMPTY_VALUE) span(class="file-label-empty","No folder selected")
      else span(class="file-label-ok", icon("circle-check"), " ", basename(p))
    })
    
    shinyFileChoose(input, "ms2_annot_file_btn", roots = roots, session = session,
                    filetypes = c("tsv","csv","txt","xlsx","xls"))
    observeEvent(input$ms2_annot_file_btn, {
      info <- parseFilePaths(roots, input$ms2_annot_file_btn)
      if (nrow(info) > 0) ms2_annot_file(as.character(info$datapath[1]))
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    output$ms2_annot_file_label <- renderUI({
      p <- ms2_annot_file()
      if (p == EMPTY_VALUE) span(class="file-label-empty","No file (optional)")
      else span(class="file-label-ok", icon("circle-check"), " ", basename(p))
    })
    
    # --- Auto-select adducts when MS1 ion mode changes ----------------------
    observeEvent(input$ms1_ion_mode, {
      pos <- input$ms1_ion_mode == "positive"
      updateCheckboxInput(session, "ms1_pos_MH",    value = pos)
      updateCheckboxInput(session, "ms1_pos_M2H",   value = pos)
      updateCheckboxInput(session, "ms1_pos_MNa",   value = pos)
      updateCheckboxInput(session, "ms1_pos_MK",    value = pos)
      updateCheckboxInput(session, "ms1_pos_MNH4",  value = pos)
      updateCheckboxInput(session, "ms1_pos_MHH2O", value = pos)
      updateCheckboxInput(session, "ms1_neg_MH",    value = !pos)
      updateCheckboxInput(session, "ms1_neg_MCl",   value = !pos)
      updateCheckboxInput(session, "ms1_neg_MFAH",  value = !pos)
      updateCheckboxInput(session, "ms1_neg_MHH2O", value = !pos)
    }, ignoreInit = TRUE)
    
    # --- Auto-select adducts when MS2 ion mode changes ----------------------
    observeEvent(input$ms2_ion_mode, {
      pos <- input$ms2_ion_mode == "positive"
      updateCheckboxInput(session, "ms2_pos_MH",    value = pos)
      updateCheckboxInput(session, "ms2_pos_M2H",   value = pos)
      updateCheckboxInput(session, "ms2_pos_MNa",   value = pos)
      updateCheckboxInput(session, "ms2_pos_MK",    value = pos)
      updateCheckboxInput(session, "ms2_pos_MNH4",  value = pos)
      updateCheckboxInput(session, "ms2_pos_MHH2O", value = pos)
      updateCheckboxInput(session, "ms2_neg_MH",    value = !pos)
      updateCheckboxInput(session, "ms2_neg_MCl",   value = !pos)
      updateCheckboxInput(session, "ms2_neg_MFAH",  value = !pos)
      updateCheckboxInput(session, "ms2_neg_MHH2O", value = !pos)
    }, ignoreInit = TRUE)
    
    # --- Helper: build adduct string ----------------------------------------
    .adducts <- function(...) {
      labels <- c(...)
      result <- paste(names(labels)[as.logical(labels)], collapse = "/")
      if (nchar(result) == 0) EMPTY_VALUE else result
    }
    
    # --- Return value --------------------------------------------------------
    reactive({
      mode <- input$annot_mode %||% "MS1"
      
      # MS1 adducts
      ms1_pos <- .adducts(
        "M+H"     = isTRUE(input$ms1_pos_MH),
        "M+2H"    = isTRUE(input$ms1_pos_M2H),
        "M+Na"    = isTRUE(input$ms1_pos_MNa),
        "M+K"     = isTRUE(input$ms1_pos_MK),
        "M+NH4"   = isTRUE(input$ms1_pos_MNH4),
        "M+H-H2O" = isTRUE(input$ms1_pos_MHH2O)
      )
      ms1_neg <- .adducts(
        "M-H"     = isTRUE(input$ms1_neg_MH),
        "M-Cl"    = isTRUE(input$ms1_neg_MCl),
        "M+FA-H"  = isTRUE(input$ms1_neg_MFAH),
        "M-H-H2O" = isTRUE(input$ms1_neg_MHH2O)
      )
      ms1_db <- .adducts(
        "hmdb"      = isTRUE(input$ms1_db_hmdb),
        "lipidmaps" = isTRUE(input$ms1_db_lipidmaps),
        "metlin"    = isTRUE(input$ms1_db_metlin),
        "kegg"      = isTRUE(input$ms1_db_kegg)
      )
      
      # MS2 adducts
      ms2_pos <- .adducts(
        "M+H"     = isTRUE(input$ms2_pos_MH),
        "M+2H"    = isTRUE(input$ms2_pos_M2H),
        "M+Na"    = isTRUE(input$ms2_pos_MNa),
        "M+K"     = isTRUE(input$ms2_pos_MK),
        "M+NH4"   = isTRUE(input$ms2_pos_MNH4),
        "M+H-H2O" = isTRUE(input$ms2_pos_MHH2O)
      )
      ms2_neg <- .adducts(
        "M-H"     = isTRUE(input$ms2_neg_MH),
        "M-Cl"    = isTRUE(input$ms2_neg_MCl),
        "M+FA-H"  = isTRUE(input$ms2_neg_MFAH),
        "M-H-H2O" = isTRUE(input$ms2_neg_MHH2O)
      )
      ms2_ms1_db <- .adducts(
        "hmdb"      = isTRUE(input$ms2_ms1_hmdb),
        "lipidmaps" = isTRUE(input$ms2_ms1_lipidmaps),
        "metlin"    = isTRUE(input$ms2_ms1_metlin),
        "kegg"      = isTRUE(input$ms2_ms1_kegg)
      )
      
      
      list(
        annot_mode = mode,
        
        # Annotated mode — only the annotation column name is needed
        annot_compound_col = input$annot_compound_col %||% "HMDB",
        
        # MS1
        ms1_ion_mode  = input$ms1_ion_mode  %||% "positive",
        ms1_ppm       = input$ms1_ppm       %||% 15L,
        ms1_input_idx = input$ms1_input_idx %||% "No",
        ms1_pos       = ms1_pos,
        ms1_neg       = ms1_neg,
        ms1_db        = ms1_db,
        ms1_annot_file = ms1_annot_file(),
        
        # MS2
        ms2_dir        = ms2_dir(),
        ms2_ion_mode   = input$ms2_ion_mode   %||% "positive",
        ms2_ppm        = input$ms2_ppm        %||% 25L,
        ms2_rt         = input$ms2_rt        %||% 10L,
        ms2_pos        = ms2_pos,
        ms2_neg        = ms2_neg,
        ms2_col_type   = input$ms2_col_type  %||% "rp", 
        ms2_db_hmdb    = input$ms2_db_hmdb    %||% "1 (priority)",
        ms2_db_mona    = input$ms2_db_mona    %||% "2 (priority)",
        ms2_db_massbank= input$ms2_db_massbank%||% "3 (priority)",
        ms2_ms1_db_priority = paste(
          input$ms2_db_hmdb     %||% "1 (priority)",
          input$ms2_db_mona     %||% "2 (priority)",
          input$ms2_db_massbank %||% "3 (priority)",
          sep = "/"
        ),
        ms2_ms1_db     = ms2_ms1_db,
        ms2_ms1_ppm    = input$ms2_ms1_ppm    %||% 15L,
        ms2_annot_idx  = input$ms2_annot_idx  %||% "No",
        ms2_annot_file = ms2_annot_file()
      )
    })
    
  })
}