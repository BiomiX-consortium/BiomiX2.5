# =============================================================================
# app.R
# BiomiX Shiny application — entry point.
#
# This file:
#   1. Sources global.R (libraries + utility functions)
#   2. Sources all modules
#   3. Defines the main UI and Server
#   4. Launches the app
#
# No business logic lives here — everything is delegated to the modules.
# =============================================================================



#Loading R enrironemnt or alternatively conda one (to be defined)
#Load renv environment:
renv::load(project = "_INSTALL")

# --- Load globals and modules ------------------------------------------------

source("global.R")
library(shinyjs)
library(shinyFiles)

source("modules/mod_startup.R")
source("modules/mod_inputs.R")
source("modules/mod_integration.R")
source("modules/mod_adv_general.R")
source("modules/mod_adv_metabolomics.R")    # replaces ms1 + ms2 separate modules
source("modules/mod_adv_mofa.R")
source("modules/mod_adv_snf_nemo.R")
source("modules/mod_adv_minttea.R")
source("modules/mod_adv_metadata_filter.R")
source("modules/mod_qc_imputation.R")


# =============================================================================
# FILESYSTEM ROOTS for shinyFiles (Linux / Docker)
#
# In Docker: BIOMIX_DATA_DIR is set by docker-compose.yml (e.g. /data)
# Outside Docker (local development): falls back to the directory where
# app.R lives, so shinyFiles browses from the app folder by default.
# =============================================================================

DATA_DIR <- Sys.getenv("BIOMIX_DATA_DIR", unset = "")

# If the env variable is not set, use the app's own directory as root.
# This works both in RStudio (getwd()) and when launched via Rscript.
if (nchar(DATA_DIR) == 0) {
  DATA_DIR <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = FALSE
  )
}

BIOMIX_ROOTS <- c(Data = DATA_DIR)


# =============================================================================
# MAIN UI
# =============================================================================

ui <- fluidPage(
  
  useShinyjs(),
  
  tags$head(
    tags$style(HTML("

      /* ---- Fonts & base colours ---------------------------------------- */
      @import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@300;400;500;600&family=Space+Mono:wght@400;700&display=swap');

      :root {
        --biomix-blue:       #4a7cbf;
        --biomix-blue-light: #e8f0fc;
        --biomix-green:      #2e9e6b;
        --biomix-warn:       #c47c1a;
        --biomix-danger:     #c0392b;
        --biomix-bg:         #f4f6fb;
        --biomix-card:       #ffffff;
        --biomix-border:     #d6dde8;
        --biomix-text:       #1e2533;
        --biomix-text-muted: #6b7a99;
        --radius:            8px;
        --shadow:            0 2px 12px rgba(74,124,191,0.10);
      }

      body {
        background: var(--biomix-bg);
        color: var(--biomix-text);
        font-family: 'DM Sans', sans-serif;
        font-size: 14px;
      }

      /* ---- Tab navigation ----------------------------------------------- */
      .nav-tabs > li > a {
        font-family: 'Space Mono', monospace;
        font-size: 12px;
        letter-spacing: 0.05em;
        color: var(--biomix-text-muted);
        border-radius: var(--radius) var(--radius) 0 0;
      }
      .nav-tabs > li.active > a {
        color: var(--biomix-blue);
        border-bottom-color: var(--biomix-card);
        font-weight: 700;
      }

      /* ---- Startup header ----------------------------------------------- */
      .startup-header {
        padding: 24px 0 16px 0;
        border-bottom: 2px solid var(--biomix-blue-light);
        margin-bottom: 24px;
      }
      .startup-title {
        font-family: 'Space Mono', monospace;
        font-size: 28px;
        font-weight: 700;
        color: var(--biomix-blue);
        letter-spacing: -0.02em;
      }
      .startup-version {
        font-size: 14px;
        color: var(--biomix-text-muted);
        margin-left: 8px;
        vertical-align: middle;
      }
      .startup-subtitle {
        color: var(--biomix-text-muted);
        font-size: 13px;
        margin-top: 4px;
      }

      /* ---- Cards -------------------------------------------------------- */
      .startup-card, .integration-panel {
        background: var(--biomix-card);
        border: 1px solid var(--biomix-border);
        border-radius: var(--radius);
        padding: 24px;
        margin-bottom: 20px;
        box-shadow: var(--shadow);
      }

      /* ---- Section labels ----------------------------------------------- */
      .section-label {
        font-family: 'Space Mono', monospace;
        font-size: 12px;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.08em;
        color: var(--biomix-blue);
        margin-bottom: 4px;
      }
      .section-label-sub {
        font-size: 12px;
        color: var(--biomix-text-muted);
        margin-bottom: 10px;
      }
      .startup-divider  { border-color: var(--biomix-blue-light); margin: 20px 0; }
      .startup-section  { margin-bottom: 4px; }

      /* ---- Status messages ---------------------------------------------- */
      .status-ok   { color: var(--biomix-green); font-size: 13px; margin-top: 6px; font-weight: 500; }
      .status-warn { color: var(--biomix-warn);  font-size: 13px; margin-top: 6px; font-weight: 500; }

      /* ---- Group labels ------------------------------------------------- */
      .group-label     { font-size: 12px; font-weight: 600; margin-bottom: 4px;
                         padding: 2px 8px; border-radius: 4px; display: inline-block; }
      .condition-label { background: #dbeafe; color: #1d4ed8; }
      .control-label   { background: #dcfce7; color: #166534; }

      /* ---- Inputs table ------------------------------------------------- */
      .inputs-header {
        background: var(--biomix-blue);
        color: white;
        font-family: 'Space Mono', monospace;
        font-size: 11px;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.06em;
        padding: 10px 15px;
        border-radius: var(--radius) var(--radius) 0 0;
        margin: 0;
      }
      .col-header { text-align: center; line-height: 1.3; }

      .input-row {
        background: var(--biomix-card);
        border: 1px solid var(--biomix-border);
        border-top: none;
        padding: 8px 15px;
        align-items: center;
        display: flex;
      }
      .input-row-alt   { background: var(--biomix-blue-light); }
      .input-row-label {
        font-family: 'Space Mono', monospace;
        font-size: 12px;
        font-weight: 700;
        color: var(--biomix-blue);
      }
      .checkbox-center { text-align: center; padding-top: 4px; }

      /* ---- Upload buttons ----------------------------------------------- */
      .btn-upload {
        background: var(--biomix-blue);
        color: white;
        border: none;
        border-radius: var(--radius);
        font-size: 12px;
        padding: 5px 14px;
        width: 100%;
        cursor: pointer;
        transition: background 0.2s;
      }
      .btn-upload:hover     { background: #3566a8; color: white; }
      .btn-upload-done {
        background: var(--biomix-green) !important;
        color: white !important;
        border: none;
        border-radius: var(--radius);
        font-size: 12px;
        padding: 5px 14px;
        width: 100%;
      }
      .upload-btn-wrapper { display: flex; flex-direction: column; gap: 4px; }
      .file-label-empty   { font-size: 11px; color: var(--biomix-text-muted); }
      .file-label-ok      { font-size: 11px; color: var(--biomix-green); font-weight: 600; }
      .file-label-path    { font-size: 10px; color: var(--biomix-text-muted); word-break: break-all; }

      /* ---- Integration panel -------------------------------------------- */
      .panel-section-title {
        font-family: 'Space Mono', monospace;
        font-size: 13px;
        font-weight: 700;
        color: var(--biomix-blue);
        margin-bottom: 16px;
        text-transform: uppercase;
        letter-spacing: 0.06em;
      }
      .integration-block       { padding: 8px 0; }
      .integration-block-label {
        font-size: 11px;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.06em;
        color: var(--biomix-text-muted);
        margin-bottom: 6px;
        background: var(--biomix-blue-light);
        padding: 3px 8px;
        border-radius: 4px;
        display: inline-block;
      }
      .param-label { font-size: 11px; color: var(--biomix-text-muted); margin-bottom: 3px; }

      /* ---- Browse buttons (shinyFiles) ---------------------------------- */
      .btn-browse {
        background: var(--biomix-blue);
        color: white !important;
        border: none;
        border-radius: var(--radius);
        font-size: 13px;
        padding: 7px 18px;
        cursor: pointer;
        transition: background 0.2s;
        white-space: nowrap;
      }
      .btn-browse:hover { background: #3566a8; }

      /* ---- Advanced modules --------------------------------------------- */
      .adv-section-header {
        display: flex;
        align-items: center;
        gap: 10px;
        margin: 20px 0 10px 0;
      }
      .adv-hr-left, .adv-hr-right {
        flex: 1;
        border: none;
        border-top: 1px solid var(--biomix-border);
        margin: 0;
      }
      .adv-section-title {
        font-family: 'Space Mono', monospace;
        font-size: 11px;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.08em;
        color: var(--biomix-blue);
        white-space: nowrap;
        background: var(--biomix-blue-light);
        padding: 3px 10px;
        border-radius: 4px;
      }
      .adv-param-block  { margin-bottom: 8px; }
      .adv-param-label  {
        font-size: 11px;
        font-weight: 600;
        color: var(--biomix-text-muted);
        margin-bottom: 3px;
        background: var(--biomix-blue-light);
        padding: 2px 7px;
        border-radius: 3px;
        display: inline-block;
      }
      .btn-browse-sm {
        background: var(--biomix-blue);
        color: white !important;
        border: none;
        border-radius: var(--radius);
        font-size: 12px;
        padding: 5px 12px;
        cursor: pointer;
        white-space: nowrap;
      }
      .btn-browse-sm:hover { background: #3566a8; }

      /* ---- Run panel ---------------------------------------------------- */
      .run-panel {
        background: var(--biomix-card);
        border: 1px solid var(--biomix-border);
        border-radius: var(--radius);
        padding: 20px 24px;
        margin-top: 20px;
        display: flex;
        align-items: center;
        gap: 20px;
        box-shadow: var(--shadow);
      }
      .btn-run {
        background: var(--biomix-blue);
        color: white !important;
        font-family: 'Space Mono', monospace;
        font-weight: 700;
        font-size: 13px;
        letter-spacing: 0.05em;
        border: none;
        border-radius: var(--radius);
        padding: 12px 32px;
        cursor: pointer;
        transition: all 0.2s;
        box-shadow: 0 4px 12px rgba(74,124,191,0.3);
      }
      .btn-run:hover {
        background: #3566a8;
        transform: translateY(-1px);
        box-shadow: 0 6px 16px rgba(74,124,191,0.4);
      }
      .run-status { font-size: 13px; color: var(--biomix-text-muted); }

    "))
  ),
  
  # Main container
  div(
    style = "max-width: 1300px; margin: 0 auto; padding: 0 20px 40px 20px;",
    
    tabsetPanel(
      id   = "main_tabs",
      type = "tabs",
      
      # -----------------------------------------------------------------------
      # TAB 1: Setup
      # -----------------------------------------------------------------------
      tabPanel(
        title = "⚙ Setup",
        value = "tab_setup",
        br(),
        
        # Startup module: metadata, groups, output folder
        startupUI("startup"),
        
        div(class = "section-label", style = "margin: 24px 0 8px 0;", "Omics Inputs"),
        
        # Inputs module: 6-row table
        inputsUI("inputs"),
        
        br(),
        
        # Integration module
        integrationUI("integration"),
        
        br(),
        
        # Run panel
        div(
          class = "run-panel",
          actionButton("btn_adv", "Advanced Options",
                       class = "btn-default", icon = icon("sliders")),
          actionButton("btn_run", "Run BiomiX",
                       class = "btn-run", icon = icon("play")),
          downloadButton("btn_download_json", "Download JSON",
                         class = "btn-default", icon = icon("download")),
          uiOutput("run_status")
        )
      ),
      
      # -----------------------------------------------------------------------
      # TAB 2: Advanced Options — inner tabset with all advanced modules
      # -----------------------------------------------------------------------
      tabPanel(
        title = "Advanced Options",
        value = "tab_advanced",
        br(),
        tabsetPanel(
          type = "tabs",
          tabPanel("General",            br(), advGeneralUI("adv_general")),
          tabPanel("Metabolomics",      br(), advMetabolomicsUI("adv_met")),
          tabPanel("MOFA & Diablo",      br(), advMofaUI("adv_mofa")),
          tabPanel("SNF / NEMO",         br(), advSnfNemoUI("adv_snf_nemo")),
          tabPanel("MintTea",            br(), advMintteaUI("adv_minttea")),
          tabPanel("Metadata Filtering", br(), advMetadataFilterUI("adv_meta_filter"))
        )
      ),
      
      # -----------------------------------------------------------------------
      # TAB 4: QC & Imputation
      # -----------------------------------------------------------------------
      tabPanel(
        title = "🔬 QC & Imputation",
        value = "tab_qc",
        br(),
        qcImputationUI("qc_imputation")
      ),
      tabPanel(
        title = "JSON Preview",
        value = "tab_json",
        br(),
        div(class = "startup-card",
            div(class = "section-label", "COMBINED_COMMANDS.json — live preview"),
            br(),
            verbatimTextOutput("json_preview")
        )
      )
    )
  )
)


# =============================================================================
# MAIN SERVER
# =============================================================================

server <- function(input, output, session) {
  
  # ---------------------------------------------------------------------------
  # Initialise all modules
  # Each module returns a reactive with its collected data.
  # ---------------------------------------------------------------------------
  
  startup_data         <- startupServer("startup",         roots = BIOMIX_ROOTS)
  inputs_data          <- inputsServer("inputs",           roots = BIOMIX_ROOTS)
  integration_data     <- integrationServer("integration")
  adv_general_data     <- advGeneralServer("adv_general",  roots = BIOMIX_ROOTS)
  adv_met_data         <- advMetabolomicsServer("adv_met", roots = BIOMIX_ROOTS)
  adv_mofa_data        <- advMofaServer("adv_mofa",        roots = BIOMIX_ROOTS)
  adv_snf_data         <- advSnfNemoServer("adv_snf_nemo",
                                           metadata_df = reactive(startup_data()$metadata_df))
  adv_minttea_data     <- advMintteaServer("adv_minttea")
  adv_meta_filter_data <- advMetadataFilterServer("adv_meta_filter",
                                                  metadata_df = reactive(startup_data()$metadata_df))
  
  # QC & Imputation module — fully independent, no inputs needed from other modules
  qcImputationServer("qc_imputation")
  
  # ---------------------------------------------------------------------------
  # "Advanced Options" button -> switch to the advanced tab
  # ---------------------------------------------------------------------------
  observeEvent(input$btn_adv, {
    updateTabsetPanel(session, "main_tabs", selected = "tab_advanced")
  })
  
  # ---------------------------------------------------------------------------
  # Central reactive: assemble the complete JSON structure.
  # This is the single point where all module outputs are combined.
  # ---------------------------------------------------------------------------
  combined_json_data <- reactive({
    sd  <- startup_data()
    id  <- inputs_data()
    ig  <- integration_data()
    ag  <- adv_general_data()
    met <- adv_met_data()
    mf  <- adv_mofa_data()
    sn  <- adv_snf_data()
    mt  <- adv_minttea_data()
    fil <- adv_meta_filter_data()
    
    adv <- list(
      # Row 1: user-configured values
      list(
        ADVANCED_OPTION_TRASCRIPTOMICS                      = as.character(ag$tx_fc),
        ADVANCED_OPTION_SUBGROUPING                         = as.character(ag$sub_zscore2),
        ADVANCED_OPTION_METABOLOMICS                        = ag$met_fc,
        ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1         = met$ms1_ion_mode,
        ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2       = met$ms1_pos,
        ADVANCED_OPTION_METHYLOMICS                         = as.character(ag$meth_fc),
        ADVANCED_OPTION_MOFA                                = mf$mofa_max_iter,
        ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2         = as.character(met$ms2_ppm),
        ADVANCED_OPTION_METABOLOMICS_ANNOTATION_GENERAL     = met$annot_mode,
        ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY    = mf$biblio_type,
        ADVANCED_OPTION_MOFA_INTERPRETATION_PATHWAY         = as.character(mf$pathway_padj),
        ADVANCED_OPTION_MOFA_INTERPRETATION_CLINICAL        = mf$clinical_numerical,
        ADVANCED_OPTION_METADATA_FILTERING_1                = fil[[1]]$column,
        ADVANCED_OPTION_METADATA_FILTERING_2                = fil[[2]]$column,
        ADVANCED_OPTION_METADATA_FILTERING_3                = fil[[3]]$column,
        ADVANCED_OPTION_METADATA_FILTERING_4                = fil[[4]]$column,
        ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_INDEX = met$ms1_input_idx,
        ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES       = met$ms1_annot_file,
        ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY          = met$ms2_dir,
        ADVANCED_OPTION_CLINIC_DATA_DIRECTORY               = mf$clinic_file,
        ADVANCED_OPTION_CLUSTERING_OPTIONS                  = ag$clust_distance,
        ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY_2  = as.character(mf$biblio_top_contrib),
        ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2       = met$ms2_pos,
        ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3       = met$ms2_ion_mode,
        ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3_INDEX = met$ms2_annot_idx,
        ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2   = met$ms2_annot_file,
        ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_4       = met$ms2_ms1_db,
        ADVANCED_OPTION_SNF_OPTIONS                         = sn$snf_neighbors,
        ADVANCED_OPTION_NEMO_OPTIONS                        = as.character(sn$nemo_neighbors),
        ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES          = sn$snf_cat_cols,
        ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS            = sn$snf_var_heatmap,
        ADVANCED_OPTION_FILE_PATH_DIABLO_DESIGN             = mf$diablo_design_file,
        ADVANCED_OPTION_DIABLO_OPTIONS                      = mf$diablo_feature_sel,
        ADVANCED_OPTION_MINTTEA_OPTIONS                     = as.character(mt$n_repeats),
        ADVANCED_OPTION_MINTTEA_FEATURES                    = as.character(mt$keep_x)
      ),
      # Row 2: minimum/alternative thresholds (fixed constants from global.R)
      COMMANDS_ADVANCED_ROW2,
      # Row 3: maximum values / extra options (fixed constants from global.R)
      COMMANDS_ADVANCED_ROW3
    )
    
    list(
      COMMANDS          = build_commands(id),
      COMMANDS_MOFA     = build_commands_mofa(ig),
      COMMANDS_ADVANCED = adv,
      COMMAND_LINE_ARGS = list(
        GROUP_1  = sd$group1,
        GROUP_2  = sd$group2,
        BASE_DIR = sd$output_dir
      ),
      DIRECTORY_INFO = list(
        METADATA_DIR = sd$metadata_path,
        OUTPUT_DIR   = sd$output_dir
      )
    )
  })
  
  # ---------------------------------------------------------------------------
  # JSON Preview tab — live formatted output
  # ---------------------------------------------------------------------------
  output$json_preview <- renderText({
    jsonlite::toJSON(combined_json_data(), pretty = TRUE,
                     auto_unbox = TRUE, na = "null")
  })
  
  # ---------------------------------------------------------------------------
  # Download button: export COMBINED_COMMANDS.json directly from the browser
  # ---------------------------------------------------------------------------
  output$btn_download_json <- downloadHandler(
    filename = function() "COMBINED_COMMANDS.json",
    content  = function(file) {
      json_string <- jsonlite::toJSON(
        combined_json_data(),
        pretty      = TRUE,
        auto_unbox  = TRUE,
        na          = "null"
      )
      writeLines(json_string, file)
    },
    contentType = "application/json"
  )
  
  # ---------------------------------------------------------------------------
  # Run button: validate → write JSON to BASE_DIR → launch BiomiX_BETA.r
  #
  # BiomiX_BETA.r expects:
  #   args[1] = GROUP_1   (condition)
  #   args[2] = GROUP_2   (control)
  #   args[3] = BASE_DIR  (app directory — where BiomiX_BETA.r lives and
  #                         where COMBINED_COMMANDS.json will be written)
  #
  # The script does setwd(BASE_DIR) and reads COMBINED_COMMANDS.json from there.
  # ---------------------------------------------------------------------------
  observeEvent(input$btn_run, {
    
    sd <- startup_data()
    
    # Validate prerequisites
    if (!sd$ready) {
      showNotification(
        "Please complete the Setup first: load metadata, select two different groups, and choose an output folder.",
        type = "warning", duration = 6
      )
      return()
    }
    
    # BASE_DIR = directory where app.R and BiomiX_BETA.r live
    # JSON must be written here so BiomiX_BETA.r can find it with setwd()
    base_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    json_path <- file.path(base_dir, "COMBINED_COMMANDS.json")
    
    # Step 1: write the JSON
    tryCatch({
      json_data <- combined_json_data()
      write_combined_json(
        commands          = json_data$COMMANDS,
        commands_mofa     = json_data$COMMANDS_MOFA,
        commands_advanced = json_data$COMMANDS_ADVANCED,
        group1            = sd$group1,
        group2            = sd$group2,
        base_dir          = base_dir,
        metadata_dir      = sd$metadata_path,
        output_dir        = sd$output_dir,
        out_path          = json_path
      )
      showNotification(
        paste0("JSON written to: ", json_path),
        type = "message", duration = 4
      )
    }, error = function(e) {
      showNotification(
        paste0("Error writing JSON: ", e$message),
        type = "error", duration = 8
      )
      return()
    })
    
    # Step 2: locate BiomiX_BETA.r (same directory as app.R)
    pipeline_script <- file.path(base_dir, "BiomiX_BETA.r")
    
    if (!file.exists(pipeline_script)) {
      showNotification(
        paste0("BiomiX_BETA.r not found in: ", base_dir,
               " — JSON was saved but pipeline was not launched."),
        type = "warning", duration = 8
      )
      return()
    }
    
    # Step 3: launch BiomiX_BETA.r via Rscript as a background process.
    # system2(..., wait = FALSE) returns immediately — the Shiny app stays
    # responsive while the pipeline runs in the background.
    #
    # Arguments passed to BiomiX_BETA.r:
    #   [1] GROUP_1   [2] GROUP_2   [3] BASE_DIR
    rscript <- if (.Platform$OS.type == "windows") {
      file.path(R.home("bin"), "Rscript.exe")
    } else {
      file.path(R.home("bin"), "Rscript")
    }
    
    system2(
      command = rscript,
      args    = c(
        shQuote(pipeline_script),
        shQuote(sd$group1),
        shQuote(sd$group2),
        shQuote(base_dir)
      ),
      wait   = TRUE,    # blocking: output streams to the R console in real time
      stdout = "",      # "" = inherit the parent R console (not captured)
      stderr = ""       # "" = same for errors
    )
    
    showNotification(
      paste0("BiomiX pipeline started! Groups: ", sd$group1, " vs ", sd$group2,
             ". Results will be saved to: ", sd$output_dir),
      type     = "message",
      duration = 10
    )
  })
  
  # ---------------------------------------------------------------------------
  # Run status indicator (shown next to the Run button)
  # ---------------------------------------------------------------------------
  output$run_status <- renderUI({
    sd <- startup_data()
    if (!sd$ready)
      div(class = "run-status", icon("circle-exclamation"), " Setup incomplete")
    else
      div(class = "status-ok",  icon("circle-check"),       " Ready to run")
  })
  
}


# =============================================================================
# LAUNCH
# =============================================================================

shinyApp(ui = ui, server = server)

