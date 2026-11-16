# =============================================================================
# mod_adv_pramigo.R
# Advanced options for the Pramigo multi-omics integration method.
#
# Pramigo parameters are organised in 4 sections (4 params each):
#   1. Training       — epoch, n_batch, lr, resume
#   2. Architecture   — n_hid, ae_n_hid, n_heads, n_layers
#   3. Graph          — topology, omic1_corr_cutoff, omic2_corr_cutoff, sample_mode
#   4. Run / Output   — leiden_resolution, leiden_min_cluster_size, cuda, checkpoint_every
#
# These 16 values are written into 4 JSON fields × 3 rows in COMMANDS_ADVANCED,
# matching the existing convention used by all other method-specific params.
# =============================================================================


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------

advPramigoUI <- function(id) {
  ns <- NS(id)

  tagList(

    # ── Section 1: Training ─────────────────────────────────────────────────
    .adv_section_header("Training"),

    fluidRow(
      column(3,
        .adv_param_block(
          "Epochs (--epoch)",
          numericInput(ns("epoch"), label = NULL, value = 100, min = 1, step = 1)
        )
      ),
      column(3,
        .adv_param_block(
          "Batch size (--n_batch)",
          numericInput(ns("n_batch"), label = NULL, value = 25, min = 1, step = 1)
        )
      ),
      column(3,
        .adv_param_block(
          "Learning rate (--lr)",
          numericInput(ns("lr"), label = NULL, value = 0.0001, min = 1e-6, step = 1e-5)
        )
      ),
      column(3,
        .adv_param_block(
          "Resume training (--resume)",
          selectInput(ns("resume"), label = NULL,
                      choices = c("off", "on"), selected = "off")
        )
      )
    ),

    # ── Section 2: Architecture ─────────────────────────────────────────────
    .adv_section_header("Architecture"),

    fluidRow(
      column(3,
        .adv_param_block(
          "Hidden units (--n_hid)",
          numericInput(ns("n_hid"), label = NULL, value = 128, min = 1, step = 1)
        )
      ),
      column(3,
        .adv_param_block(
          "AE hidden units (--ae_n_hid)",
          numericInput(ns("ae_n_hid"), label = NULL, value = 256, min = 1, step = 1)
        )
      ),
      column(3,
        .adv_param_block(
          "Attention heads (--n_heads)",
          numericInput(ns("n_heads"), label = NULL, value = 8, min = 1, step = 1)
        )
      ),
      column(3,
        .adv_param_block(
          "Transformer layers (--n_layers)",
          numericInput(ns("n_layers"), label = NULL, value = 2, min = 1, step = 1)
        )
      )
    ),

    # ── Section 3: Graph ────────────────────────────────────────────────────
    .adv_section_header("Graph"),

    fluidRow(
      column(3,
        .adv_param_block(
          "Graph topology (--topology)",
          selectInput(ns("topology"), label = NULL,
                      choices = c("full", "star"), selected = "full")
        )
      ),
      column(3,
        .adv_param_block(
          "Omic 1 corr. cutoff (--omic1_corr_cutoff)",
          numericInput(ns("omic1_corr_cutoff"), label = NULL,
                       value = 0.7, min = 0, max = 1, step = 0.05)
        )
      ),
      column(3,
        .adv_param_block(
          "Omic 2 corr. cutoff (--omic2_corr_cutoff)",
          numericInput(ns("omic2_corr_cutoff"), label = NULL,
                       value = 0.7, min = 0, max = 1, step = 0.05)
        )
      ),
      column(3,
        .adv_param_block(
          "Sample mode (--sample_mode)",
          selectInput(ns("sample_mode"), label = NULL,
                      choices = c("intersection", "union"), selected = "intersection")
        )
      )
    ),

    # ── Section 4: Run / Output ─────────────────────────────────────────────
    .adv_section_header("Run / Output"),

    fluidRow(
      column(3,
        .adv_param_block(
          "Leiden resolution (--leiden_resolution)",
          numericInput(ns("leiden_resolution"), label = NULL,
                       value = 1.5, min = 0.1, step = 0.1)
        )
      ),
      column(3,
        .adv_param_block(
          "Min cluster size (--leiden_min_cluster_size)",
          numericInput(ns("leiden_min_cluster_size"), label = NULL,
                       value = 3, min = 1, step = 1)
        )
      ),
      column(3,
             .adv_param_block(
               "Compute device (--cuda)",
               selectInput(ns("cuda"), label = NULL,
                           choices = c("GPU (0)" = "0", "CPU (1)" = "1"),
                           selected = "1")
             )
      ),
      
      column(3,
        .adv_param_block(
          "Save checkpoint every N epochs (--checkpoint_every)",
          numericInput(ns("checkpoint_every"), label = NULL, value = 1, min = 1, step = 1)
        )
      )
    ),

    br()
  )
}


# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------

advPramigoServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # -------------------------------------------------------------------------
    # Returned reactive: all 16 Pramigo parameter values
    # -------------------------------------------------------------------------
    reactive({
      list(
        epoch        = as.character(input$epoch        %||% 100),
        n_batch      = as.character(input$n_batch      %||% 25),
        lr           = as.character(input$lr           %||% 0.0001),
        resume       = input$resume                    %||% "off",
        
        n_hid        = as.character(input$n_hid        %||% 128),
        ae_n_hid     = as.character(input$ae_n_hid     %||% 256),
        n_heads      = as.character(input$n_heads      %||% 8),
        n_layers     = as.character(input$n_layers     %||% 2),
        
        topology     = input$topology                  %||% "full",
        omic1_cutoff = as.character(input$omic1_corr_cutoff %||% 0.7),
        omic2_cutoff = as.character(input$omic2_corr_cutoff %||% 0.7),
        sample_mode  = input$sample_mode               %||% "intersection",
        
        leiden_res   = as.character(input$leiden_resolution       %||% 1.5),
        leiden_min   = as.character(input$leiden_min_cluster_size %||% 3),
        cuda         = input$cuda                      %||% "1",
        checkpoint   = as.character(input$checkpoint_every        %||% 1)
      )
    })

  })
}
