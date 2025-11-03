# Kiwifruit Softening Monte Carlo — Shiny R App
# ------------------------------------------------------------
# Reimplementation of MATLAB model in R with Shiny GUI
# - Arrhenius helper (karr)
# - Triphase firmness ODE (GA/HW)
# - Monte Carlo driver mirroring MATLAB logic
# - Shiny UI to define one or multiple treatments with stepwise Temp & C2H4
# - Outputs: faceted plot (replicates + average), summary table, CSV download
# - Parameter files expected by default in ./params: para_fast.xlsx, para_med.xlsx, para_slow.xlsx
#   Sheet 1 = Green (HW), Sheet 2 = Gold (GA)
# ------------------------------------------------------------

# --- Packages ---
suppressPackageStartupMessages({
  library(shiny)
  library(deSolve)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(DT)
  library(stringr)
})

# --------- Arrhenius helper (equivalent to MATLAB karr.m) ---------
arrhenius_k <- function(kref, Ea, TempC, TrefC = 0) {
  Rgas <- 8.3143
  kref * exp((Ea / Rgas) * (1 / (TrefC + 273.15) - 1 / (TempC + 273.15)))
}

# --------- ODE core (equivalent to Triphase_firmness.m) ---------
# y = c(y1_enzyme, y2_phase1, y3_phase2, y4_ox, y5_CI)
# parms: list with F0, Ffix1, variety ("Green" or "Gold"), temp_fun, c2h4_fun
triphase_firmness <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    # Variety-dependent parameters
    Tref <- 0
    if (identical(variety, "Green")) { # HW
      Ea_Enz <- 40577
      ki <- 0.0025
      Ffix2 <- 0.1
      vmax_Eth <- 32.253
      vmax_base_ref <- 6.3451
      Km <- 2.528
      kEnz <- 0.06094
      ksref <- 0.021758
      kpref <- 0.00022615
      Ea_s <- 64492
      Ea_p <- 53699
      gamma <- 0.35611
    } else { # Gold (GA)
      Ea_Enz <- 41474
      ki <- 0.0025
      Ffix2 <- 0.1
      vmax_Eth <- 19.446
      Km <- 27.924
      kEnz <- 0.067792
      vmax_base_ref <- 5.4574
      ksref <- 0.021758
      kpref <- 0.00019547
      Ea_s <- 64493
      Ea_p <- 43944
      gamma <- 0.05611
    }
    
    # ensure monotonicity: Ffix1 <= F0
    if (Ffix1 > F0) Ffix1 <- F0
    
    # Oxidative pathway parameters
    kox_d_ref <- 0.72
    Ea_ox_d <- 86000
    kox_p <- 0.01
    
    # Interpolated environment
    T <- temp_fun(t)
    C2H4 <- c2h4_fun(t)
    
    # Temperature dependencies (Arrhenius)
    kox_d <- arrhenius_k(kox_d_ref, Ea_ox_d, T, Tref)
    ks <- arrhenius_k(ksref, Ea_s, T, Tref)
    kp <- arrhenius_k(kpref, Ea_p, T, Tref)
    vmax_base <- arrhenius_k(vmax_base_ref, Ea_Enz, T, Tref)
    vmax <- vmax_base + vmax_Eth * C2H4 / (Km + C2H4)
    
    # State derivatives
    y1 <- y1_enzyme; y2 <- y2_phase1; y3 <- y3_phase2; y4 <- y4_ox; y5 <- y5_CI
    dy1 <- kEnz * y1 * (vmax - y1) - ki * y1
    dy2 <- ks * y1 * (F0 - y2 - Ffix1)
    dy3 <- kp * y2 * (Ffix1 - y3 - Ffix2) * y1 + gamma * y5 * (Ffix1 - y3 - Ffix2)
    dy4 <- kox_p - kox_d * y4
    dy5 <- y4 * y5 * (4 - y5) # CI from 0 to 4
    
    list(c(dy1, dy2, dy3, dy4, dy5))
  })
}

# ---------- Monte Carlo driver (R version of MonteCarloFirmness.m) ----------
# inputs:
#   params_df: data.frame(E0, F0, Ffix1) sampled N rows
#   treatments: list of treatment lists, each with name, profile_df (cols: day, tempC, c2h4)
#   variety: "Green" or "Gold"
#   use_init: TRUE/FALSE; if TRUE, recenter/rescale F0 to mean F_n and sd std_n
#   F_n, std_n: numeric targets for F0 distribution if use_init
#   n_time_steps: number of output points (MATLAB used ~51)
# returns: list with $trajectories (long df) and $summary (per-treatment %<=0.8 at final day)
run_montecarlo_firmness <- function(params_df, treatments, variety = "Green",
                                    use_init = FALSE, F_n = 7, std_n = 1,
                                    n_time_steps = 51) {
  stopifnot(all(c("E0", "F0", "Ffix1") %in% names(params_df)))
  
  # If requested, re-center/scale F0s to user targets
  if (use_init) {
    m0 <- mean(params_df$F0)
    s0 <- stats::sd(params_df$F0)
    if (is.na(s0) || s0 == 0) s0 <- 1
    params_df$F0 <- (params_df$F0 - m0) * (std_n / s0) + F_n
  }
  
  # Helper to build interpolation functions for a profile
  make_env_funs <- function(profile_df) {
    # Expect columns: day (numeric), tempC (numeric), c2h4 (numeric)
    profile_df <- arrange(profile_df, day)
    # Enforce first day == 0
    if (profile_df$day[1] != 0) {
      profile_df <- bind_rows(tibble(day = 0,
                                     tempC = profile_df$tempC[1],
                                     c2h4 = profile_df$c2h4[1]),
                              profile_df)
    }
    # Last day is explicit in input, used as t_end
    t_end <- max(profile_df$day)
    
    # Build linear interpolation functions
    temp_fun <- stats::approxfun(profile_df$day, profile_df$tempC, method = "linear",
                                 rule = 2, ties = mean)
    c2h4_fun <- stats::approxfun(profile_df$day, profile_df$c2h4, method = "linear",
                                 rule = 2, ties = mean)
    list(temp_fun = temp_fun, c2h4_fun = c2h4_fun, t_end = t_end)
  }
  
  all_trajs <- list()
  summaries <- list()
  
  for (tr in treatments) {
    tr_name <- tr$name
    env <- make_env_funs(tr$profile_df)
    times <- seq(0, env$t_end, length.out = n_time_steps)
    
    # Simulate each Monte Carlo draw
    trajs <- purrr::imap_dfr(seq_len(nrow(params_df)), function(i, idx) {
      E0 <- params_df$E0[i]
      F0 <- params_df$F0[i]
      Ffix1 <- params_df$Ffix1[i]
      
      y0 <- c(y1_enzyme = E0,
              y2_phase1 = 0.01,
              y3_phase2 = 0.01,
              y4_ox = 1e-05,
              y5_CI = 1e-05)
      
      parms <- list(F0 = F0, Ffix1 = Ffix1, variety = variety,
                    temp_fun = env$temp_fun, c2h4_fun = env$c2h4_fun)
      
      ode_out <- deSolve::ode(y = y0, times = times, func = triphase_firmness, parms = parms,
                              method = "lsoda", atol = 1e-8, rtol = 1e-8)
      ode_df <- as.data.frame(ode_out)
      # Firmness trajectory (match MATLAB’s F = (F0 + 0.02) - y2 - y3)
      firmness <- (F0 + 0.02) - ode_df$y2_phase1 - ode_df$y3_phase2
      
      tibble(treatment = tr_name,
             time = ode_df$time,
             replicate = i,
             firmness = firmness)
    })
    
    # Store per-treatment
    all_trajs[[tr_name]] <- trajs
    
    # Summary: % fruits <= 0.8 kgf on final simulated day
    final_day <- max(trajs$time)
    last_rows <- filter(trajs, time == final_day)
    pct_le_0_8 <- mean(last_rows$firmness <= 0.8) * 100
    
    summaries[[tr_name]] <- tibble(
      treatment = tr_name,
      final_day = final_day,
      pct_firmness_le_0_8 = round(pct_le_0_8, 2),
      mean_final_firmness = round(mean(last_rows$firmness), 3),
      sd_final_firmness = round(stats::sd(last_rows$firmness), 3),
      n = nrow(last_rows)
    )
  }
  
  traj_df <- bind_rows(all_trajs)
  sum_df  <- bind_rows(summaries)
  list(trajectories = traj_df, summary = sum_df)
}

# -------------------------- Shiny App --------------------------
ui <- fluidPage(
  titlePanel("Kiwifruit Softening Monte Carlo (R / Shiny)"),
  sidebarLayout(
    sidebarPanel(width = 4,
                 h4("Model & Data"),
                 radioButtons("variety", "Variety", choices = c("Green (HW)" = "Green", "Gold (GA)" = "Gold"), selected = "Green"),
                 selectInput("softening", "Softening type", choices = c("fast softening", "average softening", "slow softening"), selected = "average softening"),
                 numericInput("mc_n", "Monte Carlo replicates", value = 1000, min = 10, max = 3000, step = 10),
                 checkboxInput("use_init", "Use initial firmness distribution (recenter & rescale F0)", value = FALSE),
                 fluidRow(
                   column(6, numericInput("F_n", "Target mean F0", value = 7, step = 0.1)),
                   column(6, numericInput("std_n", "Target SD F0", value = 1, step = 0.1))
                 ),
                 hr(),
                 h4("Parameter files"),
                 helpText("Place: para_fast.xlsx, para_med.xlsx, para_slow.xlsx in the folder below. Sheet 1 = Green, Sheet 2 = Gold."),
                 textInput("param_dir", "Parameter folder", value = "./params"),
                 hr(),
                 h4("Treatments (stepwise profiles)"),
                 helpText("Define one treatment per tab. For each, enter rows with Day, Temp(°C), C2H4 (ppb). Ensure a row at Day 0 and a final Day."),
                 tabsetPanel(id = "tr_tabs",
                             tabPanel("Treatment 1",
                                      textInput("tr1_name", "Treatment name", value = "T1"),
                                      DTOutput("tr1_table")
                             ),
                             tabPanel("Treatment 2",
                                      textInput("tr2_name", "Treatment name", value = "T2"),
                                      DTOutput("tr2_table")
                             ),
                             tabPanel("Treatment 3",
                                      textInput("tr3_name", "Treatment name", value = "T3"),
                                      DTOutput("tr3_table")
                             )
                 ),
                 br(),
                 actionButton("run", "Run Simulation", class = "btn btn-primary"),
                 br(), br(),
                 downloadButton("download_csv", "Download CSV results")
    ),
    mainPanel(width = 8,
              h4("Monte Carlo Firmness — Faceted by Treatment"),
              plotOutput("mc_plot", height = "550px"),
              hr(),
              h4("Summary at Final Day"),
              DTOutput("summary_table")
    )
  )
)

server <- function(input, output, session) {
  # Example starter tables
  starter_profile <- function(last_day = 12, t1 = 10, t2 = 5, t3 = 0, c1 = 10, c2 = 10, c3 = 10) {
    # Example: 0d at 10°C & 10ppb → 5d at 5°C → 10d at 0°C → 12d hold 0°C
    tibble(
      day = c(0, 5, 10, last_day),
      tempC = c(t1, t2, t3, t3),
      c2h4 = c(c1, c2, c3, c3)
    )
  }
  
  tr_tables <- reactiveValues(
    tr1 = starter_profile(),
    tr2 = starter_profile(last_day = 10, t1 = 15, t2 = 10, t3 = 2),
    tr3 = starter_profile(last_day = 7, t1 = 5, t2 = 5, t3 = 5)
  )
  
  render_tr_dt <- function(id, rv_name) {
    output[[id]] <- renderDT({
      datatable(tr_tables[[rv_name]], editable = TRUE, rownames = FALSE,
                options = list(dom = 'tip', pageLength = 5))
    })
  }
  
  render_tr_dt("tr1_table", "tr1")
  render_tr_dt("tr2_table", "tr2")
  render_tr_dt("tr3_table", "tr3")
  
  proxy_update <- function(id, data) {
    replaceData(dataTableProxy(id), data, resetPaging = FALSE, rownames = FALSE)
  }
  
  observeEvent(input$tr1_table_cell_edit, {
    info <- input$tr1_table_cell_edit
    tr_tables$tr1[info$row, info$col] <- as.numeric(info$value)
    proxy_update("tr1_table", tr_tables$tr1)
  })
  observeEvent(input$tr2_table_cell_edit, {
    info <- input$tr2_table_cell_edit
    tr_tables$tr2[info$row, info$col] <- as.numeric(info$value)
    proxy_update("tr2_table", tr_tables$tr2)
  })
  observeEvent(input$tr3_table_cell_edit, {
    info <- input$tr3_table_cell_edit
    tr_tables$tr3[info$row, info$col] <- as.numeric(info$value)
    proxy_update("tr3_table", tr_tables$tr3)
  })
  
  # Load parameter workbook based on variety & softening type
  load_params <- reactive({
    req(dir.exists(input$param_dir))
    fn <- switch(input$softening,
                 "fast softening" = "para_fast.xlsx",
                 "average softening" = "para_med.xlsx",
                 "slow softening" = "para_slow.xlsx"
    )
    path <- file.path(input$param_dir, fn)
    validate(need(file.exists(path), paste0("File not found: ", path)))
    sheet_id <- ifelse(input$variety == "Green", 1, 2)
    df <- readxl::read_excel(path, sheet = sheet_id, col_names = FALSE)
    # Expect 3 cols: E0, F0, Ffix1
    names(df) <- c("E0", "F0", "Ffix1")[seq_len(ncol(df))]
    df
  })
  
  # Run simulation
  sim_result <- eventReactive(input$run, {
    params_all <- load_params()
    N <- as.integer(input$mc_n)
    validate(need(N > 0, "Replicates must be > 0"))
    
    # sample without replacement (cap at nrow)
    N_eff <- min(N, nrow(params_all))
    idx <- sample(seq_len(nrow(params_all)), N_eff, replace = FALSE)
    params_df <- params_all[idx, c("E0", "F0", "Ffix1")] %>% as.data.frame()
    
    # Build treatments list from UI
    treatments <- list(
      list(name = input$tr1_name, profile_df = tr_tables$tr1),
      list(name = input$tr2_name, profile_df = tr_tables$tr2),
      list(name = input$tr3_name, profile_df = tr_tables$tr3)
    )
    
    run_montecarlo_firmness(
      params_df = params_df,
      treatments = treatments,
      variety = input$variety,
      use_init = isTRUE(input$use_init),
      F_n = input$F_n,
      std_n = input$std_n,
      n_time_steps = 51
    )
  })
  
  # Plot
  output$mc_plot <- renderPlot({
    res <- sim_result()
    req(res$trajectories)
    
    df <- res$trajectories
    # Average per treatment x time
    avg_df <- df %>% group_by(treatment, time) %>% summarise(mean_f = mean(firmness), .groups = 'drop')
    
    ggplot(df, aes(x = time, y = firmness, group = replicate)) +
      geom_line(alpha = 0.08, linewidth = 0.4, colour = "grey30") +
      geom_line(data = avg_df, aes(x = time, y = mean_f, group = NULL), linewidth = 1.2) +
      facet_wrap(~ treatment, scales = "free_y") +
      labs(x = "Time (days)", y = "Firmness (kgf)",
           subtitle = "Grey = individual fruit replicates; Solid line = average firmness") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })
  
  # Summary table
  output$summary_table <- renderDT({
    res <- sim_result()
    datatable(res$summary, rownames = FALSE, options = list(dom = 'tip'))
  })
  
  # Download CSV results (long format + summary in separate sheets zipped)
  output$download_csv <- downloadHandler(
    filename = function() paste0("MC_Firmness_", Sys.Date(), ".csv"),
    content = function(file) {
      res <- sim_result()
      # Single CSV with trajectories; users can pivot later
      readr::write_csv(res$trajectories, file)
    }
  )
}

shinyApp(ui, server)
