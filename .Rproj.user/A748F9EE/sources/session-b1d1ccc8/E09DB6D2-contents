# Kiwifruit Softening Monte Carlo — Shiny R App (updated)
# -----------------------------------------------------------------
# Full rewrite with robust DT editing, 1–3 treatments, add/delete rows,
# minimal 2-row profiles (Reset), and faithful MATLAB model port.
# Parameter files expected in ./params (sheet 1 = Green, sheet 2 = Gold).
# -----------------------------------------------------------------

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
  library(readr)
})

# ---------------- Arrhenius helper (karr.m) ----------------
arrhenius_k <- function(kref, Ea, TempC, TrefC = 0) {
  Rgas <- 8.3143
  kref * exp((Ea / Rgas) * (1 / (TrefC + 273.15) - 1 / (TempC + 273.15)))
}

# ---------------- Triphase ODE (Triphase_firmness.m) ----------------
# y = c(y1_enzyme, y2_phase1, y3_phase2, y4_ox, y5_CI)
triphase_firmness <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    Tref <- 0
    if (identical(variety, "Green")) { # HW
      Ea_Enz <- 40577;  ki <- 0.0025; Ffix2 <- 0.1; vmax_Eth <- 32.253
      vmax_base_ref <- 6.3451; Km <- 2.528; kEnz <- 0.06094
      ksref <- 0.021758; kpref <- 0.00022615; Ea_s <- 64492; Ea_p <- 53699
      gamma <- 0.35611
    } else { # Gold (GA)
      Ea_Enz <- 41474;  ki <- 0.0025; Ffix2 <- 0.1; vmax_Eth <- 19.446
      Km <- 27.924; kEnz <- 0.067792; vmax_base_ref <- 5.4574
      ksref <- 0.021758; kpref <- 0.00019547; Ea_s <- 64493; Ea_p <- 43944
      gamma <- 0.05611
    }
    if (Ffix1 > F0) Ffix1 <- F0
    
    # Oxidative path
    kox_d_ref <- 0.72; Ea_ox_d <- 86000; kox_p <- 0.01
    
    # Interpolated environment
    T <- temp_fun(t)
    C2H4 <- c2h4_fun(t)
    
    # Temp dependencies
    kox_d <- arrhenius_k(kox_d_ref, Ea_ox_d, T, Tref)
    ks    <- arrhenius_k(ksref,     Ea_s,   T, Tref)
    kp    <- arrhenius_k(kpref,     Ea_p,   T, Tref)
    vmax_base <- arrhenius_k(vmax_base_ref, Ea_Enz, T, Tref)
    vmax <- vmax_base + vmax_Eth * C2H4 / (Km + C2H4)
    
    y1 <- y1_enzyme; y2 <- y2_phase1; y3 <- y3_phase2; y4 <- y4_ox; y5 <- y5_CI
    dy1 <- kEnz * y1 * (vmax - y1) - ki * y1
    dy2 <- ks * y1 * (F0 - y2 - Ffix1)
    dy3 <- kp * y2 * (Ffix1 - y3 - Ffix2) * y1 + gamma * y5 * (Ffix1 - y3 - Ffix2)
    dy4 <- kox_p - kox_d * y4
    dy5 <- y4 * y5 * (4 - y5)
    
    list(c(dy1, dy2, dy3, dy4, dy5))
  })
}

# --------------- Monte Carlo driver (MATLAB parity) ---------------
run_montecarlo_firmness <- function(params_df, treatments, variety = "Green",
                                    use_init = FALSE, F_n = 7, std_n = 1,
                                    n_time_steps = 51) {
  stopifnot(all(c("E0", "F0", "Ffix1") %in% names(params_df)))
  
  if (use_init) {
    m0 <- mean(params_df$F0); s0 <- stats::sd(params_df$F0)
    if (is.na(s0) || s0 == 0) s0 <- 1
    params_df$F0 <- (params_df$F0 - m0) * (std_n / s0) + F_n
  }
  
  make_env_funs <- function(profile_df) {
    profile_df <- profile_df %>% mutate(across(everything(), as.numeric)) %>% drop_na(day, tempC, c2h4)
    profile_df <- arrange(profile_df, day)
    if (nrow(profile_df) < 2) stop("Each treatment needs ≥ 2 rows: day 0 and a final day.")
    if (profile_df$day[1] != 0) {
      profile_df <- bind_rows(tibble(day = 0, tempC = profile_df$tempC[1], c2h4 = profile_df$c2h4[1]), profile_df)
    }
    profile_df <- profile_df %>% group_by(day) %>% summarise(tempC = mean(tempC), c2h4 = mean(c2h4), .groups = 'drop')
    t_end <- max(profile_df$day)
    temp_fun <- stats::approxfun(profile_df$day, profile_df$tempC, method = "linear", rule = 2, ties = mean)
    c2h4_fun <- stats::approxfun(profile_df$day, profile_df$c2h4, method = "linear", rule = 2, ties = mean)
    list(temp_fun = temp_fun, c2h4_fun = c2h4_fun, t_end = t_end)
  }
  
  all_trajs <- list(); summaries <- list()
  for (tr in treatments) {
    tr_name <- tr$name
    env <- make_env_funs(tr$profile_df)
    times <- seq(0, env$t_end, length.out = n_time_steps)
    
    trajs <- purrr::imap_dfr(seq_len(nrow(params_df)), function(i, idx) {
      E0 <- params_df$E0[i]; F0 <- params_df$F0[i]; Ffix1 <- params_df$Ffix1[i]
      y0 <- c(y1_enzyme = E0, y2_phase1 = 0.01, y3_phase2 = 0.01, y4_ox = 1e-05, y5_CI = 1e-05)
      parms <- list(F0 = F0, Ffix1 = Ffix1, variety = variety, temp_fun = env$temp_fun, c2h4_fun = env$c2h4_fun)
      ode_df <- as.data.frame(deSolve::ode(y = y0, times = times, func = triphase_firmness, parms = parms,
                                           method = "lsoda", atol = 1e-8, rtol = 1e-8))
      firmness <- (F0 + 0.02) - ode_df$y2_phase1 - ode_df$y3_phase2
      tibble(treatment = tr_name, time = ode_df$time, replicate = i, firmness = firmness)
    })
    
    all_trajs[[tr_name]] <- trajs
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
  list(trajectories = bind_rows(all_trajs), summary = bind_rows(summaries))
}

# --------------------------- Shiny UI ---------------------------
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
                 helpText("Place para_fast.xlsx, para_med.xlsx, para_slow.xlsx in this folder. Sheet 1 = Green, Sheet 2 = Gold."),
                 textInput("param_dir", "Parameter folder", value = "./params"),
                 hr(),
                 h4("Treatments (stepwise profiles)"),
                 helpText("Minimum two rows per treatment (start @ day 0 and a final day). Add intermediate steps as needed."),
                 selectInput("treat_count", "Number of treatments", choices = c(1,2,3), selected = 1),
                 tabsetPanel(id = "tr_tabs",
                             tabPanel("Treatment 1",
                                      textInput("tr1_name", "Treatment name", value = "T1"),
                                      fluidRow(
                                        column(8, DTOutput("tr1_table")),
                                        column(4,
                                               h5("Quick setup"),
                                               numericInput("tr1_final_day", "Final day", value = 12, min = 1, step = 1),
                                               numericInput("tr1_start_temp", "Start Temp (°C)", value = 10),
                                               numericInput("tr1_start_c2h4", "Start C2H4 (ppb)", value = 10),
                                               numericInput("tr1_final_temp", "Final Temp (°C)", value = 0),
                                               numericInput("tr1_final_c2h4", "Final C2H4 (ppb)", value = 10),
                                               actionButton("tr1_reset2", "Reset to minimal 2 rows"),
                                               br(), br(),
                                               actionButton("tr1_add", "+ Add row"),
                                               actionButton("tr1_del", "− Delete last"),
                                               helpText("Keep at least two rows (start & final).")
                                        )
                                      )
                             ),
                             conditionalPanel("input.treat_count >= 2",
                                              tabPanel("Treatment 2",
                                                       textInput("tr2_name", "Treatment name", value = "T2"),
                                                       fluidRow(
                                                         column(8, DTOutput("tr2_table")),
                                                         column(4,
                                                                h5("Quick setup"),
                                                                numericInput("tr2_final_day", "Final day", value = 10, min = 1, step = 1),
                                                                numericInput("tr2_start_temp", "Start Temp (°C)", value = 15),
                                                                numericInput("tr2_start_c2h4", "Start C2H4 (ppb)", value = 10),
                                                                numericInput("tr2_final_temp", "Final Temp (°C)", value = 2),
                                                                numericInput("tr2_final_c2h4", "Final C2H4 (ppb)", value = 10),
                                                                actionButton("tr2_reset2", "Reset to minimal 2 rows"),
                                                                br(), br(),
                                                                actionButton("tr2_add", "+ Add row"),
                                                                actionButton("tr2_del", "− Delete last"),
                                                                helpText("Keep at least two rows (start & final).")
                                                         )
                                                       )
                                              )
                             ),
                             conditionalPanel("input.treat_count >= 3",
                                              tabPanel("Treatment 3",
                                                       textInput("tr3_name", "Treatment name", value = "T3"),
                                                       fluidRow(
                                                         column(8, DTOutput("tr3_table")),
                                                         column(4,
                                                                h5("Quick setup"),
                                                                numericInput("tr3_final_day", "Final day", value = 7, min = 1, step = 1),
                                                                numericInput("tr3_start_temp", "Start Temp (°C)", value = 5),
                                                                numericInput("tr3_start_c2h4", "Start C2H4 (ppb)", value = 10),
                                                                numericInput("tr3_final_temp", "Final Temp (°C)", value = 5),
                                                                numericInput("tr3_final_c2h4", "Final C2H4 (ppb)", value = 10),
                                                                actionButton("tr3_reset2", "Reset to minimal 2 rows"),
                                                                br(), br(),
                                                                actionButton("tr3_add", "+ Add row"),
                                                                actionButton("tr3_del", "− Delete last"),
                                                                helpText("Keep at least two rows (start & final).")
                                                         )
                                                       )
                                              )
                             )
                 ),
                 br(),
                 actionButton("run", "Run Simulation", class = "btn btn-primary"),
                 br(), br(),
                 downloadButton("download_csv", "Download CSV results")
    ),
    mainPanel(width = 8,
              h4("Monte Carlo Firmness — Faceted by Treatment"),
              plotOutput("mc_plot", height = "560px"),
              hr(),
              h4("Summary at Final Day"),
              DTOutput("summary_table")
    )
  )
)

# -------------------------- Shiny Server --------------------------
server <- function(input, output, session) {
  # Minimal profile builder
  minimal_profile <- function(final_day, start_temp, start_c2h4, final_temp, final_c2h4) {
    tibble(day = c(0, final_day), tempC = c(start_temp, final_temp), c2h4 = c(start_c2h4, final_c2h4))
  }
  starter_profile <- function(last_day = 12, t1 = 10, t2 = 5, t3 = 0, c1 = 10, c2 = 10, c3 = 10) {
    tibble(day = c(0, 5, 10, last_day), tempC = c(t1, t2, t3, t3), c2h4 = c(c1, c2, c3, c3))
  }
  
  tr_tables <- reactiveValues(
    tr1 = starter_profile(),
    tr2 = starter_profile(last_day = 10, t1 = 15, t2 = 10, t3 = 2),
    tr3 = starter_profile(last_day = 7, t1 = 5, t2 = 5, t3 = 5)
  )
  
  render_tr_dt <- function(id, rv_name) {
    output[[id]] <- renderDT({
      datatable(tr_tables[[rv_name]], editable = TRUE, rownames = FALSE, selection = 'none',
                options = list(dom = 'tip', pageLength = 5))
    })
  }
  render_tr_dt("tr1_table", "tr1"); render_tr_dt("tr2_table", "tr2"); render_tr_dt("tr3_table", "tr3")
  
  proxy_update <- function(id, data) {
    replaceData(dataTableProxy(id), data, resetPaging = FALSE, rownames = FALSE)
  }
  
  # Add / delete row helpers (never below 2 rows)
  add_row <- function(df) {
    if (nrow(df) < 1) return(tibble(day = 0, tempC = 0, c2h4 = 0))
    last <- tail(df, 1)
    bind_rows(df, tibble(day = as.numeric(last$day) + 1, tempC = as.numeric(last$tempC), c2h4 = as.numeric(last$c2h4)))
  }
  del_row <- function(df) { if (nrow(df) <= 2) df else df[-nrow(df), ] }
  
  observeEvent(input$tr1_add, { tr_tables$tr1 <- add_row(tr_tables$tr1); proxy_update("tr1_table", tr_tables$tr1) })
  observeEvent(input$tr1_del, { tr_tables$tr1 <- del_row(tr_tables$tr1); proxy_update("tr1_table", tr_tables$tr1) })
  observeEvent(input$tr2_add, { tr_tables$tr2 <- add_row(tr_tables$tr2); proxy_update("tr2_table", tr_tables$tr2) })
  observeEvent(input$tr2_del, { tr_tables$tr2 <- del_row(tr_tables$tr2); proxy_update("tr2_table", tr_tables$tr2) })
  observeEvent(input$tr3_add, { tr_tables$tr3 <- add_row(tr_tables$tr3); proxy_update("tr3_table", tr_tables$tr3) })
  observeEvent(input$tr3_del, { tr_tables$tr3 <- del_row(tr_tables$tr3); proxy_update("tr3_table", tr_tables$tr3) })
  
  # Minimal 2-row resets using quick-setup inputs
  observeEvent(input$tr1_reset2, {
    tr_tables$tr1 <- minimal_profile(input$tr1_final_day, input$tr1_start_temp, input$tr1_start_c2h4, input$tr1_final_temp, input$tr1_final_c2h4)
    proxy_update("tr1_table", tr_tables$tr1)
  })
  observeEvent(input$tr2_reset2, {
    tr_tables$tr2 <- minimal_profile(input$tr2_final_day, input$tr2_start_temp, input$tr2_start_c2h4, input$tr2_final_temp, input$tr2_final_c2h4)
    proxy_update("tr2_table", tr_tables$tr2)
  })
  observeEvent(input$tr3_reset2, {
    tr_tables$tr3 <- minimal_profile(input$tr3_final_day, input$tr3_start_temp, input$tr3_start_c2h4, input$tr3_final_temp, input$tr3_final_c2h4)
    proxy_update("tr3_table", tr_tables$tr3)
  })
  
  # Robust DT edits: map visible col index -> name, numeric only
  observeEvent(input$tr1_table_cell_edit, {
    info <- input$tr1_table_cell_edit
    if (is.null(info$col) || info$col == 0) return()
    col_name <- colnames(tr_tables$tr1)[info$col]
    val <- suppressWarnings(as.numeric(info$value)); if (is.na(val)) return()
    if (col_name %in% c("day","tempC","c2h4")) {
      tr_tables$tr1[info$row, col_name] <- val; proxy_update("tr1_table", tr_tables$tr1)
    }
  })
  observeEvent(input$tr2_table_cell_edit, {
    info <- input$tr2_table_cell_edit
    if (is.null(info$col) || info$col == 0) return()
    col_name <- colnames(tr_tables$tr2)[info$col]
    val <- suppressWarnings(as.numeric(info$value)); if (is.na(val)) return()
    if (col_name %in% c("day","tempC","c2h4")) {
      tr_tables$tr2[info$row, col_name] <- val; proxy_update("tr2_table", tr_tables$tr2)
    }
  })
  observeEvent(input$tr3_table_cell_edit, {
    info <- input$tr3_table_cell_edit
    if (is.null(info$col) || info$col == 0) return()
    col_name <- colnames(tr_tables$tr3)[info$col]
    val <- suppressWarnings(as.numeric(info$value)); if (is.na(val)) return()
    if (col_name %in% c("day","tempC","c2h4")) {
      tr_tables$tr3[info$row, col_name] <- val; proxy_update("tr3_table", tr_tables$tr3)
    }
  })
  
  # Load parameters from Excel based on variety & softening type
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
    names(df) <- c("E0", "F0", "Ffix1")[seq_len(ncol(df))]
    df
  })
  
  # Run simulation
  sim_result <- eventReactive(input$run, {
    params_all <- load_params()
    N <- as.integer(input$mc_n); validate(need(N > 0, "Replicates must be > 0"))
    N_eff <- min(N, nrow(params_all))
    idx <- sample(seq_len(nrow(params_all)), N_eff, replace = FALSE)
    params_df <- params_all[idx, c("E0", "F0", "Ffix1")] %>% as.data.frame()
    
    # Build treatments list per treat_count
    tr_list <- list(list(name = input$tr1_name, profile_df = tr_tables$tr1))
    if (input$treat_count >= 2) tr_list[[length(tr_list)+1]] <- list(name = input$tr2_name, profile_df = tr_tables$tr2)
    if (input$treat_count >= 3) tr_list[[length(tr_list)+1]] <- list(name = input$tr3_name, profile_df = tr_tables$tr3)
    
    run_montecarlo_firmness(
      params_df = params_df,
      treatments = tr_list,
      variety = input$variety,
      use_init = isTRUE(input$use_init),
      F_n = input$F_n,
      std_n = input$std_n,
      n_time_steps = 51
    )
  })
  
  # Plot
  output$mc_plot <- renderPlot({
    res <- sim_result(); req(res$trajectories)
    df <- res$trajectories
    avg_df <- df %>% group_by(treatment, time) %>% summarise(mean_f = mean(firmness), .groups = 'drop')
    ggplot(df, aes(x = time, y = firmness, group = replicate)) +
      geom_line(alpha = 0.08, linewidth = 0.4, colour = "grey30") +
      geom_line(data = avg_df, aes(x = time, y = mean_f, group = NULL), linewidth = 1.2) +
      facet_wrap(~ treatment, scales = "free_y") +
      labs(x = "Time (days)", y = "Firmness (kgf)", subtitle = "Grey = individual fruit replicates; Solid = average firmness") +
      theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank())
  })
  
  # Summary table
  output$summary_table <- renderDT({
    res <- sim_result(); datatable(res$summary, rownames = FALSE, options = list(dom = 'tip'))
  })
  
  # Download all trajectories (long format)
  output$download_csv <- downloadHandler(
    filename = function() paste0("MC_Firmness_", Sys.Date(), ".csv"),
    content = function(file) {
      res <- sim_result(); readr::write_csv(res$trajectories, file)
    }
  )
}

shinyApp(ui, server)
