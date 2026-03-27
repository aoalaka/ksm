# Kiwifruit Softening Monte Carlo — Shiny R App (v1.7.1)
# -----------------------------------------------------------------
# Changes from v1.7.0:
#   - Professional UI: custom flexbox sidebar (no bslib collapse toggle),
#     clean readable CSS with CSS vars, proper font sizes, polished cards
#   - No page-level scroll: viewport-filling tabs, internal panel scroll
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
  library(bslib)
  library(plotly)
  library(jsonlite)
  library(mvtnorm)
})

db_available <- requireNamespace("DBI",  quietly = TRUE) &&
                requireNamespace("odbc", quietly = TRUE)
if (db_available) suppressPackageStartupMessages({ library(DBI); library(odbc) })

# ── Helpers ─────────────────────────────────────────────────────
load_param_bank <- function(softening_type, variety) {
  params_dir <- file.path(getwd(), "params")
  if (!dir.exists(params_dir)) stop("Missing params folder: ", params_dir)
  fn <- switch(softening_type,
               "fast softening"    = "para_fast.xlsx",
               "average softening" = "para_med.xlsx",
               "slow softening"    = "para_slow.xlsx")
  path <- file.path(params_dir, fn)
  if (!file.exists(path)) stop("File not found: ", path)
  df <- readxl::read_excel(path, sheet = ifelse(variety == "Green", 1, 2), col_names = FALSE)
  names(df) <- c("E0", "F0", "Ffix1")[seq_len(ncol(df))]
  df
}

variety_db <- function(v) c("Green" = "HW", "Gold" = "GA")[v]

# Max plausible kiwifruit firmness (kgf) — anything above is a data entry error
FIRMNESS_MAX <- 12

mape_colour <- function(x) {
  dplyr::case_when(x <= 5  ~ "#27ae60",   # green — good
                   x <= 15 ~ "#e67e22",   # orange — acceptable
                   TRUE    ~ "#c0392b")   # red — poor
}

get_variety_defaults <- function(variety) {
  if (identical(variety, "Green")) {
    list(ksref = 0.021758, kpref = 0.00022615, gamma = 0.35611,
         Ea_s = 64492, Ea_p = 53699, kEnz = 0.06094, Km = 2.528,
         vmax_base_ref = 6.3451, Ea_Enz = 40577)
  } else {
    list(ksref = 0.021758, kpref = 0.00019547, gamma = 0.05611,
         Ea_s = 64493, Ea_p = 43944, kEnz = 0.067792, Km = 27.924,
         vmax_base_ref = 5.4574, Ea_Enz = 41474)
  }
}

# Build c2h4 step-function rows for a pulse application
# Returns tibble with columns day, tempC (NA), c2h4 — ready to bind_rows with a temp profile
build_eth_rows <- function(eth_df, max_day, eps = 0.01) {
  # eth_df: data.frame with columns ppm, day, duration (can be multiple rows)
  # Legacy single-value call: build_eth_rows(data.frame(ppm=100, day=2, duration=1), max_day)
  if (is.null(eth_df) || nrow(eth_df) == 0 || all(eth_df$ppm == 0) || max_day <= 0) {
    return(tibble(day = c(0, max(max_day, 1)), tempC = NA_real_, c2h4 = 0))
  }
  # Build breakpoints for each application
  all_days <- c(0, max_day)
  for (i in seq_len(nrow(eth_df))) {
    d <- eth_df$day[i]; dur <- eth_df$duration[i]
    all_days <- c(all_days,
      max(0, d - eps), min(d, max_day),
      min(d + dur, max_day), min(d + dur + eps, max_day))
  }
  all_days <- sort(unique(all_days))
  # For each timepoint, find the active ppm (last matching application wins)
  c2h4_vals <- rep(0, length(all_days))
  for (i in seq_len(nrow(eth_df))) {
    d <- eth_df$day[i]; dur <- eth_df$duration[i]; ppm <- eth_df$ppm[i]
    active <- all_days >= d & all_days <= d + dur
    c2h4_vals[active] <- pmax(c2h4_vals[active], ppm)
  }
  tibble(day = all_days, tempC = NA_real_, c2h4 = c2h4_vals)
}

# ── Arrhenius ───────────────────────────────────────────────────
arrhenius_k <- function(kref, Ea, TempC, TrefC = 0) {
  Rgas <- 8.3143
  kref * exp((Ea / Rgas) * (1 / (TrefC + 273.15) - 1 / (TempC + 273.15)))
}

# ── Triphase ODE ────────────────────────────────────────────────
triphase_firmness <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    Tref <- 0
    if (identical(variety, "Green")) {
      Ea_Enz <- 40577; ki <- 0.0025; Ffix2 <- 0.1; vmax_Eth <- 32.253
      vmax_base_ref <- 6.3451; Km <- 2.528; kEnz <- 0.06094
      ksref <- 0.021758; kpref <- 0.00022615; Ea_s <- 64492; Ea_p <- 53699; gamma <- 0.35611
    } else {
      Ea_Enz <- 41474; ki <- 0.0025; Ffix2 <- 0.1; vmax_Eth <- 19.446
      Km <- 27.924; kEnz <- 0.067792; vmax_base_ref <- 5.4574
      ksref <- 0.021758; kpref <- 0.00019547; Ea_s <- 64493; Ea_p <- 43944; gamma <- 0.05611
    }
    # Allow calibration overrides from parms
    pnm <- names(parms)
    if ("ksref" %in% pnm) ksref <- parms$ksref
    if ("kpref" %in% pnm) kpref <- parms$kpref
    if ("gamma_cal" %in% pnm) gamma <- parms$gamma_cal
    if (Ffix1 > F0) Ffix1 <- F0
    kox_d_ref <- 0.72; Ea_ox_d <- 86000; kox_p <- 0.01
    T <- temp_fun(t); C2H4 <- c2h4_fun(t)
    kox_d     <- arrhenius_k(kox_d_ref,    Ea_ox_d, T, Tref)
    ks        <- arrhenius_k(ksref,         Ea_s,   T, Tref)
    kp        <- arrhenius_k(kpref,         Ea_p,   T, Tref)
    vmax_base <- arrhenius_k(vmax_base_ref, Ea_Enz, T, Tref)
    vmax      <- vmax_base + vmax_Eth * C2H4 / (Km + C2H4)
    y1 <- y1_enzyme; y2 <- y2_phase1; y3 <- y3_phase2; y4 <- y4_ox; y5 <- y5_CI
    dy1 <- kEnz * y1 * (vmax - y1) - ki * y1
    dy2 <- ks * y1 * (F0 - y2 - Ffix1)
    dy3 <- kp * y2 * (Ffix1 - y3 - Ffix2) * y1 + gamma * y5 * (Ffix1 - y3 - Ffix2)
    dy4 <- kox_p - kox_d * y4
    dy5 <- y4 * y5 * (4 - y5)
    list(c(dy1, dy2, dy3, dy4, dy5))
  })
}

# ── Monte Carlo driver ──────────────────────────────────────────
# profile_df cols: day (required), tempC (NA ok), c2h4 (NA ok).
# Temp and c2h4 are interpolated independently from their non-NA rows.
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
    profile_df <- profile_df %>% mutate(across(c(day, tempC, c2h4), as.numeric))
    temp_rows <- profile_df %>% filter(!is.na(day), !is.na(tempC)) %>%
      group_by(day) %>% summarise(tempC = mean(tempC), .groups = "drop") %>% arrange(day)
    c2h4_rows <- profile_df %>% filter(!is.na(day), !is.na(c2h4)) %>%
      group_by(day) %>% summarise(c2h4 = mean(c2h4), .groups = "drop") %>% arrange(day)
    if (nrow(temp_rows) < 2) stop("Each treatment needs >= 2 temperature data points.")
    if (nrow(c2h4_rows) < 2) stop("Each treatment needs >= 2 ethylene data points.")
    if (temp_rows$day[1] != 0)
      temp_rows <- bind_rows(tibble(day = 0, tempC = temp_rows$tempC[1]), temp_rows)
    if (c2h4_rows$day[1] != 0)
      c2h4_rows <- bind_rows(tibble(day = 0, c2h4 = 0), c2h4_rows)
    t_end    <- max(max(temp_rows$day), max(c2h4_rows$day))
    temp_fun <- stats::approxfun(temp_rows$day, temp_rows$tempC, method = "linear", rule = 2, ties = mean)
    c2h4_fun <- stats::approxfun(c2h4_rows$day, c2h4_rows$c2h4, method = "linear", rule = 2, ties = mean)
    list(temp_fun = temp_fun, c2h4_fun = c2h4_fun, t_end = t_end)
  }
  all_trajs <- list(); summaries <- list()
  for (tr in treatments) {
    tr_name <- tr$name; env <- make_env_funs(tr$profile_df)
    times <- seq(0, env$t_end, length.out = n_time_steps)
    trajs <- purrr::imap_dfr(seq_len(nrow(params_df)), function(i, idx) {
      E0 <- params_df$E0[i]; F0 <- params_df$F0[i]; Ffix1 <- params_df$Ffix1[i]
      y0 <- c(y1_enzyme = E0, y2_phase1 = 0.01, y3_phase2 = 0.01, y4_ox = 1e-05, y5_CI = 1e-05)
      parms <- list(F0 = F0, Ffix1 = Ffix1, variety = variety,
                    temp_fun = env$temp_fun, c2h4_fun = env$c2h4_fun)
      ode_df   <- as.data.frame(deSolve::ode(y = y0, times = times, func = triphase_firmness,
                                              parms = parms, method = "lsoda", atol = 1e-8, rtol = 1e-8))
      firmness <- (F0 + 0.02) - ode_df$y2_phase1 - ode_df$y3_phase2
      tibble(treatment = tr_name, time = ode_df$time, replicate = i, firmness = firmness)
    })
    all_trajs[[tr_name]] <- trajs
    final_day <- max(trajs$time); last_rows <- dplyr::filter(trajs, time == final_day)
    pct_le_0_8 <- mean(last_rows$firmness <= 0.8) * 100
    summaries[[tr_name]] <- tibble(treatment = tr_name, final_day = final_day,
      pct_firmness_le_0_8 = round(pct_le_0_8, 2),
      mean_final_firmness = round(mean(last_rows$firmness), 3),
      sd_final_firmness   = round(stats::sd(last_rows$firmness), 3),
      n = nrow(last_rows))
  }
  list(trajectories = bind_rows(all_trajs), summary = bind_rows(summaries))
}

# ── Bayesian ODE Calibration ──────────────────────────────────────
calibration_neglogpost <- function(theta, obs_df, temp_fun, c2h4_fun, variety,
                                    prior_means, prior_sds, param_names, t_end,
                                    sigma_obs = 0.5) {
  E0    <- as.numeric(theta["E0"])
  F0    <- as.numeric(theta["F0"])
  Ffix1 <- as.numeric(theta["Ffix1"])
  if (Ffix1 >= F0) return(1e8)
  parms <- list(F0 = F0, Ffix1 = Ffix1, variety = variety,
                temp_fun = temp_fun, c2h4_fun = c2h4_fun)
  for (nm in setdiff(param_names, c("E0", "F0", "Ffix1"))) {
    if (nm == "gamma") parms[["gamma_cal"]] <- as.numeric(theta[nm])
    else parms[[nm]] <- as.numeric(theta[nm])
  }
  y0 <- c(y1_enzyme = E0, y2_phase1 = 0.01, y3_phase2 = 0.01,
          y4_ox = 1e-05, y5_CI = 1e-05)
  times <- sort(unique(c(seq(0, t_end, length.out = 51), obs_df$day)))
  times <- times[times >= 0]
  sol <- tryCatch(
    suppressWarnings(
      as.data.frame(deSolve::ode(y = y0, times = times, func = triphase_firmness,
                                  parms = parms, method = "lsoda",
                                  atol = 1e-6, rtol = 1e-6, maxsteps = 5000))),
    error = function(e) NULL)
  if (is.null(sol) || nrow(sol) < 2) return(1e8)
  pred_f_all <- (F0 + 0.02) - sol$y2_phase1 - sol$y3_phase2
  if (any(is.na(pred_f_all)) || any(!is.finite(pred_f_all))) return(1e8)
  pred_fun   <- stats::approxfun(sol$time, pred_f_all, rule = 2)
  pred_at_obs <- pred_fun(obs_df$day)
  if (any(is.na(pred_at_obs))) return(1e8)
  ll <- sum(dnorm(obs_df$firmness, pred_at_obs, sigma_obs, log = TRUE))
  theta_vals <- as.numeric(theta[param_names])
  pm_vals <- as.numeric(prior_means[param_names])
  ps_vals <- as.numeric(prior_sds[param_names])
  lp <- sum(dnorm(theta_vals, pm_vals, ps_vals, log = TRUE))
  result <- -(ll + lp)
  if (!is.finite(result)) return(1e8)
  result
}

run_calibration <- function(obs_df, temp_fun, c2h4_fun, variety,
                             param_bank, calibrate_secondary = FALSE,
                             n_posterior_samples = 500,
                             obs_init_firmness = NULL) {
  # Use observed initial firmness to inform F0 prior if available
  f0_obs <- if (!is.null(obs_init_firmness)) obs_init_firmness else mean(obs_df$firmness)
  prior_means <- c(E0 = mean(param_bank$E0),
                   F0 = f0_obs,
                   Ffix1 = mean(param_bank$Ffix1))
  prior_sds <- c(E0 = max(sd(param_bank$E0), 0.1),
                 F0 = max(sd(param_bank$F0), 1),
                 Ffix1 = max(sd(param_bank$Ffix1), 0.5))
  lower <- c(E0 = 0.001, F0 = max(0.5, f0_obs * 0.5), Ffix1 = 0.05)
  upper <- c(E0 = max(param_bank$E0) * 3,
             F0 = max(f0_obs * 1.5, max(param_bank$F0) * 1.5),
             Ffix1 = max(param_bank$Ffix1) * 3)
  param_names <- c("E0", "F0", "Ffix1")
  if (calibrate_secondary) {
    defaults <- get_variety_defaults(variety)
    for (nm in c("ksref", "kpref", "gamma")) {
      prior_means[nm] <- defaults[[nm]]
      prior_sds[nm]   <- defaults[[nm]] * 0.3
      lower[nm]       <- defaults[[nm]] * 0.05
      upper[nm]       <- defaults[[nm]] * 10
    }
    param_names <- c(param_names, "ksref", "kpref", "gamma")
  }
  # Estimate sigma from data spread
  sigma_obs <- max(sd(obs_df$firmness), 0.3)

  theta0 <- prior_means[param_names]
  t_end  <- max(obs_df$day)
  best <- list(value = Inf)

  # Try L-BFGS-B with multiple restarts
  for (restart in 1:5) {
    init <- if (restart == 1) theta0 else theta0 * runif(length(theta0), 0.8, 1.2)
    init <- pmax(init, lower[param_names])
    init <- pmin(init, upper[param_names])
    # Ensure Ffix1 < F0
    if (init["Ffix1"] >= init["F0"]) init["Ffix1"] <- init["F0"] * 0.5
    fit <- tryCatch(
      optim(init, calibration_neglogpost, method = "L-BFGS-B",
            lower = lower[param_names], upper = upper[param_names], hessian = TRUE,
            obs_df = obs_df, temp_fun = temp_fun, c2h4_fun = c2h4_fun,
            variety = variety, prior_means = prior_means, prior_sds = prior_sds,
            param_names = param_names, t_end = t_end, sigma_obs = sigma_obs),
      error = function(e) list(value = Inf))
    if (fit$value < best$value) best <- fit
  }
  # Fallback: Nelder-Mead (no gradient needed, more robust)
  if (best$value >= 1e7) {
    for (restart in 1:3) {
      init <- theta0 * runif(length(theta0), 0.8, 1.2)
      init <- pmax(init, lower[param_names])
      init <- pmin(init, upper[param_names])
      if (init["Ffix1"] >= init["F0"]) init["Ffix1"] <- init["F0"] * 0.5
      fit <- tryCatch(
        optim(init, calibration_neglogpost, method = "Nelder-Mead",
              control = list(maxit = 2000),
              obs_df = obs_df, temp_fun = temp_fun, c2h4_fun = c2h4_fun,
              variety = variety, prior_means = prior_means, prior_sds = prior_sds,
              param_names = param_names, t_end = t_end, sigma_obs = sigma_obs),
        error = function(e) list(value = Inf))
      if (fit$value < best$value) {
        # Get hessian at Nelder-Mead solution via finite differences
        best <- fit
        best$hessian <- tryCatch(
          optimHess(best$par, calibration_neglogpost,
                    obs_df = obs_df, temp_fun = temp_fun, c2h4_fun = c2h4_fun,
                    variety = variety, prior_means = prior_means, prior_sds = prior_sds,
                    param_names = param_names, t_end = t_end, sigma_obs = sigma_obs),
          error = function(e) NULL)
      }
    }
  }
  if (best$value >= 1e7) return(NULL)

  # Laplace approximation for posterior
  H <- best$hessian
  if (is.null(H)) {
    Sigma_post <- diag(prior_sds[param_names]^2)
  } else {
    Sigma_post <- tryCatch(solve(H), error = function(e) diag(prior_sds[param_names]^2))
  }
  eig <- eigen(Sigma_post, symmetric = TRUE)
  if (any(eig$values <= 0)) {
    eig$values[eig$values <= 0] <- 1e-6
    Sigma_post <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  }
  samples <- mvtnorm::rmvnorm(n_posterior_samples, mean = best$par, sigma = Sigma_post)
  colnames(samples) <- param_names
  for (nm in param_names) {
    samples[, nm] <- pmax(samples[, nm], lower[nm])
    samples[, nm] <- pmin(samples[, nm], upper[nm])
  }
  list(map_estimate = best$par, posterior_cov = Sigma_post,
       posterior_samples = as.data.frame(samples),
       prior_means = prior_means[param_names], prior_sds = prior_sds[param_names],
       convergence = best$convergence, neglogpost = best$value,
       param_names = param_names)
}

# ============================= CSS =============================

app_css <- "
/* ═══ Design tokens ═══ */
:root {
  --font: 'Inter', 'Segoe UI', system-ui, -apple-system, sans-serif;
  --c-primary: #4f46e5;
  --c-primary-h: #4338ca;
  --c-success: #059669;
  --c-success-h: #047857;
  --c-danger: #dc2626;
  --c-bg: #f8fafc;
  --c-surface: #ffffff;
  --c-border: #e2e8f0;
  --c-border-light: #f1f5f9;
  --c-text: #1e293b;
  --c-muted: #64748b;
  --c-dim: #94a3b8;
  --nav-h: 56px;
  --sidebar-w: 300px;
  --radius: 8px;
  --radius-lg: 12px;
}

body { font-family: var(--font) !important; background: var(--c-bg) !important; }

/* ── Navbar ──────────────────────────────── */
.navbar {
  background: #0f172a !important;
  min-height: var(--nav-h);
  padding: 0 1.5rem !important;
  border: none !important;
  box-shadow: 0 1px 3px rgba(0,0,0,0.12);
}
.navbar-brand {
  color: #f1f5f9 !important;
  font-weight: 700 !important;
  font-size: 0.95rem !important;
  letter-spacing: -0.01em;
}
.nav-link {
  color: rgba(248,250,252,0.55) !important;
  font-size: 0.82rem !important;
  font-weight: 500 !important;
  padding: 0.42rem 0.9rem !important;
  border-radius: 6px !important;
  margin: 0 2px;
  transition: all 0.15s ease;
}
.nav-link:hover {
  color: #e2e8f0 !important;
  background: rgba(255,255,255,0.07) !important;
}
.nav-link.active, .nav-link[aria-selected=true] {
  color: #ffffff !important;
  background: rgba(99,102,241,0.35) !important;
  font-weight: 600 !important;
}

/* ── Viewport fill ───────────────────────── */
.bslib-page-navbar > .tab-content {
  height: calc(100vh - var(--nav-h));
  overflow: hidden;
}
.bslib-page-navbar > .tab-content > .tab-pane {
  height: 100%;
  overflow: hidden;
}

/* ── Sidebar / main layout ───────────────── */
.app-row {
  display: flex;
  height: 100%;
  overflow: hidden;
}
.app-sidebar {
  width: var(--sidebar-w);
  min-width: var(--sidebar-w);
  background: var(--c-surface);
  border-right: 1px solid var(--c-border);
  overflow-y: auto;
  padding: 1.25rem;
  scrollbar-width: thin;
}
.app-sidebar::-webkit-scrollbar { width: 5px; }
.app-sidebar::-webkit-scrollbar-thumb { background: #cbd5e1; border-radius: 3px; }
.app-main {
  flex: 1;
  overflow-y: auto;
  padding: 1.25rem;
  scrollbar-width: thin;
}
.app-main::-webkit-scrollbar { width: 5px; }
.app-main::-webkit-scrollbar-thumb { background: #cbd5e1; border-radius: 3px; }

/* ── Cards ────────────────────────────────── */
.card {
  border: 1px solid var(--c-border) !important;
  border-radius: var(--radius-lg) !important;
  box-shadow: 0 1px 2px rgba(0,0,0,0.04) !important;
  background: var(--c-surface);
  margin-bottom: 1rem;
}
.kpi-card {
  border: none !important;
  box-shadow: 0 2px 8px rgba(0,0,0,0.12) !important;
  color: white !important;
}
.kpi-card .card-body * { color: inherit !important; }
.card:last-child { margin-bottom: 0; }
.card-header {
  background: var(--c-surface) !important;
  border-bottom: 1px solid var(--c-border-light) !important;
  padding: 0.75rem 1.15rem !important;
  font-size: 0.82rem !important;
  font-weight: 600 !important;
  color: var(--c-text) !important;
  border-radius: var(--radius-lg) var(--radius-lg) 0 0 !important;
}
.card-body { padding: 1rem 1.15rem !important; overflow: hidden; }

/* ── Forms ────────────────────────────────── */
label, .form-label, .control-label {
  font-size: 0.78rem !important;
  font-weight: 500 !important;
  color: var(--c-text) !important;
  margin-bottom: 0.2rem !important;
}
.form-control, .form-select {
  font-size: 0.84rem !important;
  border: 1px solid var(--c-border) !important;
  border-radius: var(--radius) !important;
  padding: 0.42rem 0.75rem !important;
  transition: border-color 0.15s, box-shadow 0.15s;
}
.form-control:focus, .form-select:focus {
  border-color: var(--c-primary) !important;
  box-shadow: 0 0 0 3px rgba(79,70,229,0.1) !important;
}
.shiny-input-container { margin-bottom: 0.65rem !important; }
.form-check-input:checked {
  background-color: var(--c-primary) !important;
  border-color: var(--c-primary) !important;
}
.radio label, .checkbox label, .form-check-label {
  font-size: 0.82rem !important;
  color: var(--c-text) !important;
}
.form-text, .shiny-help-block {
  font-size: 0.75rem !important;
  color: var(--c-dim) !important;
}

/* ── Buttons ──────────────────────────────── */
.btn {
  border-radius: var(--radius) !important;
  font-size: 0.82rem !important;
  font-weight: 500 !important;
  padding: 0.5rem 1rem !important;
  transition: all 0.15s ease;
}
.btn-primary     { background: var(--c-primary) !important; border: none !important; color: #fff !important; }
.btn-primary:hover { background: var(--c-primary-h) !important; }
.btn-success     { background: var(--c-success) !important; border: none !important; color: #fff !important; }
.btn-success:hover { background: var(--c-success-h) !important; }
.btn-outline-primary {
  border: 1px solid var(--c-primary) !important;
  color: var(--c-primary) !important;
  background: transparent !important;
}
.btn-outline-secondary {
  border: 1px solid var(--c-border) !important;
  color: var(--c-muted) !important;
  background: transparent !important;
}
.btn-outline-secondary:hover { background: var(--c-bg) !important; }
.btn.w-100 { display: flex; align-items: center; justify-content: center; gap: 0.4rem; }

/* ── Misc ─────────────────────────────────── */
hr { border-color: var(--c-border) !important; margin: 0.85rem 0 !important; opacity: 0.7; }
.alert {
  border: 1px solid transparent !important;
  font-size: 0.82rem !important;
  border-radius: var(--radius) !important;
  padding: 0.6rem 1rem !important;
}
.alert-success { background: #f0fdf4 !important; color: #166534 !important; border-color: #bbf7d0 !important; }
.alert-danger  { background: #fef2f2 !important; color: #991b1b !important; border-color: #fecaca !important; }
.alert-info    { background: #eff6ff !important; color: #1e40af !important; border-color: #bfdbfe !important; }
.alert-warning { background: #fffbeb !important; color: #92400e !important; border-color: #fde68a !important; }

/* ── Section labels ───────────────────────── */
.sec-lbl {
  font-size: 0.68rem;
  font-weight: 700;
  text-transform: uppercase;
  letter-spacing: 0.08em;
  color: var(--c-dim);
  margin: 1rem 0 0.4rem;
  display: block;
}
.sec-lbl:first-child { margin-top: 0; }

/* ── Ethylene box ─────────────────────────── */
.eth-box {
  background: var(--c-bg);
  border-radius: var(--radius);
  padding: 0.85rem;
  border: 1px solid var(--c-border);
  margin-top: 0.3rem;
  overflow: hidden;
}
.eth-box .form-group { margin-bottom: 0.5rem; }

/* ── DataTables ───────────────────────────── */
.dataTables_wrapper { font-size: 0.82rem !important; }
table.dataTable thead th {
  border-bottom: 2px solid var(--c-border) !important;
  font-weight: 600 !important;
  font-size: 0.76rem !important;
  color: var(--c-muted) !important;
  text-transform: uppercase;
  letter-spacing: 0.04em;
  padding: 0.6rem 0.75rem !important;
}
table.dataTable tbody td {
  border-bottom: 1px solid var(--c-border-light) !important;
  padding: 0.45rem 0.75rem !important;
}
table.dataTable tbody tr:hover { background: #f8fafc !important; }

/* ── Treatment tabs ───────────────────────── */
.nav-tabs { border-bottom: none !important; gap: 4px; padding: 4px 0; }
.nav-tabs .nav-link {
  color: var(--c-muted) !important;
  background: transparent !important;
  font-size: 0.82rem !important;
  font-weight: 500 !important;
  border: 1px solid transparent !important;
  padding: 0.4rem 0.9rem !important;
  margin: 0 !important;
  border-radius: 6px !important;
  transition: all 0.15s ease;
}
.nav-tabs .nav-link:hover {
  color: var(--c-text) !important;
  background: rgba(99,102,241,0.06) !important;
  border-color: #e2e8f0 !important;
}
.nav-tabs .nav-link.active {
  color: var(--c-primary) !important;
  background: rgba(99,102,241,0.08) !important;
  border-color: rgba(99,102,241,0.25) !important;
  font-weight: 600 !important;
}

/* Card header small buttons */
.card-header .btn {
  padding: 0.28rem 0.65rem !important;
  font-size: 0.75rem !important;
}

/* Content scroll helper */
.main-scroll {
  display: flex;
  flex-direction: column;
  gap: 1rem;
}
"

# ============================= UI =============================

ui <- page_navbar(
  id    = "main_nav",
  title = tags$span(
    tags$i(class = "bi bi-graph-up me-2", style = "font-size: 1rem;"),
    "Kiwifruit Firmness Simulator"
  ),
  theme = bs_theme(
    version = 5,
    bg = "#f8fafc", fg = "#1e293b",
    primary = "#4f46e5",
    secondary = "#64748b",
    success = "#059669",
    danger  = "#dc2626",
    warning = "#d97706",
    "navbar-bg" = "#0f172a"
  ),
  header = tags$head(
    tags$link(rel = "stylesheet",
      href = "https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.css"),
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "preconnect", href = "https://fonts.gstatic.com", crossorigin = ""),
    tags$link(rel = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"),
    tags$style(HTML(app_css))
  ),

  # ── Setup ─────────────────────────────────────────────────────
  nav_panel("Setup",
    div(class = "app-row",
      div(class = "app-sidebar",
        span(class = "sec-lbl", "Model & Data"),
        radioButtons("variety", "Variety",
                     choices = c("Green (HW)" = "Green", "Gold (GA)" = "Gold"),
                     selected = "Green"),
        selectInput("softening", "Softening type",
                    choices = c("fast softening", "average softening", "slow softening"),
                    selected = "average softening"),
        layout_columns(col_widths = c(7, 5),
          numericInput("mc_n", "Sample size", value = 1000, min = 10, max = 3000, step = 10),
          numericInput("F_n", "Mean F\u2080 (kgf)", value = 7, step = 0.1)
        ),
        checkboxInput("use_init", "Rescale F\u2080 to measured distribution", FALSE),
        conditionalPanel("input.use_init == true",
          numericInput("std_n", "SD", value = 1, step = 0.1)
        ),
        hr(),
        selectInput("treat_count", "Number of treatments",
                    choices = c(1, 2, 3), selected = 1),
        hr(),
        span(class = "sec-lbl", "Ethylene application"),
        checkboxInput("no_c2h4", "No ethylene (C\u2082H\u2084 = 0)", FALSE),
        div(class = "eth-box",
          conditionalPanel("input.no_c2h4 == false",
            DTOutput("eth_table", height = "auto"),
            div(class = "d-flex gap-1 mt-1",
              actionButton("eth_add", "+ Add", class = "btn btn-outline-secondary btn-sm"),
              actionButton("eth_rm",  "- Remove last", class = "btn btn-outline-secondary btn-sm")
            )
          ),
          conditionalPanel("input.no_c2h4 == true",
            tags$p(class = "text-muted mb-0", style = "font-size: 0.78rem;",
                   "Ethylene concentration set to zero for all treatments.")
          )
        ),
        hr(),
        actionButton("run", tagList(tags$i(class = "bi bi-play-fill"), "Run Simulation"),
                     class = "btn btn-primary w-100")
      ),
      div(class = "app-main",
        card(
          card_header(div(class = "d-flex justify-content-between align-items-center",
            span("Treatment Profiles"),
            tags$small(class = "text-muted", style = "font-weight: 400;",
                       "Double-click cells to edit  \u2022  Ethylene is set in the sidebar")
          )),
          card_body(uiOutput("treat_tabs"))
        )
      )
    )
  ),

  # ── Results ───────────────────────────────────────────────────
  nav_panel("Results",
    div(class = "app-main", style = "height: 100%;",
      div(class = "main-scroll",
        card(full_screen = TRUE,
          card_header(div(class = "d-flex justify-content-between align-items-center",
            span("Simulated Firmness \u2014 by Treatment"),
            downloadButton("download_csv", "Download CSV",
                           class = "btn btn-sm btn-outline-secondary")
          )),
          plotOutput("mc_plot", height = "420px")
        ),
        card(
          card_header("Summary at Final Day"),
          DTOutput("summary_table")
        )
      )
    )
  ),

  # ── Validate ──────────────────────────────────────────────────
  nav_panel("Validate",
    if (!db_available) {
      div(class = "alert alert-warning m-4",
          tags$strong("Database packages unavailable."),
          " Install ", tags$code("DBI"), ", ", tags$code("odbc"),
          " and ODBC Driver 17 for SQL Server.")
    } else {
      div(class = "app-row",
        div(class = "app-sidebar",
          uiOutput("db_status_ui"),
          actionButton("db_connect", "Connect to Database",
                       class = "btn btn-outline-primary w-100",
                       icon = shiny::icon("database")),
          hr(),
          span(class = "sec-lbl", "Filters"),
          uiOutput("val_year_ui"),
          radioButtons("val_location", "FMT Test Location",
                       choices = c("Shore", "Vessel"), selected = "Shore", inline = TRUE),
          uiOutput("val_voyage_ui"),
          uiOutput("val_hatch_ui"),
          hr(),
          span(class = "sec-lbl", "Model"),
          radioButtons("val_variety", "Variety",
                       choices = c("Green (HW)" = "Green", "Gold (GA)" = "Gold"),
                       selected = "Green"),
          selectInput("val_softening", "Softening type",
                      choices = c("fast softening", "average softening", "slow softening"),
                      selected = "average softening"),
          numericInput("val_mc_n", "Replicates per hatch",
                       value = 40, min = 10, max = 2000, step = 10),
          hr(),
          span(class = "sec-lbl", "Ethylene application"),
          checkboxInput("val_no_c2h4", "No ethylene (C\u2082H\u2084 = 0)", FALSE),
          div(class = "eth-box",
            conditionalPanel("input.val_no_c2h4 == false",
              DTOutput("val_eth_table", height = "auto"),
              div(class = "d-flex gap-1 mt-1",
                actionButton("val_eth_add", "+ Add", class = "btn btn-outline-secondary btn-sm"),
                actionButton("val_eth_rm",  "- Remove last", class = "btn btn-outline-secondary btn-sm")
              )
            ),
            conditionalPanel("input.val_no_c2h4 == true",
              tags$p(class = "text-muted mb-0", style = "font-size: 0.78rem;",
                     "Ethylene concentration set to zero.")
            )
          ),
          hr(),
          helpText(tags$i(class = "bi bi-info-circle me-1"),
                   "Initial firmness is auto-derived from each hatch's earliest FMT measurement."),
          actionButton("run_validate",
                       tagList(tags$i(class = "bi bi-check-circle"), "Run Validation"),
                       class = "btn btn-success w-100")
        ),
        div(class = "app-main",
          div(class = "main-scroll",
            uiOutput("val_kpi_cards"),
            card(
              card_header("Temperature Profiles \u2014 Cargo Condition Setpoints by Hatch"),
              plotlyOutput("val_temp_plot", height = "220px")
            ),
            card(full_screen = TRUE,
              card_header(div(class = "d-flex justify-content-between align-items-center",
                div(
                  tags$strong("Simulated vs Observed Firmness \u2014 by Hatch/Deck"),
                  br(),
                  tags$small(class = "text-muted",
                    tags$span(style = "color: #6baed6;", "\u2014"), " Sim replicates  \u2502  ",
                    tags$span(style = "color: #2c3e50; font-weight: bold;", "\u2014"), " Sim mean  \u2502  ",
                    tags$span(style = "color: #e74c3c; font-weight: bold;", "\u25cf"),
                    " Observed FMT  \u2502  ",
                    tags$span(style = "color: #e74c3c;", "\u2014"), " Obs mean  \u2502  ",
                    tags$span(style = "color: #e67e22;", "- - -"),
                    " 0.8 kgf threshold  \u2502  ",
                    "MAPE: ",
                    tags$span(style = "color: #27ae60;", "\u25a0"), " \u22645%   ",
                    tags$span(style = "color: #e67e22;", "\u25a0"), " 5\u201315%   ",
                    tags$span(style = "color: #c0392b;", "\u25a0"), " >15%")
                ),
                downloadButton("download_val_plot", "Download PNG",
                               class = "btn btn-sm btn-outline-secondary ms-3")
              )),
              plotOutput("val_plot")
            ),
            card(
              card_header(div(class = "d-flex justify-content-between align-items-center",
                span("Validation Metrics \u2014 per Hatch / Deck"),
                tags$small(class = "text-muted",
                           "% Error = (Sim \u2212 Obs) / Obs \u00d7 100.  ",
                           "Positive = model overpredicts firmness.")
              )),
              DTOutput("val_metrics")
            )
          )
        )
      )
    }
  ),

  # ── Predict ───────────────────────────────────────────────────
  nav_panel("Predict",
    if (!db_available) {
      div(class = "alert alert-warning m-4",
          tags$strong("Database packages unavailable."),
          " Install ", tags$code("DBI"), ", ", tags$code("odbc"),
          " and ODBC Driver 17 for SQL Server.")
    } else {
      div(class = "app-row",
        div(class = "app-sidebar",
          uiOutput("pred_db_status_ui"),
          hr(),
          span(class = "sec-lbl", "Filters"),
          uiOutput("pred_year_ui"),
          radioButtons("pred_location", "FMT Test Location",
                       choices = c("Shore", "Vessel"), selected = "Shore", inline = TRUE),
          uiOutput("pred_voyage_ui"),
          uiOutput("pred_hatch_ui"),
          hr(),
          span(class = "sec-lbl", "Model"),
          radioButtons("pred_variety", "Variety",
                       choices = c("Green (HW)" = "Green", "Gold (GA)" = "Gold"),
                       selected = "Green"),
          selectInput("pred_softening", "Softening type",
                      choices = c("fast softening", "average softening", "slow softening"),
                      selected = "average softening"),
          numericInput("pred_mc_n", "Replicates per hatch",
                       value = 40, min = 10, max = 2000, step = 10),
          hr(),
          span(class = "sec-lbl", "Ethylene application"),
          checkboxInput("pred_no_c2h4", "No ethylene (C\u2082H\u2084 = 0)", FALSE),
          div(class = "eth-box",
            conditionalPanel("input.pred_no_c2h4 == false",
              DTOutput("pred_eth_table", height = "auto"),
              div(class = "d-flex gap-1 mt-1",
                actionButton("pred_eth_add", "+ Add", class = "btn btn-outline-secondary btn-sm"),
                actionButton("pred_eth_rm",  "- Remove last", class = "btn btn-outline-secondary btn-sm")
              )
            ),
            conditionalPanel("input.pred_no_c2h4 == true",
              tags$p(class = "text-muted mb-0", style = "font-size: 0.78rem;",
                     "Ethylene concentration set to zero.")
            )
          ),
          hr(),
          span(class = "sec-lbl", "Future Temperature Extension"),
          helpText("Auto-filled from last known cargo setpoints.",
                   "Edit to model different future scenarios.",
                   "Day numbers are relative to each hatch's first FMT date."),
          DTOutput("pred_ext_table"),
          actionButton("pred_ext_add_row",
                       tagList(tags$i(class = "bi bi-plus-circle"), " Add Row"),
                       class = "btn btn-sm btn-outline-secondary mt-1 mb-2"),
          hr(),
          helpText(tags$i(class = "bi bi-info-circle me-1"),
                   "Initial firmness is auto-derived from each hatch's earliest FMT measurement."),
          actionButton("run_predict",
                       tagList(tags$i(class = "bi bi-arrow-right-circle"), "Run Prediction"),
                       class = "btn btn-success w-100")
        ),
        div(class = "app-main",
          div(class = "main-scroll",
            uiOutput("pred_landing_cards"),
            card(
              card_header("Temperature Profiles \u2014 Known Cargo Setpoints + Future Extension"),
              plotlyOutput("pred_temp_plot", height = "220px")
            ),
            card(full_screen = TRUE,
              card_header(div(class = "d-flex justify-content-between align-items-center",
                div(
                  tags$strong("Predicted Firmness \u2014 by Hatch/Deck"),
                  br(),
                  tags$small(class = "text-muted",
                    tags$span(style = "color: #6baed6;", "\u2014"), " Sim replicates  \u2502  ",
                    tags$span(style = "color: #2c3e50; font-weight: bold;", "\u2014"), " Sim mean  \u2502  ",
                    tags$span(style = "color: #e74c3c; font-weight: bold;", "\u25cf"),
                    " Observed FMT  \u2502  ",
                    tags$span(style = "color: #e74c3c;", "\u2014"), " Obs mean  \u2502  ",
                    tags$span(style = "color: #3498db; font-weight: bold;", "| today"),
                    "  \u2502  ",
                    tags$span(style = "color: #8e44ad; font-weight: bold;", "| landing"),
                    "  \u2502  ",
                    tags$span(style = "color: #e67e22;", "- - -"),
                    " 0.8 kgf")
                ),
                downloadButton("download_pred_plot", "Download PNG",
                               class = "btn btn-sm btn-outline-secondary ms-3")
              )),
              plotOutput("pred_plot")
            )
          )
        )
      )
    }
  ),
  nav_panel("Calibrate",
    if (!db_available) {
      div(class = "app-row", div(class = "app-main p-4",
        div(class = "alert alert-warning",
          tags$i(class = "bi bi-exclamation-triangle me-1"),
          "Database packages (DBI, odbc) are not available. Calibration requires EDW access.")))
    } else {
      div(class = "app-row",
        div(class = "app-sidebar",
          uiOutput("cal_db_status_ui"),
          hr(),
          span(class = "sec-lbl", "FILTERS"),
          uiOutput("cal_year_ui"),
          radioButtons("cal_variety", "Variety",
                       choices = c("Green (HW)" = "Green", "Gold (GA)" = "Gold"),
                       selected = "Gold", inline = TRUE),
          radioButtons("cal_location", "FMT Test Location",
                       choices = c("Shore" = "Shore", "Vessel" = "Vessel"),
                       selected = "Vessel", inline = TRUE),
          uiOutput("cal_voyage_ui"),
          uiOutput("cal_hatch_ui"),
          hr(),
          span(class = "sec-lbl", "CALIBRATION SETTINGS"),
          selectInput("cal_softening", "Softening type (for priors)",
                      choices = c("fast softening", "average softening", "slow softening"),
                      selected = "fast softening"),
          checkboxInput("cal_secondary", "Also calibrate rate constants (ksref, kpref, gamma)", FALSE),
          numericInput("cal_n_samples", "Posterior samples", value = 500, min = 100, max = 5000, step = 100),
          hr(),
          span(class = "sec-lbl", "TRAIN / VALIDATE SPLIT"),
          helpText("Calibrate on data up to cutoff day, then test prediction on later days."),
          sliderInput("cal_cutoff_day", "Train cutoff (day)", min = 1, max = 15, value = 5, step = 1),
          hr(),
          span(class = "sec-lbl", "ETHYLENE"),
          checkboxInput("cal_no_c2h4", "No ethylene (C\u2082H\u2084 = 0)", FALSE),
          div(class = "eth-box",
            conditionalPanel("input.cal_no_c2h4 == false",
              DTOutput("cal_eth_table", height = "auto"),
              div(class = "d-flex gap-1 mt-1",
                actionButton("cal_eth_add", "+ Add", class = "btn btn-outline-secondary btn-sm"),
                actionButton("cal_eth_rm", "- Remove last", class = "btn btn-outline-secondary btn-sm")
              )
            )
          ),
          hr(),
          actionButton("run_calibrate",
                        tagList(tags$i(class = "bi bi-gear"), " Run Calibration"),
                        class = "btn btn-primary w-100"),
          div(class = "mt-2",
            actionButton("cal_use_params",
                          tagList(tags$i(class = "bi bi-arrow-right-circle"), " Use for Prediction"),
                          class = "btn btn-outline-success btn-sm w-100"),
            downloadButton("cal_dl_params", "Download Parameters (CSV)",
                           class = "btn btn-outline-primary btn-sm w-100 mt-1")
          )
        ),
        div(class = "app-main",
          div(class = "main-scroll",
            uiOutput("cal_kpi_cards"),
            div(class = "card mb-3",
              div(class = "card-header fw-semibold",
                "Calibrated Fit vs Observed \u2014 per Hatch",
                tags$small(class = "text-muted ms-2",
                  "Grey region = training data | White region = validation (prediction test)")),
              div(class = "card-body p-2", plotOutput("cal_fit_plot", height = "auto"))
            ),
            div(class = "card mb-3",
              div(class = "card-header fw-semibold", "Validation Metrics \u2014 Held-out Data"),
              div(class = "card-body p-2", DTOutput("cal_val_metrics_table"))
            ),
            div(class = "card mb-3",
              div(class = "card-header fw-semibold", "Posterior Distributions"),
              div(class = "card-body p-2", plotlyOutput("cal_posterior_plot", height = "300px"))
            ),
            div(class = "card mb-3",
              div(class = "card-header fw-semibold", "Parameter Summary \u2014 per Hatch"),
              div(class = "card-body p-2", DTOutput("cal_param_table"))
            )
          )
        )
      )
    }
  ),
  nav_panel("FMT Data Viewer",
    div(class = "app-row",
      div(class = "app-sidebar",
        span(class = "sec-lbl", "Upload"),
        fileInput("fmt_files", "Choose files", accept = c(".json", ".xlsx", ".xls"),
                  multiple = TRUE),
        uiOutput("fmt_file_list"),
        actionButton("fmt_clear", "Clear all", class = "btn btn-outline-danger btn-sm w-100 mb-2"),
        uiOutput("json_table_selector"),
        hr(),
        span(class = "sec-lbl", "Target Firmness Range (kgf)"),
        div(class = "d-flex gap-2",
          numericInput("tfr_min", "Min", value = 2, min = 0, max = 12, step = 0.1, width = "50%"),
          numericInput("tfr_max", "Max", value = 2.5, min = 0, max = 12, step = 0.1, width = "50%")
        ),
        uiOutput("fmt_su_filter_ui"),
        hr(),
        span(class = "sec-lbl", "Download"),
        downloadButton("json_dl_csv",   "Download CSV",   class = "btn btn-outline-primary btn-sm w-100 mb-2"),
        downloadButton("json_dl_excel", "Download Excel", class = "btn btn-outline-primary btn-sm w-100")
      ),
      div(class = "app-main",
        tabsetPanel(id = "fmt_subtabs",
          tabPanel("Overview",
            uiOutput("json_summary"),
            uiOutput("json_meta_banner"),
            div(class = "d-flex gap-3", style = "min-height:0;",
              div(style = "flex:1; min-width:0;", uiOutput("json_firmness_panel")),
              div(style = "flex:1; min-width:0;", uiOutput("json_temp_panel"))
            ),
            DTOutput("json_dt")
          ),
          tabPanel("By Shipping Unit",
            uiOutput("fmt_su_panel")
          )
        )
      )
    )
  )
)

# =========================== Server ===========================
server <- function(input, output, session) {

  # ── Setup: treatment tables ───────────────────────────────────
  starter_profile <- function(last_day = 12, t1 = 10, t2 = 5, t3 = 0)
    tibble(day   = c(0, 5, 10, last_day),
           tempC = c(t1, t2, t3, t3))

  tr_tables <- reactiveValues(
    tr1 = starter_profile(),
    tr2 = starter_profile(last_day = 10, t1 = 15, t2 = 10, t3 = 2),
    tr3 = starter_profile(last_day = 7,  t1 = 5,  t2 = 5,  t3 = 5))

  proxy_update <- function(id, data)
    replaceData(dataTableProxy(id), data %>% mutate(`Delete` = NULL),
                resetPaging = FALSE, rownames = FALSE)

  render_tr_dt <- function(id, rv_name, del_input) {
    output[[id]] <- renderDT({
      df <- tr_tables[[rv_name]]
      if (nrow(df) < 2) { df <- tibble(day=c(0,1),tempC=c(0,0)); tr_tables[[rv_name]] <<- df }
      df %>%
        mutate(`Delete` = ifelse(dplyr::row_number() <= 2,
          '<button class="btn btn-sm btn-light" disabled><i class="bi bi-trash3"></i></button>',
          '<button class="btn btn-sm btn-link text-danger delete_btn"><i class="bi bi-trash3"></i></button>')) %>%
        datatable(editable = list(target = "cell"), rownames = FALSE,
                  selection = "none", escape = FALSE,
                  options = list(dom = "t", paging = FALSE, ordering = FALSE),
                  callback = JS(sprintf(
                    "table.on('click','button.delete_btn',function(){
   var idx=table.row($(this).closest('tr')).index()+1;
   Shiny.setInputValue('%s',idx,{priority:'event'});});", del_input)))
    })
  }
  render_tr_dt("tr1_table", "tr1", "tr1_delete_row")
  render_tr_dt("tr2_table", "tr2", "tr2_delete_row")
  render_tr_dt("tr3_table", "tr3", "tr3_delete_row")

  output$treat_tabs <- renderUI({
    make_tab <- function(lbl, name_id, val, dt_id, add_id)
      tabPanel(lbl, textInput(name_id, "Treatment name", value = val), DTOutput(dt_id),
               actionButton(add_id, tagList(tags$i(class = "bi bi-plus-circle"), " Add Row"),
                            class = "btn btn-sm btn-outline-secondary mt-2"))
    tabs <- list(make_tab("Treatment 1", "tr1_name", "T1", "tr1_table", "tr1_add_row"))
    if (input$treat_count >= 2) tabs[[2]] <- make_tab("Treatment 2","tr2_name","T2","tr2_table","tr2_add_row")
    if (input$treat_count >= 3) tabs[[3]] <- make_tab("Treatment 3","tr3_name","T3","tr3_table","tr3_add_row")
    do.call(tabsetPanel, c(list(id = "tr_tabs"), tabs))
  })

  for (tr_id in c("tr1", "tr2", "tr3")) {
    local({
      tid <- tr_id
      observeEvent(input[[paste0(tid, "_table_cell_edit")]], {
        info <- input[[paste0(tid, "_table_cell_edit")]]
        if (is.null(info$col)) return()
        col_idx  <- as.integer(info$col) + 1
        if (col_idx < 1 || col_idx > ncol(tr_tables[[tid]])) return()
        col_name <- colnames(tr_tables[[tid]])[col_idx]
        val      <- suppressWarnings(as.numeric(info$value)); if (is.na(val)) return()
        if (col_name %in% c("day", "tempC")) {
          tr_tables[[tid]][info$row, col_name] <- val
          proxy_update(paste0(tid, "_table"), tr_tables[[tid]])
        }
      })
      observeEvent(input[[paste0(tid, "_delete_row")]], {
        idx <- as.integer(input[[paste0(tid, "_delete_row")]])
        if (!is.na(idx) && nrow(tr_tables[[tid]]) > 2 && idx > 2) {
          tr_tables[[tid]] <- tr_tables[[tid]][-idx, , drop = FALSE]
          proxy_update(paste0(tid, "_table"), tr_tables[[tid]])
        }
      })
      observeEvent(input[[paste0(tid, "_add_row")]], {
        last <- tail(tr_tables[[tid]], 1)
        tr_tables[[tid]] <- bind_rows(tr_tables[[tid]],
          tibble(day = last$day + 1, tempC = last$tempC))
      })
    })
  }

  # ── Ethylene application tables (Setup / Validate / Predict) ──
  eth_default <- tibble(ppm = 100, day = 2, duration = 1)
  eth_rv      <- reactiveVal(eth_default)
  val_eth_rv  <- reactiveVal(eth_default)
  pred_eth_rv <- reactiveVal(eth_default)

  render_eth_dt <- function(df) {
    datatable(df, editable = TRUE, rownames = FALSE,
              options = list(dom = "t", ordering = FALSE, paging = FALSE,
                             scrollX = FALSE, autoWidth = FALSE,
                             columnDefs = list(list(className = "dt-center", targets = "_all"))),
              colnames = c("ppm", "Day", "Duration (d)"),
              class = "compact cell-border stripe") %>%
      formatRound(columns = c("ppm", "day", "duration"), digits = c(0, 1, 1))
  }

  output$eth_table     <- renderDT(render_eth_dt(eth_rv()))
  output$val_eth_table <- renderDT(render_eth_dt(val_eth_rv()))
  output$pred_eth_table <- renderDT(render_eth_dt(pred_eth_rv()))

  # Cell edits (with bounds check)
  edit_eth <- function(info, rv) {
    df <- rv()
    col_i <- as.integer(info$col) + 1L
    row_i <- as.integer(info$row)
    val <- suppressWarnings(as.numeric(info$value))
    if (!is.na(val) && col_i >= 1 && col_i <= ncol(df) && row_i >= 1 && row_i <= nrow(df)) {
      df[row_i, col_i] <- val; rv(df)
    }
  }
  observeEvent(input$eth_table_cell_edit,      edit_eth(input$eth_table_cell_edit, eth_rv))
  observeEvent(input$val_eth_table_cell_edit,  edit_eth(input$val_eth_table_cell_edit, val_eth_rv))
  observeEvent(input$pred_eth_table_cell_edit, edit_eth(input$pred_eth_table_cell_edit, pred_eth_rv))

  # Add / remove rows
  observeEvent(input$eth_add, {
    df <- eth_rv(); last <- tail(df, 1)
    eth_rv(bind_rows(df, tibble(ppm = last$ppm, day = last$day + 2, duration = last$duration)))
  })
  observeEvent(input$eth_rm, { df <- eth_rv(); if (nrow(df) > 1) eth_rv(df[-nrow(df), ]) })

  observeEvent(input$val_eth_add, {
    df <- val_eth_rv(); last <- tail(df, 1)
    val_eth_rv(bind_rows(df, tibble(ppm = last$ppm, day = last$day + 2, duration = last$duration)))
  })
  observeEvent(input$val_eth_rm, { df <- val_eth_rv(); if (nrow(df) > 1) val_eth_rv(df[-nrow(df), ]) })

  observeEvent(input$pred_eth_add, {
    df <- pred_eth_rv(); last <- tail(df, 1)
    pred_eth_rv(bind_rows(df, tibble(ppm = last$ppm, day = last$day + 2, duration = last$duration)))
  })
  observeEvent(input$pred_eth_rm, { df <- pred_eth_rv(); if (nrow(df) > 1) pred_eth_rv(df[-nrow(df), ]) })

  load_params <- reactive({
    validate(need(dir.exists(file.path(getwd(), "params")),
                  paste("Missing params folder:", file.path(getwd(), "params"))))
    tryCatch(load_param_bank(input$softening, input$variety),
             error = function(e) validate(need(FALSE, conditionMessage(e))))
  })

  sim_result <- eventReactive(input$run, {
    withProgress(message = "Simulating firmness...", value = 0, {
      params_all <- load_params()
      N <- as.integer(input$mc_n); validate(need(N > 0, "Replicates must be > 0"))
      N_eff <- min(N, nrow(params_all)); incProgress(0.05)
      idx <- sample(seq_len(nrow(params_all)), N_eff, replace = FALSE)
      params_df <- params_all[idx, c("E0","F0","Ffix1")] %>% as.data.frame()
      eth_df <- if (isTRUE(input$no_c2h4)) data.frame(ppm = 0, day = 0, duration = 0) else eth_rv()
      make_tr_profile <- function(temp_tbl) {
        max_day <- max(temp_tbl$day)
        bind_rows(temp_tbl %>% mutate(c2h4 = NA_real_),
                  build_eth_rows(eth_df, max_day)) %>% arrange(day)
      }
      tr_list <- list(list(name = input$tr1_name, profile_df = make_tr_profile(tr_tables$tr1)))
      if (input$treat_count >= 2) tr_list[[2]] <- list(name=input$tr2_name, profile_df=make_tr_profile(tr_tables$tr2))
      if (input$treat_count >= 3) tr_list[[3]] <- list(name=input$tr3_name, profile_df=make_tr_profile(tr_tables$tr3))
      incProgress(0.2)
      res <- run_montecarlo_firmness(params_df=params_df, treatments=tr_list, variety=input$variety,
               use_init=isTRUE(input$use_init), F_n=input$F_n, std_n=input$std_n, n_time_steps=51)
      incProgress(0.75); res
    })
  })

  observeEvent(sim_result(), { nav_select("main_nav", "Results") })

  output$mc_plot <- renderPlot({
    res <- sim_result(); req(res$trajectories)
    df <- res$trajectories
    avg_df <- df %>% group_by(treatment, time) %>% summarise(mean_f = mean(firmness), .groups = "drop")
    ggplot(df, aes(x = time, y = firmness, group = replicate)) +
      geom_line(alpha = 0.08, linewidth = 0.3, colour = "#6baed6") +
      geom_line(data = avg_df, aes(x = time, y = mean_f, group = NULL),
                linewidth = 1.3, colour = "#2c3e50") +
      facet_wrap(~ treatment) +
      ylim(0,10) +
      labs(x = "Time (days)", y = "Firmness (kgf)",
           subtitle = "Blue = individual replicates  |  Bold = mean firmness") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
  })

  output$summary_table <- renderDT({
    res <- sim_result()
    nice <- res$summary %>% dplyr::rename(
      `Treatment` = treatment, `Final day (d)` = final_day,
      `<= 0.8 kgf (%)` = pct_firmness_le_0_8,
      `Mean firmness at final day (kgf)` = mean_final_firmness,
      `SD at final day` = sd_final_firmness, `Replicates` = n)
    datatable(nice, rownames = FALSE, options = list(dom = "tip", pageLength = 5))
  })

  output$download_csv <- downloadHandler(
    filename = function() paste0("MC_Firmness_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(sim_result()$trajectories, file)
  )

  # ── Validate tab ─────────────────────────────────────────────
  if (db_available) {

    db_con_rv <- reactiveVal(NULL)
    db_status <- reactiveVal("idle")
    db_msg    <- reactiveVal("")

    observeEvent(input$db_connect, {
      notif_id <- showNotification(
        tagList(tags$strong("Connecting to EDW..."),
                " A browser window may open for Azure AD authentication."),
        duration = NULL, type = "message", closeButton = FALSE)
      tryCatch({
        con <- DBI::dbConnect(odbc::odbc(),
          UID = "abdulquadri.alaka@zespri.com",
          Driver = "ODBC Driver 17 for SQL Server",
          Server = "zapdw.database.windows.net", Database = "EDW",
          Authentication = "ActiveDirectoryInteractive")
        removeNotification(notif_id)
        db_con_rv(con); db_status("connected"); db_msg("Connected to EDW.")
        showNotification("Connected to EDW.", type = "message", duration = 3)
      }, error = function(e) {
        removeNotification(notif_id)
        db_con_rv(NULL); db_status("error"); db_msg(conditionMessage(e))
      })
    })

    output$db_status_ui <- renderUI({
      s <- db_status(); if (s == "idle") return(NULL)
      cls <- switch(s, connected = "alert alert-success py-2",
                       error     = "alert alert-danger py-2",
                       "alert alert-info py-2")
      div(class = paste(cls, "mb-2"), db_msg())
    })

    # ── Available seasons (from FMT data) ─────────────────────
    available_years <- reactive({
      req(db_con_rv())
      DBI::dbGetQuery(db_con_rv(), "
        SELECT DISTINCT year FROM [ods].[Coolchain_FMT_TestSheet_Materialise]
        WHERE year IS NOT NULL ORDER BY year DESC
      ") %>% pull(year)
    })

    output$val_year_ui <- renderUI({
      req(available_years())
      yrs <- available_years()
      cur <- as.integer(format(Sys.Date(), "%Y"))
      sel <- if (cur %in% yrs) cur else yrs[1]
      selectInput("val_year", "Season", choices = yrs, selected = sel)
    })

    output$pred_year_ui <- renderUI({
      req(available_years())
      yrs <- available_years()
      cur <- as.integer(format(Sys.Date(), "%Y"))
      sel <- if (cur %in% yrs) cur else yrs[1]
      selectInput("pred_year", "Season", choices = yrs, selected = sel)
    })

    # ── DB queries ────────────────────────────────────────────
    fmt_raw <- reactive({
      req(db_con_rv(), input$val_year)
      DBI::dbGetQuery(db_con_rv(), sprintf("
        SELECT voyage_name, variety, CAST(date AS DATE) AS date,
               fmt_testing_location, hatch_deck,
               CAST(firmness AS FLOAT) AS firmness,
               CAST(fruit_tempurature_in_hold AS FLOAT) AS fruit_tempurature_in_hold
        FROM [ods].[Coolchain_FMT_TestSheet_Materialise]
        WHERE year = %d AND firmness IS NOT NULL AND voyage_name IS NOT NULL
              AND CAST(firmness AS FLOAT) > 0 AND CAST(firmness AS FLOAT) <= 12
      ", as.integer(input$val_year))) %>%
        mutate(hatch_deck = str_replace_all(hatch_deck, "\\s+", ""), date = as.Date(date))
    })

    cargo_raw <- reactive({
      req(db_con_rv(), input$val_year)
      DBI::dbGetQuery(db_con_rv(), sprintf("
        SELECT Voyage_Name, Hatch_Code, Set_Point_Temperature_Celcius,
               CAST(CargoConditionReportDate AS DATE) AS CargoConditionReportDate
        FROM [ods].[Coolchain_CargoCondition_Report_Materialise]
        WHERE Voyage_Name IN (
          SELECT DISTINCT voyage_name
          FROM [ods].[Coolchain_FMT_TestSheet_Materialise]
          WHERE year = %d AND voyage_name IS NOT NULL
        )
      ", as.integer(input$val_year))) %>%
        mutate(Hatch_Code = str_replace_all(Hatch_Code, "\\s+", ""),
               CargoConditionReportDate = as.Date(CargoConditionReportDate))
    })

    # ── Selectors ─────────────────────────────────────────────
    output$val_voyage_ui <- renderUI({
      req(fmt_raw())
      voyages <- fmt_raw() %>%
        filter(variety == variety_db(input$val_variety),
               fmt_testing_location == input$val_location) %>%
        distinct(voyage_name) %>% arrange(voyage_name) %>% pull()
      if (!length(voyages)) return(helpText("No voyages found for selected variety/year/location."))
      selectInput("val_voyage", "Voyage", choices = voyages)
    })

    output$val_hatch_ui <- renderUI({
      req(fmt_raw(), input$val_voyage)
      hatches <- fmt_raw() %>%
        filter(voyage_name == input$val_voyage,
               variety == variety_db(input$val_variety),
               fmt_testing_location == input$val_location) %>%
        distinct(hatch_deck) %>% arrange(hatch_deck) %>% pull()
      if (!length(hatches)) return(helpText("No hatch/deck data for this voyage/location."))
      checkboxGroupInput("val_hatches", "Hatch / Deck", choices = hatches, selected = hatches)
    })

    # ── Per-hatch profiles ────────────────────────────────────
    val_profile_data <- reactive({
      req(fmt_raw(), cargo_raw(), input$val_voyage, input$val_hatches)
      var_code <- variety_db(input$val_variety)
      eth_df <- if (isTRUE(input$val_no_c2h4)) data.frame(ppm = 0, day = 0, duration = 0) else val_eth_rv()

      hatch_profiles <- list(); skipped <- character(0)
      for (hatch in input$val_hatches) {
        hatch_fmt <- fmt_raw() %>%
          filter(voyage_name == input$val_voyage, variety == var_code,
                 hatch_deck == hatch, fmt_testing_location == input$val_location,
                 !is.na(firmness))
        if (!nrow(hatch_fmt)) { skipped <- c(skipped, paste0(hatch," (no FMT)")); next }
        day0 <- min(hatch_fmt$date, na.rm = TRUE)
        # Full cargo range relative to first FMT date (may include negative days)
        temp_all <- cargo_raw() %>%
          filter(Voyage_Name == input$val_voyage, Hatch_Code == hatch,
                 !is.na(Set_Point_Temperature_Celcius)) %>%
          mutate(day = as.numeric(difftime(CargoConditionReportDate, day0, units = "days"))) %>%
          group_by(day) %>%
          summarise(tempC = mean(Set_Point_Temperature_Celcius, na.rm = TRUE), .groups = "drop") %>%
          arrange(day)
        if (!nrow(temp_all)) { skipped <- c(skipped, paste0(hatch, " (no cargo data)")); next }
        # Simulation profile: anchor day 0 to the setpoint in effect at first FMT date
        pre_fmt <- temp_all %>% filter(day <= 0) %>% tail(1)
        temp_df <- bind_rows(
          if (nrow(pre_fmt)) tibble(day = 0, tempC = pre_fmt$tempC) else tibble(),
          temp_all %>% filter(day > 0)
        ) %>% arrange(day) %>% mutate(c2h4 = NA_real_)
        if (!nrow(temp_df)) { skipped <- c(skipped, paste0(hatch, " (no cargo data)")); next }
        max_day  <- max(temp_df$day)
        hatch_profiles[[hatch]] <- list(
          day0       = day0,
          temp_all   = temp_all,   # full cargo range — used by preview plot
          profile_df = bind_rows(temp_df, build_eth_rows(eth_df, max_day)) %>% arrange(day))
      }
      if (length(hatch_profiles) == 0) {
        skip_msg <- if (length(skipped) > 0) paste0(" Skipped: ", paste(skipped, collapse = "; ")) else ""
        validate(paste0("No data for selected hatches.", skip_msg))
      }
      list(hatch_profiles = hatch_profiles, skipped = skipped)
    })

    # Temperature profile preview
    output$val_temp_plot <- renderPlotly({
      req(val_profile_data())
      pd <- val_profile_data()
      eth <- if (isTRUE(input$val_no_c2h4)) NULL else val_eth_rv()

      temp_df <- purrr::map_dfr(names(pd$hatch_profiles), function(h) {
        tf      <- pd$hatch_profiles[[h]]$temp_all
        day_min <- min(tf$day, na.rm = TRUE)
        tf %>% mutate(Hatch = h, day = day - day_min,
                      tip = paste0("Hatch: ", h, "<br>Day: ", round(day - day_min, 1), "<br>Temp: ", round(tempC, 1), "\u00b0C"))
      })

      # Build plotly directly for better tooltip control
      p <- plot_ly()
      hatches <- unique(temp_df$Hatch)
      colours <- RColorBrewer::brewer.pal(max(3, length(hatches)), "Set1")
      for (i in seq_along(hatches)) {
        hd <- temp_df %>% filter(Hatch == hatches[i])
        p <- p %>%
          add_trace(data = hd, x = ~day, y = ~tempC, type = "scatter",
                    mode = "lines+markers", line = list(shape = "hv", width = 2, color = colours[i]),
                    marker = list(size = 6, color = colours[i]),
                    name = hatches[i], text = ~tip, hoverinfo = "text")
      }
      # Ethylene markers
      if (!is.null(eth) && nrow(eth) > 0 && any(eth$ppm > 0)) {
        y_base <- min(temp_df$tempC, na.rm = TRUE) - 0.3
        em <- eth %>% filter(ppm > 0) %>%
          mutate(dur_label = ifelse(duration == 0.5, "12 h", ifelse(duration == 1, "24 h", paste0(duration, " d"))),
                 tip = paste0("\u2b25 Ethylene<br>", ppm, " ppm<br>Day ", day, "<br>Duration: ", dur_label),
                 y = y_base)
        p <- p %>%
          add_trace(data = em, x = ~day, y = ~y, type = "scatter", mode = "markers",
                    marker = list(size = 12, color = "#8e44ad", symbol = "diamond"),
                    name = "Ethylene", text = ~tip, hoverinfo = "text",
                    showlegend = TRUE)
      }
      p %>% layout(xaxis = list(title = "Days from first cargo report"),
                   yaxis = list(title = "Set Point (\u00b0C)"),
                   margin = list(t = 10, b = 40, l = 50, r = 10),
                   legend = list(orientation = "v", x = 1.02, y = 1))
    })

    # ── Simulation: one run per hatch ─────────────────────────
    val_sim <- eventReactive(input$run_validate, {
      req(val_profile_data())
      var_code <- variety_db(input$val_variety)
      pd       <- val_profile_data()
      hatch_names <- names(pd$hatch_profiles)

      withProgress(message = "Running validation...", value = 0, {
        params_all <- tryCatch(load_param_bank(input$val_softening, input$val_variety),
                               error = function(e) validate(need(FALSE, conditionMessage(e))))
        N_eff <- min(as.integer(input$val_mc_n), nrow(params_all))
        all_trajs <- list(); all_summaries <- list(); obs_list <- list()

        for (i in seq_along(hatch_names)) {
          hatch <- hatch_names[i]; hp <- pd$hatch_profiles[[hatch]]
          setProgress(i / length(hatch_names) * 0.9,
                      message = paste0("Hatch ", hatch, " (", i, "/", length(hatch_names), ")..."))
          hatch_fmt <- fmt_raw() %>%
            filter(voyage_name == input$val_voyage, variety == var_code,
                   hatch_deck == hatch, fmt_testing_location == input$val_location,
                   !is.na(firmness))
          init_f <- hatch_fmt %>% filter(date == hp$day0) %>% pull(firmness)
          F_n_h  <- if (length(init_f) && !all(is.na(init_f))) mean(init_f, na.rm=TRUE) else 7
          sd_n_h <- if (length(init_f) > 1) max(sd(init_f, na.rm=TRUE), 0.1) else 1
          idx    <- sample(seq_len(nrow(params_all)), N_eff, replace = FALSE)
          params_df <- params_all[idx, c("E0","F0","Ffix1")] %>% as.data.frame()
          res <- run_montecarlo_firmness(
            params_df = params_df,
            treatments = list(list(name = hatch, profile_df = hp$profile_df)),
            variety = input$val_variety, use_init = TRUE,
            F_n = F_n_h, std_n = sd_n_h, n_time_steps = 51)
          all_trajs[[hatch]]     <- res$trajectories
          all_summaries[[hatch]] <- res$summary
          obs_list[[hatch]] <- hatch_fmt %>%
            mutate(day = as.numeric(difftime(date, hp$day0, units = "days")),
                   treatment = hatch) %>%
            filter(day >= 0) %>%
            select(treatment, hatch_deck, date, day, fmt_testing_location, firmness)
        }
        list(sim = list(trajectories = bind_rows(all_trajs),
                        summary      = bind_rows(all_summaries)),
             obs            = bind_rows(obs_list),
             hatch_profiles = pd$hatch_profiles,
             skipped        = pd$skipped)
      })
    })

    # ── Metrics data (shared by plot labels, DT, and KPI cards) ──
    val_metrics_data <- reactive({
      req(val_sim())
      v <- val_sim(); sim_df <- v$sim$trajectories; obs_df <- v$obs
      if (!nrow(obs_df)) return(NULL)
      obs_agg <- obs_df %>%
        group_by(treatment, date, day) %>%
        summarise(obs_mean       = mean(firmness, na.rm = TRUE),
                  obs_sd         = sd(firmness,   na.rm = TRUE),
                  n_obs          = n(),
                  obs_pct_le_0_8 = mean(firmness <= 0.8, na.rm = TRUE) * 100,
                  .groups = "drop")
      obs_agg %>%
        mutate(sim_vals = purrr::map2(treatment, day, function(tr, d) {
          sub <- sim_df %>% filter(treatment == tr)
          ct  <- sub$time[which.min(abs(sub$time - d))]
          sub %>% filter(time == ct) %>%
            summarise(sim_mean = mean(firmness), sim_pct_le_0_8 = mean(firmness <= 0.8)*100,
                      .groups = "drop")
        })) %>%
        tidyr::unnest(sim_vals) %>%
        mutate(diff_kgf  = round(obs_mean - sim_mean, 2),
               pct_error = round((sim_mean - obs_mean) / obs_mean * 100, 1)) %>%
        arrange(treatment, date)
    })

    # Per-hatch MAPE summary (used by plot labels and KPI cards)
    val_hatch_mape <- reactive({
      req(val_metrics_data())
      val_metrics_data() %>%
        group_by(treatment) %>%
        summarise(mape   = round(mean(abs(pct_error), na.rm = TRUE), 1),
                  .groups = "drop") %>%
        mutate(label  = paste0("MAPE: ", mape, "%"),
               colour = mape_colour(mape))
    })

    # ── KPI cards ─────────────────────────────────────────────
    output$val_kpi_cards <- renderUI({
      req(val_hatch_mape())
      hm       <- val_hatch_mape()
      avg_mape <- round(mean(hm$mape), 1)
      best     <- hm %>% dplyr::slice_min(mape, n = 1, with_ties = FALSE)
      worst    <- hm %>% dplyr::slice_max(mape, n = 1, with_ties = FALSE)
      n_pass   <- sum(hm$mape <= 10)
      n_total  <- nrow(hm)

      kpi_card <- function(title, big_value, sub = NULL, bg = "#2c3e50") {
        div(class = "card kpi-card h-100 border-0",
            style = paste0("background:", bg, " !important; color:white; border-radius:var(--radius-lg);"),
            div(class = "card-body p-3",
                tags$p(class = "mb-1", style = "font-size:0.7rem; font-weight:600; text-transform:uppercase; letter-spacing:0.06em; opacity:0.8;", title),
                tags$p(class = "mb-0", style = "font-size:1.75rem; font-weight:700; line-height:1;", big_value),
                if (!is.null(sub)) tags$p(class = "mb-0 mt-1", style = "font-size:0.75rem; opacity:0.75;", sub)))
      }

      bg_avg  <- if (avg_mape <= 5) "#27ae60" else if (avg_mape <= 15) "#e67e22" else "#c0392b"
      bg_pass <- if (n_pass == n_total) "#27ae60" else if (n_pass >= n_total*0.75) "#e67e22" else "#c0392b"

      tagList(
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          gap = "0.75rem",
          kpi_card("Average MAPE", paste0(avg_mape, "%"),
                   sub = paste0(n_total, " hatches"), bg = bg_avg),
          kpi_card("Best Hatch",   paste0(best$mape, "%"),
                   sub = best$treatment,          bg = "#27ae60"),
          kpi_card("Worst Hatch",  paste0(worst$mape, "%"),
                   sub = worst$treatment,         bg = mape_colour(worst$mape)),
          kpi_card(paste0("Within \u226410% MAPE"),
                   paste0(n_pass, " / ", n_total),
                   sub = "hatches",               bg = bg_pass)
        ),
        tags$div(style = "height: 0.75rem;")
      )
    })

    # ── Validation plot (reactive object — reused for download) ──
    val_plot_gg <- reactive({
      req(val_sim(), val_hatch_mape())
      v      <- val_sim()
      sim_df <- v$sim$trajectories
      obs_df <- v$obs
      mape_df <- val_hatch_mape()
      avg_df <- sim_df %>% group_by(treatment, time) %>%
        summarise(mean_f = mean(firmness), .groups = "drop")
      n <- length(unique(sim_df$treatment))

      eth_sub <- if (isTRUE(input$val_no_c2h4)) {
        "  \u2022  No ethylene"
      } else {
        edf <- val_eth_rv()
        parts <- paste0(edf$ppm, " ppm d", edf$day, " \u00d7", edf$duration, "d")
        paste0("  \u2022  C\u2082H\u2084: ", paste(parts, collapse = "; "))
      }

      # Observed average line per hatch
      obs_avg_df <- if (nrow(obs_df) > 0) {
        obs_df %>% group_by(treatment, day) %>%
          summarise(mean_f = mean(firmness, na.rm = TRUE), .groups = "drop")
      } else { tibble(treatment = character(), day = numeric(), mean_f = numeric()) }

      p <- ggplot(sim_df, aes(x = time, y = firmness, group = replicate)) +
        geom_line(alpha = 0.08, linewidth = 0.3, colour = "#6baed6") +
        geom_hline(yintercept = 0.8, linetype = "dashed",
                   colour = "#e67e22", linewidth = 0.55) +
        geom_line(data = avg_df, aes(x = time, y = mean_f, group = NULL),
                  linewidth = 1.2, colour = "#2c3e50") +
        # MAPE annotation — coloured by severity, top-right of each facet
        geom_label(data = mape_df,
                   aes(x = Inf, y = Inf, label = label, colour = I(colour)),
                   hjust = 1.05, vjust = 1.3, size = 3.1, fontface = "bold",
                   fill = alpha("white", 0.82), label.size = NA,
                   label.padding = unit(0.2, "lines"), inherit.aes = FALSE) +
        facet_wrap(~ treatment, ncol = 4, scales = "free_x") +
        ylim(0,10) +
        labs(x = "Days from first FMT test (per hatch)",
             y = "Firmness (kgf)",
             title = paste("Voyage:", input$val_voyage),
             subtitle = paste0(input$val_variety, "  \u2022  ", input$val_softening, eth_sub)) +
        theme_minimal(base_size = 12) +
        theme(panel.grid.minor  = element_blank(),
              strip.text        = element_text(face = "bold"),
              plot.title        = element_text(face = "bold", size = 13),
              plot.subtitle     = element_text(colour = "grey40", size = 10),
              legend.position   = "bottom")

      if (nrow(obs_df) > 0) {
        p <- p +
          geom_point(data = obs_df,
                     aes(x = day, y = firmness, shape = fmt_testing_location),
                     colour = "#e74c3c", size = 3, stroke = 0.9,
                     inherit.aes = FALSE) +
          geom_line(data = obs_avg_df,
                    aes(x = day, y = mean_f, group = NULL),
                    linewidth = 1, colour = "#e74c3c", linetype = "solid",
                    inherit.aes = FALSE) +
          scale_shape_manual(values = c("Shore" = 16, "Vessel" = 17),
                             name = "FMT Location")
      }
      p
    })

    output$val_plot <- renderPlot({
      val_plot_gg()
    }, height = function() {
      req(val_sim())
      n <- length(unique(val_sim()$sim$trajectories$treatment))
      max(300, ceiling(n / 4) * 270)
    })

    # Download validation plot as high-res PNG
    output$download_val_plot <- downloadHandler(
      filename = function()
        paste0("validation_", input$val_voyage, "_", Sys.Date(), ".png"),
      content = function(file) {
        req(val_plot_gg())
        n     <- length(unique(val_sim()$sim$trajectories$treatment))
        ncols <- min(4, n)
        nrows <- ceiling(n / ncols)
        ggplot2::ggsave(file, plot = val_plot_gg(),
                        width  = ncols * 5,
                        height = nrows * 4,
                        dpi    = 180, bg = "white")
      }
    )

    # ── Metrics DT ───────────────────────────────────────────
    output$val_metrics <- renderDT({
      req(val_metrics_data(), val_hatch_mape())
      m    <- val_metrics_data()
      mape <- val_hatch_mape()

      mape_txt <- paste(paste0(mape$treatment, ": ", mape$mape, "%"), collapse = "   \u2502   ")
      avg_mape <- round(mean(mape$mape), 1)

      m %>%
        mutate(across(c(obs_mean, sim_mean, obs_sd), ~ round(.x, 2)),
               across(c(obs_pct_le_0_8, sim_pct_le_0_8), ~ round(.x, 1)),
               day = round(day, 1)) %>%
        dplyr::select(treatment, date, day, n_obs,
                      obs_mean, sim_mean, diff_kgf, pct_error,
                      obs_pct_le_0_8, sim_pct_le_0_8) %>%
        dplyr::rename(
          `Hatch`        = treatment,
          `Date`         = date,
          `Day`          = day,
          `N obs`        = n_obs,
          `Obs (kgf)`    = obs_mean,
          `Sim (kgf)`    = sim_mean,
          `Diff (kgf)`   = diff_kgf,
          `% Error`      = pct_error,
          `Obs ≤0.8 (%)` = obs_pct_le_0_8,
          `Sim ≤0.8 (%)` = sim_pct_le_0_8
        ) %>%
        datatable(
          rownames = FALSE,
          caption  = htmltools::tags$caption(
            style = "caption-side:bottom; text-align:left; color:#555; padding-top:8px;",
            tags$strong(paste0("Avg MAPE: ", avg_mape, "%   \u2502   ")),
            mape_txt),
          options = list(
            dom        = "tip",
            pageLength = 20,
            columnDefs = list(list(className = "dt-center", targets = "_all")))
        ) %>%
        DT::formatStyle("% Error",
          color      = DT::styleInterval(c(-15, -5, 5, 15),
                         c("#c0392b","#e67e22","#27ae60","#e67e22","#c0392b")),
          fontWeight = "bold") %>%
        DT::formatStyle(
          c("Obs \u22640.8 (%)", "Sim \u22640.8 (%)"),
          background = DT::styleColorBar(c(0, 100), "#e8f4fd"),
          backgroundSize = "100% 70%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center")
    })

    # Shared reactive: calibrated parameters from Calibrate tab
    calibrated_params_rv <- reactiveVal(NULL)

    # ── Predict tab ─────────────────────────────────────────────

    output$pred_db_status_ui <- renderUI({
      s <- db_status(); if (s == "idle") return(
        div(class = "alert alert-info py-2 mb-2",
            tags$i(class = "bi bi-database me-1"),
            "Connect to the database in the ", tags$strong("Validate"), " tab first."))
      cls <- switch(s, connected = "alert alert-success py-2",
                       error     = "alert alert-danger py-2",
                       "alert alert-info py-2")
      div(class = paste(cls, "mb-2"), db_msg())
    })

    pred_fmt_raw <- reactive({
      req(db_con_rv(), input$pred_year)
      DBI::dbGetQuery(db_con_rv(), sprintf("
        SELECT voyage_name, variety, CAST(date AS DATE) AS date,
               fmt_testing_location, hatch_deck,
               CAST(firmness AS FLOAT) AS firmness,
               CAST(fruit_tempurature_in_hold AS FLOAT) AS fruit_tempurature_in_hold
        FROM [ods].[Coolchain_FMT_TestSheet_Materialise]
        WHERE year = %d AND firmness IS NOT NULL AND voyage_name IS NOT NULL
              AND CAST(firmness AS FLOAT) > 0 AND CAST(firmness AS FLOAT) <= 12
      ", as.integer(input$pred_year))) %>%
        mutate(hatch_deck = str_replace_all(hatch_deck, "\\s+", ""), date = as.Date(date))
    })

    pred_cargo_raw <- reactive({
      req(db_con_rv(), input$pred_year)
      DBI::dbGetQuery(db_con_rv(), sprintf("
        SELECT Voyage_Name, Hatch_Code, Set_Point_Temperature_Celcius,
               CAST(CargoConditionReportDate AS DATE) AS CargoConditionReportDate
        FROM [ods].[Coolchain_CargoCondition_Report_Materialise]
        WHERE Voyage_Name IN (
          SELECT DISTINCT voyage_name
          FROM [ods].[Coolchain_FMT_TestSheet_Materialise]
          WHERE year = %d AND voyage_name IS NOT NULL
        )
      ", as.integer(input$pred_year))) %>%
        mutate(Hatch_Code = str_replace_all(Hatch_Code, "\\s+", ""),
               CargoConditionReportDate = as.Date(CargoConditionReportDate))
    })

    pred_voyage_status <- reactive({
      req(db_con_rv(), input$pred_year)
      DBI::dbGetQuery(db_con_rv(), sprintf("
        SELECT Voyage_Name, [Shipment Status] AS shipment_status
        FROM [lab_Technical].[Coolchain_Voyage_View]
        WHERE Season = %d AND Fruit_Monitored = 'Fruit Monitored'
      ", as.integer(input$pred_year)))
    })

    output$pred_voyage_ui <- renderUI({
      req(pred_fmt_raw(), pred_voyage_status())
      on_route <- pred_voyage_status() %>%
        filter(shipment_status == "On Route") %>% pull(Voyage_Name)
      voyages <- pred_fmt_raw() %>%
        filter(variety == variety_db(input$pred_variety),
               fmt_testing_location == input$pred_location,
               voyage_name %in% on_route) %>%
        distinct(voyage_name) %>% arrange(voyage_name) %>% pull()
      if (!length(voyages)) return(helpText("No On Route voyages found for selected filters."))
      selectInput("pred_voyage", "Voyage (On Route)", choices = voyages)
    })

    output$pred_hatch_ui <- renderUI({
      req(pred_fmt_raw(), input$pred_voyage)
      hatches <- pred_fmt_raw() %>%
        filter(voyage_name == input$pred_voyage,
               variety == variety_db(input$pred_variety),
               fmt_testing_location == input$pred_location) %>%
        distinct(hatch_deck) %>% arrange(hatch_deck) %>% pull()
      if (!length(hatches)) return(helpText("No hatch/deck data for this voyage/location."))
      checkboxGroupInput("pred_hatches", "Hatch / Deck", choices = hatches, selected = hatches)
    })

    # ── Editable future extension table ───────────────────────
    pred_ext_rv <- reactiveValues(table = tibble(day = c(0, 20), tempC = c(0, 0)))

    # Auto-populate extension table when voyage/hatch selection changes
    observeEvent(list(input$pred_voyage, input$pred_hatches), {
      req(pred_cargo_raw(), pred_fmt_raw(), input$pred_voyage, input$pred_hatches)
      var_code <- variety_db(input$pred_variety)
      # Find day0 per hatch (earliest FMT date) -> use min across hatches as reference
      fmt_sub <- pred_fmt_raw() %>%
        filter(voyage_name == input$pred_voyage, variety == var_code,
               hatch_deck %in% input$pred_hatches,
               fmt_testing_location == input$pred_location, !is.na(firmness))
      if (!nrow(fmt_sub)) return()
      overall_day0 <- min(fmt_sub$date, na.rm = TRUE)
      # Last cargo report per hatch
      cargo_sub <- pred_cargo_raw() %>%
        filter(Voyage_Name == input$pred_voyage, Hatch_Code %in% input$pred_hatches,
               !is.na(Set_Point_Temperature_Celcius)) %>%
        group_by(Hatch_Code) %>%
        filter(CargoConditionReportDate == max(CargoConditionReportDate)) %>%
        summarise(last_day  = as.numeric(difftime(first(CargoConditionReportDate), overall_day0, units="days")),
                  last_temp = round(mean(Set_Point_Temperature_Celcius, na.rm=TRUE), 1),
                  .groups = "drop")
      if (!nrow(cargo_sub)) return()
      ext_start <- max(0, floor(min(cargo_sub$last_day, na.rm=TRUE)))
      ext_temp  <- round(median(cargo_sub$last_temp, na.rm=TRUE), 1)
      ext_end   <- ext_start + 10   # default 10 days ahead — user can edit
      pred_ext_rv$table <- tibble(day = c(ext_start, ext_end), tempC = c(ext_temp, ext_temp))
    }, ignoreInit = TRUE)

    output$pred_ext_table <- renderDT({
      df <- pred_ext_rv$table
      df %>%
        mutate(`Delete` = ifelse(dplyr::row_number() <= 2,
          '<button class="btn btn-sm btn-light" disabled><i class="bi bi-trash3"></i></button>',
          '<button class="btn btn-sm btn-link text-danger delete_btn"><i class="bi bi-trash3"></i></button>')) %>%
        datatable(editable = list(target = "cell"), rownames = FALSE,
                  selection = "none", escape = FALSE,
                  options = list(dom = "t", paging = FALSE, ordering = FALSE),
                  callback = JS("table.on('click','button.delete_btn',function(){
  var idx=table.row($(this).closest('tr')).index()+1;
  Shiny.setInputValue('pred_ext_delete_row',idx,{priority:'event'});});"))
    })

    observeEvent(input$pred_ext_table_cell_edit, {
      info <- input$pred_ext_table_cell_edit
      if (is.null(info$col)) return()
      col_idx  <- as.integer(info$col) + 1
      if (col_idx < 1 || col_idx > ncol(pred_ext_rv$table)) return()
      col_name <- colnames(pred_ext_rv$table)[col_idx]
      val <- suppressWarnings(as.numeric(info$value)); if (is.na(val)) return()
      if (col_name %in% c("day", "tempC")) {
        pred_ext_rv$table[info$row, col_name] <- val
        replaceData(dataTableProxy("pred_ext_table"),
                    pred_ext_rv$table %>% mutate(`Delete` = NULL),
                    resetPaging = FALSE, rownames = FALSE)
      }
    })

    observeEvent(input$pred_ext_delete_row, {
      idx <- as.integer(input$pred_ext_delete_row)
      if (!is.na(idx) && nrow(pred_ext_rv$table) > 2 && idx > 2) {
        pred_ext_rv$table <- pred_ext_rv$table[-idx, , drop = FALSE]
        replaceData(dataTableProxy("pred_ext_table"),
                    pred_ext_rv$table %>% mutate(`Delete` = NULL),
                    resetPaging = FALSE, rownames = FALSE)
      }
    })

    observeEvent(input$pred_ext_add_row, {
      last <- tail(pred_ext_rv$table, 1)
      pred_ext_rv$table <- bind_rows(pred_ext_rv$table,
        tibble(day = last$day + 1, tempC = last$tempC))
    })

    # ── Per-hatch prediction profiles ────────────────────────────
    pred_profile_data <- reactive({
      req(pred_fmt_raw(), pred_cargo_raw(), input$pred_voyage, input$pred_hatches)
      var_code <- variety_db(input$pred_variety)
      ext_table <- pred_ext_rv$table
      eth_df <- if (isTRUE(input$pred_no_c2h4)) data.frame(ppm = 0, day = 0, duration = 0) else pred_eth_rv()

      hatch_profiles <- list(); skipped <- character(0)
      for (hatch in input$pred_hatches) {
        hatch_fmt <- pred_fmt_raw() %>%
          filter(voyage_name == input$pred_voyage, variety == var_code,
                 hatch_deck == hatch, fmt_testing_location == input$pred_location,
                 !is.na(firmness))
        if (!nrow(hatch_fmt)) { skipped <- c(skipped, paste0(hatch, " (no FMT)")); next }
        day0 <- min(hatch_fmt$date, na.rm = TRUE)

        # Known cargo setpoints (full history)
        temp_all <- pred_cargo_raw() %>%
          filter(Voyage_Name == input$pred_voyage, Hatch_Code == hatch,
                 !is.na(Set_Point_Temperature_Celcius)) %>%
          mutate(day = as.numeric(difftime(CargoConditionReportDate, day0, units = "days"))) %>%
          group_by(day) %>%
          summarise(tempC = mean(Set_Point_Temperature_Celcius, na.rm=TRUE), .groups="drop") %>%
          arrange(day)
        if (!nrow(temp_all)) { skipped <- c(skipped, paste0(hatch, " (no cargo data)")); next }

        last_known_day  <- max(temp_all$day)
        last_known_temp <- temp_all %>% filter(day == max(day)) %>% pull(tempC) %>% mean()
        today_day <- last_known_day   # "today" = most recent cargo report

        # Future extension rows from user table (days beyond last known)
        ext_rows <- ext_table %>%
          filter(day > last_known_day) %>%
          mutate(c2h4 = NA_real_)
        # If extension doesn't cover beyond last_known_day, extend last setpoint
        if (!nrow(ext_rows) && nrow(ext_table) > 0) {
          ext_end <- max(ext_table$day)
          if (ext_end > last_known_day)
            ext_rows <- tibble(day = ext_end, tempC = last_known_temp, c2h4 = NA_real_)
        }
        landing_day <- if (nrow(ext_rows) > 0) max(ext_rows$day) else last_known_day

        # Simulation profile: anchor day0, include pre-FMT cargo as day-0 setpoint
        pre_fmt <- temp_all %>% filter(day <= 0) %>% tail(1)
        temp_df <- bind_rows(
          if (nrow(pre_fmt)) tibble(day = 0, tempC = pre_fmt$tempC, c2h4 = NA_real_) else tibble(),
          temp_all %>% filter(day > 0) %>% mutate(c2h4 = NA_real_),
          ext_rows
        ) %>% arrange(day)
        if (!nrow(temp_df) || max(temp_df$day) <= 0) {
          skipped <- c(skipped, paste0(hatch, " (insufficient temp data)")); next
        }

        max_day <- max(temp_df$day)
        profile_df <- bind_rows(temp_df,
          build_eth_rows(eth_df, max_day)) %>% arrange(day)

        hatch_profiles[[hatch]] <- list(
          day0        = day0,
          temp_all    = temp_all,
          today_day   = today_day,
          landing_day = landing_day,
          profile_df  = profile_df)
      }
      if (length(hatch_profiles) == 0) {
        skip_msg <- if (length(skipped) > 0) paste0(" Skipped: ", paste(skipped, collapse = "; ")) else ""
        validate(paste0("No data for selected hatches.", skip_msg))
      }
      list(hatch_profiles = hatch_profiles, skipped = skipped)
    })

    # Temperature profile preview (known history + future extension + ethylene markers)
    output$pred_temp_plot <- renderPlotly({
      req(pred_profile_data())
      pd <- pred_profile_data()
      ext <- pred_ext_rv$table
      eth <- if (isTRUE(input$pred_no_c2h4)) NULL else pred_eth_rv()

      known_df <- purrr::map_dfr(names(pd$hatch_profiles), function(h) {
        tf      <- pd$hatch_profiles[[h]]$temp_all
        day_min <- min(tf$day, na.rm=TRUE)
        tf %>% mutate(Hatch = h, day = day - day_min, type = "Known")
      })
      ext_df <- purrr::map_dfr(names(pd$hatch_profiles), function(h) {
        hp      <- pd$hatch_profiles[[h]]
        day_min <- min(hp$temp_all$day, na.rm=TRUE)
        last_pt <- hp$temp_all %>% filter(day == max(day))
        ext %>% filter(day > hp$today_day) %>%
          bind_rows(tibble(day = hp$today_day, tempC = last_pt$tempC[1]), .) %>%
          mutate(Hatch = h, day = day - day_min, type = "Future (extension)")
      })
      all_df <- bind_rows(known_df, ext_df)

      p <- plot_ly()
      hatches <- unique(all_df$Hatch)
      colours <- RColorBrewer::brewer.pal(max(3, length(hatches)), "Set1")
      for (i in seq_along(hatches)) {
        kd <- all_df %>% filter(Hatch == hatches[i], type == "Known") %>%
          mutate(tip = paste0("Hatch: ", Hatch, "<br>Day: ", round(day, 1), "<br>Temp: ", round(tempC, 1), "\u00b0C"))
        ed <- all_df %>% filter(Hatch == hatches[i], type == "Future (extension)") %>%
          mutate(tip = paste0("Hatch: ", Hatch, " (ext)<br>Day: ", round(day, 1), "<br>Temp: ", round(tempC, 1), "\u00b0C"))
        if (nrow(kd) > 0) {
          p <- p %>%
            add_trace(data = kd, x = ~day, y = ~tempC, type = "scatter",
                      mode = "lines+markers", line = list(shape = "hv", width = 2, color = colours[i]),
                      marker = list(size = 6, color = colours[i]),
                      name = hatches[i], text = ~tip, hoverinfo = "text",
                      legendgroup = hatches[i])
        }
        if (nrow(ed) > 0) {
          p <- p %>%
            add_trace(data = ed, x = ~day, y = ~tempC, type = "scatter",
                      mode = "lines", line = list(shape = "hv", width = 2, dash = "dash", color = colours[i]),
                      name = paste0(hatches[i], " ext"), text = ~tip, hoverinfo = "text",
                      legendgroup = hatches[i], showlegend = FALSE)
        }
      }
      # Ethylene markers
      if (!is.null(eth) && nrow(eth) > 0 && any(eth$ppm > 0)) {
        y_base <- min(all_df$tempC, na.rm = TRUE) - 0.3
        em <- eth %>% filter(ppm > 0) %>%
          mutate(dur_label = ifelse(duration == 0.5, "12 h", ifelse(duration == 1, "24 h", paste0(duration, " d"))),
                 tip = paste0("\u2b25 Ethylene<br>", ppm, " ppm<br>Day ", day, "<br>Duration: ", dur_label),
                 y = y_base)
        p <- p %>%
          add_trace(data = em, x = ~day, y = ~y, type = "scatter", mode = "markers",
                    marker = list(size = 12, color = "#8e44ad", symbol = "diamond"),
                    name = "Ethylene", text = ~tip, hoverinfo = "text",
                    showlegend = TRUE)
      }
      p %>% layout(xaxis = list(title = "Days from first cargo report"),
                   yaxis = list(title = "Set Point (\u00b0C)"),
                   margin = list(t = 10, b = 40, l = 50, r = 10),
                   legend = list(orientation = "v", x = 1.02, y = 1))
    })

    # ── Simulation: one run per hatch ─────────────────────────
    pred_sim <- eventReactive(input$run_predict, {
      req(pred_profile_data())
      var_code    <- variety_db(input$pred_variety)
      pd          <- pred_profile_data()
      hatch_names <- names(pd$hatch_profiles)
      use_calibrated <- !is.null(calibrated_params_rv())

      withProgress(message = "Running prediction...", value = 0, {
        if (use_calibrated) {
          cal_df <- calibrated_params_rv()
          params_all <- cal_df[, c("E0", "F0", "Ffix1")] %>% as.data.frame()
        } else {
          params_all <- tryCatch(load_param_bank(input$pred_softening, input$pred_variety),
                                 error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
          if (is.null(params_all)) return(NULL)
        }
        N_eff <- min(as.integer(input$pred_mc_n), nrow(params_all))
        all_trajs <- list(); all_summaries <- list(); obs_list <- list()

        for (i in seq_along(hatch_names)) {
          hatch <- hatch_names[i]; hp <- pd$hatch_profiles[[hatch]]
          setProgress(i / length(hatch_names) * 0.9,
                      message = paste0("Hatch ", hatch, " (", i, "/", length(hatch_names), ")..."))
          hatch_fmt <- pred_fmt_raw() %>%
            filter(voyage_name == input$pred_voyage, variety == var_code,
                   hatch_deck == hatch, fmt_testing_location == input$pred_location,
                   !is.na(firmness))
          init_f <- hatch_fmt %>% filter(date == hp$day0) %>% pull(firmness)
          F_n_h  <- if (length(init_f) && !all(is.na(init_f))) mean(init_f, na.rm=TRUE) else 7
          sd_n_h <- if (length(init_f) > 1) max(sd(init_f, na.rm=TRUE), 0.1) else 1
          idx    <- sample(seq_len(nrow(params_all)), N_eff, replace = FALSE)
          params_df <- params_all[idx, c("E0","F0","Ffix1")] %>% as.data.frame()
          res <- run_montecarlo_firmness(
            params_df = params_df,
            treatments = list(list(name = hatch, profile_df = hp$profile_df)),
            variety = input$pred_variety,
            use_init = !use_calibrated,
            F_n = F_n_h, std_n = sd_n_h, n_time_steps = 51)
          all_trajs[[hatch]]     <- res$trajectories
          all_summaries[[hatch]] <- res$summary
          obs_list[[hatch]] <- hatch_fmt %>%
            mutate(day = as.numeric(difftime(date, hp$day0, units = "days")),
                   treatment = hatch) %>%
            filter(day >= 0) %>%
            select(treatment, hatch_deck, date, day, fmt_testing_location, firmness)
        }
        param_source <- if (use_calibrated) "Calibrated (Bayesian)" else "Prior (Excel)"
        list(sim = list(trajectories = bind_rows(all_trajs),
                        summary      = bind_rows(all_summaries)),
             obs            = bind_rows(obs_list),
             hatch_profiles = pd$hatch_profiles,
             skipped        = pd$skipped,
             param_source   = param_source)
      })
    })

    # ── Landing-day summary cards ──────────────────────────────
    output$pred_landing_cards <- renderUI({
      req(pred_sim())
      v <- pred_sim(); sim_df <- v$sim$trajectories
      hatch_names <- unique(sim_df$treatment)

      cards <- purrr::map(hatch_names, function(h) {
        hp <- v$hatch_profiles[[h]]
        landing <- hp$landing_day
        hatch_rows <- sim_df %>% filter(treatment == h)
        closest_t  <- hatch_rows$time[which.min(abs(hatch_rows$time - landing))]
        at_landing <- hatch_rows %>% filter(time == closest_t)
        mean_f     <- round(mean(at_landing$firmness), 2)
        sd_f       <- round(sd(at_landing$firmness), 2)
        pct_soft   <- round(mean(at_landing$firmness <= 0.8) * 100, 1)
        bg <- if (pct_soft < 10) "#27ae60" else if (pct_soft < 30) "#e67e22" else "#c0392b"
        div(class = "card kpi-card border-0 shadow-sm",
            style = paste0("background:", bg, " !important; color:white; border-radius: var(--radius-lg);"),
            div(class = "card-body p-2",
                tags$p(class = "small mb-0 fw-semibold text-uppercase opacity-75", h),
                tags$p(class = "h5 fw-bold mb-0", paste0(mean_f, " \u00b1", sd_f, " kgf")),
                tags$p(class = "small mb-0 opacity-85",
                       paste0(pct_soft, "% \u2264 0.8 kgf  \u2502  day ", round(landing,1)))))
      })
      tagList(
        div(class = "d-flex flex-wrap gap-2 mb-3", cards),
      )
    })

    # ── Prediction plot ────────────────────────────────────────
    pred_plot_gg <- reactive({
      req(pred_sim())
      v      <- pred_sim()
      sim_df <- v$sim$trajectories
      obs_df <- v$obs
      avg_df <- sim_df %>% group_by(treatment, time) %>%
        summarise(mean_f = mean(firmness), .groups="drop")

      # Today and landing day markers per hatch
      today_df   <- purrr::map_dfr(names(v$hatch_profiles), function(h)
        tibble(treatment = h, today_day   = v$hatch_profiles[[h]]$today_day))
      landing_df <- purrr::map_dfr(names(v$hatch_profiles), function(h)
        tibble(treatment = h, landing_day = v$hatch_profiles[[h]]$landing_day))

      eth_sub <- if (isTRUE(input$pred_no_c2h4)) {
        "  \u2022  No ethylene"
      } else {
        edf <- pred_eth_rv()
        parts <- paste0(edf$ppm, " ppm d", edf$day, " \u00d7", edf$duration, "d")
        paste0("  \u2022  C\u2082H\u2084: ", paste(parts, collapse = "; "))
      }

      # Observed average line per hatch
      obs_avg_df <- if (nrow(obs_df) > 0) {
        obs_df %>% group_by(treatment, day) %>%
          summarise(mean_f = mean(firmness, na.rm = TRUE), .groups = "drop")
      } else { tibble(treatment = character(), day = numeric(), mean_f = numeric()) }

      p <- ggplot(sim_df, aes(x = time, y = firmness, group = replicate)) +
        geom_line(alpha = 0.08, linewidth = 0.3, colour = "#6baed6") +
        geom_hline(yintercept = 0.8, linetype = "dashed",
                   colour = "#e67e22", linewidth = 0.55) +
        geom_vline(data = today_df,
                   aes(xintercept = today_day), colour = "#3498db",
                   linetype = "solid", linewidth = 0.7, inherit.aes = FALSE) +
        geom_vline(data = landing_df,
                   aes(xintercept = landing_day), colour = "#8e44ad",
                   linetype = "dashed", linewidth = 0.8, inherit.aes = FALSE) +
        geom_line(data = avg_df, aes(x = time, y = mean_f, group = NULL),
                  linewidth = 1.2, colour = "#2c3e50") +
        facet_wrap(~ treatment, ncol = 4, scales = "free_x") +
        ylim(0,10) +
        labs(x = "Days from first FMT test (per hatch)",
             y = "Firmness (kgf)",
             title = paste("Voyage:", input$pred_voyage),
             subtitle = paste0(input$pred_variety, "  \u2022  ", input$pred_softening, eth_sub)) +
        theme_minimal(base_size = 12) +
        theme(panel.grid.minor = element_blank(),
              strip.text       = element_text(face = "bold"),
              plot.title       = element_text(face = "bold", size = 13),
              plot.subtitle    = element_text(colour = "grey40", size = 10))

      if (nrow(obs_df) > 0) {
        p <- p +
          geom_point(data = obs_df,
                     aes(x = day, y = firmness, shape = fmt_testing_location),
                     colour = "#e74c3c", size = 3, stroke = 0.9, inherit.aes = FALSE) +
          geom_line(data = obs_avg_df,
                    aes(x = day, y = mean_f, group = NULL),
                    linewidth = 1, colour = "#e74c3c", linetype = "solid",
                    inherit.aes = FALSE) +
          scale_shape_manual(values = c("Shore" = 16, "Vessel" = 17),
                             name = "FMT Location")
      }
      p
    })

    output$pred_plot <- renderPlot({
      pred_plot_gg()
    }, height = function() {
      req(pred_sim())
      n <- length(unique(pred_sim()$sim$trajectories$treatment))
      max(300, ceiling(n / 4) * 270)
    })

    output$download_pred_plot <- downloadHandler(
      filename = function()
        paste0("prediction_", input$pred_voyage, "_", Sys.Date(), ".png"),
      content = function(file) {
        req(pred_plot_gg())
        n     <- length(unique(pred_sim()$sim$trajectories$treatment))
        ncols <- min(4, n); nrows <- ceiling(n / ncols)
        ggplot2::ggsave(file, plot = pred_plot_gg(),
                        width = ncols * 5, height = nrows * 4,
                        dpi = 180, bg = "white")
      }
    )

    # ══════════════════════════════════════════════════════════════
    # ── Calibrate tab ─────────────────────────────────────────────

    output$cal_db_status_ui <- renderUI({
      s <- db_status(); if (s == "idle") return(
        div(class = "alert alert-info py-2 mb-2",
            tags$i(class = "bi bi-database me-1"),
            "Connect to the database in the ", tags$strong("Validate"), " tab first."))
      cls <- switch(s, connected = "alert alert-success py-2",
                       error     = "alert alert-danger py-2",
                       "alert alert-info py-2")
      div(class = paste(cls, "mb-2"), db_msg())
    })

    output$cal_year_ui <- renderUI({
      req(available_years())
      selectInput("cal_year", "Season", choices = available_years(), selected = available_years()[1])
    })

    cal_fmt_raw <- reactive({
      req(db_con_rv(), input$cal_year)
      DBI::dbGetQuery(db_con_rv(), sprintf("
        SELECT voyage_name, variety, CAST(date AS DATE) AS date,
               hatch_deck, fmt_testing_location,
               CAST(firmness AS FLOAT) AS firmness,
               CAST(fruit_tempurature_in_hold AS FLOAT) AS fruit_tempurature_in_hold
        FROM [ods].[Coolchain_FMT_TestSheet_Materialise]
        WHERE year = %d AND firmness IS NOT NULL AND voyage_name IS NOT NULL
              AND CAST(firmness AS FLOAT) > 0 AND CAST(firmness AS FLOAT) <= 12
      ", as.integer(input$cal_year))) %>%
        mutate(hatch_deck = str_replace_all(hatch_deck, "\\s+", ""), date = as.Date(date))
    })

    cal_cargo_raw <- reactive({
      req(db_con_rv(), input$cal_year)
      DBI::dbGetQuery(db_con_rv(), sprintf("
        SELECT Voyage_Name, Hatch_Code, Set_Point_Temperature_Celcius,
               CAST(CargoConditionReportDate AS DATE) AS CargoConditionReportDate
        FROM [ods].[Coolchain_CargoCondition_Report_Materialise]
        WHERE Voyage_Name IN (
          SELECT DISTINCT voyage_name
          FROM [ods].[Coolchain_FMT_TestSheet_Materialise]
          WHERE year = %d AND voyage_name IS NOT NULL
        )
      ", as.integer(input$cal_year))) %>%
        mutate(Hatch_Code = str_replace_all(Hatch_Code, "\\s+", ""),
               CargoConditionReportDate = as.Date(CargoConditionReportDate))
    })

    output$cal_voyage_ui <- renderUI({
      req(cal_fmt_raw())
      voyages <- cal_fmt_raw() %>%
        filter(variety == variety_db(input$cal_variety),
               fmt_testing_location == input$cal_location) %>%
        distinct(voyage_name) %>% arrange(voyage_name) %>% pull()
      if (!length(voyages)) return(helpText("No voyages found for selected filters."))
      selectInput("cal_voyage", "Voyage", choices = voyages)
    })

    output$cal_hatch_ui <- renderUI({
      req(cal_fmt_raw(), input$cal_voyage)
      hatches <- cal_fmt_raw() %>%
        filter(voyage_name == input$cal_voyage,
               variety == variety_db(input$cal_variety),
               fmt_testing_location == input$cal_location) %>%
        distinct(hatch_deck) %>% arrange(hatch_deck) %>% pull()
      if (!length(hatches)) return(helpText("No hatch/deck data for this voyage/location."))
      checkboxGroupInput("cal_hatches", "Hatch / Deck", choices = hatches, selected = hatches)
    })

    # Ethylene table for Calibrate tab
    cal_eth_rv <- reactiveVal(tibble(ppm = 100, day = 2, duration = 1))

    output$cal_eth_table <- renderDT(
      datatable(cal_eth_rv(), editable = TRUE, rownames = FALSE,
                options = list(dom = "t", paging = FALSE, ordering = FALSE),
                class = "compact cell-border stripe"),
      server = TRUE)

    observeEvent(input$cal_eth_table_cell_edit, {
      info <- input$cal_eth_table_cell_edit
      df <- cal_eth_rv()
      col_i <- as.integer(info$col) + 1L; row_i <- as.integer(info$row)
      val <- suppressWarnings(as.numeric(info$value))
      if (!is.na(val) && col_i >= 1 && col_i <= ncol(df) && row_i >= 1 && row_i <= nrow(df)) {
        df[row_i, col_i] <- val; cal_eth_rv(df)
      }
    })
    observeEvent(input$cal_eth_add, { cal_eth_rv(bind_rows(cal_eth_rv(), tibble(ppm = 0, day = 0, duration = 1))) })
    observeEvent(input$cal_eth_rm, { df <- cal_eth_rv(); if (nrow(df) > 1) cal_eth_rv(df[-nrow(df), ]) })

    # ── Calibration run (per-hatch + train/validate split) ──────
    cal_result <- eventReactive(input$run_calibrate, {
      req(cal_fmt_raw(), cal_cargo_raw(), input$cal_voyage, input$cal_hatches)
      var_code <- variety_db(input$cal_variety)
      eth_df <- if (isTRUE(input$cal_no_c2h4)) data.frame(ppm = 0, day = 0, duration = 0) else cal_eth_rv()
      cutoff_day <- input$cal_cutoff_day

      withProgress(message = "Calibrating per hatch...", value = 0, {
        params_all <- tryCatch(load_param_bank(input$cal_softening, input$cal_variety),
                               error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
        if (is.null(params_all)) return(NULL)

        hatch_results <- list()
        hatch_sims    <- list()
        hatch_obs_all <- list()
        skipped <- c()
        n_hatches <- length(input$cal_hatches)

        for (hi in seq_along(input$cal_hatches)) {
          hatch <- input$cal_hatches[hi]
          setProgress(hi / n_hatches * 0.8,
                      message = paste0("Hatch ", hatch, " (", hi, "/", n_hatches, ")..."))

          # Get FMT data for this hatch
          hatch_fmt <- cal_fmt_raw() %>%
            filter(voyage_name == input$cal_voyage, variety == var_code,
                   hatch_deck == hatch, fmt_testing_location == input$cal_location,
                   !is.na(firmness))
          if (!nrow(hatch_fmt)) { skipped <- c(skipped, paste(hatch, "(no FMT data)")); next }
          day0 <- min(hatch_fmt$date, na.rm = TRUE)

          # Temperature profile
          temp_all <- cal_cargo_raw() %>%
            filter(Voyage_Name == input$cal_voyage, Hatch_Code == hatch,
                   !is.na(Set_Point_Temperature_Celcius)) %>%
            mutate(day = as.numeric(difftime(CargoConditionReportDate, day0, units = "days"))) %>%
            group_by(day) %>%
            summarise(tempC = mean(Set_Point_Temperature_Celcius, na.rm = TRUE), .groups = "drop") %>%
            arrange(day)
          if (!nrow(temp_all)) { skipped <- c(skipped, paste(hatch, "(no cargo data)")); next }

          pre_fmt <- temp_all %>% filter(day <= 0) %>% tail(1)
          temp_df <- bind_rows(
            if (nrow(pre_fmt)) tibble(day = 0, tempC = pre_fmt$tempC) else tibble(),
            temp_all %>% filter(day > 0)
          ) %>% arrange(day)
          if (!nrow(temp_df)) { skipped <- c(skipped, paste(hatch, "(no temp data)")); next }

          max_day <- max(temp_df$day)
          profile_df <- bind_rows(
            temp_df %>% mutate(c2h4 = NA_real_),
            build_eth_rows(eth_df, max_day)
          ) %>% arrange(day)

          # Build env functions
          pf <- profile_df %>% mutate(across(c(day, tempC, c2h4), as.numeric))
          temp_rows <- pf %>% filter(!is.na(day), !is.na(tempC)) %>%
            group_by(day) %>% summarise(tempC = mean(tempC), .groups = "drop") %>% arrange(day)
          c2h4_rows <- pf %>% filter(!is.na(day), !is.na(c2h4)) %>%
            group_by(day) %>% summarise(c2h4 = mean(c2h4), .groups = "drop") %>% arrange(day)
          if (nrow(temp_rows) < 2 || nrow(c2h4_rows) < 2) {
            skipped <- c(skipped, paste(hatch, "(insufficient profile)")); next
          }
          if (temp_rows$day[1] != 0) temp_rows <- bind_rows(tibble(day = 0, tempC = temp_rows$tempC[1]), temp_rows)
          if (c2h4_rows$day[1] != 0) c2h4_rows <- bind_rows(tibble(day = 0, c2h4 = 0), c2h4_rows)
          temp_fun <- stats::approxfun(temp_rows$day, temp_rows$tempC, method = "linear", rule = 2, ties = mean)
          c2h4_fun <- stats::approxfun(c2h4_rows$day, c2h4_rows$c2h4, method = "linear", rule = 2, ties = mean)

          # Split obs into train (<=cutoff) and validate (>cutoff)
          obs_h <- hatch_fmt %>%
            mutate(day = as.numeric(difftime(date, day0, units = "days"))) %>%
            filter(day >= 0) %>%
            select(day, firmness)
          train_obs <- obs_h %>% filter(day <= cutoff_day)
          val_obs   <- obs_h %>% filter(day > cutoff_day)

          if (nrow(train_obs) < 2) { skipped <- c(skipped, paste(hatch, "(< 2 train points)")); next }

          # Initial firmness from this hatch's earliest data
          min_day_h <- min(train_obs$day)
          obs_init_f <- mean(train_obs$firmness[train_obs$day == min_day_h], na.rm = TRUE)
          if (is.na(obs_init_f)) obs_init_f <- mean(train_obs$firmness, na.rm = TRUE)

          # Calibrate on training data only
          cal_h <- run_calibration(
            obs_df = train_obs, temp_fun = temp_fun, c2h4_fun = c2h4_fun,
            variety = input$cal_variety, param_bank = params_all,
            calibrate_secondary = isTRUE(input$cal_secondary),
            n_posterior_samples = as.integer(input$cal_n_samples),
            obs_init_firmness = obs_init_f)

          if (is.null(cal_h)) { skipped <- c(skipped, paste(hatch, "(calibration failed)")); next }

          # Run MC simulation with calibrated params over full time range
          cal_params_df <- cal_h$posterior_samples[, c("E0", "F0", "Ffix1")]
          treatments_h <- list(list(name = hatch, profile_df = pf))
          cal_sim_h <- run_montecarlo_firmness(
            params_df = cal_params_df, treatments = treatments_h,
            variety = input$cal_variety, use_init = FALSE, n_time_steps = 101)

          # Prior simulation for comparison
          N_prior <- min(nrow(cal_params_df), nrow(params_all))
          idx <- sample(seq_len(nrow(params_all)), N_prior)
          prior_sim_h <- run_montecarlo_firmness(
            params_df = params_all[idx, ], treatments = treatments_h,
            variety = input$cal_variety, use_init = TRUE,
            F_n = obs_init_f, std_n = max(sd(train_obs$firmness), 0.5),
            n_time_steps = 101)

          # Compute validation metrics on held-out data
          if (nrow(val_obs)) {
            cal_avg <- cal_sim_h$trajectories %>%
              group_by(time) %>% summarise(pred_mean = mean(firmness), .groups = "drop")
            pred_fun <- stats::approxfun(cal_avg$time, cal_avg$pred_mean, rule = 2)
            val_obs$pred <- pred_fun(val_obs$day)
            val_obs$residual <- val_obs$firmness - val_obs$pred
            val_mape <- mean(abs(val_obs$residual / val_obs$firmness), na.rm = TRUE) * 100
            val_rmse <- sqrt(mean(val_obs$residual^2, na.rm = TRUE))
            val_bias <- mean(val_obs$residual, na.rm = TRUE)
          } else {
            val_mape <- NA; val_rmse <- NA; val_bias <- NA
          }

          # Train metrics
          cal_avg_t <- cal_sim_h$trajectories %>%
            group_by(time) %>% summarise(pred_mean = mean(firmness), .groups = "drop")
          pred_fun_t <- stats::approxfun(cal_avg_t$time, cal_avg_t$pred_mean, rule = 2)
          train_obs$pred <- pred_fun_t(train_obs$day)
          train_obs$residual <- train_obs$firmness - train_obs$pred
          train_mape <- mean(abs(train_obs$residual / train_obs$firmness), na.rm = TRUE) * 100
          train_rmse <- sqrt(mean(train_obs$residual^2, na.rm = TRUE))

          hatch_results[[hatch]] <- list(
            cal = cal_h, hatch = hatch,
            train_obs = train_obs, val_obs = val_obs, all_obs = obs_h,
            profile_df = pf, max_day = max_day,
            train_mape = train_mape, train_rmse = train_rmse,
            val_mape = val_mape, val_rmse = val_rmse, val_bias = val_bias,
            n_train = nrow(train_obs), n_val = nrow(val_obs))
          hatch_sims[[hatch]] <- list(cal_sim = cal_sim_h, prior_sim = prior_sim_h)
          hatch_obs_all[[hatch]] <- obs_h %>% mutate(hatch_deck = hatch)
        }

        if (!length(hatch_results)) {
          showNotification(paste0("All hatches skipped: ", paste(skipped, collapse = ", ")),
                           type = "error", duration = 10)
          return(NULL)
        }

        setProgress(1, message = "Done!")

        list(hatch_results = hatch_results, hatch_sims = hatch_sims,
             hatch_obs_all = hatch_obs_all, param_bank = params_all,
             cutoff_day = cutoff_day, skipped = skipped)
      })
    })

    # ── KPI Cards (per-hatch summary) ────────────────────────
    output$cal_kpi_cards <- renderUI({
      req(cal_result())
      cr <- cal_result(); hr <- cr$hatch_results
      cutoff <- cr$cutoff_day

      cards <- purrr::map(names(hr), function(h) {
        r <- hr[[h]]
        val_mape <- if (!is.na(r$val_mape)) paste0(round(r$val_mape, 1), "%") else "N/A"
        train_mape <- paste0(round(r$train_mape, 1), "%")
        bg <- if (is.na(r$val_mape)) "#64748b"
              else if (r$val_mape <= 15) "#27ae60"
              else if (r$val_mape <= 30) "#e67e22"
              else "#c0392b"
        div(class = "card kpi-card border-0 shadow-sm",
            style = paste0("background:", bg, " !important; color:white; border-radius:var(--radius-lg); min-width:150px;"),
            div(class = "card-body p-2",
                tags$p(class = "small mb-0 fw-semibold text-uppercase opacity-75", paste("Hatch", h)),
                tags$p(class = "h5 fw-bold mb-0", paste0("Val MAPE: ", val_mape)),
                tags$p(class = "small mb-0 opacity-85", paste0("Train MAPE: ", train_mape)),
                tags$p(class = "small mb-0 opacity-75",
                       paste0(r$n_train, " train | ", r$n_val, " val obs"))))
      })
      skip_msg <- if (length(cr$skipped)) paste0("Skipped: ", paste(cr$skipped, collapse = ", ")) else ""
      div(class = "mb-3",
        div(class = "d-flex flex-wrap gap-2 mb-2", cards),
        tags$small(class = "text-muted",
          paste0("Cutoff: day ", cutoff, "  |  ", length(hr), " hatches calibrated"),
          if (nzchar(skip_msg)) paste0("  |  ", skip_msg) else ""))
    })

    # ── Fit vs Observed plot (ggplot, faceted per hatch) ──────
    output$cal_fit_plot <- renderPlot({
      req(cal_result())
      cr <- cal_result(); hr <- cr$hatch_results; hs <- cr$hatch_sims
      cutoff <- cr$cutoff_day

      all_cal_avg <- list(); all_pri_avg <- list(); all_obs <- list()
      for (h in names(hr)) {
        cal_traj <- hs[[h]]$cal_sim$trajectories
        pri_traj <- hs[[h]]$prior_sim$trajectories
        all_cal_avg[[h]] <- cal_traj %>%
          group_by(time) %>%
          summarise(mean_f = mean(firmness), lo = quantile(firmness, 0.025),
                    hi = quantile(firmness, 0.975), .groups = "drop") %>%
          mutate(treatment = h)
        all_pri_avg[[h]] <- pri_traj %>%
          group_by(time) %>%
          summarise(mean_f = mean(firmness), lo = quantile(firmness, 0.025),
                    hi = quantile(firmness, 0.975), .groups = "drop") %>%
          mutate(treatment = h)
        all_obs[[h]] <- hr[[h]]$all_obs %>% mutate(treatment = h)
      }
      cal_df <- bind_rows(all_cal_avg); pri_df <- bind_rows(all_pri_avg)
      obs_df <- bind_rows(all_obs)
      ncols <- min(4, length(hr))

      ggplot() +
        annotate("rect", xmin = -Inf, xmax = cutoff, ymin = -Inf, ymax = Inf,
                 fill = "#f0f0f0", alpha = 0.5) +
        geom_vline(xintercept = cutoff, colour = "#e67e22", linewidth = 0.8, linetype = "dashed") +
        geom_ribbon(data = pri_df, aes(x = time, ymin = lo, ymax = hi),
                    fill = "#bdbdbd", alpha = 0.25) +
        geom_line(data = pri_df, aes(x = time, y = mean_f),
                  colour = "#bdbdbd", linewidth = 1, linetype = "dashed") +
        geom_ribbon(data = cal_df, aes(x = time, ymin = lo, ymax = hi),
                    fill = "#4f46e5", alpha = 0.15) +
        geom_line(data = cal_df, aes(x = time, y = mean_f),
                  colour = "#4f46e5", linewidth = 1.2) +
        geom_point(data = obs_df %>% filter(day <= cutoff),
                   aes(x = day, y = firmness), colour = "#2c3e50", alpha = 0.4, size = 1.5) +
        geom_point(data = obs_df %>% filter(day > cutoff),
                   aes(x = day, y = firmness), colour = "#e74c3c", alpha = 0.5, size = 1.5) +
        geom_hline(yintercept = 0.8, colour = "#e67e22", linetype = "dashed", linewidth = 0.5) +
        facet_wrap(~ treatment, ncol = ncols) +
        coord_cartesian(ylim = c(0, 12)) +
        labs(x = "Days from first FMT", y = "Firmness (kgf)",
             subtitle = paste0("Grey = training (\u2264 day ", cutoff,
                               ")  |  White = validation  |  Blue = calibrated 95% CI  |  Grey band = prior")) +
        theme_minimal(base_size = 12) +
        theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
              plot.subtitle = element_text(colour = "grey40", size = 10))
    }, height = function() {
      req(cal_result())
      nrows <- ceiling(length(cal_result()$hatch_results) / 4)
      max(300, nrows * 280)
    })

    # ── Validation metrics table ──────────────────────────────
    output$cal_val_metrics_table <- renderDT({
      req(cal_result())
      hr <- cal_result()$hatch_results
      tbl <- purrr::map_dfr(names(hr), function(h) {
        r <- hr[[h]]
        tibble(Hatch = h, `Train Obs` = r$n_train, `Val Obs` = r$n_val,
               `Train MAPE (%)` = round(r$train_mape, 1),
               `Val MAPE (%)` = if (!is.na(r$val_mape)) round(r$val_mape, 1) else NA,
               `Val RMSE (kgf)` = if (!is.na(r$val_rmse)) round(r$val_rmse, 2) else NA,
               `Val Bias (kgf)` = if (!is.na(r$val_bias)) round(r$val_bias, 2) else NA)
      })
      datatable(tbl, rownames = FALSE, options = list(dom = "t", paging = FALSE),
                class = "compact cell-border stripe") %>%
        formatStyle("Val MAPE (%)",
          backgroundColor = styleInterval(c(15, 30), c("#d4edda", "#fff3cd", "#f8d7da")))
    })

    # ── Posterior density plot (all hatches overlaid) ──────────
    output$cal_posterior_plot <- renderPlotly({
      req(cal_result())
      hr <- cal_result()$hatch_results
      pn <- c("E0", "F0", "Ffix1")
      hatch_cols <- c("#4f46e5", "#e74c3c", "#27ae60", "#e67e22",
                      "#8e44ad", "#3498db", "#2c3e50", "#d4a017")

      p <- plotly::subplot(
        lapply(seq_along(pn), function(i) {
          nm <- pn[i]; pp <- plot_ly()
          for (j in seq_along(names(hr))) {
            h <- names(hr)[j]; cal_h <- hr[[h]]$cal
            samp <- cal_h$posterior_samples[[nm]]
            d <- density(samp)
            col <- hatch_cols[((j - 1) %% length(hatch_cols)) + 1]
            pp <- pp %>% add_lines(x = d$x, y = d$y, line = list(color = col, width = 1.5),
              name = if (i == 1) h else NA, showlegend = (i == 1), legendgroup = h)
          }
          first_cal <- hr[[names(hr)[1]]]$cal
          pr_m <- first_cal$prior_means[nm]; pr_s <- first_cal$prior_sds[nm]
          x_rng <- seq(pr_m - 3*pr_s, pr_m + 3*pr_s, length.out = 200)
          pp <- pp %>% add_lines(x = x_rng, y = dnorm(x_rng, pr_m, pr_s),
            line = list(color = "#bdbdbd", width = 1.5, dash = "dash"),
            name = if (i == 1) "Prior" else NA, showlegend = (i == 1), legendgroup = "prior")
          pp %>% layout(xaxis = list(title = nm),
                        yaxis = list(title = if (i == 1) "Density" else ""))
        }),
        nrows = 1, shareY = FALSE, titleX = TRUE, titleY = TRUE
      ) %>% layout(legend = list(orientation = "h", y = -0.2), margin = list(t = 10))
      p
    })

    # ── Parameter summary table (per hatch) ───────────────────
    output$cal_param_table <- renderDT({
      req(cal_result())
      hr <- cal_result()$hatch_results
      pn <- c("E0", "F0", "Ffix1")
      tbl <- purrr::map_dfr(names(hr), function(h) {
        cal_h <- hr[[h]]$cal
        tibble(Hatch = h, Parameter = pn,
               `MAP Estimate` = round(cal_h$map_estimate[pn], 5),
               `Posterior Mean` = round(colMeans(cal_h$posterior_samples[pn]), 5),
               `Posterior SD` = round(sapply(cal_h$posterior_samples[pn], sd), 5),
               `95% CI Lower` = round(sapply(cal_h$posterior_samples[pn], quantile, 0.025), 5),
               `95% CI Upper` = round(sapply(cal_h$posterior_samples[pn], quantile, 0.975), 5))
      })
      datatable(tbl, rownames = FALSE, options = list(dom = "t", paging = FALSE, pageLength = 50),
                class = "compact cell-border stripe")
    })

    # ── Download calibrated params (all hatches) ──────────────
    output$cal_dl_params <- downloadHandler(
      filename = function() paste0("calibrated_params_", Sys.Date(), ".csv"),
      content = function(file) {
        req(cal_result())
        hr <- cal_result()$hatch_results
        all_params <- purrr::map_dfr(names(hr), function(h) {
          hr[[h]]$cal$posterior_samples %>% mutate(hatch = h)
        })
        write.csv(all_params, file, row.names = FALSE)
      }
    )

    # ── Bridge to Predict tab ─────────────────────────────────
    observeEvent(input$cal_use_params, {
      req(cal_result())
      # Combine posterior samples from all hatches
      hr <- cal_result()$hatch_results
      all_samples <- purrr::map_dfr(names(hr), function(h) hr[[h]]$cal$posterior_samples)
      calibrated_params_rv(all_samples)
      showNotification("Calibrated parameters loaded. Switch to Predict tab to use them.",
                       type = "message", duration = 5)
    })

  } # end if (db_available)

  # ══════════════════════════════════════════════════════════════
  # FMT Data Viewer tab
  # ══════════════════════════════════════════════════════════════
  json_data <- reactiveVal(NULL)    # named list of data frames
  json_selected <- reactiveVal(NULL)
  fmt_uploaded_files <- reactiveVal(character(0))  # track file names

  # ── Parse a single JSON file into named list of data frames ──
  parse_json_file <- function(filepath) {
    raw <- jsonlite::fromJSON(filepath, flatten = TRUE)
    tbls <- list()
    for (nm in names(raw)) {
      v <- raw[[nm]]
      if (is.data.frame(v) && nrow(v) > 0) {
        list_cols <- sapply(v, is.list)
        if (any(list_cols)) v <- v[, !list_cols, drop = FALSE]
        if (ncol(v) > 0) tbls[[nm]] <- v
      } else if (is.list(v) && length(v) > 0) {
        tryCatch({
          df <- as.data.frame(v)
          list_cols <- sapply(df, is.list)
          if (any(list_cols)) df <- df[, !list_cols, drop = FALSE]
          if (nrow(df) > 0 && ncol(df) > 0) tbls[[nm]] <- df
        }, error = function(e) NULL)
      }
    }
    tbls
  }

  # ── Parse a single Excel file into named list of data frames ──
  parse_excel_file <- function(filepath, filename) {
    tbls <- list()
    sheets <- readxl::excel_sheets(filepath)
    # Prefer "data export sheet" if it exists
    export_sheet <- sheets[tolower(sheets) %in% c("data export sheet", "sheet1")]
    other_sheets <- setdiff(sheets, c(export_sheet, "Settings", "Data Check"))

    for (s in sheets) {
      tryCatch({
        df <- suppressWarnings(readxl::read_excel(filepath, sheet = s, col_names = TRUE))
        # Skip sheets with mostly unnamed columns (template sheets like Lib. 1, etc.)
        named_pct <- mean(!grepl("^\\.\\.\\.\\d+$", names(df)))
        if (nrow(df) > 0 && ncol(df) > 1 && named_pct >= 0.3) {
          # Clean column names
          names(df) <- make.names(names(df), unique = TRUE)
          tag <- paste0(tools::file_path_sans_ext(filename), " :: ", s)
          tbls[[tag]] <- df
        }
      }, error = function(e) NULL)
    }
    tbls
  }

  # ── Handle file uploads (multiple, accumulating) ──
  observeEvent(input$fmt_files, {
    req(input$fmt_files)
    files <- input$fmt_files
    n_files <- nrow(files)
    existing_tbls <- json_data() %||% list()
    new_names <- c()

    withProgress(message = "Loading FMT data", value = 0, {
      for (i in seq_len(n_files)) {
        fname <- files$name[i]
        fpath <- files$datapath[i]
        ext <- tolower(tools::file_ext(fname))

        incProgress(0, detail = paste0("Reading ", fname, " (", i, "/", n_files, ")"))

        tryCatch({
          new_tbls <- if (ext == "json") {
            parse_json_file(fpath)
          } else if (ext %in% c("xlsx", "xls")) {
            parse_excel_file(fpath, fname)
          } else {
            showNotification(paste("Unsupported file:", fname), type = "warning")
            list()
          }

          incProgress(0.5 / n_files, detail = paste0("Merging ", fname))

          # Merge: if same table name exists, append rows; otherwise add new
          for (nm in names(new_tbls)) {
            if (nm %in% names(existing_tbls)) {
              common_cols <- intersect(names(existing_tbls[[nm]]), names(new_tbls[[nm]]))
              if (length(common_cols) > 0) {
                existing_tbls[[nm]] <- bind_rows(
                  existing_tbls[[nm]][, common_cols, drop = FALSE],
                  new_tbls[[nm]][, common_cols, drop = FALSE])
              }
            } else {
              existing_tbls[[nm]] <- new_tbls[[nm]]
            }
          }
          new_names <- c(new_names, fname)
          incProgress(0.5 / n_files, detail = paste0(fname, " done"))
        }, error = function(e) {
          showNotification(paste("Error reading", fname, ":", e$message), type = "error")
          incProgress(1 / n_files)
        })
      }

      incProgress(0, detail = "Building views...")
    })

    json_data(existing_tbls)
    fmt_uploaded_files(c(fmt_uploaded_files(), new_names))
    if (length(existing_tbls)) json_selected("__ALL__")

    total_rows <- sum(sapply(existing_tbls, nrow))
    showNotification(
      paste0(length(new_names), " file(s) loaded — ", length(existing_tbls),
             " tables, ", format(total_rows, big.mark = ","), " total rows"),
      type = "message", duration = 4)
  })

  # ── Clear all uploaded data ──
  observeEvent(input$fmt_clear, {
    json_data(NULL)
    json_selected(NULL)
    fmt_uploaded_files(character(0))
  })

  # ── Show uploaded file list ──
  output$fmt_file_list <- renderUI({
    fnames <- fmt_uploaded_files()
    if (!length(fnames)) return(NULL)
    tags$div(style = "margin:6px 0;",
      tags$p(style = "font-size:11px; color:#64748b; margin:0 0 4px;",
        paste0(length(fnames), " file(s) loaded:")),
      lapply(fnames, function(f) {
        ext <- tolower(tools::file_ext(f))
        bg <- if (ext == "json") "#dbeafe" else "#d1fae5"
        fg <- if (ext == "json") "#1e40af" else "#065f46"
        tags$span(style = paste0(
          "display:inline-block; padding:2px 8px; margin:1px 2px; border-radius:12px;",
          "font-size:11px; font-weight:500; background:", bg, "; color:", fg, ";"), f)
      })
    )
  })

  output$json_table_selector <- renderUI({
    req(json_data())
    nms <- names(json_data())
    if (!length(nms)) return(helpText("No tabular data found."))
    choices <- c("All files (combined)" = "__ALL__", setNames(nms, nms))
    tagList(
      selectInput("json_table_pick", "Table",
                  choices = choices, selected = "__ALL__"),
      tags$p(style = "font-size:11px; color:#94a3b8; margin:-8px 0 8px; line-height:1.3;",
        "Select a specific table to view individual file distributions and data.")
    )
  })

  observeEvent(input$json_table_pick, { json_selected(input$json_table_pick) })

  current_json_df <- reactive({
    req(json_data(), json_selected())
    sel <- json_selected()
    tbls <- json_data()
    if (sel == "__ALL__") {
      # Show the largest table by default
      biggest <- names(which.max(sapply(tbls, nrow)))
      df <- tbls[[biggest]]
    } else {
      req(sel %in% names(tbls))
      df <- tbls[[sel]]
    }
    req(is.data.frame(df))
    # Drop columns that are all NULL/NA or nested lists
    keep <- sapply(df, function(col) !all(is.na(col)) && !is.list(col))
    df[, keep, drop = FALSE]
  })

  output$json_summary <- renderUI({
    req(current_json_df())
    df <- current_json_df()
    div(class = "d-flex gap-3 mb-2",
      tags$span(class = "badge bg-primary", paste0(nrow(df), " rows")),
      tags$span(class = "badge bg-secondary", paste0(ncol(df), " columns")),
      tags$span(class = "badge bg-info text-dark", json_selected())
    )
  })

  output$json_meta_banner <- renderUI({
    meta <- json_header_meta()
    if (is.null(meta)) return(NULL)

    # Pill builder
    pill <- function(text, bg = "#f1f5f9", fg = "#334155") {
      tags$span(style = paste0(
        "display:inline-block; padding:3px 10px; margin:2px 3px; border-radius:20px;",
        "font-size:12px; font-weight:500; background:", bg, "; color:", fg, ";"), text)
    }
    # Label for a group
    grp_label <- function(text) {
      tags$span(style = "font-size:11px; font-weight:600; color:#94a3b8; text-transform:uppercase; letter-spacing:0.5px; margin-right:4px;", text)
    }
    # Row wrapper
    grp_row <- function(label, pills) {
      div(style = "display:flex; align-items:center; flex-wrap:wrap; gap:2px; margin-bottom:4px;",
          grp_label(label), pills)
    }

    rows <- list()
    # Year + Voyage on one row
    key_pills <- list()
    if (length(meta$calendarYear))
      key_pills <- c(key_pills, lapply(meta$calendarYear, function(y) pill(y, "#dbeafe", "#1e40af")))
    if (length(meta$voyageName))
      key_pills <- c(key_pills, lapply(meta$voyageName, function(v) pill(v, "#ede9fe", "#5b21b6")))
    if (length(key_pills)) rows <- c(rows, list(grp_row("Voyage", key_pills)))

    na_pill <- pill("\u2014", "#f1f5f9", "#cbd5e1")

    # Libraries
    if (length(meta$libraryCodes)) {
      lib_pills <- lapply(meta$libraryCodes, function(l) pill(l, "#d1fae5", "#065f46"))
    } else {
      lib_pills <- list(na_pill)
    }
    rows <- c(rows, list(grp_row("Libraries", lib_pills)))

    # Shipping Units — reflect the sidebar filter selection
    all_sus <- meta$shippingUnits
    su_sel <- input$fmt_su_filter
    if (length(all_sus)) {
      su_pills <- lapply(all_sus, function(s) {
        if (!is.null(su_sel) && !(s %in% su_sel))
          pill(s, "#f1f5f9", "#cbd5e1")
        else
          pill(s, "#fef3c7", "#92400e")
      })
    } else {
      su_pills <- list(na_pill)
    }
    rows <- c(rows, list(grp_row("Shipping Units", su_pills)))

    # Assessed dates
    if (length(meta$assessedDates)) {
      unique_days <- unique(substr(meta$assessedDates, 1, 10))
      unique_times <- unique(substr(meta$assessedDates, 12, 16))
      date_pills <- lapply(unique_days, function(d) pill(d, "#e0e7ff", "#3730a3"))
      time_pills <- lapply(unique_times, function(t) pill(t, "#f1f5f9", "#475569"))
      rows <- c(rows, list(grp_row("Date", date_pills)))
      if (length(unique_times) > 1) rows <- c(rows, list(grp_row("Times", time_pills)))
    } else {
      rows <- c(rows, list(grp_row("Date", list(na_pill))))
    }

    div(style = "background:white; border:1px solid #e2e8f0; border-radius:8px; padding:10px 14px; margin-bottom:12px;",
        rows)
  })

  output$json_dt <- renderDT({
    req(current_json_df())
    datatable(current_json_df(), rownames = FALSE, filter = "top",
              options = list(scrollX = TRUE, pageLength = 25,
                             dom = "lfrtip", autoWidth = FALSE),
              class = "compact stripe cell-border hover")
  })

  output$json_dl_csv <- downloadHandler(
    filename = function() paste0(json_selected(), "_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(current_json_df(), file)
  )

  output$json_dl_excel <- downloadHandler(
    filename = function() paste0(json_selected(), "_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      if (requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(current_json_df(), file)
      } else {
        readr::write_csv(current_json_df(), file)
      }
    }
  )

  # ── Firmness / temp distribution ──────────────────────────
  # Helper: natural sort (Lib. 1, Lib. 2, ... Lib. 9, Lib. 10)
  nat_sort <- function(x) {
    nums <- suppressWarnings(as.numeric(gsub("[^0-9]", "", x)))
    x[order(ifelse(is.na(nums), Inf, nums), x)]
  }

  # Helper: find a column by fuzzy name matching
  find_col <- function(df, patterns) {
    nms <- names(df)
    for (p in patterns) {
      hit <- grep(p, nms, ignore.case = TRUE, perl = TRUE, value = TRUE)
      if (length(hit)) return(hit[1])
    }
    NA_character_
  }

  # Helper: extract metadata (year, voyage, libraries, shipping units, dates)
  json_header_meta <- reactive({
    req(json_data())
    tbls <- json_data()
    sel <- json_selected()
    if (!is.null(sel) && sel != "__ALL__") tbls <- tbls[intersect(sel, names(tbls))]

    cal_years <- character(0); voy_names <- character(0)
    assessed_dates <- character(0); library_codes <- character(0)
    shipping_units <- character(0); headers <- NULL

    # Strategy 1: JSON-style headers table (calendarYear, voyageName)
    for (nm in names(tbls)) {
      df <- tbls[[nm]]
      if (is.data.frame(df) && all(c("calendarYear", "voyageName") %in% names(df))) {
        headers <- df
        cal_years <- unique(na.omit(as.character(df$calendarYear)))
        voy_names <- unique(na.omit(as.character(df$voyageName)))
        break
      }
    }

    # Strategy 2: Excel-style columns (Year, Voyage, Library, etc.)
    for (nm in names(tbls)) {
      df <- tbls[[nm]]
      if (!is.data.frame(df)) next
      yr_col <- find_col(df, c("^Year$"))
      voy_col <- find_col(df, c("^Voyage$", "^Voyage\\."))
      if (!is.na(yr_col) && !length(cal_years))
        cal_years <- unique(na.omit(as.character(df[[yr_col]])))
      if (!is.na(voy_col) && !length(voy_names))
        voy_names <- unique(na.omit(as.character(df[[voy_col]])))
      # Libraries from Excel
      lib_col <- find_col(df, c("^Library$", "^Library\\."))
      if (!is.na(lib_col) && !length(library_codes))
        library_codes <- nat_sort(unique(na.omit(as.character(df[[lib_col]]))))
      # Shipping units from Excel
      su_col <- find_col(df, c("^Cooling\\.Section", "^Hold\\.Deck$", "^shippingUnit$"))
      if (!is.na(su_col) && !length(shipping_units))
        shipping_units <- sort(unique(na.omit(as.character(df[[su_col]]))))
      # Date from Excel
      date_col <- find_col(df, c("^Date$", "^Time\\.Assessed"))
      if (!is.na(date_col) && !length(assessed_dates)) {
        raw <- as.character(df[[date_col]])
        assessed_dates <- sort(unique(na.omit(substr(raw, 1, 16))))
      }
    }

    # Strategy 3: JSON libraries table for dates/codes
    libs <- NULL
    for (nm in names(tbls)) {
      if (grepl("Librar", nm, ignore.case = TRUE)) { libs <- tbls[[nm]]; break }
    }
    if (is.null(libs)) {
      for (nm in names(tbls)) {
        df <- tbls[[nm]]
        if (is.data.frame(df) && "assessedDatetime" %in% names(df)) { libs <- df; break }
      }
    }
    if (!is.null(libs) && is.data.frame(libs)) {
      if ("assessedDatetime" %in% names(libs) && !length(assessed_dates)) {
        raw <- as.character(libs$assessedDatetime)
        assessed_dates <- sort(unique(na.omit(substr(raw, 1, 16))))
      }
      if ("libraryCode" %in% names(libs) && !length(library_codes))
        library_codes <- nat_sort(unique(na.omit(as.character(libs$libraryCode))))
      if ("shippingUnit" %in% names(libs) && !length(shipping_units))
        shipping_units <- sort(unique(na.omit(as.character(libs$shippingUnit))))
    }

    # Fallback: dates from JSON headers
    if (!length(assessed_dates) && !is.null(headers)) {
      date_col <- intersect(c("testDate", "testDateNzt", "testDateVlt"), names(headers))[1]
      if (!is.na(date_col)) {
        raw <- as.character(headers[[date_col]])
        assessed_dates <- sort(unique(na.omit(substr(raw, 1, 16))))
      }
    }

    if (!length(cal_years) && !length(voy_names)) return(NULL)

    list(
      calendarYear   = cal_years,
      voyageName     = voy_names,
      assessedDates  = assessed_dates,
      libraryCodes   = library_codes,
      shippingUnits  = shipping_units,
      headers        = headers
    )
  })

  # ── Available shipping units (for filter) ───────────────────
  fmt_available_sus <- reactive({
    req(json_data())
    tbls <- json_data()
    sus <- c()
    for (nm in names(tbls)) {
      df <- tbls[[nm]]
      if (!is.data.frame(df) || nrow(df) == 0) next
      su_col <- find_col(df, c("Cooling\\.Section", "Hold\\.Deck$", "Hold\\.Deck\\.", "^shippingUnit$"))
      if (!is.na(su_col)) sus <- c(sus, as.character(df[[su_col]]))
    }
    sus <- unique(sus[!is.na(sus) & sus != ""])
    nat_sort(sus)
  })

  output$fmt_su_filter_ui <- renderUI({
    sus <- fmt_available_sus()
    if (!length(sus)) return(NULL)
    tagList(
      hr(),
      span(class = "sec-lbl", "Filter by Shipping Unit"),
      checkboxGroupInput("fmt_su_filter", NULL, choices = sus, selected = sus)
    )
  })

  json_firmness <- reactive({
    req(json_data())
    tbls <- json_data()
    sel <- json_selected()
    if (!is.null(sel) && sel != "__ALL__") tbls <- tbls[intersect(sel, names(tbls))]
    all_results <- list()

    for (nm in names(tbls)) {
      df <- tbls[[nm]]
      if (!is.data.frame(df) || nrow(df) == 0) next

      # Format 1: JSON-style (firmnessA + firmnessB)
      if (all(c("firmnessA", "firmnessB") %in% names(df))) {
        res <- df %>%
          mutate(firmnessA = suppressWarnings(as.numeric(firmnessA)),
                 firmnessB = suppressWarnings(as.numeric(firmnessB)),
                 avg_firmness = (firmnessA + firmnessB) / 2) %>%
          filter(!is.na(avg_firmness), avg_firmness > 0, avg_firmness <= FIRMNESS_MAX)
        # Try to join variety + shipping unit from pallets
        for (pnm in names(tbls)) {
          pdf <- tbls[[pnm]]
          if (is.data.frame(pdf) && "offlineId" %in% names(pdf) && "palletOfflineId" %in% names(res)) {
            cols_to_pick <- c(palletOfflineId = "offlineId")
            if ("varietyCode" %in% names(pdf)) cols_to_pick <- c(cols_to_pick, varietyCode = "varietyCode")
            su_col <- find_col(pdf, c("^shippingUnit$"))
            if (!is.na(su_col)) cols_to_pick <- c(cols_to_pick, shipping_unit = su_col)
            pal_sel <- pdf %>% select(all_of(cols_to_pick)) %>%
              mutate(across(everything(), as.character)) %>% distinct()
            res <- res %>% mutate(palletOfflineId = as.character(palletOfflineId)) %>%
              left_join(pal_sel, by = "palletOfflineId")
            break
          }
        }
        if (nrow(res)) all_results <- c(all_results, list(res))
      }

      # Format 2: Excel-style (Firmness column, single value per fruit)
      firm_col <- find_col(df, c("^Firmness$", "^Firmness\\.(?!A|B)", "^Firmness\\.RAW$"))
      if (!is.na(firm_col) && !"firmnessA" %in% names(df)) {
        res <- df %>%
          mutate(avg_firmness = suppressWarnings(as.numeric(.data[[firm_col]]))) %>%
          filter(!is.na(avg_firmness), avg_firmness > 0, avg_firmness <= FIRMNESS_MAX)
        # Map variety
        var_col <- find_col(df, c("^Variety$", "^Variety\\.Group$", "^varietyCode$"))
        if (!is.na(var_col)) res$varietyCode <- as.character(res[[var_col]])
        # Map shipping unit
        su_col <- find_col(df, c("Cooling\\.Section", "Hold\\.Deck$", "Hold\\.Deck\\.", "shippingUnit"))
        if (!is.na(su_col)) res$shipping_unit <- as.character(res[[su_col]])
        res$source_file <- nm
        if (nrow(res)) all_results <- c(all_results, list(res))
      }
    }

    if (!length(all_results)) return(NULL)

    # Combine all results — keep common analysis columns + shipping_unit
    combined <- bind_rows(lapply(all_results, function(r) {
      out <- tibble(avg_firmness = r$avg_firmness)
      if ("varietyCode" %in% names(r)) out$varietyCode <- r$varietyCode
      if ("shipping_unit" %in% names(r)) out$shipping_unit <- r$shipping_unit
      out
    }))

    # Apply shipping unit filter from sidebar
    su_sel <- input$fmt_su_filter
    if (!is.null(su_sel) && "shipping_unit" %in% names(combined)) {
      combined <- combined %>% filter(is.na(shipping_unit) | shipping_unit %in% su_sel)
    }

    combined
  })

  output$json_firmness_panel <- renderUI({
    df <- json_firmness()
    if (is.null(df) || nrow(df) == 0) return(NULL)

    vals <- df$avg_firmness
    n     <- length(vals)
    mn    <- round(mean(vals), 2)
    med   <- round(median(vals), 2)
    sdev  <- round(sd(vals), 2)
    lo    <- round(min(vals), 2)
    hi    <- round(max(vals), 2)

    # TFR stats
    tfr_lo <- input$tfr_min %||% 2.5
    tfr_hi <- input$tfr_max %||% 6.0
    pct_below  <- round(100 * mean(vals < tfr_lo), 1)
    pct_within <- round(100 * mean(vals >= tfr_lo & vals <= tfr_hi), 1)
    pct_above  <- round(100 * mean(vals > tfr_hi), 1)
    pct_lt1    <- round(100 * mean(vals < 1), 1)

    # Stat cards — clean, no left-border accent
    stat_card <- function(label, value, colour = "#334155") {
      div(style = paste0("flex:1; min-width:65px; padding:5px 8px; background:#f8fafc;",
                          "border-radius:6px; border:1px solid #e2e8f0;"),
          tags$p(style = "margin:0; font-size:10px; color:#64748b; text-transform:uppercase; letter-spacing:0.3px; line-height:1.2;", label),
          tags$p(style = paste0("margin:1px 0 0; font-size:13px; font-weight:700; color:", colour, "; line-height:1.2;"), value))
    }

    div(class = "card mb-2",
      div(class = "card-header py-2 fw-semibold d-flex justify-content-between align-items-center",
        tags$span(style = "font-size:13px;", "Firmness Distribution"),
        tags$small(class = "text-muted fw-normal", paste0("n = ", n))),
      div(class = "card-body p-2",
        div(class = "d-flex gap-1 flex-wrap mb-1",
          stat_card("Mean", paste0(mn, " kgf"), "#2563eb"),
          stat_card("Median", paste0(med, " kgf"), "#7c3aed"),
          stat_card("SD", sdev, "#64748b"),
          stat_card("Range", paste0(lo, "\u2013", hi), "#0891b2")
        ),
        div(class = "d-flex gap-1 flex-wrap mb-1",
          stat_card("<1 kgf", paste0(pct_lt1, "%"),
                    if (pct_lt1 > 5) "#dc2626" else "#64748b"),
          stat_card(paste0("<", tfr_lo), paste0(pct_below, "%"),
                    if (pct_below > 20) "#dc2626" else "#64748b"),
          stat_card(paste0("In TFR (", tfr_lo, "\u2013", tfr_hi, ")"), paste0(pct_within, "%"),
                    if (pct_within >= 80) "#16a34a" else if (pct_within >= 50) "#ea580c" else "#dc2626"),
          stat_card(paste0(">", tfr_hi), paste0(pct_above, "%"),
                    if (pct_above > 20) "#dc2626" else "#64748b")
        ),
        plotlyOutput("json_firmness_plot", height = "240px")
      )
    )
  })

  output$json_firmness_plot <- renderPlotly({
    df <- json_firmness()
    req(df, nrow(df) > 0)

    vals <- df$avg_firmness
    mn <- mean(vals); med <- median(vals)
    has_variety <- "varietyCode" %in% names(df) && length(unique(df$varietyCode)) > 1
    variety_cols <- c("#2563eb", "#dc2626", "#16a34a", "#ea580c", "#7c3aed", "#0891b2")

    font_cfg <- list(family = "Inter, Segoe UI, system-ui, sans-serif")

    tfr_lo <- input$tfr_min %||% 2.5
    tfr_hi <- input$tfr_max %||% 6.0

    # Colour bars by TFR zone: below = warm red, within = slate blue, above = amber
    bar_col_below  <- "rgba(239,68,68,0.45)"    # red-500
    bar_col_within <- "rgba(59,130,246,0.45)"    # blue-500
    bar_col_above  <- "rgba(245,158,11,0.45)"    # amber-500
    bdr_below  <- "#ef4444"
    bdr_within <- "#3b82f6"
    bdr_above  <- "#f59e0b"

    # Pre-compute bins as % of total, colour by TFR zone
    bin_size <- 0.5
    breaks <- seq(0, 12, by = bin_size)
    bin_mid <- head(breaks, -1) + bin_size / 2
    counts <- as.numeric(table(cut(vals, breaks, right = FALSE, include.lowest = TRUE)))
    pcts <- round(100 * counts / length(vals), 2)

    # Assign zone colour per bin
    zone_col <- ifelse(bin_mid < tfr_lo, bar_col_below,
                  ifelse(bin_mid <= tfr_hi, bar_col_within, bar_col_above))
    zone_bdr <- ifelse(bin_mid < tfr_lo, bdr_below,
                  ifelse(bin_mid <= tfr_hi, bdr_within, bdr_above))
    zone_lbl <- ifelse(bin_mid < tfr_lo, paste0("< ", tfr_lo),
                  ifelse(bin_mid <= tfr_hi, "In TFR", paste0("> ", tfr_hi)))

    p <- plot_ly(x = bin_mid, y = pcts, type = "bar", width = bin_size * 0.9,
                 marker = list(color = zone_col, line = list(color = zone_bdr, width = 0.8)),
                 text = paste0(pcts, "%"), hoverinfo = "text+x",
                 showlegend = FALSE)

    p %>% layout(
      font = font_cfg,
      xaxis = list(title = list(text = "Average Firmness (kgf)", font = list(size = 12, color = "#334155")),
                   gridcolor = "#f1f5f9", zeroline = FALSE, tickfont = list(size = 10, color = "#64748b")),
      yaxis = list(title = list(text = "Percentage (%)", font = list(size = 12, color = "#334155")),
                   gridcolor = "#f1f5f9", zeroline = FALSE, tickfont = list(size = 10, color = "#64748b")),
      plot_bgcolor = "#ffffff", paper_bgcolor = "#ffffff",
      shapes = list(
        # TFR boundary lines — dark green, distinct from bar colours
        list(type = "line", x0 = tfr_lo, x1 = tfr_lo, y0 = 0, y1 = 1, yref = "paper",
             line = list(color = "#15803d", width = 1.8, dash = "dot")),
        list(type = "line", x0 = tfr_hi, x1 = tfr_hi, y0 = 0, y1 = 1, yref = "paper",
             line = list(color = "#15803d", width = 1.8, dash = "dot")),
        # Mean & Median — thin, unobtrusive
        list(type = "line", x0 = mn, x1 = mn, y0 = 0, y1 = 1, yref = "paper",
             line = list(color = "#1e293b", width = 1.5, dash = "dash")),
        list(type = "line", x0 = med, x1 = med, y0 = 0, y1 = 1, yref = "paper",
             line = list(color = "#6b21a8", width = 1.5, dash = "dashdot"))
      ),
      annotations = list(
        list(x = mn, y = 1.02, yref = "paper", text = paste0("<b>Mean</b> ", round(mn, 2)),
             showarrow = FALSE, yanchor = "bottom", font = list(color = "#1e293b", size = 10)),
        list(x = med, y = 0.94, yref = "paper", text = paste0("<b>Med</b> ", round(med, 2)),
             showarrow = FALSE, yanchor = "bottom", font = list(color = "#6b21a8", size = 10)),
        list(x = tfr_lo, y = 0.88, yref = "paper", text = paste0("<b>", tfr_lo, "</b>"),
             showarrow = FALSE, xanchor = "right", font = list(color = "#15803d", size = 9)),
        list(x = tfr_hi, y = 0.88, yref = "paper", text = paste0("<b>", tfr_hi, "</b>"),
             showarrow = FALSE, xanchor = "left", font = list(color = "#15803d", size = 9))
      ),
      legend = list(orientation = "h", y = -0.18, font = list(size = 10)),
      margin = list(t = 30, b = 40, l = 50, r = 20)
    ) %>% plotly::config(displayModeBar = FALSE)
  })

  # ── Temperature Distribution ─────────────────────────
  json_temp_data <- reactive({
    req(json_data())
    tbls <- json_data()
    sel <- json_selected()
    if (!is.null(sel) && sel != "__ALL__") tbls <- tbls[intersect(sel, names(tbls))]
    result <- tibble()

    for (nm in names(tbls)) {
      df <- tbls[[nm]]
      if (!is.data.frame(df) || nrow(df) == 0) next

      # Determine shipping unit for this table's rows
      su_vec <- rep(NA_character_, nrow(df))
      su_col <- find_col(df, c("Cooling\\.Section", "Hold\\.Deck$", "Hold\\.Deck\\.", "shippingUnit"))
      if (!is.na(su_col)) su_vec <- as.character(df[[su_col]])
      # JSON: join pallets for SU
      if (is.na(su_col) && "palletOfflineId" %in% names(df)) {
        for (pnm in names(tbls)) {
          pdf <- tbls[[pnm]]
          if (!is.data.frame(pdf)) next
          su_col2 <- find_col(pdf, c("^shippingUnit$"))
          if (!is.na(su_col2) && "offlineId" %in% names(pdf)) {
            idx <- match(df$palletOfflineId, pdf$offlineId)
            su_vec <- as.character(pdf[[su_col2]][idx])
            break
          }
        }
      }

      # JSON-style columns
      if ("fruitInHoldTemperature" %in% names(df)) {
        vals <- suppressWarnings(as.numeric(df$fruitInHoldTemperature))
        keep <- !is.na(vals) & vals > -20 & vals < 60
        if (any(keep)) result <- bind_rows(result, tibble(temp = vals[keep], type = "In Hold", shipping_unit = su_vec[keep]))
      }
      if ("fruitAtTestTemperature" %in% names(df)) {
        vals <- suppressWarnings(as.numeric(df$fruitAtTestTemperature))
        keep <- !is.na(vals) & vals > -20 & vals < 60
        if (any(keep)) result <- bind_rows(result, tibble(temp = vals[keep], type = "At Test", shipping_unit = su_vec[keep]))
      }

      # Excel-style columns (Fruit Tempurature in Hold, Fruit Temperature at Test)
      hold_col <- find_col(df, c("Fruit.Tempurature.in.Hold", "Fruit.Temperature.in.Hold",
                                  "fruitInHoldTemperature"))
      test_col <- find_col(df, c("Fruit.Temperature.at.Test", "fruitAtTestTemperature"))
      if (!is.na(hold_col) && hold_col != "fruitInHoldTemperature") {
        vals <- suppressWarnings(as.numeric(df[[hold_col]]))
        keep <- !is.na(vals) & vals > -20 & vals < 60
        if (any(keep)) result <- bind_rows(result, tibble(temp = vals[keep], type = "In Hold", shipping_unit = su_vec[keep]))
      }
      if (!is.na(test_col) && test_col != "fruitAtTestTemperature") {
        vals <- suppressWarnings(as.numeric(df[[test_col]]))
        keep <- !is.na(vals) & vals > -20 & vals < 60
        if (any(keep)) result <- bind_rows(result, tibble(temp = vals[keep], type = "At Test", shipping_unit = su_vec[keep]))
      }
    }

    req(nrow(result) > 0)

    # Apply shipping unit filter from sidebar
    su_sel <- input$fmt_su_filter
    if (!is.null(su_sel) && "shipping_unit" %in% names(result)) {
      result <- result %>% filter(is.na(shipping_unit) | shipping_unit %in% su_sel)
    }

    result
  })

  output$json_temp_panel <- renderUI({
    df <- json_temp_data()
    req(df, nrow(df) > 0)

    stats_by_type <- df %>%
      group_by(type) %>%
      summarise(n = n(), mn = round(mean(temp), 1), med = round(median(temp), 1),
                sdev = round(sd(temp), 1), lo = round(min(temp), 1), hi = round(max(temp), 1),
                .groups = "drop")

    # Teal for Hold, warm amber for Test — professional, distinct
    type_colours <- c("In Hold" = "#0f766e", "At Test" = "#b45309")

    stat_card <- function(label, value, colour = "#334155") {
      div(style = paste0("flex:1; min-width:65px; padding:5px 8px; background:#f8fafc;",
                          "border-radius:6px; border:1px solid #e2e8f0;"),
          tags$p(style = "margin:0; font-size:10px; color:#64748b; text-transform:uppercase; letter-spacing:0.3px; line-height:1.2;", label),
          tags$p(style = paste0("margin:1px 0 0; font-size:13px; font-weight:700; color:", colour, "; line-height:1.2;"), value))
    }

    total_n <- nrow(df)
    cards_list <- list()
    for (i in seq_len(nrow(stats_by_type))) {
      row <- stats_by_type[i, ]
      tc <- if (row$type %in% names(type_colours)) type_colours[row$type] else "#64748b"
      cards_list <- c(cards_list, list(
        stat_card(paste0(row$type, " Mean"), paste0(row$mn, "\u00b0C"), tc),
        stat_card(paste0(row$type, " Med"), paste0(row$med, "\u00b0C"), tc),
        stat_card(paste0(row$type, " SD"), as.character(row$sdev), "#475569"),
        stat_card(paste0(row$type, " Range"), paste0(row$lo, "\u2013", row$hi), "#475569")
      ))
    }

    div(class = "card mb-2",
      div(class = "card-header py-2 fw-semibold d-flex justify-content-between align-items-center",
        tags$span(style = "font-size:13px;", "Temperature Distribution"),
        tags$small(class = "text-muted fw-normal", paste0("n = ", total_n))),
      div(class = "card-body p-2",
        div(class = "d-flex gap-1 flex-wrap mb-1", cards_list),
        plotlyOutput("json_temp_plot", height = "260px")
      )
    )
  })

  output$json_temp_plot <- renderPlotly({
    df <- json_temp_data()
    req(df, nrow(df) > 0)

    types <- unique(df$type)
    colours <- c("In Hold" = "#0f766e", "At Test" = "#b45309")
    fill_colours <- c("In Hold" = "rgba(15,118,110,0.30)", "At Test" = "rgba(180,83,9,0.30)")
    font_cfg <- list(family = "Inter, Segoe UI, system-ui, sans-serif")

    p <- plot_ly(alpha = 0.65) %>% layout(barmode = "overlay")
    for (tp in types) {
      sub <- df %>% filter(type == tp)
      col <- if (tp %in% names(colours)) colours[tp] else "#64748b"
      fcol <- if (tp %in% names(fill_colours)) fill_colours[tp] else "rgba(100,116,139,0.30)"
      p <- p %>% add_histogram(
        x = sub$temp, name = tp, histnorm = "percent",
        xbins = list(start = 0, end = 12, size = 0.5),
        marker = list(color = fcol, line = list(color = col, width = 1.2)))
    }

    # Mean & median lines per type
    shapes <- list(); annotations <- list()
    y_pos <- 1.02
    for (tp in types) {
      sub <- df %>% filter(type == tp)
      mn <- mean(sub$temp); med <- median(sub$temp)
      col <- if (tp %in% names(colours)) colours[tp] else "#64748b"
      shapes <- c(shapes, list(
        list(type = "line", x0 = mn, x1 = mn, y0 = 0, y1 = 1, yref = "paper",
             line = list(color = col, width = 2, dash = "dash")),
        list(type = "line", x0 = med, x1 = med, y0 = 0, y1 = 1, yref = "paper",
             line = list(color = col, width = 1.5, dash = "dot"))))
      annotations <- c(annotations, list(
        list(x = mn, y = y_pos, yref = "paper",
             text = paste0("<b>", tp, " Mean</b> ", round(mn, 1), "\u00b0C"),
             showarrow = FALSE, yanchor = "bottom",
             font = list(color = col, size = 11)),
        list(x = med, y = y_pos - 0.08, yref = "paper",
             text = paste0("<b>", tp, " Median</b> ", round(med, 1), "\u00b0C"),
             showarrow = FALSE, yanchor = "bottom",
             font = list(color = col, size = 10))))
      y_pos <- y_pos - 0.18
    }

    p %>% layout(
      font = font_cfg,
      xaxis = list(title = list(text = "Temperature (\u00b0C)", font = list(size = 13, color = "#334155")),
                   gridcolor = "#f1f5f9", zeroline = FALSE, tickfont = list(size = 11, color = "#64748b")),
      yaxis = list(title = list(text = "Percentage (%)", font = list(size = 13, color = "#334155")),
                   gridcolor = "#f1f5f9", zeroline = FALSE, tickfont = list(size = 11, color = "#64748b")),
      plot_bgcolor = "#ffffff", paper_bgcolor = "#ffffff",
      shapes = shapes, annotations = annotations,
      legend = list(orientation = "h", y = -0.18, font = list(size = 10)),
      margin = list(t = 30, b = 40, l = 50, r = 20)) %>%
      plotly::config(displayModeBar = FALSE)
  })

  # ── Firmness by Shipping Unit (faceted) ─────────────────────
  # Reactive: firmness with shipping unit attached
  fmt_firmness_by_su <- reactive({
    req(json_data())
    tbls <- json_data()
    sel <- json_selected()
    if (!is.null(sel) && sel != "__ALL__") tbls <- tbls[intersect(sel, names(tbls))]
    all_results <- list()

    for (nm in names(tbls)) {
      df <- tbls[[nm]]
      if (!is.data.frame(df) || nrow(df) == 0) next

      # Excel-style: has Firmness + Cooling.Section / Hold.Deck
      firm_col <- find_col(df, c("^Firmness$", "^Firmness\\.(?!A|B)", "^Firmness\\.RAW$"))
      su_col <- find_col(df, c("Cooling\\.Section", "Hold\\.Deck$", "Hold\\.Deck\\.", "shippingUnit"))
      if (!is.na(firm_col) && !is.na(su_col)) {
        res <- df %>%
          mutate(avg_firmness = suppressWarnings(as.numeric(.data[[firm_col]])),
                 shipping_unit = as.character(.data[[su_col]])) %>%
          filter(!is.na(avg_firmness), avg_firmness > 0, avg_firmness <= FIRMNESS_MAX,
                 !is.na(shipping_unit), shipping_unit != "")
        if (nrow(res)) all_results <- c(all_results, list(res[, c("avg_firmness", "shipping_unit")]))
      }

      # JSON-style: firmnessA/B — need to join through pallets to get shipping unit
      if (all(c("firmnessA", "firmnessB") %in% names(df))) {
        res <- df %>%
          mutate(firmnessA = suppressWarnings(as.numeric(firmnessA)),
                 firmnessB = suppressWarnings(as.numeric(firmnessB)),
                 avg_firmness = (firmnessA + firmnessB) / 2) %>%
          filter(!is.na(avg_firmness), avg_firmness > 0, avg_firmness <= FIRMNESS_MAX)
        # Try to get shipping unit from pallets table (which has shippingUnit or equivalent)
        for (pnm in names(tbls)) {
          pdf <- tbls[[pnm]]
          if (!is.data.frame(pdf)) next
          su_col2 <- find_col(pdf, c("^shippingUnit$"))
          if (!is.na(su_col2) && "offlineId" %in% names(pdf) && "palletOfflineId" %in% names(res)) {
            pal_su <- pdf %>%
              select(palletOfflineId = offlineId, shipping_unit = !!sym(su_col2)) %>%
              mutate(shipping_unit = as.character(shipping_unit)) %>%
              distinct()
            res <- res %>% left_join(pal_su, by = "palletOfflineId") %>%
              filter(!is.na(shipping_unit), shipping_unit != "")
            break
          }
        }
        if ("shipping_unit" %in% names(res) && nrow(res))
          all_results <- c(all_results, list(res[, c("avg_firmness", "shipping_unit")]))
      }
    }

    if (!length(all_results)) return(NULL)
    bind_rows(all_results)
  })

  # Reactive: temperature with shipping unit attached
  fmt_temp_by_su <- reactive({
    req(json_data())
    tbls <- json_data()
    sel <- json_selected()
    if (!is.null(sel) && sel != "__ALL__") tbls <- tbls[intersect(sel, names(tbls))]
    result <- tibble()

    for (nm in names(tbls)) {
      df <- tbls[[nm]]
      if (!is.data.frame(df) || nrow(df) == 0) next

      # Determine shipping unit per row
      su_vec <- rep(NA_character_, nrow(df))
      su_col <- find_col(df, c("Cooling\\.Section", "Hold\\.Deck$", "Hold\\.Deck\\.", "shippingUnit"))
      if (!is.na(su_col)) su_vec <- as.character(df[[su_col]])
      if (is.na(su_col) && "palletOfflineId" %in% names(df)) {
        for (pnm in names(tbls)) {
          pdf <- tbls[[pnm]]
          if (!is.data.frame(pdf)) next
          su_col2 <- find_col(pdf, c("^shippingUnit$"))
          if (!is.na(su_col2) && "offlineId" %in% names(pdf)) {
            idx <- match(df$palletOfflineId, pdf$offlineId)
            su_vec <- as.character(pdf[[su_col2]][idx])
            break
          }
        }
      }

      # JSON-style
      if ("fruitInHoldTemperature" %in% names(df)) {
        vals <- suppressWarnings(as.numeric(df$fruitInHoldTemperature))
        keep <- !is.na(vals) & vals > -20 & vals < 60 & !is.na(su_vec) & su_vec != ""
        if (any(keep)) result <- bind_rows(result, tibble(temp = vals[keep], type = "In Hold", shipping_unit = su_vec[keep]))
      }
      if ("fruitAtTestTemperature" %in% names(df)) {
        vals <- suppressWarnings(as.numeric(df$fruitAtTestTemperature))
        keep <- !is.na(vals) & vals > -20 & vals < 60 & !is.na(su_vec) & su_vec != ""
        if (any(keep)) result <- bind_rows(result, tibble(temp = vals[keep], type = "At Test", shipping_unit = su_vec[keep]))
      }

      # Excel-style
      hold_col <- find_col(df, c("Fruit.Tempurature.in.Hold", "Fruit.Temperature.in.Hold", "fruitInHoldTemperature"))
      test_col <- find_col(df, c("Fruit.Temperature.at.Test", "fruitAtTestTemperature"))
      if (!is.na(hold_col) && hold_col != "fruitInHoldTemperature") {
        vals <- suppressWarnings(as.numeric(df[[hold_col]]))
        keep <- !is.na(vals) & vals > -20 & vals < 60 & !is.na(su_vec) & su_vec != ""
        if (any(keep)) result <- bind_rows(result, tibble(temp = vals[keep], type = "In Hold", shipping_unit = su_vec[keep]))
      }
      if (!is.na(test_col) && test_col != "fruitAtTestTemperature") {
        vals <- suppressWarnings(as.numeric(df[[test_col]]))
        keep <- !is.na(vals) & vals > -20 & vals < 60 & !is.na(su_vec) & su_vec != ""
        if (any(keep)) result <- bind_rows(result, tibble(temp = vals[keep], type = "At Test", shipping_unit = su_vec[keep]))
      }
    }

    if (!nrow(result)) return(NULL)
    result
  })

  output$fmt_su_panel <- renderUI({
    df <- fmt_firmness_by_su()
    if (is.null(df) || nrow(df) == 0) {
      return(div(style = "padding:40px; text-align:center; color:#94a3b8;",
        tags$p(style = "font-size:16px;", "No shipping unit data available."),
        tags$p(style = "font-size:13px;", "Upload Excel FMT files with Cooling Section / Hold Deck columns.")))
    }

    sus <- nat_sort(unique(df$shipping_unit))
    tfr_lo <- input$tfr_min %||% 2.5
    tfr_hi <- input$tfr_max %||% 6.0

    stat_card <- function(label, value, colour = "#334155") {
      div(style = paste0("flex:1; min-width:55px; padding:4px 6px; background:#f8fafc;",
                          "border-radius:6px; border:1px solid #e2e8f0;"),
          tags$p(style = "margin:0; font-size:9px; color:#64748b; text-transform:uppercase; letter-spacing:0.3px; line-height:1.2;", label),
          tags$p(style = paste0("margin:1px 0 0; font-size:12px; font-weight:700; color:", colour, "; line-height:1.2;"), value))
    }

    # Get temp data by SU
    temp_df <- fmt_temp_by_su()

    # Build a card per shipping unit
    cards <- lapply(sus, function(su) {
      sub <- df %>% filter(shipping_unit == su)
      vals <- sub$avg_firmness
      n <- length(vals)
      mn <- round(mean(vals), 2)
      med <- round(median(vals), 2)
      pct_below  <- round(100 * mean(vals < tfr_lo), 1)
      pct_within <- round(100 * mean(vals >= tfr_lo & vals <= tfr_hi), 1)
      pct_above  <- round(100 * mean(vals > tfr_hi), 1)
      pct_lt1    <- round(100 * mean(vals < 1), 1)

      # Temperature stats for this SU
      temp_cards <- list()
      if (!is.null(temp_df) && nrow(temp_df) > 0) {
        su_temp <- temp_df %>% filter(shipping_unit == su)
        if (nrow(su_temp) > 0) {
          for (tp in c("In Hold", "At Test")) {
            tp_vals <- su_temp %>% filter(type == tp) %>% pull(temp)
            if (length(tp_vals) > 0) {
              tc <- if (tp == "In Hold") "#0f766e" else "#b45309"
              temp_cards <- c(temp_cards, list(
                stat_card(paste0(tp, " Avg"), paste0(round(mean(tp_vals), 1), "\u00b0C"), tc),
                stat_card(paste0(tp, " Range"), paste0(round(min(tp_vals), 1), "\u2013", round(max(tp_vals), 1)), tc)
              ))
            }
          }
        }
      }

      div(style = "flex:1; min-width:280px; max-width:50%;",
        div(class = "card mb-2",
          div(class = "card-header py-1 fw-semibold d-flex justify-content-between align-items-center",
            tags$span(style = "font-size:13px;", su),
            tags$small(class = "text-muted fw-normal", paste0("n = ", n))),
          div(class = "card-body p-2",
            div(class = "d-flex gap-1 flex-wrap mb-1",
              stat_card("Mean", paste0(mn, " kgf"), "#1e293b"),
              stat_card("Median", paste0(med, " kgf"), "#6b21a8"),
              stat_card("<1", paste0(pct_lt1, "%"), if (pct_lt1 > 5) "#dc2626" else "#64748b"),
              stat_card(paste0("<", tfr_lo), paste0(pct_below, "%"),
                        if (pct_below > 20) "#dc2626" else "#64748b"),
              stat_card("In TFR", paste0(pct_within, "%"),
                        if (pct_within >= 80) "#16a34a" else if (pct_within >= 50) "#ea580c" else "#dc2626"),
              stat_card(paste0(">", tfr_hi), paste0(pct_above, "%"),
                        if (pct_above > 20) "#dc2626" else "#64748b")
            ),
            if (length(temp_cards)) div(class = "d-flex gap-1 flex-wrap mb-1", temp_cards),
            plotlyOutput(paste0("su_plot_", gsub("[^A-Za-z0-9]", "_", su)), height = "200px")
          )
        )
      )
    })

    # Render individual plots
    lapply(sus, function(su) {
      plot_id <- paste0("su_plot_", gsub("[^A-Za-z0-9]", "_", su))
      local({
        my_su <- su
        output[[plot_id]] <- renderPlotly({
          df_all <- fmt_firmness_by_su()
          req(df_all)
          sub <- df_all %>% filter(shipping_unit == my_su)
          vals <- sub$avg_firmness
          mn <- mean(vals); med <- median(vals)
          tfr_lo <- input$tfr_min %||% 2.5
          tfr_hi <- input$tfr_max %||% 6.0
          font_cfg <- list(family = "Inter, Segoe UI, system-ui, sans-serif")

          # Pre-compute bins as % of total, colour by TFR zone
          bin_size <- 0.5
          breaks <- seq(0, 12, by = bin_size)
          bin_mid <- head(breaks, -1) + bin_size / 2
          counts <- as.numeric(table(cut(vals, breaks, right = FALSE, include.lowest = TRUE)))
          pcts <- round(100 * counts / length(vals), 2)

          zone_col <- ifelse(bin_mid < tfr_lo, "rgba(239,68,68,0.45)",
                        ifelse(bin_mid <= tfr_hi, "rgba(59,130,246,0.45)", "rgba(245,158,11,0.45)"))
          zone_bdr <- ifelse(bin_mid < tfr_lo, "#ef4444",
                        ifelse(bin_mid <= tfr_hi, "#3b82f6", "#f59e0b"))

          p <- plot_ly(x = bin_mid, y = pcts, type = "bar", width = bin_size * 0.9,
                       marker = list(color = zone_col, line = list(color = zone_bdr, width = 0.8)),
                       text = paste0(pcts, "%"), hoverinfo = "text+x", showlegend = FALSE)

          p %>% layout(
            font = font_cfg,
            xaxis = list(title = "", gridcolor = "#f1f5f9", zeroline = FALSE,
                         tickfont = list(size = 9, color = "#64748b")),
            yaxis = list(title = list(text = "%", font = list(size = 10, color = "#64748b")),
                         gridcolor = "#f1f5f9", zeroline = FALSE,
                         tickfont = list(size = 9, color = "#64748b")),
            plot_bgcolor = "#ffffff", paper_bgcolor = "#ffffff",
            shapes = list(
              list(type = "line", x0 = tfr_lo, x1 = tfr_lo, y0 = 0, y1 = 1, yref = "paper",
                   line = list(color = "#15803d", width = 1.2, dash = "dot")),
              list(type = "line", x0 = tfr_hi, x1 = tfr_hi, y0 = 0, y1 = 1, yref = "paper",
                   line = list(color = "#15803d", width = 1.2, dash = "dot")),
              list(type = "line", x0 = mn, x1 = mn, y0 = 0, y1 = 1, yref = "paper",
                   line = list(color = "#1e293b", width = 1.2, dash = "dash")),
              list(type = "line", x0 = med, x1 = med, y0 = 0, y1 = 1, yref = "paper",
                   line = list(color = "#6b21a8", width = 1.2, dash = "dashdot"))
            ),
            margin = list(t = 5, b = 25, l = 35, r = 15)) %>%
          plotly::config(displayModeBar = FALSE)
        })
      })
    })

    div(style = "padding:8px;",
      tags$h6(style = "color:#334155; margin-bottom:8px;", "Firmness Distribution by Shipping Unit"),
      div(class = "d-flex flex-wrap gap-2", cards)
    )
  })
}

shinyApp(ui, server)
