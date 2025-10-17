
# ========== 0) Packages ==========
pkgs <- c("survival","survminer","ggplot2","dplyr","scales","ragg")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) {
  install.packages(p, repos = "https://cloud.r-project.org")
}
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(scales)
agg_png <- ragg::agg_png 

# ========== 1) Utils & data ==========
source("utils.R")  

df <- read_dirt("DirtSlurper3100.csv")
cat("Dataset:", nrow(df), "rows\n")

# ========== 2) Helpers ==========
normalize_event <- function(x) {
  if (is.logical(x)) {
    y <- x
  } else if (is.numeric(x)) {
    y <- x != 0
  } else {
    xc <- toupper(trimws(as.character(x)))
    y <- xc %in% c("1","TRUE","T","Y","YES","DAMAGE","FAIL","FAILED","EVENT")
  }
  y[is.na(y)] <- FALSE
  y
}

aft_to_weibull <- function(fit) {
  stopifnot(inherits(fit, "survreg"))
  list(shape = as.numeric(1/fit$scale),
       scale = as.numeric(exp(unname(coef(fit)[1]))))
}

weibull_quantile_ci <- function(fit, p, alpha = 0.05) {
  pr <- predict(fit, type = "quantile", p = 1 - p, se.fit = TRUE)
  t  <- as.numeric(pr$fit)
  se <- as.numeric(pr$se.fit)
  z  <- qnorm(1 - alpha/2)
  c(est = t, lcl = pmax(0, t - z*se), ucl = t + z*se)
}

km_diag_points <- function(sf) {
  S <- sf$surv; tt <- sf$time
  keep <- which(S > 0 & S < 1 & tt > 0)
  if (!length(keep)) return(NULL)
  data.frame(x = log(tt[keep]), y = log(-log(S[keep])))
}

# ========== 3) Core runner ==========
run_weibull <- function(component, clock = "calendar") {
  cat("\n====================================================\n")
  cat("Component:", toupper(component), " | Clock:", clock, "\n")
  cat("====================================================\n")

  data <- prepare_surv_data(df, component, clock) %>%
          filter(!is.na(time) & time > 0)

  if (!"event" %in% names(data)) stop("ERROR: 'event' column missing after prepare_surv_data().")
  data$event <- normalize_event(data$event)

  if (nrow(data) == 0) { warning("No usable rows for ", component); return(invisible(NULL)) }

  ev <- sum(data$event) 
  if (ev == 0L) { warning("All censored (no failures) for ", component, " — skipping."); return(invisible(NULL)) }

  surv_obj <- create_surv(data)

  # ---- Weibull AFT ----
  fit_wb <- survreg(surv_obj ~ 1, data = data, dist = "weibull")
  pars   <- aft_to_weibull(fit_wb)

  # Key quantiles + CI
  q50 <- weibull_quantile_ci(fit_wb, p = 0.50)
  q90 <- weibull_quantile_ci(fit_wb, p = 0.90)

  # Exponential comparison
  fit_exp <- try(survreg(surv_obj ~ 1, data = data, dist = "exponential"), silent = TRUE)
  AIC_wb  <- AIC(fit_wb)
  AIC_exp <- if (inherits(fit_exp, "survreg")) AIC(fit_exp) else NA_real_
  dAIC    <- if (is.na(AIC_exp)) NA_real_ else (AIC_exp - AIC_wb)

  # ---- KM (nonparametric) ----
  sf_km <- survfit(Surv(time, event) ~ 1, data = data)

  # ---- Parametric Weibull curve ----
  t_max  <- max(data$time, na.rm = TRUE)
  t_grid <- seq(0, t_max, length.out = 400)
  S_wb   <- exp(- (t_grid / pars$scale)^pars$shape)
  wb_df  <- data.frame(time = t_grid, surv = S_wb)

  # ---- Plot: KM + Weibull overlay ----
  plt <- ggsurvplot(
    sf_km,
    data = data,
    conf.int = TRUE,
    risk.table = FALSE,   
    ggtheme = theme_minimal(),
    title = paste0("Weibull vs KM — ", toupper(component)),
    xlab = ifelse(clock=="calendar","Time (days)","Usage (hours)"),
    ylab = "Survival probability S(t)"
  )
  plt$plot <- plt$plot +
    geom_line(data = wb_df, aes(x = time, y = surv), linewidth = 1, color = "firebrick3") +
    annotate("text", x = t_max * 0.75, y = 0.1, label = "Weibull fit", color = "firebrick3", hjust = 0)

  # ---- Diagnostic: Weibull probability plot ----
  diag_df <- km_diag_points(sf_km)
  if (!is.null(diag_df)) {
    r2 <- summary(lm(y ~ x, data = diag_df))$r.squared
    diag_plot <- ggplot(diag_df, aes(x, y)) +
      geom_point(alpha = 0.6) + geom_smooth(method = "lm", se = FALSE) +
      labs(title = paste0("Weibull probability plot — ", toupper(component),
                          sprintf("  (R^2 = %.3f)", r2)),
           x = "log(time)", y = "log(-log(Ŝ(t)))") +
      theme_minimal()
  } else {
    diag_plot <- ggplot() + theme_void() + labs(title = "Diagnostic unavailable (insufficient KM points)")
  }

  # ---- Save PNGs ----
  ggsave(paste0("WB_", component, "_", clock, ".png"),
         plt$plot, width = 7, height = 5, dpi = 300, device = agg_png)
  ggsave(paste0("WB_diag_", component, "_", clock, ".png"),
         diag_plot, width = 6.5, height = 5, dpi = 300, device = agg_png)

  res <- list(
    component           = as.character(component),
    clock               = as.character(clock),
    N                   = as.integer(nrow(data)),
    events              = as.integer(ev),
    event_pct           = as.numeric(round(100 * ev / nrow(data), 1)),
    shape_k             = as.numeric(pars$shape),
    scale_lambda        = as.numeric(pars$scale),
    median_est          = as.numeric(q50["est"]),
    median_lcl          = as.numeric(q50["lcl"]),
    median_ucl          = as.numeric(q50["ucl"]),
    L10_est             = as.numeric(q90["est"]),
    L10_lcl             = as.numeric(q90["lcl"]),
    L10_ucl             = as.numeric(q90["ucl"]),
    AIC_Weibull         = as.numeric(AIC_wb),
    AIC_Exponential     = as.numeric(AIC_exp),
    dAIC_Exp_minus_WB   = as.numeric(dAIC),
    logLik_Weibull      = as.numeric(logLik(fit_wb))
  )
  out <- as.data.frame(res, stringsAsFactors = FALSE, optional = TRUE)
  row.names(out) <- NULL  
  write.csv(out, paste0("WB_summary_", component, "_", clock, ".csv"), row.names = FALSE)

  
}

# ========== 4) Run all three components ==========
for (comp in c("battery","ir","impact")) run_weibull(comp, clock = "calendar")