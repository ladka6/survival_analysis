
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
  as.data.frame(nd)
}

# koef table
safe_coef_table <- function(fit){
  cf_raw <- tryCatch(summary(fit)$coefficients, error=function(e) NULL)
  if (is.null(cf_raw)) return(data.frame(term=character(), est=numeric(), se=numeric(), lcl=numeric(), ucl=numeric()))
  cf <- as.data.frame(cf_raw, stringsAsFactors=FALSE)
  if (nrow(cf)==0 || ncol(cf)==0) return(data.frame(term=character(), est=numeric(), se=numeric(), lcl=numeric(), ucl=numeric()))
  cn <- tolower(gsub("[^a-z0-9]+","", colnames(cf))); colnames(cf) <- cn
  selcol <- function(cands){ x <- intersect(cands, colnames(cf)); if (length(x)) x[1] else NA_character_ }
  est <- selcol(c("est","estimate","coef","beta")); se <- selcol(c("se","stderr","secoef"))
  lcl <- selcol(c("lcl","lower","lower95","l95"));  ucl <- selcol(c("ucl","upper","upper95","u95"))
  if (is.na(est)) return(data.frame(term=character(), est=numeric(), se=numeric(), lcl=numeric(), ucl=numeric()))
  ton <- function(z) suppressWarnings(as.numeric(z))
  out <- data.frame(term=rownames(cf),
                    est=ton(cf[[est]]),
                    se =if(!is.na(se))  ton(cf[[se]])  else NA_real_,
                    lcl=if(!is.na(lcl)) ton(cf[[lcl]]) else NA_real_,
                    ucl=if(!is.na(ucl)) ton(cf[[ucl]]) else NA_real_,
                    stringsAsFactors=FALSE)
  rownames(out) <- NULL; out
}

tidy_time_ratios <- function(fit){
  cf <- safe_coef_table(fit); if (nrow(cf)==0) return(cf)
  keep <- c("(Intercept)","intercept","pets","carpet_score","total_usage_time")
  cf2 <- cf[ cf$term %in% keep | grepl("pets|carpet|usage|total_usage", cf$term), , drop=FALSE ]
  if (nrow(cf2)==0) return(cf2)
  cf2$time_ratio <- exp(cf2$est); cf2$lower_tr <- exp(cf2$lcl); cf2$upper_tr <- exp(cf2$ucl)
  dplyr::select(cf2, term, est, se, time_ratio, lower_tr, upper_tr)
}

get_L10 <- function(fit, newdata=NULL){
  s <- summary(fit, type="quantile", quantiles=0.10, newdata=newdata)
  ton <- function(z) suppressWarnings(as.numeric(z))
  if (is.list(s)) c(L10=ton(s[[1]]$est), L10_lcl=ton(s[[1]]$lcl), L10_ucl=ton(s[[1]]$ucl))
  else            c(L10=ton(s[,"est"]),  L10_lcl=ton(s[,"lcl"]),  L10_ucl=ton(s[,"ucl"]))
}

plot_model_only <- function(fit, d, comp, label){
  if (!PLOT_MODEL_ONLY) return(invisible(NULL))
  newd <- build_newdata(d)
  tgrid <- seq(0, quantile(d$time, 0.995, na.rm=TRUE), length.out=200)
  cur   <- summary(fit, type="survival", newdata=newd, t=tgrid)[[1]]
  df    <- data.frame(time=cur$time, surv=cur$est)
  p <- ggplot(df, aes(time, surv)) + geom_line(linewidth=0.7) +
       labs(title=paste0(toupper(comp), ": ", label, " (model-only)"),
            x="Time", y="S(t)") + theme_minimal()
  file <- file.path(OUT_DIR, paste0("para_Model_only_", label, "_", comp, ".png"))
  ggsave(file, p, width=6.5, height=4.2, dpi=160)
  file
}

read_existing_km <- function(component){
  if (!USE_EXISTING_KM) return(NULL)
  for (f in c(paste0("KM_summary_",component,"_calendar.csv"),
              paste0("KM_summary_",component,".csv")))
    if (file.exists(f)) return(tryCatch(readr::read_csv(f, show_col_types=FALSE), error=function(e) NULL))
  NULL
}

# ---------- data cleansing ----------
sanitize_df <- function(df){
  keep <- intersect(c("time","event","pets","carpet_score","total_usage_time"), names(df))
  df <- df[, keep, drop=FALSE]
  for (nm in intersect(c("carpet_score","total_usage_time"), names(df)))
    df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  df$event <- as.integer(df$event)
  df <- df[!is.na(df$time) & df$time >= 0, , drop=FALSE]
  df$time[df$time==0] <- .Machine$double.eps
  # drop zero-variance covariates
  zvars <- sapply(df[, setdiff(colnames(df), c("time","event")), drop=FALSE],
                  function(x){ if(is.numeric(x)) sd(x, na.rm=TRUE) else if(is.logical(x)) length(unique(x[!is.na(x)])) else if(is.factor(x)) length(levels(x)) else NA })
  dropcols <- names(zvars)[ (!is.na(zvars)) & (zvars==0 | zvars==1) ]
  if (length(dropcols)) df <- df[, setdiff(colnames(df), dropcols), drop=FALSE]
  # drop rows with non-finite numerics
  numc <- intersect(c("carpet_score","total_usage_time"), names(df))
  if (length(numc)){
    ok <- rep(TRUE, nrow(df))
    for (nm in numc) ok <- ok & (is.na(df[[nm]]) | is.finite(df[[nm]]))
    df <- df[ok, , drop=FALSE]
  }
  df
}


fit_exp_safe <- function(form, df){
  tryCatch(flexsurvreg(formula=form, data=df, dist="exponential"), error=function(e) NULL)
}

fit_wb_safe <- function(form, df){
  form0 <- as.formula("Surv(time,event) ~ 1")
  fit0  <- tryCatch(flexsurvreg(formula=form0, data=df, dist="weibull"), error=function(e) NULL)
  inits <- NULL
  if (!is.null(fit0)) {
    par0 <- coef(fit0)
    inits <- c(par0[ c("shape","scale") ], rep(0, length(attr(terms(form), "term.labels"))))
  }
  fit1 <- tryCatch(flexsurvreg(formula=form, data=df, dist="weibull",
                               inits=inits, hessian=FALSE), error=function(e) NULL)
  if (!is.null(fit1)) return(fit1)
  tryCatch(flexsurvreg(formula=form0, data=df, dist="weibull"), error=function(e) NULL)
}

OUT_DIR <- "."
DATA_CANDIDATES <- c("DirtSlurper3100.csv","DirtSlurper3100 (1).csv","DirtSlurper3100_clean.csv")
USE_EXISTING_KM <- TRUE
PLOT_MODEL_ONLY <- TRUE

  ev <- sum(data$event) 
  if (ev == 0L) { warning("All censored (no failures) for ", component, " — skipping."); return(invisible(NULL)) }

  surv_obj <- create_surv(data)

  fit_exp <- fit_exp_safe(form, df)
  fit_wb  <- if (n_ev >= 3) fit_wb_safe(form, df) else NULL

  AIC_wb  <- if (!is.null(fit_wb)) AIC(fit_wb) else NA_real_
  AIC_exp <- if (!is.null(fit_exp)) AIC(fit_exp) else NA_real_

  aic_tbl <- data.frame(component=comp, clock=clock,
                        AIC_Weibull=AIC_wb, AIC_Exponential=AIC_exp,
                        dAIC_Exp_minus_WB = if (is.na(AIC_wb) || is.na(AIC_exp)) NA_real_ else AIC_exp - AIC_wb)
  write.csv(aic_tbl, file.path(OUT_DIR, paste0("para_AIC_",comp,"_",clock,".csv")), row.names=FALSE)

  tr_out <- data.frame()
  if (!is.null(fit_wb))  tr_out <- dplyr::bind_rows(tr_out, tidy_time_ratios(fit_wb)  |> dplyr::mutate(model="Weibull",    component=comp, clock=clock))
  if (!is.null(fit_exp)) tr_out <- dplyr::bind_rows(tr_out, tidy_time_ratios(fit_exp) |> dplyr::mutate(model="Exponential", component=comp, clock=clock))
  write.csv(tr_out, file.path(OUT_DIR, paste0("para_time_ratios_",comp,"_",clock,".csv")), row.names=FALSE)

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

for (comp in c("battery","ir","impact")) run_component(comp, clock="calendar")
cat("\nParametric finished — utils.R pipeline reused (v4). Outputs prefixed with para_.\n")
