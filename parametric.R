##########################################
# parametric.R — Parametric Survival (AFT) ONLY — v4 (defansif) + para_ prefix
# - All outputs (CSV/PNG) are prefixed with "para_" so you can tell they come from Parametric
##########################################

suppressPackageStartupMessages({
  if (!require("dplyr"))     install.packages("dplyr")
  if (!require("readr"))     install.packages("readr")
  if (!require("flexsurv"))  install.packages("flexsurv")
  if (!require("survival"))  install.packages("survival")
  if (!require("ggplot2"))   install.packages("ggplot2")
})

library(dplyr); library(readr); library(flexsurv); library(survival); library(ggplot2)

# === Minimal KM + Olasılık Grafikleri (tek dosya) ===
.safe_clip <- function(p, eps=1e-9) pmin(pmax(p, eps), 1 - eps)

# df (time,event) -> KM zaman serisi (time, surv)
.km_series <- function(df){
  stopifnot(all(c("time","event") %in% names(df)))
  fit <- survfit(Surv(time, event) ~ 1, data = df)
  d <- data.frame(time = fit$time, surv = fit$surv)
  d <- d[is.finite(d$time) & d$time > 0 & is.finite(d$surv), , drop=FALSE]
  d
}

.plot_weibull_just <- function(km_ts, comp, clock, out_dir="."){
  if (nrow(km_ts) < 3) return(invisible(NULL))
  d <- within(km_ts, { surv <- .safe_clip(surv); X <- log(time); Y <- log(-log(surv)) })
  p <- ggplot(d, aes(X, Y)) +
    geom_point(alpha=0.8) +
    geom_smooth(method="lm", se=FALSE, linewidth=0.7) +
    labs(title=paste0(toupper(comp),": Weibull justification"),
         x="log(t)", y="log(-log S(t))") +
    theme_minimal()
  ggsave(file.path(out_dir, paste0("para_Weibull_just_", comp, "_", clock, ".png")),
         p, width=7, height=4.6, dpi=160)
}

.plot_exp_just <- function(km_ts, comp, clock, out_dir="."){
  if (nrow(km_ts) < 3) return(invisible(NULL))
  d <- within(km_ts, { surv <- .safe_clip(surv) })
  p <- ggplot(d, aes(time, log(surv))) +
    geom_point(alpha=0.8) +
    geom_smooth(method="lm", se=FALSE, linewidth=0.7) +
    labs(title=paste0(toupper(comp),": Exponential justification"),
         x="Time", y="log S(t)") +
    theme_minimal()
  ggsave(file.path(out_dir, paste0("para_Exp_just_", comp, "_", clock, ".png")),
         p, width=7, height=4.6, dpi=160)
}

# Opsiyonel PP-plot (parametrik fit karşılaştırması)
.plot_pp <- function(fit, df, comp, label, clock, out_dir="."){
  d <- df[df$event==1 & is.finite(df$time) & df$time>0, , drop=FALSE]
  if (nrow(d) < 3) return(invisible(NULL))
  d <- d[order(d$time), ]
  n <- nrow(d)
  d$F_emp <- (seq_len(n) - 0.5) / n
  s <- summary(fit, type="cdf", t=d$time)[[1]]
  d$F_mod <- .safe_clip(s$est)
  p <- ggplot(d, aes(F_emp, F_mod)) +
    geom_point(alpha=0.85) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    labs(title=paste0(toupper(comp),": PP-Plot (",label,")"),
         x="Empirical CDF", y="Model CDF") +
    theme_minimal()
  ggsave(file.path(out_dir, paste0("para_PP_", label, "_", comp, "_", clock, ".png")),
         p, width=6.4, height=6.0, dpi=160)
}

# --- KM vs Weibull overlay (KM step in BLUE, Weibull model in RED) ---
plot_km_vs_weibull <- function(fit_wb, df, comp, clock="calendar", out_dir="."){
  if (is.null(fit_wb)) return(invisible(NULL))

  # KM step
  km <- survfit(Surv(time, event) ~ 1, data = df)
  km_df <- data.frame(time = c(0, km$time), surv = c(1, km$surv))

  # Weibull model on same grid (or fallback grid)
  newd  <- build_newdata(df)
  tgrid <- sort(unique(km_df$time))
  if (length(tgrid) < 2) tgrid <- seq(0, quantile(df$time, 0.995, na.rm=TRUE), length.out=200)
  cur   <- summary(fit_wb, type="survival", newdata=newd, t=tgrid)[[1]]
  wb_df <- data.frame(time = cur$time, surv = cur$est)

  p <- ggplot() +
    geom_step(data=km_df, aes(time, surv, color="KM"), linewidth=1.1) +           # BLUE
    geom_line(data=wb_df, aes(time, surv, color="Weibull model"), linewidth=1.2) +# RED
    scale_color_manual(values = c("KM"="#377EB8", "Weibull model"="#E41A1C"), name = "Curve") +
    labs(title = paste0(toupper(comp), ": KM vs Weibull (", clock, ")"),
         x = "Time", y = "S(t)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right")
  ggsave(file.path(out_dir, paste0("para_KM_vs_Weibull_", comp, "_", clock, ".png")),
         p, width=7.2, height=4.6, dpi=220)
}

# Tek çağrıyla tüm grafikleri üret
emit_km_probability_plots <- function(df, comp, clock="calendar", out_dir=".", fit_wb=NULL, fit_exp=NULL){
  km_ts <- .km_series(df)
  .plot_weibull_just(km_ts, comp, clock, out_dir)
  .plot_exp_just(km_ts, comp, clock, out_dir)
  if (!is.null(fit_wb))  .plot_pp(fit_wb,  df, comp, "Weibull",    clock, out_dir)
  if (!is.null(fit_exp)) .plot_pp(fit_exp, df, comp, "Exponential", clock, out_dir)
}
# === /Minimal KM + Olasılık Grafikleri ===

OUT_DIR <- "."
DATA_CANDIDATES <- c("DirtSlurper3100.csv","DirtSlurper3100 (1).csv","DirtSlurper3100_clean.csv")
USE_EXISTING_KM <- TRUE
PLOT_MODEL_ONLY <- TRUE

# ---------- utils.R (KM ile aynı) ----------
source("utils.R")   # read_dirt(), prepare_surv_data()

load_data_like_km <- function(){
  for(p in DATA_CANDIDATES) if (file.exists(p)) {
    df <- read_dirt(p); cat(sprintf("Loaded %d rows × %d columns from %s.\n", nrow(df), ncol(df), p)); return(df)
  }
  csvs <- list.files(pattern="\\.csv$", ignore.case=TRUE)
  if (length(csvs)>=1){ df <- read_dirt(csvs[1]); cat(sprintf("Loaded %d rows × %d columns from %s.\n", nrow(df), ncol(df), csvs[1])); return(df) }
  stop("CSV not found.")
}

# ---------- yardımcılar ----------
make_formula <- function(df){
  cand <- c("pets","carpet_score","total_usage_time")
  rhs  <- cand[cand %in% names(df)]
  if (length(rhs)==0) as.formula("Surv(time,event) ~ 1")
  else as.formula(paste0("Surv(time,event) ~ ", paste(rhs, collapse=" + ")))
}

# tip-uyumlu newdata
build_newdata <- function(df){
  nd <- list()
  for (nm in intersect(c("pets","carpet_score","total_usage_time"), names(df))){
    x <- df[[nm]]
    if (is.logical(x)){
      mode_val <- as.logical(names(sort(table(x), decreasing=TRUE))[1]); nd[[nm]] <- mode_val
    } else if (is.factor(x)){
      mode_val <- names(sort(table(x), decreasing=TRUE))[1]; nd[[nm]] <- factor(mode_val, levels=levels(x))
    } else if (is.numeric(x)){
      nd[[nm]] <- mean(x, na.rm=TRUE)
    } else {
      mode_val <- names(sort(table(x), decreasing=TRUE))[1]; nd[[nm]] <- mode_val
    }
  }
  as.data.frame(nd)
}

# koef tablosu güvenli okuma
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

# ---------- veri temizliği ----------
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

# ---------- güvenli fit yardımcıları ----------
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

# ---------- ana koşucu ----------
run_component <- function(comp, clock="calendar"){
  raw <- load_data_like_km()
  df  <- prepare_surv_data(raw, comp, clock=clock)
  df  <- sanitize_df(df)

  n_ev <- sum(df$event, na.rm=TRUE)
  form <- make_formula(df)

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

  nd <- build_newdata(df)
  L10_wb  <- if (!is.null(fit_wb))  get_L10(fit_wb,  newdata=nd) else c(L10=NA, L10_lcl=NA, L10_ucl=NA)
  L10_exp <- if (!is.null(fit_exp)) get_L10(fit_exp, newdata=nd) else c(L10=NA, L10_lcl=NA, L10_ucl=NA)

  km_ref <- read_existing_km(comp); km_L10 <- NA_real_; km_note <- "KM summary not found"
  if (!is.null(km_ref)) {
    nmz <- tolower(names(km_ref))
    if ("l10" %in% nmz)      km_L10 <- suppressWarnings(as.numeric(km_ref[[which(nmz=="l10")]][1]))
    if ("l10_days" %in% nmz) km_L10 <- suppressWarnings(as.numeric(km_ref[[which(nmz=="l10_days")]][1]))
    km_note <- "KM summary loaded"
  }

  sum_tbl <- data.frame(
    component=comp, clock=clock,
    logLik_Weibull     = if (!is.null(fit_wb))  as.numeric(logLik(fit_wb))  else NA_real_,
    logLik_Exponential = if (!is.null(fit_exp)) as.numeric(logLik(fit_exp)) else NA_real_,
    AIC_Weibull=AIC_wb, AIC_Exponential=AIC_exp,
    dAIC_Exp_minus_WB = if (is.na(AIC_wb) || is.na(AIC_exp)) NA_real_ else AIC_exp - AIC_wb,
    L10_Weibull=L10_wb["L10"],  L10_Weibull_lcl=L10_wb["L10_lcl"], L10_Weibull_ucl=L10_wb["L10_ucl"],
    L10_Exponential=L10_exp["L10"], L10_Exp_lcl=L10_exp["L10_lcl"], L10_Exp_ucl=L10_exp["L10_ucl"],
    KM_L10_ref=km_L10, KM_note=km_note,
    n_events = n_ev, n_rows = nrow(df)
  )
  write.csv(sum_tbl, file.path(OUT_DIR, paste0("para_summary_",comp,"_",clock,".csv")), row.names=FALSE)

  if (!is.null(fit_wb))  plot_model_only(fit_wb,  df, comp, paste0("Weibull_",clock))
  if (!is.null(fit_exp)) plot_model_only(fit_exp, df, comp, paste0("Exponential_",clock))

  # --- Probability (Justification) Plots from KM + PP vs Models ---
  try({
    emit_km_probability_plots(
      df      = df,           # run_component içindeki veri (time,event içermeli)
      comp    = comp,         # "battery" | "ir" | "impact"
      clock   = clock,        # "calendar" vb.
      out_dir = OUT_DIR,      # mevcut çıktı klasörün
      fit_wb  = fit_wb,       # varsa; yoksa NULL
      fit_exp = fit_exp       # varsa; yoksa NULL
    )
  }, silent = TRUE)

  cat("\n===", toupper(comp), "(", clock, ") ===\n")
  print(aic_tbl); print(utils::head(tr_out)); print(sum_tbl)

  # KM vs Weibull overlay (only Weibull)
  if (!is.null(fit_wb)) plot_km_vs_weibull(fit_wb, df, comp, clock, OUT_DIR)
}

for (comp in c("battery","ir","impact")) run_component(comp, clock="calendar")
cat("\nParametric finished — utils.R pipeline reused (v4). Outputs prefixed with para_.\n")
