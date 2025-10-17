# ------------------------------
# 0. Libraries & Utilities
# ------------------------------
if (!require("survminer")) install.packages("survminer")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("gridExtra")) install.packages("gridExtra")

library(survminer)
library(ggplot2)
library(survival)
library(gridExtra)

source("utils.R")

# ------------------------------
# 1. Load and Prepare Data
# ------------------------------
df <- read_dirt("DirtSlurper3100.csv")
component_list <- c("battery", "ir", "impact")

# ------------------------------
# 2. Diagnostic Function
# ------------------------------
run_diagnostics <- function(component, clock = "calendar") {
  message("\n===== Component: ", toupper(component), " =====")
  
  # --- Prepare data ---
  data <- prepare_surv_data(df, component, clock)
  surv_obj <- create_surv(data)
  
  # --- Fit Kaplan–Meier ---
  km_fit <- survfit(surv_obj ~ pets, data = data)
  
  # ------------------------------
  # A. Kaplan–Meier and Log–log Plots
  # ------------------------------
  p1 <- ggsurvplot(
    km_fit,
    data = data,
    conf.int = TRUE,
    risk.table = TRUE,
    palette = c("#1f77b4", "#ff7f0e"),
    ggtheme = theme_minimal(),
    title = paste("Kaplan–Meier Curve –", toupper(component))
  )
  
  p2 <- ggsurvplot(
    km_fit,
    data = data,
    fun = "cloglog",  # log(-log(S(t))) vs log(t)
    palette = c("#1f77b4", "#ff7f0e"),
    ggtheme = theme_minimal(),
    title = paste("Log–log Survival Plot –", toupper(component))
  )
  
  # ------------------------------
  # B. Cox Model and PH Diagnostics
  # ------------------------------
  cox_fit <- coxph(surv_obj ~ pets + carpet_score, data = data)
  schoenfeld <- cox.zph(cox_fit)
  global_p <- round(schoenfeld$table["GLOBAL", "p"], 4)
  message("Global PH test p-value: ", global_p)
  
  p3 <- ggcoxzph(schoenfeld) +
    ggtitle(paste("Schoenfeld Residuals –", toupper(component))) +
    theme_minimal()
  
  # ------------------------------
  # C. Cox–Snell Residuals (Global Fit Check)
  # ------------------------------
  fit_data <- model.frame(cox_fit)  # ensures same subset
  cs_resid <- -residuals(cox_fit, type = "martingale")
  
  cs_fit <- survfit(Surv(cs_resid, fit_data$event) ~ 1)
  cs_data <- data.frame(
    x = cs_fit$time,
    y = -log(cs_fit$surv)
  )
  
  p4 <- ggplot(cs_data, aes(x = x, y = y)) +
    geom_line(color = "#2ca02c", linewidth = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    labs(
      title = paste("Cox–Snell Residuals –", toupper(component)),
      x = "Cox–Snell residuals",
      y = "Cumulative hazard"
    ) +
    theme_minimal()
  
  # ------------------------------
  # D. Weibull AFT Model + QQ Plot
  # ------------------------------
  aft_fit <- survreg(Surv(time, event) ~ pets + carpet_score, data = data, dist = "weibull")
  
  res <- residuals(aft_fit, type = "working")
  res_df <- data.frame(res = res[!is.na(res)])  
  
  p5 <- ggplot(res_df, aes(sample = res)) +
    stat_qq(color = "#9467bd") +
    stat_qq_line(linetype = "dashed", color = "gray40") +
    labs(
      title = paste("Weibull QQ Plot –", toupper(component)),
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal()
  
  # ------------------------------
  # 3. Save Plots
  # ------------------------------
  dir.create("figures/diagnostics", showWarnings = FALSE, recursive = TRUE)
  
  ggsave(paste0("figures/diagnostics/KM_", component, ".png"),
         p1$plot, width = 7, height = 5, dpi = 300)
  ggsave(paste0("figures/diagnostics/loglog_", component, ".png"),
         p2$plot, width = 7, height = 5, dpi = 300)
  ggsave(paste0("figures/diagnostics/schoenfeld_", component, ".png"),
         p3, width = 7, height = 5, dpi = 300)
  ggsave(paste0("figures/diagnostics/coxsnell_", component, ".png"),
         p4, width = 7, height = 5, dpi = 300)
  ggsave(paste0("figures/diagnostics/weibullQQ_", component, ".png"),
         p5, width = 7, height = 5, dpi = 300)
  
}

# ------------------------------
# 4. Run Diagnostics for Each Component
# ------------------------------
for (comp in component_list) {
  run_diagnostics(comp)
}

