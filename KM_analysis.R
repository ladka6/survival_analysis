
# ------------------------------
# 0. Libraries & Utils
# ------------------------------
if (!require("survminer")) install.packages("survminer")
if (!require("ggplot2")) install.packages("ggplot2")
library(survminer)
library(ggplot2)

# Import shared utilities
source("utils.R")

# ------------------------------
# 1. Load Data
# ------------------------------
df <- read_dirt("DirtSlurper3100.csv")

# Quick overview
cat("\nDataset contains", nrow(df), "observations.\n")

# ------------------------------
# 2. Kaplan–Meier Function
# ------------------------------
run_km <- function(component, clock = "calendar", stratify = FALSE) {
  cat("\n====================================================\n")
  cat("Component:", toupper(component), " | Clock:", clock, "\n")
  cat("====================================================\n")
  
  data <- prepare_surv_data(df, component, clock)
  surv_obj <- create_surv(data)
  
  fit <- if (stratify) {
    survfit(surv_obj ~ pets, data = data)
  } else {
    survfit(surv_obj ~ 1, data = data)
  }
  
  # --- Plot ---
  ggsurv <- ggsurvplot(
    fit,
    data = data,
    conf.int = TRUE,
    risk.table = TRUE,
    ggtheme = theme_minimal(),
    palette = if (stratify) c("steelblue", "firebrick3") else "firebrick3",
    title = paste("Kaplan–Meier:", toupper(component)),
    xlab = ifelse(clock == "calendar", "Time (days)", "Usage (hours)"),
    ylab = "Survival Probability S(t)",
    pval = stratify, 
    pval.method = TRUE
  )
  
  fname <- paste0("KM_", component, "_", clock, ifelse(stratify, "_pets", ""), ".png")
  ggsave(fname, print(ggsurv), width = 7, height = 5.5, dpi = 300)
  print(ggsurv)
  
  if (stratify) {
    cat("\n--- Log-Rank Test ---\n")
    cat(sprintf("Testing for difference in %s survival curves (With Pets vs. No Pets).\n\n", toupper(component)))
    
    log_rank_test <- survdiff(surv_obj ~ pets, data = data)
    print(log_rank_test)
  }
  
  events <- sum(data$event)
  med <- summary(fit)$table["median"]
  quant <- quantile(fit, probs = 0.10)$quantile
  
  cat("\n--- KM Summary ---\n")
  cat(sprintf("N = %d | Events = %d (%.1f%%)\n", nrow(data), events, 100 * events / nrow(data)))
  cat(sprintf("Median survival: %.1f days\n", med))
  cat(sprintf("L10 (S=0.90): %.1f days\n", quant))
  
  write.csv(summary(fit)$table, paste0("KM_summary_", component, "_", clock, ifelse(stratify, "_pets", ""), ".csv"), row.names = TRUE)
}

# ------------------------------
# 3. Run Analyses
# ------------------------------
for (comp in c("battery", "ir", "impact")) run_km(comp)
run_km("battery")
run_km("battery", stratify = TRUE)

run_km("ir")
run_km("ir", stratify = TRUE)

run_km("impact")
run_km("impact", stratify = TRUE)
