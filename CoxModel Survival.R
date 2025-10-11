###############################################################
# Cox Proportional Hazards Modeling - DirtSlurper3100 Project #
###############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(survival)
  library(survminer)
  library(broom)
})

# Auto-detect where real data starts (skips README lines)
all_lines <- read_lines("DirtSlurper3100.csv")
hdr_line  <- which(grepl("^Registration date,Total usage time", all_lines))[1]

raw <- read_csv("DirtSlurper3100.csv", skip = hdr_line - 1, show_col_types = FALSE)

# Data preparation 
df <- raw %>%
  mutate(
    reg_date     = dmy(`Registration date`),
    failure_date = suppressWarnings(dmy(na_if(`Failure date`, "---"))),
    obs_date     = if_else(!is.na(failure_date),
                           failure_date,
                           as.Date("2019-12-31")),
    time_days    = as.numeric(obs_date - reg_date),
    pets   = factor(Pets, levels = c("NO", "YES")),
    carpet = as.integer(`Carpet score`)
  ) %>%
  select(time_days, pets, carpet,
         `Battery status`, `IR status`, `Impact status`)

# prepare per-component dataset
mk_component <- function(data, status_col) {
  data %>%
    mutate(event = as.integer(.data[[status_col]] == "Damage")) %>%
    select(time_days, event, pets, carpet) %>%
    filter(!is.na(time_days), !is.na(event))
}

battery_df <- mk_component(df, "Battery status")
ir_df      <- mk_component(df, "IR status")
impact_df  <- mk_component(df, "Impact status")

# Cox model
run_cox <- function(data, label) {
  cat("\n====================\n", label, "\n====================\n")
  fit <- coxph(Surv(time_days, event) ~ pets + carpet, data = data)
  print(summary(fit))
  
  # PH assumption test
  ph <- cox.zph(fit)
  print(ph)
  
  dir.create("results_cox", showWarnings = FALSE)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    rename(HR = estimate, CI_low = conf.low, CI_high = conf.high, p = p.value) %>%
    write.csv(file = paste0("results_cox/HR_", label, ".csv"), row.names = FALSE)
  
  capture.output(summary(fit),
                 file = paste0("results_cox/Summary_", label, ".txt"))
  capture.output(ph,
                 file = paste0("results_cox/PHtest_", label, ".txt"))
  
  # survival curves
  newdat <- expand.grid(pets = levels(data$pets),
                        carpet = mean(data$carpet, na.rm = TRUE))
  sf <- survfit(fit, newdata = newdat)
  p <- ggsurvplot(sf, data = data, conf.int = TRUE,
                  legend.title = "Pets",
                  ggtheme = theme_minimal(),
                  title = paste("Cox PH -", label))
  ggsave(paste0("results_cox/Cox_", label, ".png"), p$plot, width = 7, height = 5, dpi = 300)
}

# run models
run_cox(battery_df, "Battery")
run_cox(ir_df, "IR")
run_cox(impact_df, "Impact")

cat("\nCox PH models complete\n")