##########################################
# survival_utils.R
# Common utility functions for survival analysis
# (Kaplan–Meier, Cox, Log-rank)
# Author: Ege Erdal
# Course: Survival Analysis for Data Scientists (2AMS11)
##########################################

# ------------------------------
# Required Packages
# ------------------------------
if (!require("dplyr")) install.packages("dplyr")
if (!require("stringr")) install.packages("stringr")
if (!require("survival")) install.packages("survival")

library(dplyr)
library(stringr)
library(survival)

# ------------------------------------------------------------
# 1. Read & Clean the DirtSlurper Dataset
# ------------------------------------------------------------
read_dirt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header_idx <- tail(grep("Registration date", lines, ignore.case = TRUE), 1)
  if (length(header_idx) == 0) stop("Header line not found in file.")
  
  dat <- read.csv(
    text = paste(lines[header_idx:length(lines)], collapse = "\n"),
    na.strings = c("---", "N/A", "NA", ""),
    stringsAsFactors = FALSE
  )
  
  # Normalize column names
  names(dat) <- tolower(gsub("[\\.\\s]+", "_", names(dat)))
  
  # Fix known corrupted names
  name_map <- c(
    "regi_tration_date" = "registration_date",
    "total_u_age_time"  = "total_usage_time",
    "pet_"              = "pets",
    "carpet_core"       = "carpet_score",
    "battery_tatu_"     = "battery_status",
    "impact_tatu_"      = "impact_status",
    "ir_tatu_"          = "ir_status"
  )
  for (old in names(name_map)) {
    if (old %in% names(dat)) {
      names(dat)[names(dat) == old] <- name_map[[old]]
    }
  }
  
  # Convert date columns
  if ("registration_date" %in% names(dat)) {
    dat$registration_date <- as.Date(dat$registration_date, "%d%b%Y")
  }
  if ("failure_date" %in% names(dat)) {
    dat$failure_date <- as.Date(dat$failure_date, "%d%b%Y")
  }
  
  # Logical and numeric coercions
  if ("pets" %in% names(dat)) dat$pets <- toupper(trimws(dat$pets)) == "YES"
  if ("sent_for_repair" %in% names(dat)) dat$sent_for_repair <- toupper(trimws(dat$sent_for_repair)) == "YES"
  if ("total_usage_time" %in% names(dat)) dat$total_usage_time <- as.numeric(dat$total_usage_time)
  if ("carpet_score" %in% names(dat)) dat$carpet_score <- as.numeric(dat$carpet_score)
  
  message("✅ Loaded ", nrow(dat), " rows × ", ncol(dat), " columns.")
  return(dat)
}

# ------------------------------------------------------------
# 2. Prepare Survival-Ready Data
# ------------------------------------------------------------
prepare_surv_data <- function(df, component, clock = "calendar", censor_date = as.Date("2019-12-31")) {
  status_col <- grep(paste0("^", component, "_status$"), names(df), ignore.case = TRUE, value = TRUE)
  if (length(status_col) == 0) stop("Status column not found for ", component)
  
  df$event <- toupper(trimws(df[[status_col]])) == "DAMAGE"
  
  if (clock == "calendar") {
    df$end_date <- ifelse(
      !is.na(df$failure_date),
      pmin(df$failure_date, censor_date),
      censor_date
    )
    df$end_date <- as.Date(df$end_date, origin = "1970-01-01")
    df$time <- as.numeric(df$end_date - df$registration_date)
  } else {
    df$time <- as.numeric(df$total_usage_time)
  }
  
  df <- df %>% filter(!is.na(time) & time >= 0)
  return(df)
}

# ------------------------------------------------------------
# 3. Create a Surv Object
# ------------------------------------------------------------
create_surv <- function(df) {
  return(Surv(df$time, df$event))
}
