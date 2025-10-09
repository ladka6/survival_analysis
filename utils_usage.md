# Survival Utilities (utils.R)

## Overview
This module provides reusable helper functions for survival analysis of the **DirtSlurper3100** dataset.  
It standardizes data loading, cleaning, and right-censoring for all models:
- Kaplan–Meier (nonparametric)
- Cox Proportional Hazards
- Log-Rank test

Using these utilities ensures consistent time/event definitions and reproducibility across analyses.

---

## 🧩 Functions

### `read_dirt(path)`
**Purpose:**  
Reads the DirtSlurper dataset from a `.csv` file, locates the correct header, normalizes column names, fixes encoding, and parses dates.

**Returns:**  
A clean `data.frame` with standardized variables:
- `registration_date`, `failure_date`
- `total_usage_time`, `carpet_score`
- `battery_status`, `impact_status`, `ir_status`
- `pets`, `sent_for_repair`

---

### `prepare_surv_data(df, component, clock = "calendar", censor_date = "2019-12-31")`
**Purpose:**  
Constructs the time-to-event dataset for a given component (`"battery"`, `"ir"`, or `"impact"`).  
Automatically applies right-censoring at the specified date.

**Arguments:**
- `df`: Data frame from `read_dirt()`
- `component`: Component name
- `clock`: `"calendar"` (days) or `"usage"` (hours)
- `censor_date`: Optional override of censoring date

**Returns:**  
A filtered `data.frame` with:
- `time`: duration until failure/censor
- `event`: binary indicator (1 = failure, 0 = censored)

---

### `create_surv(df)`
**Purpose:**  
Creates a `Surv` object from the prepared data.  
Used directly by KM, Cox, and log-rank models.

**Example:**
```r
df <- read_dirt("DirtSlurper3100 (1).csv")
battery_df <- prepare_surv_data(df, "battery")
surv_obj <- create_surv(battery_df)

# Kaplan–Meier
fit_km <- survfit(surv_obj ~ 1, data = battery_df)

# Cox Regression
fit_cox <- coxph(surv_obj ~ pets + carpet_score, data = battery_df)

# Log-Rank Test
survdiff(Surv(time, event) ~ pets, data = battery_df)
