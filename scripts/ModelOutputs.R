################################################################################
# Script: ModelOutputs.R
# Author: Isaac Harris
#
# Purpose:
#   Helper functions to turn mixed models (glmmTMB / lme4) into nicely formatted
#   gt tables, and export them (e.g. as .docx) for reporting / supplements.
#
# Requirements:
#   - Models already fitted in the environment:
#       P0_age, F1_age, F3_age
#       P0_LRS, F1_LRS, F3_LRS
#       P0_fit, F1_fit, F3_fit
#       P0_sur, F1_sur, F3_sur
#       transgen_fit, transgen_LRS
################################################################################

library(here)
library(gt)
library(broom.mixed)
library(car)
library(glmmTMB)
library(lme4)

################################################################################
# Helper: Mixed model table (GLMM / LMM)                                       #
################################################################################
make_mixed_model_table <- function(
    model,
    title = "Mixed Model Summary",
    subtitle = NULL,
    anova_type = 2,
    pairwise_target = NULL
) {
  stopifnot(requireNamespace("gt",          quietly = TRUE))
  stopifnot(requireNamespace("broom.mixed", quietly = TRUE))
  stopifnot(requireNamespace("car",         quietly = TRUE))
  stopifnot(requireNamespace("emmeans",     quietly = TRUE))
  
  is_lme4    <- inherits(model, c("lmerMod", "glmerMod"))
  is_glmmTMB <- inherits(model, "glmmTMB")
  
  ## --- 1. Formulas -------------------------------------------------
  model_formula <- tryCatch(
    paste(deparse(formula(model)), collapse = " "),
    error = function(e) ""
  )
  if (is_glmmTMB) {
    model_formula <- paste(deparse(formula(model, component = "cond")), collapse = " ")
  }
  
  ## --- 2. Fixed effects with CI ------------------------------------
  fx <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE, conf.level = 0.95)
  if (!"df" %in% names(fx)) fx$df <- NA_real_
  
  fx$p.value_fmt <- ifelse(fx$p.value < 0.001, "<0.001", sprintf("%.3f", fx$p.value))
  fx$ci_fmt      <- sprintf("[%0.2f, %0.2f]", fx$conf.low, fx$conf.high)
  
  # Format numeric columns
  fx$estimate  <- round(fx$estimate, 2)
  fx$std.error <- round(fx$std.error, 2)
  fx$statistic <- round(fx$statistic, 2)
  fx$df_fmt    <- ifelse(is.na(fx$df), "", as.character(round(fx$df, 1)))
  
  predictors_block <- data.frame(
    Block        = "Predictors",
    Term         = fx$term,
    Estimate     = as.character(fx$estimate),
    StdError     = as.character(fx$std.error),
    CI           = fx$ci_fmt,
    df           = fx$df_fmt,
    Statistic    = as.character(fx$statistic),
    p            = fx$p.value_fmt,
    Random_Value = NA_character_,
    stringsAsFactors = FALSE
  )
  
  ## --- 3. Random effects -------------------------------------------
  rand_rows_var <- NULL
  rand_rows_n   <- NULL
  residual_row  <- NULL
  
  if (is_lme4) {
    rv <- try(broom.mixed::tidy(model, effects = "ran_vars"), silent = TRUE)
    if (!inherits(rv, "try-error") && nrow(rv) > 0) {
      rand_rows_var <- data.frame(
        Term  = paste0("τ00 ", rv$group, ifelse(rv$term == "(Intercept)", "", paste0(":", rv$term))),
        Value = rv$estimate, stringsAsFactors = FALSE)
    }
    re_struct <- try(lme4::getME(model, "flist"), silent = TRUE)
    if (!inherits(re_struct, "try-error") && length(re_struct) > 0) {
      n_levels <- vapply(re_struct, nlevels, integer(1))
      rand_rows_n <- data.frame(Term = paste0("N ", names(n_levels)), Value = as.numeric(n_levels), stringsAsFactors = FALSE)
    }
    if (inherits(model, "lmerMod")) {
      residual_row <- data.frame(Term = "Residual", Value = sigma(model)^2, stringsAsFactors = FALSE)
    }
  } else if (is_glmmTMB) {
    vc_df <- try(as.data.frame(VarCorr(model)), silent = TRUE)
    if (!inherits(vc_df, "try-error") && nrow(vc_df) > 0) {
      rand_rows_var <- data.frame(
        Term = paste0("τ00 ", vc_df$grp, ifelse(vc_df$var1 == "(Intercept)", "", paste0(":", vc_df$var1))),
        Value = vc_df$vcov, stringsAsFactors = FALSE)
    }
    rand_terms <- lme4::findbars(formula(model, component = "cond"))
    if (length(rand_terms) > 0) {
      grp_names <- unique(vapply(rand_terms, function(x) as.character(x[[3]]), character(1)))
      # FIX: Removed the stray '+' that was here
      n_list <- lapply(grp_names, function(g) length(unique(model$frame[[g]])))
      rand_rows_n <- data.frame(Term = paste0("N ", grp_names), Value = as.numeric(unlist(n_list)), stringsAsFactors = FALSE)
    }
    if (identical(model$modelInfo$family$family, "gaussian")) {
      residual_row <- data.frame(Term = "Residual", Value = sigma(model)^2, stringsAsFactors = FALSE)
    }
  }
  
  obs_row <- data.frame(Term = "Observations", Value = stats::nobs(model), stringsAsFactors = FALSE)
  random_disp <- do.call(rbind, Filter(Negate(is.null), list(rand_rows_var, rand_rows_n, residual_row, obs_row)))
  
  random_block <- data.frame(
    Block        = "Random Effects",
    Term         = random_disp$Term,
    Estimate     = NA_character_,
    StdError     = NA_character_,
    CI           = NA_character_,
    df           = NA_character_,
    Statistic    = NA_character_,
    p            = NA_character_,
    Random_Value = as.character(round(random_disp$Value, 2)),
    stringsAsFactors = FALSE
  )
  
  ## --- 4. ANOVA block ----------------------------------------------
  aov_tab <- try(car::Anova(model, type = anova_type), silent = TRUE)
  anova_block <- NULL
  if (!inherits(aov_tab, "try-error")) {
    aov_df <- as.data.frame(aov_tab)
    aov_df$Term <- rownames(aov_df)
    stat_col <- if ("F value" %in% names(aov_df)) "F value" else "Chisq"
    p_col    <- if ("Pr(>F)" %in% names(aov_df)) "Pr(>F)" else "Pr(>Chisq)"
    
    anova_block <- data.frame(
      Block        = "ANOVA",
      Term         = aov_df$Term,
      Estimate     = NA_character_,
      StdError     = NA_character_,
      CI           = NA_character_,
      df           = as.character(round(aov_df$Df, 1)),
      Statistic    = paste0(gsub(" value", "", stat_col), " = ", round(aov_df[[stat_col]], 2)),
      p            = ifelse(aov_df[[p_col]] < 0.001, "<0.001", sprintf("%.3f", aov_df[[p_col]])),
      Random_Value = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  ## --- 5. NEW: EMMEANS / CONTRASTS -----------
  contrast_block <- NULL
  
  if (!is.null(pairwise_target)) {
    # 1. Calculate Emmeans & Pairs
    emm <- emmeans::emmeans(model, specs = pairwise_target)
    
    # 2. Robust Summary (Get Estimates + P-values in one table)
    emm_pairs <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey")
    s_pairs   <- as.data.frame(summary(emm_pairs, infer = c(TRUE, TRUE)))
    
    if (nrow(s_pairs) > 0) {
      # 3. Dynamic Column Finding (Prevent 0-length errors)
      
      # ESTIMATE: Try 'estimate' (linear) or 'ratio' (log/logit)
      est_col_name <- intersect(names(s_pairs), c("estimate", "ratio", "emmean"))
      est_vals <- if(length(est_col_name) > 0) s_pairs[[est_col_name[1]]] else s_pairs[[2]]
      
      # STATISTIC: Find any column ending in .ratio (t.ratio, z.ratio, etc)
      stat_col_name <- grep("\\.ratio$", names(s_pairs), value = TRUE)
      if (length(stat_col_name) > 0) {
        stat_vals <- s_pairs[[stat_col_name[1]]]
      } else {
        stat_vals <- rep(NA, nrow(s_pairs)) # Safe fallback
      }
      
      # DF: Handle missing DF
      df_vals <- if("df" %in% names(s_pairs)) s_pairs$df else rep(NA, nrow(s_pairs))
      
      # 4. Format
      est_fmt  <- sprintf("%.2f", est_vals)
      se_fmt   <- sprintf("%.2f", s_pairs$SE)
      ci_fmt   <- sprintf("[%0.2f, %0.2f]", s_pairs$lower.CL, s_pairs$upper.CL)
      stat_fmt <- ifelse(is.na(stat_vals), "", sprintf("%.2f", stat_vals))
      df_fmt   <- ifelse(is.na(df_vals) | is.infinite(df_vals), "", as.character(round(df_vals, 1)))
      p_fmt    <- ifelse(s_pairs$p.value < 0.001, "<0.001", sprintf("%.3f", s_pairs$p.value))
      
      contrast_block <- data.frame(
        Block        = "Post-hoc Contrasts",
        Term         = as.character(s_pairs$contrast),
        Estimate     = est_fmt,
        StdError     = se_fmt,
        CI           = NA,
        df           = df_fmt,
        Statistic    = stat_fmt,
        p            = p_fmt,
        Random_Value = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  
  ## --- 6. Combine and Table ----------------------------------------
  combo <- rbind(predictors_block, random_block, anova_block, contrast_block)
  
  tbl <- gt::gt(combo) |>
    gt::cols_label(
      Block = "", Term = "", Estimate = "Estimate", StdError = "Std. Error",
      CI = "95% CI", df = "df", Statistic = "Statistic", p = "p", Random_Value = "Random Effects"
    ) |>
    gt::tab_spanner(label = "Fixed Effects", columns = c(Estimate, StdError, CI, df, Statistic, p)) |>
    gt::tab_style(style = gt::cell_text(weight = "bold"), locations = gt::cells_body(columns = "Block")) |>
    gt::fmt_missing(columns = gt::everything(), missing_text = "") |>
    gt::tab_options(table.font.size = gt::px(12), data_row.padding = gt::px(2)) |>
    gt::tab_header(title = title, subtitle = if(is.null(subtitle)) model_formula else subtitle)
  
  return(tbl)
}


################################################################################
# Helper: Survival model table                                                 #
################################################################################

# This version is essentially the same layout, just with slightly different
# grouping-factor handling to deal nicely with nested random effects like Plate/Worm.

make_survival_table <- function(
    model,
    title = "Mixed Model Summary",
    subtitle = NULL,
    anova_type = 2  # 2 = Type II, 3 = Type III
) {
  stopifnot(requireNamespace("gt", quietly = TRUE))
  stopifnot(requireNamespace("broom.mixed", quietly = TRUE))
  stopifnot(requireNamespace("car", quietly = TRUE))
  
  is_lme4 <- inherits(model, c("lmerMod", "glmerMod"))
  is_glmmTMB <- inherits(model, "glmmTMB")
  
  ## --- helper: pull grouping factor names from random effects ---
  extract_grouping_factors <- function(formula_obj) {
    bars <- lme4::findbars(formula_obj)
    if (length(bars) == 0) return(character(0))
    
    raw_terms <- lapply(bars, function(b) b[[3]])
    groups <- unlist(lapply(raw_terms, function(expr) {
      ch <- as.character(expr)
      grp_raw <- paste(ch, collapse = "")
      parts <- strsplit(grp_raw, "/", fixed = TRUE)[[1]]
      
      if (length(parts) == 1) return(parts)
      acc <- parts[1]
      out <- acc
      if (length(parts) > 1) {
        for (i in 2:length(parts)) {
          acc <- paste0(acc, ":", parts[i])
          out <- c(out, acc)
        }
      }
      out
    }))
    unique(groups)
  }
  
  ## --- 1. Formulas ---
  model_formula <- tryCatch(
    paste(deparse(formula(model)), collapse = " "),
    error = function(e) ""
  )
  
  disp_formula <- NULL
  zi_formula <- NULL
  
  if (is_glmmTMB) {
    model_formula <- paste(deparse(formula(model, component = "cond")), collapse = " ")
    disp_formula <- tryCatch(paste(deparse(formula(model, component = "disp")), collapse = " "), error = function(e) NULL)
    zi_formula <- tryCatch(paste(deparse(formula(model, component = "zi")), collapse = " "), error = function(e) NULL)
    if (!is.null(disp_formula) && disp_formula == "~1") disp_formula <- NULL
    if (!is.null(zi_formula) && zi_formula == "~0") zi_formula <- NULL
  }
  
  ## --- 2. Fixed effects ---
  fx <- broom.mixed::tidy(model, effects = "fixed")
  if (!"df" %in% names(fx)) fx$df <- NA_real_
  fx$p.value_fmt <- ifelse(fx$p.value < 0.001, "<0.001", sprintf("%.3f", fx$p.value))
  
  fixed_block <- data.frame(
    Block = "Predictors",
    Term = fx$term,
    Estimate = fx$estimate,
    StdError = fx$std.error,
    df = fx$df,
    Statistic = fx$statistic,
    p = fx$p.value_fmt,
    stringsAsFactors = FALSE
  )
  
  ## --- 3. Random effects ---
  rand_rows_var <- NULL
  rand_rows_n <- NULL
  residual_row <- NULL
  
  if (is_lme4) {
    rv <- try(broom.mixed::tidy(model, effects = "ran_vars"), silent = TRUE)
    if (!inherits(rv, "try-error") && nrow(rv) > 0) {
      rand_rows_var <- data.frame(
        Term = paste0("τ00 ", rv$group, ifelse(rv$term == "(Intercept)", "", paste0(":", rv$term))),
        Value = rv$estimate,
        stringsAsFactors = FALSE
      )
    }
    
    re_struct <- try(lme4::getME(model, "flist"), silent = TRUE)
    if (!inherits(re_struct, "try-error") && length(re_struct) > 0) {
      n_levels <- vapply(re_struct, nlevels, integer(1))
      rand_rows_n <- data.frame(
        Term = paste0("N ", names(n_levels)),
        Value = as.numeric(n_levels),
        stringsAsFactors = FALSE
      )
    }
    
    if (inherits(model, "lmerMod")) {
      sigma2 <- sigma(model)^2
      residual_row <- data.frame(Term = "Residual", Value = sigma2, stringsAsFactors = FALSE)
    }
  } else if (is_glmmTMB) {
    vc_df <- try(as.data.frame(VarCorr(model)), silent = TRUE)
    if (!inherits(vc_df, "try-error") && nrow(vc_df) > 0) {
      rand_rows_var <- data.frame(
        Term = paste0("τ00 ", vc_df$grp, ifelse(vc_df$var1 == "(Intercept)", "", paste0(":", vc_df$var1))),
        Value = vc_df$vcov,
        stringsAsFactors = FALSE
      )
      rand_rows_var <- rand_rows_var[!duplicated(rand_rows_var$Term), , drop = FALSE]
    }
    
    grp_names <- extract_grouping_factors(formula(model, component = "cond"))
    if (length(grp_names) > 0) {
      mf <- model$frame
      get_group_vec <- function(gname) {
        if (gname %in% names(mf)) return(mf[[gname]])
        if (grepl(":", gname, fixed = TRUE)) {
          parts <- strsplit(gname, ":", fixed = TRUE)[[1]]
          if (all(parts %in% names(mf))) return(interaction(mf[parts], drop = TRUE))
        }
        rep(NA, nrow(mf))
      }
      n_list <- lapply(grp_names, function(g) {
        gv <- get_group_vec(g)
        length(unique(gv[!is.na(gv)]))
      })
      rand_rows_n <- data.frame(
        Term = paste0("N ", grp_names),
        Value = as.numeric(unlist(n_list)),
        stringsAsFactors = FALSE
      )
    }
    
    fam <- model$modelInfo$family$family
    if (identical(fam, "gaussian")) {
      sigma2 <- sigma(model)^2
      residual_row <- data.frame(Term = "Residual", Value = sigma2, stringsAsFactors = FALSE)
    }
  }
  
  obs_row <- data.frame(Term = "Observations", Value = stats::nobs(model), stringsAsFactors = FALSE)
  
  random_disp <- do.call(
    rbind,
    Filter(Negate(is.null), list(rand_rows_var, rand_rows_n, residual_row, obs_row))
  )
  
  ## --- 4. ANOVA ---
  aov_tab <- try(car::Anova(model, type = anova_type), silent = TRUE)
  anova_block <- NULL
  
  if (!inherits(aov_tab, "try-error")) {
    aov_df <- as.data.frame(aov_tab)
    aov_df$Term <- rownames(aov_df)
    
    if ("F value" %in% names(aov_df)) {
      stat_label <- "F"
      anova_disp <- data.frame(
        Term = aov_df$Term,
        df = aov_df$Df,
        TestStat = aov_df[["F value"]],
        p = aov_df[["Pr(>F)"]],
        stringsAsFactors = FALSE
      )
    } else if ("Chisq" %in% names(aov_df)) {
      stat_label <- "Chisq"
      anova_disp <- data.frame(
        Term = aov_df$Term,
        df = aov_df$Df,
        TestStat = aov_df[["Chisq"]],
        p = aov_df[["Pr(>Chisq)"]],
        stringsAsFactors = FALSE
      )
    } else {
      stat_label <- "Stat"
      p_col <- grep("^Pr", names(aov_df), value = TRUE)
      anova_disp <- data.frame(
        Term = aov_df$Term,
        df = if ("Df" %in% names(aov_df)) aov_df$Df else NA,
        TestStat = NA,
        p = if (length(p_col) > 0) aov_df[[p_col[1]]] else NA,
        stringsAsFactors = FALSE
      )
    }
    
    anova_disp$p <- ifelse(
      is.na(anova_disp$p), "",
      ifelse(anova_disp$p < 0.001, "<0.001", sprintf("%.3f", anova_disp$p))
    )
    anova_disp$TestStat <- round(anova_disp$TestStat, 2)
    anova_disp$df <- round(anova_disp$df, 2)
    
    anova_block <- data.frame(
      Block = "ANOVA",
      Term = anova_disp$Term,
      Estimate = NA,
      StdError = NA,
      CI = NA,
      df = anova_disp$df,
      Statistic = paste0(stat_label, " = ", anova_disp$TestStat),
      p = anova_disp$p,
      Random_Value = NA,
      stringsAsFactors = FALSE
    )
  }
  
  ## --- 5. Format numbers ---
  round_if_num <- function(x, digits = 2) if (is.numeric(x)) round(x, digits) else x
  
  fixed_block$Estimate <- round_if_num(fixed_block$Estimate)
  fixed_block$StdError <- round_if_num(fixed_block$StdError)
  fixed_block$Statistic <- round_if_num(fixed_block$Statistic)
  fixed_block$df <- ifelse(is.na(fixed_block$df), "", round(fixed_block$df, 1))
  
  if (!is.null(random_disp) && nrow(random_disp) > 0) {
    random_disp$Value <- ifelse(is.numeric(random_disp$Value), round(random_disp$Value, 2), random_disp$Value)
  } else {
    random_disp <- data.frame(Term = "Observations", Value = stats::nobs(model), stringsAsFactors = FALSE)
  }
  
  random_block <- data.frame(
    Block = "Random Effects",
    Term = random_disp$Term,
    Estimate = NA,
    StdError = NA,
    CI = NA,
    df = NA,
    Statistic = NA,
    p = NA,
    Random_Value = random_disp$Value,
    stringsAsFactors = FALSE
  )
  
  predictors_block <- data.frame(
    Block = "Predictors",
    CI = NA,
    Term = fixed_block$Term,
    Estimate = fixed_block$Estimate,
    StdError = fixed_block$StdError,
    df = fixed_block$df,
    Statistic = fixed_block$Statistic,
    p = fixed_block$p,
    Random_Value = NA,
    stringsAsFactors = FALSE
  )
  
  combo <- rbind(predictors_block, random_block, anova_block)
  
  ## --- 6. Build gt table ---
  tbl <- gt::gt(combo) |>
    gt::cols_label(
      Block = "",
      Term = "",
      Estimate = "Estimate",
      StdError = "Std. Error",
      df = "df",
      Statistic = "Statistic",
      p = "p",
      Random_Value = "Random Effects"
    ) |>
    gt::tab_spanner(
      label = "Fixed Effects",
      columns = c(Estimate, StdError, df, Statistic, p)
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(columns = "Block")
    ) |>
    gt::fmt_missing(columns = gt::everything(), missing_text = "") |>
    gt::tab_options(
      table.font.size = gt::px(12),
      data_row.padding = gt::px(2)
    )
  
  ## --- 7. Header ---
  if (is.null(subtitle)) {
    subtitle_text <- paste0(
      "Model formula: ", model_formula,
      if (!is.null(disp_formula)) paste0(" | Dispersion: ", disp_formula) else "",
      if (!is.null(zi_formula)) paste0(" | Zero-inflation: ", zi_formula) else ""
    )
  } else {
    subtitle_text <- subtitle
  }
  
  tbl <- tbl |> gt::tab_header(title = title, subtitle = subtitle_text)
  return(tbl)
}

################################################################################
# Export all model tables                                                      #
################################################################################

# Ensure output directory exists
out_dir <- here("model_output")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# P0, F1, F3 – Age-specific reproduction
P0ASR <- make_mixed_model_table(P0_age, title = "P0 Age Specific Reproduction")
gt::gtsave(P0ASR, file = file.path(out_dir, "P0ASR.docx"))

F1ASR <- make_mixed_model_table(F1_age, title = "F1 Age Specific Reproduction")
gt::gtsave(F1ASR, file = file.path(out_dir, "F1ASR.docx"))

F3ASR <- make_mixed_model_table(F3_age, title = "F3 Age Specific Reproduction")
gt::gtsave(F3ASR, file = file.path(out_dir, "F3ASR.docx"))

# Lifetime Reproductive Success (LRS)
P0LRS <- make_mixed_model_table(P0_LRS, title = "P0 Lifetime Reproductive Success", pairwise_target = ~ Treatment)
gt::gtsave(P0LRS, file = file.path(out_dir, "P0LRS.docx"))

F1LRS <- make_mixed_model_table(F1_LRS, title = "F1 Lifetime Reproductive Success", pairwise_target = ~ Treatment)
gt::gtsave(F1LRS, file = file.path(out_dir, "F1LRS.docx"))

F3LRS <- make_mixed_model_table(F3_LRS, title = "F3 Lifetime Reproductive Success", pairwise_target = ~ Treatment)
gt::gtsave(F3LRS, file = file.path(out_dir, "F3LRS.docx"))

# Rate-sensitive fitness

P0fit <- make_mixed_model_table(P0_fit, title = "P0 Rate Sensitive Fitness", pairwise_target = ~ Treatment)
gt::gtsave(P0fit, file = file.path(out_dir, "P0fit.docx"))

F1fit <- make_mixed_model_table(F1_fit, title = "F1 Rate Sensitive Fitness", pairwise_target = ~ Treatment)
gt::gtsave(F1fit, file = file.path(out_dir, "F1fit.docx"))

F3fit <- make_mixed_model_table(F3_fit, title = "F3 Rate Sensitive Fitness", pairwise_target = ~ Treatment)
gt::gtsave(F3fit, file = file.path(out_dir, "F3fit.docx"))

# Survival models
P0sur <- make_survival_table(P0_sur, title = "P0 Survival")
gt::gtsave(P0sur, file = file.path(out_dir, "P0sur.docx"))

F1sur <- make_survival_table(F1_Sur, title = "F1 Survival")
gt::gtsave(F1sur, file = file.path(out_dir, "F1sur.docx"))

F3sur <- make_survival_table(F3_surv, title = "F3 Survival")
gt::gtsave(F3sur, file = file.path(out_dir, "F3sur.docx"))

# Transgen models
transgenfit <- make_mixed_model_table(transgen_fit, title = "transgenerational Rate Sensitive Fitness")
gt::gtsave(transgenfit, file = file.path(out_dir, "transgenfit.docx"))

transgenrepro <- make_mixed_model_table(transgen_LRS, title = "transgenerational LRS")
gt::gtsave(transgenrepro, file = file.path(out_dir, "transgenLRS.docx"))
