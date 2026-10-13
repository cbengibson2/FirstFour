rm(list = ls())

# ============================================================
# unified_analysis.R
#
# Single-pass computation of all three selection criteria:
#   1. Greedy forward AIC selection
#   2. Univariate AIC rank  (free from greedy step-1 candidates)
#   3. Take-one-out from the 4-term baseline model  (TOO-4)
#   4. Take-one-out from the full 12-term model     (TOO-12)
#
# Each network generates one pe object at one n_nonevents value,
# shared across all four analyses.  A single term_groups list is
# the sole source of truth: its keys are the labels stored in
# every output object; its values are the column names in pe.
#
# Outputs
#   stability_results.Rdata          per-network results (greedy + TOO)
#   greedy_results.Rdata             greedy_results list, same schema as
#                                    greedy_results7 (for dw_figures.R etc.)
#   paper/aic_sensitivity_table.tex
#   paper/aic_sensitivity_plot.pdf
#   paper/greedy_aic_curve.pdf
#   paper/greedy_aic_curve_vis.pdf / _val.pdf
#   paper/greedy_vis_heatmap.pdf
#   paper/greedy_valence_heatmap.pdf
#   paper/freq_in_top4.pdf
# ============================================================

library(llrem)
library(survival)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================
# CONFIGURATION  -- edit these for your machine
# ============================================================

PROJ_DIR    <- "/path/to/Stability"       # project root
DATA_DIR    <- "/path/to/DataP"           # root of network data files
N_NONEVENTS <- 120L                       # non-events per step (single value for all analyses)
N_CORES     <- parallel::detectCores() - 1L

# ============================================================
# TERM GROUPS  -- single source of truth
#
# Key   = label stored in every output object (path$added,
#         column headers, table rows, figure y-axis)
# Value = exact column name that llrem writes into pe
#
# Using key == value wherever possible to eliminate any
# label/column translation.
# ============================================================

term_groups <- list(
  Imitation      = "Imitation",
  DyadRecencyLog = "DyadRecencyLog",   # forward-dyad log-time recency
  OutDegSnd      = "OutDegSnd",
  Repetition     = "Repetition",
  SndInitiation  = "SndInitiation",
  RRecInitiation = "RRecInitiation",
  RecRecencyLog  = "RecRecencyLog",
  RecRecencyRank = "RecRecencyRank",
  WeightedRec    = "WeightedRec",
  IOviolation    = "IOviolation",
  NOTPSnd        = "NOTPSnd",
  NITPSnd        = "NITPSnd"
)

term_labels  <- names(term_groups)
term_cols    <- unlist(term_groups, use.names = FALSE)
n_terms      <- length(term_groups)
baseline     <- c("SndInitiation", "IOviolation", "Repetition", "RecRecencyRank")

# ============================================================
# HELPER: locate data file on this machine
#
# greedy_results store the original path under $df.  Strip the
# basename and search the local DATA_DIR subdirectory tree.
# Only needed when re-analysing from a saved greedy_results;
# the main analysis below iterates DATA_DIR directly.
# ============================================================

resolve_path <- function(stored_path) {
  bn   <- basename(stored_path)
  stem <- sub("[.]Rdatadf[.]Rdata$", "", bn, ignore.case = TRUE)
  stem <- sub("[.]rdatadf[.]Rdata$", "", stem, ignore.case = TRUE)
  cands <- list.files(DATA_DIR, pattern = paste0("^", stem, "[.]"),
                      recursive = TRUE, full.names = TRUE,
                      ignore.case = TRUE)
  if (length(cands) == 0L) NA_character_ else cands[1L]
}

# ============================================================
# HELPER: load events from a data file
#
# Supports two formats:
#   poss_edges  -- full risk-set frame; filter occurred == 1
#   edgelist    -- pre-filtered event list
#
# Returns a data.frame(time, sender, receiver) sorted by time,
# and a list(cat, visible) metadata, or NULL on failure.
# ============================================================

load_events <- function(dat_path) {
  e <- new.env()
  tryCatch(load(dat_path, envir = e), error = function(err) NULL)

  if (exists("poss_edges", envir = e)) {
    pe_obj <- e$poss_edges
    el <- pe_obj[pe_obj$occurred == 1L,
                 c("time2", "sender", "reciever")]   # note original typo
    colnames(el) <- c("time", "sender", "receiver")
    cat_val <- if (!is.null(pe_obj$cat))     pe_obj$cat[1]     else NA_character_
    vis_val <- if (!is.null(pe_obj$visible)) pe_obj$visible[1] else NA_character_
  } else if (exists("edgelist", envir = e)) {
    el <- as.data.frame(e$edgelist)[, 1:3]
    colnames(el) <- c("time", "sender", "receiver")
    cat_val <- NA_character_
    vis_val <- NA_character_
  } else {
    return(NULL)
  }

  el <- el[complete.cases(el), ]
  el <- el[order(el$time), ]
  if (nrow(el) < 10L) return(NULL)

  # Path-based category/visibility inference (fallback when not in file)
  p <- tolower(dat_path)
  if (is.na(cat_val) || cat_val == "") {
    cat_val <- dplyr::case_when(
      grepl("gang",              p) ~ "Gang Violence",
      grepl("ftf|glance",        p) ~ "Glances",
      grepl("war|conflict",      p) ~ "International Conflict",
      grepl("crime",             p) ~ "Phone (Crime)",
      grepl("phone|dod",         p) ~ "Phone (DoD)",
      grepl("email",             p) ~ "Emails",
      grepl("fb|uci|\\bpm\\b",   p) ~ "PMs",
      grepl("sms",               p) ~ "SMS",
      grepl("twitter",           p) ~ "Twitter",
      grepl("class",             p) ~ "Classroom",
      TRUE                          ~ tools::file_path_sans_ext(basename(dat_path))
    )
  }
  if (is.na(vis_val) || vis_val == "") {
    vis_val <- if (cat_val %in% c("Gang Violence",
                                   "International Conflict",
                                   "Phone (DoD)")) "Y" else "N"
  }
  # Phone (DoD) override: stored value may be wrong in source files
  if (cat_val == "Phone (DoD)") vis_val <- "Y"

  list(el = el, cat = cat_val, visible = vis_val)
}

# ============================================================
# HELPER: build pe from an event list
# ============================================================

build_pe <- function(el, n_nonevents) {
  prep         <- prepare_edgelist(el)
  all_nodes    <- unique(c(prep$sender, prep$receiver))
  base_riskset <- llrem:::get_base_riskset(NULL, all_nodes)
  pe <- llrem:::remstats_generate(prep,
                                  n            = length(all_nodes),
                                  n_nonevents  = n_nonevents,
                                  base_riskset = base_riskset)
  llrem:::add_survival_cols(pe)
}

# ============================================================
# HELPER: fit a Cox model and return AIC (NA on failure)
# ============================================================

fit_aic <- function(pe, cols) {
  f <- as.formula(paste0(
    "Surv(time1, occurred) ~ ",
    paste(cols, collapse = " + "),
    " + strata(time2)"
  ))
  m <- tryCatch(survival::coxph(f, data = pe), error = function(e) NULL)
  if (is.null(m)) NA_real_ else AIC(m)
}

# ============================================================
# HELPER: greedy forward selection (returns path list)
# ============================================================

greedy_forward <- function(pe, term_groups) {
  null_fit <- tryCatch(
    survival::coxph(Surv(time1, occurred) ~ strata(time2), data = pe),
    error = function(e) NULL
  )
  if (is.null(null_fit)) return(NULL)

  null_aic     <- AIC(null_fit)
  current_aic  <- null_aic
  current_cols <- character(0)
  remaining    <- term_groups
  path         <- vector("list", length(term_groups))

  for (step in seq_along(term_groups)) {
    candidates <- vapply(names(remaining), function(lbl) {
      cols <- c(current_cols, remaining[[lbl]])
      fit_aic(pe, cols)
    }, numeric(1))

    if (all(is.na(candidates))) {
      path[[step]] <- list(step = step, added = NA_character_,
                           aic = NA_real_, delta_aic = NA_real_,
                           prop_delta_aic = NA_real_,
                           candidates = candidates)
      break
    }

    best_lbl <- names(which.min(candidates))
    best_aic <- candidates[[best_lbl]]

    path[[step]] <- list(
      step           = step,
      added          = best_lbl,
      aic            = best_aic,
      delta_aic      = best_aic - current_aic,
      prop_delta_aic = (best_aic - current_aic) / abs(current_aic),
      candidates     = candidates
    )

    current_cols <- c(current_cols, remaining[[best_lbl]])
    remaining    <- remaining[names(remaining) != best_lbl]
    current_aic  <- best_aic
  }

  path
}

# ============================================================
# PER-NETWORK WORKER
#
# Generates pe ONCE per network, then runs:
#   1. Greedy forward selection
#   2. TOO from 4-term baseline model   (if all 4 baseline in top-4)
#   3. TOO from full 12-term model
#
# Univariate ranks are read from path[[1]]$candidates (free).
# ============================================================

analyze_one_network <- function(dat_path, term_groups, term_labels,
                                term_cols, baseline, n_nonevents) {
  na4  <- setNames(rep(NA_real_, length(baseline)),    baseline)
  na12 <- setNames(rep(NA_real_, length(term_labels)), term_labels)
  null_result <- list(
    df = dat_path, cat = NA_character_, visible = NA_character_,
    n = NA_integer_, events = NA_integer_,
    path = NULL, too4 = na4, too12 = na12
  )

  # 1. Load events
  ev <- tryCatch(load_events(dat_path), error = function(e) NULL)
  if (is.null(ev)) return(null_result)

  # 2. Build pe (one call, shared by all analyses)
  pe <- tryCatch(build_pe(ev$el, n_nonevents), error = function(e) NULL)
  if (is.null(pe)) return(null_result)

  # 3. Greedy forward selection
  path <- tryCatch(greedy_forward(pe, term_groups), error = function(e) NULL)
  if (is.null(path)) return(null_result)

  # 4. TOO from 4-term baseline model
  top4_labels <- sapply(path[1:4], function(s) s$added)
  top4_labels <- top4_labels[!is.na(top4_labels)]
  too4 <- if (length(top4_labels) == 4L && all(baseline %in% top4_labels)) {
    top4_cols <- unlist(term_groups[top4_labels], use.names = FALSE)
    aic4 <- fit_aic(pe, top4_cols)
    if (is.na(aic4)) na4 else {
      sapply(setNames(baseline, baseline), function(lbl) {
        drop_cols <- top4_cols[top4_cols != term_groups[[lbl]]]
        fit_aic(pe, drop_cols) - aic4
      })
    }
  } else na4

  # 5. TOO from full 12-term model
  aic12 <- fit_aic(pe, term_cols)
  too12 <- if (is.na(aic12)) na12 else {
    setNames(
      sapply(seq_along(term_labels), function(i)
        fit_aic(pe, term_cols[-i]) - aic12),
      term_labels
    )
  }

  list(
    df      = dat_path,
    cat     = ev$cat,
    visible = ev$visible,
    n       = max(c(ev$el$sender, ev$el$receiver)),
    events  = nrow(ev$el),
    path    = path,
    too4    = too4,
    too12   = too12
  )
}

# ============================================================
# PARALLEL RUN
# ============================================================

datafiles <- list.files(DATA_DIR, pattern = "\\.[Rr]data$",
                        recursive = TRUE, full.names = TRUE)
datafiles <- datafiles[!grepl("finch", datafiles, ignore.case = TRUE)]
cat(sprintf("Found %d data files in %s\n", length(datafiles), DATA_DIR))

cl <- parallel::makeCluster(N_CORES)
parallel::clusterExport(cl,
  c("term_groups", "term_labels", "term_cols", "baseline",
    "n_terms", "N_NONEVENTS",
    "load_events", "build_pe", "fit_aic", "greedy_forward",
    "analyze_one_network", "DATA_DIR"),
  envir = environment())
parallel::clusterCall(cl, function(pkgs)
  lapply(pkgs, library, character.only = TRUE),
  c("llrem", "survival", "dplyr"))

cat(sprintf("Running analysis on %d cores (n_nonevents = %d)...\n",
            N_CORES, N_NONEVENTS))

stability_results <- tryCatch(
  parallel::parLapply(cl, datafiles, function(dp)
    analyze_one_network(dp, term_groups, term_labels,
                        term_cols, baseline, N_NONEVENTS)),
  finally = parallel::stopCluster(cl)
)

# Drop failures
ok <- !sapply(stability_results, function(x) is.null(x$path))
stability_results <- stability_results[ok]
cat(sprintf("%d / %d networks successfully processed.\n",
            length(stability_results), length(datafiles)))

# ============================================================
# SAVE
# ============================================================

PROJ_DIR_OUT <- PROJ_DIR   # keep reference before any reassignment
save(stability_results, term_groups, term_labels, term_cols, baseline,
     file = file.path(PROJ_DIR_OUT, "stability_results.Rdata"))
cat("Saved: stability_results.Rdata\n")

# Also save a backward-compatible greedy_results object for dw_figures.R
greedy_results <- lapply(stability_results, function(x)
  list(df = x$df, cat = x$cat, visible = x$visible,
       n = x$n, events = x$events, path = x$path))
save(greedy_results,
     file = file.path(PROJ_DIR_OUT, "greedy_results.Rdata"))
cat("Saved: greedy_results.Rdata  (backward-compatible)\n")

# ============================================================
# BUILD ANALYSIS MATRICES
# ============================================================

cats   <- sapply(stability_results, function(x) x$cat)
vis    <- sapply(stability_results, function(x) x$visible)
vis[cats == "Phone (DoD)"] <- "Y"
vis_bin   <- ifelse(vis == "Y", "High", "Low")
val_cats  <- c("Gang Violence", "International Conflict")
val_group <- ifelse(cats %in% val_cats, "Conflictual", "Affiliative")
domains   <- sort(unique(cats))
n_dom     <- length(domains)
n_nets    <- length(stability_results)
n_steps   <- n_terms

order_mat <- do.call(rbind, lapply(stability_results, function(x)
  sapply(x$path, function(s) s$added)))

step_prop_mat <- do.call(rbind, lapply(stability_results, function(x)
  sapply(x$path, function(s) as.numeric(s$prop_delta_aic))))

cum_prop_mat <- t(apply(step_prop_mat, 1, cumsum))

too4_mat <- do.call(rbind, lapply(stability_results, function(x) x$too4))
too12_mat <- do.call(rbind, lapply(stability_results, function(x) x$too12))

# ============================================================
# DOMAIN-WEIGHTED HELPERS
# ============================================================

domain_idx <- setNames(
  lapply(domains, function(d) which(cats == d)),
  domains)

dw_top4 <- function(dom_idx_list) {
  dm <- do.call(rbind, lapply(dom_idx_list, function(idx)
    sapply(term_labels, function(lbl)
      mean(apply(order_mat[idx, 1:4, drop = FALSE], 1,
                 function(r) lbl %in% r), na.rm = TRUE))))
  colMeans(dm)
}

dw_step_stats <- function(dom_idx_list) {
  mat <- do.call(rbind, lapply(dom_idx_list, function(idx)
    colMeans(step_prop_mat[idx, , drop = FALSE], na.rm = TRUE)))
  list(mean = colMeans(mat),
       se   = apply(mat, 2, sd) / sqrt(nrow(mat)))
}

# ============================================================
# ANALYSIS 1: UNIVARIATE RANK  (from greedy step-1 candidates)
# ============================================================

univ_ranks_mat <- do.call(rbind, lapply(stability_results, function(x) {
  cands <- x$path[[1]]$candidates
  cands <- cands[!is.na(cands)]
  rnks  <- rank(cands, ties.method = "min")
  out   <- setNames(rep(NA_real_, n_terms), term_labels)
  out[names(rnks)] <- rnks
  out
}))

# Domain-weighted univariate top-4
dw_univ_top4 <- colMeans(
  do.call(rbind, lapply(domain_idx, function(idx)
    colMeans(univ_ranks_mat[idx, , drop = FALSE] <= 4L, na.rm = TRUE))),
  na.rm = TRUE)

cat("\n=== UNIVARIATE TOP-4 (domain-weighted) ===\n")
print(round(sort(dw_univ_top4, decreasing = TRUE), 3))

# ============================================================
# ANALYSIS 2: GREEDY FREQ-IN-TOP-4 (domain-weighted)
# ============================================================

dw_greedy_top4 <- dw_top4(domain_idx)

cat("\n=== GREEDY FREQ-IN-TOP-4 (domain-weighted) ===\n")
print(round(sort(dw_greedy_top4, decreasing = TRUE), 3))

# ============================================================
# ANALYSIS 3: TOO FROM 4-TERM BASELINE MODEL
# ============================================================

mean_too4 <- colMeans(too4_mat, na.rm = TRUE)
n_too4    <- sum(rowSums(!is.na(too4_mat)) > 0)

cat(sprintf("\n=== TOO-4 LOSS (n = %d networks where all baseline in top-4) ===\n",
            n_too4))
print(round(sort(mean_too4, decreasing = TRUE), 2))

# ============================================================
# ANALYSIS 4: TOO FROM FULL 12-TERM MODEL
# ============================================================

mean_too12 <- colMeans(too12_mat, na.rm = TRUE)
n_too12    <- sum(!is.na(too12_mat[, 1]))

cat(sprintf("\n=== TOO-12 LOSS (n = %d networks) ===\n", n_too12))
print(round(sort(mean_too12, decreasing = TRUE), 2))

# ============================================================
# DOMAIN-WEIGHTED PER-STEP AIC
# ============================================================

dw_st <- dw_step_stats(domain_idx)
dw_mean <- dw_st$mean
dw_se   <- dw_st$se
dw_cum  <- cumsum(dw_mean)

cat("\n=== DOMAIN-WEIGHTED PER-STEP AIC ===\n")
for (s in seq_len(n_steps))
  cat(sprintf("  Step %2d: mean=%+.6f  SE=%.6f  cum=%+.6f\n",
              s, dw_mean[s], dw_se[s], dw_cum[s]))
cat(sprintf("\nPeak: %.6f at step %d\n", min(dw_cum), which.min(dw_cum)))
cat(sprintf("Step-4 recovery: %.1f%%\n", 100 * dw_cum[4] / min(dw_cum)))

# ============================================================
# COMBINED COMPARISON TABLE (console)
# ============================================================

hi_domains <- sort(unique(cats[vis_bin   == "High"]))
lo_domains <- sort(unique(cats[vis_bin   == "Low"]))

hi_top4 <- dw_top4(setNames(lapply(hi_domains, function(d) which(cats == d)), hi_domains))
lo_top4 <- dw_top4(setNames(lapply(lo_domains, function(d) which(cats == d)), lo_domains))

comparison <- data.frame(
  Term         = term_labels,
  Greedy_dw    = round(dw_greedy_top4[term_labels], 3),
  Univ_dw      = round(dw_univ_top4[term_labels],   3),
  Greedy_hi    = round(hi_top4[term_labels],         3),
  Greedy_lo    = round(lo_top4[term_labels],         3),
  TOO4_loss    = round(mean_too4[term_labels],        2),
  TOO12_loss   = round(mean_too12[term_labels],       2),
  row.names    = NULL
)
comparison <- comparison[order(-comparison$Greedy_dw), ]

cat("\n=== COMBINED COMPARISON (domain-weighted) ===\n")
print(comparison, row.names = FALSE)

# ============================================================
# LATEX TABLE  (aic_sensitivity_table.tex)
# ============================================================

PAPER_DIR <- file.path(PROJ_DIR_OUT, "paper")

above <- comparison[comparison$Term %in% baseline, ]
below <- comparison[!comparison$Term %in% baseline, ]

tex <- c(
  "\\begin{table}[H]",
  "\\centering", "\\small",
  paste0(
    "\\caption{Four selection criteria for the robust endogenous baseline.",
    " \\textit{Greedy top-4}: domain-weighted proportion of domains in which",
    " each term appears in greedy steps~1--4 ($N=", n_nets, "$ networks,",
    " $", n_dom, "$ domains).",
    " \\textit{Univariate top-4}: domain-weighted proportion of domains in",
    " which the term ranks in the top~4 by step-1 univariate AIC (no path",
    " dependence).",
    " \\textit{TOO-4 loss}: mean AIC increase when the term is dropped from",
    " the joint 4-term baseline model ($n=", n_too4, "$ networks where all",
    " four baseline terms were in the greedy top-4; baseline terms only).",
    " \\textit{TOO-12 loss}: mean AIC increase when the term is dropped from",
    " the full 12-term model ($n=", n_too12, "$).",
    " All analyses use $n_{\\text{nonevents}}=", N_NONEVENTS, "$.}"
  ),
  "\\label{tab:aic_sensitivity}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  paste("\\textbf{Term} & \\textbf{Greedy top-4} & \\textbf{Univ.\\ top-4}",
        "& \\textbf{TOO-4 loss} & \\textbf{TOO-12 loss} \\\\"),
  "\\midrule"
)

fmt_row <- function(r) {
  too4  <- if (is.na(r$TOO4_loss))  "---" else sprintf("%.1f", r$TOO4_loss)
  too12 <- if (is.na(r$TOO12_loss)) "---" else sprintf("%.1f", r$TOO12_loss)
  sprintf("  %s & %.0f\\%% & %.0f\\%% & %s & %s \\\\",
          r$Term, r$Greedy_dw * 100, r$Univ_dw * 100, too4, too12)
}

for (i in seq_len(nrow(above))) tex <- c(tex, fmt_row(above[i, ]))
tex <- c(tex, "\\midrule")
for (i in seq_len(nrow(below))) tex <- c(tex, fmt_row(below[i, ]))
tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")

writeLines(tex, file.path(PAPER_DIR, "aic_sensitivity_table.tex"))
cat("Saved: paper/aic_sensitivity_table.tex\n")

# ============================================================
# FIGURES
# ============================================================

# -- Shared themes -------------------------------------------

curve_theme <- theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

heatmap_theme <- theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 8),
        legend.position = "right")

# -- Term ordering for heatmap y-axis (by mean selection step) ---

all_mss <- rowMeans(
  do.call(rbind, lapply(domain_idx, function(idx) {
    fm <- matrix(0, nrow = n_terms, ncol = n_steps,
                 dimnames = list(term_labels, NULL))
    for (i in idx)
      for (s in seq_len(n_steps)) {
        g <- stability_results[[i]]$path[[s]]$added
        if (!is.na(g) && g %in% term_labels)
          fm[g, s] <- fm[g, s] + 1
      }
    fm <- fm / length(idx)
    rowSums(fm * matrix(seq_len(n_steps), nrow = n_terms,
                        ncol = n_steps, byrow = TRUE))
  })),
  na.rm = TRUE)
term_ord <- names(sort(all_mss))

# -- Figure: Overall elbow curve (base R) --------------------

ov_st <- dw_step_stats(domain_idx)
xs <- seq_len(n_steps)

pdf(file.path(PAPER_DIR, "greedy_aic_curve.pdf"), height = 5, width = 8)
par(mar = c(4, 4.5, 3, 1))
plot(xs, ov_st$mean,
     type = "b", pch = 19,
     xlab = "Step",
     ylab = expression("Proportional " * Delta * "AIC (per step)"),
     main = "Greedy Forward Selection: Marginal AIC Gain per Step (domain-weighted)",
     xaxt = "n",
     ylim = range(c(ov_st$mean - 1.96 * ov_st$se,
                    ov_st$mean + 1.96 * ov_st$se, 0), na.rm = TRUE))
polygon(c(xs, rev(xs)),
        c(ov_st$mean + 1.96 * ov_st$se, rev(ov_st$mean - 1.96 * ov_st$se)),
        col = adjustcolor("steelblue", alpha.f = 0.2), border = NA)
axis(1, at = xs, labels = xs)
abline(h = 0, lty = 2, col = "gray60")
dev.off()
cat("Saved: paper/greedy_aic_curve.pdf\n")

# -- Helpers for split curve and heatmap figures -------------

make_curve_data_dw <- function(dom_idx_list, label) {
  st <- dw_step_stats(dom_idx_list)
  data.frame(step = seq_len(n_steps), mean = st$mean, se = st$se, group = label)
}

make_curve_plot <- function(df) {
  ggplot(df, aes(x = step, y = mean, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                alpha = 0.15, colour = NA) +
    geom_line() + geom_point(size = 1.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
    annotate("rect", xmin = 4.5, xmax = n_steps + 0.5,
             ymin = -Inf, ymax = Inf, fill = "gray90", alpha = 0.4) +
    scale_x_continuous(breaks = seq_len(n_steps)) +
    labs(x = "Greedy step",
         y = expression("Proportional " * Delta * "AIC"),
         colour = NULL, fill = NULL) +
    curve_theme
}

make_long_heatmap_dw <- function(dom_idx_list, facet_label) {
  # Domain-weighted frequency matrix: term x step
  n_doms <- length(dom_idx_list)
  fm_sum <- matrix(0, nrow = n_terms, ncol = n_steps,
                   dimnames = list(term_labels, paste0("Step", seq_len(n_steps))))
  for (idx in dom_idx_list) {
    fm <- matrix(0, nrow = n_terms, ncol = n_steps,
                 dimnames = list(term_labels, paste0("Step", seq_len(n_steps))))
    for (i in idx)
      for (s in seq_len(n_steps)) {
        g <- stability_results[[i]]$path[[s]]$added
        if (!is.na(g) && g %in% term_labels)
          fm[g, s] <- fm[g, s] + 1
      }
    fm_sum <- fm_sum + fm / length(idx)
  }
  fm <- fm_sum / n_doms
  df <- as.data.frame(fm)
  df$term <- rownames(fm)
  pivot_longer(df, -term, names_to = "Step", values_to = "freq") |>
    mutate(
      step_num = as.integer(sub("Step", "", Step)),
      group    = facet_label,
      label    = ifelse(freq >= 0.05, sprintf("%.0f", freq * 100), "")
    )
}

# -- Split curve and heatmap figures -------------------------

hi_idx_list <- lapply(hi_domains, function(d) which(cats == d))
lo_idx_list <- lapply(lo_domains, function(d) which(cats == d))
af_domains  <- sort(unique(cats[val_group == "Affiliative"]))
co_domains  <- sort(unique(cats[val_group == "Conflictual"]))
af_idx_list <- lapply(af_domains, function(d) which(cats == d))
co_idx_list <- lapply(co_domains, function(d) which(cats == d))

lbl_hi <- sprintf("High-visibility (%d domains)", length(hi_domains))
lbl_lo <- sprintf("Low-visibility (%d domains)",  length(lo_domains))
lbl_af <- sprintf("Affiliative (%d domains)",     length(af_domains))
lbl_co <- sprintf("Conflictual (%d domains)",     length(co_domains))

df_vis_curve <- bind_rows(
  make_curve_data_dw(hi_idx_list, lbl_hi),
  make_curve_data_dw(lo_idx_list, lbl_lo))
df_val_curve <- bind_rows(
  make_curve_data_dw(af_idx_list, lbl_af),
  make_curve_data_dw(co_idx_list, lbl_co))

ggsave(file.path(PAPER_DIR, "greedy_aic_curve_vis.pdf"),
       make_curve_plot(df_vis_curve), width = 6, height = 4)
ggsave(file.path(PAPER_DIR, "greedy_aic_curve_val.pdf"),
       make_curve_plot(df_val_curve), width = 6, height = 4)
cat("Saved: paper/greedy_aic_curve_vis.pdf, paper/greedy_aic_curve_val.pdf\n")

for (split in c("vis", "val")) {
  if (split == "vis") {
    df_split <- bind_rows(
      make_long_heatmap_dw(hi_idx_list, lbl_hi),
      make_long_heatmap_dw(lo_idx_list, lbl_lo))
    outfile <- file.path(PAPER_DIR, "greedy_vis_heatmap.pdf")
  } else {
    df_split <- bind_rows(
      make_long_heatmap_dw(af_idx_list, lbl_af),
      make_long_heatmap_dw(co_idx_list, lbl_co))
    outfile <- file.path(PAPER_DIR, "greedy_valence_heatmap.pdf")
  }
  df_split <- mutate(df_split, term = factor(term, levels = term_ord))
  p <- ggplot(df_split, aes(x = step_num, y = term, fill = freq)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 2.6, colour = "white", fontface = "bold") +
    facet_wrap(~group, ncol = 2) +
    scale_fill_gradient(low = "white", high = "#2171b5", limits = c(0, 1),
                        name = "Prop.\ndomains") +
    scale_x_continuous(breaks = seq_len(n_steps), expand = c(0, 0)) +
    labs(x = "Greedy step", y = NULL) +
    heatmap_theme
  ggsave(outfile, p, width = 10, height = 4.5)
}
cat("Saved: paper/greedy_vis_heatmap.pdf, paper/greedy_valence_heatmap.pdf\n")

# -- freq_in_top4 dot plot -----------------------------------

lbl_ov <- sprintf("Overall (%d domains)", n_dom)
af_top4 <- dw_top4(af_idx_list)
co_top4 <- dw_top4(co_idx_list)

top4_long <- data.frame(
  term    = rep(term_labels, 4),
  freq    = c(dw_greedy_top4, hi_top4, lo_top4, dw_greedy_top4),  # overall reused
  group   = rep(c(lbl_ov, lbl_hi, lbl_lo, lbl_ov), each = n_terms)
) |>
  filter(group != lbl_ov | !duplicated(term)) |>
  mutate(term = factor(term,
    levels = names(sort(dw_greedy_top4[term_labels])))) |>
  mutate(group = factor(group, levels = c(lbl_ov, lbl_hi, lbl_lo)))

# Rebuild cleanly
top4_df <- data.frame(
  term    = term_labels,
  overall = as.numeric(dw_greedy_top4[term_labels]),
  high    = as.numeric(hi_top4[term_labels]),
  low     = as.numeric(lo_top4[term_labels])
) |>
  arrange(desc(overall)) |>
  mutate(term = factor(term, levels = rev(term)))

top4_long2 <- top4_df |>
  pivot_longer(c(overall, high, low), names_to = "group", values_to = "freq") |>
  mutate(group = recode(group,
    overall = lbl_ov,
    high    = lbl_hi,
    low     = lbl_lo))

p_top4 <- ggplot(top4_long2, aes(x = freq, y = term, colour = group, shape = group)) +
  geom_point(size = 3.2) +
  scale_colour_manual(
    values = setNames(c("black", "#2171b5", "#cb181d"), c(lbl_ov, lbl_hi, lbl_lo)),
    name = NULL) +
  scale_shape_manual(
    values = setNames(c(18, 16, 17), c(lbl_ov, lbl_hi, lbl_lo)),
    name = NULL) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "gray60") +
  scale_x_continuous(labels = scales::percent_format(),
                     limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(x = "Proportion of domains selected in top 4", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(PAPER_DIR, "freq_in_top4.pdf"), p_top4, width = 7, height = 4.5)
cat("Saved: paper/freq_in_top4.pdf\n")

# -- Two-panel sensitivity figure ----------------------------

term_order <- rev(comparison$Term)

freq_df <- comparison |>
  mutate(highlight = ifelse(Term %in% baseline, "Baseline", "Other"),
         Term = factor(Term, levels = term_order)) |>
  pivot_longer(cols = c(Greedy_dw, Univ_dw),
               names_to = "criterion", values_to = "value") |>
  mutate(criterion = recode(criterion,
    Greedy_dw = "Greedy top-4 (domain-wtd)",
    Univ_dw   = "Univariate top-4 (domain-wtd)"))

p_freq <- ggplot(freq_df, aes(x = value, y = Term,
                               colour = highlight, shape = criterion)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "gray70") +
  geom_point(size = 2.8, alpha = 0.85,
             position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(Baseline = "#2171b5", Other = "gray55"),
                      name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent,
                     name = "Proportion of domains") +
  labs(y = NULL, title = "(a) Selection frequency") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        legend.box = "vertical")

too_df <- comparison |>
  mutate(highlight = ifelse(Term %in% baseline, "Baseline", "Other"),
         Term = factor(Term, levels = term_order)) |>
  pivot_longer(cols = c(TOO4_loss, TOO12_loss),
               names_to = "model", values_to = "loss") |>
  filter(!is.na(loss)) |>
  mutate(model = recode(model,
    TOO4_loss  = "Drop from 4-term model",
    TOO12_loss = "Drop from 12-term model"))

p_too <- ggplot(too_df, aes(x = loss, y = Term,
                              colour = highlight, shape = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray70") +
  geom_point(size = 2.8, alpha = 0.85,
             position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(Baseline = "#2171b5", Other = "gray55"),
                      name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  scale_x_continuous(name = "Mean AIC increase when term dropped") +
  labs(y = NULL, title = "(b) Take-one-out AIC loss") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        legend.box = "vertical")

p_combined <- gridExtra::arrangeGrob(p_freq, p_too, ncol = 2)
ggsave(file.path(PAPER_DIR, "aic_sensitivity_plot.pdf"),
       p_combined, width = 11, height = 5.5)
cat("Saved: paper/aic_sensitivity_plot.pdf\n")

cat("\nDone.\n")
cat(sprintf("  %d networks / %d domains\n", n_nets, n_dom))
cat(sprintf("  n_nonevents = %d\n", N_NONEVENTS))
cat(sprintf("  term label == column name throughout (no translation)\n"))
