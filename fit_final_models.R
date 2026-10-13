rm(list = ls())

# ============================================================
# fit_final_models.R
#
# For each of the 217 networks in greedy_results7.Rdata, refit
# the "final" model (last greedy step where prop_delta_aic < 0),
# extract per-term coefficients and SEs, and aggregate into a
# summary table by visibility group for the paper appendix.
#
# fb* networks are skipped (raw files not available locally).
# Outputs: final_model_coefs.Rdata + appendix_coef_table.tex
# ============================================================

library(llrem)
library(survival)
library(parallel)

PROJ_DIR <- "/path/to/FirstFour"       # edit to your local project root

# ============================================================
# term_groups: label -> column name in pe (same as greedy_opt2.R)
# ============================================================
term_groups <- list(
  Imitation      = "Imitation",
  DyadRecency    = "DyadRecencyLog",
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

# ============================================================
# Map stored Windows basename -> local raw file path
# ============================================================
resolve_local_path <- function(bn) {
  stem <- sub("[.]Rdatadf[.]Rdata$", "", bn, ignore.case = TRUE)
  stem <- sub("[.]rdatadf[.]Rdata$", "", stem, ignore.case = TRUE)

  if (grepl("^eu[0-9]", stem))
    return(sprintf("~/Downloads/DataP/Email/%s.Rdata", stem))
  if (stem == "TSOF")
    return("~/Downloads/DataP/Email/TSOF.Rdata")
  if (grepl("^gangnet", stem))
    return(sprintf("~/Downloads/DataP/Gangs/%s.Rdata", stem))
  if (grepl("^network", stem))
    return(sprintf("~/Downloads/DataP/ftf/%s.Rdata", stem))
  if (grepl("^ucimsg", stem))
    return(sprintf("~/Downloads/DataP/PMs/%s.Rdata", stem))
  if (grepl("^fb", stem))
    return(NA_character_)  # not available locally
  if (grepl("^sms", stem))
    return(sprintf("~/Downloads/DataP/SMS/%s.Rdata", stem))
  if (stem == "twitter_boston_relevent")
    return("~/Downloads/DataP/Twitter/twitter_boston_relevent.rdata")
  if (stem == "conflict")
    return("~/Downloads/DataP/War/conflict.Rdata")
  if (stem == "crimephone")
    return("~/Projects/Endogeneity/DataP/Phone/crimephone.Rdata")
  if (stem == "Phone")
    return("~/Downloads/DataP/Phone/Phone.Rdata")
  if (stem == "McFarland")
    return("~/Downloads/DataP/Classroom/McFarland.Rdata")
  return(NA_character_)
}

# ============================================================
# Load raw edgelist from any supported format;
# return data.frame with cols: time, sender, receiver
# ============================================================
load_edgelist <- function(path) {
  e <- new.env()
  load(path_expand(path), envir = e)
  el <- e$edgelist
  if (is.matrix(el)) el <- as.data.frame(el)
  # Standardise: first three columns = time, sender, receiver
  el <- el[, 1:3]
  colnames(el) <- c("time", "sender", "receiver")
  el <- el[order(el$time), ]
  el
}

path_expand <- function(p) path.expand(p)

# ============================================================
# Determine the "final step" for one network's path:
# last step where prop_delta_aic < 0
# ============================================================
get_final_step <- function(path) {
  pdas <- sapply(path, function(s) s$prop_delta_aic)
  good <- which(!is.na(pdas) & pdas < 0)
  if (length(good) == 0L) 0L else max(good)
}

# ============================================================
# Fit the final model for one network element from greedy_results
# Returns: list(df, cat, visible, final_step, n_terms, coef, se)
#          or NULL on any failure
# ============================================================
fit_one <- function(elem, term_groups) {
  bn        <- basename(elem$df)
  local_path <- resolve_local_path(bn)

  if (is.na(local_path) || !file.exists(path.expand(local_path)))
    return(NULL)

  final_step <- get_final_step(elem$path)
  if (final_step == 0L) return(NULL)

  # Extract the term labels added in steps 1..final_step
  term_labels <- sapply(elem$path[seq_len(final_step)], function(s) s$added)

  # Map labels to column names in pe
  term_cols <- unlist(term_groups[term_labels], use.names = FALSE)
  if (any(is.na(term_cols))) return(NULL)

  # Load and prepare edgelist
  raw_el <- tryCatch(load_edgelist(local_path), error = function(e) NULL)
  if (is.null(raw_el)) return(NULL)

  prep_el <- tryCatch(prepare_edgelist(raw_el), error = function(e) NULL)
  if (is.null(prep_el)) return(NULL)

  all_nodes    <- unique(c(prep_el$sender, prep_el$receiver))
  base_riskset <- llrem:::get_base_riskset(NULL, all_nodes)
  pe <- tryCatch(
    llrem:::remstats_generate(prep_el,
                              n            = length(all_nodes),
                              n_nonevents  = 200L,
                              base_riskset = base_riskset),
    error = function(e) NULL
  )
  if (is.null(pe)) return(NULL)
  pe <- llrem:::add_survival_cols(pe)

  # Build formula
  f <- as.formula(paste0(
    "Surv(time1, occurred) ~ ",
    paste(term_cols, collapse = " + "),
    " + strata(time2)"
  ))

  fit <- tryCatch(coxph(f, data = pe), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  cs  <- coef(summary(fit))
  # Map col names back to labels
  col2label <- setNames(names(term_groups), unlist(term_groups))
  row_labels <- col2label[rownames(cs)]

  list(
    df         = elem$df,
    cat        = elem$cat,
    visible    = elem$visible,
    final_step = final_step,
    coef       = setNames(cs[, "coef"],       row_labels),
    se         = setNames(cs[, "se(coef)"],   row_labels)
  )
}

# ============================================================
# Run fits in parallel
# ============================================================
load(file.path(PROJ_DIR, "greedy_results7.Rdata"))

cats    <- sapply(greedy_results, function(x) x$cat)
vis     <- sapply(greedy_results, function(x) x$visible)
vis[cats == "Phone (DoD)"] <- "Y"
vis_bin <- ifelse(vis == "Y", "High", "Low")

n_cores <- min(detectCores() - 1L, 8L)
cl <- makeCluster(n_cores)
clusterExport(cl, c("term_groups", "resolve_local_path", "load_edgelist",
                    "path_expand", "get_final_step", "fit_one"),
              envir = environment())
clusterCall(cl, function() {
  library(llrem)
  library(survival)
})

cat("Fitting final models for", length(greedy_results), "networks...\n")
fits <- tryCatch(
  parLapply(cl, greedy_results, function(elem) fit_one(elem, term_groups)),
  finally = stopCluster(cl)
)

# Remove NULLs
n_total   <- length(fits)
fits      <- fits[!sapply(fits, is.null)]
vis_fits  <- vis_bin[!sapply(fits[seq_len(n_total)], is.null)]
cat(sprintf("%d / %d networks successfully fit.\n", length(fits), n_total))

# Rebuild vis_bin for the retained fits
fit_vis <- sapply(fits, function(x) {
  c_val <- x$cat
  v_val <- x$visible
  if (c_val == "Phone (DoD)") v_val <- "Y"
  ifelse(v_val == "Y", "High", "Low")
})

save(fits, fit_vis, file = file.path(PROJ_DIR, "final_model_coefs.Rdata"))
cat("Saved to final_model_coefs.Rdata\n")

# ============================================================
# Aggregate: mean, SD, N by term × vis group
# ============================================================

all_terms <- c("SndInitiation", "IOviolation", "Repetition", "RecRecencyRank",
               "RecRecencyLog", "RRecInitiation", "WeightedRec",
               "Imitation", "OutDegSnd", "DyadRecency", "NOTPSnd", "NITPSnd")

# Collect coef vectors
all_coefs <- do.call(rbind, lapply(seq_along(fits), function(i) {
  x <- fits[[i]]
  v <- fit_vis[i]
  row <- c(vis = v, setNames(rep(NA_real_, length(all_terms)), all_terms))
  for (nm in names(x$coef)) {
    if (nm %in% all_terms) row[nm] <- x$coef[nm]
  }
  row
}))
all_coefs <- as.data.frame(all_coefs, stringsAsFactors = FALSE)
for (t in all_terms) all_coefs[[t]] <- as.numeric(all_coefs[[t]])

summarise_term <- function(x) {
  x <- x[!is.na(x)]
  n    <- length(x)
  mn   <- mean(x)
  sd_  <- sd(x)
  tstat <- if (n > 1) mn / (sd_ / sqrt(n)) else NA
  stars <- if (!is.na(tstat)) {
    if (abs(tstat) >= 3.29) "***"
    else if (abs(tstat) >= 2.58) "**"
    else if (abs(tstat) >= 1.96) "*"
    else ""
  } else ""
  list(n = n, mean = mn, sd = sd_, tstat = tstat, stars = stars,
       cell = sprintf("%.3f (%.3f)%s", mn, sd_, stars))
}

groups <- list(
  hi  = all_coefs$vis == "High",
  lo  = all_coefs$vis == "Low",
  all = rep(TRUE, nrow(all_coefs))
)

rows <- lapply(all_terms, function(t) {
  hi_s  <- summarise_term(all_coefs[[t]][groups$hi])
  lo_s  <- summarise_term(all_coefs[[t]][groups$lo])
  all_s <- summarise_term(all_coefs[[t]][groups$all])
  data.frame(
    Term     = t,
    N_hi     = hi_s$n,
    N_lo     = lo_s$n,
    N_all    = all_s$n,
    Hi       = hi_s$cell,
    Lo       = lo_s$cell,
    Overall  = all_s$cell,
    stringsAsFactors = FALSE
  )
})
tab <- do.call(rbind, rows)
tab <- tab[tab$N_all > 0, ]  # drop terms not appearing in any final model

cat("\n=== Summary coefficient table (final models) ===\n\n")
print(tab, row.names = FALSE)

# ============================================================
# Write LaTeX table
# ============================================================

# Term group structure matching the paper
term_group_labels <- list(
  "Forward-dyad history"            = c("SndInitiation", "Repetition"),
  "Reverse-dyad history"            = c("RecRecencyRank", "RecRecencyLog",
                                        "RRecInitiation", "WeightedRec"),
  "Popularity, sender activity, triadic" = c("Imitation", "OutDegSnd",
                                              "DyadRecency", "NOTPSnd", "NITPSnd"),
  "Sequential ordering"             = c("IOviolation")
)

n_hi_total  <- sum(fit_vis == "High")
n_lo_total  <- sum(fit_vis == "Low")
n_all_total <- length(fit_vis)

tex_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  sprintf("\\caption{Mean (SD) coefficients in network-specific final greedy models, by visibility group. $N$ gives the number of networks whose final model included each term. High-visibility: $N=%d$; low-visibility: $N=%d$; total: $N=%d$ (7 Facebook PM networks excluded: raw data unavailable). Stars denote one-sample $t$-test vs.\\ 0: $^*p<.05$, $^{**}p<.01$, $^{***}p<.001$.}",
            n_hi_total, n_lo_total, n_all_total),
  "\\label{tab:final_coefs}",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  sprintf("Term & $N$ (High) & High-vis & $N$ (Low) & Low-vis & Overall ($N=%d$) \\\\",
          n_all_total),
  "\\midrule"
)

for (grp_name in names(term_group_labels)) {
  tex_lines <- c(tex_lines,
    sprintf("\\multicolumn{6}{l}{\\textit{%s}} \\\\", grp_name))
  for (t in term_group_labels[[grp_name]]) {
    row_idx <- which(tab$Term == t)
    if (length(row_idx) == 0) next
    r <- tab[row_idx, ]
    hi_cell  <- if (r$N_hi  > 0) sprintf("%s", r$Hi)      else "---"
    lo_cell  <- if (r$N_lo  > 0) sprintf("%s", r$Lo)      else "---"
    all_cell <- if (r$N_all > 0) sprintf("%s", r$Overall)  else "---"
    tex_lines <- c(tex_lines,
      sprintf("\\quad %s & %d & %s & %d & %s & %s \\\\",
              t, r$N_hi, hi_cell, r$N_lo, lo_cell, all_cell))
  }
}

tex_lines <- c(tex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

tex_path <- file.path(PROJ_DIR, "paper", "appendix_coef_table.tex")
writeLines(tex_lines, tex_path)
cat("\nLaTeX table written to:", tex_path, "\n")
