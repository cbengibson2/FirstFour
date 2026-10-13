load("final_model_coefs.Rdata")

cats <- sapply(fits, function(x) x$cat)
vis  <- sapply(fits, function(x) x$visible)
vis[cats == "Phone (DoD)"] <- "Y"
fit_vis <- ifelse(vis == "Y", "High", "Low")

all_term_labels <- c(
  "SndInitiation", "Repetition",
  "RecRecencyRank", "RecRecencyLog", "RRecInitiation", "WeightedRec",
  "Imitation", "OutDegSnd", "DyadRecency", "NOTPSnd", "NITPSnd",
  "IOviolation"
)

build_mat <- function(slot) {
  do.call(rbind, lapply(fits, function(x) {
    row <- setNames(rep(NA_real_, length(all_term_labels)), all_term_labels)
    for (nm in names(x[[slot]]))
      if (nm %in% all_term_labels) row[nm] <- x[[slot]][nm]
    row
  }))
}
coef_mat <- build_mat("coef")
se_mat   <- build_mat("se")

summarise_term_se <- function(coefs, ses) {
  ok  <- !is.na(coefs) & !is.na(ses)
  n   <- sum(ok)
  if (n == 0) return(list(n = 0, cell = "---"))
  mn  <- mean(coefs[ok])
  mse <- mean(ses[ok])
  t   <- mn / mse
  stars <- if (abs(t) >= 3.29) "***" else if (abs(t) >= 2.58) "**" else
            if (abs(t) >= 1.96) "*" else ""
  list(n = n, cell = sprintf("%.3f (%.3f)%s", mn, mse, stars))
}

tab <- do.call(rbind, lapply(all_term_labels, function(t) {
  hi <- summarise_term_se(coef_mat[fit_vis == "High", t], se_mat[fit_vis == "High", t])
  lo <- summarise_term_se(coef_mat[fit_vis == "Low",  t], se_mat[fit_vis == "Low",  t])
  ov <- summarise_term_se(coef_mat[, t], se_mat[, t])
  data.frame(Term = t, N_hi = hi$n, Hi = hi$cell,
             N_lo = lo$n, Lo = lo$cell,
             N_all = ov$n, Overall = ov$cell,
             stringsAsFactors = FALSE)
}))
tab <- tab[tab$N_all > 0, ]

n_hi  <- sum(fit_vis == "High")
n_lo  <- sum(fit_vis == "Low")
n_all <- length(fit_vis)

term_groups_table <- list(
  "Forward-dyad history"                 = c("SndInitiation", "Repetition"),
  "Reverse-dyad history"                 = c("RecRecencyRank", "RecRecencyLog",
                                              "RRecInitiation", "WeightedRec"),
  "Popularity, sender activity, triadic" = c("Imitation", "OutDegSnd",
                                              "DyadRecency", "NOTPSnd", "NITPSnd"),
  "Sequential ordering"                  = "IOviolation"
)

caption <- paste0(
  "Mean coefficient and mean model SE (in parentheses) across",
  " network-specific final greedy models, by visibility group.",
  " $N$ is the number of networks whose final model included each term.",
  " Stars indicate mean~coef / mean~SE $\\geq 1.96$ ($^*$), $2.58$ ($^{**}$),",
  " or $3.29$ ($^{***}$), reflecting average within-model precision.",
  " Terms with $N \\leq 3$ appear in very few networks;",
  " their entries are included for completeness but should not be",
  " treated as evidence of cross-network generalizability."
)

lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  paste0("\\caption{", caption, "}"),
  "\\label{tab:final_coefs}",
  "\\begin{tabular}{lcrcrr}",
  "\\toprule",
  sprintf("& \\multicolumn{2}{c}{High-vis ($n=%d$)} & \\multicolumn{2}{c}{Low-vis ($n=%d$)} & Overall ($n=%d$) \\\\",
          n_hi, n_lo, n_all),
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Term & Mean (SE) & $N$ & Mean (SE) & $N$ & Mean (SE) \\\\",
  "\\midrule"
)

for (grp_name in names(term_groups_table)) {
  lines <- c(lines,
    sprintf("\\multicolumn{6}{l}{\\textit{%s}} \\\\", grp_name))
  for (t in term_groups_table[[grp_name]]) {
    r <- tab[tab$Term == t, ]
    if (nrow(r) == 0) next
    lines <- c(lines,
      sprintf("\\quad %s & %s & %d & %s & %d & %s \\\\",
              t, r$Hi, r$N_hi, r$Lo, r$N_lo, r$Overall))
  }
}

lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")

writeLines(lines, "paper/appendix_coef_table.tex")
cat("Written to paper/appendix_coef_table.tex\n")
