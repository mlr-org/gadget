#!/usr/bin/env Rscript
# Summarize GADGET vs effector benchmark results and generate figures.
# Run: Rscript simulation/summarize_benchmark.R [--indir DIR] [--figdir DIR]
#                                                [--fixed-N N] [--fixed-D D]

args <- commandArgs(trailingOnly = TRUE)
indir   <- "simulation/results/benchmark"
figdir  <- "simulation/results/figures"
fixed_N <- 1000L
fixed_D <- 10L
i <- 1
while (i <= length(args)) {
  if      (args[i] == "--indir"   && i < length(args)) { indir   <- args[i+1]; i <- i+2 }
  else if (args[i] == "--figdir"  && i < length(args)) { figdir  <- args[i+1]; i <- i+2 }
  else if (args[i] == "--fixed-N" && i < length(args)) { fixed_N <- as.integer(args[i+1]); i <- i+2 }
  else if (args[i] == "--fixed-D" && i < length(args)) { fixed_D <- as.integer(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(ggplot2)
library(scales)
if (!requireNamespace("gridExtra", quietly = TRUE)) stop("Install gridExtra")
library(gridExtra)

# ---------------------------------------------------------------------------
# Load CSVs
# ---------------------------------------------------------------------------

load_csv <- function(filename, model_type) {
  f <- file.path(indir, filename)
  if (!file.exists(f)) return(data.table())
  dt <- fread(f)
  dt[, model_type := model_type]
  dt
}

dt <- rbindlist(list(
  load_csv("benchmark_gadget_rf.csv",   "rf"),
  load_csv("benchmark_gadget_toy.csv",  "toy"),
  load_csv("benchmark_effector_rf.csv", "rf"),
  load_csv("benchmark_effector_toy.csv","toy")
), use.names = TRUE, fill = TRUE)

if (nrow(dt) == 0) stop("No benchmark data found. Run benchmark scripts first.")
if ("status" %in% names(dt)) {
  n_err <- sum(dt[["status"]] == "error", na.rm = TRUE)
  if (is.finite(n_err) && n_err > 0L) {
    message(
      "Warning: ", n_err,
      " gadget benchmark row(s) have status=error (NA times); see error_message in raw CSV."
    )
  }
}
if (!("impl" %in% names(dt))) dt[, impl := ""]
dt[is.na(impl) | impl == "", impl := fifelse(package == "gadget", "default", package)]
dt[, n_grid      := as.integer(fifelse(is.na(n_grid)      | n_grid      == "", NA_character_, as.character(n_grid)))]
dt[, n_intervals := as.integer(fifelse(is.na(n_intervals) | n_intervals == "", NA_character_, as.character(n_intervals)))]

# ---------------------------------------------------------------------------
# Predict baseline subtraction (RF global methods only)
# ---------------------------------------------------------------------------

load_baseline <- function(filename) {
  f <- file.path(indir, filename)
  if (!file.exists(f)) return(data.table(package=character(), N=integer(), D=integer(), predict_time_mean=numeric()))
  fread(f)[, .(package, N, D, predict_time_mean)]
}
base <- rbind(load_baseline("predict_baseline_gadget.csv"),
              load_baseline("predict_baseline_effector.csv"))
dt <- merge(dt, base, by = c("package","N","D"), all.x = TRUE)
dt[, predict_time_mean := fifelse(is.na(predict_time_mean), 0, predict_time_mean)]
dt[, method_time := fifelse(
  model_type == "rf" & method %in% c("global_pdp","global_ale"),
  time_sec - predict_time_mean, time_sec
)]

# ---------------------------------------------------------------------------
# Aggregate
# ---------------------------------------------------------------------------

# Aggregate within each sub-experiment first to avoid pooling (N=fixed_N, D=fixed_D, res=default)
# being counted multiple times (it would otherwise appear in all three sub-experiments).
if ("sub_experiment" %in% names(dt)) {
  agg_sub <- dt[, .(
    time_mean        = mean(time_sec,    na.rm = TRUE),
    time_sd          = sd(time_sec,      na.rm = TRUE),
    method_time_mean = mean(method_time, na.rm = TRUE),
    method_time_sd   = sd(method_time,  na.rm = TRUE),
    n_rep = .N
  ), by = .(package, impl, method, sub_experiment, model_type, N, D, n_grid, n_intervals)]
  # Then take the mean of sub-experiment means for the same (N, D, res) cell
  agg <- agg_sub[, .(
    time_mean        = mean(time_mean,        na.rm = TRUE),
    time_sd          = mean(time_sd,          na.rm = TRUE),
    method_time_mean = mean(method_time_mean, na.rm = TRUE),
    method_time_sd   = mean(method_time_sd,   na.rm = TRUE),
    n_rep            = sum(n_rep)
  ), by = .(package, impl, method, model_type, N, D, n_grid, n_intervals)]
} else {
  agg <- dt[, .(
    time_mean        = mean(time_sec,    na.rm = TRUE),
    time_sd          = sd(time_sec,      na.rm = TRUE),
    method_time_mean = mean(method_time, na.rm = TRUE),
    method_time_sd   = sd(method_time,  na.rm = TRUE),
    n_rep = .N
  ), by = .(package, impl, method, model_type, N, D, n_grid, n_intervals)]
}

write.csv(agg, file.path(indir, "summary.csv"), row.names = FALSE)
message("Written: ", file.path(indir, "summary.csv"))

# Derived columns
agg[, method_base := sub("_(split|total)$", "", method)]

agg[, pkg := fifelse(
  package == "gadget" & impl == "effect-cpp", "effect-cpp",
  fifelse(package == "gadget" & impl == "effect-r", "effect-r",
  fifelse(package == "effector", "effector", paste0(package, "-", impl)))
)]
agg[package == "gadget"   & grepl("_split$", method), pkg := "gadget-split"]
agg[package == "gadget"   & grepl("_total$", method), pkg := "gadget-total"]
agg[package == "effector" & grepl("_split$", method), pkg := "effector-split"]
agg[package == "effector" & grepl("_total$", method), pkg := "effector-total"]

# ---------------------------------------------------------------------------
# Fixed colour palette (consistent across all panels)
# ---------------------------------------------------------------------------

PKG_COLOURS <- c(
  "effect-cpp"    = "#1f77b4",
  "effect-r"      = "#aec7e8",
  "effector"      = "#2ca02c",
  "gadget-split"  = "#1f77b4",
  "gadget-total"  = "#9467bd",
  "effector-split"= "#2ca02c",
  "effector-total"= "#d62728"
)

# ---------------------------------------------------------------------------
# Panel builder
# ---------------------------------------------------------------------------

mkplot <- function(sub, xvar, xlab, yvar = "method_time_mean", all_pkgs) {
  sub <- sub[!is.na(sub[[xvar]]) & sub[[xvar]] != ""]
  if (nrow(sub) == 0 || length(unique(sub[[xvar]])) < 2) return(NULL)
  xv <- as.numeric(sub[[xvar]])
  sub[[xvar]] <- xv
  if (nrow(sub[!is.na(sub[[yvar]])]) == 0) return(NULL)
  sub  <- as.data.frame(sub)
  sub  <- sub[order(sub$pkg, sub[[xvar]]), , drop = FALSE]
  sdcol <- if (yvar == "method_time_mean") "method_time_sd" else "time_sd"
  sub[[sdcol]][is.na(sub[[sdcol]])] <- 0

  use_log_x <- identical(xvar, "N") && all(xv > 0, na.rm = TRUE)
  if (use_log_x) {
    sub$.xplot <- log10(sub[[xvar]])
    x_aes <- aes(x = .xplot, y = !!rlang::sym(yvar),
                 color = pkg, fill = pkg, group = pkg)
  } else {
    x_aes <- aes(x = !!rlang::sym(xvar), y = !!rlang::sym(yvar),
                 color = pkg, fill = pkg, group = pkg)
  }

  pal_use  <- PKG_COLOURS[all_pkgs]
  pal_fill <- alpha(pal_use, 0.25)

  p <- ggplot(sub, x_aes) +
    geom_ribbon(aes(ymin = pmax(.data[[yvar]] - .data[[sdcol]], 0),
                    ymax = .data[[yvar]] + .data[[sdcol]]),
                alpha = 0.22, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_color_manual(values = pal_use,  drop = FALSE, breaks = all_pkgs) +
    scale_fill_manual( values = pal_fill, drop = FALSE, guide  = "none") +
    theme_bw(base_size = 10) +
    theme(
      legend.position   = "bottom",
      legend.title      = element_blank(),
      legend.key.size   = unit(0.5, "lines"),
      legend.text       = element_text(size = 8),
      plot.title        = element_text(size = 9, face = "bold", hjust = 0.5),
      panel.grid.minor  = element_blank()
    ) +
    labs(y = "Runtime (s)", x = xlab)

  if (use_log_x) {
    n_brks   <- sort(unique(xv))
    log_brks <- log10(n_brks)
    p <- p + scale_x_continuous(
      breaks = log_brks,
      labels = function(b) comma(round(10^b)),
      expand = expansion(mult = c(0.12, 0.14))
    )
  } else if (is.numeric(xv) && max(xv, na.rm = TRUE) > 1000) {
    p <- p + scale_x_continuous(labels = comma)
  }
  p
}

# ---------------------------------------------------------------------------
# Shared-legend extraction helper
# ---------------------------------------------------------------------------

extract_legend <- function(p) {
  g     <- ggplotGrob(p)
  legs  <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(legs) == 0) return(NULL)
  g$grobs[[legs[1]]]
}

# ---------------------------------------------------------------------------
# Build figures: one per scope (global / regional)
# ---------------------------------------------------------------------------

default_res <- 20L

fig_defs <- list(
  # Figure 1: global feature-effect computation (effect-cpp / effect-r / effector)
  list(
    scope = "global", title = "Global Methods",
    pkg_filter = NULL,  # all pkgs
    methods = list(
      list(method = "global_pdp", label = "PDP", res_var = "n_grid",      res_label = "Grid size"),
      list(method = "global_ale", label = "ALE", res_var = "n_intervals",  res_label = "Intervals")
    )
  ),
  # Figure 2: regional split only (gadget-split vs effector-split)
  # Shows the pure tree-fitting cost on both sides; total excluded here.
  list(
    scope = "regional", title = "Regional Methods - Split (tree fitting only)",
    pkg_filter = c("gadget-split", "effector-split"),
    methods = list(
      list(method = "regional_pdp", label = "PDP", res_var = "n_grid",     res_label = "Grid size"),
      list(method = "regional_ale", label = "ALE", res_var = "n_intervals", res_label = "Intervals")
    )
  ),
  # Figure 3: full pipeline total (gadget-total vs effector-total)
  list(
    scope = "total", title = "Regional Methods - Total (global + split)",
    pkg_filter = c("gadget-total", "effector-total"),
    methods = list(
      list(method = "regional_pdp", label = "PDP", res_var = "n_grid",     res_label = "Grid size"),
      list(method = "regional_ale", label = "ALE", res_var = "n_intervals", res_label = "Intervals")
    )
  )
)

model_types  <- c("rf",       "toy")
model_labels <- c(rf = "RF model", toy = "Toy model")

for (fd in fig_defs) {
  # Collect pkg labels for this figure (optionally filtered)
  scope_data <- agg[method_base %in% sapply(fd$methods, `[[`, "method")]
  if (!is.null(fd$pkg_filter)) scope_data <- scope_data[pkg %in% fd$pkg_filter]
  all_pkgs   <- if (!is.null(fd$pkg_filter)) fd$pkg_filter else sort(unique(scope_data$pkg))

  panels <- list()   # (nrow * 3) list of ggplot / grob

  for (md in fd$methods) {
    for (mt in model_types) {
      sub  <- agg[method_base == md$method & model_type == mt]
      if (!is.null(fd$pkg_filter)) sub <- sub[pkg %in% fd$pkg_filter]
      if (nrow(sub) == 0) next

      is_pdp  <- grepl("pdp", md$method)
      res_col <- if (is_pdp) "n_grid" else "n_intervals"

      sub_n <- sub[D == fixed_D & !is.na(get(res_col)) & get(res_col) == default_res]
      sub_d <- sub[N == fixed_N & !is.na(get(res_col)) & get(res_col) == default_res]
      sub_r <- sub[N == fixed_N & D == fixed_D & !is.na(get(res_col))]

      row_label <- sprintf("%s / %s", md$label, model_labels[[mt]])

      p_n <- mkplot(sub_n, "N",      "Sample size N", all_pkgs = all_pkgs)
      p_d <- mkplot(sub_d, "D",      "Features D",    all_pkgs = all_pkgs)
      p_r <- mkplot(sub_r, res_col,  md$res_label,    all_pkgs = all_pkgs)

      add_title <- function(p, suffix) {
        if (!is.null(p)) p + ggtitle(paste0(row_label, " - ", suffix)) else NULL
      }
      p_n <- add_title(p_n, "vs N")
      p_d <- add_title(p_d, "vs D")
      p_r <- add_title(p_r, paste0("vs ", md$res_label))

      blank <- ggplot() + theme_void()
      panels <- c(panels, list(
        if (!is.null(p_n)) p_n else blank,
        if (!is.null(p_d)) p_d else blank,
        if (!is.null(p_r)) p_r else blank
      ))
    }
  }

  nrow_panels <- length(panels) / 3
  if (nrow_panels == 0) next

  # Extract shared legend from the first non-null panel
  first_p  <- Filter(Negate(is.null), panels)[[1]]
  leg_grob <- extract_legend(first_p)

  # Strip legends from all panels (they share one at the bottom)
  panels_no_leg <- lapply(panels, function(p) {
    p + theme(legend.position = "none")
  })

  # Assemble: title | grid of panels | shared legend
  title_grob <- grid::textGrob(
    fd$title,
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )

  panel_grid <- gridExtra::arrangeGrob(
    grobs = panels_no_leg,
    ncol  = 3,
    padding = unit(0.3, "lines")
  )

  panel_h_cm <- 3.2 * nrow_panels * 2.54  # approx in cm (will be normalised)
  leg_h_cm   <- 1.2
  ttl_h_cm   <- 0.6

  if (!is.null(leg_grob)) {
    combined <- gridExtra::arrangeGrob(
      title_grob, panel_grid, leg_grob,
      ncol    = 1,
      heights = c(ttl_h_cm, panel_h_cm, leg_h_cm)
    )
  } else {
    combined <- gridExtra::arrangeGrob(
      title_grob, panel_grid,
      ncol    = 1,
      heights = c(ttl_h_cm, panel_h_cm)
    )
  }

  fname <- sprintf("%s_methods.png", fd$scope)
  ggsave(
    file.path(figdir, fname), combined,
    width = 15, height = 3.2 * nrow_panels + 1.2, dpi = 150
  )
  message("Written: ", file.path(figdir, fname))
}

message("Done.")
