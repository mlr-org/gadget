#' Build per-feature ALE panels (mean curve only)
#'
#' Helper used by \code{plot_tree_ale()} to generate per-feature ALE mean
#' panels (optionally with overlaid observation points).
#'
#' @param curves (`list()`) \cr
#'   Output of \code{prepare_plot_data_ale} for a node.
#' @param color_ale (`character(1)`) \cr
#'   Color for ALE curves.
#' @param target_feature_name (`character(1)`) \cr
#'   Target column name; used as the y-axis label (same convention as PD plots).
#' @param mean_center (`logical(1)`) \cr
#'   Whether ALE curves are mean-centered; controls legend text (\code{"Mean centered ALE"} vs
#'   \code{"ALE"}), matching PD ICE/PDP labeling.
#' @param ymin,ymax (`numeric(1)` or `NULL`) \cr
#'   Y-axis limits.
#' @param show_point (`logical(1)`) \cr
#'   Whether to add observation points.
#' @param point_values (`list()` or `NULL`) \cr
#'   Per-feature data.frames with \code{x}, \code{y}; used when \code{show_point = TRUE}.
#' @param x_limits (`list()` or `NULL`) \cr
#'   Per-feature x-axis: numeric \code{c(xmin, xmax)} or character (level order).
#'
#' @return (`list()`) \cr
#'   Named list of ggplot objects per feature.
#' @keywords internal
plot_regional_ale = function(curves, color_ale = "lightcoral", target_feature_name,
  mean_center = TRUE, ymin = NULL, ymax = NULL, show_point = FALSE,
  point_values = NULL, x_limits = NULL) {
  checkmate::assert_character(target_feature_name, len = 1, .var.name = "target_feature_name")
  checkmate::assert_flag(mean_center, .var.name = "mean_center")
  if (!length(curves)) {
    return(list())
  }

  ale_label = if (mean_center) "Mean centered ALE" else "ALE"
  color_vals = color_ale
  names(color_vals) = ale_label

  feats = names(curves)
  out = lapply(feats, function(feat) {
    mean_dt = curves[[feat]]$mean_effect
    is_numeric = "x_grid" %in% names(mean_dt) && is.numeric(mean_dt$x_grid)
    xlim_feat = if (!is.null(x_limits)) x_limits[[feat]] else NULL

    has_points = show_point && !is.null(point_values) && !is.null(point_values[[feat]])
    obs_df = if (has_points) point_values[[feat]] else NULL

    mean_dt$series = ale_label

    if (is_numeric) {
      valid = !is.na(mean_dt$d_l) & !is.na(mean_dt$x_grid)
      n_x = length(unique(mean_dt$x_grid[valid]))
      has_line = n_x >= 2L

      p_ale = ggplot2::ggplot(mean_dt, ggplot2::aes(x = x_grid, y = d_l, color = series)) +
        ggplot2::theme_bw() +
        ggplot2::xlab(feat) +
        ggplot2::ylab(target_feature_name)

      if (has_line) {
        p_ale = p_ale + ggplot2::geom_line(linewidth = 1.2)
      } else {
        p_ale = p_ale + ggplot2::geom_point(size = 1.5)
      }

      if (!is.null(obs_df)) {
        p_ale = p_ale +
          ggplot2::geom_point(
            data = obs_df,
            ggplot2::aes(x = x, y = y),
            inherit.aes = FALSE,
            alpha = 0.3, size = 0.8, color = "black"
          )
      }
    } else {
      if (!is.null(xlim_feat)) {
        level_order = xlim_feat
      } else {
        level_order = as.character(mean_dt$x_grid)
      }
      mean_dt$level = factor(mean_dt$x_grid, levels = level_order)

      valid = !is.na(mean_dt$d_l) & !is.na(mean_dt$level)
      n_x = length(unique(mean_dt$level[valid]))
      has_line = n_x >= 2L

      p_ale = ggplot2::ggplot(mean_dt, ggplot2::aes(x = level, y = d_l, group = 1, color = series)) +
        ggplot2::theme_bw() +
        ggplot2::xlab(feat) +
        ggplot2::ylab(target_feature_name)

      if (has_line) {
        p_ale = p_ale + ggplot2::geom_line(linewidth = 1.2)
      }
      p_ale = p_ale + ggplot2::geom_point(size = 2)

      if (!is.null(obs_df)) {
        obs_df$level = factor(obs_df$x, levels = level_order)
        p_ale = p_ale +
          ggplot2::geom_point(
            data = obs_df,
            ggplot2::aes(x = level, y = y),
            inherit.aes = FALSE,
            alpha = 0.3, size = 1.0, color = "black",
            position = ggplot2::position_jitter(width = 0.12, height = 0)
          )
      }
    }

    p_ale = p_ale +
      ggplot2::scale_color_manual(values = color_vals, name = NULL) +
      ggplot2::labs(color = NULL) +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 8),
        legend.key.size = ggplot2::unit(0.6, "lines"),
        axis.title = ggplot2::element_text(size = 9),
        axis.text = ggplot2::element_text(size = 9)
      )

    leg_theme = if (show_point) {
      ggplot2::theme(
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.box.background = ggplot2::element_rect(fill = "white", color = "grey65", linewidth = 0.2),
        legend.margin = ggplot2::margin(3, 3, 3, 3)
      )
    } else {
      ggplot2::theme(
        legend.position = "inside",
        legend.position.inside = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.background = ggplot2::element_rect(fill = NA, color = NA),
        legend.box.background = ggplot2::element_rect(fill = NA, color = "grey", linewidth = 0.1)
      )
    }
    p_ale = p_ale + leg_theme

    coord_args = list()
    if (is_numeric && !is.null(xlim_feat) && length(xlim_feat) == 2L) {
      coord_args$xlim = xlim_feat
    }
    if (!is.null(ymin) && !is.null(ymax)) {
      ylim_top = ymax
      if (show_point && is.finite(ymin) && is.finite(ymax)) {
        rng = ymax - ymin
        ylim_top = ymax + max(rng, 1e-6) * 0.12
      }
      coord_args$ylim = c(ymin, ylim_top)
    }
    if (length(coord_args)) {
      p_ale = p_ale + do.call(ggplot2::coord_cartesian, coord_args)
    }

    p_ale
  })

  names(out) = feats
  out
}
