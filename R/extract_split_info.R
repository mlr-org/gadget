#' Extract Split Information from Tree Structure
#'
#' Given tree (depth-list of Node objects) and optional split_benchmark:
#' flattens nodes; builds one row per node (depth, id, n_obs, split_feature,
#' split_value, int_imp, int_imp_j, etc.); merges timing if split_benchmark has
#' node_id/depth.
#' Returns data frame.
#'
#' @param tree (`list()`) \cr
#'   Depth-based list of Node objects (from \code{convert_tree_to_list}).
#' @param split_benchmark (`data.frame()` or `list()` or `NULL`) \cr
#'   Optional split timing info with columns \code{node_id} and \code{depth}.
#'
#' @return (`data.frame()`) \cr
#'   One row per node: split feature/value, statistics, and (if available) timing.
#'
#' @details
#' This function is used internally by the GadgetTree framework to extract
#' and summarize the structure and statistics of effect-based decision trees.
#' It is useful for interpretation, reporting, and benchmarking.
#'
#' @keywords internal
extract_split_info = function(tree, split_benchmark = NULL) {
  # Collect all importance$imp_j names
  all_intimp_names = unique(unlist(
    lapply(tree, function(depth) {
      lapply(depth, function(node) {
        if (!is.null(node$importance) && !is.null(node$importance$imp_j)) names(node$importance$imp_j) else NULL
      })
    })
  ))
  all_intimp_names = all_intimp_names[!is.na(all_intimp_names)]

  rows = lapply(unlist(tree, recursive = FALSE), function(node) {
    if (is.null(node)) {
      return(NULL)
    }
    n_obs = if (is.null(node$subset_idx) || length(node$subset_idx) == 0) 0 else length(node$subset_idx)
    is_final = isTRUE(node$improvement_met) | isTRUE(node$stop_criterion_met) |
      is.null(node$children) || (is.list(node$children) && all(mlr3misc::map_lgl(node$children, is.null)))
    row = data.frame(
      depth = as.integer(node$depth),
      id = as.integer(node$id),
      n_obs = as.integer(n_obs),
      node_type = if (is.null(node$id)) NA_character_ else
        if (node$id == 1) "root" else if (node$id %% 2 == 0) "left" else "right",
      split_feature = if (is.null(node$split)) NA_character_ else as.character(node$split$feature),
      split_value = if (is.null(node$split)) NA else node$split$value,
      node_objective = if (is.null(node$objective$value)) NA else node$objective$value,
      int_imp = if (is.null(node$importance)) NA else round(node$importance$imp, 2),
      split_feature_parent = if (is.null(node$parent) || is.null(node$parent$split_feature))
        NA_character_ else as.character(node$parent$split_feature),
      split_value_parent = if (is.null(node$parent) || is.null(node$parent$split_value))
        NA else node$parent$split_value,
      objective_value_parent = if (is.null(node$parent) || is.null(node$parent$objective_value))
        NA else node$parent$objective_value,
      int_imp_parent = if (is.null(node$parent) || is.null(node$parent$int_imp)) NA else round(node$parent$int_imp, 2),
      is_final = is_final,
      stringsAsFactors = FALSE
    )
    # Ensure all int_imp.* fields exist
    for (nm in all_intimp_names) {
      row[[paste0("int_imp_", nm)]] =
        if (!is.null(node$importance$imp_j) && nm %in% names(node$importance$imp_j))
          round(node$importance$imp_j[[nm]], 2) else NA
    }
    row
  })
  df_split = do.call(rbind, rows)
  rownames(df_split) = NULL

  # Field order: place int_imp.* fields after int_imp
  all_cols = names(df_split)
  intimp_pos = which(all_cols == "int_imp")
  intimp_dot_cols = grep("^int_imp_", all_cols, value = TRUE)
  cols_wo_dot = setdiff(all_cols, intimp_dot_cols)
  new_order = append(cols_wo_dot, intimp_dot_cols, after = intimp_pos)
  df_split = df_split[, new_order]

  # Merge split_benchmark info if provided
  if (!is.null(split_benchmark) && length(split_benchmark) > 0) {
    # Ensure split_benchmark is a data.frame
    if (is.list(split_benchmark) && !is.data.frame(split_benchmark)) {
      split_benchmark = do.call(rbind, lapply(split_benchmark, as.data.frame))
    }
    # Check if split_benchmark is valid data frame with rows
    if (is.data.frame(split_benchmark) && nrow(split_benchmark) > 0) {
      df_split = merge(df_split, split_benchmark,
        by.x = c("id", "depth"), by.y = c("node_id", "depth"),
        all.x = TRUE, sort = TRUE)
    }
  }
  df_split
}
