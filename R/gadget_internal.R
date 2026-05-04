# Internal helpers and C++ symbols (not for direct use).

risk_from_stats = function(n, s1, s2) ifelse(n <= 1L, 0.0, s2 - (s1 * s1) / n)

assert_ale_effect_list = function(Y, var_name = "Y") {
  required_cols = c("row_id", "interval_index", "d_l", "int_n", "int_s1", "int_s2")
  checkmate::assert_list(Y, min.len = 1, .var.name = var_name)
  checkmate::assert_true(all(mlr3misc::map_lgl(Y, is.data.frame)), .var.name = var_name)
  missing_cols = mlr3misc::map(Y, function(dt) setdiff(required_cols, colnames(dt)))
  bad = Filter(function(x) length(x) > 0, missing_cols)
  if (length(bad) > 0) {
    cli::cli_abort(
      "{.arg {var_name}} elements are missing required columns: {.val {unique(unlist(bad))}}."
    )
  }
}

default_predict_fun = function(model, data) {
  predict_newdata_fast_dispatch(model, data)
}

has_predict_method = function(model, method) {
  !is.null(model[[method]]) && is.function(model[[method]])
}

# Coerce mlr3 Prediction / matrix / vector outputs to a numeric vector (one value per row).
# For multiclass probability PD/ICE on a specific class, pass a custom `predict_fun` in PdStrategy.
# When `expected_n` is set, abort if length differs (including empty vs nonempty rows).
extract_numeric_prediction = function(pred, expected_n = NULL) {
  out = if (inherits(pred, "Prediction")) {
    if (!is.null(pred$response)) {
      as.numeric(pred$response)
    } else if (!is.null(pred$prob)) {
      if (ncol(pred$prob) > 1L) {
        cli::cli_warn(
          "Multiclass model detected: using first class probability column ({colnames(pred$prob)[1L]}) as prediction. Pass {.arg predict_fun} to override."
        )
      }
      as.numeric(pred$prob[, 1L])
    } else {
      cli::cli_abort(
        "{.cls Prediction} object has neither {.field response} nor {.field prob}; cannot extract numeric predictions."
      )
    }
  } else if (is.list(pred) && !is.null(pred$response)) {
    as.numeric(pred$response)
  } else if (is.list(pred) && !is.null(pred$predictions)) {
    as.numeric(pred$predictions)
  } else if (is.data.frame(pred)) {
    as.numeric(pred[[1L]])
  } else if (is.matrix(pred)) {
    as.numeric(pred[, 1L])
  } else {
    as.numeric(pred)
  }

  if (!is.null(expected_n)) {
    if (length(out) != expected_n) {
      cli::cli_abort(
        "Prediction length mismatch: got {length(out)} value{?s} but expected {expected_n} (one per row of {.arg newdata})."
      )
    }
  }
  out
}

# Prediction dispatch priority:
# 1) mlr3 fast path (predict_newdata_fast, requires data.table),
# 2) mlr3 standard path (predict_newdata),
# 3) generic stats::predict fallback.
predict_newdata_fast_dispatch = function(model, newdata) {
  n_rows = nrow(newdata)
  if (has_predict_method(model, "predict_newdata_fast")) {
    newdata_fast = if (data.table::is.data.table(newdata)) newdata else data.table::as.data.table(newdata)
    pred = model$predict_newdata_fast(newdata_fast)
    return(extract_numeric_prediction(pred, expected_n = n_rows))
  }
  if (has_predict_method(model, "predict_newdata")) {
    pred = model$predict_newdata(newdata)
    return(extract_numeric_prediction(pred, expected_n = n_rows))
  }
  pred = tryCatch(
    stats::predict(model, newdata = newdata, type = "response"),
    error = function(e1) {
      tryCatch(
        stats::predict(model, data = newdata, type = "response"),
        error = function(e2) {
          cli::cli_abort(c(
            "Prediction failed for both {.code stats::predict(..., newdata =)} and {.code stats::predict(..., data =)}.",
            "i" = paste0("newdata: ", conditionMessage(e1)),
            "i" = paste0("data: ", conditionMessage(e2))
          ))
        }
      )
    }
  )
  extract_numeric_prediction(pred, expected_n = n_rows)
}

# Wrap each line of a node label (split by \\n) for ggraph tree plots.
wrap_tree_label = function(text, width = 34L) {
  if (length(text) == 0L) {
    return("")
  }
  text = text[[1L]]
  if (is.na(text) || !nzchar(text)) {
    return(text)
  }
  lines = strsplit(text, "\n", fixed = TRUE)[[1L]]
  wrapped = vapply(lines, function(line) {
    if (!nzchar(line)) {
      return(line)
    }
    paste(strwrap(line, width = width, exdent = 2L), collapse = "\n")
  }, character(1L))
  paste(wrapped, collapse = "\n")
}

#' Internal C++ helpers and package symbols
#'
#' Functions and symbols used internally by the package. Not intended for direct use.
#'
#' @name gadget_internal
#' @aliases ale_sweep_cpp calculate_ale_heterogeneity_list_cpp
#'   calculate_ale_heterogeneity_single_cpp re_mean_center_ice_cpp
#'   search_best_split_cpp default_predict_fun predict_newdata_fast_dispatch
#'   extract_numeric_prediction has_predict_method cpp_pd_stack_newdata
#'   risk_from_stats assert_ale_effect_list
#'   d_l interval_index level x x_grid x_left x_right y
#' @keywords internal
NULL
