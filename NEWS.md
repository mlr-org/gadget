# gadget 0.1.0

- AleStrategy categorical splits now apply ordered-prefix partitions consistently in fitted trees and display category sets for those splits (no issue).
- AleStrategy no longer exposes the retired `with_stab` option; ALE split search now always uses the fast bias-corrected self-gain ranking objective (no issue).
- calculate_ale_fast() now normalizes custom predict_fun outputs like the R ALE path and errors on prediction length mismatches (no issue).
- AleStrategy now uses only the selected split's objective rows when multiple ALE split candidates tie (no issue).
- extract_split_info() keeps categorical split level sets out of the default summary table (no issue).
- plot_tree_pd() now displays categorical split conditions as category sets instead of equality labels (no issue).
- plot_tree_pd() now names returned node plots with actual tree node ids instead of depth-local positions (no issue).
- prepare_split_data_pd() now respects feature_set for PD effects independently from split_feature candidates (no issue).
