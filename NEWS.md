# xplaineff 0.1.0

- Package renamed to xplaineff for CRAN submission (no issue).
- AleStrategy categorical splits now apply ordered-prefix partitions consistently in fitted trees and display category sets for those splits (no issue).
- AleStrategy now accepts a bare prediction function as `model` in the default ALE prediction path (no issue).
- AleStrategy split search now uses the bias-corrected self-gain ranking objective for ALE self-feature splits (no issue).
- AleStrategy now uses only the selected split's objective rows when multiple ALE split candidates tie (no issue).
- calculate_ale() and calculate_ale_fast() now preserve fractional interval bounds for integer features and restore shared prediction scratch data between features (no issue).
- calculate_ale_fast() now normalizes custom predict_fun outputs like the R ALE path and errors on prediction length mismatches (no issue).
- compute_ice_r() now preserves fractional grid values for cached integer features in the PD R backend (no issue).
- extract_split_info() keeps categorical split level sets out of the default summary table (no issue).
- plot_tree_pd() now displays categorical split conditions as category sets instead of equality labels (no issue).
- plot_tree_pd() now names returned node plots with actual tree node ids instead of depth-local positions (no issue).
- prepare_split_data_pd() now infers effect features from precomputed PD/ICE results when feature_set is omitted and still uses all non-target columns as split candidates (no issue).
