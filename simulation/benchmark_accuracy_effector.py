#!/usr/bin/env python3
"""Accuracy benchmark (effector): split recovery vs simulation/efficiency_benchmark_plan.md §9.

Requires shared CSVs from: Rscript simulation/generate_accuracy_data.R
Run: python simulation/benchmark_accuracy_effector.py [--datadir DIR] [--outdir DIR]
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd

try:
    import effector
    from effector.space_partitioning import Best
except ImportError:
    sys.stderr.write("effector not installed. Run: pip install effector\n")
    sys.exit(1)


N_VEC = [200, 500, 1000, 5000]
D_VEC = [5, 10, 20]
VARIANTS = ["num_0", "num_05", "cat"]
N_SEEDS = 30
FOI = 0
N_GRID = 20
N_INTERVALS = 20


def parse_int_vec(s: str) -> list:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def parse_str_vec(s: str) -> list:
    return [x.strip() for x in s.split(",") if x.strip()]


def true_threshold(variant: str) -> float | None:
    if variant == "num_0":
        return 0.0
    if variant == "num_05":
        return 0.5
    return None


def load_csv(datadir: str, N: int, D: int, variant: str, seed: int) -> tuple[pd.DataFrame, np.ndarray]:
    fn = os.path.join(datadir, "acc_N{}_D{}_{}_seed{}.csv".format(N, D, variant, seed))
    if not os.path.isfile(fn):
        sys.stderr.write("Missing dataset {} — run: Rscript simulation/generate_accuracy_data.R\n".format(fn))
        sys.exit(1)
    df = pd.read_csv(fn)
    if variant == "cat":
        # Character A/B -> 0/1 in [:,2] for numpy-only X while preserving two levels for utils.get_feature_types
        xv = df["x3"].astype(str).values
        x3n = np.where(xv == "B", 1.0, 0.0)
        df_num = df.copy()
        df_num["x3"] = x3n
        X = df_num[["x{}".format(i + 1) for i in range(D)]].values.astype(np.float64)
    else:
        X = df[["x{}".format(i + 1) for i in range(D)]].values.astype(np.float64)
    return df, X


def toy_predict(variant: str):

    def pred(Xm: np.ndarray) -> np.ndarray:
        x1, x2, x3 = Xm[:, 0], Xm[:, 1], Xm[:, 2]
        if variant == "cat":
            m = x3 >= 0.5
        else:
            th = 0.5 if variant == "num_05" else 0.0
            m = x3 > th
        return 5.0 * x1 + 5.0 * x2 + np.where(m, 10.0 * x1 - 10.0 * x2, 0.0)

    return pred


def axis_limits_from_X(X: np.ndarray) -> np.ndarray:
    """Shape (2, D): row0 mins, row1 maxs — effector convention [lower, upper] per feature."""
    lo = X.min(axis=0)
    hi = X.max(axis=0)
    # small padding for numerical stability
    pad = 1e-6 * (hi - lo + 1e-9)
    return np.vstack([lo - pad, hi + pad])


def true_left_mask(X: np.ndarray, variant: str) -> np.ndarray:
    x3 = X[:, 2]
    if variant == "cat":
        return x3 < 0.5
    th = true_threshold(variant)
    return x3 <= th + 1e-12


def first_split_props(part) -> tuple:
    """Return (foc_index, foc_split_position, foc_type). Compatible with effector versions
    that expose either `important_splits` (list of dicts) or `splits_tree` (Tree nodes)."""
    if part is None:
        return None, None, None
    imp = getattr(part, "important_splits", None)
    if isinstance(imp, list) and len(imp) > 0:
        s0 = imp[0]
        return s0["foc_index"], s0["foc_split_position"], s0["foc_type"]
    tree = getattr(part, "splits_tree", None)
    if tree is not None and len(getattr(tree, "nodes", [])) > 1:
        for n in tree.nodes:
            info = getattr(n, "info", None) or {}
            if info.get("level") == 1 and "foc_index" in info:
                return info["foc_index"], info["foc_split_position"], info.get("foc_type")
    return None, None, None


def split_feat_hit(foc_index: int | None, variant: str) -> bool:
    return foc_index == 2


def split_point_mae(foc_index, position, foc_type, variant: str) -> float | None:
    if variant == "cat":
        return None
    if foc_index != 2 or position is None:
        return None
    th = true_threshold(variant)
    return abs(float(position) - float(th))


def node_acc(X: np.ndarray, variant: str, foc_index, position, foc_type) -> float:
    if foc_index is None:
        return float("nan")
    tl = true_left_mask(X, variant)
    if foc_type == "cat":
        pl = X[:, foc_index] == position
        p2 = X[:, foc_index] != position
    else:
        pl = X[:, foc_index] < position
        p2 = X[:, foc_index] >= position
    a1 = np.mean(pl == tl)
    a2 = np.mean(p2 == tl)
    return float(max(a1, a2))


def run_regional_pdp(X: np.ndarray, predict, axis_limits: np.ndarray, feature_names: list, spp: Best):
    m = effector.RegionalPDP(
        data=X,
        model=predict,
        axis_limits=axis_limits,
        nof_instances="all",
        feature_names=feature_names,
    )
    m.fit(
        features=[FOI],
        candidate_conditioning_features="all",
        space_partitioner=spp,
    )
    return m.partitioners["feature_{}".format(FOI)]


def run_regional_ale(X: np.ndarray, predict, axis_limits: np.ndarray, feature_names: list, spp: Best):
    bm = effector.axis_partitioning.Fixed(nof_bins=N_INTERVALS, min_points_per_bin=0)
    m = effector.RegionalALE(
        data=X,
        model=predict,
        axis_limits=axis_limits,
        nof_instances="all",
        feature_names=feature_names,
    )
    m.fit(
        features=[FOI],
        candidate_conditioning_features="all",
        space_partitioner=spp,
        binning_method=bm,
        points_for_mean_heterogeneity=N_INTERVALS,
    )
    return m.partitioners["feature_{}".format(FOI)]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--datadir", default="simulation/data/accuracy")
    p.add_argument("--outdir", default="simulation/results/accuracy")
    p.add_argument("--n-seeds", type=int, default=N_SEEDS)
    p.add_argument("--N-vec", default=None)
    p.add_argument("--D-vec", default=None)
    p.add_argument("--variants", default=None)
    args = p.parse_args()

    n_vec = parse_int_vec(args.N_vec) if args.N_vec else N_VEC
    d_vec = parse_int_vec(args.D_vec) if args.D_vec else D_VEC
    variants = parse_str_vec(args.variants) if args.variants else VARIANTS
    n_seeds = args.n_seeds

    os.makedirs(args.outdir, exist_ok=True)
    rows = []
    spp = Best(
        max_depth=1,
        min_samples_leaf=50,
        min_heterogeneity_decrease_pcg=0.01,
        numerical_features_grid_size=40,
    )

    for variant in variants:
        predict = toy_predict(variant)
        for N in n_vec:
            for D in d_vec:
                feat_names = ["x{}".format(i + 1) for i in range(D)]
                for s in range(1, n_seeds + 1):
                    seed = 1000 + s
                    _, X = load_csv(args.datadir, N, D, variant, seed)
                    al = axis_limits_from_X(X)
                    part_pdp = run_regional_pdp(X, predict, al, feat_names, spp)
                    fi, pos, ft = first_split_props(part_pdp)
                    hit = split_feat_hit(fi, variant)
                    mae = split_point_mae(fi, pos, ft, variant)
                    nacc = node_acc(X, variant, fi, pos, ft)
                    rows.append(
                        {
                            "package": "effector",
                            "method": "effector_rpdp",
                            "variant": variant,
                            "N": N,
                            "D": D,
                            "seed": seed,
                            "split_feat_correct": hit,
                            "split_pt_error": mae if mae is not None else "",
                            "node_acc": nacc,
                            "effect_mse_node1": "",
                            "effect_mse_node2": "",
                        }
                    )

                    part_ale = run_regional_ale(X, predict, al, feat_names, spp)
                    fi, pos, ft = first_split_props(part_ale)
                    hit = split_feat_hit(fi, variant)
                    mae = split_point_mae(fi, pos, ft, variant)
                    nacc = node_acc(X, variant, fi, pos, ft)
                    rows.append(
                        {
                            "package": "effector",
                            "method": "effector_rale",
                            "variant": variant,
                            "N": N,
                            "D": D,
                            "seed": seed,
                            "split_feat_correct": hit,
                            "split_pt_error": mae if mae is not None else "",
                            "node_acc": nacc,
                            "effect_mse_node1": "",
                            "effect_mse_node2": "",
                        }
                    )

    outp = os.path.join(args.outdir, "accuracy_effector.csv")
    with open(outp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("Written: {}".format(outp))


if __name__ == "__main__":
    main()
