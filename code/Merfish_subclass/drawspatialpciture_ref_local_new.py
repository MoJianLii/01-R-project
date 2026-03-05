# -*- coding: utf-8 -*-
"""
Mouse1 strict1to1 (exclude layer1) reference-style drawing script.

Notes
1) Keep original reference draw script untouched.
2) First two panels follow reference plotting style:
   - scatter marker '.'
   - point size 0.5
   - xlim(0,11), ylim(11,0), equal axis, no ticks
3) Data source for first two panels uses local base_dir:
   E:\\zaw\\2602\\cell_metadata1
4) Third panel overlays strict 1to1 pair regions on wholebrain (E:\\zaw\\Total).
"""

from __future__ import print_function

import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ===================== Paths =====================
BASE_DIR = r"E:\zaw\2602\cell_metadata1"  # user requested local reference source
WHOLEBRAIN_DIR = r"E:\zaw\Total"

BASE_OUT_DIR = r"E:\zaw\2511\mouseMerfish_zhuang_subclass\ws0.4_ss0.02\Mouse1_strict1to1_PLOTS_ROI_MATCHED"
TABLES_DIR = os.path.join(BASE_OUT_DIR, "tables")
PAIRS_FILE = os.path.join(TABLES_DIR, "Mouse1_strict_pairs_excludeLayer1.tsv")

CLUSTER_FILE_OVER3 = r"E:\zaw\2511\mouseMerfish_zhuang_subclass\ws0.4_ss0.02\mouse_subclass_cluster_total_over3.txt"
CLUSTER_FILE_FALLBACK = r"E:\zaw\2511\mouseMerfish_zhuang_subclass\ws0.4_ss0.02\mouse_subclass_cluster_total.txt"

OUT_DIR = os.path.join(BASE_OUT_DIR, "mouse1_excludeLayer1_refstyle_from_local")
OUT_COMBINED = os.path.join(OUT_DIR, "combined")
OUT_P1 = os.path.join(OUT_DIR, "panel1_subclass")
OUT_P2 = os.path.join(OUT_DIR, "panel2_substructure")
OUT_P3 = os.path.join(OUT_DIR, "panel3_overlay")
OUT_P3_ANN = os.path.join(OUT_DIR, "panel3_overlay_with_cluster_stats")
OUT_COMBINED_ANN = os.path.join(OUT_DIR, "combined_with_cluster_stats")


# ===================== Params =====================
KEEP_LAYERS = ["2/3", "4", "5", "6a", "6b"]
ID_COLS = ["Glut_Neruon_cell_ids", "GABA_Neruon_cell_ids", "Non_Neruon_cell_ids"]
TOP_SUBSTRUCTURE = 28
TOP_SUBCLASS = 35
MAX_BG_POINTS = 280000
MAX_HI_POINTS = 320000

# Optional quick test: set env DRAW_MAX_PAIRS=20
DRAW_MAX_PAIRS = os.environ.get("DRAW_MAX_PAIRS", "").strip()
DRAW_MAX_PAIRS = int(DRAW_MAX_PAIRS) if DRAW_MAX_PAIRS.isdigit() else None


def ensure_dirs(paths):
    for p in paths:
        if not os.path.exists(p):
            os.makedirs(p)


def safe_token(x):
    s = "NA" if x is None else str(x)
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s if s else "NA"


def normalize_layer(x):
    s = str(x).strip().lower()
    s = re.sub(r"^layer\s*", "", s)
    s = re.sub(r"\s+", "", s)
    if s in ["23", "2-3", "2_3"]:
        s = "2/3"
    return s


def pick_cluster_file():
    if os.path.exists(CLUSTER_FILE_OVER3):
        return CLUSTER_FILE_OVER3
    if os.path.exists(CLUSTER_FILE_FALLBACK):
        return CLUSTER_FILE_FALLBACK
    raise IOError("No cluster_total file found (over3/fallback).")


def parse_ids(raw):
    if raw is None:
        return []
    try:
        if np.isnan(raw):
            return []
    except Exception:
        pass
    txt = str(raw)
    if (not txt) or (txt.lower() == "nan"):
        return []
    return re.findall(r"\d+", txt)


def load_pairs():
    if not os.path.exists(PAIRS_FILE):
        raise IOError("Pairs file not found: {}".format(PAIRS_FILE))
    df = pd.read_csv(PAIRS_FILE, sep="\t", dtype=str, low_memory=False)
    req = ["pair_id", "slide_use", "label1", "label2"]
    miss = [c for c in req if c not in df.columns]
    if miss:
        raise ValueError("Missing required columns in pair file: {}".format(miss))

    df["slide_use"] = df["slide_use"].astype(str).str.strip()
    df = df[df["slide_use"].str.startswith("C57BL6J-1.", na=False)].copy()

    if "layer_use" in df.columns:
        df["layer_use"] = df["layer_use"].map(normalize_layer)
    elif "layer_f" in df.columns:
        df["layer_use"] = df["layer_f"].map(normalize_layer)
    else:
        df["layer_use"] = ""

    df = df[df["layer_use"].isin(KEEP_LAYERS)].copy()
    df["pair_id"] = pd.to_numeric(df["pair_id"], errors="coerce")
    df = df.dropna(subset=["pair_id"])
    df["pair_id"] = df["pair_id"].astype(int)
    df = df.sort_values(["slide_use", "pair_id"]).drop_duplicates(
        subset=["pair_id", "slide_use", "label1", "label2"]
    )
    if DRAW_MAX_PAIRS is not None:
        df = df.head(DRAW_MAX_PAIRS).copy()
    return df.reset_index(drop=True)


def build_label_to_ids(cluster_file, labels):
    head = pd.read_csv(cluster_file, sep="\t", nrows=0)
    needed = set(["slide", "layer", "subclass", "merge_regions"] + ID_COLS)
    miss = [c for c in needed if c not in head.columns]
    if miss:
        raise ValueError("cluster_total missing columns: {}".format(miss))

    usecols = list(needed) + (["label"] if "label" in head.columns else [])
    clu = pd.read_csv(cluster_file, sep="\t", dtype=str, usecols=usecols, low_memory=False)
    if "label" not in clu.columns:
        clu["label"] = (
            clu["slide"].astype(str).str.strip()
            + "_"
            + clu["layer"].astype(str).str.strip()
            + "_"
            + clu["subclass"].astype(str).str.strip()
            + "_"
            + clu["merge_regions"].astype(str).str.strip()
        )
    clu["label"] = clu["label"].astype(str).str.strip()
    label_set = set([str(x).strip() for x in labels if str(x).strip()])
    clu = clu[clu["label"].isin(label_set)].copy()

    out = {}
    for _, row in clu.iterrows():
        ids = []
        for c in ID_COLS:
            ids.extend(parse_ids(row.get(c)))
        out[row["label"]] = list(dict.fromkeys(ids))
    return out


def _pick_first_col(columns, candidates):
    for c in candidates:
        if c in columns:
            return c
    return None


def _safe_int(v):
    try:
        return int(float(v))
    except Exception:
        return None


def build_label_to_cluster_stats(cluster_file, labels):
    """
    Build label -> {glu_num, gaba_num, ei_ratio}
    """
    head = pd.read_csv(cluster_file, sep="\t", nrows=0)
    cols = head.columns.tolist()

    needed = set(["slide", "layer", "subclass", "merge_regions"] + ID_COLS)
    miss = [c for c in needed if c not in cols]
    if miss:
        raise ValueError("cluster_total missing columns for stats: {}".format(miss))

    glu_num_col = _pick_first_col(cols, [
        "Glut_Neruon_cell_ids_num", "Glut_Neuron_cell_ids_num", "GLU_cell_num", "cluster_GLU_cell_num"
    ])
    gaba_num_col = _pick_first_col(cols, [
        "GABA_Neruon_cell_ids_num", "GABA_Neuron_cell_ids_num", "GABA_cell_num", "cluster_GABA_cell_num"
    ])

    usecols = list(needed) + (["label"] if "label" in cols else [])
    if glu_num_col:
        usecols.append(glu_num_col)
    if gaba_num_col:
        usecols.append(gaba_num_col)
    usecols = list(dict.fromkeys(usecols))

    clu = pd.read_csv(cluster_file, sep="\t", dtype=str, usecols=usecols, low_memory=False)
    if "label" not in clu.columns:
        clu["label"] = (
            clu["slide"].astype(str).str.strip()
            + "_"
            + clu["layer"].astype(str).str.strip()
            + "_"
            + clu["subclass"].astype(str).str.strip()
            + "_"
            + clu["merge_regions"].astype(str).str.strip()
        )
    clu["label"] = clu["label"].astype(str).str.strip()
    label_set = set([str(x).strip() for x in labels if str(x).strip()])
    clu = clu[clu["label"].isin(label_set)].copy()

    out = {}
    for _, row in clu.iterrows():
        glu_n = _safe_int(row.get(glu_num_col)) if glu_num_col else None
        gaba_n = _safe_int(row.get(gaba_num_col)) if gaba_num_col else None
        if glu_n is None:
            glu_n = len(parse_ids(row.get("Glut_Neruon_cell_ids")))
        if gaba_n is None:
            gaba_n = len(parse_ids(row.get("GABA_Neruon_cell_ids")))
        ei = (float(glu_n) / float(gaba_n)) if (gaba_n is not None and gaba_n > 0) else np.nan
        out[row["label"]] = {
            "glu_num": int(glu_n) if glu_n is not None else 0,
            "gaba_num": int(gaba_n) if gaba_n is not None else 0,
            "ei_ratio": ei,
        }
    return out


def load_local_cell_metadata(base_dir):
    csv_files = sorted(
        [f for f in os.listdir(base_dir) if f.endswith(".csv") and f.lower().startswith("cell_metadata")]
    )
    if not csv_files:
        raise IOError("No cell_metadata*.csv found in {}".format(base_dir))

    cell = {}
    for i, f in enumerate(csv_files, start=1):
        key = "Zhuang-ABCA-{}".format(i)
        fp = os.path.join(base_dir, f)
        df = pd.read_csv(fp, dtype={"cell_label": str}, low_memory=False)
        req = ["cell_label", "brain_section_label", "x", "y"]
        miss = [c for c in req if c not in df.columns]
        if miss:
            raise ValueError("{} missing columns {}".format(fp, miss))
        df["cell_label"] = df["cell_label"].astype(str).str.strip()
        df["brain_section_label"] = df["brain_section_label"].astype(str).str.strip()
        df["x"] = pd.to_numeric(df["x"], errors="coerce")
        df["y"] = pd.to_numeric(df["y"], errors="coerce")
        df = df.dropna(subset=["x", "y"]).copy()
        cell[key] = df
        print("{}: Number of cells = {}, Number of sections = {}".format(
            key, len(df), df["brain_section_label"].nunique()
        ))
    return cell


_WHOLEBRAIN_CACHE = {}


def load_wholebrain_slide(slide):
    if slide in _WHOLEBRAIN_CACHE:
        return _WHOLEBRAIN_CACHE[slide]

    fp1 = os.path.join(WHOLEBRAIN_DIR, "{}.txt".format(slide))
    fp2 = os.path.join(WHOLEBRAIN_DIR, "{}.tsv".format(slide))
    if os.path.exists(fp1):
        fp = fp1
    elif os.path.exists(fp2):
        fp = fp2
    else:
        raise IOError("Wholebrain file missing for {}".format(slide))

    cols = pd.read_csv(fp, sep="\t", nrows=0).columns.tolist()
    req = ["cell_label", "x", "y"]
    miss = [c for c in req if c not in cols]
    if miss:
        raise ValueError("{} missing required columns {}".format(fp, miss))
    usecols = req + [c for c in ["class", "subclass", "ccf_region_name"] if c in cols]
    df = pd.read_csv(fp, sep="\t", dtype={"cell_label": str}, usecols=usecols, low_memory=False)
    df["cell_label"] = df["cell_label"].astype(str).str.strip()
    df["x"] = pd.to_numeric(df["x"], errors="coerce")
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    df = df.dropna(subset=["x", "y"]).copy()
    if "class" not in df.columns:
        df["class"] = "Unknown"
    if "subclass" not in df.columns:
        df["subclass"] = "Unknown"
    if "ccf_region_name" not in df.columns:
        df["ccf_region_name"] = "Unknown"
    _WHOLEBRAIN_CACHE[slide] = df
    return df


def c57_to_zhuang(slide):
    # C57BL6J-1.079 -> Zhuang-ABCA-1.079
    m = re.match(r"^C57BL6J-(\d+)\.(\d+)$", str(slide))
    if not m:
        return None, None
    mouse_id = m.group(1)
    sec_id = m.group(2)
    dkey = "Zhuang-ABCA-{}".format(mouse_id)
    sec = "{}.{}".format(dkey, sec_id)
    return dkey, sec


def make_discrete_colors(values, cmap_name="tab20"):
    uniq = sorted(pd.Series(values).fillna("Unknown").astype(str).unique().tolist())
    n = max(1, len(uniq))
    cmap = plt.get_cmap(cmap_name, n)
    return dict((k, cmap(i)) for i, k in enumerate(uniq))


# ---------- reference-like plotting primitives ----------
def subplot_section(ax, xx, yy, cc=None):
    ax.scatter(xx, yy, s=0.5, color=cc, marker=".")
    ax.set_ylim(11, 0)
    ax.set_xlim(0, 11)
    ax.axis("equal")
    ax.set_xticks([])
    ax.set_yticks([])


def sample_df(df, nmax, seed=1):
    if len(df) <= nmax:
        return df
    return df.sample(n=nmax, random_state=seed)


def draw_panel1_subclass(ax, section_df, title):
    # panel1 switched to subclass annotation
    tmp = section_df.copy()
    top = tmp["subclass"].astype(str).value_counts().head(TOP_SUBCLASS).index
    tmp["subclass_show"] = np.where(tmp["subclass"].astype(str).isin(top), tmp["subclass"], "Other")
    color_map = make_discrete_colors(tmp["subclass_show"].values, cmap_name="tab20")
    cc = tmp["subclass_show"].astype(str).map(color_map).values
    subplot_section(ax, tmp["x"].values, tmp["y"].values, cc=cc)
    ax.set_title(title, fontsize=12)


def draw_panel2_substructure(ax, section_df, title):
    tmp = section_df.copy()
    # avoid too many legends/colors: keep top regions, others -> Other
    top = tmp["ccf_region_name"].astype(str).value_counts().head(TOP_SUBSTRUCTURE).index
    tmp["sub_show"] = np.where(tmp["ccf_region_name"].astype(str).isin(top), tmp["ccf_region_name"], "Other")
    color_map = make_discrete_colors(tmp["sub_show"].values, cmap_name="turbo")
    cc = tmp["sub_show"].astype(str).map(color_map).values
    subplot_section(ax, tmp["x"].values, tmp["y"].values, cc=cc)
    ax.set_title(title, fontsize=12)


def draw_panel3_overlay(ax, whole_df, ids_glu, ids_gab, title, anno_text=None):
    bg = sample_df(whole_df[["x", "y"]].copy(), MAX_BG_POINTS, seed=1)
    # slightly darker background for better contrast
    ax.scatter(bg["x"], bg["y"], s=8, color="#7A7A7A", marker=".", alpha=1.0, linewidths=0)

    id_union = set(ids_glu).union(set(ids_gab))
    hi = whole_df[whole_df["cell_label"].isin(id_union)][["cell_label", "x", "y"]].copy()
    hi = sample_df(hi, MAX_HI_POINTS, seed=1)

    set_glu = set(ids_glu)
    set_gab = set(ids_gab)
    hi["grp"] = np.where(
        hi["cell_label"].isin(set_glu) & hi["cell_label"].isin(set_gab),
        "Overlap",
        np.where(hi["cell_label"].isin(set_gab), "GabaOnly", "GlutOnly"),
    )
    # deepen highlight colors to make ROI stand out
    colors = {"GlutOnly": "#FF0000", "Overlap": "#0033FF", "GabaOnly": "#B3131B"}
    for g in ["GlutOnly", "Overlap", "GabaOnly"]:
        part = hi[hi["grp"] == g]
        if len(part) > 0:
            ax.scatter(
                part["x"], part["y"],
                s=10.0, color=colors[g], marker="o",
                alpha=1.0, linewidths=0, zorder=3
            )

    ax.set_ylim(11, 0)
    ax.set_xlim(0, 11)
    ax.axis("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=12)
    if anno_text:
        ax.text(
            0.02, 0.98, anno_text,
            transform=ax.transAxes, ha="left", va="top",
            fontsize=9, color="black",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#666666", alpha=0.82)
        )

    return {
        "n_glu_only_drawn": int((hi["grp"] == "GlutOnly").sum()),
        "n_overlap_drawn": int((hi["grp"] == "Overlap").sum()),
        "n_gab_only_drawn": int((hi["grp"] == "GabaOnly").sum()),
    }


def main():
    ensure_dirs([OUT_DIR, OUT_COMBINED, OUT_P1, OUT_P2, OUT_P3, OUT_P3_ANN, OUT_COMBINED_ANN, TABLES_DIR])

    pairs = load_pairs()
    if pairs.empty:
        raise RuntimeError("No mouse1 strict pairs after excluding layer1.")

    cluster_file = pick_cluster_file()
    label_to_ids = build_label_to_ids(cluster_file, set(pairs["label1"]).union(set(pairs["label2"])))
    label_to_stats = build_label_to_cluster_stats(cluster_file, set(pairs["label1"]).union(set(pairs["label2"])))
    local_cell = load_local_cell_metadata(BASE_DIR)

    logs = []
    for _, row in pairs.iterrows():
        pair_id = int(row["pair_id"])
        slide = str(row["slide_use"])
        layer = str(row.get("layer_use", "NA"))
        label1 = str(row["label1"])
        label2 = str(row["label2"])
        ov = str(row.get("ov01", "NA"))

        dkey, sec = c57_to_zhuang(slide)
        if dkey is None or dkey not in local_cell:
            logs.append({
                "pair_id": pair_id, "slide": slide, "layer_use": layer,
                "label1": label1, "label2": label2, "status": "NoLocalDatasetKey"
            })
            continue

        try:
            whole_df = load_wholebrain_slide(slide)
        except Exception as e:
            logs.append({
                "pair_id": pair_id, "slide": slide, "layer_use": layer,
                "label1": label1, "label2": label2, "status": "WholebrainLoadError: {}".format(e)
            })
            continue

        sec_df = local_cell[dkey]
        sec_df = sec_df[sec_df["brain_section_label"] == sec][["cell_label", "x", "y"]].copy()
        if sec_df.empty:
            logs.append({
                "pair_id": pair_id, "slide": slide, "layer_use": layer,
                "label1": label1, "label2": label2, "status": "NoSectionInLocalMeta"
            })
            continue

        # enrich first two panels with subclass/substructure by cell_label
        sec_df = sec_df.merge(
            whole_df[["cell_label", "subclass", "ccf_region_name"]] if "subclass" in whole_df.columns else whole_df[["cell_label", "ccf_region_name"]],
            on="cell_label", how="left"
        )
        if "subclass" not in sec_df.columns:
            sec_df["subclass"] = "Unknown"
        sec_df["subclass"] = sec_df["subclass"].fillna("Unknown")
        sec_df["ccf_region_name"] = sec_df["ccf_region_name"].fillna("Unknown")

        ids_glu = label_to_ids.get(label1, [])
        ids_gab = label_to_ids.get(label2, [])
        if len(ids_glu) == 0 or len(ids_gab) == 0:
            logs.append({
                "pair_id": pair_id, "slide": slide, "layer_use": layer,
                "label1": label1, "label2": label2, "status": "MissingAllIDs",
                "n_glu_ids": len(ids_glu), "n_gab_ids": len(ids_gab)
            })
            continue

        s1 = label_to_stats.get(label1, {"glu_num": 0, "gaba_num": 0, "ei_ratio": np.nan})
        s2 = label_to_stats.get(label2, {"glu_num": 0, "gaba_num": 0, "ei_ratio": np.nan})
        c1_ei = "NA" if np.isnan(s1["ei_ratio"]) else "{:.2f}".format(s1["ei_ratio"])
        c2_ei = "NA" if np.isnan(s2["ei_ratio"]) else "{:.2f}".format(s2["ei_ratio"])
        anno_text = (
            "cluster1  GLU:{}, GABA:{}, EI:{}\n"
            "cluster2  GLU:{}, GABA:{}, EI:{}"
        ).format(s1["glu_num"], s1["gaba_num"], c1_ei, s2["glu_num"], s2["gaba_num"], c2_ei)

        stem = "{}_pair{:04d}_layer{}".format(safe_token(slide), pair_id, safe_token(layer))
        try:
            fig, ax = plt.subplots(1, 3, figsize=(18.8, 6.6))
            draw_panel1_subclass(
                ax[0], sec_df,
                "{} Subclass\nCell Number: {}".format(sec, len(sec_df))
            )
            draw_panel2_substructure(
                ax[1], sec_df,
                "{} Substructure\nCell Number: {}".format(sec, len(sec_df))
            )
            ov_stat = draw_panel3_overlay(
                ax[2], whole_df, ids_glu, ids_gab,
                "{} layer {} cluster1:{} cluster2:{}".format(
                    slide, layer, len(ids_glu), len(ids_gab)
                ),
                anno_text=anno_text
            )
            fig.tight_layout()
            fig.savefig(os.path.join(OUT_COMBINED_ANN, "{}_3panel_refstyle_withStats.png".format(stem)), dpi=350, facecolor="white")
            fig.savefig(os.path.join(OUT_COMBINED_ANN, "{}_3panel_refstyle_withStats.pdf".format(stem)), dpi=350, facecolor="white")
            plt.close(fig)

            # also save three individual panels
            p1, a1 = plt.subplots(1, 1, figsize=(6.2, 6.2))
            draw_panel1_subclass(a1, sec_df, "{} Subclass".format(sec))
            p1.tight_layout()
            p1.savefig(os.path.join(OUT_P1, "{}_P1_subclass.png".format(stem)), dpi=350, facecolor="white")
            plt.close(p1)

            p2, a2 = plt.subplots(1, 1, figsize=(6.2, 6.2))
            draw_panel2_substructure(a2, sec_df, "{} Substructure".format(sec))
            p2.tight_layout()
            p2.savefig(os.path.join(OUT_P2, "{}_P2_substructure.png".format(stem)), dpi=350, facecolor="white")
            plt.close(p2)

            p3, a3 = plt.subplots(1, 1, figsize=(6.2, 6.2))
            draw_panel3_overlay(
                a3, whole_df, ids_glu, ids_gab,
                "{} Pair {}".format(slide, pair_id),
                anno_text=anno_text
            )
            p3.tight_layout()
            p3.savefig(os.path.join(OUT_P3_ANN, "{}_P3_overlay_withStats.png".format(stem)), dpi=350, facecolor="white")
            plt.close(p3)

            logs.append({
                "pair_id": pair_id,
                "slide": slide,
                "zhuang_section": sec,
                "layer_use": layer,
                "label1": label1,
                "label2": label2,
                "ov01": ov,
                "n_section_cells": int(len(sec_df)),
                "n_total_cells_slide": int(len(whole_df)),
                "n_glu_ids": int(len(ids_glu)),
                "n_gab_ids": int(len(ids_gab)),
                "cluster1_glu_num": int(s1["glu_num"]),
                "cluster1_gaba_num": int(s1["gaba_num"]),
                "cluster1_ei_ratio": c1_ei,
                "cluster2_glu_num": int(s2["glu_num"]),
                "cluster2_gaba_num": int(s2["gaba_num"]),
                "cluster2_ei_ratio": c2_ei,
                "n_glu_only_drawn": ov_stat["n_glu_only_drawn"],
                "n_overlap_drawn": ov_stat["n_overlap_drawn"],
                "n_gab_only_drawn": ov_stat["n_gab_only_drawn"],
                "status": "OK",
                "file_stem": stem,
            })
        except Exception as e:
            logs.append({
                "pair_id": pair_id, "slide": slide, "layer_use": layer,
                "label1": label1, "label2": label2,
                "status": "RenderError: {}".format(e)
            })

    log_df = pd.DataFrame(logs)
    log_path = os.path.join(TABLES_DIR, "Mouse1_excludeLayer1_refstyle_local_draw_log.tsv")
    log_df.to_csv(log_path, sep="\t", index=False)

    if len(log_df) > 0:
        log_df["ok"] = log_df["status"].eq("OK")
        summary = log_df.groupby("layer_use", dropna=False).agg(
            n_pairs=("pair_id", "count"),
            n_ok=("ok", "sum"),
        ).reset_index()
        summary["pct_ok"] = np.where(summary["n_pairs"] > 0, 100.0 * summary["n_ok"] / summary["n_pairs"], np.nan)
    else:
        summary = pd.DataFrame(columns=["layer_use", "n_pairs", "n_ok", "pct_ok"])

    summary_path = os.path.join(TABLES_DIR, "Mouse1_excludeLayer1_refstyle_local_summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)

    print("Done. Output dir: {}".format(OUT_DIR))
    print("Draw log: {}".format(log_path))
    print("Summary: {}".format(summary_path))


if __name__ == "__main__":
    main()
