# Anne Nguyen (an875)
# Last Update: 02/03/2026 
# Script for classifying molecules from update_extractingMolecules_MutilContig_V7_SingleSampleNew_tosent_fixed.py's OUTPUT .CSV FILE into 4 main categories
"""
Molecule Classification into 4 categories (Normal/Fused x With/Without Telomere).

Key Logic Choices:
- Do NOT drop rows where Contig_Site == 0 and Contig_Position == 0:
  because those rows still represent real molecule labels (including molecule ends and telomere-channel labels)
- For telomere-positive molecules, choose the LabelChannel=1 label closest to the telomeric end of the molecule
- For telomere-negative molecules, estimate telomere position using a MEAN offset learned per contig
  from telomere-positive molecules: (telomere - anchor) in Qmap_position and siteID
- Determine telomeric molecule end using:
  (A) which contig end is telomeric (depends on chromosome arm + contig orientation)
  (B) molecule alignment Ori in the CSV (per molecule), which indicates forward/reverse mapping to contig

Outputs:
1) One summary per contig: 4 categories with count + percent among classified molecules.
2) Per molecule output: Molecule_ID with x and y + debugging columns and excluded reason if applicable.
"""

from __future__ import annotations
import math
from typing import Dict, List, Optional, Tuple
import pandas as pd

# ========== INPUTS & CONSTANTS ==========
INPUT_CSV = "Result_sgTelo_target_Mol_data.csv" # Output from update_extractingMolecules_MutilContig_V7_SingleSampleNew_tosent_fixed.py

TELOMERE_CHANNEL = 1 # Define Telomere Channel 
# Fusion Threshold: x distance (bp) between telomeric molecule end and telomere position (real/estimated)
THRESHOLD_BP = 10_000

# IMPORTANT: contig_orientation here is the CONTIG orientation (+/-) relative to chromosome p->q, not the per-molecule Ori column in the CSV
CONTIG_META = [
    {"contig_id": 1041, "arm": "q", "contig_orientation": "-", "target_site": 26},
]

OUT_SUMMARY_PER_CONTIG = "molecule_classification_summary_per_contig.csv"
OUT_PER_MOLECULE_XY = "molecule_classification_per_molecule.csv"
OUT_PER_MOLECULE_DEBUG = "molecule_classification_per_molecule_debug.csv"

# ========== INPUT VALIDATION ==========
REQUIRED_COLS = [
    "Molecule ID",
    "Contig_ID",
    "Ori",  # per-molecule alignment orientation (+/-) in the input CSV
    "Qmap_position",
    "LabelChannel",
    "siteID",
    "Contig_Site",
    "Contig_Position",
]

# ========== CATEGORY DEFINITIONS ==========
MAIN_CATEGORIES = [
    "Normal_With_Telomere",
    "Fused_With_Telomere",
    "Normal_Without_Telomere",
    "Fused_Without_Telomere",
]


# ========== LOAD & PREP DATA ==========
def load_csv(path: str) -> pd.DataFrame:
    """
    Purpose: Load input CSV and coerce types

    Notes:
    - We keep ALL rows, including Contig_Site==0, because they represent real molecule labels
    - We only drop rows missing essential identifiers (Molecule ID / Contig_ID / Qmap_position / siteID / Ori)
    """
    df = pd.read_csv(path)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Coerce numeric fields (becomes NaN if conversion fails)
    num_cols = [
        "Molecule ID",
        "Contig_ID",
        "Qmap_position",
        "LabelChannel",
        "siteID",
        "Contig_Site",
        "Contig_Position",
    ]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Drop rows missing the essentials
    df = df.dropna(subset=["Molecule ID", "Contig_ID", "Ori", "Qmap_position", "LabelChannel", "siteID"])

    # Make types consistent
    df["Molecule ID"] = df["Molecule ID"].astype(int)
    df["Contig_ID"] = df["Contig_ID"].astype(int)
    df["LabelChannel"] = df["LabelChannel"].astype(int)
    df["siteID"] = df["siteID"].astype(int)
    df["Ori"] = df["Ori"].astype(str).str.strip()

    return df


# ========== DEFINE TELOMERIC-END LOGIC ==========
def telomeric_contig_end(arm: str, contig_orientation: str) -> str:
    """
    Purpose: Determine which contig_site end corresponds to the telomere-end of a chromosome arm

    Returns:
      'low'  => telomere is at low contig_site numbers
      'high' => telomere is at high contig_site numbers

    Rules:
      If contig is '+' (stored p->q):
        p-arm telomere = low contig_site
        q-arm telomere = high contig_site

      If contig is '-' (stored q->p, reversed):
        p-arm telomere = high contig_site
        q-arm telomere = low contig_site
    """
    arm = arm.strip().lower()
    contig_orientation = contig_orientation.strip()

    if arm not in ("p", "q"):
        raise ValueError(f"arm must be 'p' or 'q', got {arm}")
    if contig_orientation not in ("+", "-"):
        raise ValueError(f"contig_orientation must be '+' or '-', got {contig_orientation}")

    if contig_orientation == "+":
        return "low" if arm == "p" else "high"
    else:
        return "high" if arm == "p" else "low"


def molecule_telomeric_end(
    telo_contig_end: str,
    molecule_ori: str,
    qmap_min: float,
    qmap_max: float,
    siteid_min: int,
    siteid_max: int,
) -> Tuple[float, int, str]:
    """
    Purpose: Determine which physical end of the MOLECULE is the telomeric end

    Inputs:
      - telo_contig_end: 'low' or 'high' (from telomeric_contig_end)
      - molecule_ori: value from input CSV column Ori ('+' or '-')
        IMPORTANT: This is per-molecule alignment orientation (NOT contig orientation)

    Interpretation of molecule Ori in this dataset:
      - Ori = '+' : contig sites increase as Qmap_position increases
      - Ori = '-' : contig sites decrease as Qmap_position increases

    Mapping table:
      telo_end=high & Ori='+' -> telomeric molecule end is qmap_max
      telo_end=high & Ori='-' -> telomeric molecule end is qmap_min
      telo_end=low  & Ori='+' -> telomeric molecule end is qmap_min
      telo_end=low  & Ori='-' -> telomeric molecule end is qmap_max

    Returns:
      (telomeric_end_qpos, telomeric_end_siteID, telomeric_end_side_str)
      side_str is 'qmap_min' or 'qmap_max' (useful debug column)
    """
    molecule_ori = str(molecule_ori).strip()
    if molecule_ori not in ("+", "-"):
        raise ValueError(f"Molecule Ori must be '+' or '-', got {molecule_ori}")
    if telo_contig_end not in ("low", "high"):
        raise ValueError(f"telo_contig_end must be 'low' or 'high', got {telo_contig_end}")

    if telo_contig_end == "high":
        if molecule_ori == "+":
            return qmap_max, siteid_max, "qmap_max"
        else:
            return qmap_min, siteid_min, "qmap_min"
    else:  # telo_contig_end == "low"
        if molecule_ori == "+":
            return qmap_min, siteid_min, "qmap_min"
        else:
            return qmap_max, siteid_max, "qmap_max"


# ========== MOLECULE-LEVEL PICKING ==========
def get_molecule_ends(mol_df: pd.DataFrame) -> Tuple[float, float, int, int]:
    """
    Purpose: Get physical molecule ends from ALL rows (including unmapped labels)

    Returns:
      qmap_min, qmap_max, siteID_at_qmap_min, siteID_at_qmap_max

    Note: If multiple rows share the same qmap_min or qmap_max, we choose min/max siteID among them for stable behavior
    """
    qmin = float(mol_df["Qmap_position"].min())
    qmax = float(mol_df["Qmap_position"].max())

    siteid_min = int(mol_df.loc[mol_df["Qmap_position"] == qmin, "siteID"].min())
    siteid_max = int(mol_df.loc[mol_df["Qmap_position"] == qmax, "siteID"].max())

    return qmin, qmax, siteid_min, siteid_max


def get_mapped_rows(mol_df: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose: Get mapped rows used for anchor selection
    A row is considered mapped if Contig_Site > 0 and Contig_Position > 0.

    Note: We do NOT delete unmapped rows globally; we only use this filtered view for anchor selection
    """
    m = mol_df.dropna(subset=["Contig_Site", "Contig_Position"]).copy()
    m = m[(m["Contig_Site"] > 0) & (m["Contig_Position"] > 0)]
    return m


def pick_anchor(mapped_df: pd.DataFrame, telo_contig_end: str) -> Optional[Tuple[int, float, int]]:
    """
    Purpose: Choose the anchor landmark on the molecule aka the mapped contig label nearest the telomeric contig end

    If telo_contig_end == 'low'  => anchor is the MIN mapped Contig_Site observed in this molecule.
    If telo_contig_end == 'high' => anchor is the MAX mapped Contig_Site observed in this molecule.

    Returns:
      (anchor_contig_site, anchor_qmap_position, anchor_siteID) or None if no mapped rows exist.

    Debug Note:
    Some molecules may have multiple rows mapping to the same Contig_Site
    -> We choose a stable representative: closest to the median Qmap_position among those rows
    """
    if mapped_df.empty:
        return None

    if telo_contig_end == "low":
        anchor_site = int(mapped_df["Contig_Site"].min())
    else:
        anchor_site = int(mapped_df["Contig_Site"].max())

    sub = mapped_df[mapped_df["Contig_Site"] == anchor_site].copy()
    if sub.empty:
        return None

    q_med = float(sub["Qmap_position"].median())
    sub["__dist__"] = (sub["Qmap_position"] - q_med).abs()
    best = sub.loc[sub["__dist__"].idxmin()]

    return anchor_site, float(best["Qmap_position"]), int(best["siteID"])


def pick_telomere_label_closest_to_end(
    mol_df: pd.DataFrame,
    telomeric_mol_end_qpos: float,
) -> Optional[Tuple[float, int]]:
    """
    Purpose: Choose the telomere label (LabelChannel=1) that is closest to the telomeric end of the molecule

    Returns:
      (tel_qmap_position, tel_siteID) or None if no telomere labels exist.
    """
    tel = mol_df[mol_df["LabelChannel"] == TELOMERE_CHANNEL].copy()
    if tel.empty:
        return None

    tel["__dist__"] = (tel["Qmap_position"] - telomeric_mol_end_qpos).abs()
    best = tel.loc[tel["__dist__"].idxmin()]
    return float(best["Qmap_position"]), int(best["siteID"])


# ========== MEAN OFFSETS PER CONTIG ==========
def safe_mean(values: List[float]) -> Optional[float]:
    vals = [
        v for v in values
        if v is not None and not (isinstance(v, float) and math.isnan(v))
    ]
    if not vals:
        return None
    return float(sum(vals) / len(vals))


def compute_mean_offsets_per_contig(per_molecule_rows: List[Dict]) -> Tuple[Optional[float], Optional[float], int]:
    """
    Purpose: Compute mean offsets (telomere - anchor) per contig using ONLY
      - molecules with a REAL telomere label (Telomere_Source == 'real')
      - molecules with an anchor (Anchor_Qmap_Position_bp not None)

    Offsets:
      mean_Δqpos = mean(tel_qpos - anchor_qpos)
      mean_Δsite = mean(tel_siteID - anchor_siteID)

    Returns:
      (mean_Δqpos, mean_Δsite, n_used)
    """
    dq: List[float] = []
    ds: List[float] = []

    for r in per_molecule_rows:
        if r.get("Telomere_Source") != "real":
            continue

        aq = r.get("Anchor_Qmap_Position_bp")
        asid = r.get("Anchor_SiteID")
        tq = r.get("Telomere_Qmap_Position_Used_bp")
        tsid = r.get("Telomere_SiteID_Used")

        if aq is None or asid is None or tq is None or tsid is None:
            continue

        dq.append(tq - aq)
        ds.append(float(tsid - asid))

    return safe_mean(dq), safe_mean(ds), len(dq)


def estimate_telomere_from_anchor(
    anchor_qpos: float,
    anchor_siteid: int,
    mean_dq: float,
    mean_ds: Optional[float],
) -> Tuple[float, int]:
    """
    Purpose: Estimate telomere position (qmap, siteID) when no telomere label exists

    tel_qpos_est = anchor_qpos + mean_Δqpos
    tel_site_est = anchor_siteid + mean_Δsite (rounded)
    """
    tel_qpos_est = float(anchor_qpos + mean_dq)
    if mean_ds is None:
        tel_site_est = int(anchor_siteid)
    else:
        tel_site_est = int(round(anchor_siteid + mean_ds))
    return tel_qpos_est, tel_site_est


# =========== CLASSIFICATION & CALCULATIONS ==========
def compute_x_y(
    mol_end_qpos_tel: float,
    mol_end_siteid_tel: int,
    tel_qpos_used: float,
    tel_siteid_used: int,
) -> Tuple[float, int]:
    """
    x: bp distance between telomeric molecule end and telomere position (real/estimated)
    y: label count difference (siteID distance)
    """
    x = float(abs(mol_end_qpos_tel - tel_qpos_used))
    y = int(abs(mol_end_siteid_tel - tel_siteid_used))
    return x, y


def category_from_x(telomere_present: bool, x_bp: float, threshold_bp: int) -> str:
    fused = x_bp > threshold_bp
    if telomere_present:
        return "Fused_With_Telomere" if fused else "Normal_With_Telomere"
    else:
        return "Fused_Without_Telomere" if fused else "Normal_Without_Telomere"


# ========== PER-CONTIG PROCESSING ==========
def process_contig(
    df: pd.DataFrame,
    contig_meta: Dict,
    threshold_bp: int,
) -> Tuple[List[Dict], Optional[float], Optional[float], int]:
    """
    Purpose: Process a single contig

    Returns:
      per_rows: list of per-molecule dicts (includes category or excluded reason)
      mean_dq, mean_ds, n_used: calibration stats per contig
    """
    contig_id = int(contig_meta["contig_id"])
    arm = str(contig_meta["arm"])
    contig_orientation = str(contig_meta["contig_orientation"]).strip()

    df_c = df[df["Contig_ID"] == contig_id].copy()
    if df_c.empty:
        return [], None, None, 0

    telo_end = telomeric_contig_end(arm, contig_orientation)

    # First Pass: compute molecule telomeric end, pick anchor, and if telomere exists pick it (real)
    per_rows: List[Dict] = []

    for mol_id, mol_df in df_c.groupby("Molecule ID", sort=False):
        mol_id = int(mol_id)

        # Ori appears constant per molecule in sample input; if not, use mode/most common
        ori = str(mol_df["Ori"].iloc[0]).strip()

        qmin, qmax, siteid_min, siteid_max = get_molecule_ends(mol_df)

        mol_end_qpos_tel, mol_end_siteid_tel, mol_tel_end_side = molecule_telomeric_end(
            telo_end, ori, qmin, qmax, siteid_min, siteid_max
        )

        tel_pick = pick_telomere_label_closest_to_end(mol_df, mol_end_qpos_tel)
        tel_present = tel_pick is not None

        mapped_df = get_mapped_rows(mol_df)
        anchor = pick_anchor(mapped_df, telo_end)

        # Base row with debug context (helps you explain & troubleshoot)
        row: Dict = {
            # IDs / metadata
            "Contig_ID": contig_id,
            "Chrom_Arm": arm,
            "Contig_Orientation": contig_orientation,
            "Telomeric_Contig_End": telo_end,
            "Molecule_ID": mol_id,
            "Molecule_Ori": ori,

            # Debug: physical molecule ends
            "Qmap_Min_bp": qmin,
            "Qmap_Max_bp": qmax,
            "SiteID_at_Qmap_Min": siteid_min,
            "SiteID_at_Qmap_Max": siteid_max,

            # Debug: which molecule end is telomeric for this molecule
            "Telomeric_Molecule_End_Side": mol_tel_end_side,  # 'qmap_min' or 'qmap_max'
            "Telomeric_Molecule_End_Qmap_bp": mol_end_qpos_tel,
            "Telomeric_Molecule_End_SiteID": mol_end_siteid_tel,

            # Anchor info (may be None)
            "Anchor_Contig_Site": anchor[0] if anchor else None,
            "Anchor_Qmap_Position_bp": anchor[1] if anchor else None,
            "Anchor_SiteID": anchor[2] if anchor else None,

            # Telomere info (real if present; estimated later if absent)
            "Telomere_Present": bool(tel_present),
            "Telomere_Source": "real" if tel_present else "estimated",
            "Telomere_Qmap_Position_Used_bp": tel_pick[0] if tel_present else None,
            "Telomere_SiteID_Used": tel_pick[1] if tel_present else None,

            # Outputs
            "x_bp": None,
            "y_labels": None,
            "Category": None,
            "Excluded_Reason": None,
        }

        # If telomere is present, we can compute x/y immediately (anchor not required)
        # (Anchor is only required for estimation when telomere is absent)
        if tel_present:
            x, y = compute_x_y(
                mol_end_qpos_tel, mol_end_siteid_tel,
                row["Telomere_Qmap_Position_Used_bp"], row["Telomere_SiteID_Used"]
            )
            row["x_bp"] = x
            row["y_labels"] = y
            row["Category"] = category_from_x(True, x, threshold_bp)
            row["Telomere_Source"] = "real"
        else:
            # Telomere absent: if we also have no anchor, we cannot estimate telomere position honestly
            if anchor is None:
                row["Telomere_Source"] = "none"
                row["Excluded_Reason"] = "no_mapped_anchor_for_estimation"

        per_rows.append(row)

    # Calibration: mean offsets per contig using telomere-positive molecules that have anchors
    mean_dq, mean_ds, n_used = compute_mean_offsets_per_contig(per_rows)

    # Second Pass: estimate telomere for telomere-negative molecules (that are not already excluded),then compute x/y and classify
    for row in per_rows:
        if row["Telomere_Source"] != "estimated":
            continue  # already real or excluded

        if mean_dq is None:
            row["Telomere_Source"] = "none"
            row["Excluded_Reason"] = "no_telomere_positive_molecules_for_mean"
            continue

        aq = row["Anchor_Qmap_Position_bp"]
        asid = row["Anchor_SiteID"]
        if aq is None or asid is None:
            row["Telomere_Source"] = "none"
            row["Excluded_Reason"] = "no_anchor_for_estimation"
            continue

        tel_q_est, tel_sid_est = estimate_telomere_from_anchor(aq, int(asid), mean_dq, mean_ds)

        row["Telomere_Qmap_Position_Used_bp"] = tel_q_est
        row["Telomere_SiteID_Used"] = tel_sid_est
        row["Telomere_Source"] = "estimated"

        x, y = compute_x_y(
            row["Telomeric_Molecule_End_Qmap_bp"], int(row["Telomeric_Molecule_End_SiteID"]),
            tel_q_est, tel_sid_est
        )
        row["x_bp"] = x
        row["y_labels"] = y
        row["Category"] = category_from_x(False, x, threshold_bp)

    return per_rows, mean_dq, mean_ds, n_used


# ========== OUTPUT BUILDERS ============
def summary_for_one_contig(per_df_contig: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose: Produce a 4-row summary for one contig: Contig_ID, Category, Count, Percent

    Percent is computed among classified molecules only (those in MAIN_CATEGORIES)
    Excluded molecules do not affect the denominator
    """
    contig_id = int(per_df_contig["Contig_ID"].iloc[0])

    df_main = per_df_contig[per_df_contig["Category"].isin(MAIN_CATEGORIES)].copy()
    total_classified = len(df_main)

    rows: List[Dict] = []
    for cat in MAIN_CATEGORIES:
        n = int((df_main["Category"] == cat).sum())
        pct = (100.0 * n / total_classified) if total_classified > 0 else 0.0
        rows.append({
            "Contig_ID": contig_id,
            "Category": cat,
            "Count": n,
            "Percent": round(pct, 2),
        })

    return pd.DataFrame(rows)


def make_summary_per_contig(per_molecule_df: pd.DataFrame) -> pd.DataFrame:
    summaries: List[pd.DataFrame] = []
    for contig_id, df_c in per_molecule_df.groupby("Contig_ID", sort=True):
        summaries.append(summary_for_one_contig(df_c))
    if summaries:
        return pd.concat(summaries, ignore_index=True)
    return pd.DataFrame(columns=["Contig_ID", "Category", "Count", "Percent"])


def make_per_molecule_xy_short(per_molecule_df: pd.DataFrame) -> pd.DataFrame:
    """
    SHORT VERSION per-molecule XY output: main outputs only (+ excluded reason)
    """
    cols = [
        "Contig_ID",
        "Molecule_ID",
        "Category",
        "x_bp",
        "y_labels",
        "Telomere_Source",
        "Excluded_Reason",
    ]
    return per_molecule_df[cols].copy()


def make_per_molecule_debug(per_molecule_df: pd.DataFrame) -> pd.DataFrame:
    """
    DEBUG VERSION per-molecule output: includes orientation/anchor/telomere/end details
    """
    # Keep everything (already structured), but you can reorder for readability:
    debug_cols = [
        "Contig_ID", "Molecule_ID",
        "Category", "x_bp", "y_labels", "Telomere_Source", "Excluded_Reason",
        "Chrom_Arm", "Contig_Orientation", "Telomeric_Contig_End", "Molecule_Ori",
        "Qmap_Min_bp", "Qmap_Max_bp",
        "Telomeric_Molecule_End_Side", "Telomeric_Molecule_End_Qmap_bp", "Telomeric_Molecule_End_SiteID",
        "Telomere_Present", "Telomere_Qmap_Position_Used_bp", "Telomere_SiteID_Used",
        "Anchor_Contig_Site", "Anchor_Qmap_Position_bp", "Anchor_SiteID",
    ]
    debug_cols = [c for c in debug_cols if c in per_molecule_df.columns]
    return per_molecule_df[debug_cols].copy()


# =========== MAIN ===========
def main() -> None:
    df = load_csv(INPUT_CSV)

    all_rows: List[Dict] = []

    for meta in CONTIG_META:
        rows, mean_dq, mean_ds, n_used = process_contig(df, meta, THRESHOLD_BP)
        all_rows.extend(rows)

        # Console logging helps debug quickly while developing
        cid = meta["contig_id"]
        print(f"\nContig {cid} calibration (mean offsets from telomere-positive molecules w/ anchors):")
        print(f"  mean_Δqpos (bp) = {mean_dq}")
        print(f"  mean_Δsite (labels) = {mean_ds}")
        print(f"  n_used = {n_used}")

    per_df = pd.DataFrame(all_rows)

    # Build outputs
    summary_df = make_summary_per_contig(per_df)
    per_xy_short_df = make_per_molecule_xy_short(per_df)
    per_debug_df = make_per_molecule_debug(per_df)

    # Write outputs
    summary_df.to_csv(OUT_SUMMARY_PER_CONTIG, index=False)
    per_xy_short_df.to_csv(OUT_PER_MOLECULE_XY, index=False)
    per_debug_df.to_csv(OUT_PER_MOLECULE_DEBUG, index=False)

    summary_df.to_csv(OUT_SUMMARY_PER_CONTIG, index=False)
    per_xy_short_df.to_csv(OUT_PER_MOLECULE_XY, index=False)
    per_debug_df.to_csv(OUT_PER_MOLECULE_DEBUG, index=False)

    # Print summary to console per contig
    print("\n=== Summary per contig (4 categories; percentages among classified molecules) ===")
    if not summary_df.empty:
        for contig_id, sdf in summary_df.groupby("Contig_ID", sort=True):
            print(f"\nContig {contig_id}")
            print(sdf.to_string(index=False))
    else:
        print("(No summary rows.)")


if __name__ == "__main__":
    main()
