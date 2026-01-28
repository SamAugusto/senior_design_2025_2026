# Anne Nguyen (an875)
# Last Update: 01/28/2026 
# Script for classifying molecules from update_extractingMolecules_MutilContig_V7_SingleSampleNew_tosent_fixed.py's OUTPUT .CSV FILE into 4 main categories

import pandas as pd
from typing import Dict, List, Tuple, Optional

# ----- INPUTS & CONSTANTS -----
INPUT_CSV = "Result_sgTelo_target_Mol_data.csv"

# Can populate this list from the Chromsome Sample text file later
# For now, we hardcode example that matches the INPUT_CSV file we are looking at
# - arm: "p" or "q" chromosome arm
# - target_site: landmark/target contig site ID that was manually picked using Bionano for a specific contig
# I believe in the text file, the corresponding columns of interest are: "contig", "p or q", and "contig_site ID"
TARGETS = [
    {"contig_id": 1041, "arm": "q", "target_site": 26},  # Example
]

# Assumption (for now): contig assembled in p -> q direction
ASSUME_CONTIG_P_TO_Q = True

# Define constants
TELOMERE_CHANNEL = 1
THRESHOLD_BP = 10_000  # 10 kb
OUTPUT_PREFIX = "classification_outputs"


def telomeric_side_from_arm(arm: str) -> str:
    """
    Under the assumption contig is oriented p -> q:
      • p-arm telomere is on the LOW contig-site end
      • q-arm telomere is on the HIGH contig-site end
    """
    arm = arm.lower().strip()
    if arm not in {"p", "q"}:
        raise ValueError(f"arm must be 'p' or 'q', got: {arm}")
    if not ASSUME_CONTIG_P_TO_Q:
        # If we later stop assuming p->q, we should set telomeric_side manually
        raise ValueError("ASSUME_CONTIG_P_TO_Q is False; telomeric side must be determined another way.")
    return "low" if arm == "p" else "high"


def load_and_prepare_csv(input_csv: str) -> pd.DataFrame:
    """
    - Load CSV and coerce key columns to numeric
    - Filter to include all must-have columns for analysis
    """
    df = pd.read_csv(input_csv)

    for col in ["Molecule ID", "Contig_ID", "LabelChannel", "Qmap_position", "Contig_Site", "Contig_Position"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Must-have fields for grouping and distance calculations
    df = df.dropna(subset=["Molecule ID", "Contig_ID", "LabelChannel", "Qmap_position"])
    df["Molecule ID"] = df["Molecule ID"].astype(int)
    df["Contig_ID"] = df["Contig_ID"].astype(int)
    return df


def mapped_rows(mol_df: pd.DataFrame) -> pd.DataFrame:
    """
    - Filter rows with meaningful contig mapping ONLY to use for landmark coverage + choosing contig-site end
    - Does NOT include rows with Contig_Position & Contig_Site = 0 since that means molecule with those sites did not map to contig
    Used for:
      • confirming target_site coverage
      • computing landmark_qpos
      • computing alignment-end qpos (min/max Contig_Site depending on telomeric side)
    """
    m = mol_df.dropna(subset=["Contig_Position", "Contig_Site"]).copy()
    m = m[(m["Contig_Position"] > 0) & (m["Contig_Site"] > 0)]
    return m


# ----- Main Classification Logic -----
def classify_one_molecule(mol_df: pd.DataFrame, # All rows in one specific molecule (one Molecule ID)
                          target_site: int,     # Chosen landmark/target contig site
                          telomeric_side: str,  # Either "high" or "low" Contig Sites: (low end) ------contig------ (high end)
                          telomere_channel: int,# Usually LabelChannel = 1 
                          threshold_bp: int     # Threshold (bp) to determine if molecule is normal or fused
                          ) -> Tuple[str,             # category
                                     Optional[float], # distance_bp
                                     bool,            # telomere_present
                                     Optional[float], # landmark_qpos
                                     Optional[float], # end_qpos
                                     Optional[float], # telomere_qpos_used (only if telomere_present)
                                     ]:    
    """
    Returns:
      category (str): 4 main categories or Unclassified
      distance_bp (float | None): the distance used for threshold comparison
      telomere_present (bool): whether telomere label was present
      landmark_qpos (None if Unclassified)
      end_qpos (None if Unclassified)
      telomere_qpos_used (None if no telomere labels present)
    """
    m = mapped_rows(mol_df)

    # Determine if molecule cover the target/landmark contig site
    if m.empty or not (m["Contig_Site"] == target_site).any():
        # Telomere present is still detectable, but molecule is unclassifiable for this target
        tel_present = (mol_df["LabelChannel"] == telomere_channel).any()
        return "Unclassified", None, bool(tel_present), None, None, None

    # Find the position of landmark/target contig site on molecule (Qmap_position at target site)
    # This gives where the landmark site lies on this molecule in molecule coordinates
    landmark_qpos = float(m.loc[m["Contig_Site"] == target_site, "Qmap_position"].median())

    # Determine alignment-end on telomeric side using contig-site end
    # This gives where the molecule’s alignment reaches the telomeric end of the contig, in molecule coordinates
    telomeric_side = telomeric_side.lower().strip()
    if telomeric_side == "high":
        end_site = m["Contig_Site"].max()
    elif telomeric_side == "low":
        end_site = m["Contig_Site"].min()
    else:
        raise ValueError("telomeric_side must be 'high' or 'low'")

    end_qpos = float(m.loc[m["Contig_Site"] == end_site, "Qmap_position"].median())

    # Determine if telomere label exists (may be unmapped to contig, but have Qmap_position)
    tel_rows = mol_df[mol_df["LabelChannel"] == telomere_channel]
    tel_present = len(tel_rows) > 0

    # If telomere label EXISTS
    if tel_present:
        # Compute distance from alignment end to closest telomere label (Qmap coords)
        deltas = (tel_rows["Qmap_position"] - end_qpos).abs()
        min_idx = deltas.idxmin()
        tel_qpos_used = float(tel_rows.loc[min_idx, "Qmap_position"])
        distance_bp = float(deltas.loc[min_idx])
        fused = distance_bp > threshold_bp
        category = "Fused_With_Telomere" if fused else "Normal_With_Telomere"
        return category, distance_bp, True, landmark_qpos, end_qpos, tel_qpos_used
    # If telomere label DOES NOT EXIST
    else:
        # Compute distance from alignment end to landmark/target contig site (Qmap coords)
        tel_qpos_used = None
        distance_bp = float(abs(end_qpos - landmark_qpos))
        fused = distance_bp > threshold_bp
        category = "Fused_Without_Telomere" if fused else "Normal_Without_Telomere"
        return category, distance_bp, False, landmark_qpos, end_qpos, tel_qpos_used
    

# ----- Main Output Files Formatting -----
def classify_contig_with_reports(df: pd.DataFrame,
                                 contig_id: int,
                                 arm: str,
                                 target_site: int,
                                 telomere_channel: int,
                                 threshold_bp: int,
                                 output_prefix: str
                                 ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    For a single contig:
      - classify each molecule
      - produce per-molecule report dataframe
      - produce summary dataframe
      - write both dataframe to CSV

    Returns:
      per_molecule_df, summary_df
    """
    telomeric_side = telomeric_side_from_arm(arm) # Determine telomeric side from arm information already provided
    df_c = df[df["Contig_ID"] == int(contig_id)].copy() # Filter to select one contig at a time for analysis

    # Per-molecule report rows
    report_rows: List[Dict] = []

    # Loop over molecules in this contig
    for mol_id, mol_df in df_c.groupby("Molecule ID", sort=False):
        mol_id_int = int(mol_id)
        category, dist_bp, tel_present, landmark_qpos, end_qpos, tel_qpos_used = classify_one_molecule(
            mol_df=mol_df,
            target_site=target_site,
            telomeric_side=telomeric_side,
            telomere_channel=telomere_channel,
            threshold_bp=threshold_bp,
        )

        report_rows.append(
            {
                "Contig_ID": contig_id,
                "Chrom_Arm": arm,
                "Telomeric_Side_Assumed": telomeric_side,
                "Target_Contig_Site": target_site,
                "Molecule_ID": mol_id_int,
                "Telomere_Present": tel_present,
                "Category": category,
                # DEBUG: key molecule-space positions
                "Target_Site_Qmap_Position_bp": landmark_qpos,
                "Alignment_End_Qmap_Position_bp": end_qpos,
                "Telomere_Qmap_Position_Used_bp": tel_qpos_used,  # None if no telomere
                "Threshold_bp": threshold_bp,
                "Calculated_Distance_bp_For_Threshold_Comparison": dist_bp
            }
        )

    per_molecule_df = pd.DataFrame(report_rows)

    # Summary counts + percents (exclude Unclassified from denominator)
    main_cats = [
        "Fused_With_Telomere",
        "Normal_With_Telomere",
        "Fused_Without_Telomere",
        "Normal_Without_Telomere",
    ]
    counts = per_molecule_df["Category"].value_counts()
    total_main = int(sum(counts.get(c, 0) for c in main_cats))

    summary_rows = []
    for c in main_cats:
        n = int(counts.get(c, 0))
        pct = (100.0 * n / total_main) if total_main > 0 else 0.0
        summary_rows.append({"Contig_ID": contig_id, "Category": c, "Count": n, "Percent": round(pct, 2)})

    # Add unclassified count
    summary_rows.append(
        {
            "Contig_ID": contig_id,
            "Category": "Unclassified",
            "Count": int(counts.get("Unclassified", 0)),
            "Percent": "N/A",
        }
    )

    summary_df = pd.DataFrame(summary_rows)

    # Write outputs
    per_molecule_path = f"{output_prefix}_contig{contig_id}_per_molecule.csv"
    summary_path = f"{output_prefix}_contig{contig_id}_summary.csv"

    per_molecule_df.to_csv(per_molecule_path, index=False)
    summary_df.to_csv(summary_path, index=False)

    print(f"Saved per-molecule report: {per_molecule_path}")
    print(f"Saved summary:             {summary_path}")

    return per_molecule_df, summary_df


def main() -> None:
    df = load_and_prepare_csv(INPUT_CSV)

    # Track which contigs exist in the data
    contigs_in_data = set(df["Contig_ID"].unique().tolist())
    targets_contigs = set(int(t["contig_id"]) for t in TARGETS)

    # Warn if data contains contigs we aren't analyzing
    extra = sorted(list(contigs_in_data - targets_contigs))
    if extra:
        print(f"WARNING: input CSV contains contigs not listed in TARGETS (skipping): {extra}")

    # Process each target contig and write outputs
    all_per_molecule = []
    all_summary = []

    for t in TARGETS:
        contig_id = int(t["contig_id"])
        arm = str(t["arm"])
        target_site = int(t["target_site"])

        if contig_id not in contigs_in_data:
            print(f"WARNING: contig {contig_id} not found in input CSV. Skipping.")
            continue

        print(f"\n===== Processing contig {contig_id} (arm={arm}, target_site={target_site}) =====")
        per_mol_df, summary_df = classify_contig_with_reports(
            df=df,
            contig_id=contig_id,
            arm=arm,
            target_site=target_site,
            telomere_channel=TELOMERE_CHANNEL,
            threshold_bp=THRESHOLD_BP,
            output_prefix=OUTPUT_PREFIX,
        )

        print(summary_df.to_string(index=False))

        all_per_molecule.append(per_mol_df)
        all_summary.append(summary_df)


if __name__ == "__main__":
    main()