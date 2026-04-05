import pandas as pd
import os

# ── Config ────────────────────────────────────────────────────────────────────
INPUT_FILE  = "../data/full_filtered_q_sv_support_molecules_discussion.xlsx"
OUTPUT_FILE = "../analysis/2nd_milestone_output.txt"
CONTINUITY_THRESHOLD_BP = 1000   # max gap between consecutive segments to be "continuous"
SHEET_NAME  = "raw_data"          # sheet that holds the alignment records
# ─────────────────────────────────────────────────────────────────────────────


def load_data(path: str) -> pd.DataFrame:
    """Load the xlsx and return only the columns we need."""
    df = pd.read_excel(path, sheet_name=SHEET_NAME)
    required = {"QueryID", "RefID", "RefStartCoord", "RefStopCoord"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing expected columns: {missing}")
    return df


def select_molecules(df: pd.DataFrame) -> list:
    """
    Ask the user whether to use automatic selection (QueryIDs that appear
    more than once) or to enter IDs manually.
    Returns a list of QueryID values to analyse.
    """
    print("\n=== Molecule Selection ===")
    print("  [1] Automatic – keep only QueryIDs that appear more than once")
    print("  [2] Manual   – enter the QueryIDs you want to analyse")
    choice = input("Choose an option (1 or 2): ").strip()

    if choice == "1":
        counts = df["QueryID"].value_counts()
        selected = counts[counts > 1].index.tolist()
        print(f"\n  → {len(selected)} molecules selected automatically.")
        return selected

    elif choice == "2":
        raw = input("Enter QueryIDs separated by spaces: ").strip()
        entered = [int(x) for x in raw.split() if x.strip()]
        # Validate that they exist in the file
        valid   = [qid for qid in entered if qid in df["QueryID"].values]
        invalid = [qid for qid in entered if qid not in df["QueryID"].values]
        if invalid:
            print(f"  ⚠ The following IDs were not found and will be skipped: {invalid}")
        print(f"\n  → {len(valid)} molecules selected.")
        return valid

    else:
        print("  Invalid choice. Defaulting to automatic selection.")
        counts = df["QueryID"].value_counts()
        return counts[counts > 1].index.tolist()


def classify_molecule(rows: pd.DataFrame) -> str:
    """
    Given all alignment rows for a single QueryID, decide:
      - 'Not Fused'      – all segments are continuous (gap ≤ threshold)
      - 'Fused (Intra)'  – not continuous, but all rows share the same RefID
      - 'Fused (Inter)'  – not continuous AND at least one RefID differs
    """
    # Sort segments by their start coordinate on the reference
    rows = rows.sort_values("RefStartCoord").reset_index(drop=True)

    # Check continuity between consecutive segments
    is_continuous = True
    for i in range(1, len(rows)):
        gap = rows.loc[i, "RefStartCoord"] - rows.loc[i - 1, "RefStopCoord"]
        if gap > CONTINUITY_THRESHOLD_BP:
            is_continuous = False
            break

    if is_continuous:
        return "Not Fused"

    # It is fused – determine intra vs inter
    if rows["RefID"].nunique() == 1:
        return "Fused (Intra)"
    else:
        return "Fused (Inter)"


def write_output(results: dict, path: str):
    """
    Write a bare-bones tab-separated file:
      header row: Fused (Inter)\tFused (Intra)\tNot Fused
      one molecule ID per row, empty string if a column has fewer entries.
    """
    inter = sorted(results.get("Fused (Inter)", []))
    intra = sorted(results.get("Fused (Intra)", []))
    not_f = sorted(results.get("Not Fused",    []))

    # Pad shorter lists with empty strings so rows line up
    max_len = max(len(inter), len(intra), len(not_f), 1)
    inter += [""] * (max_len - len(inter))
    intra += [""] * (max_len - len(intra))
    not_f += [""] * (max_len - len(not_f))

    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "w") as fh:
        fh.write("Fused (Inter)\tFused (Intra)\tNot Fused\n")
        for a, b, c in zip(inter, intra, not_f):
            fh.write(f"{a}\t{b}\t{c}\n")

    print(f"\n  Results written to: {path}")


def main():
    print("Loading data …")
    df = load_data(INPUT_FILE)
    print(f"  {len(df)} rows loaded, {df['QueryID'].nunique()} unique QueryIDs.")

    selected_ids = select_molecules(df)
    if not selected_ids:
        print("No molecules to analyse. Exiting.")
        return

    # Filter to selected molecules only
    df_sel = df[df["QueryID"].isin(selected_ids)]

    # Classify each molecule
    print("\nClassifying molecules …")
    results = {"Fused (Inter)": [], "Fused (Intra)": [], "Not Fused": []}

    for qid, group in df_sel.groupby("QueryID"):
        category = classify_molecule(group)
        results[category].append(qid)

    # Print summary to console
    print("\n=== Classification Summary ===")
    for cat, ids in results.items():
        print(f"  {cat:<20}: {len(ids)} molecules")

    write_output(results, OUTPUT_FILE)


if __name__ == "__main__":
    main()
