import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import os
import csv
import re
import pandas as pd
from collections import defaultdict
import numpy as np


# ===== Analysis Functions (from main.py) =====

def selecting_molecules(df):
    unique_df_min = df.drop_duplicates(subset=["Molecule ID"], keep="first")
    unique_df_max = df.drop_duplicates(subset=["Molecule ID"], keep="last")
    molec_id = dict(zip(
        unique_df_min["Molecule ID"],
        zip(
            unique_df_min["Qmap_position"], unique_df_min["siteID"],
            unique_df_max["Qmap_position"], unique_df_max["siteID"],
            unique_df_min["Ori"]
        )
    ))
    return molec_id


def mk_categories():
    categories = ["fused_telomere", "not_fused_telomere_(normal)", "not_fused_no_telomere", "fused_no_telomere"]
    bins = defaultdict(set)
    for cat in categories:
        bins[cat]
    return bins


def finding_averages(df, contig, molec_id):
    general_dist_array = []
    contig_gap_dist = defaultdict(list)

    for molecule in molec_id.keys():
        mask = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == contig)
        try:
            idx = df.index[mask][0]
            current_pos = df.iloc[idx, 5]
            red_mask_before = df.iloc[idx - 1, 6]
            red_mask_after = df.iloc[idx + 1, 6]

            if red_mask_before == 1 or red_mask_after == 1:
                if red_mask_before == 1:
                    neighbor_pos = df.iloc[idx - 1, 5]
                    general_dist_array.append(current_pos - neighbor_pos)
                if red_mask_after == 1:
                    neighbor_pos = df.iloc[idx + 1, 5]
                    general_dist_array.append(current_pos - neighbor_pos)

            for i in range(0, 11):
                if contig + 5 - i == contig:
                    continue
                mask_i = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == contig + 5 - i)
                if not df[mask_i].empty:
                    gap_distance = abs(df[mask].iloc[0, 5] - df[mask_i].iloc[0, 5])
                    contig_gap_dist[contig + 5 - i].append(gap_distance)
        except IndexError:
            continue

    general_dist_avg = float(np.average(general_dist_array)) if general_dist_array else 0.0
    contig_gap_dist_avg = {k: float(np.average(v)) for k, v in contig_gap_dist.items()}
    return general_dist_avg, contig_gap_dist_avg


def classify_molecules(df, bins, molec_id, contig, label_avg, gap_avg, chrom_orientation):
    for molecule in molec_id.keys():
        mask = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == contig)

        try:
            idx = df.index[mask][0]
        except IndexError:
            # Try nearby contigs ±5
            check = True
            contig_cp_lst = [contig - 1, contig + 1, contig - 2, contig + 2,
                             contig - 3, contig + 3, contig - 4, contig + 4,
                             contig - 5, contig + 5]
            contig_cp_idx = -1
            while check:
                if contig_cp_idx >= len(contig_cp_lst) - 1:
                    break
                contig_cp_idx += 1
                contig_cp = contig_cp_lst[contig_cp_idx]
                mask = (df["Molecule ID"] == molecule) & (df["Contig_Site"] == contig_cp)
                try:
                    idx = df.index[mask][0]
                    red_mask_before = df.iloc[idx - 1, 6]
                    red_mask_after = df.iloc[idx + 1, 6]
                    check = False

                    if red_mask_before == 1 or red_mask_after == 1:
                        if red_mask_before == 1:
                            qmap_pos = df.iloc[idx - 1, 5]
                            row_pos = df.iloc[idx - 1, 7]
                        else:
                            qmap_pos = df.iloc[idx + 1, 5]
                            row_pos = df.iloc[idx + 1, 7]
                        distance_kb, row_distance = _calc_distances(
                            molec_id[molecule], qmap_pos, row_pos, chrom_orientation)
                        _add_to_bin(bins, molecule, distance_kb, row_distance, has_label=True)
                    else:
                        qmap_pos = df[mask].iloc[0, 5] + label_avg + gap_avg.get(contig_cp, 0)
                        row_pos = df[mask].iloc[0, 7]
                        distance_kb, row_distance = _calc_distances(
                            molec_id[molecule], qmap_pos, row_pos, chrom_orientation)
                        _add_to_bin(bins, molecule, distance_kb, row_distance, has_label=False)
                except IndexError:
                    continue
            continue

        red_mask_before = df.iloc[idx - 1, 6]
        red_mask_after = df.iloc[idx + 1, 6]

        if red_mask_before == 1 or red_mask_after == 1:
            if red_mask_before == 1:
                qmap_pos = df.iloc[idx - 1, 5]
                row_pos = df.iloc[idx - 1, 7]
            else:
                qmap_pos = df.iloc[idx + 1, 5]
                row_pos = df.iloc[idx + 1, 7]
            distance_kb, row_distance = _calc_distances(
                molec_id[molecule], qmap_pos, row_pos, chrom_orientation)
            _add_to_bin(bins, molecule, distance_kb, row_distance, has_label=True)
        else:
            qmap_pos = df[mask].iloc[0, 5] + label_avg
            row_pos = df[mask].iloc[0, 7]
            distance_kb, row_distance = _calc_distances(
                molec_id[molecule], qmap_pos, row_pos, chrom_orientation)
            _add_to_bin(bins, molecule, distance_kb, row_distance, has_label=False)

    return bins


def _calc_distances(mol_data, qmap_pos, row_pos, chrom_orientation):
    ori = mol_data[-1]
    if "p" in chrom_orientation and ori == "+":
        distance_kb = abs(mol_data[0] - qmap_pos)
        row_distance = abs(mol_data[1] - row_pos)
    elif "p" in chrom_orientation and ori == "-":
        distance_kb = abs(mol_data[2] - qmap_pos)
        row_distance = abs(mol_data[3] - row_pos)
    elif "q" in chrom_orientation and ori == "+":
        distance_kb = abs(mol_data[2] - qmap_pos)
        row_distance = abs(mol_data[3] - row_pos)
    elif "q" in chrom_orientation and ori == "-":
        distance_kb = abs(mol_data[0] - qmap_pos)
        row_distance = abs(mol_data[1] - row_pos)
    else:
        raise ValueError(f"Unexpected orientation combination: {chrom_orientation}, {ori}")
    return distance_kb, row_distance


def _add_to_bin(bins, molecule, distance_kb, row_distance, has_label):
    entry = f"{molecule}_({distance_kb},{row_distance})"
    if has_label:
        if distance_kb >= 10_000:
            bins["fused_telomere"].add(entry)
        else:
            bins["not_fused_telomere_(normal)"].add(entry)
    else:
        if distance_kb >= 10_000:
            bins["fused_no_telomere"].add(entry)
        else:
            bins["not_fused_no_telomere"].add(entry)


def parse_bin_entry(entry):
    """Parse 'moleculeID_(distance_kb, row_distance)' -> (mol_id, distance_kb, row_distance)."""
    match = re.match(r'^(.+)_\(([^,]+),\s*([^)]+)\)$', str(entry))
    if match:
        return match.group(1), float(match.group(2)), float(match.group(3))
    return str(entry), 0.0, 0.0


# ===== Category display names =====
CATEGORIES = [
    ("fused_telomere",             "Fused Telomere"),
    ("not_fused_telomere_(normal)", "Not Fused Telomere (Normal)"),
    ("not_fused_no_telomere",      "Not Fused No Telomere"),
    ("fused_no_telomere",          "Fused No Telomere"),
]


# ===== GUI =====

class TelomereApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Automated Detection of Telomeric Structural Variations")
        self.root.geometry("1400x900")

        self.excel_path = tk.StringVar()
        self.chromosome_var = tk.StringVar()
        self.orientation_var = tk.StringVar(value="p")
        self.contig_var = tk.StringVar()

        self.bins = None
        self.label_avg = None

        self._create_widgets()

    def _create_widgets(self):
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=3)
        main_frame.rowconfigure(1, weight=1)

        # Title
        ttk.Label(
            main_frame,
            text="Automated Detection of Telomeric Structural Variations",
            font=('Arial', 16, 'bold')
        ).grid(row=0, column=0, columnspan=2, pady=(0, 12))

        # ── Left panel ──────────────────────────────────────────────────────
        left = ttk.Frame(main_frame)
        left.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(0, 12))
        left.columnconfigure(0, weight=1)

        # Inputs
        inp = ttk.LabelFrame(left, text="Analysis Inputs", padding="10")
        inp.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        inp.columnconfigure(1, weight=1)

        # Excel file row
        ttk.Label(inp, text="Excel File:").grid(row=0, column=0, sticky=tk.W, pady=4)
        file_row = ttk.Frame(inp)
        file_row.grid(row=0, column=1, sticky=(tk.W, tk.E), pady=4)
        file_row.columnconfigure(0, weight=1)
        ttk.Entry(file_row, textvariable=self.excel_path).grid(
            row=0, column=0, sticky=(tk.W, tk.E), padx=(0, 5))
        ttk.Button(file_row, text="Browse", command=self._browse_excel, width=8).grid(row=0, column=1)

        # Chromosome
        ttk.Label(inp, text="Chromosome:").grid(row=1, column=0, sticky=tk.W, pady=4)
        ttk.Entry(inp, textvariable=self.chromosome_var, width=14).grid(
            row=1, column=1, sticky=tk.W, pady=4)

        # Orientation
        ttk.Label(inp, text="Orientation:").grid(row=2, column=0, sticky=tk.W, pady=4)
        ori_frame = ttk.Frame(inp)
        ori_frame.grid(row=2, column=1, sticky=tk.W, pady=4)
        ttk.Radiobutton(ori_frame, text="p arm", variable=self.orientation_var, value="p").grid(
            row=0, column=0, padx=(0, 12))
        ttk.Radiobutton(ori_frame, text="q arm", variable=self.orientation_var, value="q").grid(
            row=0, column=1)

        # Contig
        ttk.Label(inp, text="Contig Site:").grid(row=3, column=0, sticky=tk.W, pady=4)
        ttk.Entry(inp, textvariable=self.contig_var, width=14).grid(
            row=3, column=1, sticky=tk.W, pady=4)

        ttk.Button(inp, text="Run Analysis", command=self._run_analysis).grid(
            row=4, column=0, columnspan=2, pady=(12, 0), ipadx=20, ipady=4)

        # Summary
        summ = ttk.LabelFrame(left, text="Summary", padding="10")
        summ.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        summ.columnconfigure(1, weight=1)

        self._summ_labels = {}
        rows = [
            ("total",                        "Total Molecules"),
            ("pct_fusion",                   "Total % Fusion"),
            ("sep", None),
            ("fused_telomere",               "Fused Telomere"),
            ("not_fused_telomere_(normal)",  "Not Fused Telomere (Normal)"),
            ("not_fused_no_telomere",        "Not Fused No Telomere"),
            ("fused_no_telomere",            "Fused No Telomere"),
            ("sep2", None),
            ("label_avg",                    "Label Distance Avg (kb)"),
        ]
        r = 0
        for key, label in rows:
            if label is None:
                ttk.Separator(summ, orient=tk.HORIZONTAL).grid(
                    row=r, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=4)
                r += 1
                continue
            ttk.Label(summ, text=label + ":", font=('Arial', 9, 'bold')).grid(
                row=r, column=0, sticky=tk.W, pady=2)
            val = ttk.Label(summ, text="—", font=('Arial', 9))
            val.grid(row=r, column=1, sticky=tk.W, padx=(10, 0), pady=2)
            self._summ_labels[key] = val
            r += 1

        # Export button
        ttk.Button(left, text="Export Full Data to CSV", command=self._export_csv).grid(
            row=2, column=0, pady=8, ipadx=10, ipady=4)

        # ── Right panel: tabbed results ──────────────────────────────────────
        right = ttk.Frame(main_frame)
        right.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N, tk.S))
        right.columnconfigure(0, weight=1)
        right.rowconfigure(0, weight=1)

        results_box = ttk.LabelFrame(right, text="Classification Results", padding="10")
        results_box.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        results_box.columnconfigure(0, weight=1)
        results_box.rowconfigure(0, weight=1)

        self.notebook = ttk.Notebook(results_box)
        self.notebook.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        self._trees = {}
        self._count_labels = {}
        for cat_key, cat_display in CATEGORIES:
            tab = ttk.Frame(self.notebook, padding="5")
            self.notebook.add(tab, text=cat_display)
            tab.columnconfigure(0, weight=1)
            tab.rowconfigure(1, weight=1)

            count_lbl = ttk.Label(tab, text="Count: —", font=('Arial', 9, 'bold'))
            count_lbl.grid(row=0, column=0, sticky=tk.W, pady=(0, 4))
            self._count_labels[cat_key] = count_lbl

            tree_frame = ttk.Frame(tab)
            tree_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
            tree_frame.columnconfigure(0, weight=1)
            tree_frame.rowconfigure(0, weight=1)

            tree = ttk.Treeview(
                tree_frame,
                columns=("molecule", "distance_kb", "row_distance"),
                show="headings",
                height=22,
            )
            tree.heading("molecule",      text="Molecule ID")
            tree.heading("distance_kb",   text="Distance (kb)")
            tree.heading("row_distance",  text="Row Distance")
            tree.column("molecule",      width=220, anchor=tk.W)
            tree.column("distance_kb",   width=150, anchor=tk.E)
            tree.column("row_distance",  width=150, anchor=tk.E)
            tree.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

            vsb = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL, command=tree.yview)
            vsb.grid(row=0, column=1, sticky=(tk.N, tk.S))
            tree.configure(yscrollcommand=vsb.set)

            hsb = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL, command=tree.xview)
            hsb.grid(row=1, column=0, sticky=(tk.W, tk.E))
            tree.configure(xscrollcommand=hsb.set)

            self._trees[cat_key] = tree

        # Status bar
        self._status_var = tk.StringVar(value="Ready — select an Excel file and enter parameters")
        ttk.Label(self.root, textvariable=self._status_var, relief=tk.SUNKEN, anchor=tk.W).grid(
            row=1, column=0, sticky=(tk.W, tk.E), padx=2, pady=(2, 0))

    # ── Event handlers ───────────────────────────────────────────────────────

    def _browse_excel(self):
        path = filedialog.askopenfilename(
            title="Select Excel file",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        if path:
            self.excel_path.set(path)
            self._set_status(f"File selected: {os.path.basename(path)}")

    def _run_analysis(self):
        excel_file  = self.excel_path.get().strip()
        chromosome  = self.chromosome_var.get().strip()
        orientation = self.orientation_var.get().strip()
        contig_str  = self.contig_var.get().strip()

        if not excel_file:
            messagebox.showwarning("Missing Input", "Please select an Excel file.")
            return
        if not chromosome:
            messagebox.showwarning("Missing Input", "Please enter a chromosome (e.g. 4).")
            return
        if not contig_str:
            messagebox.showwarning("Missing Input", "Please enter a contig site number.")
            return
        try:
            contig = int(contig_str)
        except ValueError:
            messagebox.showerror("Invalid Input", "Contig site must be an integer.")
            return

        sheet_name = chromosome + orientation
        self._set_status(f"Running analysis — sheet '{sheet_name}', contig {contig}…")
        self.root.update_idletasks()

        try:
            df = pd.read_excel(excel_file, sheet_name=sheet_name)
            molec_id    = selecting_molecules(df)
            empty_bins  = mk_categories()
            label_avg, gap_avg = finding_averages(df, contig, molec_id)
            bins = classify_molecules(
                df, empty_bins, molec_id, contig,
                label_avg, gap_avg, orientation + chromosome
            )
            self.bins      = bins
            self.label_avg = label_avg
            self._populate_results(bins, label_avg)
            self._set_status("Analysis complete.")
        except Exception as exc:
            messagebox.showerror("Error", f"Analysis failed:\n{exc}")
            self._set_status("Analysis failed — see error dialog.")

    def _populate_results(self, bins, label_avg):
        total = sum(len(v) for v in bins.values())

        for cat_key, _ in CATEGORIES:
            tree    = self._trees[cat_key]
            entries = bins.get(cat_key, set())

            tree.delete(*tree.get_children())
            for entry in sorted(entries):
                mol_id, dist_kb, row_dist = parse_bin_entry(entry)
                tree.insert("", tk.END, values=(
                    mol_id,
                    f"{dist_kb:,.2f}",
                    f"{row_dist:,.2f}",
                ))

            count = len(entries)
            pct   = round(100 * count / total, 3) if total > 0 else 0.0
            self._count_labels[cat_key].config(text=f"Count: {count}  ({pct}%)")

        # Summary panel
        self._summ_labels["total"].config(text=str(total))

        fused    = len(bins["fused_telomere"]) + len(bins["fused_no_telomere"])
        pct_fuse = round(100 * fused / total, 3) if total > 0 else 0.0
        self._summ_labels["pct_fusion"].config(text=f"{pct_fuse}%")

        for cat_key, _ in CATEGORIES:
            count = len(bins.get(cat_key, set()))
            pct   = round(100 * count / total, 3) if total > 0 else 0.0
            self._summ_labels[cat_key].config(text=f"{count}  ({pct}%)")

        self._summ_labels["label_avg"].config(text=f"{label_avg:,.2f}")

    def _export_csv(self):
        if self.bins is None:
            messagebox.showwarning("No Data", "Run the analysis first before exporting.")
            return

        filename = filedialog.asksaveasfilename(
            title="Save full data as CSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not filename:
            return

        try:
            with open(filename, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(["Molecule ID", "Category", "Distance (kb)", "Row Distance"])
                for cat_key, cat_display in CATEGORIES:
                    for entry in sorted(self.bins.get(cat_key, set())):
                        mol_id, dist_kb, row_dist = parse_bin_entry(entry)
                        writer.writerow([mol_id, cat_display, f"{dist_kb:.4f}", f"{row_dist:.4f}"])

            messagebox.showinfo("Exported", f"Data saved to:\n{filename}")
            self._set_status(f"Exported: {os.path.basename(filename)}")
        except Exception as exc:
            messagebox.showerror("Export Error", f"Could not save file:\n{exc}")

    def _set_status(self, message):
        self._status_var.set(message)


def main():
    root = tk.Tk()
    TelomereApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
