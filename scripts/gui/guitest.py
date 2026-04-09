import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import os
import csv
import re
import pandas as pd

from main import selecting_molecules, mk_categories, finding_averages, classify_molecules


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
        self.sheet_var = tk.StringVar()
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

        # Sheet picker
        ttk.Label(inp, text="Sheet:").grid(row=1, column=0, sticky=tk.W, pady=4)
        self._sheet_combo = ttk.Combobox(inp, textvariable=self.sheet_var,
                                         state="readonly", width=20)
        self._sheet_combo.grid(row=1, column=1, sticky=tk.W, pady=4)
        self._sheet_combo.bind("<<ComboboxSelected>>", self._on_sheet_selected)

        # Chromosome
        ttk.Label(inp, text="Chromosome:").grid(row=2, column=0, sticky=tk.W, pady=4)
        ttk.Entry(inp, textvariable=self.chromosome_var, width=14).grid(
            row=2, column=1, sticky=tk.W, pady=4)

        # Orientation
        ttk.Label(inp, text="Orientation:").grid(row=3, column=0, sticky=tk.W, pady=4)
        ori_frame = ttk.Frame(inp)
        ori_frame.grid(row=3, column=1, sticky=tk.W, pady=4)
        ttk.Radiobutton(ori_frame, text="p arm", variable=self.orientation_var, value="p").grid(
            row=0, column=0, padx=(0, 12))
        ttk.Radiobutton(ori_frame, text="q arm", variable=self.orientation_var, value="q").grid(
            row=0, column=1)

        # Contig
        ttk.Label(inp, text="Contig Site:").grid(row=4, column=0, sticky=tk.W, pady=4)
        ttk.Entry(inp, textvariable=self.contig_var, width=14).grid(
            row=4, column=1, sticky=tk.W, pady=4)

        ttk.Button(inp, text="Run Analysis", command=self._run_analysis).grid(
            row=5, column=0, columnspan=2, pady=(12, 0), ipadx=20, ipady=4)

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

    def _on_sheet_selected(self, _event=None):
        """Auto-fill chromosome and orientation from sheet names like '4p+', '13q-'."""
        sheet = self.sheet_var.get()
        m = re.match(r'^(\d+)([pq])[+-]?$', sheet)
        if m:
            self.chromosome_var.set(m.group(1))
            self.orientation_var.set(m.group(2))

    def _browse_excel(self):
        path = filedialog.askopenfilename(
            title="Select Excel file",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        if not path:
            return
        self.excel_path.set(path)
        try:
            xl = pd.ExcelFile(path)
            sheets = xl.sheet_names
            self._sheet_combo["values"] = sheets
            self.sheet_var.set(sheets[0])
            self._set_status(f"Loaded: {os.path.basename(path)}  ({len(sheets)} sheet(s))")
        except Exception as exc:
            messagebox.showerror("Error", f"Could not read Excel file:\n{exc}")

    def _run_analysis(self):
        excel_file  = self.excel_path.get().strip()
        sheet_name  = self.sheet_var.get().strip()
        chromosome  = self.chromosome_var.get().strip()
        orientation = self.orientation_var.get().strip()
        contig_str  = self.contig_var.get().strip()

        if not excel_file:
            messagebox.showwarning("Missing Input", "Please select an Excel file.")
            return
        if not sheet_name:
            messagebox.showwarning("Missing Input", "Please select a sheet.")
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

