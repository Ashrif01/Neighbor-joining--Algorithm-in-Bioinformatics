# Neighbor-Joining (NJ) Phylogenetic Tree Construction

This project implements the **Neighbor-Joining (NJ)** algorithm in Python for reconstructing phylogenetic trees from distance matrices in **PHYLIP** format. It compares the performance and output accuracy against **QuickTree**, a standard NJ tool written in C.

## 📁 Project Structure

```
neighbor-joining-bioinformatics/
├── code/
│   ├── nj.py                  # Custom NJ algorithm implementation
│   └── rfdist.py              # RF-distance comparison script (from Project 4)
│
├── data/
│   ├── unique_distance_matrices/  # 14 input .phy files
│   ├── nj_outputs/                # NJ-generated trees
│   └── quicktree_outputs/        # QuickTree-generated trees
│
├── results/
│   ├── final_report_table.csv     # Runtime + speedup comparison table
│   ├── nj_runtime_results.csv     # NJ runtimes
│   ├── quicktree_runtime.csv      # QuickTree runtimes
│   └── rfdist.csv                 # RF-distance results
│
├── docs/
│   ├── NJ_Project_Report.pdf      # Final report (2-page, as required)
│   └── NJ_Project_Report.docx
│
├── README.md
└── LICENSE
```

---

## 🚀 How to Run the NJ Code

```bash
python code/nj.py path/to/input.phy path/to/output.nwk
```

To visualize a generated tree (optional):

```python
from Bio import Phylo
tree = Phylo.read("output.nwk", "newick")
Phylo.draw(tree)
```

---

## 🧪 Experiments

- 14 distance matrices from `unique_distance_matrices/` were processed
- Runtime of NJ vs QuickTree was measured
- RF-distances were computed using `rfdist.py` comparing topology similarity

---

## 📊 Report Results

Final report includes:
- Correct output tree from `example_slide4.phy`
- Python NJ implementation details
- System info + timing method
- Full results table: QuickTree vs NJ runtime + speedups
- RF-distance comparison for validation

---

## ⚠️ Notes

- RapidNJ was excluded due to platform compatibility issues (MacOS / Windows).
- QuickTree comparison was done using PowerShell timing (`Measure-Command`).

---

## 👨‍💻 Author

**Al Ashrif Bin Ahamed**  
Aarhus University  
Spring 2025 – Algorithm in Bioinformatics