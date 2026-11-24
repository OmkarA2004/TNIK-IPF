from pathlib import Path
import csv
from datetime import date

ROOT = Path.cwd()
RESULTS = ROOT / "09_results"
MD_RESULTS = RESULTS / "molecular_dynamics"
OUTPUT = RESULTS / "TNIK_IPF_PUBLICATION_READINESS_AUDIT.md"

required_files = [
    "09_results/TNIK_IPF_MASTER_RESULTS_REPORT.md",
    "09_results/TNIK_IPF_MASTER_EVIDENCE_INDEX.csv",
    "09_results/TNIK_ML_final_metrics_table.csv",
    "09_results/TNIK_docking_final_evidence_table.csv",
    "09_results/molecular_dynamics/TNIK_MD_final_evidence_table.csv",
    "09_results/TNIK_IPF_integrated_results_summary.md",
    "09_results/TNIK_IPF_evidence_inventory.md",
]

lines = []
lines.append("# TNIK–IPF Publication Readiness Audit")
lines.append("")
lines.append(f"Audit date: {date.today().isoformat()}")
lines.append("")
lines.append("## 1. Existing evidence files")
lines.append("")
lines.append("| File | Status |")
lines.append("|---|---|")

for item in required_files:
    path = ROOT / item
    status = "FOUND" if path.exists() else "MISSING"
    lines.append(f"| `{item}` | {status} |")

lines.append("")
lines.append("## 2. Evidence index audit")
lines.append("")

index_path = ROOT / "09_results/TNIK_IPF_MASTER_EVIDENCE_INDEX.csv"

if index_path.exists():
    with open(index_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    lines.append(f"- Evidence-index entries: **{len(rows)}**")

    missing = []
    for row in rows:
        values = list(row.values())
        possible_paths = [
            value.strip() for value in values
            if value and (
                "/" in value or value.endswith((".csv", ".md", ".png",
                                                 ".pdb", ".pdbqt", ".xtc"))
            )
        ]

        for value in possible_paths:
            candidate = ROOT / value
            if not candidate.exists() and value not in missing:
                missing.append(value)

    if missing:
        lines.append(f"- Paths requiring manual review: **{len(missing)}**")
        for item in missing:
            lines.append(f"  - `{item}`")
    else:
        lines.append("- No missing paths detected by this basic path check.")
else:
    lines.append("- Evidence index is missing; this check could not be performed.")

lines.append("")
lines.append("## 3. Existing manuscript figures")
lines.append("")
figure_dirs = [
    ROOT / "09_results",
    ROOT / "06_docking/interaction_analysis/comparative_analysis",
    ROOT / "06_docking/molecular_dynamics/analysis",
]

figures = []
for folder in figure_dirs:
    if folder.exists():
        figures.extend(folder.rglob("*.png"))

figures = sorted(set(figures))

if figures:
    lines.append(f"PNG figures found: **{len(figures)}**")
    for figure in figures:
        lines.append(f"- `{figure.relative_to(ROOT)}`")
else:
    lines.append("No PNG figures were found in the checked locations.")

lines.append("")
lines.append("## 4. Results that can form the current paper")
lines.append("")
lines.append("The following evidence has already been generated and can be considered "
             "for manuscript preparation, subject to final verification:")
lines.append("")
lines.append("1. TNIK ChEMBL bioactivity data curation and pIC50 dataset.")
lines.append("2. RDKit Morgan fingerprint generation.")
lines.append("3. Random Forest baseline model and held-out test evaluation.")
lines.append("4. Random-split and scaffold-split cross-validation.")
lines.append("5. High-confidence subset evaluation.")
lines.append("6. Y-randomisation and applicability-domain assessment.")
lines.append("7. TNIK 6RA7 structure preparation and JWK redocking.")
lines.append("8. Virtual screening and candidate prioritisation.")
lines.append("9. PLIP interaction analysis, subject to manual interaction review.")
lines.append("10. Preliminary restrained MD analyses for Ponatinib and Tepotinib.")

lines.append("")
lines.append("## 5. Work requiring attention before submission")
lines.append("")
lines.append("- Verify all manuscript numbers against the final CSV tables.")
lines.append("- Explain the weaker and variable scaffold-validation performance.")
lines.append("- Confirm the docking benchmark is not described as virtual-screening enrichment.")
lines.append("- Manually review questionable PLIP interaction assignments.")
lines.append("- Clearly disclose the short MD timescales and heavy-atom restraints.")
lines.append("- Check the exact provenance and assay context of the Ponatinib Kd observation.")
lines.append("- Check whether a credible known-active/inactive docking benchmark can be assembled.")
lines.append("- Confirm figure labels, resolution, captions and consistency.")
lines.append("- Check all references and software versions.")

lines.append("")
lines.append("## 6. Important scientific limitations")
lines.append("")
lines.append("- No experimental TNIK validation is available.")
lines.append("- Docking scores do not establish experimental binding affinity or inhibition.")
lines.append("- The Ponatinib cell-lysate Kinobeads Kd observation is not a purified-TNIK IC50 result.")
lines.append("- The current short restrained MD simulations do not establish long-term unrestrained stability.")
lines.append("- Computational prioritisation does not establish IPF therapeutic efficacy.")
lines.append("")
lines.append("## 7. Recommended manuscript framing")
lines.append("")
lines.append("Present this work as an exploratory computational drug-repurposing and "
             "candidate-prioritisation study targeting TNIK in the context of IPF.")
lines.append("")
lines.append("Do not describe Ponatinib or Tepotinib as experimentally confirmed TNIK inhibitors.")
lines.append("")
lines.append("## 8. Next actions")
lines.append("")
lines.append("1. Review this audit and resolve missing or inconsistent evidence.")
lines.append("2. Finalise the manuscript figure and table selection.")
lines.append("3. Prepare the Materials and Methods section from the actual scripts and parameters.")
lines.append("4. Draft Results using only verified outputs.")
lines.append("5. Draft Discussion, limitations and conclusion.")
lines.append("6. Select a student journal or conference and check its author guidelines.")
lines.append("")
lines.append("## 9. Audit interpretation")
lines.append("")
lines.append("This report is a file-based checklist, not an independent scientific validation. "
             "A FOUND status confirms file presence only, not correctness of its scientific contents.")

OUTPUT.write_text("\n".join(lines), encoding="utf-8")

print("\nPUBLICATION READINESS AUDIT CREATED")
print("-----------------------------------")
print(f"Output: {OUTPUT}")
print(f"Required evidence files found: "
      f"{sum((ROOT / x).exists() for x in required_files)}/{len(required_files)}")
print(f"Figures found: {len(figures)}")
print(f"Evidence-index entries: {len(rows) if index_path.exists() else 'Not available'}")
print("\nPlease open the generated Markdown report and share the Terminal output.")
