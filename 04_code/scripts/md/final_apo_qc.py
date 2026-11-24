from pathlib import Path

files = {
    "NPT checkpoint": "06_docking/molecular_dynamics/parameters/TNIK_npt_crescale.cpt",
    "NPT trajectory": "06_docking/molecular_dynamics/parameters/TNIK_npt_crescale.xtc",
    "NPT topology": "06_docking/molecular_dynamics/parameters/TNIK_npt_600ps.tpr",
    "RMSD": "06_docking/molecular_dynamics/analysis/TNIK_apo_backbone_RMSD.xvg",
    "Radius of gyration": "06_docking/molecular_dynamics/analysis/TNIK_apo_Rg.xvg",
    "RMSF": "06_docking/molecular_dynamics/analysis/TNIK_apo_RMSF.xvg",
    "DSSP": "06_docking/molecular_dynamics/analysis/TNIK_apo_secondary_structure.xvg",
    "Hydrogen bonds": "06_docking/molecular_dynamics/analysis/TNIK_apo_hbond_number.xvg",
    "DSSP plot": "09_results/molecular_dynamics/TNIK_apo_DSSP_composition.png",
    "Hydrogen-bond plot": "09_results/molecular_dynamics/TNIK_apo_hbond_timeseries.png",
    "RMSD plot": "09_results/molecular_dynamics/TNIK_apo_RMSD.png",
    "Rg plot": "09_results/molecular_dynamics/TNIK_apo_Rg.png",
    "RMSF plot": "09_results/molecular_dynamics/TNIK_apo_RMSF_original_numbering.png",
}

print("APO-TNIK FINAL FILE CHECK")
print("=" * 45)

missing = []

for name, path in files.items():
    p = Path(path)
    status = "FOUND" if p.is_file() and p.stat().st_size > 0 else "MISSING"
    print(f"{name:22s}: {status}")
    if status == "MISSING":
        missing.append(name)

print("\nQC result:", "ALL FILES FOUND" if not missing else "FILES MISSING")
