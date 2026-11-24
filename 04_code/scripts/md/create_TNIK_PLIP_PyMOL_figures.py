from pathlib import Path
import subprocess

ROOT = Path.cwd()

BASE = ROOT / "06_docking/interaction_analysis"
COMPLEX_DIR = BASE / "PLIP_complexes"
RESULT_DIR = BASE / "PLIP_results"
OUT = BASE / "comparative_analysis/3D_figures"

OUT.mkdir(parents=True, exist_ok=True)

# Residues confirmed by the PLIP XML reports and validated against
# the original TNIK PDB residue numbering.
systems = {
    "Ponatinib": {
        "complex": COMPLEX_DIR / "ponatinib_seed11_complex.pdb",
        "residues": [31, 36, 39, 54, 107, 115, 170, 171],
        "session": OUT / "TNIK_Ponatinib_PLIP.pse",
        "image": OUT / "TNIK_Ponatinib_PLIP.png",
        "script": OUT / "ponatinib_figure.pml",
    },
    "Tepotinib": {
        "complex": COMPLEX_DIR / "tepotinib_seed22_complex.pdb",
        "residues": [29, 31, 32, 33, 36, 39, 52, 54, 69, 108, 160, 170, 171],
        "session": OUT / "TNIK_Tepotinib_PLIP.pse",
        "image": OUT / "TNIK_Tepotinib_PLIP.png",
        "script": OUT / "tepotinib_figure.pml",
    },
}

pymol = Path("/opt/homebrew/bin/pymol")

if not pymol.exists():
    raise FileNotFoundError(f"PyMOL executable not found: {pymol}")

for ligand_name, data in systems.items():
    complex_file = data["complex"]

    if not complex_file.exists():
        raise FileNotFoundError(f"Missing complex: {complex_file}")

    residue_selection = "+".join(map(str, data["residues"]))

    commands = f"""
reinitialize
load {complex_file}, complex
remove solvent

hide everything
show cartoon, chain A
color gray80, chain A

select ligand, resn MOL
show sticks, ligand
color cyan, ligand
color red, ligand and elem O
color blue, ligand and elem N
color yellow, ligand and elem S
color forest, ligand and elem F

select pocket, chain A and resi {residue_selection}
show sticks, pocket
color orange, pocket

set cartoon_transparency, 0.15
set stick_radius, 0.19
set sphere_scale, 0.25
set ray_opaque_background, 1
set antialias, 2
set ray_trace_mode, 1
set label_size, 18
set label_color, black
set dash_width, 2.5
set bg_rgb, [1, 1, 1]

orient ligand
zoom ligand, 12

# Keep the ligand and highlighted binding-pocket residues visible.
deselect
hide labels
set orthoscopic, 1

ray 2400, 1800
png {data["image"]}, dpi=300

save {data["session"]}
quit
"""

    data["script"].write_text(commands)

    print(f"\nGenerating {ligand_name} figure...")

    result = subprocess.run(
        [str(pymol), "-cq", str(data["script"])],
        capture_output=True,
        text=True,
        timeout=300,
    )

    if result.returncode != 0:
        print(result.stdout[-4000:])
        print(result.stderr[-4000:])
        raise RuntimeError(f"PyMOL failed for {ligand_name}")

    if not data["image"].exists() or not data["session"].exists():
        print(result.stdout[-4000:])
        print(result.stderr[-4000:])
        raise RuntimeError(f"Expected PyMOL output missing for {ligand_name}")

    print("  Image:", data["image"].name)
    print("  Session:", data["session"].name)
    print("  Image size:", data["image"].stat().st_size, "bytes")
    print("  Session size:", data["session"].stat().st_size, "bytes")

print("\nBoth PyMOL figures generated successfully.")
print("Output directory:", OUT)
