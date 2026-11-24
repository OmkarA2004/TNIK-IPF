from pathlib import Path
import subprocess

ROOT = Path.cwd()
BASE = ROOT / "06_docking/interaction_analysis"
COMPLEX_DIR = BASE / "PLIP_complexes"
OUT = BASE / "comparative_analysis/3D_figures/closeup"

OUT.mkdir(parents=True, exist_ok=True)

PYMOL = Path("/opt/homebrew/bin/pymol")

systems = {
    "Ponatinib": {
        "complex": COMPLEX_DIR / "ponatinib_seed11_complex.pdb",
        "residues": [31, 36, 39, 54, 107, 115, 170, 171],
    },
    "Tepotinib": {
        "complex": COMPLEX_DIR / "tepotinib_seed22_complex.pdb",
        "residues": [29, 31, 32, 33, 36, 39, 52, 54, 69, 108, 160, 170, 171],
    },
}

if not PYMOL.exists():
    raise FileNotFoundError(f"PyMOL not found: {PYMOL}")

for name, data in systems.items():
    pdb = data["complex"]

    if not pdb.exists():
        raise FileNotFoundError(f"Complex file missing: {pdb}")

    residues = "+".join(map(str, data["residues"]))
    safe_name = name.lower()

    image = OUT / f"TNIK_{safe_name}_PLIP_closeup.png"
    session = OUT / f"TNIK_{safe_name}_PLIP_closeup.pse"
    pml = OUT / f"{safe_name}_closeup.pml"

    commands = f"""
reinitialize
load {pdb}, complex
remove solvent

hide everything
show cartoon, chain A
color gray80, chain A
set cartoon_transparency, 0.25

select ligand, resn MOL
show sticks, ligand
set stick_radius, 0.24, ligand
color cyan, ligand
color red, ligand and elem O
color blue, ligand and elem N
color yellow, ligand and elem S
color forest, ligand and elem F

select pocket, chain A and resi {residues}
select pocket_sidechains, pocket and not name N+CA+C+O
show sticks, pocket_sidechains
set stick_radius, 0.18, pocket_sidechains
color orange, pocket_sidechains

# Label the alpha-carbon position of each highlighted residue.
label pocket and name CA, "%s%s" % (resn, resi)
set label_size, 20
set label_color, black
set label_outline_color, white
set label_position, [0, 0, 2]

bg_color white
set orthoscopic, 1
set antialias, 2
set ray_opaque_background, 1
set ray_trace_mode, 1
set depth_cue, 0

orient ligand
zoom (ligand or pocket), 5

deselect
ray 2000, 1600
png {image}, dpi=300
save {session}
quit
"""

    pml.write_text(commands)

    print(f"\nCreating close-up figure: {name}")

    result = subprocess.run(
        [str(PYMOL), "-cq", str(pml)],
        capture_output=True,
        text=True,
        timeout=300,
    )

    if result.returncode != 0:
        print(result.stdout[-3000:])
        print(result.stderr[-3000:])
        raise RuntimeError(f"PyMOL failed for {name}")

    if not image.exists() or not session.exists():
        print(result.stdout[-3000:])
        print(result.stderr[-3000:])
        raise RuntimeError(f"Output missing for {name}")

    print("Image:", image.name, f"({image.stat().st_size:,} bytes)")
    print("Session:", session.name)

print("\n===== CLOSE-UP FIGURES GENERATED =====")
for file in sorted(OUT.iterdir()):
    if file.is_file():
        print(file.name, f"{file.stat().st_size:,} bytes")
