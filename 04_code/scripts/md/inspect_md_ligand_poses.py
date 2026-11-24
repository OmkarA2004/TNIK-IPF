from pathlib import Path
import re

folder = Path("06_docking/reproducibility")

for ligand in ["ponatinib", "tepotinib"]:
    print(f"\n{ligand.upper()}")
    print("-" * 45)

    for seed in [11, 22, 33, 44, 55]:
        path = folder / f"{ligand}_seed_{seed}_out.pdbqt"

        if not path.exists():
            print(f"Seed {seed}: MISSING")
            continue

        lines = path.read_text().splitlines()
        score = "Not found"

        for line in lines:
            if "REMARK VINA RESULT:" in line:
                match = re.search(
                    r"REMARK VINA RESULT:\s*([-0-9.]+)", line
                )
                if match:
                    score = match.group(1)
                break

        atoms = sum(
            line.startswith(("ATOM", "HETATM")) for line in lines
        )

        print(f"Seed {seed}: score={score} kcal/mol; atoms={atoms}")

