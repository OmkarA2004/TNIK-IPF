from pathlib import Path
from collections import Counter

source = Path("02_raw_data/pdb_original/6RA7.pdb")
output = Path("06_docking/molecular_dynamics/structures/TNIK_6RA7_protein_raw.pdb")
report = Path("06_docking/molecular_dynamics/logs/protein_extraction_report.txt")

standard = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL"
}

kept = []
residues = set()
altloc_count = Counter()

for line in source.read_text().splitlines():
    if not line.startswith(("ATOM  ", "HETATM")):
        continue

    if line[21] != "A":
        continue

    resname = line[17:20].strip()
    altloc = line[16]

    if resname not in standard | {"TPO"}:
        continue

    if altloc not in (" ", "A"):
        continue

    if altloc != " ":
        altloc_count[resname] += 1

    kept.append(line[:16] + " " + line[17:])
    residues.add((resname, line[22:27].strip()))

kept.extend(["TER", "END"])
output.write_text("\n".join(kept) + "\n")

counts = Counter(r[0] for r in residues)
report.parent.mkdir(parents=True, exist_ok=True)
report.write_text(
    f"Source: {source}\nOutput: {output}\n"
    f"Unique residues: {len(residues)}\n"
    f"Residue counts: {dict(counts)}\n"
    f"TPO residues: {sorted(r for r in residues if r[0]=='TPO')}\n"
    f"Alternate conformers selected: {dict(altloc_count)}\n"
)

print("Protein extraction completed.")
print("Local report saved:", report)
print("TPO residues retained:", sorted(r for r in residues if r[0] == "TPO"))
