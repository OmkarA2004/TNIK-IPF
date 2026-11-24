load 07_explainability/plip/Tepotinib_analogue_55/TNIK_analogue_55_complex_LIG.pdb, complex

hide everything
remove solvent

show cartoon, polymer.protein
color lightblue, polymer.protein

select ligand, resn LIG
show sticks, ligand
color orange, ligand

select contact_residues, polymer.protein within 4.5 of ligand
show sticks, contact_residues
color cyan, contact_residues

select hb_residues, resi 54+106+108 and polymer.protein
show sticks, hb_residues
color yellow, hb_residues

distance hb_lys54, (ligand), (resi 54 and polymer.protein), 3.5
distance hb_glu106, (ligand), (resi 106 and polymer.protein), 3.5
distance hb_cys108, (ligand), (resi 108 and polymer.protein), 3.5

set dash_color, yellow
set dash_width, 2
set ray_opaque_background, off
bg_color white

orient ligand
zoom ligand, 12

set label_size, 18
set ray_shadows, 0

ray 1800, 1400
png 09_results/figures/TNIK_analogue_55_3D_interactions.png, dpi=300

save 09_results/figures/TNIK_analogue_55_3D_interactions.pse
