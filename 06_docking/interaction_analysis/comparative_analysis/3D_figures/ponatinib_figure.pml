
reinitialize
load /Users/omkardaravattu/Desktop/TNIK_IPF_Project/06_docking/interaction_analysis/PLIP_complexes/ponatinib_seed11_complex.pdb, complex
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

select pocket, chain A and resi 31+36+39+54+107+115+170+171
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
png /Users/omkardaravattu/Desktop/TNIK_IPF_Project/06_docking/interaction_analysis/comparative_analysis/3D_figures/TNIK_Ponatinib_PLIP.png, dpi=300

save /Users/omkardaravattu/Desktop/TNIK_IPF_Project/06_docking/interaction_analysis/comparative_analysis/3D_figures/TNIK_Ponatinib_PLIP.pse
quit
