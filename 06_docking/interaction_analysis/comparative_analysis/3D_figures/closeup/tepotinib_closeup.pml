
reinitialize
load /Users/omkardaravattu/Desktop/TNIK_IPF_Project/06_docking/interaction_analysis/PLIP_complexes/tepotinib_seed22_complex.pdb, complex
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

select pocket, chain A and resi 29+31+32+33+36+39+52+54+69+108+160+170+171
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
png /Users/omkardaravattu/Desktop/TNIK_IPF_Project/06_docking/interaction_analysis/comparative_analysis/3D_figures/closeup/TNIK_tepotinib_PLIP_closeup.png, dpi=300
save /Users/omkardaravattu/Desktop/TNIK_IPF_Project/06_docking/interaction_analysis/comparative_analysis/3D_figures/closeup/TNIK_tepotinib_PLIP_closeup.pse
quit
