
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pathlib import Path

root = Path.home() / "Desktop" / "TNIK_IPF_Project"
out = root / "10_report" / "final_ppt"
out.mkdir(parents=True, exist_ok=True)

prs = Presentation()
prs.slide_width = Inches(13.33)
prs.slide_height = Inches(7.5)

BLUE = RGBColor(31, 78, 121)
DARK = RGBColor(40, 40, 40)

def add_slide(title, bullets, image=None):
    slide = prs.slides.add_slide(prs.slide_layouts[6])

    header = slide.shapes.add_textbox(
        Inches(0.55), Inches(0.35), Inches(12.2), Inches(0.65)
    )
    p = header.text_frame.paragraphs[0]
    p.text = title
    p.font.size = Pt(27)
    p.font.bold = True
    p.font.color.rgb = BLUE

    body_width = 6.0 if image else 12.0
    body = slide.shapes.add_textbox(
        Inches(0.7), Inches(1.35), Inches(body_width), Inches(5.5)
    )
    tf = body.text_frame
    tf.word_wrap = True

    for i, item in enumerate(bullets):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = item
        p.font.size = Pt(19)
        p.font.color.rgb = DARK
        p.space_after = Pt(15)

    if image and image.is_file():
        slide.shapes.add_picture(
            str(image), Inches(7.0), Inches(1.45),
            width=Inches(5.7)
        )

    footer = slide.shapes.add_textbox(
        Inches(0.6), Inches(7.12), Inches(12), Inches(0.2)
    )
    fp = footer.text_frame.paragraphs[0]
    fp.text = "TNIK–IPF Computational Study | Supervisor Review"
    fp.font.size = Pt(9)
    fp.font.color.rgb = RGBColor(120, 120, 120)

# 1. Title
slide = prs.slides.add_slide(prs.slide_layouts[6])
box = slide.shapes.add_textbox(
    Inches(1), Inches(1.6), Inches(11.3), Inches(3)
)
tf = box.text_frame
p = tf.paragraphs[0]
p.text = "Computational Prioritisation of Approved Drugs for TNIK"
p.font.size = Pt(34)
p.font.bold = True
p.font.color.rgb = BLUE
p.alignment = PP_ALIGN.CENTER

p = tf.add_paragraph()
p.text = "In the Context of Idiopathic Pulmonary Fibrosis"
p.font.size = Pt(25)
p.alignment = PP_ALIGN.CENTER

p = tf.add_paragraph()
p.text = "\nProject progress and supervisor review\n10 November 2026"
p.font.size = Pt(18)
p.alignment = PP_ALIGN.CENTER

# 2. Aim and workflow
add_slide("1. Research Aim and Computational Workflow", [
    "Aim: prioritise approved drugs for further investigation against TNIK in IPF.",
    "ChEMBL bioactivity retrieval and data preparation.",
    "Morgan fingerprints and Random Forest regression.",
    "SHAP interpretation and applicability-domain assessment.",
    "Approved-drug docking, interaction analysis and preliminary MD equilibration."
])

# 3. Dataset
add_slide("2. TNIK Bioactivity Dataset", [
    "2,177 ChEMBL activity records retrieved.",
    "554 exact IC50 measurements selected for pIC50 analysis.",
    "453 unique compounds after deduplication.",
    "2,048-bit Morgan fingerprints generated.",
    "Train/test split: 362 training and 91 test compounds."
])

# 4. Model performance
add_slide("3. Random Forest Model Performance", [
    "Test R² = 0.4914.",
    "Test RMSE = 0.7976.",
    "Test MAE = 0.5917.",
    "Five-fold random CV mean R² = 0.5439 ± 0.0931.",
    "Scaffold validation was weaker, indicating limited generalisation to unseen chemical scaffolds."
], root / "09_results/figures/TNIK_RF_actual_vs_predicted.png")

# 5. Explainability
add_slide("4. Model Interpretation", [
    "SHAP feature-importance and summary analyses generated.",
    "Applicability-domain assessment completed.",
    "86 of 91 test compounds were inside the assessed domain.",
    "Five compounds were outside; their higher prediction error requires cautious interpretation.",
    "Fingerprint features do not independently establish molecular mechanisms."
], root / "09_results/figures/TNIK_SHAP_summary.png")

# 6. Docking
add_slide("5. Approved-Drug Docking Prioritisation", [
    "TNIK crystal structure: PDB 6RA7.",
    "JWK redocking heavy-atom RMSD: 0.782 Å.",
    "Shortlisted: Ponatinib, Tepotinib, Entrectinib and Nilotinib.",
    "Docking scores ranged from −10.523 to −11.149 kcal/mol.",
    "Scores are computational predictions, not proof of TNIK inhibition."
], root / "09_results/figures/TNIK_four_drug_comparative_interactions.png")

# 7. Experimental evidence and MD
add_slide("6. Experimental Affinity and Molecular Dynamics", [
    "Ponatinib has a reported TNIK Kd of 385 nM in a Kinobeads cell-lysate competition assay.",
    "Tepotinib, Entrectinib and Nilotinib have weak or censored affinity records.",
    "MD included 500 ps NVT and 1 ns NPT restrained equilibration.",
    "These simulations were restrained equilibration, not production MD.",
    "Docking, affinity and MD evidence must be interpreted separately."
], root / "06_docking/molecular_dynamics/analysis/npt_stability/plots/TNIK_three_complexes_backbone_RMSD.png")

# 8. Limitations and future work
add_slide("7. Limitations and Next Steps", [
    "Moderate model performance and weaker scaffold validation.",
    "Docking and short restrained MD do not establish biological efficacy.",
    "ChemBERTa development and model comparison remain pending.",
    "Analogue generation and ADMET assessment remain pending.",
    "Next: supervisor feedback, final evidence review and experimental validation."
])

# Save
file = out / "TNIK_IPF_Supervisor_Presentation.pptx"
prs.save(file)
print(f"Presentation created: {file}")
print(f"Slides: {len(prs.slides)}")

