import streamlit as st
import pandas as pd
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import Draw

from src.prediction import predict_compound

# --------------------------------------------------
# PAGE CONFIGURATION
# --------------------------------------------------

st.set_page_config(
    page_title="TNIK-Predict | AI Drug Discovery",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --------------------------------------------------
# DESIGN SYSTEM
# --------------------------------------------------

st.markdown("""

<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Manrope:wght@400;500;600;700;800&display=swap');

:root {
    --navy: #10243A;
    --navy2: #183A55;
    --teal: #00A99D;
    --teal-light: #E6F7F4;
    --muted: #718096;
    --border: #E7EDF2;
    --background: #F5F8FA;
    --white: #FFFFFF;
}

html, body, [class*="css"] {
    font-family: 'DM Sans', sans-serif;
}

.stApp {
    background: var(--background);
}

.block-container {
    max-width: 1500px;
    padding: 2.2rem 3rem 4rem;
}

h1, h2, h3, h4 {
    font-family: 'Manrope', sans-serif !important;
    color: var(--navy) !important;
    letter-spacing: -0.7px;
}

h1 { font-weight: 800 !important; }
h2, h3 { font-weight: 700 !important; }

p, label, .stMarkdown {
    color: #43556A;
}

/* SIDEBAR */
[data-testid="stSidebar"] {
    background: #10243A;
    border-right: 1px solid #263F55;
}

[data-testid="stSidebar"] * {
    color: #F3F8FC !important;
}

[data-testid="stSidebar"] [data-testid="stRadio"] label {
    padding: 12px 14px;
    border-radius: 10px;
    margin-bottom: 5px;
    transition: background .2s ease;
}

[data-testid="stSidebar"] [data-testid="stRadio"] label:hover {
    background: #24445C;
}

/* HERO */
.hero {
    background: linear-gradient(115deg, #10243A 0%, #16445C 58%, #008F89 100%);
    padding: 48px 48px;
    border-radius: 22px;
    color: white;
    margin-bottom: 32px;
    box-shadow: 0 16px 40px rgba(16,36,58,.13);
    position: relative;
    overflow: hidden;
}

.hero h1 {
    color: white !important;
    font-size: 42px !important;
    font-weight: 800;
    line-height: 1.15;
    margin-bottom: 12px;
}

.hero p {
    color: #D7E9EE;
    font-size: 16px;
    max-width: 700px;
    line-height: 1.7;
}

.eyebrow {
    color: #65E1D2;
    font-size: 11px;
    font-weight: 800;
    letter-spacing: 2.2px;
    text-transform: uppercase;
}

/* SECTION LABELS */
.section-label {
    color: #009D93;
    font-size: 11px;
    font-weight: 800;
    letter-spacing: 1.8px;
    text-transform: uppercase;
    margin-bottom: 7px;
}

/* PANELS AND CARDS */
.panel, [data-testid="stVerticalBlockBorderWrapper"] {
    background: white;
    border: 1px solid var(--border);
    border-radius: 17px;
    box-shadow: 0 5px 22px rgba(16,36,58,.045);
}

.panel {
    padding: 25px;
    margin-bottom: 20px;
}

.metric-card, .model-card {
    background: white;
    border: 1px solid var(--border);
    border-radius: 16px;
    padding: 23px;
    min-height: 125px;
    box-shadow: 0 5px 20px rgba(16,36,58,.04);
    transition: transform .2s ease, box-shadow .2s ease;
}

.metric-card:hover, .model-card:hover {
    transform: translateY(-3px);
    box-shadow: 0 10px 25px rgba(16,36,58,.09);
}

.metric-title, .model-name {
    color: #718096;
    font-size: 11px;
    font-weight: 700;
    letter-spacing: .8px;
    text-transform: uppercase;
}

.metric-value, .model-value {
    color: var(--navy);
    font-family: 'Manrope', sans-serif;
    font-size: 29px;
    font-weight: 800;
    margin-top: 10px;
}

.model-card {
    min-height: 155px;
}

/* INPUTS */
.stTextInput input, .stTextArea textarea {
    background: #FFFFFF !important;
    border: 1px solid #DCE5EC !important;
    border-radius: 11px !important;
    padding: 13px 15px !important;
    color: var(--navy) !important;
}

.stTextInput input:focus {
    border-color: var(--teal) !important;
    box-shadow: 0 0 0 3px rgba(0,169,157,.12) !important;
}

/* BUTTONS */
.stButton > button, .stDownloadButton > button {
    background: linear-gradient(110deg, #008F89, #00B3A4);
    color: white !important;
    border: none;
    border-radius: 10px;
    padding: 12px 20px;
    font-weight: 700;
    min-height: 44px;
    box-shadow: 0 5px 13px rgba(0,143,137,.17);
    transition: all .2s ease;
}

.stButton > button:hover, .stDownloadButton > button:hover {
    background: #087E7B;
    border: none;
    transform: translateY(-1px);
}

/* CHARTS AND DATA */
[data-testid="stDataFrame"] {
    border: 1px solid var(--border);
    border-radius: 12px;
    overflow: hidden;
}

[data-testid="stMetric"] {
    background: white;
    border: 1px solid var(--border);
    border-radius: 13px;
    padding: 15px;
}

/* ALERTS */
[data-testid="stAlert"] {
    border-radius: 12px;
    border: 1px solid var(--border);
}

/* TABS */
.stTabs [data-baseweb="tab-list"] {
    gap: 20px;
    border-bottom: 1px solid var(--border);
}

.stTabs [data-baseweb="tab"] {
    color: #718096;
    font-weight: 700;
}

.stTabs [aria-selected="true"] {
    color: #009D93 !important;
}

/* DIVIDERS AND FOOTER */
hr {
    border-color: var(--border);
}

footer {
    visibility: hidden;
}

@media (max-width: 768px) {
    .block-container {
        padding: 1.2rem 1rem 2rem;
    }

    .hero {
        padding: 28px 23px;
    }

    .hero h1 {
        font-size: 30px !important;
    }
}
</style>

""", unsafe_allow_html=True)


# --------------------------------------------------
# SIDEBAR
# --------------------------------------------------

with st.sidebar:
    st.markdown("""
    <div style="padding:18px 5px 25px 5px;">
        <div style="font-size:30px;">🧬</div>
        <div style="font-family:Manrope;font-size:23px;font-weight:800;">
            TNIK-Predict
        </div>
        <div style="font-size:11px;color:#AFC8D5;">
            AI-ASSISTED DRUG DISCOVERY
        </div>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("---")

    page = st.radio(
        "NAVIGATION",
        [
            "Dashboard",
            "Compound Screening",
            "Model Information",
            "About the Project"
        ],
        label_visibility="visible"
    )

    st.markdown("---")

    st.markdown("""
    <div style="font-size:12px;color:#B8CBD5;">
        <b>RESEARCH PLATFORM</b><br><br>
        Target: TNIK<br>
        Research area: Idiopathic Pulmonary Fibrosis<br>
        Models: 4 machine-learning models<br><br>
        <span style="color:#79D8C9;">● Local application</span>
    </div>
    """, unsafe_allow_html=True)


# --------------------------------------------------
# HERO
# --------------------------------------------------

st.markdown("""
<div class="hero">
    <div class="eyebrow">COMPUTATIONAL DRUG DISCOVERY PLATFORM</div>
    <h1>TNIK-Predict</h1>
    <p>
        AI-assisted molecular screening and activity prediction
        for investigating potential TNIK inhibitors in
        idiopathic pulmonary fibrosis research.
    </p>
</div>
""", unsafe_allow_html=True)


# --------------------------------------------------
# DASHBOARD
# --------------------------------------------------

if page == "Dashboard":

    st.markdown('<div class="section-label">PLATFORM OVERVIEW</div>',
                unsafe_allow_html=True)
    st.header("Welcome to TNIK-Predict")

    st.write(
        "Explore molecular structures, estimate TNIK-related activity "
        "and compare predictions from four machine-learning models."
    )

    a, b, c, d = st.columns(4)

    dashboard_metrics = [
        ("TARGET", "TNIK", "Target protein"),
        ("ML MODELS", "04", "Trained models"),
        ("FINGERPRINT", "2048", "Morgan fingerprint bits"),
        ("REFERENCE SET", "453", "Unique compounds")
    ]

    for col, item in zip([a, b, c, d], dashboard_metrics):
        with col:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-title">{item[0]}</div>
                <div class="metric-value">{item[1]}</div>
                <div style="color:#718096;font-size:12px;">{item[2]}</div>
            </div>
            """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    left, right = st.columns([1.2, 1])

    with left:
        st.markdown('<div class="panel">', unsafe_allow_html=True)
        st.subheader("Compound screening")
        st.write(
            "Submit a molecular SMILES representation to generate "
            "predictions and chemical descriptors."
        )
        if st.button("Open compound screening", key="dashboard_screen"):
            st.session_state["navigation_request"] = "Compound Screening"
            st.info("Select Compound Screening from the sidebar.")
        st.markdown('</div>', unsafe_allow_html=True)

    with right:
        st.markdown('<div class="panel">', unsafe_allow_html=True)
        st.subheader("Scientific workflow")
        st.write("01  Molecular input")
        st.write("02  Morgan fingerprint generation")
        st.write("03  Machine-learning prediction")
        st.write("04  Similarity analysis")
        st.write("05  Results interpretation")
        st.markdown('</div>', unsafe_allow_html=True)


# --------------------------------------------------
# COMPOUND SCREENING
# --------------------------------------------------

elif page == "Compound Screening":

    st.markdown('<div class="section-label">MOLECULAR ANALYSIS</div>',
                unsafe_allow_html=True)
    st.header("Compound screening")
    st.write(
        "Enter a valid SMILES representation to analyse a compound "
        "using the trained TNIK prediction models."
    )

    with st.container(border=True):
        st.markdown("### Molecular input")

        smiles = st.text_input(
            "Compound SMILES",
            placeholder="Example: CC(=O)Oc1ccccc1C(=O)O",
            help="Enter a valid molecular SMILES string."
        )

        run_prediction = st.button(
            "Run molecular prediction  →",
            type="primary",
            use_container_width=True
        )

    if run_prediction:

        if not smiles.strip():
            st.warning("Please enter a compound SMILES.")

        else:
            try:
                with st.spinner(
                    "Generating molecular fingerprints and running models..."
                ):
                    result = predict_compound(smiles)

                st.session_state["last_result"] = result
                st.session_state["last_smiles"] = smiles

            except Exception as error:
                st.error(f"Unable to process compound: {error}")

    if "last_result" in st.session_state:

        result = st.session_state["last_result"]
        current_smiles = st.session_state["last_smiles"]

        molecule = Chem.MolFromSmiles(current_smiles)

        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown('<div class="section-label">ANALYSIS RESULTS</div>',
                    unsafe_allow_html=True)
        st.header("Molecular analysis report")

        left, right = st.columns([1, 1.4])

        with left:
            st.markdown('<div class="panel">', unsafe_allow_html=True)
            st.subheader("Chemical structure")
            st.image(
                Draw.MolToImage(molecule, size=(480, 350)),
                use_container_width=True
            )
            st.caption("Canonical SMILES")
            st.code(result["canonical_smiles"])
            st.markdown('</div>', unsafe_allow_html=True)

        with right:
            st.markdown('<div class="panel">', unsafe_allow_html=True)
            st.subheader("Molecular properties")

            properties = result["properties"]

            property_columns = st.columns(2)

            for index, (name, value) in enumerate(properties.items()):
                with property_columns[index % 2]:
                    st.markdown(f"""
                    <div class="metric-card" style="margin-bottom:12px;">
                        <div class="metric-title">{name}</div>
                        <div class="metric-value" style="font-size:22px;">
                            {value}
                        </div>
                    </div>
                    """, unsafe_allow_html=True)

            st.markdown('</div>', unsafe_allow_html=True)

        st.markdown('<div class="section-label">MODEL OUTPUT</div>',
                    unsafe_allow_html=True)
        st.header("Machine-learning predictions")

        predictions = result["predictions"]

        model_columns = st.columns(4)

        for col, (name, values) in zip(model_columns, predictions.items()):
            with col:
                st.markdown(f"""
                <div class="model-card">
                    <div class="model-name">{name}</div>
                    <div class="model-value">
                        {values["predicted_pIC50"]:.3f}
                    </div>
                    <div class="model-unit">Predicted pIC50</div>
                    <div style="color:#078C87;font-size:13px;
                                font-weight:700;margin-top:7px;">
                        {values["estimated_IC50_nM"]:.2f} nM estimated
                    </div>
                </div>
                """, unsafe_allow_html=True)

        st.markdown("<br>", unsafe_allow_html=True)

        chart_data = pd.DataFrame([
            {
                "Model": name,
                "Predicted pIC50": values["predicted_pIC50"]
            }
            for name, values in predictions.items()
        ])

        st.markdown('<div class="panel">', unsafe_allow_html=True)
        st.subheader("Model comparison")

        fig = px.bar(
            chart_data,
            x="Model",
            y="Predicted pIC50",
            text_auto=".3f",
            color="Predicted pIC50",
            color_continuous_scale=["#A8DCD6", "#078C87", "#10263F"]
        )

        fig.update_layout(
            template="plotly_white",
            height=350,
            showlegend=False,
            coloraxis_showscale=False,
            margin=dict(l=10, r=10, t=25, b=10),
            font=dict(family="DM Sans", color="#10263F"),
            yaxis_title="Predicted pIC50",
            xaxis_title=""
        )

        st.plotly_chart(fig, use_container_width=True)
        st.markdown('</div>', unsafe_allow_html=True)

        st.markdown('<div class="section-label">CHEMICAL SPACE</div>',
                    unsafe_allow_html=True)
        st.header("Chemical similarity analysis")

        similarity = result["similarity"]
        similarity_value = similarity["maximum_similarity"]

        sim1, sim2 = st.columns(2)

        with sim1:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-title">Maximum Tanimoto similarity</div>
                <div class="metric-value">{similarity_value:.3f}</div>
                <div style="color:#718096;font-size:12px;">
                    Compared with the reference dataset
                </div>
            </div>
            """, unsafe_allow_html=True)

        with sim2:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-title">Closest reference compound</div>
                <div class="metric-value" style="font-size:21px;">
                    {similarity["closest_compound"]}
                </div>
                <div style="color:#718096;font-size:12px;">
                    ChEMBL compound identifier
                </div>
            </div>
            """, unsafe_allow_html=True)

        if similarity_value < 0.30:
            st.warning(
                "Low similarity to the reference compounds. "
                "The prediction may involve extrapolation beyond "
                "familiar chemical structures. Interpret cautiously. "
                "This is a preliminary similarity flag, not a validated "
                "applicability-domain classification."
            )
        else:
            st.info(
                "Similarity has been calculated. It does not by itself "
                "establish prediction reliability or applicability."
            )

        # --------------------------------------------------
        # DOCKING AND PLIP INTERACTION EVIDENCE
        # --------------------------------------------------
        from rdkit import Chem

        st.markdown('<div class="section-label">STRUCTURAL EVIDENCE</div>',
                    unsafe_allow_html=True)
        st.header("Docking and protein–ligand interactions")

        evidence = None
        submitted_molecule = Chem.MolFromSmiles(
            result["compound"]["canonical_smiles"]
        )

        for candidate in ["ponatinib", "tepotinib"]:
            candidate_evidence = get_compound_evidence(candidate)
            plip = candidate_evidence.get("plip")

            if plip and plip.get("ligand_smiles"):
                reference_molecule = Chem.MolFromSmiles(
                    plip["ligand_smiles"]
                )
                if (submitted_molecule is not None
                        and reference_molecule is not None
                        and Chem.MolToSmiles(submitted_molecule)
                        == Chem.MolToSmiles(reference_molecule)):
                    evidence = candidate_evidence
                    break

        if evidence and evidence.get("plip"):
            plip = evidence["plip"]
            docking = evidence.get("docking")

            if docking:
                st.metric(
                    "AutoDock Vina score",
                    f'{docking["vina_score_kcal_mol"]:.3f} kcal/mol'
                )
                st.caption(
                    "Computational docking score; not an experimental "
                    "binding affinity."
                )

            display_plip_evidence(candidate, plip)

        else:
            st.info(
                "No compound-specific PLIP interaction report is currently "
                "available for this submitted molecule. Docking evidence "
                "is shown only for compounds with matching project reports."
            )

        st.markdown('<div class="panel">', unsafe_allow_html=True)
        st.subheader("Download results")

        report_df = pd.DataFrame([
            {
                "Model": name,
                "Predicted pIC50": values["predicted_pIC50"],
                "Estimated IC50 (nM)": values["estimated_IC50_nM"]
            }
            for name, values in predictions.items()
        ])

        st.download_button(
            "Download prediction data (CSV)",
            data=report_df.to_csv(index=False),
            file_name="TNIK_Predict_results.csv",
            mime="text/csv",
            use_container_width=True
        )

        st.markdown('</div>', unsafe_allow_html=True)

        st.info(
            "Scientific interpretation: pIC50 and IC50 values displayed "
            "here are model-derived estimates, not experimental results. "
            "Predicted TNIK activity does not establish biological efficacy "
            "against idiopathic pulmonary fibrosis."
        )


# --------------------------------------------------
# MODEL INFORMATION
# --------------------------------------------------

elif page == "Model Information":

    st.markdown('<div class="section-label">METHODOLOGY</div>',
                unsafe_allow_html=True)
    st.header("Machine-learning models")

    st.write(
        "TNIK-Predict currently integrates four trained regression "
        "models developed using molecular fingerprints."
    )

    model_info = pd.DataFrame([
        ["Random Forest", "Ensemble decision trees", "2048-bit Morgan"],
        ["SVR", "Support vector regression", "2048-bit Morgan"],
        ["Gradient Boosting", "Boosted decision trees", "2048-bit Morgan"],
        ["Extra Trees", "Randomized ensemble trees", "2048-bit Morgan"]
    ], columns=["Model", "Method", "Input features"])

    st.dataframe(model_info, hide_index=True, use_container_width=True)

    st.subheader("Molecular representation")

    st.markdown("""
    - Fingerprint type: Morgan fingerprint
    - Radius: 2
    - Fingerprint size: 2048 bits
    - Training reference: 453 unique compounds
    - Biological activity endpoint: pIC50
    - Data source: ChEMBL-derived TNIK activity data
    """)

    st.warning(
        "Model predictions are intended for computational research "
        "and hypothesis generation, not clinical decision-making."
    )


# --------------------------------------------------
# ABOUT
# --------------------------------------------------

elif page == "About the Project":

    st.markdown('<div class="section-label">ABOUT</div>',
                unsafe_allow_html=True)
    st.header("About TNIK-Predict")

    st.write("""
    TNIK-Predict is a computational drug-discovery research platform
    designed to investigate potential TNIK inhibitors in the context
    of idiopathic pulmonary fibrosis (IPF).

    The platform combines molecular structure processing, chemical
    fingerprints, trained machine-learning models and chemical
    similarity analysis.
    """)

    st.subheader("Research objective")

    st.write(
        "To support early-stage computational screening and prioritisation "
        "of compounds for further investigation of TNIK-related activity."
    )

    st.subheader("Scientific limitations")

    st.write("""
    1. Predicted activity is not experimentally confirmed activity.
    2. Model estimates may be less reliable for unfamiliar chemical structures.
    3. TNIK inhibition does not automatically imply therapeutic benefit in IPF.
    4. Experimental validation is required to establish biological activity.
    """)


# --------------------------------------------------
# FOOTER
# --------------------------------------------------

st.markdown("""
<div class="footer">
    <b>TNIK-Predict</b> · AI-assisted computational drug discovery<br>
    Research prototype · Not intended for clinical decisions<br>
    © TNIK-Predict Research Project
</div>
""", unsafe_allow_html=True)
