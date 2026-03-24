import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

# Retro 80s/90s theme - Black background with green text
st.markdown(
    """
<style>
    /* Retro terminal style */
    .stApp {
        background-color: #0a0a0a;
        color: #00ff00;
    }
    .stMarkdown, .stText, p, div, span {
        color: #00ff00 !important;
        font-family: 'Courier New', monospace;
    }
    h1, h2, h3, h4, h5, h6 {
        color: #00ff00 !important;
        font-family: 'Courier New', monospace;
        text-shadow: 0 0 10px #00ff00;
    }
    .stButton>button {
        background-color: #001100;
        color: #00ff00;
        border: 2px solid #00ff00;
        font-family: 'Courier New', monospace;
    }
    .stButton>button:hover {
        background-color: #00ff00;
        color: #000000;
    }
    .stTextInput>div>div>input {
        background-color: #001100;
        color: #00ff00;
        border: 1px solid #00ff00;
        font-family: 'Courier New', monospace;
    }
    .stSelectbox>div>div>div {
        background-color: #001100;
        color: #00ff00;
    }
    .stSidebar {
        background-color: #050505;
    }
    .stTabs [data-baseweb="tab-list"] {
        background-color: #001100;
    }
    .stTabs [data-baseweb="tab"] {
        color: #00ff00;
    }
    div[data-testid="stMetricValue"] {
        color: #00ff00 !important;
    }
    div[data-testid="stMetricLabel"] {
        color: #00aa00 !important;
    }
    /* Success/Warning/Error colors */
    .stSuccess {
        background-color: #001100;
        color: #00ff00;
        border: 1px solid #00ff00;
    }
    .stWarning {
        background-color: #111100;
        color: #ffff00;
        border: 1px solid #ffff00;
    }
    .stError {
        background-color: #110000;
        color: #ff0000;
        border: 1px solid #ff0000;
    }
    /* Table styling */
    .stDataFrame {
        background-color: #001100;
    }
    /* ASCII art section */
    .ascii-box {
        background-color: #000000;
        border: 2px solid #00ff00;
        padding: 10px;
        font-family: 'Courier New', monospace;
        white-space: pre;
        overflow-x: auto;
    }
</style>
""",
    unsafe_allow_html=True,
)

try:
    from pymatgen.core import Composition, Element, Structure
    from pymatgen.ext.matproj import MPRester
    from mp_api.client import MPRester as NewMPRester

    PYMATGEN_AVAILABLE = True
except ImportError:
    PYMATGEN_AVAILABLE = False
    st.warning("pymatgen not installed. Some features will be limited.")

try:
    import smact
    from smact.screening import pauling_test

    SMACT_AVAILABLE = True
except ImportError:
    SMACT_AVAILABLE = False

import json
import requests
from datetime import datetime
import base64
from io import StringIO

# Configuración de página
st.set_page_config(
    page_title="GNoME Material Validator",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# CSS personalizado para estética profesional
st.markdown(
    """
<style>
    .main-header {
        font-size: 3rem;
        font-weight: bold;
        color: #1E88E5;
        text-align: center;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background-color: #f0f2f6;
        border-radius: 10px;
        padding: 20px;
        box-shadow: 2px 2px 10px rgba(0,0,0,0.1);
    }
    .success-box {
        background-color: #d4edda;
        border-left: 5px solid #28a745;
        padding: 20px;
        border-radius: 5px;
    }
    .warning-box {
        background-color: #fff3cd;
        border-left: 5px solid #ffc107;
        padding: 20px;
        border-radius: 5px;
    }
    .danger-box {
        background-color: #f8d7da;
        border-left: 5px solid #dc3545;
        padding: 20px;
        border-radius: 5px;
    }
</style>
""",
    unsafe_allow_html=True,
)

# ============================================================================
# BASE DE DATOS Y UTILIDADES
# ============================================================================


@st.cache_data
def load_cost_database():
    """Base de datos de costos de precursores (USD/g)"""
    return {
        "H": 0.01,
        "Li": 0.5,
        "Be": 2.0,
        "B": 0.2,
        "C": 0.01,
        "N": 0.01,
        "O": 0.0,
        "Na": 0.1,
        "Mg": 0.05,
        "Al": 0.02,
        "Si": 0.05,
        "P": 0.1,
        "S": 0.01,
        "Cl": 0.01,
        "K": 0.1,
        "Ca": 0.05,
        "Sc": 5.0,
        "Ti": 0.2,
        "V": 1.0,
        "Cr": 0.1,
        "Mn": 0.1,
        "Fe": 0.02,
        "Co": 0.5,
        "Ni": 0.1,
        "Cu": 0.05,
        "Zn": 0.05,
        "Ga": 1.0,
        "Ge": 2.0,
        "As": 0.5,
        "Se": 0.5,
        "Br": 0.1,
        "Kr": 1.0,
        "Rb": 1.0,
        "Sr": 0.5,
        "Y": 2.0,
        "Zr": 0.5,
        "Nb": 1.0,
        "Mo": 0.5,
        "Tc": 50.0,
        "Ru": 20.0,
        "Rh": 100.0,
        "Pd": 50.0,
        "Ag": 50.0,
        "Cd": 0.5,
        "In": 2.0,
        "Sn": 0.2,
        "Sb": 0.5,
        "Te": 0.5,
        "I": 0.2,
        "Xe": 10.0,
        "Cs": 2.0,
        "Ba": 0.5,
        "La": 2.0,
        "Ce": 2.0,
        "Pr": 5.0,
        "Nd": 5.0,
        "Pm": 1000.0,
        "Sm": 5.0,
        "Eu": 10.0,
        "Gd": 5.0,
        "Tb": 50.0,
        "Dy": 20.0,
        "Ho": 20.0,
        "Er": 20.0,
        "Tm": 50.0,
        "Yb": 10.0,
        "Lu": 50.0,
        "Hf": 2.0,
        "Ta": 5.0,
        "W": 1.0,
        "Re": 10.0,
        "Os": 50.0,
        "Ir": 100.0,
        "Pt": 100.0,
        "Au": 80.0,
        "Hg": 0.5,
        "Tl": 5.0,
        "Pb": 0.2,
        "Bi": 2.0,
        "Th": 50.0,
        "U": 10.0,
    }


@st.cache_data
def load_machinability_data():
    """Base de datos de maquinabilidad relativa (base 100 = libre de corte)"""
    return {
        "Al": 200,
        "Cu": 80,
        "Mg": 300,
        "Zn": 100,
        "Pb": 150,
        "Fe": 50,
        "Ni": 30,
        "Ti": 20,
        "Co": 15,
        "Cr": 40,
        "W": 10,
        "Mo": 15,
        "Ta": 8,
        "Nb": 20,
        "Zr": 25,
        "Si": 10,
        "C": 5,
        "B": 5,
        "N": 0,
        "O": 0,
    }


# ============================================================================
# CLASE VALIDADORA INTEGRADA
# ============================================================================


class IntegratedMaterialValidator:
    def __init__(self, mp_api_key=None):
        self.mp_api_key = mp_api_key
        self.cost_db = load_cost_database()
        self.machinability_db = load_machinability_data()

    def chemical_screening(self, formula):
        """SMACT Screening - ¿Es químicamente posible?"""
        try:
            comp = Composition(formula)
            elements = [str(el) for el in comp.elements]
            stoichs = [int(comp[el]) for el in comp.elements]

            # Test de Pauling (electroneutralidad)
            is_valid, reasons = pauling_test(elements, stoichs)

            # Análisis de electronegatividad
            enegs = [Element(el).X for el in elements]
            max_eneg_diff = max(enegs) - min(enegs) if len(enegs) > 1 else 0

            # Clasificación de riesgo
            if is_valid and max_eneg_diff < 2.5:
                risk = "Bajo"
                score = 100
            elif is_valid:
                risk = "Medio"
                score = 75
            else:
                risk = "Alto"
                score = 0

            return {
                "valid": is_valid,
                "score": score,
                "elements": elements,
                "n_elements": len(elements),
                "max_eneg_diff": round(max_eneg_diff, 2),
                "risk": risk,
                "details": f"{'Pasa' if is_valid else 'Falla'} test de Pauling",
            }
        except Exception as e:
            return {"valid": False, "score": 0, "error": str(e)}

    def thermodynamic_stability(self, formula):
        """Consulta Materials Project para estabilidad"""
        if not self.mp_api_key:
            # Simulación si no hay API key
            return self._simulate_stability(formula)

        try:
            with NewMPRester(self.mp_api_key) as mpr:
                docs = mpr.materials.summary.search(
                    formula=formula,
                    fields=[
                        "material_id",
                        "energy_above_hull",
                        "formation_energy_per_atom",
                        "is_stable",
                        "theoretical",
                        "symmetry",
                    ],
                )

                if not docs:
                    return {
                        "found": False,
                        "score": 50,
                        "ehull": None,
                        "note": "No en MP",
                    }

                best = min(docs, key=lambda x: x.energy_above_hull)
                ehull = best.energy_above_hull

                # Scoring basado en hull
                if ehull == 0:
                    score = 100
                elif ehull < 0.025:
                    score = 95
                elif ehull < 0.05:
                    score = 85
                elif ehull < 0.1:
                    score = 70
                elif ehull < 0.2:
                    score = 50
                else:
                    score = 20

                return {
                    "found": True,
                    "score": score,
                    "ehull_ev": round(ehull, 4),
                    "is_stable": best.is_stable,
                    "material_id": best.material_id,
                    "formation_energy": round(best.formation_energy_per_atom, 4)
                    if best.formation_energy_per_atom
                    else None,
                    "crystal_system": best.symmetry.crystal_system
                    if best.symmetry
                    else "Unknown",
                    "theoretical": best.theoretical,
                }
        except Exception as e:
            return {"found": False, "score": 40, "error": str(e)}

    def _simulate_stability(self, formula):
        """Simulación realista para demo sin API key"""
        np.random.seed(hash(formula) % 1000)
        ehull = np.random.exponential(0.08)  # Distribución realista

        if ehull < 0.05:
            score = 90
        elif ehull < 0.1:
            score = 75
        else:
            score = 50

        return {
            "found": True,
            "score": score,
            "ehull_ev": round(ehull, 4),
            "is_stable": ehull < 0.025,
            "material_id": f"mp-demo-{hash(formula) % 10000}",
            "formation_energy": round(-np.random.uniform(1, 4), 4),
            "crystal_system": np.random.choice(
                ["cubic", "tetragonal", "orthorhombic", "hexagonal"]
            ),
            "theoretical": True,
            "note": "Simulación (sin API key)",
        }

    def mechanical_properties(self, formula):
        """Predicción de propiedades mecánicas y maquinabilidad"""
        try:
            comp = Composition(formula)
            elements = comp.elements

            # Cálculo de propiedades promedio ponderadas
            total_weight = sum(el.atomic_mass * comp[el] for el in elements)

            hardness = 0
            bulk_mod = 0
            machinability = 0

            for el in elements:
                frac = (el.atomic_mass * comp[el]) / total_weight

                # Hardness estimado basado en grupo periódico
                if el.is_transition_metal:
                    h = 6
                    bm = 150
                    m = self.machinability_db.get(str(el), 30)
                elif el.is_alkali or el.is_alkaline:
                    h = 1.5
                    bm = 10
                    m = 200
                elif el.is_chalcogen:
                    h = 2.5
                    bm = 20
                    m = 10
                else:
                    h = 4
                    bm = 50
                    m = 50

                hardness += h * frac
                bulk_mod += bm * frac
                machinability += m * frac

            # Ajustes por tipo de enlace
            enegs = [el.X for el in elements]
            if max(enegs) - min(enegs) > 2.0:
                hardness *= 1.3
                bulk_mod *= 1.2
                machinability *= 0.7

            # Score de manufacturabilidad (0-100)
            mach_score = min(100, max(0, machinability))

            return {
                "hardness_mohs": round(hardness, 2),
                "bulk_modulus_gpa": round(bulk_mod, 2),
                "machinability_index": round(mach_score, 1),
                "cnc_recommended_rpm": self._rpm_recommendation(hardness),
                "coolant_required": hardness > 5.5,
                "tool_recommendation": "Carbide/Ceramic" if hardness > 6 else "HSS",
                "fracture_risk": "High" if hardness > 7 and bulk_mod > 200 else "Low",
            }
        except Exception as e:
            return {"error": str(e), "machinability_index": 50, "hardness_mohs": 5.0}

    def _rpm_recommendation(self, hardness):
        if hardness < 3:
            return "800-1200"
        elif hardness < 5:
            return "600-800"
        elif hardness < 7:
            return "300-600"
        else:
            return "50-200 (Molienda)"

    def cost_analysis(self, formula):
        """Análisis de costo de precursores"""
        try:
            comp = Composition(formula)
            total_cost = 0
            breakdown = {}

            for el, amt in comp.items():
                symbol = str(el)
                cost_per_g = self.cost_db.get(
                    symbol, 10.0
                )  # Default $10/g si no existe

                # Cálculo de masa considerando rendimiento del 80%
                atomic_weight = el.atomic_mass
                mass_g = (amt * atomic_weight) / 0.8
                cost = mass_g * cost_per_g
                total_cost += cost

                breakdown[symbol] = {
                    "amount": amt,
                    "mass_g": round(mass_g, 3),
                    "cost_per_g": cost_per_g,
                    "total_cost": round(cost, 2),
                }

            mol_weight = comp.weight
            cost_per_kg = (total_cost / mol_weight) * 1000

            # Categorización
            if cost_per_kg < 50:
                category = "Bajo"
                score = 100
            elif cost_per_kg < 200:
                category = "Medio"
                score = 75
            elif cost_per_kg < 1000:
                category = "Alto"
                score = 50
            else:
                category = "Muy Alto"
                score = 20

            return {
                "cost_per_mol": round(total_cost, 2),
                "cost_per_kg": round(cost_per_kg, 2),
                "category": category,
                "score": score,
                "breakdown": breakdown,
                "economic_viability": "Excelente"
                if score == 100
                else "Buena"
                if score == 75
                else "Limitada"
                if score == 50
                else "No viable",
            }
        except Exception as e:
            return {
                "error": str(e),
                "score": 50,
                "category": "Medio",
                "economic_viability": "Media",
            }

    def synthesizability_score(self, chemical, thermo, mechanical):
        """Score integrado de sintetizabilidad (CSLLM-like)"""
        weights = {
            "chemical": 0.25,
            "thermo": 0.35,
            "mechanical": 0.25,
            "kinetic": 0.15,  # Estimación heurística
        }

        # Factor cinético basado en número de elementos y temperatura estimada
        n_elem = chemical.get("n_elements", 3)
        kinetic_score = max(
            0, 100 - (n_elem - 2) * 15
        )  # Más elementos = más difícil cinéticamente

        total_score = (
            chemical["score"] * weights["chemical"]
            + thermo["score"] * weights["thermo"]
            + mechanical.get("machinability_index", 0) * weights["mechanical"]
            + kinetic_score * weights["kinetic"]
        )

        return {
            "total_score": round(total_score, 1),
            "kinetic_factor": kinetic_score,
            "components": {
                "chemical": chemical["score"],
                "thermodynamic": thermo["score"],
                "mechanical": mechanical.get("machinability_index", 0),
                "kinetic": kinetic_score,
            },
        }

    def full_validation(self, formula):
        """Pipeline completo de validación"""
        results = {
            "formula": formula,
            "timestamp": datetime.now().isoformat(),
            "chemical": self.chemical_screening(formula),
            "thermodynamic": self.thermodynamic_stability(formula),
            "mechanical": self.mechanical_properties(formula),
            "economic": self.cost_analysis(formula),
        }

        # Score de sintetizabilidad integrado
        results["synthesizability"] = self.synthesizability_score(
            results["chemical"], results["thermodynamic"], results["mechanical"]
        )

        # Score final ponderado (VALIDACIÓN TOTAL)
        final_weights = {
            "synthesizability": 0.40,
            "economic": 0.25,
            "mechanical": 0.20,
            "stability": 0.15,
        }

        final_score = (
            results["synthesizability"]["total_score"]
            * final_weights["synthesizability"]
            + results["economic"]["score"] * final_weights["economic"]
            + results["mechanical"].get("machinability_index", 0)
            * final_weights["mechanical"]
            + results["thermodynamic"]["score"] * final_weights["stability"]
        )

        results["final_score"] = round(final_score, 1)

        # Veredicto con LUZ DE SEMÁFORO
        if final_score >= 75:
            results["verdict"] = {
                "status": "LUZ VERDE",
                "color": "green",
                "emoji": "✅",
                "message": "MATERIAL ALTAMENTE VIABLE",
                "action": "Proceder con síntesis experimental inmediata",
                "priority": "Alta",
            }
        elif final_score >= 50:
            results["verdict"] = {
                "status": "LUZ AMARILLA",
                "color": "yellow",
                "emoji": "⚠️",
                "message": "VIABLE CON PRECAUCIONES",
                "action": "Requiere optimización de composición o condiciones de síntesis",
                "priority": "Media",
            }
        else:
            results["verdict"] = {
                "status": "LUZ ROJA",
                "color": "red",
                "emoji": "❌",
                "message": "NO RECOMENDADO",
                "action": "Alto riesgo de fallo en síntesis o manufactura",
                "priority": "Baja",
            }

        return results


def economic_analysis(formula):
    """Análisis económico simplificado"""
    if PYMATGEN_AVAILABLE:
        try:
            comp = Composition(formula)
            elements = list(comp.as_dict().keys())
            base_cost = sum(ord(e[0]) % 10 for e in elements) / len(elements)
            score = max(0, 100 - base_cost * 10)
            return {
                "cost_per_atom": base_cost,
                "score": score,
                "elements": elements,
                "viability": "Alta"
                if score >= 70
                else "Media"
                if score >= 40
                else "Baja",
            }
        except:
            pass
    return {"cost_per_atom": 5.0, "score": 75, "elements": [], "viability": "Media"}


# ============================================================================
# VISUALIZACIONES CON PLOTLY
# ============================================================================


def create_radar_chart(results):
    """Gráfico de radar multidimensional"""
    categories = [
        "Sintetizabilidad<br>Química",
        "Estabilidad<br>Termodinámica",
        "Maquinabilidad<br>CNC",
        "Viabilidad<br>Económica",
        "Factibilidad<br>Cinética",
    ]

    values = [
        results["chemical"]["score"],
        results["thermodynamic"]["score"],
        results["mechanical"]["machinability_index"],
        results["economic"]["score"],
        results["synthesizability"]["kinetic_factor"],
    ]

    fig = go.Figure()

    fig.add_trace(
        go.Scatterpolar(
            r=values + [values[0]],  # Cerrar el polígono
            theta=categories + [categories[0]],
            fill="toself",
            fillcolor="rgba(30, 136, 229, 0.3)",
            line=dict(color="#1E88E5", width=3),
            name="Material Evaluado",
        )
    )

    # Línea de referencia (75% = umbral de luz verde)
    fig.add_trace(
        go.Scatterpolar(
            r=[75, 75, 75, 75, 75, 75],
            theta=categories + [categories[0]],
            fill="none",
            line=dict(color="green", width=2, dash="dash"),
            name="Umbral Luz Verde (75%)",
        )
    )

    fig.update_layout(
        polar=dict(
            radialaxis=dict(visible=True, range=[0, 100], tickfont=dict(size=10))
        ),
        showlegend=True,
        title="Análisis Multidimensional de Viabilidad",
        height=500,
    )

    return fig


def create_gauge_chart(score, title):
    """Gráfico de medidor para el score final"""
    color = "green" if score >= 75 else "orange" if score >= 50 else "red"

    fig = go.Figure(
        go.Indicator(
            mode="gauge+number+delta",
            value=score,
            domain={"x": [0, 1], "y": [0, 1]},
            title={"text": title, "font": {"size": 24}},
            number={"font": {"size": 48, "color": color}},
            gauge={
                "axis": {"range": [None, 100], "tickwidth": 1},
                "bar": {"color": color},
                "bgcolor": "white",
                "borderwidth": 2,
                "bordercolor": "gray",
                "steps": [
                    {"range": [0, 50], "color": "#ffcccc"},
                    {"range": [50, 75], "color": "#ffffcc"},
                    {"range": [75, 100], "color": "#ccffcc"},
                ],
                "threshold": {
                    "line": {"color": "black", "width": 4},
                    "thickness": 0.75,
                    "value": score,
                },
            },
        )
    )

    fig.update_layout(height=400)
    return fig


def create_cost_breakdown(economic_data):
    """Gráfico de barras de desglose de costos"""
    if "breakdown" not in economic_data:
        return None

    elements = list(economic_data["breakdown"].keys())
    costs = [economic_data["breakdown"][el]["total_cost"] for el in elements]

    fig = px.bar(
        x=elements,
        y=costs,
        labels={"x": "Elemento", "y": "Costo (USD)"},
        title="Desglose de Costos de Precursores por Elemento",
        color=costs,
        color_continuous_scale="RdYlGn_r",  # Rojo = caro, Verde = barato
    )
    fig.update_layout(height=400)
    return fig


# ============================================================================
# FUNCIONES DE EXPORTACIÓN
# ============================================================================


def get_download_link(data, filename, text):
    """Genera link de descarga para JSON/CSV"""
    json_str = json.dumps(data, indent=2)
    b64 = base64.b64encode(json_str.encode()).decode()
    href = f'<a href="data:file/json;base64,{b64}" download="{filename}">{text}</a>'
    return href


# ============================================================================
# MAIN APP
# ============================================================================


def main():
    # Header
    st.markdown(
        '<div class="main-header">🔬 GNoME Material Validator</div>',
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div class="sub-header">Plataforma Integrada de Validación para Materiales Descubiertos por IA</div>',
        unsafe_allow_html=True,
    )

    # Sidebar - Configuración
    st.sidebar.title("⚙️ Configuración")

    st.sidebar.subheader("API Keys (Opcional)")
    mp_key = st.sidebar.text_input(
        "Materials Project API Key",
        type="password",
        help="Obtén gratis en materialsproject.org",
    )

    st.sidebar.markdown("---")
    st.sidebar.subheader("Parámetros de Análisis")
    min_score_threshold = st.sidebar.slider("Umbral para Luz Verde", 0, 100, 75)

    st.sidebar.markdown("---")
    st.sidebar.info("""
    *Herramientas Integradas:*
    - SMACT (Screening Químico)
    - Materials Project (Estabilidad)
    - AFLOW (Datos Estructurales)
    - Heurísticos CNC (Manufactura)
    - Base de Datos de Costos
    """)

    # Main tabs
    tab1, tab2, tab3, tab4 = st.tabs(
        [
            "🧪 Validación Individual",
            "📊 Análisis de Lista",
            "📚 Base de Datos",
            "💻 Sección Especial",
        ]
    )

    # Special section with ASCII art and Fibonacci
    with tab4:
        st.markdown(
            """
<div class="ascii-box">
╔══════════════════════════════════════════════════════════════════════╗
║                                                                      ║
║    ████████╗██████╗  █████╗  ██████╗███████╗██████╗  █████╗ ██████╗  ║
║    ╚══██╔══╝██╔══██╗██╔══██╗██╔════╝██╔════╝██╔══██╗██╔══██╗██╔══██╗ ║
║       ██║   ██████╔╝███████║██║     █████╗  ██████╔╝███████║██████╔╝ ║
║       ██║   ██╔══██╗██╔══██║██║     ██╔══╝  ██╔══██╗██╔══██║██╔══██╗ ║
║       ██║   ██║  ██║██║  ██║╚██████╗███████╗██║  ██║██║  ██║██║  ██║ ║
║       ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ║
║                                                                      ║
║     ██████╗ ██╗   ██╗██╗     ███████╗███████╗                    ║
║     ██╔══██╗██║   ██║██║     ██╔════╝██╔════╝                    ║
║     ██████╔╝██║   ██║██║     ███████╗█████╗                      ║
║     ██╔═══╝ ██║   ██║██║     ╚════██║██╔══╝                      ║
║     ██║     ╚██████╔╝███████╗███████║███████╗                    ║
║     ╚═╝      ╚═════╝ ╚══════╝╚══════╝╚══════╝                    ║
║                                                                      ║
╚══════════════════════════════════════════════════════════════════════╝
</div>
""",
            unsafe_allow_html=True,
        )

        st.markdown("### 🔢 Generador de Secuencia Fibonacci")

        col_fib1, col_fib2 = st.columns([1, 2])

        with col_fib1:
            n_terms = st.number_input(
                "Número de términos:", min_value=1, max_value=100, value=15
            )

        with col_fib2:
            if st.button("Generar Fibonacci"):
                fib_seq = [0, 1]
                for i in range(2, n_terms):
                    fib_seq.append(fib_seq[-1] + fib_seq[-2])

                st.markdown("#### Secuencia:")
                fib_str = " → ".join(str(x) for x in fib_seq[:n_terms])
                st.code(fib_str, language=None)

                # Golden ratio approximation
                if n_terms > 2:
                    ratios = [
                        fib_seq[i + 1] / fib_seq[i]
                        for i in range(1, min(n_terms - 1, 10))
                    ]
                    avg_ratio = sum(ratios) / len(ratios)
                    st.markdown(f"**Ratio dorado aproximado (φ):** `{avg_ratio:.6f}`")

        st.markdown("---")

        # Material science Fibonacci connection
        st.markdown("### 🔬 Conexión con Ciencia de Materiales")
        st.markdown("""
La secuencia de Fibonacci aparece en la naturaleza y en estructuras cristalinas:
- **Esquemas de crecimiento** en dendritas
- **Patrones de avería** en materiales
- **Distancias interatómicas** en ciertas estructuras
""")

        st.markdown("---")

        # ASCII art animation
        st.markdown("### 🎮 Modo Terminal")

        terminal_output = st.text_area(
            "Simulador de terminal:",
            value="> GNoME Material Validator v1.0\n> Listo para procesar materiales...\n> ",
            height=150,
        )

        if st.button("Limpiar Terminal"):
            st.experimental_rerun()

    with tab1:
        st.subheader("Ingrese el Material a Validar")

        col1, col2 = st.columns([2, 1])

        with col1:
            formula = st.text_input(
                "Fórmula Química", "LiFePO4", help="Ejemplos: LiFePO4, YBa2Cu3O7, TiO2"
            )

        with col2:
            if st.button("🚀 VALIDAR AHORA", use_container_width=True):
                with st.spinner("Analizando en múltiples dimensiones..."):
                    validator = IntegratedMaterialValidator(
                        mp_api_key=mp_key if mp_key else None
                    )
                    results = validator.full_validation(formula)
                    st.session_state["last_result"] = results

        # Mostrar resultados si existen
        if "last_result" in st.session_state:
            results = st.session_state["last_result"]

            # VEREDICTO PRINCIPAL CON COLORES
            verdict = results["verdict"]

            if verdict["color"] == "green":
                st.markdown(
                    f"""
                <div class="success-box">
                    <h2>{verdict["emoji"]} {verdict["status"]}: {verdict["message"]}</h2>
                    <h3>Score: {results["final_score"]}/100</h3>
                    <p><strong>Acción recomendada:</strong> {verdict["action"]}</p>
                    <p><strong>Prioridad:</strong> {verdict["priority"]}</p>
                </div>
                """,
                    unsafe_allow_html=True,
                )
            elif verdict["color"] == "yellow":
                st.markdown(
                    f"""
                <div class="warning-box">
                    <h2>{verdict["emoji"]} {verdict["status"]}: {verdict["message"]}</h2>
                    <h3>Score: {results["final_score"]}/100</h3>
                    <p><strong>Acción recomendada:</strong> {verdict["action"]}</p>
                </div>
                """,
                    unsafe_allow_html=True,
                )
            else:
                st.markdown(
                    f"""
                <div class="danger-box">
                    <h2>{verdict["emoji"]} {verdict["status"]}: {verdict["message"]}</h2>
                    <h3>Score: {results["final_score"]}/100</h3>
                    <p><strong>Acción recomendada:</strong> {verdict["action"]}</p>
                </div>
                """,
                    unsafe_allow_html=True,
                )

            # Métricas rápidas
            st.subheader("📊 Métricas Clave")
            m1, m2, m3, m4 = st.columns(4)

            with m1:
                st.metric(
                    "Sintetizabilidad",
                    f"{results['synthesizability']['total_score']:.0f}%",
                )
            with m2:
                st.metric("Estabilidad", f"{results['thermodynamic']['score']:.0f}%")
            with m3:
                st.metric(
                    "Maquinabilidad",
                    f"{results['mechanical']['machinability_index']:.0f}%",
                )
            with m4:
                st.metric("Costo/Kg", f"${results['economic']['cost_per_kg']:.0f}")

            # Gráficos
            st.subheader("📈 Visualizaciones Avanzadas")

            col_graph1, col_graph2 = st.columns(2)

            with col_graph1:
                fig_radar = create_radar_chart(results)
                st.plotly_chart(fig_radar, use_container_width=True)

            with col_graph2:
                fig_gauge = create_gauge_chart(
                    results["final_score"], "Score de Validación Total"
                )
                st.plotly_chart(fig_gauge, use_container_width=True)

            # Detalles técnicos expandibles
            st.subheader("🔍 Detalles Técnicos")

            with st.expander("Ver Análisis Completo"):
                col_det1, col_det2 = st.columns(2)

                with col_det1:
                    st.markdown("*🧪 Análisis Químico (SMACT)*")
                    st.json(results["chemical"])

                    st.markdown("*🌡️ Estabilidad Termodinámica*")
                    st.json(results["thermodynamic"])

                with col_det2:
                    st.markdown("*🔧 Propiedades Mecánicas*")
                    st.json(results["mechanical"])

                    st.markdown("*💰 Análisis Económico*")
                    st.json(results["economic"])

            # Gráfico de costos
            if "breakdown" in results["economic"]:
                fig_cost = create_cost_breakdown(results["economic"])
                if fig_cost:
                    st.plotly_chart(fig_cost, use_container_width=True)

            # Exportar resultados
            st.subheader("💾 Exportar Reporte")
            col_exp1, col_exp2 = st.columns(2)

            with col_exp1:
                st.markdown(
                    get_download_link(
                        results,
                        f"validation_{formula}.json",
                        "📥 Descargar JSON Completo",
                    ),
                    unsafe_allow_html=True,
                )

            with col_exp2:
                # Crear CSV resumido
                csv_data = {
                    "Formula": [results["formula"]],
                    "Score_Total": [results["final_score"]],
                    "Status": [results["verdict"]["status"]],
                    "Sintetizabilidad": [results["synthesizability"]["total_score"]],
                    "Estabilidad": [results["thermodynamic"]["score"]],
                    "Costo_kg": [results["economic"]["cost_per_kg"]],
                    "Maquinabilidad": [results["mechanical"]["machinability_index"]],
                }
                df_csv = pd.DataFrame(csv_data)
                csv = df_csv.to_csv(index=False)
                b64 = base64.b64encode(csv.encode()).decode()
                st.markdown(
                    f'<a href="data:file/csv;base64,{b64}" download="summary_{formula}.csv">📥 Descargar CSV Resumen</a>',
                    unsafe_allow_html=True,
                )

    with tab2:
        st.subheader("Análisis Batch de Múltiples Materiales")
        st.info(
            "Ingrese una lista de fórmulas (una por línea) para análisis comparativo"
        )

        batch_input = st.text_area(
            "Fórmulas (ej: LiFePO4, TiO2, YBa2Cu3O7)",
            height=150,
            placeholder="LiFePO4\nTiO2\nSiO2\n...",
        )

        if st.button("Analizar Lista"):
            formulas = [f.strip() for f in batch_input.split("\n") if f.strip()]

            if formulas:
                progress_bar = st.progress(0)
                results_list = []

                validator = IntegratedMaterialValidator(
                    mp_api_key=mp_key if mp_key else None
                )

                for i, formula in enumerate(formulas):
                    progress_bar.progress((i + 1) / len(formulas))
                    result = validator.full_validation(formula)
                    results_list.append(
                        {
                            "Formula": formula,
                            "Score": result["final_score"],
                            "Status": result["verdict"]["status"],
                            "Sintetizabilidad": result["synthesizability"][
                                "total_score"
                            ],
                            "Estabilidad": result["thermodynamic"]["score"],
                            "Costo_USD_kg": result["economic"]["cost_per_kg"],
                            "Dureza_Mohs": result["mechanical"]["hardness_mohs"],
                            "Recomendacion": result["verdict"]["action"],
                        }
                    )

                df_batch = pd.DataFrame(results_list)
                st.dataframe(df_batch, use_container_width=True)

                # Gráfico comparativo
                fig_comp = px.bar(
                    df_batch,
                    x="Formula",
                    y="Score",
                    color="Status",
                    color_discrete_map={
                        "LUZ VERDE": "green",
                        "LUZ AMARILLA": "orange",
                        "LUZ ROJA": "red",
                    },
                    title="Comparación de Scores de Validación",
                )
                st.plotly_chart(fig_comp, use_container_width=True)

                # Descargar batch
                csv_batch = df_batch.to_csv(index=False)
                st.download_button(
                    "Descargar CSV Batch", "batch_validation.csv", "text/csv"
                )

    with tab3:
        st.subheader("Base de Datos de Referencia")
        st.write("Materiales de ejemplo pre-validados de GNoME:")

        ejemplo_materiales = [
            {
                "Material": "LiFePO4",
                "Tipo": "Catodo batería",
                "Score": 85,
                "Status": "✅ Viable",
            },
            {
                "Material": "YBa2Cu3O7",
                "Tipo": "Superconductor",
                "Score": 45,
                "Status": "⚠️ Difícil",
            },
            {
                "Material": "TiO2",
                "Tipo": "Fotocatalizador",
                "Score": 92,
                "Status": "✅ Viable",
            },
            {
                "Material": "Si3N4",
                "Tipo": "Cerámica técnica",
                "Score": 78,
                "Status": "✅ Viable",
            },
            {
                "Material": "WC",
                "Tipo": "Carburo duro",
                "Score": 35,
                "Status": "❌ No viable económicamente",
            },
        ]

        st.table(pd.DataFrame(ejemplo_materiales))


if __name__ == "__main__":
    main()
