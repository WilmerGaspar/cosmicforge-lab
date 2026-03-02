"""
CosmicForge Lab - Personal Edition
Diseño de Materiales Inspirado en Firmas Astrofísicas
Versión 3.1 - Edición Personal Mejorada

Mejoras:
- Letras legibles y claras
- PDF con formato correcto
- Proceso de producción detallado
- Guía de síntesis paso a paso
"""

import streamlit as st
import json
import sys
import os
import numpy as np
from pathlib import Path
from datetime import datetime
import io
import base64

# Add modules directory to path
MODULES_PATH = Path(__file__).parent / "modules"
sys.path.insert(0, str(MODULES_PATH))

from modules.physics_calculator import PhysicsCalculator
from modules.chemistry_engine import ChemistryEngine
from modules.file_generators import FileGenerator
from modules.visualizer import Visualizer
from modules.file_parser import FileParser

# Plotly para 3D
try:
    import plotly.graph_objects as go
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except:
    PLOTLY_AVAILABLE = False

# PDF generation con FPDF (más confiable)
try:
    from fpdf import FPDF
    PDF_AVAILABLE = True
except:
    try:
        from reportlab.lib.pagesizes import letter, A4
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
        from reportlab.lib import colors
        PDF_AVAILABLE = True
        USE_REPORTLAB = True
    except:
        PDF_AVAILABLE = False
        USE_REPORTLAB = False

# ============================================================================
# PAGE CONFIGURATION
# ============================================================================

st.set_page_config(
    page_title="CosmicForge Lab",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================================
# SESSION STATE
# ============================================================================

if 'astro_data' not in st.session_state:
    st.session_state.astro_data = None
if 'physical_props' not in st.session_state:
    st.session_state.physical_props = None
if 'recipe' not in st.session_state:
    st.session_state.recipe = None
if 'raw_json' not in st.session_state:
    st.session_state.raw_json = None

# ============================================================================
# DATABASES
# ============================================================================

# Ejemplos de objetos astrofísicos
ASTROPHYSICAL_EXAMPLES = {
    "Orion Nebula M42": {
        "fractal_dimension": 1.654, "criticality_score": 0.722, "entropy": 0.019,
        "anisotropy": 0.329, "turbulence_beta": 2.278, "lyapunov_max": -0.227,
        "mode": "balanced", "type": "Nebulosa de Emisión", "distance": "1,344 ly"
    },
    "Crab Nebula M1": {
        "fractal_dimension": 1.82, "criticality_score": 0.89, "entropy": 0.032,
        "anisotropy": 0.456, "turbulence_beta": 2.45, "lyapunov_max": -0.15,
        "mode": "turbulent", "type": "Remanente Supernova", "distance": "6,500 ly"
    },
    "Ring Nebula M57": {
        "fractal_dimension": 1.45, "criticality_score": 0.55, "entropy": 0.012,
        "anisotropy": 0.28, "turbulence_beta": 1.95, "lyapunov_max": -0.35,
        "mode": "stable", "type": "Nebulosa Planetaria", "distance": "2,283 ly"
    },
    "Triangulum Galaxy M33": {
        "fractal_dimension": 1.93, "criticality_score": 0.78, "entropy": 0.00000964,
        "anisotropy": 0.23, "turbulence_beta": 1.84, "lyapunov_max": 0.11,
        "mode": "balanced", "type": "Galaxia Espiral", "distance": "2.73M ly"
    },
    "Andromeda Galaxy M31": {
        "fractal_dimension": 1.88, "criticality_score": 0.85, "entropy": 0.000012,
        "anisotropy": 0.19, "turbulence_beta": 1.92, "lyapunov_max": 0.08,
        "mode": "balanced", "type": "Galaxia Espiral", "distance": "2.537M ly"
    },
    "Eagle Nebula M16": {
        "fractal_dimension": 1.78, "criticality_score": 0.81, "entropy": 0.028,
        "anisotropy": 0.42, "turbulence_beta": 2.35, "lyapunov_max": -0.18,
        "mode": "turbulent", "type": "Nebulosa de Emisión", "distance": "7,000 ly"
    },
    "Tarantula Nebula": {
        "fractal_dimension": 1.91, "criticality_score": 0.93, "entropy": 0.045,
        "anisotropy": 0.58, "turbulence_beta": 2.68, "lyapunov_max": 0.15,
        "mode": "turbulent", "type": "Región HII", "distance": "160,000 ly"
    },
    "Pillars of Creation": {
        "fractal_dimension": 1.85, "criticality_score": 0.88, "entropy": 0.035,
        "anisotropy": 0.52, "turbulence_beta": 2.55, "lyapunov_max": -0.12,
        "mode": "turbulent", "type": "Nube Molecular", "distance": "6,500 ly"
    },
    "Whirlpool Galaxy M51": {
        "fractal_dimension": 1.95, "criticality_score": 0.91, "entropy": 0.000015,
        "anisotropy": 0.22, "turbulence_beta": 1.98, "lyapunov_max": 0.12,
        "mode": "turbulent", "type": "Galaxia Interactuante", "distance": "23M ly"
    },
    "Helix Nebula NGC 7293": {
        "fractal_dimension": 1.48, "criticality_score": 0.58, "entropy": 0.011,
        "anisotropy": 0.31, "turbulence_beta": 2.02, "lyapunov_max": -0.32,
        "mode": "stable", "type": "Nebulosa Planetaria", "distance": "655 ly"
    },
}

# Sistema de aleaciones completo
ALLOY_SYSTEMS = {
    "Ti-Al": {
        "name": "Titanio-Aluminio (Gamma-TiAl)",
        "ratio": "Ti:Al = 1:1",
        "melting_point": "1460°C",
        "density": "3.76 g/cm³",
        "applications": ["Turbinas aeronáuticas", "Componentes automotrices", "Implantes médicos"],
        "synthesis": [
            "1. Fundir Ti puro (99.9%) en horno de arco eléctrico bajo argón",
            "2. Añadir Al puro (99.99%) en proporción 1:1 atómica",
            "3. Homogeneizar a 1400°C por 4 horas",
            "4. Enfriar lentamente (5°C/min) hasta 1000°C",
            "5. Tratamiento térmico a 900°C por 2 horas",
            "6. Enfriar al aire"
        ],
        "equipment": ["Horno de arco eléctrico", "Atmósfera de argón", "Crisol de grafito", "Horno de tratamiento térmico"],
        "precursors": ["Ti esponja (99.9%)", "Al lingote (99.99%)"],
        "safety": ["Usar guantes refractarios", "Atmósfera inerte obligatoria", "Ventilación adecuada"]
    },
    "Ti-Al-V (TA6V)": {
        "name": "Titanio TA6V (Ti-6Al-4V)",
        "ratio": "Ti:Al:V = 90:6:4",
        "melting_point": "1650°C",
        "density": "4.43 g/cm³",
        "applications": ["Aeroespacial", "Implantes ortopédicos", "Industria química"],
        "synthesis": [
            "1. Preparar aleación maestra Ti-Al (60:40) previamente",
            "2. Fundir Ti esponja en horno de arco bajo vacío",
            "3. Añadir aleación maestra Ti-Al y V metálico",
            "4. Fundir y voltear 3-5 veces para homogeneizar",
            "5. Colar en molde de grafito precalentado",
            "6. Forjar a 950°C para eliminar segregación",
            "7. Tratamiento térmico: solución a 950°C + envejecimiento a 540°C"
        ],
        "equipment": ["Horno de arco bajo vacío", "Crisol de cobre refrigerado", "Prensa de forja", "Horno de tratamiento"],
        "precursors": ["Ti esponja grado aeronáutico", "Al lingote", "V metálico"],
        "safety": ["Vacío alto (<10⁻³ mbar)", "Protección radiación UV", "Evitar contaminación con O₂, N₂, H₂"]
    },
    "Fe-Ni (Invar)": {
        "name": "Invar (Fe-36Ni)",
        "ratio": "Fe:Ni = 64:36",
        "melting_point": "1425°C",
        "density": "8.12 g/cm³",
        "applications": ["Instrumentos de precisión", "Relojes", "Sellos térmicos"],
        "synthesis": [
            "1. Fundir Fe electrolítico en horno de inducción",
            "2. Añadir Ni electrolítico gradualmente",
            "3. Mantener a 1500°C por 30 min bajo argón",
            "4. Colar en lingotera precalentada",
            "5. Laminar en caliente a 1000°C",
            "6. Recocer a 850°C por 1 hora",
            "7. Enfriar lentamente en horno"
        ],
        "equipment": ["Horno de inducción", "Atmósfera inerte", "Laminador", "Horno de recocido"],
        "precursors": ["Fe electrolítico", "Ni electrolítico (99.9%)"],
        "safety": ["Evitar oxidación", "Controlar temperatura exacta"]
    },
    "Cu-Ni (Cuproníquel)": {
        "name": "Cuproníquel 70/30",
        "ratio": "Cu:Ni = 70:30",
        "melting_point": "1170°C",
        "density": "8.94 g/cm³",
        "applications": ["Componentes marinos", "Monedas", "Intercambiadores de calor"],
        "synthesis": [
            "1. Fundir Cu electrolítico en horno de inducción",
            "2. Añadir Ni gradualmente bajo fundente bórax",
            "3. Desoxidar con Cu-P (0.02%)",
            "4. Colar a 1250°C en molde metálico",
            "5. Laminar en caliente a 800°C",
            "6. Recocer a 650°C por 2 horas"
        ],
        "equipment": ["Horno de inducción", "Fundente bórax", "Laminador"],
        "precursors": ["Cu cátodo (99.99%)", "Ni electrolítico"],
        "safety": ["Fundente para evitar oxidación", "Ventilación"]
    },
    "Al-Cu (Duraluminio)": {
        "name": "Duraluminio 2024",
        "ratio": "Al:Cu:Mg:Mn = 93.5:4.4:1.5:0.6",
        "melting_point": "660°C",
        "density": "2.78 g/cm³",
        "applications": ["Estructuras aeronáuticas", "Componentes automotrices", "Herramientas"],
        "synthesis": [
            "1. Fundir Al primario en crisol de grafito",
            "2. Añadir Cu, Mg y Mn como aleaciones maestras",
            "3. Degasear con cloro o nitrógeno",
            "4. Refinar con fundente salino",
            "5. Colar a 700°C en molde metálico",
            "6. Tratamiento T6: solución 495°C + temple + envejecimiento 190°C"
        ],
        "equipment": ["Horno de resistencia", "Crisol grafito", "Molde metálico", "Horno de tratamiento"],
        "precursors": ["Al primario", "Cu electrolítico", "Mg lingote", "Mn metálico"],
        "safety": ["Evitar humedad (explosión)", "Degaseado obligatorio"]
    },
    "TiO2 (Titania)": {
        "name": "Óxido de Titanio (TiO₂)",
        "ratio": "Ti:O = 1:2",
        "melting_point": "1843°C",
        "density": "4.23 g/cm³",
        "applications": ["Pigmentos", "Fotocatálisis", "Celdas solares", "Recubrimientos"],
        "synthesis": [
            "MÉTODO SOL-GEL:",
            "1. Disolver precursor Ti (isopropóxido de Ti) en etanol",
            "2. Añadir agua destilada lentamente con agitación",
            "3. Ajustar pH a 2 con HNO₃",
            "4. Envejecer gel 24 horas",
            "5. Secar a 100°C por 12 horas",
            "6. Calcinación: 400°C (anatasa) o 800°C (rutilo)",
            "",
            "MÉTODO HIDROTERMAL:",
            "1. Mezclar TiCl₄ con NaOH 10M",
            "2. Transferir a autoclave",
            "3. Calentar a 180°C por 12 horas",
            "4. Lavar con agua hasta pH neutro",
            "5. Secar a 80°C"
        ],
        "equipment": ["Matraz reacción", "Agitador magnético", "Horno mufla", "Autoclave (opcional)"],
        "precursors": ["Isopropóxido de Ti", "TiCl₄", "Etanol", "HNO₃"],
        "safety": ["Guantes y gafas", "Campana extractora", "TiCl₄ es corrosivo"]
    },
    "ZnO (Óxido de Zinc)": {
        "name": "Óxido de Zinc (ZnO)",
        "ratio": "Zn:O = 1:1",
        "melting_point": "1975°C",
        "density": "5.61 g/cm³",
        "applications": ["Protectores solares", "Sensores de gas", "Varistores", "Catalizadores"],
        "synthesis": [
            "MÉTODO DE PRECIPITACIÓN:",
            "1. Disolver Zn(NO₃)₂ en agua desionizada (0.5M)",
            "2. Preparar solución NaOH 1M",
            "3. Añadir NaOH gota a gota con agitación vigorosa",
            "4. Mantener pH 10-11",
            "5. Envejecer precipitado 2 horas",
            "6. Filtrar y lavar con agua/etanol",
            "7. Secar a 80°C por 12 horas",
            "8. Calcinación: 400-600°C por 4 horas",
            "",
            "MÉTODO SOLVOTÉRMICO:",
            "1. Disolver acetato de Zn en metanol",
            "2. Añadir NaOH en metanol",
            "3. Calentar a 60°C por 2 horas",
            "4. Centrifugar y lavar",
            "5. Secar a 60°C"
        ],
        "equipment": ["Matraz", "Agitador", "Horno mufla", "Centrífuga"],
        "precursors": ["Zn(NO₃)₂", "NaOH", "Acetato de Zn", "Metanol"],
        "safety": ["Guantes", "Evitar inhalación de polvo"]
    },
    "Fe₂O₃ (Hematita)": {
        "name": "Óxido de Hierro (Hematita α-Fe₂O₃)",
        "ratio": "Fe:O = 2:3",
        "melting_point": "1565°C",
        "density": "5.26 g/cm³",
        "applications": ["Pigmentos", "Catálisis", "Sensores magnéticos", "Electrodos"],
        "synthesis": [
            "MÉTODO DE PRECIPITACIÓN:",
            "1. Disolver FeCl₃·6H₂O en agua (0.2M)",
            "2. Añadir NH₄OH hasta pH 8-9",
            "3. Envejecer precipitado 24 horas",
            "4. Filtrar y lavar con agua caliente",
            "5. Secar a 100°C",
            "6. Calcinación: 500-700°C por 3 horas",
            "",
            "MÉTODO HIDROTERMAL:",
            "1. Mezclar FeCl₃ y NaOH (relación 1:3)",
            "2. Autoclave a 180°C por 12 horas",
            "3. Lavar hasta pH neutro",
            "4. Secar a 80°C"
        ],
        "equipment": ["Matraz", "Agitador", "Horno mufla", "Autoclave"],
        "precursors": ["FeCl₃·6H₂O", "Fe(NO₃)₃", "NH₄OH", "NaOH"],
        "safety": ["Manchas la piel", "Usar guantes"]
    },
}

# Elementos y compatibilidades
ELEMENTS_DATABASE = {
    "Ti": {"name": "Titanio", "atomic_num": 22, "mw": 47.867, "mp": 1668, "compatible": ["Al", "V", "Fe", "Ni", "O", "N", "C"]},
    "Al": {"name": "Aluminio", "atomic_num": 13, "mw": 26.982, "mp": 660, "compatible": ["Ti", "Cu", "Mg", "Si", "Zn", "O"]},
    "Fe": {"name": "Hierro", "atomic_num": 26, "mw": 55.845, "mp": 1538, "compatible": ["Ni", "Cr", "Co", "Mn", "C", "O"]},
    "Zn": {"name": "Zinc", "atomic_num": 30, "mw": 65.38, "mp": 420, "compatible": ["Cu", "Al", "O", "S"]},
    "Cu": {"name": "Cobre", "atomic_num": 29, "mw": 63.546, "mp": 1085, "compatible": ["Ni", "Zn", "Al", "Sn", "O"]},
    "Ni": {"name": "Níquel", "atomic_num": 28, "mw": 58.693, "mp": 1455, "compatible": ["Fe", "Co", "Cu", "Cr", "Ti", "O"]},
    "Co": {"name": "Cobalto", "atomic_num": 27, "mw": 58.933, "mp": 1495, "compatible": ["Ni", "Fe", "Mn", "O"]},
    "Mn": {"name": "Manganeso", "atomic_num": 25, "mw": 54.938, "mp": 1246, "compatible": ["Fe", "Co", "O"]},
    "Ag": {"name": "Plata", "atomic_num": 47, "mw": 107.868, "mp": 962, "compatible": ["Cu", "Au", "Pd", "O"]},
    "Au": {"name": "Oro", "atomic_num": 79, "mw": 196.967, "mp": 1064, "compatible": ["Ag", "Pt", "Pd", "Cu"]},
    "Pt": {"name": "Platino", "atomic_num": 78, "mw": 195.084, "mp": 1768, "compatible": ["Pd", "Au", "Rh", "Ir", "O"]},
    "Pd": {"name": "Paladio", "atomic_num": 46, "mw": 106.42, "mp": 1555, "compatible": ["Pt", "Au", "Ag", "Ni", "O"]},
}

# Citas de Manin
MANIN_QUOTES = [
    "La matemática es la parte más pura de la cultura humana, y sin embargo está profundamente conectada con la realidad física.",
    "La demostrabilidad es un concepto que trasciende la lógica formal: es la búsqueda de certeza en un universo de incertidumbre.",
    "Entre lo demostrable y lo indemostrable existe una zona fronteriza donde nace la creatividad científica.",
    "La física teórica y las matemáticas comparten un misterio común: ¿por qué sus abstracciones describen tan bien la realidad?",
    "La verdad matemática tiene un estatus ontológico único: existe independientemente de nuestra capacidad para demostrarla.",
]

# ============================================================================
# CUSTOM CSS - LETRAS CLARAS Y LEGIBLES
# ============================================================================

st.markdown("""
<style>
    /* Fuente base */
    html, body, [class*="css"] {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif !important;
        font-size: 16px !important;
        line-height: 1.6 !important;
    }
    
    /* Títulos principales */
    h1 {
        font-size: 2.5rem !important;
        font-weight: 700 !important;
        color: #1565C0 !important;
        text-align: center !important;
        margin-bottom: 0.5rem !important;
    }
    
    h2 {
        font-size: 1.8rem !important;
        font-weight: 600 !important;
        color: #1976D2 !important;
        margin-top: 1.5rem !important;
    }
    
    h3 {
        font-size: 1.4rem !important;
        font-weight: 600 !important;
        color: #1E88E5 !important;
    }
    
    h4 {
        font-size: 1.2rem !important;
        font-weight: 600 !important;
        color: #42A5F5 !important;
    }
    
    /* Párrafos y texto */
    p, div, span, label {
        font-size: 15px !important;
        color: #333333 !important;
    }
    
    /* Métricas */
    [data-testid="stMetricValue"] {
        font-size: 1.8rem !important;
        font-weight: 700 !important;
        color: #1565C0 !important;
    }
    
    [data-testid="stMetricLabel"] {
        font-size: 0.95rem !important;
        font-weight: 500 !important;
        color: #666666 !important;
    }
    
    /* Cards */
    .info-card {
        background: linear-gradient(135deg, #E3F2FD 0%, #BBDEFB 100%);
        padding: 1.5rem;
        border-radius: 12px;
        border-left: 5px solid #1976D2;
        margin: 1rem 0;
    }
    
    .success-card {
        background: linear-gradient(135deg, #E8F5E9 0%, #C8E6C9 100%);
        padding: 1.5rem;
        border-radius: 12px;
        border-left: 5px solid #43A047;
        margin: 1rem 0;
    }
    
    .warning-card {
        background: linear-gradient(135deg, #FFF8E1 0%, #FFECB3 100%);
        padding: 1.5rem;
        border-radius: 12px;
        border-left: 5px solid #FFA000;
        margin: 1rem 0;
    }
    
    .quote-card {
        background: linear-gradient(135deg, #FCE4EC 0%, #F8BBD0 100%);
        padding: 1.5rem;
        border-radius: 12px;
        border-left: 5px solid #E91E63;
        margin: 1rem 0;
        font-style: italic;
    }
    
    .production-card {
        background: #FAFAFA;
        padding: 1.5rem;
        border-radius: 12px;
        border: 2px solid #1976D2;
        margin: 1rem 0;
    }
    
    .step-box {
        background: #F5F5F5;
        padding: 1rem;
        border-radius: 8px;
        margin: 0.5rem 0;
        border-left: 4px solid #42A5F5;
    }
    
    /* Tablas */
    table {
        font-size: 14px !important;
    }
    
    th {
        background-color: #1976D2 !important;
        color: white !important;
        font-weight: 600 !important;
    }
    
    td {
        padding: 10px !important;
    }
    
    /* Botones */
    .stButton button {
        font-size: 16px !important;
        font-weight: 600 !important;
        padding: 12px 24px !important;
    }
    
    /* Selectbox y inputs */
    .stSelectbox label, .stTextInput label, .stNumberInput label {
        font-size: 15px !important;
        font-weight: 500 !important;
    }
    
    /* Sidebar */
    [data-testid="stSidebar"] {
        font-size: 14px !important;
    }
    
    [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2 {
        font-size: 1.3rem !important;
    }
    
    /* Expander */
    .streamlit-expanderHeader {
        font-size: 15px !important;
        font-weight: 500 !important;
    }
    
    /* Código */
    code {
        font-family: 'Consolas', 'Monaco', monospace !important;
        font-size: 13px !important;
        background: #F5F5F5 !important;
        padding: 2px 6px !important;
        border-radius: 4px !important;
    }
    
    /* Footer */
    .footer {
        text-align: center;
        padding: 2rem;
        color: #666;
        font-size: 14px;
        border-top: 1px solid #E0E0E0;
        margin-top: 2rem;
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# PARSER ANOMALY DETECTOR
# ============================================================================

def parse_anomaly_detector_json(raw_data: dict) -> dict:
    """Convierte JSON de Anomaly Detector al formato CosmicForge"""
    try:
        plugin_results = raw_data.get('plugin_results', {})
        fractal_base = plugin_results.get('fractal_base', {})
        kolmogorov = plugin_results.get('kolmogorov_1941', {})
        lyapunov = plugin_results.get('lyapunov_stability', {})
        renormalization = plugin_results.get('renormalization_group', {})
        anisotropy_data = plugin_results.get('anisotropy', {})
        entropy_data = plugin_results.get('entropy', {})
        
        return {
            'object_name': raw_data.get('filename', 'Unknown').replace('.jpg', '').replace('.png', ''),
            'mode': raw_data.get('mode', 'balanced'),
            'global_score': raw_data.get('global_score', 0.5),
            'fractal_dimension': fractal_base.get('d0', 1.5),
            'criticality_score': renormalization.get('criticality_score', 0.5),
            'entropy': entropy_data.get('normalized_entropy', 0.01),
            'anisotropy': anisotropy_data.get('anisotropy_index', 0.3),
            'turbulence_beta': kolmogorov.get('beta', 2.0),
            'lyapunov_max': lyapunov.get('max_lyapunov', -0.1),
            'source': 'Anomaly Detector'
        }
    except Exception as e:
        st.error(f"Error: {e}")
        return None

# ============================================================================
# PDF GENERATOR - VERSIÓN MEJORADA
# ============================================================================

def create_pdf_report(astro_data, physical_props, recipe, production_guide):
    """Crea PDF legible y bien formateado"""
    
    try:
        from fpdf import FPDF
        
        class PDF(FPDF):
            def header(self):
                self.set_font('Arial', 'B', 12)
                self.cell(0, 10, 'CosmicForge Lab - Informe Tecnico', 0, 1, 'C')
                self.ln(5)
            
            def footer(self):
                self.set_y(-15)
                self.set_font('Arial', 'I', 8)
                self.cell(0, 10, f'Pagina {self.page_no()}', 0, 0, 'C')
        
        pdf = PDF()
        pdf.add_page()
        pdf.set_auto_page_break(auto=True, margin=15)
        
        # Título
        pdf.set_font('Arial', 'B', 20)
        pdf.cell(0, 15, 'COSMICFORGE LAB', 0, 1, 'C')
        pdf.set_font('Arial', '', 12)
        pdf.cell(0, 8, 'Informe de Diseno de Material', 0, 1, 'C')
        pdf.ln(10)
        
        # Información general
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'INFORMACION GENERAL', 0, 1)
        pdf.set_font('Arial', '', 11)
        pdf.cell(0, 8, f"Objeto Astrofisico: {astro_data.get('object_name', 'Unknown')}", 0, 1)
        pdf.cell(0, 8, f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M')}", 0, 1)
        pdf.cell(0, 8, f"Material: {recipe.get('metal', 'Ti')}-O ({recipe.get('material_type', 'Oxido')})", 0, 1)
        pdf.ln(5)
        
        # Propiedades
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'PROPIEDADES DEL MATERIAL', 0, 1)
        pdf.set_font('Arial', '', 11)
        
        props = [
            f"Porosidad: {physical_props.get('porosity', 0):.2%}",
            f"Densidad: {physical_props.get('density', 0):.3f} g/cm3",
            f"Conductividad Termica: {physical_props.get('thermal_conductivity', 0):.2f} W/mK",
            f"Modulo Elastico: {physical_props.get('elastic_modulus', 0):.2f} GPa",
            f"Area Superficial: {physical_props.get('surface_area', 0):.1f} m2/g",
            f"Band Gap: {physical_props.get('band_gap', 0):.3f} eV",
            f"Calidad: {physical_props.get('quality_score', 0):.1%}"
        ]
        
        for prop in props:
            pdf.cell(0, 7, f"  - {prop}", 0, 1)
        pdf.ln(5)
        
        # Proceso de producción
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'PROCESO DE PRODUCCION', 0, 1)
        pdf.set_font('Arial', '', 10)
        
        if production_guide:
            for step in production_guide.get('synthesis', []):
                if step.strip():
                    # Limpiar caracteres especiales
                    step_clean = step.encode('latin-1', 'replace').decode('latin-1')
                    pdf.multi_cell(0, 6, f"  {step_clean}")
        pdf.ln(5)
        
        # Equipamiento
        pdf.set_font('Arial', 'B', 12)
        pdf.cell(0, 8, 'Equipamiento Necesario:', 0, 1)
        pdf.set_font('Arial', '', 10)
        for eq in production_guide.get('equipment', []):
            pdf.cell(0, 6, f"  - {eq}", 0, 1)
        pdf.ln(3)
        
        # Precursores
        pdf.set_font('Arial', 'B', 12)
        pdf.cell(0, 8, 'Precursores:', 0, 1)
        pdf.set_font('Arial', '', 10)
        for prec in production_guide.get('precursors', []):
            pdf.cell(0, 6, f"  - {prec}", 0, 1)
        pdf.ln(3)
        
        # Seguridad
        pdf.set_font('Arial', 'B', 12)
        pdf.cell(0, 8, 'Seguridad:', 0, 1)
        pdf.set_font('Arial', '', 10)
        for safe in production_guide.get('safety', []):
            pdf.cell(0, 6, f"  - {safe}", 0, 1)
        
        # Nueva página para conclusión
        pdf.add_page()
        
        # Conclusión filosófica
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'REFLEXION CIENTIFICA', 0, 1)
        
        quote = MANIN_QUOTES[hash(astro_data.get('object_name', '')) % len(MANIN_QUOTES)]
        pdf.set_font('Arial', 'I', 11)
        pdf.multi_cell(0, 7, f'"{quote}"')
        pdf.set_font('Arial', '', 10)
        pdf.cell(0, 8, '- Yuri I. Manin, "Lo demostrable e indemostrable"', 0, 1, 'R')
        pdf.ln(10)
        
        # Análisis
        pdf.set_font('Arial', 'B', 12)
        pdf.cell(0, 10, 'ANALISIS EPISTEMOLOGICO', 0, 1)
        pdf.set_font('Arial', '', 10)
        
        analysis = f"""El material predicho a partir de {astro_data.get('object_name', 'este objeto astrofisico')} representa una hipotesis cientifica derivada de analogias matematicas entre estructuras cosmicas y propiedades materiales.

Nivel de certeza matematica: {'Alto' if physical_props.get('quality_score', 0) > 0.6 else 'Moderado'} - Las ecuaciones estan formalmente derivadas.

Nivel de certeza fisica: Requiere validacion - La traduccion astrofisica a material es una hipotesis por verificar experimentalmente.

Proximos pasos recomendados:
1. Simulacion DFT con VASP o Quantum ESPRESSO
2. Validacion experimental de propiedades
3. Comparacion con bases de datos (Materials Project, AFLOW)"""
        
        pdf.multi_cell(0, 6, analysis)
        
        # Footer
        pdf.ln(20)
        pdf.set_font('Arial', 'I', 9)
        pdf.cell(0, 8, 'Generado por CosmicForge Lab', 0, 1, 'C')
        
        return pdf.output(dest='S').encode('latin-1')
        
    except Exception as e:
        # Fallback: crear texto simple
        report_text = f"""
COSMICFORGE LAB - INFORME TECNICO
================================

Objeto Astrofisico: {astro_data.get('object_name', 'Unknown')}
Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M')}
Material: {recipe.get('metal', 'Ti')}-O

PROPIEDADES:
- Porosidad: {physical_props.get('porosity', 0):.2%}
- Densidad: {physical_props.get('density', 0):.3f} g/cm3
- Conductividad: {physical_props.get('thermal_conductivity', 0):.2f} W/mK
- Modulo Elastico: {physical_props.get('elastic_modulus', 0):.2f} GPa

PROCESO DE PRODUCCION:
{chr(10).join(production_guide.get('synthesis', []))}

EQUIPAMIENTO:
{chr(10).join(['- ' + e for e in production_guide.get('equipment', [])])}

SEGURIDAD:
{chr(10).join(['- ' + s for s in production_guide.get('safety', [])])}

---
Generado por CosmicForge Lab
"""
        return report_text.encode('utf-8')

# ============================================================================
# VISUALIZACIÓN 3D
# ============================================================================

def create_3d_structure(porosity, density, metal="Ti"):
    if not PLOTLY_AVAILABLE:
        return None
    
    np.random.seed(42)
    n_atoms = int(200 * (1 - porosity))
    n_metal = int(n_atoms * 0.4)
    n_oxygen = n_atoms - n_metal
    
    fig = go.Figure()
    
    # Metal atoms
    fig.add_trace(go.Scatter3d(
        x=np.random.randn(n_metal) * 5,
        y=np.random.randn(n_metal) * 5,
        z=np.random.randn(n_metal) * 5,
        mode='markers',
        marker=dict(size=10, color='#E53935', opacity=0.8),
        name=f'{metal} (Metal)'
    ))
    
    # Oxygen atoms
    fig.add_trace(go.Scatter3d(
        x=np.random.randn(n_oxygen) * 5,
        y=np.random.randn(n_oxygen) * 5,
        z=np.random.randn(n_oxygen) * 5,
        mode='markers',
        marker=dict(size=6, color='#43A047', opacity=0.7),
        name='O (Oxigeno)'
    ))
    
    fig.update_layout(
        title=f'Estructura Cristalina 3D - Porosidad: {porosity:.1%}',
        scene=dict(xaxis_title='X (A)', yaxis_title='Y (A)', zaxis_title='Z (A)'),
        height=450
    )
    
    return fig

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def validate_astrophysical_data(data):
    errors = []
    required = ['fractal_dimension', 'criticality_score', 'entropy', 'anisotropy', 'turbulence_beta']
    for field in required:
        if field not in data:
            errors.append(f"Campo faltante: {field}")
    return len(errors) == 0, errors

def display_metrics(data):
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Dimension Fractal", f"{data.get('fractal_dimension', 0):.4f}")
        st.metric("Criticalidad", f"{data.get('criticality_score', 0):.4f}")
    with col2:
        st.metric("Entropia", f"{data.get('entropy', 0):.6f}")
        st.metric("Anisotropia", f"{data.get('anisotropy', 0):.4f}")
    with col3:
        st.metric("Turbulencia B", f"{data.get('turbulence_beta', 0):.4f}")
        st.metric("Lyapunov max", f"{data.get('lyapunov_max', 0):.6f}")
    with col4:
        st.metric("Objeto", str(data.get('object_name', 'N/A'))[:15])
        st.metric("Fuente", str(data.get('source', 'Manual'))[:10])

# ============================================================================
# MAIN APPLICATION
# ============================================================================

# Header
st.title("CosmicForge Lab")
st.markdown("<p style='text-align:center; color:#666; font-size:18px;'>Diseno de Materiales Inspirado en Firmas Astrofisicas</p>", unsafe_allow_html=True)

# Sidebar
st.sidebar.header("Importar Datos")

input_method = st.sidebar.radio("Metodo:", ["Ejemplos Astrofisicos", "JSON Anomaly Detector", "Editor Manual"])

if input_method == "Ejemplos Astrofisicos":
    example = st.sidebar.selectbox("Selecciona:", list(ASTROPHYSICAL_EXAMPLES.keys()))
    if st.sidebar.button("Cargar Ejemplo", type="primary"):
        data = ASTROPHYSICAL_EXAMPLES[example].copy()
        data['object_name'] = example
        data['source'] = 'Base de Datos'
        st.session_state.astro_data = data
        st.sidebar.success(f"Cargado: {example}")
    
    if example in ASTROPHYSICAL_EXAMPLES:
        ex = ASTROPHYSICAL_EXAMPLES[example]
        st.sidebar.info(f"Tipo: {ex.get('type', 'N/A')}\nDistancia: {ex.get('distance', 'N/A')}")

elif input_method == "JSON Anomaly Detector":
    uploaded = st.sidebar.file_uploader("Archivo JSON", type=['json'])
    if uploaded:
        try:
            raw = json.loads(uploaded.read().decode('utf-8'))
            if 'plugin_results' in raw:
                st.session_state.astro_data = parse_anomaly_detector_json(raw)
                st.sidebar.success("JSON Anomaly Detector cargado")
            else:
                is_valid, _ = validate_astrophysical_data(raw)
                if is_valid:
                    st.session_state.astro_data = raw
                    st.sidebar.success("JSON cargado")
        except Exception as e:
            st.sidebar.error(f"Error: {e}")

elif input_method == "Editor Manual":
    name = st.sidebar.text_input("Nombre", "Objeto Personalizado")
    fd = st.sidebar.slider("Dimension Fractal", 0.0, 3.0, 1.5)
    cs = st.sidebar.slider("Criticalidad", 0.0, 1.0, 0.5)
    en = st.sidebar.number_input("Entropia", 0.0, 1.0, 0.01, format="%.6f")
    an = st.sidebar.slider("Anisotropia", 0.0, 1.0, 0.3)
    tb = st.sidebar.slider("Turbulencia", 0.0, 5.0, 2.0)
    lm = st.sidebar.number_input("Lyapunov", -1.0, 1.0, -0.1)
    
    if st.sidebar.button("Aplicar", type="primary"):
        st.session_state.astro_data = {
            'object_name': name, 'fractal_dimension': fd, 'criticality_score': cs,
            'entropy': en, 'anisotropy': an, 'turbulence_beta': tb, 'lyapunov_max': lm,
            'mode': 'custom', 'source': 'Editor'
        }
        st.sidebar.success("Parametros aplicados")

# Material selection
st.sidebar.header("Configuracion del Material")

material_type = st.sidebar.selectbox("Tipo:", ["Oxido", "Metal", "Ceramica", "Nanoparticula"])
metal = st.sidebar.selectbox("Metal:", list(ELEMENTS_DATABASE.keys()))

# Aleación
st.sidebar.subheader("Sistema de Aleacion")
use_alloy = st.sidebar.checkbox("Usar aleacion predefinida")

selected_alloy_key = None
if use_alloy:
    alloy_options = [k for k in ALLOY_SYSTEMS.keys() if metal in k]
    if alloy_options:
        selected_alloy_key = st.sidebar.selectbox("Aleacion:", alloy_options)
        if selected_alloy_key:
            alloy_info = ALLOY_SYSTEMS[selected_alloy_key]
            st.sidebar.info(f"{alloy_info['name']}\nRatio: {alloy_info['ratio']}\nMP: {alloy_info['melting_point']}")

# Limpiar
if st.session_state.astro_data and st.sidebar.button("Limpiar Datos"):
    st.session_state.astro_data = None
    st.session_state.physical_props = None
    st.session_state.recipe = None
    st.rerun()

# Main content
if st.session_state.astro_data:
    astro_data = st.session_state.astro_data
    
    st.header("Firma Astrofisica Detectada")
    display_metrics(astro_data)
    
    # Elementos compatibles
    st.markdown("---")
    st.subheader("Elementos Compatibles")
    if metal in ELEMENTS_DATABASE:
        compat = ELEMENTS_DATABASE[metal]['compatible']
        st.write(f"**{ELEMENTS_DATABASE[metal]['name']}** es compatible con: **{', '.join(compat)}**")
    
    # Generar
    if st.button("Generar Receta de Material", type="primary"):
        with st.spinner("Calculando..."):
            try:
                physics = PhysicsCalculator()
                chem = ChemistryEngine()
                
                physical_props = physics.calculate_all_properties(astro_data)
                recipe = chem.generate_complete_recipe(metal=metal, material_type=material_type.lower(), astrophysical_data=astro_data)
                
                st.session_state.physical_props = physical_props
                st.session_state.recipe = recipe
                st.success("Calculos completados!")
                st.rerun()
            except Exception as e:
                st.error(f"Error: {e}")

# Results
if st.session_state.physical_props and st.session_state.recipe:
    pp = st.session_state.physical_props
    rc = st.session_state.recipe
    ad = st.session_state.astro_data
    
    # Obtener guía de producción
    production_guide = ALLOY_SYSTEMS.get(selected_alloy_key) if selected_alloy_key else None
    
    # Si no hay guía de aleación, crear guía genérica
    if not production_guide:
        production_guide = {
            'name': f"{metal}-Oxido",
            'synthesis': [
                f"1. Preparar solucion de precursor de {metal} (0.5M) en agua desionizada",
                "2. Ajustar pH con NaOH o NH4OH segun requerimientos",
                f"3. Calentar a {int(rc.get('reaction_conditions', {}).get('temperature_C', 400))}C",
                f"4. Mantener por {rc.get('reaction_conditions', {}).get('reaction_time_hours', 2):.1f} horas",
                "5. Filtrar y lavar precipitado con agua/etanol",
                "6. Secar a 100C por 12 horas",
                "7. Calcinacion a 400-600C por 4 horas",
                "8. Caracterizar por XRD, SEM, BET"
            ],
            'equipment': ["Matraz de reaccion", "Agitador magnetico", "Horno mufla", "Filtro", "Estufa de secado"],
            'precursors': [f"Nitrato de {metal}", "NaOH", "Agua desionizada", "Etanol"],
            'safety': ["Guantes y gafas", "Campana extractora", "Evitar inhalacion de polvos"]
        }
    
    tabs = st.tabs(["Propiedades", "Sintesis", "Produccion Detallada", "Estructura 3D", "Archivos", "PDF", "Conclusion"])
    
    with tabs[0]:
        st.subheader("Propiedades del Material")
        c1, c2, c3, c4 = st.columns(4)
        with c1:
            st.metric("Porosidad", f"{pp.get('porosity', 0):.1%}")
            st.metric("Densidad", f"{pp.get('density', 0):.3f} g/cm3")
        with c2:
            st.metric("Conductividad Termica", f"{pp.get('thermal_conductivity', 0):.2f} W/mK")
            st.metric("Modulo Elastico", f"{pp.get('elastic_modulus', 0):.2f} GPa")
        with c3:
            st.metric("Area Superficial", f"{pp.get('surface_area', 0):.1f} m2/g")
            st.metric("Band Gap", f"{pp.get('band_gap', 0):.3f} eV")
        with c4:
            st.metric("Calidad", f"{pp.get('quality_score', 0):.1%}")
            st.metric("Energia Activacion", f"{pp.get('activation_energy', 0):.4f} eV")
        
        # Materials Project links
        st.markdown("### Validacion con Bases de Datos")
        formula = f"{metal}O2"
        st.markdown(f"""
        <div class="info-card">
        <strong>Consulta tu material:</strong><br>
        - <a href="https://materialsproject.org/materials?search={formula}" target="_blank">Materials Project - {formula}</a><br>
        - <a href="https://www.aflowlib.org/?search={formula}" target="_blank">AFLOW Library</a><br>
        - <a href="http://oqmd.org/materials?search={formula}" target="_blank">OQMD Database</a>
        </div>
        """, unsafe_allow_html=True)
    
    with tabs[1]:
        st.subheader("Receta de Sintesis")
        cond = rc.get('reaction_conditions', {})
        
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("#### Condiciones de Reaccion")
            st.write(f"**Temperatura:** {cond.get('temperature_C', 0):.0f} C")
            st.write(f"**Tiempo:** {cond.get('reaction_time_hours', 0):.2f} horas")
            st.write(f"**pH:** {cond.get('pH', 7):.1f}")
        
        with c2:
            st.markdown("#### Precursores")
            for p in rc.get('precursors', []):
                st.write(f"- {p}")
        
        st.markdown("#### Procedimiento General")
        for i, step in enumerate(rc.get('step_by_step', []), 1):
            st.write(f"**{i}.** {step}")
    
    with tabs[2]:
        st.subheader("Guia de Produccion Detallada")
        
        if selected_alloy_key:
            st.markdown(f"<div class='info-card'><h3>{production_guide['name']}</h3><p><strong>Ratio:</strong> {production_guide['ratio']} | <strong>Punto de Fusion:</strong> {production_guide['melting_point']} | <strong>Densidad:</strong> {production_guide['density']}</p></div>", unsafe_allow_html=True)
        
        st.markdown("#### Aplicaciones")
        apps = production_guide.get('applications', ['Ciencia de materiales', 'Investigacion'])
        for app in apps:
            st.write(f"- {app}")
        
        st.markdown("#### Proceso de Sintesis Paso a Paso")
        for step in production_guide.get('synthesis', []):
            if step.strip():
                st.markdown(f"<div class='step-box'>{step}</div>", unsafe_allow_html=True)
        
        c1, c2, c3 = st.columns(3)
        with c1:
            st.markdown("#### Equipamiento")
            for eq in production_guide.get('equipment', []):
                st.write(f"- {eq}")
        
        with c2:
            st.markdown("#### Precursores")
            for prec in production_guide.get('precursors', []):
                st.write(f"- {prec}")
        
        with c3:
            st.markdown("#### Seguridad")
            st.markdown("<div class='warning-card'>", unsafe_allow_html=True)
            for safe in production_guide.get('safety', []):
                st.write(f"- {safe}")
            st.markdown("</div>", unsafe_allow_html=True)
    
    with tabs[3]:
        st.subheader("Visualizacion 3D")
        if PLOTLY_AVAILABLE:
            fig = create_3d_structure(pp.get('porosity', 0.3), pp.get('density', 2.0), metal)
            if fig:
                st.plotly_chart(fig, use_container_width=True)
                st.info("Arrastra para rotar | Scroll para zoom")
        else:
            st.warning("Instala plotly: pip install plotly")
    
    with tabs[4]:
        st.subheader("Archivos para Simuladores")
        fg = FileGenerator()
        
        c1, c2, c3 = st.columns(3)
        with c1:
            poscar = fg.generate_poscar(pp, rc, 'oxide')
            st.download_button("POSCAR (VASP/QE)", poscar, "POSCAR.vasp", "text/plain")
        with c2:
            lammps = fg.generate_lammps_input(pp, rc, 'oxide')
            st.download_button("LAMMPS", lammps, "input.lmp", "text/plain")
        with c3:
            csv = fg.generate_properties_csv(pp, rc)
            st.download_button("CSV", csv, "properties.csv", "text/csv")
    
    with tabs[5]:
        st.subheader("Generar Reporte PDF")
        
        if st.button("Crear PDF", type="primary"):
            pdf_bytes = create_pdf_report(ad, pp, rc, production_guide)
            
            if isinstance(pdf_bytes, bytes):
                st.download_button(
                    "Descargar PDF",
                    pdf_bytes,
                    f"CosmicForge_{ad.get('object_name', 'material')}.pdf",
                    "application/pdf"
                )
            else:
                st.download_button(
                    "Descargar TXT (fallback)",
                    pdf_bytes,
                    f"CosmicForge_{ad.get('object_name', 'material')}.txt",
                    "text/plain"
                )
    
    with tabs[6]:
        st.subheader("Conclusion Cientifica")
        
        quote = MANIN_QUOTES[hash(ad.get('object_name', '')) % len(MANIN_QUOTES)]
        st.markdown(f"""
        <div class="quote-card">
        <em>"{quote}"</em><br><br>
        <strong>- Yuri I. Manin, "Lo demostrable e indemostrable"</strong>
        </div>
        """, unsafe_allow_html=True)
        
        quality = pp.get('quality_score', 0.5)
        st.markdown(f"""
        <div class="production-card">
        <h4>Analisis Epistemologico</h4>
        <p>La prediccion del material <strong>{metal}-O</strong> a partir de <strong>{ad.get('object_name', 'este objeto')}</strong> 
        ejemplifica la tension entre lo demostrable (las ecuaciones matematicas) y lo indemostrable (la validez fisica ultima).</p>
        
        <p><strong>Nivel de certeza matematica:</strong> {'Alto' if quality > 0.6 else 'Moderado'} - Ecuaciones formalmente derivadas.</p>
        <p><strong>Nivel de certeza fisica:</strong> Requiere validacion experimental.</p>
        
        <h4>Proximos Pasos:</h4>
        <ol>
        <li>Simulacion DFT con VASP o Quantum ESPRESSO</li>
        <li>Validacion experimental en laboratorio</li>
        <li>Comparacion con Materials Project</li>
        </ol>
        </div>
        """, unsafe_allow_html=True)

else:
    st.info("Carga un ejemplo o archivo JSON para comenzar")

# Footer
st.markdown("---")
st.markdown("<div class='footer'>CosmicForge Lab v3.1 | Edicion Personal</div>", unsafe_allow_html=True)
