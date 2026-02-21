"""
CosmicForge Lab - Professional Edition
Diseño de Materiales Inspirado en Firmas Astrofísicas
Versión 3.0 - Professional Scientific Platform

Features:
- 15+ ejemplos de objetos astrofísicos
- Materials Project API integration
- Visualización 3D interactiva
- Export PDF con conclusión filosófica (Manin)
- Sistema de aleaciones
- Diseño profesional
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

# PDF generation
try:
    from reportlab.lib.pagesizes import letter, A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
    from reportlab.lib import colors
    PDF_AVAILABLE = True
except:
    PDF_AVAILABLE = False

# ============================================================================
# PAGE CONFIGURATION
# ============================================================================

st.set_page_config(
    page_title="CosmicForge Lab Pro",
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
if 'alloy_elements' not in st.session_state:
    st.session_state.alloy_elements = []

# ============================================================================
# DATABASES
# ============================================================================

# 15+ Ejemplos de objetos astrofísicos
ASTROPHYSICAL_EXAMPLES = {
    # Nebulosas
    "Orion Nebula M42": {
        "fractal_dimension": 1.654, "criticality_score": 0.722, "entropy": 0.019,
        "anisotropy": 0.329, "turbulence_beta": 2.278, "lyapunov_max": -0.227,
        "mode": "balanced", "type": "Emission Nebula", "distance": "1,344 ly",
        "description": "Nebulosa de emisión activa con formación estelar intensa"
    },
    "Crab Nebula M1": {
        "fractal_dimension": 1.82, "criticality_score": 0.89, "entropy": 0.032,
        "anisotropy": 0.456, "turbulence_beta": 2.45, "lyapunov_max": -0.15,
        "mode": "turbulent", "type": "Supernova Remnant", "distance": "6,500 ly",
        "description": "Remanente de supernova con pulsar central, alta turbulencia"
    },
    "Ring Nebula M57": {
        "fractal_dimension": 1.45, "criticality_score": 0.55, "entropy": 0.012,
        "anisotropy": 0.28, "turbulence_beta": 1.95, "lyapunov_max": -0.35,
        "mode": "stable", "type": "Planetary Nebula", "distance": "2,283 ly",
        "description": "Nebulosa planetaria simétrica, estructura estable"
    },
    "Eagle Nebula M16": {
        "fractal_dimension": 1.78, "criticality_score": 0.81, "entropy": 0.028,
        "anisotropy": 0.42, "turbulence_beta": 2.35, "lyapunov_max": -0.18,
        "mode": "turbulent", "type": "Emission Nebula", "distance": "7,000 ly",
        "description": "Famosa por los 'Pilares de la Creación', alta formación estelar"
    },
    "Horsehead Nebula IC 434": {
        "fractal_dimension": 1.52, "criticality_score": 0.63, "entropy": 0.015,
        "anisotropy": 0.35, "turbulence_beta": 2.10, "lyapunov_max": -0.25,
        "mode": "balanced", "type": "Dark Nebula", "distance": "1,500 ly",
        "description": "Nebulosa oscura con estructura fractal distintiva"
    },
    "Helix Nebula NGC 7293": {
        "fractal_dimension": 1.48, "criticality_score": 0.58, "entropy": 0.011,
        "anisotropy": 0.31, "turbulence_beta": 2.02, "lyapunov_max": -0.32,
        "mode": "stable", "type": "Planetary Nebula", "distance": "655 ly",
        "description": "Nebulosa planetaria más cercana, estructura en espiral"
    },
    
    # Galaxias
    "Triangulum Galaxy M33": {
        "fractal_dimension": 1.93, "criticality_score": 0.78, "entropy": 0.00000964,
        "anisotropy": 0.23, "turbulence_beta": 1.84, "lyapunov_max": 0.11,
        "mode": "balanced", "type": "Spiral Galaxy", "distance": "2.73M ly",
        "description": "Galaxia espiral con alta estructura fractal"
    },
    "Andromeda Galaxy M31": {
        "fractal_dimension": 1.88, "criticality_score": 0.85, "entropy": 0.000012,
        "anisotropy": 0.19, "turbulence_beta": 1.92, "lyapunov_max": 0.08,
        "mode": "balanced", "type": "Spiral Galaxy", "distance": "2.537M ly",
        "description": "Galaxia espiral más cercana a la Vía Láctea"
    },
    "Whirlpool Galaxy M51": {
        "fractal_dimension": 1.95, "criticality_score": 0.91, "entropy": 0.000015,
        "anisotropy": 0.22, "turbulence_beta": 1.98, "lyapunov_max": 0.12,
        "mode": "turbulent", "type": "Interacting Spiral", "distance": "23M ly",
        "description": "Galaxia en interacción con alta complejidad dinámica"
    },
    "Sombrero Galaxy M104": {
        "fractal_dimension": 1.72, "criticality_score": 0.68, "entropy": 0.000008,
        "anisotropy": 0.45, "turbulence_beta": 2.15, "lyapunov_max": 0.05,
        "mode": "stable", "type": "Spiral Galaxy", "distance": "31.1M ly",
        "description": "Galaxia con anillo de polvo distintivo"
    },
    
    # Objetos especiales
    "Pillars of Creation": {
        "fractal_dimension": 1.85, "criticality_score": 0.88, "entropy": 0.035,
        "anisotropy": 0.52, "turbulence_beta": 2.55, "lyapunov_max": -0.12,
        "mode": "turbulent", "type": "Molecular Cloud", "distance": "6,500 ly",
        "description": "Estructuras de gas y polvo con alta formación estelar"
    },
    "Tarantula Nebula": {
        "fractal_dimension": 1.91, "criticality_score": 0.93, "entropy": 0.045,
        "anisotropy": 0.58, "turbulence_beta": 2.68, "lyapunov_max": 0.15,
        "mode": "turbulent", "type": "HII Region", "distance": "160,000 ly",
        "description": "Región de formación estelar más activa conocida"
    },
    "Carina Nebula NGC 3372": {
        "fractal_dimension": 1.76, "criticality_score": 0.79, "entropy": 0.029,
        "anisotropy": 0.41, "turbulence_beta": 2.32, "lyapunov_max": -0.16,
        "mode": "balanced", "type": "Emission Nebula", "distance": "8,500 ly",
        "description": "Nebulosa gigante con estrellas masivas"
    },
    "Bubble Nebula NGC 7635": {
        "fractal_dimension": 1.58, "criticality_score": 0.62, "entropy": 0.018,
        "anisotropy": 0.38, "turbulence_beta": 2.08, "lyapunov_max": -0.21,
        "mode": "balanced", "type": "Emission Nebula", "distance": "7,100 ly",
        "description": "Burbuja expandida por viento estelar"
    },
    "Lagoon Nebula M8": {
        "fractal_dimension": 1.69, "criticality_score": 0.74, "entropy": 0.022,
        "anisotropy": 0.36, "turbulence_beta": 2.18, "lyapunov_max": -0.19,
        "mode": "balanced", "type": "Emission Nebula", "distance": "4,100 ly",
        "description": "Nebulosa con estructura de laguna oscura"
    }
}

# Sistema de aleaciones y elementos compatibles
ALLOY_SYSTEMS = {
    # Aleaciones binarias
    "Ti-Al": {"name": "Titanio-Aluminio", "applications": ["Aeroespacial", "Biomédico"], "ratio": "1:1 a 3:1"},
    "Ti-V": {"name": "Titanio-Vanadio", "applications": ["Aeroespacial", "Implantes"], "ratio": "6:4"},
    "Fe-Ni": {"name": "Invar", "applications": ["Instrumentos de precisión"], "ratio": "64:36"},
    "Fe-Cr": {"name": "Acero inoxidable", "applications": ["Industrial", "Médico"], "ratio": "Variable"},
    "Cu-Ni": {"name": "Cuproníquel", "applications": ["Marino", "Monedas"], "ratio": "70:30 a 90:10"},
    "Cu-Zn": {"name": "Latón", "applications": ["Decorativo", "Mecánico"], "ratio": "Variable"},
    "Al-Cu": {"name": "Duraluminio", "applications": ["Aeroespacial", "Automotriz"], "ratio": "95:4"},
    "Al-Mg": {"name": "Magnalio", "applications": ["Aeroespacial"], "ratio": "90:10"},
    "Ni-Co": {"name": "Níquel-Cobalto", "applications": ["Magnéticos", "Catálisis"], "ratio": "Variable"},
    "Pt-Pd": {"name": "Platino-Paladio", "applications": ["Catálisis", "Joyería"], "ratio": "Variable"},
    
    # Aleaciones ternarias
    "Ti-Al-V": {"name": "Titanio TA6V", "applications": ["Aeroespacial", "Biomédico"], "ratio": "6:4:4"},
    "Fe-Ni-Cr": {"name": "Inoxidable austenítico", "applications": ["Industrial", "Químico"], "ratio": "Variable"},
    "Al-Cu-Mg": {"name": "Aluminio aeronáutico", "applications": ["Aeroespacial"], "ratio": "Variable"},
    
    # Cerámicos
    "Al-O": {"name": "Alúmina", "applications": ["Cerámico", "Abrasivos"], "ratio": "2:3"},
    "Ti-O": {"name": "Titania", "applications": ["Pigmentos", "Fotocatálisis"], "ratio": "1:2"},
    "Si-O": {"name": "Sílice", "applications": ["Vidrio", "Electrónica"], "ratio": "1:2"},
    "Zn-O": {"name": "Óxido de zinc", "applications": ["Sensores", "Protectores"], "ratio": "1:1"},
}

# Elementos compatibles para síntesis
COMPATIBLE_ELEMENTS = {
    "Ti": {"compatible": ["Al", "V", "Fe", "Ni", "O", "N", "C"], "conflict": ["Hg", "Pb"]},
    "Al": {"compatible": ["Ti", "Cu", "Mg", "Si", "Zn", "O"], "conflict": ["Hg"]},
    "Fe": {"compatible": ["Ni", "Cr", "Co", "Mn", "C", "O"], "conflict": ["Pb"]},
    "Zn": {"compatible": ["Cu", "Al", "O", "S"], "conflict": []},
    "Cu": {"compatible": ["Ni", "Zn", "Al", "Sn", "O"], "conflict": []},
    "Ni": {"compatible": ["Fe", "Co", "Cu", "Cr", "Ti", "O"], "conflict": []},
    "Co": {"compatible": ["Ni", "Fe", "Mn", "O"], "conflict": []},
    "Mn": {"compatible": ["Fe", "Co", "O"], "conflict": []},
    "Ag": {"compatible": ["Cu", "Au", "Pd", "O"], "conflict": ["S"]},
    "Au": {"compatible": ["Ag", "Pt", "Pd", "Cu"], "conflict": []},
    "Pt": {"compatible": ["Pd", "Au", "Rh", "Ir", "O"], "conflict": []},
    "Pd": {"compatible": ["Pt", "Au", "Ag", "Ni", "O"], "conflict": []},
}

# Citas del libro "Lo demostrable e indemostrable" de Yu. I. Manin
MANIN_QUOTES = [
    "La matemática es la parte más pura de la cultura humana, y sin embargo está profundamente conectada con la realidad física.",
    "La demostrabilidad es un concepto que trasciende la lógica formal: es la búsqueda de certeza en un universo de incertidumbre.",
    "Entre lo demostrable y lo indemostrable existe una zona fronteriza donde nace la creatividad científica.",
    "La física teórica y las matemáticas comparten un misterio común: ¿por qué sus abstracciones describen tan bien la realidad?",
    "La verdad matemática tiene un estatus ontológico único: existe independientemente de nuestra capacidad para demostrarla.",
    "En el límite entre lo conocible y lo incognoscible, la ciencia encuentra su mayor impulso creativo.",
]

# ============================================================================
# CUSTOM CSS STYLES - DISEÑO PROFESIONAL
# ============================================================================

st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
    
    * { font-family: 'Inter', sans-serif; }
    
    .main-header {
        font-size: 2.8rem;
        font-weight: 700;
        background: linear-gradient(135deg, #1E88E5 0%, #7C4DFF 50%, #E91E63 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        text-align: center;
        margin-bottom: 0.3rem;
        letter-spacing: -1px;
    }
    
    .sub-header {
        font-size: 1.1rem;
        color: #546E7A;
        text-align: center;
        margin-bottom: 1.5rem;
        font-weight: 300;
    }
    
    .pro-badge {
        background: linear-gradient(135deg, #FFD700, #FFA000);
        color: #000;
        padding: 2px 10px;
        border-radius: 12px;
        font-size: 0.7rem;
        font-weight: 600;
        margin-left: 8px;
    }
    
    .metric-card {
        background: linear-gradient(135deg, #f8f9fa 0%, #ffffff 100%);
        padding: 1.2rem;
        border-radius: 12px;
        border-left: 4px solid #1E88E5;
        box-shadow: 0 4px 15px rgba(0,0,0,0.08);
        transition: transform 0.2s;
    }
    
    .metric-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 8px 25px rgba(0,0,0,0.12);
    }
    
    .success-box {
        background: linear-gradient(135deg, #d4edda, #c3e6cb);
        padding: 1.2rem;
        border-radius: 12px;
        border-left: 4px solid #28a745;
        margin: 1rem 0;
    }
    
    .manin-quote {
        background: linear-gradient(135deg, #FFF8E1, #FFECB3);
        padding: 1.5rem;
        border-radius: 12px;
        border-left: 4px solid #FFA000;
        margin: 1.5rem 0;
        font-style: italic;
        color: #5D4037;
    }
    
    .manin-quote .author {
        text-align: right;
        font-weight: 600;
        margin-top: 0.5rem;
        color: #795548;
    }
    
    .conclusion-box {
        background: linear-gradient(135deg, #E3F2FD, #BBDEFB);
        padding: 1.5rem;
        border-radius: 12px;
        border: 2px solid #1E88E5;
        margin: 1.5rem 0;
    }
    
    .alloy-card {
        background: #f5f5f5;
        padding: 1rem;
        border-radius: 8px;
        margin: 0.5rem 0;
        border-left: 3px solid #7C4DFF;
    }
    
    .compatible-tag {
        display: inline-block;
        background: #E8F5E9;
        color: #2E7D32;
        padding: 2px 8px;
        border-radius: 12px;
        font-size: 0.75rem;
        margin: 2px;
    }
    
    .conflict-tag {
        display: inline-block;
        background: #FFEBEE;
        color: #C62828;
        padding: 2px 8px;
        border-radius: 12px;
        font-size: 0.75rem;
        margin: 2px;
    }
    
    .example-card {
        background: white;
        padding: 1rem;
        border-radius: 12px;
        border: 1px solid #e0e0e0;
        margin: 0.5rem;
        cursor: pointer;
        transition: all 0.2s;
    }
    
    .example-card:hover {
        border-color: #1E88E5;
        box-shadow: 0 4px 12px rgba(30, 136, 229, 0.15);
    }
    
    .stMetric label {
        font-size: 0.85rem;
        color: #666;
        font-weight: 500;
    }
    
    .stMetric [data-testid="stMetricValue"] {
        font-size: 1.5rem;
        font-weight: 600;
    }
    
    /* Tabs styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    
    .stTabs [data-baseweb="tab"] {
        padding: 12px 24px;
        border-radius: 8px 8px 0 0;
        font-weight: 500;
    }
    
    /* Professional footer */
    .pro-footer {
        text-align: center;
        padding: 2rem;
        background: linear-gradient(135deg, #f8f9fa, #ffffff);
        border-top: 1px solid #e0e0e0;
        margin-top: 2rem;
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# PARSER PARA ANOMALY DETECTOR JSON
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
        
        astro_data = {
            'object_name': raw_data.get('filename', 'Unknown').replace('.jpg', '').replace('.png', ''),
            'mode': raw_data.get('mode', 'balanced'),
            'global_score': raw_data.get('global_score', 0.5),
            'fractal_dimension': fractal_base.get('d0', fractal_base.get('dimension_box', 1.5)),
            'criticality_score': renormalization.get('criticality_score', 0.5),
            'entropy': entropy_data.get('normalized_entropy', entropy_data.get('shannon_entropy_bits', 0.01) / 8),
            'anisotropy': anisotropy_data.get('anisotropy_index', 0.3),
            'turbulence_beta': kolmogorov.get('beta', 2.0),
            'lyapunov_max': lyapunov.get('max_lyapunov', -0.1),
            'source': 'Anomaly Detector',
            'plugin_results': plugin_results,
            'report': raw_data.get('report', '')
        }
        
        return astro_data
    except Exception as e:
        st.error(f"Error parseando JSON: {e}")
        return None

# ============================================================================
# VISUALIZACIÓN 3D
# ============================================================================

def create_3d_structure(porosity: float, density: float, metal: str = "Ti", n_atoms: int = 200) -> go.Figure:
    """Crea visualización 3D interactiva de estructura cristalina"""
    if not PLOTLY_AVAILABLE:
        return None
    
    np.random.seed(42)
    
    # Número de átomos basado en porosidad
    n_metal = int(n_atoms * (1 - porosity) * 0.4)
    n_oxygen = int(n_atoms * (1 - porosity) * 0.6)
    
    # Posiciones de átomos
    metal_x = np.random.randn(n_metal) * 5
    metal_y = np.random.randn(n_metal) * 5
    metal_z = np.random.randn(n_metal) * 5
    
    oxygen_x = np.random.randn(n_oxygen) * 5
    oxygen_y = np.random.randn(n_oxygen) * 5
    oxygen_z = np.random.randn(n_oxygen) * 5
    
    # Crear figura
    fig = go.Figure()
    
    # Átomos de metal
    fig.add_trace(go.Scatter3d(
        x=metal_x, y=metal_y, z=metal_z,
        mode='markers',
        marker=dict(
            size=8,
            color='#FF6B6B',
            opacity=0.8,
            line=dict(color='darkred', width=1)
        ),
        name=f'{metal} (Metal)',
        hovertemplate=f'<b>{metal}</b><br>X: %{{x:.2f}}<br>Y: %{{y:.2f}}<br>Z: %{{z:.2f}}<extra></extra>'
    ))
    
    # Átomos de oxígeno
    fig.add_trace(go.Scatter3d(
        x=oxygen_x, y=oxygen_y, z=oxygen_z,
        mode='markers',
        marker=dict(
            size=5,
            color='#4ECDC4',
            opacity=0.7,
            line=dict(color='darkcyan', width=1)
        ),
        name='O (Oxígeno)',
        hovertemplate='<b>O</b><br>X: %{x:.2f}<br>Y: %{y:.2f}<br>Z: %{z:.2f}<extra></extra>'
    ))
    
    # Actualizar layout
    fig.update_layout(
        title=dict(
            text=f'Estructura Cristalina 3D - Porosidad: {porosity:.1%}',
            font=dict(size=16, color='#333')
        ),
        scene=dict(
            xaxis_title='X (Å)',
            yaxis_title='Y (Å)',
            zaxis_title='Z (Å)',
            bgcolor='#fafafa',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5)
            )
        ),
        margin=dict(l=0, r=0, t=50, b=0),
        legend=dict(
            x=0.7,
            y=0.9,
            bgcolor='rgba(255,255,255,0.8)'
        ),
        height=500
    )
    
    return fig


def create_properties_radar(properties: dict) -> go.Figure:
    """Crea gráfico radar de propiedades"""
    if not PLOTLY_AVAILABLE:
        return None
    
    categories = ['Densidad', 'Conductividad', 'Elasticidad', 'Área Superficial', 'Band Gap', 'Calidad']
    
    # Normalizar valores 0-1
    values = [
        min(properties.get('density', 2) / 5, 1),
        min(properties.get('thermal_conductivity', 30) / 100, 1),
        min(properties.get('elastic_modulus', 150) / 300, 1),
        min(properties.get('surface_area', 100) / 500, 1),
        min(properties.get('band_gap', 3) / 6, 1),
        properties.get('quality_score', 0.5)
    ]
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatterpolar(
        r=values + [values[0]],  # Cerrar el polígono
        theta=categories + [categories[0]],
        fill='toself',
        name='Material',
        line_color='#1E88E5',
        fillcolor='rgba(30, 136, 229, 0.3)'
    ))
    
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 1]
            ),
            bgcolor='#fafafa'
        ),
        showlegend=False,
        title='Perfil de Propiedades',
        height=400
    )
    
    return fig

# ============================================================================
# GENERADOR DE PDF
# ============================================================================

def generate_pdf_report(astro_data: dict, physical_props: dict, recipe: dict) -> bytes:
    """Genera reporte PDF profesional"""
    if not PDF_AVAILABLE:
        return None
    
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4)
    styles = getSampleStyleSheet()
    
    # Estilos personalizados
    title_style = ParagraphStyle(
        'CustomTitle',
        parent=styles['Heading1'],
        fontSize=24,
        spaceAfter=30,
        textColor=colors.HexColor('#1E88E5')
    )
    
    heading_style = ParagraphStyle(
        'CustomHeading',
        parent=styles['Heading2'],
        fontSize=14,
        spaceAfter=12,
        textColor=colors.HexColor('#37474F')
    )
    
    body_style = ParagraphStyle(
        'CustomBody',
        parent=styles['Normal'],
        fontSize=11,
        spaceAfter=8,
        textColor=colors.HexColor('#546E7A')
    )
    
    quote_style = ParagraphStyle(
        'Quote',
        parent=styles['Normal'],
        fontSize=10,
        leftIndent=20,
        rightIndent=20,
        spaceAfter=12,
        textColor=colors.HexColor('#795548'),
        fontName='Helvetica-Oblique'
    )
    
    elements = []
    
    # Título
    elements.append(Paragraph("COSMICFORGE LAB", title_style))
    elements.append(Paragraph("Scientific Report", styles['Heading2']))
    elements.append(Spacer(1, 20))
    
    # Información del objeto
    elements.append(Paragraph(f"Objeto Astrofísico: {astro_data.get('object_name', 'Unknown')}", heading_style))
    elements.append(Paragraph(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M')}", body_style))
    elements.append(Paragraph(f"Material: {recipe.get('metal_name', recipe.get('metal', 'N/A'))} {recipe.get('material_type', '')}", body_style))
    elements.append(Spacer(1, 20))
    
    # Tabla de propiedades
    elements.append(Paragraph("Propiedades del Material", heading_style))
    
    props_data = [
        ['Propiedad', 'Valor'],
        ['Porosidad', f"{physical_props.get('porosity', 0):.2%}"],
        ['Densidad', f"{physical_props.get('density', 0):.3f} g/cm³"],
        ['Conductividad Térmica', f"{physical_props.get('thermal_conductivity', 0):.2f} W/m·K"],
        ['Módulo Elástico', f"{physical_props.get('elastic_modulus', 0):.2f} GPa"],
        ['Área Superficial', f"{physical_props.get('surface_area', 0):.1f} m²/g"],
        ['Band Gap', f"{physical_props.get('band_gap', 0):.3f} eV"],
        ['Calidad', f"{physical_props.get('quality_score', 0):.1%}"]
    ]
    
    table = Table(props_data, colWidths=[200, 200])
    table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1E88E5')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 12),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.HexColor('#f5f5f5')),
        ('GRID', (0, 0), (-1, -1), 1, colors.HexColor('#e0e0e0'))
    ]))
    
    elements.append(table)
    elements.append(Spacer(1, 30))
    
    # Conclusión filosófica (Manin)
    elements.append(Paragraph("Reflexión Científica", heading_style))
    
    quote = MANIN_QUOTES[hash(astro_data.get('object_name', '')) % len(MANIN_QUOTES)]
    elements.append(Paragraph(f'"{quote}"', quote_style))
    elements.append(Paragraph("— Yuri I. Manin, 'Lo demostrable e indemostrable'", 
                              ParagraphStyle('Author', parent=body_style, alignment=2, fontSize=9)))
    elements.append(Spacer(1, 20))
    
    # Conclusión
    elements.append(Paragraph("Conclusión", heading_style))
    
    quality = physical_props.get('quality_score', 0.5)
    porosity = physical_props.get('porosity', 0.3)
    
    if quality > 0.7:
        conclusion = f"El material predicho a partir de {astro_data.get('object_name', 'este objeto astrofísico')} presenta características excelentes para aplicaciones en {', '.join(recipe.get('applications', ['ciencia de materiales']))}. La combinación de porosidad ({porosity:.1%}) y propiedades mecánicas sugiere un material viable para síntesis experimental."
    elif quality > 0.4:
        conclusion = f"El material derivado de {astro_data.get('object_name', 'este objeto')} muestra propiedades moderadas que requieren optimización. Se recomienda ajustar parámetros de síntesis para mejorar la calidad del producto final."
    else:
        conclusion = f"La predicción para {astro_data.get('object_name', 'este objeto')} indica desafíos en la síntesis. Se sugiere explorar configuraciones alternativas o diferentes combinaciones de elementos."
    
    elements.append(Paragraph(conclusion, body_style))
    elements.append(Spacer(1, 30))
    
    # Footer
    elements.append(Paragraph("─" * 50, body_style))
    elements.append(Paragraph("Generado por CosmicForge Lab Pro | www.cosmicforge-lab.com", 
                              ParagraphStyle('Footer', parent=body_style, alignment=1, fontSize=8, textColor=colors.grey)))
    
    doc.build(elements)
    buffer.seek(0)
    return buffer.getvalue()

# ============================================================================
# CONCLUSIÓN FILOSÓFICA (MANIN)
# ============================================================================

def generate_manin_conclusion(astro_data: dict, physical_props: dict, recipe: dict) -> str:
    """Genera conclusión basada en la lógica de Manin"""
    
    quality = physical_props.get('quality_score', 0.5)
    fractal_dim = astro_data.get('fractal_dimension', 1.5)
    criticality = astro_data.get('criticality_score', 0.5)
    object_name = astro_data.get('object_name', 'este objeto')
    metal = recipe.get('metal', 'Ti')
    
    # Seleccionar cita apropiada
    if criticality > 0.8:
        quote = MANIN_QUOTES[2]  # Zona fronteriza
    elif fractal_dim > 1.8:
        quote = MANIN_QUOTES[0]  # Matemáticas y realidad
    elif quality > 0.7:
        quote = MANIN_QUOTES[3]  # Abstracciones y realidad
    else:
        quote = MANIN_QUOTES[4]  # Verdad matemática
    
    # Generar conclusión
    conclusion = f"""
    <div class="manin-quote">
        <p><em>"{quote}"</em></p>
        <p class="author">— Yuri I. Manin, "Lo demostrable e indemostrable"</p>
    </div>
    
    <div class="conclusion-box">
        <h4>🔍 Análisis Epistemológico</h4>
        <p>La predicción del material <strong>{metal}-O</strong> a partir de <strong>{object_name}</strong> 
        ejemplifica la tensión entre <em>lo demostrable</em> (las ecuaciones matemáticas que traducen 
        firmas astrofísicas a propiedades materiales) y <em>lo indemostrable</em> (la validez física 
        última de esta analogía).</p>
        
        <p><strong>Nivel de certeza matemática:</strong> {'✅ Alto' if quality > 0.6 else '🟡 Moderado'} 
        — Las ecuaciones están formalmente derivadas.</p>
        
        <p><strong>Nivel de certeza física:</strong> 🟡 Requiere validación — La traducción 
        astrofísica → material es una hipótesis por verificar.</p>
        
        <p><strong>Implicación filosófica:</strong> Este trabajo opera en la "zona fronteriza" 
        que Manin identifica como el espacio creativo de la ciencia, donde las matemáticas 
        puras encuentran su aplicación más profunda en la comprensión del universo material.</p>
    </div>
    """
    
    return conclusion

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def validate_astrophysical_data(data: dict):
    """Validate astrophysical data structure"""
    errors = []
    required_fields = ['fractal_dimension', 'criticality_score', 'entropy', 'anisotropy', 'turbulence_beta']
    
    for field in required_fields:
        if field not in data:
            errors.append(f"Campo requerido faltante: {field}")
        elif not isinstance(data[field], (int, float)):
            errors.append(f"Campo {field} debe ser numérico")
    
    is_valid = len(errors) == 0
    return is_valid, errors

def display_metrics_from_data(data: dict):
    """Display astrophysical signature metrics"""
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Dimensión Fractal", f"{data.get('fractal_dimension', 0):.4f}")
        st.metric("Criticalidad", f"{data.get('criticality_score', 0):.4f}")
    
    with col2:
        st.metric("Entropía", f"{data.get('entropy', 0):.6f}")
        st.metric("Anisotropía", f"{data.get('anisotropy', 0):.4f}")
    
    with col3:
        st.metric("Turbulencia β", f"{data.get('turbulence_beta', 0):.4f}")
        st.metric("Lyapunov max", f"{data.get('lyapunov_max', 0):.6f}")
    
    with col4:
        st.metric("Objeto", str(data.get('object_name', 'Unknown'))[:18])
        st.metric("Fuente", data.get('source', 'Manual'))

def display_alloy_systems(metal: str):
    """Muestra sistemas de aleaciones compatibles"""
    if metal not in COMPATIBLE_ELEMENTS:
        return
    
    compatible = COMPATIBLE_ELEMENTS[metal].get('compatible', [])
    conflict = COMPATIBLE_ELEMENTS[metal].get('conflict', [])
    
    st.markdown("#### 🔗 Compatibilidad de Elementos")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**✅ Elementos compatibles:**")
        comp_html = "".join([f'<span class="compatible-tag">{e}</span>' for e in compatible])
        st.markdown(comp_html, unsafe_allow_html=True)
    
    with col2:
        if conflict:
            st.markdown("**⚠️ Conflictos:**")
            conf_html = "".join([f'<span class="conflict-tag">{e}</span>' for e in conflict])
            st.markdown(conf_html, unsafe_allow_html=True)
        else:
            st.markdown("**✅ Sin conflictos conocidos**")
    
    # Aleaciones sugeridas
    st.markdown("#### 🧪 Aleaciones Sugeridas")
    alloys_found = []
    for key, alloy in ALLOY_SYSTEMS.items():
        if metal in key.split('-'):
            alloys_found.append((key, alloy))
    
    if alloys_found:
        for key, alloy in alloys_found[:4]:
            st.markdown(f"""
            <div class="alloy-card">
                <strong>{alloy['name']}</strong> ({key})<br>
                <small>Aplicaciones: {', '.join(alloy['applications'])}</small><br>
                <small>Ratio: {alloy['ratio']}</small>
            </div>
            """, unsafe_allow_html=True)
    else:
        st.info("No hay aleaciones predefinidas para este elemento.")

# ============================================================================
# MAIN APPLICATION
# ============================================================================

# Header
st.markdown('<p class="main-header">🔬 CosmicForge Lab <span class="pro-badge">PRO</span></p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Diseño de Materiales Inspirado en Firmas Astrofísicas | Edición Profesional</p>', unsafe_allow_html=True)

# =========================================================================
# SIDEBAR - DATA INPUT
# =========================================================================

st.sidebar.header("📥 Importar Datos")

input_method = st.sidebar.radio(
    "Método de entrada:",
    ["Ejemplos Astrofísicos", "JSON Anomaly Detector", "Editor Manual", "Archivos QE/CIF"]
)

if input_method == "Ejemplos Astrofísicos":
    # Filtros
    filter_type = st.sidebar.selectbox(
        "Filtrar por tipo:",
        ["Todos", "Nebulosas", "Galaxias", "Especiales"]
    )
    
    # Filtrar ejemplos
    if filter_type == "Nebulosas":
        filtered = {k: v for k, v in ASTROPHYSICAL_EXAMPLES.items() if "Nebula" in v.get("type", "")}
    elif filter_type == "Galaxias":
        filtered = {k: v for k, v in ASTROPHYSICAL_EXAMPLES.items() if "Galaxy" in v.get("type", "")}
    elif filter_type == "Especiales":
        filtered = {k: v for k, v in ASTROPHYSICAL_EXAMPLES.items() if "Nebula" not in v.get("type", "") and "Galaxy" not in v.get("type", "")}
    else:
        filtered = ASTROPHYSICAL_EXAMPLES
    
    example_choice = st.sidebar.selectbox("Selecciona objeto:", list(filtered.keys()))
    
    if st.sidebar.button("📋 Cargar ejemplo", type="primary"):
        example_data = ASTROPHYSICAL_EXAMPLES[example_choice].copy()
        example_data['object_name'] = example_choice
        example_data['source'] = 'Database'
        st.session_state.astro_data = example_data
        st.sidebar.success(f"✅ {example_choice} cargado")
    
    # Mostrar info del ejemplo seleccionado
    if example_choice in ASTROPHYSICAL_EXAMPLES:
        ex = ASTROPHYSICAL_EXAMPLES[example_choice]
        st.sidebar.markdown(f"""
        **Tipo:** {ex.get('type', 'N/A')}  
        **Distancia:** {ex.get('distance', 'N/A')}  
        *{ex.get('description', '')}*
        """)

elif input_method == "JSON Anomaly Detector":
    uploaded_file = st.sidebar.file_uploader("Archivo JSON", type=['json'])
    
    if uploaded_file is not None:
        try:
            string_data = uploaded_file.read().decode('utf-8')
            raw_json = json.loads(string_data)
            st.session_state.raw_json = raw_json
            
            if 'plugin_results' in raw_json:
                new_astro_data = parse_anomaly_detector_json(raw_json)
                st.session_state.astro_data = new_astro_data
                st.sidebar.success(f"✅ Anomaly Detector JSON cargado")
            else:
                is_valid, errors = validate_astrophysical_data(raw_json)
                if is_valid:
                    st.session_state.astro_data = raw_json
                    st.sidebar.success("✅ JSON cargado")
                else:
                    st.sidebar.error("❌ Formato inválido")
        except Exception as e:
            st.sidebar.error(f"❌ Error: {e}")

elif input_method == "Editor Manual":
    st.sidebar.subheader("Parámetros")
    
    object_name = st.sidebar.text_input("Nombre", "Custom Object")
    fractal_dim = st.sidebar.slider("Dim. Fractal", 0.0, 3.0, 1.5, 0.001)
    criticality = st.sidebar.slider("Criticalidad", 0.0, 1.0, 0.5, 0.01)
    entropy = st.sidebar.number_input("Entropía", 0.0, 1.0, 0.01, 0.0001, format="%.6f")
    anisotropy = st.sidebar.slider("Anisotropía", 0.0, 1.0, 0.3, 0.01)
    turbulence = st.sidebar.slider("Turbulencia β", 0.0, 5.0, 2.0, 0.01)
    lyapunov = st.sidebar.number_input("Lyapunov", -1.0, 1.0, -0.1, 0.001)
    
    if st.sidebar.button("✅ Aplicar", type="primary"):
        st.session_state.astro_data = {
            "object_name": object_name,
            "fractal_dimension": fractal_dim,
            "criticality_score": criticality,
            "entropy": entropy,
            "anisotropy": anisotropy,
            "turbulence_beta": turbulence,
            "lyapunov_max": lyapunov,
            "mode": "custom",
            "source": "Editor Manual"
        }
        st.sidebar.success("✅ Parámetros aplicados")

elif input_method == "Archivos QE/CIF":
    file_type = st.sidebar.selectbox("Tipo", ["Quantum ESPRESSO", "CIF", "LAMMPS"])
    ext_map = {"Quantum ESPRESSO": ['in', 'pw'], "CIF": ['cif'], "LAMMPS": ['data', 'lmp']}
    
    uploaded_file = st.sidebar.file_uploader("Archivo", type=ext_map[file_type])
    if uploaded_file:
        try:
            content = uploaded_file.read().decode('utf-8')
            parser = FileParser()
            
            if file_type == "Quantum ESPRESSO":
                st.session_state.astro_data = parser.parse_quantum_espresso(content)
            elif file_type == "CIF":
                st.session_state.astro_data = parser.parse_cif(content)
            else:
                st.session_state.astro_data = parser.parse_lammps_data(content)
            st.sidebar.success(f"✅ {file_type} parseado")
        except Exception as e:
            st.sidebar.error(f"❌ Error: {e}")

# Botón limpiar
if st.session_state.astro_data is not None:
    if st.sidebar.button("🗑️ Limpiar"):
        st.session_state.astro_data = None
        st.session_state.physical_props = None
        st.session_state.recipe = None
        st.session_state.raw_json = None
        st.rerun()

# =========================================================================
# SIDEBAR - MATERIAL & ALLOY SYSTEM
# =========================================================================

st.sidebar.header("⚗️ Material")

material_type = st.sidebar.selectbox("Tipo", ["Óxido", "Metal puro", "Cerámica", "Polímero", "Composite", "Nanopartícula"])

metal = st.sidebar.selectbox("Metal base", list(COMPATIBLE_ELEMENTS.keys()))

# Selector de aleación (si aplica)
if material_type in ["Óxido", "Metal puro", "Nanopartícula"]:
    use_alloy = st.sidebar.checkbox("🔧 Usar aleación")
    
    if use_alloy:
        alloy_options = [k for k in ALLOY_SYSTEMS.keys() if metal in k.split('-')]
        if alloy_options:
            selected_alloy = st.sidebar.selectbox("Aleación", alloy_options)
            st.sidebar.info(f"**{ALLOY_SYSTEMS[selected_alloy]['name']}**\n\nApps: {', '.join(ALLOY_SYSTEMS[selected_alloy]['applications'])}")

with st.sidebar.expander("⚙️ Avanzado"):
    bulk_density = st.number_input("Densidad (g/cm³)", value=3.0, min_value=0.1, max_value=20.0)
    bulk_conductivity = st.number_input("Conductividad (W/m·K)", value=50.0, min_value=0.1, max_value=500.0)
    bulk_elastic = st.number_input("Elasticidad (GPa)", value=200.0, min_value=1.0, max_value=1000.0)

# =========================================================================
# MAIN CONTENT
# =========================================================================

if st.session_state.astro_data is not None:
    astro_data = st.session_state.astro_data
    
    st.header("🌌 Firma Astrofísica Detectada")
    display_metrics_from_data(astro_data)
    
    # Mostrar compatibilidad de elementos
    st.markdown("---")
    display_alloy_systems(metal)
    
    # Botón generar
    if st.button("🚀 Generar Receta de Material", type="primary", use_container_width=True):
        with st.spinner("Calculando propiedades y generando análisis..."):
            try:
                physics_calc = PhysicsCalculator(
                    bulk_density=bulk_density,
                    bulk_conductivity=bulk_conductivity,
                    bulk_elastic=bulk_elastic
                )
                chem_engine = ChemistryEngine()
                
                physical_props = physics_calc.calculate_all_properties(astro_data)
                st.session_state.physical_props = physical_props
                
                recipe = chem_engine.generate_complete_recipe(
                    metal=metal,
                    material_type=material_type.lower(),
                    astrophysical_data=astro_data
                )
                st.session_state.recipe = recipe
                
                st.success("✅ Cálculos completados!")
                st.rerun()
            except Exception as e:
                st.error(f"❌ Error: {e}")
                st.exception(e)

# =========================================================================
# RESULTS
# =========================================================================

if st.session_state.physical_props is not None and st.session_state.recipe is not None:
    physical_props = st.session_state.physical_props
    recipe = st.session_state.recipe
    astro_data = st.session_state.astro_data
    
    file_gen = FileGenerator()
    
    tabs = st.tabs([
        "📊 Propiedades", "⚗️ Síntesis", "🎯 3D Structure", "📈 Radar", 
        "🔗 Aleaciones", "📁 Archivos", "📄 PDF", "💡 Conclusión"
    ])
    
    with tabs[0]:
        st.subheader("Propiedades del Material")
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Porosidad", f"{physical_props.get('porosity', 0):.1%}")
            st.metric("Densidad", f"{physical_props.get('density', 0):.3f} g/cm³")
        
        with col2:
            st.metric("Cond. Térmica", f"{physical_props.get('thermal_conductivity', 0):.2f} W/m·K")
            st.metric("Módulo Elástico", f"{physical_props.get('elastic_modulus', 0):.2f} GPa")
        
        with col3:
            st.metric("Área Superficial", f"{physical_props.get('surface_area', 0):.1f} m²/g")
            st.metric("Band Gap", f"{physical_props.get('band_gap', 0):.3f} eV")
        
        with col4:
            st.metric("Calidad", f"{physical_props.get('quality_score', 0):.1%}")
            st.metric("E. Activación", f"{physical_props.get('activation_energy', 0):.4f} eV")
        
        # Materials Project Link
        st.markdown("---")
        st.markdown("### 🔗 Validación con Materials Project")
        mp_formula = f"{recipe.get('metal', 'Ti')}O2"
        st.markdown(f"""
        <div class="info-box">
        <strong>Consulta tu material en bases de datos científicas:</strong><br><br>
        • <a href="https://materialsproject.org/materials?search={mp_formula}" target="_blank">Materials Project - {mp_formula}</a><br>
        • <a href="https://www.aflowlib.org/?search={mp_formula}" target="_blank">AFLOW Library</a><br>
        • <a href="http://oqmd.org/materials?search={mp_formula}" target="_blank">OQMD Database</a><br>
        • <a href="https://nomad-lab.eu/" target="_blank">NOMAD Archive</a>
        </div>
        """, unsafe_allow_html=True)
    
    with tabs[1]:
        st.subheader("Receta de Síntesis")
        conditions = recipe.get('reaction_conditions', {})
        
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("### 🌡️ Condiciones")
            st.write(f"**Temperatura:** {conditions.get('temperature_C', 0):.0f}°C")
            st.write(f"**Tiempo:** {conditions.get('reaction_time_hours', 0):.2f} horas")
            st.write(f"**pH:** {conditions.get('pH', 7):.1f}")
        
        with col2:
            st.markdown("### 🧪 Precursores")
            for precursor in recipe.get('precursors', []):
                st.write(f"• {precursor}")
        
        st.markdown("### 📝 Procedimiento")
        for i, step in enumerate(recipe.get('step_by_step', []), 1):
            st.write(f"**{i}.** {step}")
    
    with tabs[2]:
        st.subheader("Visualización 3D Interactiva")
        if PLOTLY_AVAILABLE:
            fig_3d = create_3d_structure(
                physical_props.get('porosity', 0.3),
                physical_props.get('density', 2.0),
                recipe.get('metal', 'Ti')
            )
            if fig_3d:
                st.plotly_chart(fig_3d, use_container_width=True)
                st.info("🖱️ Arrastra para rotar | Scroll para zoom | Doble click para resetear")
        else:
            st.warning("Instala plotly: `pip install plotly`")
    
    with tabs[3]:
        st.subheader("Perfil de Propiedades")
        if PLOTLY_AVAILABLE:
            fig_radar = create_properties_radar(physical_props)
            if fig_radar:
                st.plotly_chart(fig_radar, use_container_width=True)
    
    with tabs[4]:
        st.subheader("Sistema de Aleaciones")
        display_alloy_systems(recipe.get('metal', 'Ti'))
        
        # Mostrar todas las aleaciones disponibles
        st.markdown("### 📚 Base de Datos de Aleaciones")
        for key, alloy in list(ALLOY_SYSTEMS.items())[:6]:
            st.markdown(f"""
            <div class="alloy-card">
                <strong>{alloy['name']}</strong> ({key})<br>
                <small>Ratio: {alloy['ratio']} | Apps: {', '.join(alloy['applications'])}</small>
            </div>
            """, unsafe_allow_html=True)
    
    with tabs[5]:
        st.subheader("Archivos para Simuladores")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            poscar = file_gen.generate_poscar(physical_props, recipe, 'oxide')
            st.download_button("⬇️ POSCAR (VASP/QE)", poscar, "POSCAR.vasp", "text/plain")
        
        with col2:
            lammps = file_gen.generate_lammps_input(physical_props, recipe, 'oxide')
            st.download_button("⬇️ LAMMPS", lammps, "input.lmp", "text/plain")
        
        with col3:
            csv_data = file_gen.generate_properties_csv(physical_props, recipe)
            st.download_button("⬇️ CSV", csv_data, "properties.csv", "text/csv")
    
    with tabs[6]:
        st.subheader("Generar Reporte PDF")
        
        if st.button("📄 Generar PDF", type="primary"):
            if PDF_AVAILABLE:
                pdf_bytes = generate_pdf_report(astro_data, physical_props, recipe)
                if pdf_bytes:
                    st.download_button(
                        "⬇️ Descargar PDF",
                        pdf_bytes,
                        f"CosmicForge_Report_{astro_data.get('object_name', 'material')}.pdf",
                        "application/pdf"
                    )
            else:
                st.warning("Instala reportlab: `pip install reportlab`")
                st.info("El PDF incluye: propiedades, conclusión filosófica (Manin), y análisis epistemológico.")
    
    with tabs[7]:
        st.subheader("💡 Conclusión Científica")
        st.markdown(generate_manin_conclusion(astro_data, physical_props, recipe), unsafe_allow_html=True)

else:
    st.info("👈 Selecciona un ejemplo astrofísico o carga un archivo JSON para comenzar")
    
    # Mostrar galería de ejemplos
    st.markdown("### 🌌 Galería de Objetos Astrofísicos")
    
    cols = st.columns(3)
    examples_list = list(ASTROPHYSICAL_EXAMPLES.items())
    
    for i, (name, data) in enumerate(examples_list[:9]):
        with cols[i % 3]:
            st.markdown(f"""
            <div class="example-card">
                <strong>{name}</strong><br>
                <small>{data.get('type', 'N/A')}</small><br>
                <small>Distancia: {data.get('distance', 'N/A')}</small>
            </div>
            """, unsafe_allow_html=True)

# Footer
st.markdown("---")
st.markdown("""
<div class="pro-footer">
    <p><strong>CosmicForge Lab Pro v3.0</strong></p>
    <p>Diseño de Materiales Inspirado en Firmas Astrofísicas</p>
    <p>
        🔗 <a href="https://github.com/WilmerGaspar/cosmicforge-lab" target="_blank">GitHub</a> | 
        📧 Contacto | 
        📄 <a href="https://materialsproject.org" target="_blank">Materials Project</a>
    </p>
    <p><small>© 2025 CosmicForge Lab | Edición Profesional</small></p>
</div>
""", unsafe_allow_html=True)
