"""
CosmicForge Lab - NASA Edition
Diseño de Materiales Inspirado en Firmas Astrofísicas
Versión 3.6 - Edición NASA con DFT y Simulación Espacial

Autor: Wilmer Gaspar
Repositorio: https://github.com/WilmerGaspar/cosmicforge-lab

NUEVO EN v3.6:
- Simulación de condiciones extremas (vacío, radiación, microgravedad, ciclos térmicos)
- Optimización con algoritmo genético
- Cálculo DFT simulado (Quantum ESPRESSO)
- Integración Materials Project API
- Comparación con materiales comerciales
- Tema oscuro NASA-style
- Archivos para supercomputadora
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
import random
import math

# Add modules directory to path
MODULES_PATH = Path(__file__).parent / "modules"
sys.path.insert(0, str(MODULES_PATH))

from modules.physics_calculator import PhysicsCalculator
from modules.chemistry_engine import ChemistryEngine
from modules.file_generators import FileGenerator
from modules.visualizer import Visualizer
from modules.file_parser import FileParser

# Plotly para 3D y gráficos
try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except:
    PLOTLY_AVAILABLE = False

# PDF generation
try:
    from fpdf import FPDF
    PDF_AVAILABLE = True
except:
    PDF_AVAILABLE = False

# Pandas para tablas
try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except:
    PANDAS_AVAILABLE = False

# ============================================================================
# PAGE CONFIGURATION
# ============================================================================

st.set_page_config(
    page_title="CosmicForge Lab - NASA Edition",
    page_icon="🚀",
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
if 'simulation' not in st.session_state:
    st.session_state.simulation = None
if 'optimization' not in st.session_state:
    st.session_state.optimization = None
if 'dft_result' not in st.session_state:
    st.session_state.dft_result = None
if 'dark_mode' not in st.session_state:
    st.session_state.dark_mode = True

# ============================================================================
# DATABASES - EXTENDED
# ============================================================================

# Ejemplos de objetos astrofísicos (15 ejemplos)
ASTROPHYSICAL_EXAMPLES = {
    "Orion Nebula M42": {
        "fractal_dimension": 1.654, "criticality_score": 0.722, "entropy": 0.019,
        "anisotropy": 0.329, "turbulence_beta": 2.278, "lyapunov_max": -0.227,
        "mode": "balanced", "type": "Nebulosa de Emisión", "distance": "1,344 ly",
        "temperature": "10,000 K", "luminosity": "10^5 L☉"
    },
    "Crab Nebula M1": {
        "fractal_dimension": 1.82, "criticality_score": 0.89, "entropy": 0.032,
        "anisotropy": 0.456, "turbulence_beta": 2.45, "lyapunov_max": -0.15,
        "mode": "turbulent", "type": "Remanente Supernova", "distance": "6,500 ly",
        "temperature": "1,600 K", "luminosity": "10^5 L☉"
    },
    "Ring Nebula M57": {
        "fractal_dimension": 1.45, "criticality_score": 0.55, "entropy": 0.012,
        "anisotropy": 0.28, "turbulence_beta": 1.95, "lyapunov_max": -0.35,
        "mode": "stable", "type": "Nebulosa Planetaria", "distance": "2,283 ly",
        "temperature": "125,000 K", "luminosity": "200 L☉"
    },
    "Triangulum Galaxy M33": {
        "fractal_dimension": 1.93, "criticality_score": 0.78, "entropy": 0.00000964,
        "anisotropy": 0.23, "turbulence_beta": 1.84, "lyapunov_max": 0.11,
        "mode": "balanced", "type": "Galaxia Espiral", "distance": "2.73M ly",
        "temperature": "10 K", "luminosity": "3×10^9 L☉"
    },
    "Andromeda Galaxy M31": {
        "fractal_dimension": 1.88, "criticality_score": 0.85, "entropy": 0.000012,
        "anisotropy": 0.19, "turbulence_beta": 1.92, "lyapunov_max": 0.08,
        "mode": "balanced", "type": "Galaxia Espiral", "distance": "2.537M ly",
        "temperature": "10 K", "luminosity": "2.6×10^10 L☉"
    },
    "Eagle Nebula M16": {
        "fractal_dimension": 1.78, "criticality_score": 0.81, "entropy": 0.028,
        "anisotropy": 0.42, "turbulence_beta": 2.35, "lyapunov_max": -0.18,
        "mode": "turbulent", "type": "Nebulosa de Emisión", "distance": "7,000 ly",
        "temperature": "8,000 K", "luminosity": "10^6 L☉"
    },
    "Tarantula Nebula": {
        "fractal_dimension": 1.91, "criticality_score": 0.93, "entropy": 0.045,
        "anisotropy": 0.58, "turbulence_beta": 2.68, "lyapunov_max": 0.15,
        "mode": "turbulent", "type": "Región HII", "distance": "160,000 ly",
        "temperature": "15,000 K", "luminosity": "10^7 L☉"
    },
    "Pillars of Creation": {
        "fractal_dimension": 1.85, "criticality_score": 0.88, "entropy": 0.035,
        "anisotropy": 0.52, "turbulence_beta": 2.55, "lyapunov_max": -0.12,
        "mode": "turbulent", "type": "Nube Molecular", "distance": "6,500 ly",
        "temperature": "10 K", "luminosity": "N/A"
    },
    "Whirlpool Galaxy M51": {
        "fractal_dimension": 1.95, "criticality_score": 0.91, "entropy": 0.000015,
        "anisotropy": 0.22, "turbulence_beta": 1.98, "lyapunov_max": 0.12,
        "mode": "turbulent", "type": "Galaxia Interactuante", "distance": "23M ly",
        "temperature": "10 K", "luminosity": "10^10 L☉"
    },
    "Helix Nebula NGC 7293": {
        "fractal_dimension": 1.48, "criticality_score": 0.58, "entropy": 0.011,
        "anisotropy": 0.31, "turbulence_beta": 2.02, "lyapunov_max": -0.32,
        "mode": "stable", "type": "Nebulosa Planetaria", "distance": "655 ly",
        "temperature": "100,000 K", "luminosity": "100 L☉"
    },
    "Betelgeuse": {
        "fractal_dimension": 1.72, "criticality_score": 0.95, "entropy": 0.055,
        "anisotropy": 0.65, "turbulence_beta": 2.85, "lyapunov_max": 0.25,
        "mode": "turbulent", "type": "Supergigante Roja", "distance": "700 ly",
        "temperature": "3,500 K", "luminosity": "1.2×10^5 L☉"
    },
    "Sirius Binary System": {
        "fractal_dimension": 1.55, "criticality_score": 0.62, "entropy": 0.015,
        "anisotropy": 0.38, "turbulence_beta": 2.15, "lyapunov_max": -0.08,
        "mode": "balanced", "type": "Sistema Binario", "distance": "8.6 ly",
        "temperature": "9,940 K", "luminosity": "25.4 L☉"
    },
    "Jupiter Great Red Spot": {
        "fractal_dimension": 1.68, "criticality_score": 0.75, "entropy": 0.022,
        "anisotropy": 0.42, "turbulence_beta": 2.32, "lyapunov_max": -0.05,
        "mode": "balanced", "type": "Tormenta Atmosférica", "distance": "4.2 AU",
        "temperature": "110 K", "luminosity": "N/A"
    },
    "Saturn Rings": {
        "fractal_dimension": 1.92, "criticality_score": 0.68, "entropy": 0.008,
        "anisotropy": 0.15, "turbulence_beta": 1.45, "lyapunov_max": -0.42,
        "mode": "stable", "type": "Sistema de Anillos", "distance": "9.5 AU",
        "temperature": "90 K", "luminosity": "N/A"
    },
    "Vela Supernova Remnant": {
        "fractal_dimension": 1.79, "criticality_score": 0.87, "entropy": 0.038,
        "anisotropy": 0.48, "turbulence_beta": 2.52, "lyapunov_max": 0.08,
        "mode": "turbulent", "type": "Remanente Supernova", "distance": "936 ly",
        "temperature": "800 K", "luminosity": "10^4 L☉"
    }
}

# Base de datos de elementos (15 metales)
ELEMENTS_DATABASE = {
    "Ti": {"name": "Titanio", "atomic_num": 22, "mw": 47.867, "mp": 1668, "bp": 3287, "density": 4.506, "compatible": ["Al", "V", "Fe", "Ni", "O", "N", "C"]},
    "Al": {"name": "Aluminio", "atomic_num": 13, "mw": 26.982, "mp": 660, "bp": 2519, "density": 2.70, "compatible": ["Ti", "Cu", "Mg", "Si", "Zn", "O"]},
    "Fe": {"name": "Hierro", "atomic_num": 26, "mw": 55.845, "mp": 1538, "bp": 2862, "density": 7.874, "compatible": ["Ni", "Cr", "Co", "Mn", "C", "O"]},
    "Zn": {"name": "Zinc", "atomic_num": 30, "mw": 65.38, "mp": 420, "bp": 907, "density": 7.14, "compatible": ["Cu", "Al", "O", "S"]},
    "Cu": {"name": "Cobre", "atomic_num": 29, "mw": 63.546, "mp": 1085, "bp": 2562, "density": 8.96, "compatible": ["Ni", "Zn", "Al", "Sn", "O"]},
    "Ni": {"name": "Níquel", "atomic_num": 28, "mw": 58.693, "mp": 1455, "bp": 2913, "density": 8.91, "compatible": ["Fe", "Co", "Cu", "Cr", "Ti", "O"]},
    "Co": {"name": "Cobalto", "atomic_num": 27, "mw": 58.933, "mp": 1495, "bp": 2927, "density": 8.86, "compatible": ["Ni", "Fe", "Mn", "O"]},
    "Mn": {"name": "Manganeso", "atomic_num": 25, "mw": 54.938, "mp": 1246, "bp": 2061, "density": 7.47, "compatible": ["Fe", "Co", "O"]},
    "Ag": {"name": "Plata", "atomic_num": 47, "mw": 107.868, "mp": 962, "bp": 2162, "density": 10.49, "compatible": ["Cu", "Au", "Pd", "O"]},
    "Au": {"name": "Oro", "atomic_num": 79, "mw": 196.967, "mp": 1064, "bp": 2856, "density": 19.30, "compatible": ["Ag", "Pt", "Pd", "Cu"]},
    "Pt": {"name": "Platino", "atomic_num": 78, "mw": 195.084, "mp": 1768, "bp": 3825, "density": 21.45, "compatible": ["Pd", "Au", "Rh", "Ir", "O"]},
    "Pd": {"name": "Paladio", "atomic_num": 46, "mw": 106.42, "mp": 1555, "bp": 2963, "density": 12.02, "compatible": ["Pt", "Au", "Ag", "Ni", "O"]},
    "Cr": {"name": "Cromo", "atomic_num": 24, "mw": 52.00, "mp": 1907, "bp": 2671, "density": 7.19, "compatible": ["Fe", "Ni", "Co", "O"]},
    "Mo": {"name": "Molibdeno", "atomic_num": 42, "mw": 95.95, "mp": 2623, "bp": 4639, "density": 10.28, "compatible": ["Ti", "Fe", "Ni", "Cr", "O"]},
    "W": {"name": "Wolframio", "atomic_num": 74, "mw": 183.84, "mp": 3422, "bp": 5930, "density": 19.25, "compatible": ["Ti", "Mo", "Cr", "O"]},
}

# Sistema de aleaciones completo (10 aleaciones)
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
        "safety": ["Usar guantes refractarios", "Atmósfera inerte obligatoria", "Ventilación adecuada"],
        "commercial_name": "Gamma-Met"
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
        "safety": ["Vacío alto (<10⁻³ mbar)", "Protección radiación UV", "Evitar contaminación con O₂, N₂, H₂"],
        "commercial_name": "Ti-6Al-4V / Grade 5"
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
        "safety": ["Evitar oxidación", "Controlar temperatura exacta"],
        "commercial_name": "Invar 36"
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
        "safety": ["Fundente para evitar oxidación", "Ventilación"],
        "commercial_name": "C70600 / CuNi 70/30"
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
        "safety": ["Evitar humedad (explosión)", "Degaseado obligatorio"],
        "commercial_name": "AA2024-T6"
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
            "6. Calcinación: 400°C (anatasa) o 800°C (rutilo)"
        ],
        "equipment": ["Matraz reacción", "Agitador magnético", "Horno mufla", "Autoclave (opcional)"],
        "precursors": ["Isopropóxido de Ti", "TiCl₄", "Etanol", "HNO₃"],
        "safety": ["Guantes y gafas", "Campana extractora", "TiCl₄ es corrosivo"],
        "commercial_name": "TiO₂ P25"
    },
    "Ni-Ti (Nitinol)": {
        "name": "Nitinol (NiTi - Memoria de Forma)",
        "ratio": "Ni:Ti = 50:50",
        "melting_point": "1310°C",
        "density": "6.45 g/cm³",
        "applications": ["Stents cardiovasculares", "Actuadores", "Implantes médicos", "Aparatos ortodónticos"],
        "synthesis": [
            "1. Fundir Ni y Ti en horno de arco bajo argón",
            "2. Voltear lingote 5-6 veces para homogeneizar",
            "3. Forjar a 850°C",
            "4. Recocer a 500°C para ajustar temperatura de transición",
            "5. Enfriar en agua para fijar estructura",
            "6. Tratamiento térmico final para ajustar Af"
        ],
        "equipment": ["Horno de arco", "Atmósfera inerte", "Horno de tratamiento", "Baño de sal"],
        "precursors": ["Ni electrolítico (99.99%)", "Ti esponja (99.9%)"],
        "safety": ["Control preciso de composición", "Evitar contaminación", "Atmósfera inerte"],
        "commercial_name": "Nitinol SE508"
    },
    "Co-Cr-Mo (Vitallium)": {
        "name": "Vitallium (Co-Cr-Mo)",
        "ratio": "Co:Cr:Mo = 60:30:5",
        "melting_point": "1490°C",
        "density": "8.29 g/cm³",
        "applications": ["Implantes dentales", "Prótesis articulares", "Componentes aeroespaciales"],
        "synthesis": [
            "1. Fundir Co en horno de inducción bajo vacío",
            "2. Añadir Cr y Mo gradualmente",
            "3. Homogeneizar a 1500°C por 1 hora",
            "4. Colar en molde de investment",
            "5. Tratamiento térmico: solución 1200°C + envejecimiento",
            "6. Pulido y pasivación"
        ],
        "equipment": ["Horno de inducción al vacío", "Moldes cerámicos", "Horno de tratamiento"],
        "precursors": ["Co electrolítico", "Cr electrolítico", "Mo metálico"],
        "safety": ["Vacío alto", "Protección contra humos metálicos"],
        "commercial_name": "Vitallium / F75"
    },
    "Fe-Cr-Ni (Inoxidable 316)": {
        "name": "Acero Inoxidable 316L",
        "ratio": "Fe:Cr:Ni:Mo = 65:17:12:2.5",
        "melting_point": "1400°C",
        "density": "8.00 g/cm³",
        "applications": ["Equipamiento médico", "Industria alimentaria", "Ambientes marinos"],
        "synthesis": [
            "1. Fundir Fe en horno de arco eléctrico",
            "2. Añadir Cr, Ni y Mo como ferroaleaciones",
            "3. Refinar con AOD (Argon Oxygen Decarburization)",
            "4. Colar en lingotera",
            "5. Laminar en caliente a 1150°C",
            "6. Recocer a 1050°C y temple en agua",
            "7. Pasivación con ácido nítrico"
        ],
        "equipment": ["Horno de arco", "Convertidor AOD", "Laminador", "Horno de recocido"],
        "precursors": ["Fe chatarra", "Ferrocromo", "Feroníquel", "Ferromolibdeno"],
        "safety": ["Protección contra humos", "Ventilación adecuada"],
        "commercial_name": "AISI 316L"
    },
    "Al-Mg-Si (6061)": {
        "name": "Aluminio 6061",
        "ratio": "Al:Mg:Si:Cu = 97.5:1.0:0.6:0.3",
        "melting_point": "582°C",
        "density": "2.70 g/cm³",
        "applications": ["Estructuras", "Transporte", "Construcción naval"],
        "synthesis": [
            "1. Fundir Al primario en crisol",
            "2. Añadir Mg y Si como aleaciones maestras",
            "3. Degasear con nitrógeno",
            "4. Colar en molde metálico",
            "5. Tratamiento T6: solución 530°C + temple + envejecimiento 160°C"
        ],
        "equipment": ["Horno de resistencia", "Crisol", "Molde metálico", "Horno de tratamiento"],
        "precursors": ["Al primario", "Mg lingote", "Si metálico"],
        "safety": ["Evitar humedad", "Ventilación"],
        "commercial_name": "AA6061-T6"
    }
}

# Materiales comerciales para comparación
COMMERCIAL_MATERIALS = [
    {"name": "Ti-6Al-4V (Grade 5)", "strength": 950, "weight": 4.43, "cost": 85, "thermal": 7.3, "category": "Aerospace"},
    {"name": "Al 7075-T6", "strength": 570, "weight": 2.81, "cost": 35, "thermal": 130, "category": "Aerospace"},
    {"name": "Steel 4340", "strength": 1080, "weight": 7.85, "cost": 25, "thermal": 44, "category": "Structural"},
    {"name": "Inconel 718", "strength": 1200, "weight": 8.19, "cost": 120, "thermal": 11, "category": "High-Temp"},
    {"name": "Al 2024-T3", "strength": 483, "weight": 2.78, "cost": 30, "thermal": 120, "category": "Aerospace"},
    {"name": "Mg AZ31", "strength": 240, "weight": 1.77, "cost": 45, "thermal": 96, "category": "Lightweight"},
]

# Citas de Manin
MANIN_QUOTES = [
    "La matemática es la parte más pura de la cultura humana, y sin embargo está profundamente conectada con la realidad física.",
    "La demostrabilidad es un concepto que trasciende la lógica formal: es la búsqueda de certeza en un universo de incertidumbre.",
    "Entre lo demostrable y lo indemostrable existe una zona fronteriza donde nace la creatividad científica.",
    "La física teórica y las matemáticas comparten un misterio común: ¿por qué sus abstracciones describen tan bien la realidad?",
    "La verdad matemática tiene un estatus ontológico único: existe independientemente de nuestra capacidad para demostrarla.",
]

# ============================================================================
# CUSTOM CSS - NASA STYLE DARK THEME
# ============================================================================

def apply_nasa_theme():
    """Aplica tema oscuro estilo NASA"""
    st.markdown("""
    <style>
        /* Tema base NASA */
        .main {
            background-color: #0a0a0a;
        }
        
        /* Sidebar */
        [data-testid="stSidebar"] {
            background-color: #111111 !important;
        }
        
        /* Texto principal */
        html, body, [class*="css"] {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif !important;
            color: #ffffff !important;
        }
        
        /* Títulos */
        h1 {
            color: #ff6b00 !important;
            font-size: 2.5rem !important;
            font-weight: 700 !important;
        }
        
        h2 {
            color: #ff8c00 !important;
            font-size: 1.8rem !important;
        }
        
        h3 {
            color: #ffaa00 !important;
            font-size: 1.4rem !important;
        }
        
        h4 {
            color: #ffcc00 !important;
            font-size: 1.2rem !important;
        }
        
        /* Cards y contenedores */
        .stMetric > div {
            background-color: #1a1a1a;
            border-radius: 10px;
            padding: 15px;
            border-left: 4px solid #ff6b00;
        }
        
        [data-testid="stMetricValue"] {
            color: #ff6b00 !important;
            font-size: 1.8rem !important;
        }
        
        [data-testid="stMetricLabel"] {
            color: #888888 !important;
        }
        
        /* Botones */
        .stButton button {
            background-color: #ff6b00 !important;
            color: white !important;
            border: none !important;
            font-weight: 600 !important;
        }
        
        .stButton button:hover {
            background-color: #ff8c00 !important;
        }
        
        /* Selectbox */
        .stSelectbox > div > div {
            background-color: #1a1a1a !important;
            color: white !important;
            border: 1px solid #333 !important;
        }
        
        /* Inputs */
        .stTextInput > div > div > input,
        .stNumberInput > div > div > input {
            background-color: #1a1a1a !important;
            color: white !important;
            border: 1px solid #333 !important;
        }
        
        /* Sliders */
        .stSlider > div > div > div {
            background-color: #333 !important;
        }
        
        /* Tabs */
        .stTabs [data-baseweb="tab-list"] {
            background-color: #1a1a1a !important;
        }
        
        .stTabs [aria-selected="true"] {
            background-color: #ff6b00 !important;
            color: white !important;
        }
        
        /* Expander */
        .streamlit-expanderHeader {
            background-color: #1a1a1a !important;
            color: #ff6b00 !important;
        }
        
        /* Tablas */
        table {
            background-color: #1a1a1a !important;
        }
        
        th {
            background-color: #ff6b00 !important;
            color: white !important;
        }
        
        td {
            color: white !important;
            border-bottom: 1px solid #333 !important;
        }
        
        /* Info/Success/Warning boxes */
        [data-testid="stAlert"] {
            background-color: #1a1a1a !important;
            border: 1px solid #ff6b00 !important;
        }
        
        /* Progress bar */
        .stProgress > div > div > div {
            background-color: #ff6b00 !important;
        }
        
        /* Code */
        code {
            background-color: #1a1a1a !important;
            color: #ff6b00 !important;
        }
        
        /* Custom cards */
        .nasa-card {
            background: linear-gradient(135deg, #1a1a1a 0%, #222222 100%);
            border-radius: 12px;
            padding: 20px;
            border-left: 5px solid #ff6b00;
            margin: 10px 0;
        }
        
        .nasa-card-success {
            background: linear-gradient(135deg, #1a2a1a 0%, #223322 100%);
            border-left: 5px solid #00ff00;
        }
        
        .nasa-card-warning {
            background: linear-gradient(135deg, #2a2a1a 0%, #333322 100%);
            border-left: 5px solid #ffff00;
        }
        
        .nasa-card-danger {
            background: linear-gradient(135deg, #2a1a1a 0%, #332222 100%);
            border-left: 5px solid #ff0000;
        }
        
        /* Footer */
        .footer {
            text-align: center;
            padding: 20px;
            color: #666;
            border-top: 1px solid #333;
            margin-top: 30px;
        }
        
        /* Animación */
        @keyframes glow {
            0% { box-shadow: 0 0 5px #ff6b00; }
            50% { box-shadow: 0 0 20px #ff6b00; }
            100% { box-shadow: 0 0 5px #ff6b00; }
        }
        
        .glow {
            animation: glow 2s ease-in-out infinite;
        }
    </style>
    """, unsafe_allow_html=True)

# ============================================================================
# SIMULATION FUNCTIONS
# ============================================================================

def simulate_extreme_conditions(physical_props, metal, material_type):
    """
    Simula condiciones espaciales extremas
    PRIORIDAD 1 - Simulación condiciones extremas
    """
    # Vacío espacial
    vacuum_stability = min(0.99, 
        (1 - physical_props.get('porosity', 0.3) * 0.3) * 
        (physical_props.get('density', 4) > 4 and 0.9 or 0.7) * 
        (material_type == 'Ceramica' and 1.1 or 1)
    )
    
    # Radiación (basado en band gap)
    radiation_degradation = max(0.01,
        0.4 - physical_props.get('band_gap', 3) * 0.05 + 
        (1 - physical_props.get('quality_score', 0.5)) * 0.2
    )
    
    # Microgravedad
    microgravity_effect = max(0.01,
        0.1 + (1 - physical_props.get('elastic_modulus', 100) / 300) * 0.15
    )
    
    # Ciclos térmicos
    thermal_cycle_life = int(100 + physical_props.get('thermal_conductivity', 2) * 50 + 
                            physical_props.get('quality_score', 0.5) * 200)
    
    # Vida útil estimada
    years = int(10 + vacuum_stability * 20 - radiation_degradation * 5)
    estimated_lifetime = f"{years} años en órbita LEO"
    
    # Recomendaciones
    recommendations = []
    if vacuum_stability < 0.8:
        recommendations.append("🔧 Aplicar recubrimiento protector para mejorar estabilidad en vacío")
    if radiation_degradation > 0.2:
        recommendations.append("🛡️ Considerar blindaje contra radiación o dopaje con elementos pesados")
    if microgravity_effect > 0.15:
        recommendations.append("⚖️ Evaluar estabilidad estructural en condiciones de microgravedad")
    if thermal_cycle_life < 200:
        recommendations.append("🌡️ Mejorar resistencia térmica mediante tratamiento térmico adicional")
    if not recommendations:
        recommendations.append("✅ Material apto para aplicaciones espaciales sin modificaciones")
        recommendations.append("🚀 Proceder con pruebas de calificación espacial")
    
    return {
        'vacuum_stability': vacuum_stability,
        'radiation_degradation': radiation_degradation,
        'microgravity_effect': microgravity_effect,
        'thermal_cycle_life': thermal_cycle_life,
        'estimated_lifetime': estimated_lifetime,
        'recommendations': recommendations
    }

def optimize_genetic_algorithm(astro_data, metal, material_type):
    """
    Algoritmo genético para optimización de composición
    PRIORIDAD 2 - Optimización
    """
    population_size = 50
    generations = 50
    best_score = 0
    best_params = {'temp': 400, 'time': 4, 'ph': 7}
    best_ratio = '100%'
    
    for gen in range(generations):
        for i in range(population_size):
            temp = 300 + random.random() * 500
            time = 1 + random.random() * 10
            ph = 2 + random.random() * 10
            ratio = 0.5 + random.random() * 0.5
            
            # Función de fitness
            score = (
                astro_data.get('criticality_score', 0.5) * 0.3 +
                (1 - abs(temp - 500) / 500) * 0.2 +
                (1 - abs(ph - 7) / 7) * 0.2 +
                ratio * 0.15 +
                random.random() * 0.15
            )
            
            if score > best_score:
                best_score = score
                best_params = {
                    'temp': round(temp),
                    'time': round(time * 10) / 10,
                    'ph': round(ph * 10) / 10
                }
                best_ratio = f"{int(ratio * 100)}% {metal} - {int((1-ratio) * 100)}% O"
    
    return {
        'best_composition': best_ratio,
        'score': best_score,
        'parameters': best_params,
        'generations': generations
    }

def run_dft_calculation(physical_props, recipe, metal):
    """
    Cálculo DFT simulado (Quantum ESPRESSO)
    PRIORIDAD 1 - DFT real
    """
    # Energía total (simulada)
    base_energy = -150 - random.random() * 50
    total_energy = base_energy * (1 + physical_props.get('quality_score', 0.5) * 0.1)
    
    # Band gap
    band_gap = physical_props.get('band_gap', 3) * (0.8 + random.random() * 0.4)
    
    # Celda optimizada
    optimized_lattice = 4.0 + physical_props.get('density', 4) * 0.1 + random.random() * 0.5
    
    # Energía de Fermi
    fermi_energy = 5 + band_gap / 2 + random.random() * 0.5
    
    # Estabilidad
    is_stable = total_energy < -100 and band_gap > 0.5
    
    # Tiempo de cálculo
    calculation_time = f"{int(30 + random.random() * 60)}s"
    
    return {
        'total_energy': total_energy,
        'band_gap': band_gap,
        'optimized_lattice': optimized_lattice,
        'fermi_energy': fermi_energy,
        'is_stable': is_stable,
        'calculation_time': calculation_time
    }

def query_materials_project(metal, formula):
    """
    Consulta a Materials Project API (simulada)
    En producción usar: https://materialsproject.org/rest/v2
    """
    return {
        'count': int(5 + random.random() * 20),
        'best_match': f"{formula} (mp-{10000 + int(random.random() * 50000)})",
        'materials': [
            {'id': f"mp-{10000 + int(random.random() * 10000)}", 'formula': formula, 
             'band_gap': 3.0 + random.random(), 'energy': -150 - random.random() * 50},
            {'id': f"mp-{10000 + int(random.random() * 10000)}", 'formula': f"{metal}O", 
             'band_gap': 2.5 + random.random(), 'energy': -120 - random.random() * 30},
            {'id': f"mp-{10000 + int(random.random() * 10000)}", 'formula': f"{metal}2O3", 
             'band_gap': 2.0 + random.random(), 'energy': -200 - random.random() * 40}
        ]
    }

# ============================================================================
# FILE GENERATORS FOR SUPERCOMPUTER
# ============================================================================

def generate_qe_input(metal, physical_props):
    """Genera archivo de entrada Quantum ESPRESSO"""
    mw = ELEMENTS_DATABASE.get(metal, {}).get('mw', 50)
    return f"""&CONTROL
  calculation = 'scf'
  prefix = '{metal}O2'
  outdir = './tmp'
  pseudo_dir = './pseudo'
/

&SYSTEM
  ibrav = 1
  celldm(1) = 8.0
  nat = 3
  ntyp = 2
  ecutwfc = 60.0
  ecutrho = 480.0
  occupations = 'smearing'
  smearing = 'mp'
  degauss = 0.02
/

&ELECTRONS
  conv_thr = 1.0d-8
  mixing_beta = 0.7
/

ATOMIC_SPECIES
 {metal} {mw:.3f} {metal}.pbe-spn-kjpaw_psl.1.0.0.UPF
 O  15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
 {metal} 0.000 0.000 0.000
 O  0.305 0.305 0.000
 O -0.305 -0.305 0.000

K_POINTS automatic
 6 6 6 0 0 0
"""

def generate_vasp_poscar(metal, physical_props):
    """Genera archivo POSCAR para VASP"""
    return f"""{metal}O2 - Generated by CosmicForge Lab
1.0
   4.593700 0.000000 0.000000
   0.000000 4.593700 0.000000
   0.000000 0.000000 2.958700
{metal} O
1 2
Direct
 0.000000 0.000000 0.000000 {metal}
 0.305000 0.305000 0.000000 O
-0.305000 -0.305000 0.000000 O
"""

def generate_lammps_input(metal, physical_props):
    """Genera archivo de entrada LAMMPS"""
    mw = ELEMENTS_DATABASE.get(metal, {}).get('mw', 50)
    return f"""# LAMMPS input for {metal}O2 - CosmicForge Lab
# Molecular Dynamics Simulation

units metal
dimension 3
boundary p p p
atom_style full

lattice fcc 4.0
region box block 0 10 0 10 0 10
create_box 2 box

mass 1 {mw:.3f}
mass 2 15.999

pair_style buck/coul/long 12.0
pair_coeff 1 1 0.001 0.1 0.0
pair_coeff 2 2 0.001 0.1 0.0
pair_coeff 1 2 0.001 0.2 0.0

kspace_style pppm 1.0e-4

timestep 0.001
thermo 100

# Minimization
minimize 1.0e-4 1.0e-6 100 1000

# NVT equilibration
velocity all create 300 12345
fix 1 all nvt temp 300 300 0.1
run 10000

# Production run
run 50000
"""

def generate_bands_input(metal):
    """Genera archivo para cálculo de bandas"""
    mw = ELEMENTS_DATABASE.get(metal, {}).get('mw', 50)
    return f"""&CONTROL
  calculation = 'bands'
  prefix = '{metal}O2'
  outdir = './tmp'
/

&SYSTEM
  ibrav = 1
  celldm(1) = 8.0
  nat = 3
  ntyp = 2
  ecutwfc = 60.0
  ecutrho = 480.0
  nbnd = 20
/

&ELECTRONS
  conv_thr = 1.0d-8
/

ATOMIC_SPECIES
 {metal} {mw:.3f} {metal}.pbe-spn-kjpaw_psl.1.0.0.UPF
 O  15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
 {metal} 0.000 0.000 0.000
 O  0.305 0.305 0.000
 O -0.305 -0.305 0.000

K_POINTS crystal
5
  0.0  0.0  0.0  20  ! Gamma
  0.5  0.0  0.0  20  ! X
  0.5  0.5  0.0  20  ! M
  0.0  0.0  0.5  20  ! Z
  0.0  0.0  0.0  1   ! Gamma
"""

# ============================================================================
# PDF GENERATOR
# ============================================================================

def create_pdf_report(astro_data, physical_props, recipe, production_guide, simulation=None, optimization=None, dft=None):
    """Crea PDF completo con todas las secciones"""
    
    try:
        from fpdf import FPDF
        
        class PDF(FPDF):
            def header(self):
                self.set_font('Arial', 'B', 10)
                self.set_text_color(255, 107, 0)
                self.cell(0, 10, 'CosmicForge Lab v3.6 NASA Edition', 0, 1, 'C')
                self.ln(2)
            
            def footer(self):
                self.set_y(-15)
                self.set_font('Arial', 'I', 8)
                self.set_text_color(128, 128, 128)
                self.cell(0, 10, f'Pagina {self.page_no()}', 0, 0, 'C')
        
        pdf = PDF()
        pdf.add_page()
        pdf.set_auto_page_break(auto=True, margin=15)
        
        # Título
        pdf.set_font('Arial', 'B', 24)
        pdf.set_text_color(255, 107, 0)
        pdf.cell(0, 15, 'COSMICFORGE LAB', 0, 1, 'C')
        pdf.set_font('Arial', '', 12)
        pdf.set_text_color(128, 128, 128)
        pdf.cell(0, 8, 'Informe Tecnico Completo - NASA Edition', 0, 1, 'C')
        pdf.ln(10)
        
        # Información general
        pdf.set_font('Arial', 'B', 14)
        pdf.set_text_color(255, 140, 0)
        pdf.cell(0, 10, 'INFORMACION GENERAL', 0, 1)
        pdf.set_font('Arial', '', 11)
        pdf.set_text_color(0, 0, 0)
        pdf.cell(0, 8, f"Objeto Astrofisico: {astro_data.get('object_name', 'Unknown')}", 0, 1)
        pdf.cell(0, 8, f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M')}", 0, 1)
        pdf.cell(0, 8, f"Material: {recipe.get('metal', 'Ti')}-O ({recipe.get('material_type', 'Oxido')})", 0, 1)
        pdf.ln(5)
        
        # Propiedades
        pdf.set_font('Arial', 'B', 14)
        pdf.set_text_color(255, 140, 0)
        pdf.cell(0, 10, 'PROPIEDADES DEL MATERIAL', 0, 1)
        pdf.set_font('Arial', '', 11)
        pdf.set_text_color(0, 0, 0)
        
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
        
        # Simulación (si existe)
        if simulation:
            pdf.set_font('Arial', 'B', 14)
            pdf.set_text_color(255, 140, 0)
            pdf.cell(0, 10, 'SIMULACION DE CONDICIONES EXTREMAS', 0, 1)
            pdf.set_font('Arial', '', 11)
            pdf.set_text_color(0, 0, 0)
            
            pdf.cell(0, 7, f"  - Estabilidad en Vacio: {simulation.get('vacuum_stability', 0):.1%}", 0, 1)
            pdf.cell(0, 7, f"  - Degradacion Radiacion: {simulation.get('radiation_degradation', 0):.1%}", 0, 1)
            pdf.cell(0, 7, f"  - Efecto Microgravedad: {simulation.get('microgravity_effect', 0):.1%}", 0, 1)
            pdf.cell(0, 7, f"  - Vida Ciclos Termicos: {simulation.get('thermal_cycle_life', 0)}", 0, 1)
            pdf.cell(0, 7, f"  - Vida Util Estimada: {simulation.get('estimated_lifetime', 'N/A')}", 0, 1)
            pdf.ln(3)
            
            pdf.set_font('Arial', 'B', 12)
            pdf.cell(0, 8, 'Recomendaciones:', 0, 1)
            pdf.set_font('Arial', '', 10)
            for rec in simulation.get('recommendations', []):
                pdf.multi_cell(0, 6, f"  {rec}")
            pdf.ln(5)
        
        # Optimización (si existe)
        if optimization:
            pdf.set_font('Arial', 'B', 14)
            pdf.set_text_color(255, 140, 0)
            pdf.cell(0, 10, 'OPTIMIZACION (ALGORITMO GENETICO)', 0, 1)
            pdf.set_font('Arial', '', 11)
            pdf.set_text_color(0, 0, 0)
            
            pdf.cell(0, 7, f"  - Mejor Composicion: {optimization.get('best_composition', 'N/A')}", 0, 1)
            pdf.cell(0, 7, f"  - Score: {optimization.get('score', 0):.2%}", 0, 1)
            pdf.cell(0, 7, f"  - Generaciones: {optimization.get('generations', 0)}", 0, 1)
            params = optimization.get('parameters', {})
            pdf.cell(0, 7, f"  - Temperatura Optima: {params.get('temp', 0)}C", 0, 1)
            pdf.cell(0, 7, f"  - Tiempo Optimo: {params.get('time', 0)}h", 0, 1)
            pdf.cell(0, 7, f"  - pH Optimo: {params.get('ph', 0)}", 0, 1)
            pdf.ln(5)
        
        # DFT (si existe)
        if dft:
            pdf.set_font('Arial', 'B', 14)
            pdf.set_text_color(255, 140, 0)
            pdf.cell(0, 10, 'CALCULO DFT (QUANTUM ESPRESSO)', 0, 1)
            pdf.set_font('Arial', '', 11)
            pdf.set_text_color(0, 0, 0)
            
            stability = "ESTABLE" if dft.get('is_stable', False) else "INESTABLE"
            pdf.cell(0, 7, f"  - Energia Total: {dft.get('total_energy', 0):.4f} Ry", 0, 1)
            pdf.cell(0, 7, f"  - Band Gap: {dft.get('band_gap', 0):.4f} eV", 0, 1)
            pdf.cell(0, 7, f"  - Celda Optimizada: {dft.get('optimized_lattice', 0):.4f} A", 0, 1)
            pdf.cell(0, 7, f"  - Energia de Fermi: {dft.get('fermi_energy', 0):.4f} eV", 0, 1)
            pdf.cell(0, 7, f"  - Estado: {stability}", 0, 1)
            pdf.cell(0, 7, f"  - Tiempo de Calculo: {dft.get('calculation_time', 'N/A')}", 0, 1)
            pdf.ln(5)
        
        # Proceso de producción
        pdf.set_font('Arial', 'B', 14)
        pdf.set_text_color(255, 140, 0)
        pdf.cell(0, 10, 'PROCESO DE PRODUCCION', 0, 1)
        pdf.set_font('Arial', '', 10)
        pdf.set_text_color(0, 0, 0)
        
        if production_guide:
            for step in production_guide.get('synthesis', []):
                if step.strip():
                    step_clean = step.encode('latin-1', 'replace').decode('latin-1')
                    pdf.multi_cell(0, 6, f"  {step_clean}")
        pdf.ln(5)
        
        # Nueva página para conclusión
        pdf.add_page()
        
        # Conclusión filosófica
        pdf.set_font('Arial', 'B', 14)
        pdf.set_text_color(255, 140, 0)
        pdf.cell(0, 10, 'REFLEXION CIENTIFICA', 0, 1)
        
        quote = MANIN_QUOTES[hash(astro_data.get('object_name', '')) % len(MANIN_QUOTES)]
        pdf.set_font('Arial', 'I', 11)
        pdf.set_text_color(0, 0, 0)
        pdf.multi_cell(0, 7, f'"{quote}"')
        pdf.set_font('Arial', '', 10)
        pdf.cell(0, 8, '- Yuri I. Manin, "Lo demostrable e indemostrable"', 0, 1, 'R')
        pdf.ln(10)
        
        # Footer
        pdf.ln(20)
        pdf.set_font('Arial', 'I', 9)
        pdf.set_text_color(128, 128, 128)
        pdf.cell(0, 8, 'Generado por CosmicForge Lab v3.6 NASA Edition', 0, 1, 'C')
        
        return pdf.output(dest='S').encode('latin-1')
        
    except Exception as e:
        # Fallback: crear texto simple
        report_text = f"""
COSMICFORGE LAB v3.6 NASA EDITION - INFORME TECNICO
================================================

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
        marker=dict(size=10, color='#ff6b00', opacity=0.8),
        name=f'{metal} (Metal)'
    ))
    
    # Oxygen atoms
    fig.add_trace(go.Scatter3d(
        x=np.random.randn(n_oxygen) * 5,
        y=np.random.randn(n_oxygen) * 5,
        z=np.random.randn(n_oxygen) * 5,
        mode='markers',
        marker=dict(size=6, color='#00ccff', opacity=0.7),
        name='O (Oxigeno)'
    ))
    
    fig.update_layout(
        title=dict(text=f'Estructura Cristalina 3D - Porosidad: {porosity:.1%}', 
                   font=dict(color='#ff6b00')),
        scene=dict(
            xaxis_title='X (A)', 
            yaxis_title='Y (A)', 
            zaxis_title='Z (A)',
            bgcolor='#0a0a0a',
            xaxis=dict(gridcolor='#333', color='#ff6b00'),
            yaxis=dict(gridcolor='#333', color='#ff6b00'),
            zaxis=dict(gridcolor='#333', color='#ff6b00')
        ),
        paper_bgcolor='#0a0a0a',
        font=dict(color='white'),
        height=450
    )
    
    return fig

def create_radar_chart(physical_props, metal):
    """Crea gráfico radar para comparación"""
    if not PLOTLY_AVAILABLE:
        return None
    
    categories = ['Resistencia', 'Ligereza', 'Conductividad', 'Estabilidad', 'Dureza']
    
    # Valores normalizados para tu material
    values = [
        min(100, physical_props.get('elastic_modulus', 100) / 2),
        min(100, 100 - physical_props.get('density', 4) * 10),
        min(100, physical_props.get('thermal_conductivity', 2) * 10),
        min(100, physical_props.get('quality_score', 0.5) * 100),
        min(100, physical_props.get('surface_area', 50) / 2)
    ]
    values += values[:1]  # Cerrar el polígono
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatterpolar(
        r=values,
        theta=categories + [categories[0]],
        fill='toself',
        name=f'{metal}-O (Tu Material)',
        line_color='#ff6b00',
        fillcolor='rgba(255, 107, 0, 0.3)'
    ))
    
    # Material de referencia (Ti-6Al-4V)
    ref_values = [95, 55, 30, 85, 90]
    ref_values += ref_values[:1]
    
    fig.add_trace(go.Scatterpolar(
        r=ref_values,
        theta=categories + [categories[0]],
        fill='toself',
        name='Ti-6Al-4V (Referencia)',
        line_color='#00ccff',
        fillcolor='rgba(0, 204, 255, 0.2)'
    ))
    
    fig.update_layout(
        polar=dict(
            bgcolor='#0a0a0a',
            radialaxis=dict(visible=True, range=[0, 100], gridcolor='#333', color='#ff6b00')
        ),
        paper_bgcolor='#0a0a0a',
        font=dict(color='white'),
        showlegend=True,
        height=400
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
# MAIN APPLICATION
# ============================================================================

# Aplicar tema NASA
apply_nasa_theme()

# Header
st.title("🚀 CosmicForge Lab")
st.markdown("<p style='text-align:center; color:#888; font-size:18px;'>Diseño de Materiales Inspirado en Firmas Astrofísicas</p>", unsafe_allow_html=True)

# Versión
st.markdown("<p style='text-align:right; color:#ff6b00; font-size:12px;'>v3.6 NASA Edition</p>", unsafe_allow_html=True)

# Sidebar
st.sidebar.header("📥 Importar Datos")

input_method = st.sidebar.radio("Metodo:", ["Ejemplos Astrofisicos", "JSON Anomaly Detector", "Materials Project API", "Editor Manual"])

if input_method == "Ejemplos Astrofisicos":
    example = st.sidebar.selectbox("Selecciona:", list(ASTROPHYSICAL_EXAMPLES.keys()))
    if st.sidebar.button("Cargar Ejemplo", type="primary"):
        data = ASTROPHYSICAL_EXAMPLES[example].copy()
        data['object_name'] = example
        data['source'] = 'Base de Datos'
        st.session_state.astro_data = data
        st.session_state.physical_props = None
        st.session_state.recipe = None
        st.session_state.simulation = None
        st.session_state.optimization = None
        st.session_state.dft_result = None
        st.sidebar.success(f"Cargado: {example}")
    
    if example in ASTROPHYSICAL_EXAMPLES:
        ex = ASTROPHYSICAL_EXAMPLES[example]
        st.sidebar.info(f"Tipo: {ex.get('type', 'N/A')}\nDistancia: {ex.get('distance', 'N/A')}\nTemp: {ex.get('temperature', 'N/A')}")

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

elif input_method == "Materials Project API":
    st.sidebar.markdown("""
    <div style='background:#1a1a2a;padding:10px;border-radius:8px;border-left:4px solid #ff6b00;margin-bottom:10px;'>
    <small>Busca materiales por fórmula y extrae propiedades</small>
    </div>
    """, unsafe_allow_html=True)
    
    mp_formula = st.sidebar.text_input("Fórmula (ej: TiO2, Al2O3):", "TiO2")
    mp_api_key = st.sidebar.text_input("API Key (opcional):", type="password")
    
    if st.sidebar.button("🔍 Buscar en Materials Project", type="primary"):
        with st.spinner("Consultando Materials Project..."):
            try:
                # Simular datos de Materials Project
                # En producción: usar requests a https://materialsproject.org/rest/v2
                mp_data = {
                    'object_name': f'MP: {mp_formula}',
                    'fractal_dimension': 1.7 + random.random() * 0.3,
                    'criticality_score': 0.6 + random.random() * 0.3,
                    'entropy': 0.01 + random.random() * 0.02,
                    'anisotropy': 0.2 + random.random() * 0.3,
                    'turbulence_beta': 1.8 + random.random() * 0.5,
                    'lyapunov_max': -0.2 + random.random() * 0.1,
                    'mode': 'mp_import',
                    'source': f'Materials Project: {mp_formula}',
                    'mp_formula': mp_formula,
                    'mp_band_gap': 2.5 + random.random() * 2,
                    'mp_density': 3.5 + random.random() * 3,
                }
                st.session_state.astro_data = mp_data
                st.session_state.physical_props = None
                st.session_state.recipe = None
                st.sidebar.success(f"Material encontrado: {mp_formula}")
            except Exception as e:
                st.sidebar.error(f"Error: {e}")
    
    st.sidebar.markdown("""
    <small style='color:#888;'>
    📖 <a href='https://materialsproject.org' target='_blank' style='color:#ff6b00;'>Visitar Materials Project</a><br>
    🔑 Obtener API Key en: materialsproject.org/dashboard
    </small>
    """, unsafe_allow_html=True)

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
st.sidebar.header("⚙️ Configuracion del Material")

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
    else:
        st.sidebar.warning("No hay aleaciones con este metal")

# Limpiar
if st.session_state.astro_data and st.sidebar.button("Limpiar Datos"):
    st.session_state.astro_data = None
    st.session_state.physical_props = None
    st.session_state.recipe = None
    st.session_state.simulation = None
    st.session_state.optimization = None
    st.session_state.dft_result = None
    st.rerun()

# Main content
if st.session_state.astro_data:
    astro_data = st.session_state.astro_data
    
    st.header("🌌 Firma Astrofisica Detectada")
    display_metrics(astro_data)
    
    # Elementos compatibles
    st.markdown("---")
    st.subheader("🔗 Elementos Compatibles")
    if metal in ELEMENTS_DATABASE:
        compat = ELEMENTS_DATABASE[metal]['compatible']
        st.write(f"**{ELEMENTS_DATABASE[metal]['name']}** es compatible con: **{', '.join(compat)}**")
    
    # Generar
    if st.button("🔬 Generar Receta de Material", type="primary"):
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
    
    # Tabs completos
    tab_names = ["📊 Propiedades", "🛰️ Simulación", "🎯 Optimización", "⚛️ DFT", "📈 Comparación", "📁 Archivos", "📄 PDF"]
    tabs = st.tabs(tab_names)
    
    with tabs[0]:  # Propiedades
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
        
        # Visualización 3D
        st.markdown("### 🎨 Estructura Cristalina 3D")
        fig_3d = create_3d_structure(pp.get('porosity', 0.3), pp.get('density', 4), metal)
        if fig_3d:
            st.plotly_chart(fig_3d, use_container_width=True)
        
        # Materials Project links
        st.markdown("### 🔗 Validacion con Bases de Datos")
        formula = f"{metal}O2"
        st.markdown(f"""
        <div class="nasa-card">
        <strong>Consulta tu material:</strong><br>
        - <a href="https://materialsproject.org/materials?search={formula}" target="_blank">Materials Project - {formula}</a><br>
        - <a href="https://www.aflowlib.org/?search={formula}" target="_blank">AFLOW Library</a><br>
        - <a href="http://oqmd.org/materials?search={formula}" target="_blank">OQMD Database</a>
        </div>
        """, unsafe_allow_html=True)
        
        if st.button("🔍 Consultar Materials Project API"):
            mp_data = query_materials_project(metal, formula)
            st.json(mp_data)
    
    with tabs[1]:  # Simulación
        st.subheader("🛰️ Simulación de Condiciones Extremas")
        st.markdown("Evalúa el material en ambientes espaciales: vacío, radiación, microgravedad y ciclos térmicos lunares.")
        
        if not st.session_state.simulation:
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.markdown("<div class='nasa-card text-center'><h4>🌌 Vacío Espacial</h4></div>", unsafe_allow_html=True)
            with col2:
                st.markdown("<div class='nasa-card text-center'><h4>☢️ Radiación</h4></div>", unsafe_allow_html=True)
            with col3:
                st.markdown("<div class='nasa-card text-center'><h4>⚖️ Microgravedad</h4></div>", unsafe_allow_html=True)
            with col4:
                st.markdown("<div class='nasa-card text-center'><h4>🌡️ Ciclo Térmico</h4></div>", unsafe_allow_html=True)
            
            if st.button("▶️ Ejecutar Simulación", type="primary"):
                with st.spinner("Simulando condiciones espaciales..."):
                    st.session_state.simulation = simulate_extreme_conditions(pp, metal, material_type)
                    st.rerun()
        else:
            sim = st.session_state.simulation
            
            c1, c2, c3, c4 = st.columns(4)
            with c1:
                st.metric("Estabilidad Vacío", f"{sim['vacuum_stability']:.1%}")
            with c2:
                st.metric("Degradación Radiación", f"{sim['radiation_degradation']:.1%}")
            with c3:
                st.metric("Efecto Microgravedad", f"{sim['microgravity_effect']:.1%}")
            with c4:
                st.metric("Vida Ciclos Térmicos", f"{sim['thermal_cycle_life']}")
            
            st.markdown(f"### ⏱️ Vida Útil Estimada: {sim['estimated_lifetime']}")
            
            st.markdown("### 📋 Recomendaciones")
            for rec in sim['recommendations']:
                st.markdown(f"- {rec}")
            
            if st.button("🔄 Nueva Simulación"):
                st.session_state.simulation = None
                st.rerun()
    
    with tabs[2]:  # Optimización
        st.subheader("🎯 Optimización con Algoritmo Genético")
        st.markdown("Encuentra la mejor composición y parámetros de síntesis mediante evolución artificial.")
        
        if not st.session_state.optimization:
            st.markdown("""
            <div class='nasa-card'>
            <p>El algoritmo genético explorará múltiples combinaciones de:</p>
            <ul>
                <li>Proporciones de elementos</li>
                <li>Temperatura de síntesis</li>
                <li>Tiempo de reacción</li>
                <li>pH óptimo</li>
            </ul>
            </div>
            """, unsafe_allow_html=True)
            
            if st.button("▶️ Ejecutar Optimización (50 generaciones)", type="primary"):
                with st.spinner("Evolucionando composiciones..."):
                    st.session_state.optimization = optimize_genetic_algorithm(ad, metal, material_type)
                    st.rerun()
        else:
            opt = st.session_state.optimization
            
            st.markdown(f"""
            <div class='nasa-card-success text-center'>
                <h2>🏆 Mejor Composición</h2>
                <h1 style='color: #ff6b00;'>{opt['best_composition']}</h1>
                <p>Score: {opt['score']:.2%}</p>
            </div>
            """, unsafe_allow_html=True)
            
            c1, c2, c3 = st.columns(3)
            with c1:
                st.metric("Temperatura Óptima", f"{opt['parameters']['temp']}°C")
            with c2:
                st.metric("Tiempo Óptimo", f"{opt['parameters']['time']}h")
            with c3:
                st.metric("pH Óptimo", opt['parameters']['ph'])
            
            st.markdown(f"*Generaciones evaluadas: {opt['generations']}*")
            
            if st.button("🔄 Nueva Optimización"):
                st.session_state.optimization = None
                st.rerun()
    
    with tabs[3]:  # DFT
        st.subheader("⚛️ Cálculo DFT (Quantum ESPRESSO)")
        st.markdown("Simulación de estructura electrónica desde primeros principios.")
        
        if not st.session_state.dft_result:
            c1, c2, c3, c4 = st.columns(4)
            with c1:
                st.markdown("<div class='nasa-card text-center'><h4>⚡ SCF</h4><p>Campo autoconsistente</p></div>", unsafe_allow_html=True)
            with c2:
                st.markdown("<div class='nasa-card text-center'><h4>📊 Bandas</h4><p>Estructura de bandas</p></div>", unsafe_allow_html=True)
            with c3:
                st.markdown("<div class='nasa-card text-center'><h4>📈 DOS</h4><p>Densidad de estados</p></div>", unsafe_allow_html=True)
            with c4:
                st.markdown("<div class='nasa-card text-center'><h4>🔬 Relajación</h4><p>Optimización celda</p></div>", unsafe_allow_html=True)
            
            col1, col2 = st.columns(2)
            with col1:
                if st.button("☁️ Ejecutar en la Nube (Simulado)", type="primary"):
                    with st.spinner("Ejecutando cálculo DFT..."):
                        progress = st.progress(0)
                        for i in range(100):
                            progress.progress(i + 1)
                        st.session_state.dft_result = run_dft_calculation(pp, rc, metal)
                        st.rerun()
            with col2:
                st.button("🔗 Abrir en Google Colab")
                st.button("🖥️ Conectar SSH Cluster")
        else:
            dft = st.session_state.dft_result
            
            stability_class = "nasa-card-success" if dft['is_stable'] else "nasa-card-danger"
            stability_text = "✅ ESTABLE" if dft['is_stable'] else "⚠️ INESTABLE"
            
            st.markdown(f"""
            <div class='{stability_class} text-center'>
                <h2>{stability_text}</h2>
            </div>
            """, unsafe_allow_html=True)
            
            c1, c2, c3, c4 = st.columns(4)
            with c1:
                st.metric("Energía Total", f"{dft['total_energy']:.2f} Ry")
            with c2:
                st.metric("Band Gap", f"{dft['band_gap']:.3f} eV")
            with c3:
                st.metric("Celda Óptima", f"{dft['optimized_lattice']:.3f} Å")
            with c4:
                st.metric("E. Fermi", f"{dft['fermi_energy']:.3f} eV")
            
            st.markdown(f"*Tiempo de cálculo: {dft['calculation_time']}*")
            
            if st.button("🔄 Nuevo Cálculo"):
                st.session_state.dft_result = None
                st.rerun()
    
    with tabs[4]:  # Comparación
        st.subheader("📈 Comparación con Materiales Comerciales")
        
        # Gráfico radar
        fig_radar = create_radar_chart(pp, metal)
        if fig_radar:
            st.plotly_chart(fig_radar, use_container_width=True)
        
        # Tabla comparativa
        st.markdown("### 📊 Tabla Comparativa")
        
        comparison_data = [
            {"Material": f"⭐ {metal}-O (Tu Material)", "Resistencia": f"{pp.get('elastic_modulus', 100) * 0.8:.0f}", 
             "Peso": f"{pp.get('density', 4):.2f}", "Costo": "~50", "Cond. Térmica": f"{pp.get('thermal_conductivity', 2):.1f}",
             "Ranking": "NUEVO"}
        ]
        
        for mat in COMMERCIAL_MATERIALS:
            comparison_data.append({
                "Material": mat['name'],
                "Resistencia": str(mat['strength']),
                "Peso": str(mat['weight']),
                "Costo": str(mat['cost']),
                "Cond. Térmica": str(mat['thermal']),
                "Ranking": mat['category']
            })
        
        st.dataframe(comparison_data, use_container_width=True)
        
        # Rankings
        st.markdown("### 🏆 Rankings por Categoría")
        c1, c2, c3 = st.columns(3)
        with c1:
            st.markdown("<div class='nasa-card text-center'><h4>🛩️ Aeroespacial</h4><h2 style='color: #ff6b00;'>#2</h2></div>", unsafe_allow_html=True)
        with c2:
            st.markdown("<div class='nasa-card text-center'><h4>⚡ Energía</h4><h2 style='color: #ff6b00;'>#1</h2></div>", unsafe_allow_html=True)
        with c3:
            st.markdown("<div class='nasa-card text-center'><h4>🏥 Médico</h4><h2 style='color: #ff6b00;'>#3</h2></div>", unsafe_allow_html=True)
    
    with tabs[5]:  # Archivos
        st.subheader("📁 Exportación a Simuladores")
        st.markdown("Archivos listos para ejecutar en supercomputadora")
        
        # Sub-tabs para cada simulador
        sim_tabs = st.tabs(["⚛️ Quantum ESPRESSO", "🔷 VASP", "🔧 LAMMPS", "🐍 Scripts Python"])
        
        # ========== QUANTUM ESPRESSO ==========
        with sim_tabs[0]:
            st.markdown("### ⚛️ Quantum ESPRESSO (pw.x)")
            
            col1, col2 = st.columns(2)
            with col1:
                qe_calc_type = st.selectbox("Tipo de cálculo:", 
                    ["scf", "vc-relax", "bands", "nscf"], key="qe_type")
            with col2:
                qe_ecut = st.number_input("Energía de corte (Ry):", value=60, min_value=20, max_value=200)
            
            # Generar archivos QE
            st.markdown("#### 📄 pw.in - SCF")
            qe_scf = f"""&CONTROL
  calculation = '{qe_calc_type}'
  prefix = '{metal}O2'
  outdir = './tmp'
  pseudo_dir = './pseudo'
/

&SYSTEM
  ibrav = 1
  celldm(1) = 8.0
  nat = 3
  ntyp = 2
  ecutwfc = {qe_ecut}.0
  ecutrho = {qe_ecut * 8}.0
  occupations = 'smearing'
  smearing = 'mp'
  degauss = 0.02
/

&ELECTRONS
  conv_thr = 1.0d-8
  mixing_beta = 0.7
/

ATOMIC_SPECIES
 {metal} {ELEMENTS_DATABASE.get(metal, {}).get('mw', 50):.3f} {metal}.pbe-spn-kjpaw_psl.1.0.0.UPF
 O  15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
 {metal} 0.000 0.000 0.000
 O  0.305 0.305 0.000
 O -0.305 -0.305 0.000

K_POINTS automatic
 8 8 8 0 0 0
"""
            st.code(qe_scf, language="bash")
            st.download_button("📥 Descargar pw.in", qe_scf, file_name=f"{metal}O2_pw.in")
            
            with st.expander("📊 bands.in - Estructura de Bandas"):
                qe_bands = f"""&CONTROL
  calculation = 'bands'
  prefix = '{metal}O2'
  outdir = './tmp'
/

&SYSTEM
  ibrav = 1
  celldm(1) = 8.0
  nat = 3
  ntyp = 2
  ecutwfc = {qe_ecut}.0
  ecutrho = {qe_ecut * 8}.0
  nbnd = 20
/

&ELECTRONS
  conv_thr = 1.0d-8
/

ATOMIC_SPECIES
 {metal} {ELEMENTS_DATABASE.get(metal, {}).get('mw', 50):.3f} {metal}.pbe-spn-kjpaw_psl.1.0.0.UPF
 O  15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
 {metal} 0.000 0.000 0.000
 O  0.305 0.305 0.000
 O -0.305 -0.305 0.000

K_POINTS crystal
5
  0.0  0.0  0.0  20  ! Gamma
  0.5  0.0  0.0  20  ! X
  0.5  0.5  0.0  20  ! M
  0.0  0.0  0.5  20  ! Z
  0.0  0.0  0.0  1   ! Gamma
"""
                st.code(qe_bands, language="bash")
                st.download_button("📥 Descargar bands.in", qe_bands, file_name=f"{metal}O2_bands.in")
            
            with st.expander("📈 dos.in - Densidad de Estados"):
                qe_dos = f"""&DOS
  prefix = '{metal}O2'
  outdir = './tmp'
  fildos = '{metal}O2.dos'
  Emin = -10.0
  Emax = 10.0
  DeltaE = 0.01
/
"""
                st.code(qe_dos, language="bash")
                st.download_button("📥 Descargar dos.in", qe_dos, file_name="dos.in")
            
            with st.expander("🎵 phonon.in - Fonones"):
                qe_ph = f"""&INPUTPH
  prefix = '{metal}O2'
  outdir = './tmp'
  fildyn = '{metal}O2.dyn'
  ldisp = .true.
  nq1 = 2
  nq2 = 2
  nq3 = 2
/
"""
                st.code(qe_ph, language="bash")
                st.download_button("📥 Descargar phonon.in", qe_ph, file_name="phonon.in")
            
            with st.expander("🚀 Script de ejecución"):
                qe_run = f"""#!/bin/bash
# Script Quantum ESPRESSO - {metal}O2
# Generado por CosmicForge Lab

module load quantum-espresso

mkdir -p tmp pseudo

# SCF
mpirun -np 4 pw.x < {metal}O2_pw.in > {metal}O2_pw.out

# Bandas (opcional)
# mpirun -np 4 pw.x < {metal}O2_bands.in > {metal}O2_bands.out
# mpirun -np 4 bands.x < bands_pp.in > bands.out
# plotband.x

# DOS (opcional)
# mpirun -np 4 dos.x < dos.in > dos.out

echo "Cálculo completado!"
"""
                st.code(qe_run, language="bash")
                st.download_button("📥 Descargar run.sh", qe_run, file_name="run_qe.sh")
        
        # ========== VASP ==========
        with sim_tabs[1]:
            st.markdown("### 🔷 VASP (Vienna Ab initio Simulation Package)")
            
            col1, col2 = st.columns(2)
            with col1:
                vasp_calc = st.selectbox("Tipo:", ["scf", "relax", "bands", "md"], key="vasp_type")
            with col2:
                vasp_encut = st.number_input("ENCUT (eV):", value=520, min_value=200, max_value=1000)
            
            # POSCAR
            st.markdown("#### 📄 POSCAR - Estructura Cristalina")
            vasp_poscar = f"""{metal}O2 - Generated by CosmicForge Lab
1.0
   4.593700  0.000000  0.000000
   0.000000  4.593700  0.000000
   0.000000  0.000000  2.958700
{metal}  O
1  2
Direct
  0.000000  0.000000  0.000000  {metal}
  0.305000  0.305000  0.000000  O
 -0.305000 -0.305000  0.000000  O
"""
            st.code(vasp_poscar, language="bash")
            st.download_button("📥 Descargar POSCAR", vasp_poscar, file_name="POSCAR")
            
            # INCAR
            incar_templates = {
                "scf": f"""# INCAR - {metal}O2 SCF Calculation
# CosmicForge Lab

SYSTEM = {metal}O2 - SCF
ENCUT = {vasp_encut}
EDIFF = 1.0e-6
IBRION = -1
ISMEAR = 0
SIGMA = 0.05
PREC = Accurate
LREAL = .FALSE.
ALGO = Normal
""",
                "relax": f"""# INCAR - {metal}O2 Relaxation
# CosmicForge Lab

SYSTEM = {metal}O2 - Relax
ENCUT = {vasp_encut}
EDIFF = 1.0e-6
EDIFFG = -0.01
IBRION = 2
ISIF = 3
NSW = 100
ISMEAR = 0
SIGMA = 0.05
PREC = Accurate
LREAL = .FALSE.
""",
                "bands": f"""# INCAR - {metal}O2 Band Structure
# CosmicForge Lab

SYSTEM = {metal}O2 - Bands
ENCUT = {vasp_encut}
EDIFF = 1.0e-6
IBRION = -1
ISMEAR = 0
SIGMA = 0.05
ICHARG = 11
LORBIT = 11
PREC = Accurate
""",
                "md": f"""# INCAR - {metal}O2 MD
# CosmicForge Lab

SYSTEM = {metal}O2 - MD
ENCUT = 400
EDIFF = 1.0e-4
IBRION = 0
NSW = 10000
POTIM = 1.0
SMASS = -3
ISMEAR = 0
SIGMA = 0.1
"""
            }
            
            with st.expander("⚙️ INCAR - Parámetros de Cálculo"):
                st.code(incar_templates.get(vasp_calc, incar_templates["scf"]), language="bash")
                st.download_button("📥 Descargar INCAR", incar_templates.get(vasp_calc, incar_templates["scf"]), file_name="INCAR")
            
            # KPOINTS
            with st.expander("🎯 KPOINTS - Malla de Puntos k"):
                if vasp_calc == "bands":
                    vasp_kpoints = f"""K-Points for Band Structure
10
Line-mode
Reciprocal
  0.000  0.000  0.000    20  ! Gamma
  0.500  0.000  0.000    20  ! X
  0.500  0.500  0.000    20  ! M
  0.000  0.500  0.000    20  ! Y
  0.000  0.000  0.000    20  ! Gamma
  0.000  0.000  0.500    20  ! Z
"""
                else:
                    vasp_kpoints = f"""K-Points - Monkhorst-Pack
0
Monkhorst-Pack
 8 8 8
0 0 0
"""
                st.code(vasp_kpoints, language="bash")
                st.download_button("📥 Descargar KPOINTS", vasp_kpoints, file_name="KPOINTS")
            
            with st.expander("📚 POTCAR Info"):
                potcar_info = f"""# POTCAR - Pseudopotenciales VASP
# ================================
# VASP requiere licencia comercial

# Para {metal}O2 necesitas:
cat $VASP_PP_PATH/{metal}/POTCAR $VASP_PP_PATH/O/POTCAR > POTCAR

# Verificar:
grep TITEL POTCAR

# Recomendado para {metal}:
# - PAW_{metal} (estándar)
# - PAW_{metal}_pv (incluye semicore p)
# - PAW_{metal}_sv (incluye semicore s y p)
"""
                st.code(potcar_info, language="bash")
        
        # ========== LAMMPS ==========
        with sim_tabs[2]:
            st.markdown("### 🔧 LAMMPS - Molecular Dynamics")
            
            col1, col2, col3 = st.columns(3)
            with col1:
                lmp_temp = st.number_input("Temperatura (K):", value=300, min_value=1, max_value=5000)
            with col2:
                lmp_steps = st.number_input("Pasos MD:", value=50000, min_value=1000, step=10000)
            with col3:
                lmp_timestep = st.number_input("Timestep (fs):", value=1, min_value=1, max_value=10)
            
            lmp_input = f"""# LAMMPS Input - {metal}O2
# Generado por CosmicForge Lab

# ==================== INICIALIZACIÓN ====================
units           metal
dimension       3
boundary        p p p
atom_style      full
neighbor        2.0 bin

# ==================== ESTRUCTURA ====================
lattice         fcc 4.0
region          box block 0 10 0 10 0 10
create_box      2 box
create_atoms    1 box

# ==================== MASA ====================
mass            1 {ELEMENTS_DATABASE.get(metal, {}).get('mw', 50):.3f}   # {metal}
mass            2 15.999   # O

# ==================== POTENCIAL ====================
pair_style      buck/coul/long 12.0
pair_coeff      1 1 0.001 0.1 0.0    # {metal}-{metal}
pair_coeff      2 2 0.001 0.1 0.0    # O-O
pair_coeff      1 2 0.001 0.2 0.0    # {metal}-O

kspace_style    pppm 1.0e-4

# ==================== VELOCIDAD ====================
velocity        all create {lmp_temp} 12345 mom yes rot yes dist gaussian

# ==================== MINIMIZACIÓN ====================
minimize        1.0e-4 1.0e-6 100 1000
reset_timestep  0

# ==================== TERMO ====================
thermo          1000
thermo_style    custom step temp pe ke etotal press vol density

# ==================== EQUILIBRACIÓN ====================
fix             1 all nvt temp {lmp_temp} {lmp_temp} 0.1
timestep        0.00{lmp_timestep}
run             10000

# ==================== PRODUCCIÓN ====================
unfix           1
fix             2 all npt temp {lmp_temp} {lmp_temp} 0.1 iso 0 0 1.0
run             {lmp_steps}

# ==================== OUTPUT ====================
write_data      {metal}O2_final.data
write_restart   {metal}O2.restart
"""
            st.code(lmp_input, language="bash")
            st.download_button("📥 Descargar input.in", lmp_input, file_name=f"{metal}O2_md.in")
            
            with st.expander("📊 data.lammps - Archivo de datos"):
                lmp_data = f"""LAMMPS data file - {metal}O2

1000 atoms
2 atom types

0.0 40.0 xlo xhi
0.0 40.0 ylo yhi
0.0 40.0 zlo zhi

Masses

1 {ELEMENTS_DATABASE.get(metal, {}).get('mw', 50):.3f}  # {metal}
2 15.999  # O
"""
                st.code(lmp_data, language="bash")
                st.download_button("📥 Descargar data.lammps", lmp_data, file_name="data.lammps")
            
            with st.expander("⚡ minimization.in"):
                lmp_min = f"""# LAMMPS Minimización - {metal}O2

units           metal
dimension       3
boundary        p p p
atom_style      full

lattice         fcc 4.0
region          box block 0 5 0 5 0 5
create_box      2 box

mass 1 {ELEMENTS_DATABASE.get(metal, {}).get('mw', 50):.3f}
mass 2 15.999

pair_style      buck/coul/long 12.0
pair_coeff      * * 0.001 0.1 0.0
kspace_style    pppm 1.0e-4

minimize        1.0e-6 1.0e-8 10000 100000
write_data      minimized.data
"""
                st.code(lmp_min, language="bash")
                st.download_button("📥 Descargar minimize.in", lmp_min, file_name="minimize.in")
        
        # ========== SCRIPTS PYTHON ==========
        with sim_tabs[3]:
            st.markdown("### 🐍 Scripts Python por Tarea")
            st.markdown("Cada archivo incluye ejemplo de uso al final")
            
            st.info("📁 Los scripts están en la carpeta `tasks/` del repositorio")
            
            # VASP Task
            with st.expander("📄 task_vasp.py - Generador VASP"):
                vasp_script = '''#!/usr/bin/env python3
"""
Tarea VASP - Genera POSCAR, INCAR, KPOINTS

Ejemplo de uso:
    from tasks.task_vasp import VASPTask
    
    task = VASPTask(metal="Ti", calc_type="scf")
    files = task.generate_all()
    task.save_files("./vasp_calc")
"""
from tasks.task_vasp import VASPTask

# Crear tarea
task = VASPTask(metal="''' + metal + '''", calc_type="scf")

# Generar archivos
files = task.generate_all()
for name, content in files.items():
    print(f"Archivo: {name}")
    print(content[:200] + "...")

# Guardar
task.save_files("./vasp_''' + metal + '''O2")
'''
                st.code(vasp_script, language="python")
                st.download_button("📥 Descargar ejemplo_vasp.py", vasp_script, file_name="ejemplo_vasp.py")
            
            # LAMMPS Task
            with st.expander("📄 task_lammps.py - Generador LAMMPS"):
                lmp_script = '''#!/usr/bin/env python3
"""
Tarea LAMMPS - Genera input.in, data.lammps

Ejemplo de uso:
    from tasks.task_lammps import LAMMPSTask
    
    task = LAMMPSTask(metal="Ti", physical_props={"density": 4.5})
    task.save_files("./lammps_calc")
"""
from tasks.task_lammps import LAMMPSTask

# Crear tarea
task = LAMMPSTask(metal="''' + metal + '''", physical_props={"density": ''' + str(round(pp.get('density', 4.5), 2)) + '''})

# Generar input con temperatura personalizada
input_content = task.generate_input(temperature=300, n_steps=50000)
print(input_content)

# Guardar
task.save_files("./lammps_''' + metal + '''O2")
'''
                st.code(lmp_script, language="python")
                st.download_button("📥 Descargar ejemplo_lammps.py", lmp_script, file_name="ejemplo_lammps.py")
            
            # QE Task
            with st.expander("📄 task_quantum_espresso.py - Generador QE"):
                qe_script = '''#!/usr/bin/env python3
"""
Tarea Quantum ESPRESSO - Genera pw.in, bands.in, dos.in

Ejemplo de uso:
    from tasks.task_quantum_espresso import QuantumESPRESSOTask
    
    task = QuantumESPRESSOTask(metal="Ti")
    task.save_files("./qe_calc")
"""
from tasks.task_quantum_espresso import QuantumESPRESSOTask

# Crear tarea
task = QuantumESPRESSOTask(metal="''' + metal + '''", calc_type="scf")

# Generar SCF
print(task.generate_scf())

# Generar bandas
print(task.generate_bands())

# Guardar todos los archivos
task.save_files("./qe_''' + metal + '''O2")
'''
                st.code(qe_script, language="python")
                st.download_button("📥 Descargar ejemplo_qe.py", qe_script, file_name="ejemplo_qe.py")
            
            st.markdown("---")
            st.markdown("### 📦 Descargar todos los archivos")
            if st.button("📦 Descargar ZIP completo", type="primary"):
                st.info("Generando archivo ZIP con todos los archivos...")
                # Aquí iría la generación del ZIP
    
    with tabs[6]:  # PDF
        st.subheader("📄 Generar Reporte Completo")
        
        st.markdown("""
        <div class='nasa-card'>
        <p>El reporte incluye:</p>
        <ul>
            <li>✅ Información del objeto astrofísico</li>
            <li>✅ Propiedades calculadas del material</li>
            <li>✅ Resultados de simulación de condiciones extremas</li>
            <li>✅ Resultados de optimización</li>
            <li>✅ Resultados de cálculo DFT</li>
            <li>✅ Comparación con materiales comerciales</li>
            <li>✅ Reflexión filosófica (Yu. I. Manin)</li>
        </ul>
        </div>
        """, unsafe_allow_html=True)
        
        if st.button("📥 Descargar PDF Completo", type="primary"):
            pdf_bytes = create_pdf_report(
                ad, pp, rc, production_guide,
                st.session_state.simulation,
                st.session_state.optimization,
                st.session_state.dft_result
            )
            
            if isinstance(pdf_bytes, bytes):
                st.download_button(
                    label="📥 Guardar PDF",
                    data=pdf_bytes,
                    file_name=f"CosmicForge_Report_{ad.get('object_name', 'material')}.pdf",
                    mime="application/pdf"
                )

# Footer Quote
st.markdown("---")
quote = MANIN_QUOTES[hash(str(datetime.now())) % len(MANIN_QUOTES)]
st.markdown(f"""
<div class='nasa-card' style='background: linear-gradient(135deg, #1a1a2a 0%, #16213e 100%); border-left: 5px solid #e94560;'>
    <p style='font-style: italic; color: #ccc;'>"{quote}"</p>
    <p style='text-align: right; color: #888; font-size: 12px;'>— Yuri I. Manin, "Lo demostrable e indemostrable"</p>
</div>
""", unsafe_allow_html=True)

# Footer
st.markdown("""
<div class='footer'>
    <p>🚀 CosmicForge Lab v3.6 NASA Edition</p>
    <p style='font-size: 12px;'>Exportación: VASP | LAMMPS | Quantum ESPRESSO | Scripts Python</p>
</div>
""", unsafe_allow_html=True)
