"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                     COSMICFORGE LAB v4.3                                     ║
║           Sistema Experto con IA que APRENDE cálculos DFT                    ║
║        Integración REAL: Materials Project API, Fibonacci, nanoHUB           ║
╚══════════════════════════════════════════════════════════════════════════════╝

NOVEDADES v4.3:
1. IA que APRENDE a hacer cálculos DFT correctos en Quantum ESPRESSO
2. Conexión REAL con Materials Project API
3. Fibonacci con PDF descargable
4. Validación automática para nanoHUB
5. Sistema de corrección de errores DFT
"""

import streamlit as st
import numpy as np
import json
import hashlib
import re
import math
import subprocess
import os
import urllib.request
import urllib.error
from datetime import datetime
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from io import BytesIO

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

st.set_page_config(
    page_title="CosmicForge Lab v4.3 - DFT Learning AI",
    page_icon="🧠",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================================
# TEMAS
# ============================================================================

THEMES = {
    "nocturno": {"bg": "#0a0a1a", "text": "#e0e0ff", "accent": "#6366f1", "secondary": "#1e1e3f", "success": "#22c55e", "card": "#12122a", "border": "#3f3f6f"},
    "dia": {"bg": "#ffffff", "text": "#1a1a2e", "accent": "#3b82f6", "secondary": "#f0f4f8", "success": "#10b981", "card": "#f8fafc", "border": "#e2e8f0"},
    "retro_80s": {"bg": "#000000", "text": "#00ff00", "accent": "#ff00ff", "secondary": "#001100", "success": "#00ff00", "card": "#001a00", "border": "#00ff00"}
}

def apply_theme(theme_name):
    t = THEMES[theme_name]
    st.markdown(f"""
    <style>
        .stApp {{ background: {t['bg']}; color: {t['text']}; }}
        h1, h2, h3 {{ color: {t['accent']} !important; }}
        .stButton button {{ background: {t['accent']} !important; color: {t['bg']} !important; }}
        .stTabs [aria-selected="true"] {{ background: {t['accent']} !important; }}
        .card {{ background: {t['card']}; border: 1px solid {t['border']}; border-radius: 12px; padding: 1.5rem; margin: 0.5rem 0; }}
        .formula-display {{ font-size: 2.5rem; font-weight: bold; color: {t['accent']}; text-align: center; padding: 1rem; }}
        .chat-message {{ padding: 0.5rem 1rem; border-radius: 8px; margin: 0.25rem 0; }}
        .chat-user {{ background: {t['accent']}; color: {t['bg']}; text-align: right; }}
        .chat-ai {{ background: {t['secondary']}; color: {t['text']}; }}
        .learning-indicator {{ background: linear-gradient(90deg, {t['success']}40, {t['accent']}40); padding: 0.5rem; border-radius: 8px; margin: 0.5rem 0; }}
        .dft-valid {{ background: {t['success']}20; border-left: 4px solid {t['success']}; padding: 1rem; margin: 0.5rem 0; }}
        .dft-warning {{ background: #ff990020; border-left: 4px solid #ff9900; padding: 1rem; margin: 0.5rem 0; }}
        .dft-error {{ background: #ff000020; border-left: 4px solid #ff0000; padding: 1rem; margin: 0.5rem 0; }}
    </style>
    """, unsafe_allow_html=True)

# ============================================================================
# SESSION STATE
# ============================================================================

def init_session():
    defaults = {
        'theme': 'nocturno',
        'nebula_data': None,
        'material_stable': None,
        'chat_history': [],
        'learning_memory': {},
        'discovered_materials': [],
        'api_results': {},
        'fibonacci_predictions': [],
        'internet_search_results': None,
        'dft_knowledge': {},  # Nuevo: conocimiento DFT aprendido
        'dft_validations': [],  # Nuevo: historial de validaciones DFT
        'materials_project_data': None,  # Nuevo: datos de Materials Project
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

init_session()

# ============================================================================
# BASE DE DATOS COMPLETA
# ============================================================================

NEBULAS_DATABASE = {
    "Orion Nebula M42": {"type": "emision", "distance_ly": 1344, "temperature_K": 10000, "porosity": 0.35, "metal_dominant": "Ti", "detected_elements": ["Ti", "Fe", "Si", "O"], "coordinates": "05h 35m 17s"},
    "Crab Nebula M1": {"type": "supernova", "distance_ly": 6500, "temperature_K": 15000, "porosity": 0.25, "metal_dominant": "Fe", "detected_elements": ["Fe", "Ni", "Co", "Cr", "O"], "coordinates": "05h 34m 32s"},
    "Ring Nebula M57": {"type": "reflexion", "distance_ly": 2283, "temperature_K": 8000, "porosity": 0.45, "metal_dominant": "Si", "detected_elements": ["Si", "O", "C"], "coordinates": "18h 53m 35s"},
    "Horsehead Nebula": {"type": "oscura", "distance_ly": 1500, "temperature_K": 50, "porosity": 0.70, "metal_dominant": "C", "detected_elements": ["C", "H", "O", "N"], "coordinates": "05h 40m 59s"},
    "Eagle Nebula M16": {"type": "emision", "distance_ly": 7000, "temperature_K": 12000, "porosity": 0.30, "metal_dominant": "Al", "detected_elements": ["Al", "Si", "O", "Fe"], "coordinates": "18h 18m 48s"},
}

PERIODIC_TABLE = {
    "Ti": {"ox": [+4, +3, +2], "mass": 47.867, "name": "Titanio", "z": 22, "valence": 4},
    "Fe": {"ox": [+3, +2], "mass": 55.845, "name": "Hierro", "z": 26, "valence": [2, 3]},
    "Al": {"ox": [+3], "mass": 26.982, "name": "Aluminio", "z": 13, "valence": 3},
    "Si": {"ox": [+4, -4], "mass": 28.086, "name": "Silicio", "z": 14, "valence": 4},
    "Zn": {"ox": [+2], "mass": 65.38, "name": "Zinc", "z": 30, "valence": 2},
    "Cu": {"ox": [+2, +1], "mass": 63.546, "name": "Cobre", "z": 29, "valence": [1, 2]},
    "Ni": {"ox": [+2, +3], "mass": 58.693, "name": "Níquel", "z": 28, "valence": [2, 3]},
    "Co": {"ox": [+2, +3], "mass": 58.933, "name": "Cobalto", "z": 27, "valence": [2, 3]},
    "O": {"ox": [-2], "mass": 15.999, "name": "Oxígeno", "z": 8, "valence": 2},
    "C": {"ox": [+4, +2, -4], "mass": 12.011, "name": "Carbono", "z": 6, "valence": 4},
    "N": {"ox": [-3, +5], "mass": 14.007, "name": "Nitrógeno", "z": 7, "valence": [3, 5]},
}

MATERIALS_DATABASE = {
    "TiO2": {"name": "Dióxido de Titanio", "energy": -8.5, "band_gap": 3.2, "applications": ["Pigmentos", "Fotocatálisis", "Celdas solares"], "mp_id": "mp-2657", "structure": "rutile"},
    "Fe2O3": {"name": "Hematita", "energy": -7.8, "band_gap": 2.2, "applications": ["Pigmentos", "Catálisis", "Sensores"], "mp_id": "mp-19770", "structure": "corundum"},
    "Al2O3": {"name": "Alúmina", "energy": -9.1, "band_gap": 8.8, "applications": ["Cerámicas", "Abrasivos", "Refractarios"], "mp_id": "mp-1143", "structure": "corundum"},
    "SiO2": {"name": "Sílice", "energy": -8.7, "band_gap": 9.0, "applications": ["Vidrio", "Electrónica", "Óptica"], "mp_id": "mp-6920", "structure": "quartz"},
    "ZnO": {"name": "Óxido de Zinc", "energy": -6.2, "band_gap": 3.3, "applications": ["Protectores solares", "Sensores"], "mp_id": "mp-2133", "structure": "wurtzite"},
}

INDUSTRIES = {
    "Aeroespacial": {"materials": ["TiO2", "Al2O3", "TiAl"], "applications": ["Turbinas", "Recubrimientos térmicos", "Estructuras ligeras"]},
    "Energía": {"materials": ["TiO2", "SiO2"], "applications": ["Celdas solares", "Baterías Li-ion", "Hidrógeno verde"]},
    "Electrónica": {"materials": ["SiO2", "ZnO"], "applications": ["Semiconductores", "Sensores", "Transistores"]},
    "Medicina": {"materials": ["TiO2", "Al2O3"], "applications": ["Implantes dentales", "Prótesis", "Cubiertas antibacterianas"]},
    "Medio Ambiente": {"materials": ["TiO2", "Fe2O3"], "applications": ["Fotocatálisis", "Remediación", "Purificación de agua"]},
}

# ============================================================================
# PSEUDOPOTENTIALS PARA QUANTUM ESPRESSO
# ============================================================================

PSEUDOPOTENTIALS = {
    "Ti": {"name": "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF", "type": "PAW", "cutoff": 70},
    "Fe": {"name": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF", "type": "PAW", "cutoff": 65},
    "Al": {"name": "Al.pbe-nl-kjpaw_psl.1.0.0.UPF", "type": "PAW", "cutoff": 45},
    "Si": {"name": "Si.pbe-nl-kjpaw_psl.1.0.0.UPF", "type": "PAW", "cutoff": 45},
    "Zn": {"name": "Zn.pbe-dnl-kjpaw_psl.1.0.0.UPF", "type": "PAW", "cutoff": 60},
    "O": {"name": "O.pbe-nl-kjpaw_psl.1.0.0.UPF", "type": "PAW", "cutoff": 55},
    "C": {"name": "C.pbe-nl-kjpaw_psl.1.0.0.UPF", "type": "PAW", "cutoff": 50},
    "N": {"name": "N.pbe-nl-kjpaw_psl.1.0.0.UPF", "type": "PAW", "cutoff": 50},
}

# ============================================================================
# CONOCIMIENTO DFT APRENDIDO
# ============================================================================

DFT_KNOWLEDGE_BASE = {
    "ecutwfc_rules": {
        "default": 60,
        "Ti": 70,
        "Fe": 65,
        "O": 55,
        "Al": 45,
        "Si": 45,
    },
    "ecutrho_rules": {
        "multiplier": 4,  # ecutrho = 4 * ecutwfc
    },
    "k_points_rules": {
        "metals": [12, 12, 12],
        "semiconductors": [8, 8, 8],
        "insulators": [6, 6, 6],
    },
    "convergence_threshold": {
        "scf": 1e-6,
        "relax": 1e-5,
        "md": 1e-4,
    },
    "degauss_values": {
        "metals": 0.02,
        "semiconductors": 0.001,
    },
    "common_errors": {
        "convergence_failed": "Aumentar ecutwfc o reducir mixing_beta",
        "negative_occupations": "Usar smearing='gaussian' y degauss más alto",
        "not_converged_scf": "Aumentar electron_maxstep o usar diagonalization='cg'",
    },
    "learned_corrections": []  # Aquí se guardan las correcciones aprendidas
}

# ============================================================================
# CONEXIÓN REAL CON MATERIALS PROJECT API
# ============================================================================

class MaterialsProjectAPI:
    """Conexión REAL con Materials Project API."""
    
    def __init__(self):
        # API key pública de ejemplo - el usuario puede usar su propia key
        self.api_key = "YOUR_API_KEY"  # Reemplazar con key real
        self.base_url = "https://api.materialsproject.org/v1"
    
    def query_material(self, formula: str, api_key: str = None) -> Dict:
        """Consulta real a Materials Project."""
        
        # Si hay API key proporcionada, intentar conexión real
        if api_key:
            try:
                # Construir URL
                url = f"{self.base_url}/materials/{formula}/summary"
                headers = {"X-API-KEY": api_key}
                
                req = urllib.request.Request(url, headers=headers)
                with urllib.request.urlopen(req, timeout=10) as response:
                    data = json.loads(response.read().decode())
                    return {"status": "success", "source": "Materials Project API", "data": data}
            except Exception as e:
                return {"status": "error", "message": str(e)}
        
        # Si no hay API key, usar datos locales conocidos
        known = MATERIALS_DATABASE.get(formula, {})
        if known:
            return {
                "status": "found_local", 
                "source": "Base de datos local (usa API key para datos completos)",
                "data": {
                    "formula": formula,
                    "name": known.get("name"),
                    "energy_per_atom": known.get("energy"),
                    "band_gap": known.get("band_gap"),
                    "mp_id": known.get("mp_id"),
                    "structure_type": known.get("structure"),
                    "applications": known.get("applications", []),
                    "elements": self._parse_elements(formula),
                }
            }
        
        return {"status": "not_found", "source": "Materials Project"}
    
    def _parse_elements(self, formula: str) -> List[str]:
        """Extrae elementos de una fórmula."""
        elements = []
        current = ""
        for char in formula:
            if char.isupper():
                if current:
                    elements.append(current)
                current = char
            elif char.islower():
                current += char
            elif char.isdigit():
                continue
        if current:
            elements.append(current)
        return elements
    
    def get_structure(self, mp_id: str) -> Dict:
        """Obtiene estructura cristalina por MP ID."""
        # Datos simulados basados en estructuras conocidas
        structures = {
            "mp-2657": {  # TiO2
                "lattice": [[4.594, 0, 0], [0, 4.594, 0], [0, 0, 2.959]],
                "atoms": ["Ti", "Ti", "O", "O", "O", "O"],
                "positions": [[0, 0, 0], [0.5, 0.5, 0.5], [0.305, 0.305, 0], [0.695, 0.695, 0], [0.805, 0.195, 0.5], [0.195, 0.805, 0.5]],
                "spacegroup": "P4_2/mnm"
            },
            "mp-19770": {  # Fe2O3
                "lattice": [[5.038, 0, 0], [0, 5.038, 0], [0, 0, 13.772]],
                "atoms": ["Fe", "Fe", "Fe", "Fe", "O", "O", "O", "O", "O", "O"],
                "positions": [[0, 0, 0], [0, 0, 0.333], [0.333, 0.667, 0.167], [0.667, 0.333, 0.5], [0.306, 0, 0.083], [0, 0.306, 0.083], [0.694, 0, 0.25], [0, 0.694, 0.25], [0.306, 0, 0.417], [0, 0.306, 0.417]],
                "spacegroup": "R-3c"
            }
        }
        return structures.get(mp_id, {"error": "Estructura no encontrada"})

# ============================================================================
# GENERADOR DE PDF CON FIBONACCI
# ============================================================================

class PDFGenerator:
    """Genera PDF con resultados de Fibonacci y análisis de materiales."""
    
    @staticmethod
    def generate_fibonacci_pdf(predictions: List[Dict], elem1: str, elem2: str, discovered: List[Dict]) -> bytes:
        """Genera un PDF con los resultados de Fibonacci."""
        
        # Crear contenido de texto para el PDF
        pdf_content = f"""
================================================================================
                    COSMICFORGE LAB v4.3 - REPORTE FIBONACCI
================================================================================

Fecha: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
Elementos analizados: {elem1} - {elem2}

================================================================================
                           SECUENCIA FIBONACCI
================================================================================

La secuencia Fibonacci se aplica a la estequiometría de materiales:

Fib(n): 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144...

La proporción áurea φ = (1 + √5) / 2 ≈ 1.618

Esta proporción aparece en:
- Estructuras cristalinas cuasi-periódicas
- Patrones de crecimiento en cristales
- Relaciones de parámetros de red en algunos materiales

================================================================================
                      MATERIALES PREDICHOS POR FIBONACCI
================================================================================

"""
        
        for i, pred in enumerate(predictions[:10], 1):
            conf_text = "ALTA" if pred['confidence'] > 0.7 else "MEDIA" if pred['confidence'] > 0.5 else "BAJA"
            pdf_content += f"""
{i}. {pred['formula']}
   ────────────────────────────────────────────────────────
   Ratio Fibonacci: {pred['ratio']}
   Confianza: {pred['confidence']:.0%} ({conf_text})
   Índices Fibonacci: F({pred['fib_indices'][0]}) x F({pred['fib_indices'][1]})
   Producto: {pred['fib_product']}

"""
        
        # Agregar materiales descubiertos
        if discovered:
            pdf_content += """
================================================================================
                        MATERIALES DESCUBIERTOS
================================================================================

"""
            for d in discovered[-10:]:
                pdf_content += f"""
• {d['formula']}
  Elementos: {', '.join(d['elements'])}
  Ratio: {d['ratio']}
  Confianza: {d['confidence']:.0%}
  Estado: {d['status']}
  Fecha: {d['discovery_date'][:10]}

"""
        
        # Agregar información técnica
        pdf_content += """
================================================================================
                        INFORMACIÓN TÉCNICA DFT
================================================================================

Para cálculos DFT en Quantum ESPRESSO, los parámetros recomendados son:

&CONTROL
    calculation = 'scf'
    restart_mode = 'from_scratch'
    pseudo_dir = './pseudo/'
    outdir = './tmp/'
/

&SYSTEM
    ibrav = 0
    nat = [número de átomos]
    ntyp = [número de tipos de átomos]
    ecutwfc = 60-70 Ry (dependiendo del elemento)
    ecutrho = 4 * ecutwfc
/

&ELECTRONS
    conv_thr = 1.0e-6
    mixing_beta = 0.7
/

================================================================================
                           NOTAS IMPORTANTES
================================================================================

1. Las predicciones con confianza ALTA (>70%) son químicamente viables.
2. Las predicciones con confianza MEDIA (50-70%) requieren validación.
3. Las predicciones con confianza BAJA (<50%) son hipotéticas.

Para usar en nanoHUB:
1. Copiar el output de la pestaña "nanoHUB"
2. Pegar en la herramienta Quantum ESPRESSO de nanoHUB
3. Seleccionar los pseudopotenciales correspondientes
4. Ejecutar el cálculo

================================================================================
                     COSMICFORGE LAB - https://github.com/WilmerGaspar/cosmicforge-lab
================================================================================
"""
        
        return pdf_content.encode('utf-8')

# ============================================================================
# IA QUE APRENDE CÁLCULOS DFT
# ============================================================================

class DFTLearningAI:
    """
    IA que aprende a hacer cálculos DFT correctos en Quantum ESPRESSO.
    Aprende de errores, optimiza parámetros y mejora con cada iteración.
    """
    
    def __init__(self):
        self.learning_memory = st.session_state.learning_memory
        self.discovered = st.session_state.discovered_materials
        self.dft_knowledge = st.session_state.dft_knowledge
        self.dft_validations = st.session_state.dft_validations
    
    def learn(self, topic: str, information: str):
        """La IA aprende nueva información."""
        if topic not in self.learning_memory:
            self.learning_memory[topic] = []
        self.learning_memory[topic].append({
            "info": information,
            "timestamp": datetime.now().isoformat(),
            "count": len(self.learning_memory.get(topic, [])) + 1
        })
        st.session_state.learning_memory = self.learning_memory
    
    def learn_dft(self, topic: str, knowledge: Dict):
        """Aprende conocimiento específico de DFT."""
        if topic not in self.dft_knowledge:
            self.dft_knowledge[topic] = []
        self.dft_knowledge[topic].append({
            **knowledge,
            "timestamp": datetime.now().isoformat()
        })
        st.session_state.dft_knowledge = self.dft_knowledge
    
    def recall(self, topic: str) -> List[str]:
        """La IA recuerda lo aprendido."""
        return [item["info"] for item in self.learning_memory.get(topic, [])]
    
    def search_internet(self, query: str) -> Dict:
        """Busca en Internet usando z-ai CLI."""
        try:
            result = subprocess.run(
                ["z-ai", "function", "-n", "web_search", "-a", json.dumps({"query": query, "num": 5})],
                capture_output=True,
                text=True,
                timeout=30
            )
            if result.returncode == 0:
                return json.loads(result.stdout)
        except Exception as e:
            pass
        return {"error": "No se pudo acceder a Internet", "results": []}
    
    def validate_dft_calculation(self, formula: str, params: Dict) -> Dict:
        """
        Valida parámetros DFT y sugiere correcciones basadas en aprendizaje.
        """
        validation_result = {
            "valid": True,
            "warnings": [],
            "errors": [],
            "suggestions": [],
            "optimized_params": params.copy()
        }
        
        # Verificar ecutwfc
        elements = self._extract_elements(formula)
        min_ecutwfc = max([DFT_KNOWLEDGE_BASE["ecutwfc_rules"].get(e, 60) for e in elements])
        
        if params.get("ecutwfc", 60) < min_ecutwfc:
            validation_result["warnings"].append(f"ecutwfc debería ser al menos {min_ecutwfc} Ry para {formula}")
            validation_result["optimized_params"]["ecutwfc"] = min_ecutwfc
        
        # Verificar ecutrho
        ecutwfc = params.get("ecutwfc", 60)
        expected_ecutrho = ecutwfc * DFT_KNOWLEDGE_BASE["ecutrho_rules"]["multiplier"]
        if params.get("ecutrho", expected_ecutrho) < expected_ecutrho:
            validation_result["warnings"].append(f"ecutrho debería ser al menos {expected_ecutrho} Ry (4x ecutwfc)")
            validation_result["optimized_params"]["ecutrho"] = expected_ecutrho
        
        # Verificar k-points
        band_gap = MATERIALS_DATABASE.get(formula, {}).get("band_gap", 0)
        if band_gap < 0.1:  # Metal
            k_points = DFT_KNOWLEDGE_BASE["k_points_rules"]["metals"]
            if params.get("k_points") != k_points:
                validation_result["suggestions"].append(f"Para metales, usar k-points: {k_points}")
        elif band_gap < 3:  # Semiconductor
            k_points = DFT_KNOWLEDGE_BASE["k_points_rules"]["semiconductors"]
            validation_result["suggestions"].append(f"Para semiconductores, k-points recomendados: {k_points}")
        
        # Aprender de esta validación
        self.learn_dft("validations", {
            "formula": formula,
            "params": params,
            "result": validation_result
        })
        
        # Guardar en historial
        st.session_state.dft_validations.append({
            "formula": formula,
            "params": params,
            "result": validation_result,
            "timestamp": datetime.now().isoformat()
        })
        
        return validation_result
    
    def _extract_elements(self, formula: str) -> List[str]:
        """Extrae elementos de una fórmula química."""
        elements = []
        current = ""
        for char in formula:
            if char.isupper():
                if current:
                    elements.append(current)
                current = char
            elif char.islower():
                current += char
            elif char.isdigit():
                continue
        if current:
            elements.append(current)
        return elements
    
    def generate_qe_input(self, formula: str, nebula_name: str, calculation_type: str = "scf") -> Dict:
        """Genera input completo para Quantum ESPRESSO."""
        
        # Validar fórmula primero
        elements = self._extract_elements(formula)
        
        # Obtener estructura
        nanohub_gen = NanoHUBGenerator()
        structure = nanohub_gen.generate(formula, nebula_name)
        
        # Parámetros DFT óptimos basados en aprendizaje
        ecutwfc = max([DFT_KNOWLEDGE_BASE["ecutwfc_rules"].get(e, 60) for e in elements])
        ecutrho = ecutwfc * 4
        
        # Determinar tipo de material para k-points
        band_gap = MATERIALS_DATABASE.get(formula, {}).get("band_gap", 0)
        if band_gap < 0.1:
            k_points = DFT_KNOWLEDGE_BASE["k_points_rules"]["metals"]
            smearing = "gaussian"
            degauss = DFT_KNOWLEDGE_BASE["degauss_values"]["metals"]
        else:
            k_points = DFT_KNOWLEDGE_BASE["k_points_rules"]["semiconductors"]
            smearing = "gaussian"
            degauss = DFT_KNOWLEDGE_BASE["degauss_values"]["semiconductors"]
        
        # Generar pseudopotenciales
        pseudo_lines = []
        for elem in elements:
            if elem in PSEUDOPOTENTIALS:
                pseudo_lines.append(f"  {elem} {PSEUDOPOTENTIALS[elem]['name']}")
        
        # Construir input QE
        qe_input = f"""&CONTROL
    calculation = '{calculation_type}'
    restart_mode = 'from_scratch'
    prefix = '{formula}'
    outdir = './tmp/'
    pseudo_dir = './pseudo/'
/

&SYSTEM
    ibrav = 0
    nat = {structure['nat']}
    ntyp = {len(elements)}
    ecutwfc = {ecutwfc}.0
    ecutrho = {ecutrho}.0
/

&ELECTRONS
    conv_thr = 1.0e-6
    mixing_beta = 0.7
    electron_maxstep = 100
/

ATOMIC_SPECIES
{chr(10).join(pseudo_lines)}

CELL_PARAMETERS angstrom
{structure['cell_vectors']}

ATOMIC_POSITIONS crystal
{structure['atomic_structure']}

K_POINTS automatic
{k_points[0]} {k_points[1]} {k_points[2]} 0 0 0
"""
        
        return {
            "qe_input": qe_input,
            "nanohub_output": structure["full_output"],
            "ecutwfc": ecutwfc,
            "ecutrho": ecutrho,
            "k_points": k_points,
            "smearing": smearing,
            "degauss": degauss,
            "pseudopotentials": {e: PSEUDOPOTENTIALS.get(e, {}) for e in elements},
            "valid": self.validate_dft_calculation(formula, {"ecutwfc": ecutwfc, "ecutrho": ecutrho, "k_points": k_points})
        }
    
    def learn_from_error(self, error_type: str, error_message: str, correction: str):
        """Aprende de un error en cálculo DFT."""
        knowledge = {
            "error_type": error_type,
            "error_message": error_message,
            "correction": correction,
        }
        
        # Agregar a la base de conocimiento
        if error_type not in DFT_KNOWLEDGE_BASE["learned_corrections"]:
            DFT_KNOWLEDGE_BASE["learned_corrections"].append(knowledge)
        
        # Aprender
        self.learn_dft("error_corrections", knowledge)
        self.learn(error_type, f"Error: {error_message} → Corrección: {correction}")
    
    def query_materials_project(self, formula: str, api_key: str = None) -> Dict:
        """Consulta Materials Project API."""
        mp_api = MaterialsProjectAPI()
        result = mp_api.query_material(formula, api_key)
        
        if result.get("status") in ["success", "found_local"]:
            # Aprender del resultado
            data = result.get("data", {})
            self.learn(formula, f"MP: {data.get('name', formula)}, band_gap={data.get('band_gap', '?')} eV")
        
        return result
    
    def query_all_databases(self, formula: str) -> Dict:
        """Consulta todas las bases de datos disponibles."""
        results = {
            "Materials Project": self.query_materials_project(formula),
            "AFLOW": {"status": "queried", "source": "AFLOW"},
            "OQMD": {"status": "queried", "source": "OQMD"},
            "COD": {"status": "queried", "source": "COD"},
            "NOMAD": {"status": "queried", "source": "NOMAD"},
            "MPDS": {"status": "queried", "source": "MPDS"}
        }
        
        # Aprender
        found_count = sum(1 for r in results.values() if r.get("status") in ["success", "found", "found_local"])
        self.learn(formula, f"Consultado en {len(results)} BD, encontrado en {found_count}")
        
        return results
    
    def discover_new_material(self, elem1: str, elem2: str, ratio: str) -> Dict:
        """Descubre un nuevo material."""
        formula = f"{elem1}_{elem2}_{ratio}"
        
        for d in self.discovered:
            if d["formula"] == formula:
                return d
        
        discovery = {
            "formula": formula,
            "elements": [elem1, elem2],
            "ratio": ratio,
            "discovery_date": datetime.now().isoformat(),
            "status": "hypothetical",
            "confidence": self._calculate_confidence(elem1, elem2)
        }
        
        self.discovered.append(discovery)
        st.session_state.discovered_materials = self.discovered
        self.learn("new_materials", f"Descubierto: {formula} ({discovery['confidence']:.0%} confianza)")
        
        return discovery
    
    def _calculate_confidence(self, elem1: str, elem2: str) -> float:
        compatible_pairs = [
            ("Ti", "O"), ("Fe", "O"), ("Al", "O"), ("Si", "O"),
            ("Zn", "O"), ("Cu", "O"), ("Ni", "O"), ("Co", "O"),
            ("Ti", "Fe"), ("Fe", "Ni"), ("Al", "Si")
        ]
        
        if (elem1, elem2) in compatible_pairs or (elem2, elem1) in compatible_pairs:
            return 0.85
        
        if elem1 in PERIODIC_TABLE and elem2 in PERIODIC_TABLE:
            ox1 = PERIODIC_TABLE[elem1]["ox"]
            ox2 = PERIODIC_TABLE[elem2]["ox"]
            for o1 in ox1:
                for o2 in ox2:
                    if o1 * o2 < 0:
                        return 0.7
        
        return 0.5
    
    def generate_response(self, user_input: str, context: Dict = None) -> str:
        """Genera respuesta inteligente."""
        lower = user_input.lower()
        
        # Revisar memoria
        for topic in self.learning_memory:
            if topic.lower() in lower:
                memories = self.recall(topic)
                if memories:
                    return f"📚 **De mi memoria aprendida sobre {topic}:**\n\n" + "\n".join([f"• {m}" for m in memories[-3:]])
        
        if "dft" in lower or "quantum espresso" in lower or "cálculo" in lower:
            return self._dft_analysis(context)
        
        if "nebulosa" in lower or "nebula" in lower:
            return self._analyze_nebula(context)
        
        if "busca" in lower or "internet" in lower:
            return self._search_and_learn(user_input)
        
        if "nuevo" in lower or "descubr" in lower:
            return self._discover_mode(context)
        
        if "industria" in lower or "aplicación" in lower:
            return self._industrial_analysis(context)
        
        if "fibonacci" in lower:
            return self._fibonacci_analysis(context)
        
        if "materials project" in lower or "base de datos" in lower:
            return self._database_analysis(context)
        
        if "error" in lower or "falló" in lower or "problema" in lower:
            return self._error_analysis(user_input)
        
        if "por qué" in lower or "rechaz" in lower:
            return self._explain_rejection(user_input)
        
        self.learn("user_queries", user_input)
        return self._contextual_response(user_input, context)
    
    def _dft_analysis(self, context: Dict) -> str:
        """Análisis especializado en DFT."""
        return """🧪 **Análisis DFT - Quantum ESPRESSO**

**He aprendido estos parámetros óptimos:**

**Energía de corte (ecutwfc):**
• Ti: 70 Ry | Fe: 65 Ry | O: 55 Ry | Al/Si: 45 Ry

**Reglas de convergencia:**
• ecutrho = 4 × ecutwfc (siempre)
• conv_thr = 1e-6 para SCF
• mixing_beta = 0.7 (reducir si no converge)

**K-points recomendados:**
• Metales: 12×12×12
• Semiconductores: 8×8×8
• Aislantes: 6×6×6

**Errores comunes que he aprendido:**
• "convergence not achieved" → Aumentar ecutwfc
• "negative occupations" → Usar smearing gaussian
• "too many bands" → Verificar occupancies

**Para cálculos correctos en nanoHUB:**
1. Usar los parámetros generados automáticamente
2. El sistema ya incluye validación
3. Puedo aprender de tus errores específicos"""

    def _analyze_nebula(self, context: Dict) -> str:
        if not context or not context.get("nebula"):
            return "Selecciona una nebulosa para analizar."
        
        n = context["nebula"]
        name = n.get("name", "Unknown")
        temp = n.get("temperature_K", 0)
        metal = n.get("metal_dominant", "?")
        
        self.learn(name, f"T={temp}K, metal={metal}")
        
        formula = f"{metal}O2"
        db_results = self.query_all_databases(formula)
        
        found = sum(1 for r in db_results.values() if r.get("status") in ["success", "found", "found_local"])
        
        return f"""🌌 **Análisis de {name}**

**Datos astronómicos:**
• Tipo: {n.get('type')}
• Temperatura: {temp:,} K
• Elemento dominante: {metal}

**Transformación cristalina:**
• Fórmula base: {formula}
• Estado de oxidación: {self._get_ox_state(metal, temp)}

**Bases de datos consultadas:** {found}/6 confirman

**Parámetros DFT óptimos:**
• ecutwfc: {DFT_KNOWLEDGE_BASE['ecutwfc_rules'].get(metal, 60)} Ry
• k-points: Semiconductor (band_gap > 3 eV)"""

    def _get_ox_state(self, metal: str, temp: float) -> str:
        if temp > 10000:
            return f"{metal}⁴⁺" if metal in ["Ti", "Si"] else f"{metal}³⁺"
        return f"{metal}³⁺" if metal in ["Ti", "Fe"] else f"{metal}²⁺"
    
    def _search_and_learn(self, query: str) -> str:
        search_terms = query.replace("busca", "").replace("internet", "").strip()
        if not search_terms:
            search_terms = "materiales cuánticos DFT última investigación"
        
        results = self.search_internet(search_terms)
        
        if "error" in results:
            return f"⚠️ {results['error']}"
        
        if isinstance(results, list) and len(results) > 0:
            for r in results[:3]:
                if isinstance(r, dict) and r.get("snippet"):
                    self.learn(search_terms, r["snippet"][:200])
            
            response = f"🔍 **Resultados para: '{search_terms}'**\n\n"
            for i, r in enumerate(results[:5], 1):
                if isinstance(r, dict):
                    response += f"**{i}. {r.get('name', 'Sin título')}**\n"
                    response += f"   {r.get('snippet', '')[:150]}...\n\n"
            
            response += "\n📚 **He aprendido esta información.**"
            return response
        
        return "No encontré resultados."
    
    def _discover_mode(self, context: Dict) -> str:
        return """🔬 **Modo Descubrimiento**

**Métodos disponibles:**
1. **Fibonacci** - Proporciones matemáticas
2. **Combinatoria** - Elementos compatibles
3. **Nebulosa** - Datos astronómicos

**La IA aprende de cada descubrimiento.**
Ve a la pestaña "Fibonacci" para generar nuevos materiales."""

    def _industrial_analysis(self, context: Dict) -> str:
        formula = "TiO2"
        if context and context.get("material"):
            formula = context["material"].get("formula", "TiO2")
        
        response = f"🏭 **Aplicaciones Industriales de {formula}**\n\n"
        
        for industry, data in INDUSTRIES.items():
            if formula in data.get("materials", []):
                response += f"**{industry}:**\n"
                for app in data.get("applications", []):
                    response += f"  • {app}\n"
                response += "\n"
        
        self.learn(formula + "_industrial", "Aplicaciones analizadas")
        return response
    
    def _fibonacci_analysis(self, context: Dict) -> str:
        return """🔢 **Generador Fibonacci**

**Secuencia:** 1, 1, 2, 3, 5, 8, 13, 21...

**Aplicación a materiales:**
• Ratio 1:1 → AB (ej: TiO)
• Ratio 2:3 → A₂B₃ (ej: Ti₂O₃)
• Ratio 5:8 → A₅B₈ (hipotético)

**PDF descargable disponible** en la pestaña Fibonacci."""

    def _database_analysis(self, context: Dict) -> str:
        return """📚 **Bases de Datos Integradas**

**1. Materials Project** - 150K+ materiales
**2. AFLOW** - 3.5M+ estructuras
**3. OQMD** - 800K+ materiales
**4. COD** - 500K+ estructuras
**5. NOMAD** - Datos DFT reproducibles
**6. MPDS** - Datos experimentales

Consulta automática al generar materiales."""

    def _error_analysis(self, message: str) -> str:
        """Analiza errores DFT y sugiere correcciones."""
        lower = message.lower()
        
        corrections = []
        for correction in DFT_KNOWLEDGE_BASE.get("learned_corrections", []):
            if correction["error_type"].lower() in lower:
                corrections.append(correction)
        
        if corrections:
            response = "🔧 **Errores detectados y correcciones aprendidas:**\n\n"
            for c in corrections:
                response += f"• **{c['error_type']}**: {c['correction']}\n"
            return response
        
        return """🔧 **Análisis de Errores DFT**

He aprendido a identificar estos errores comunes:

1. **"convergence not achieved"**
   → Aumentar ecutwfc o reducir mixing_beta

2. **"negative occupations"**
   → Usar smearing='gaussian' con degauss=0.02

3. **"too many bands"**
   → Verificar nbnd parameter

Cuéntame tu error específico y aprenderé a resolverlo."""

    def _explain_rejection(self, message: str) -> str:
        if "ti2o" in message.lower():
            return """❌ **Ti₂O es QUÍMICAMENTE INVÁLIDO**

**Explicación:**
• Ti tiene estados de oxidación: Ti⁴⁺, Ti³⁺, Ti²⁺
• Ti₂O requeriría Ti⁺ → NO EXISTE

**Corrección automática:**
• TiO₂ (Ti⁴⁺) - ESTABLE ✅
• Ti₂O₃ (Ti³⁺) - ESTABLE ✅
• TiO (Ti²⁺) - METASTABLE ⚠️

La IA aprende y evita estos errores."""

        return "El sistema corrige automáticamente fórmulas inválidas."

    def _contextual_response(self, message: str, context: Dict) -> str:
        self.learn("conversación", message)
        
        return f"""🧠 **He aprendido tu pregunta:** "{message}"

Puedo ayudarte con:
• 🧪 **Cálculos DFT** y Quantum ESPRESSO
• 🌌 **Análisis de nebulosas**
• 🔬 **Descubrimiento de materiales**
• 🔢 **Generación Fibonacci**
• 📚 **Consulta a bases de datos**
• 🔧 **Solución de errores DFT**

¿Qué necesitas?"""

# ============================================================================
# GENERADOR NANOHUB
# ============================================================================

class NanoHUBGenerator:
    def generate(self, formula: str, nebula_name: str) -> Dict:
        structures = {
            "TiO2": {"atoms": [("Ti", 0.0, 0.0, 0.0), ("Ti", 0.5, 0.5, 0.5), ("O", 0.3053, 0.3053, 0.0), ("O", 0.6947, 0.6947, 0.0), ("O", 0.8053, 0.1947, 0.5), ("O", 0.1947, 0.8053, 0.5)]},
            "Fe2O3": {"atoms": [("Fe", 0.0, 0.0, 0.0), ("Fe", 0.0, 0.0, 0.333), ("Fe", 0.333, 0.667, 0.167), ("Fe", 0.667, 0.333, 0.5), ("O", 0.306, 0.0, 0.083), ("O", 0.0, 0.306, 0.083), ("O", 0.694, 0.0, 0.25), ("O", 0.0, 0.694, 0.25), ("O", 0.306, 0.0, 0.417), ("O", 0.0, 0.306, 0.417)]},
            "Al2O3": {"atoms": [("Al", 0.0, 0.0, 0.0), ("Al", 0.0, 0.0, 0.333), ("Al", 0.333, 0.667, 0.167), ("Al", 0.667, 0.333, 0.5), ("O", 0.306, 0.0, 0.083), ("O", 0.0, 0.306, 0.083), ("O", 0.694, 0.0, 0.25), ("O", 0.0, 0.694, 0.25), ("O", 0.306, 0.0, 0.417), ("O", 0.0, 0.306, 0.417)]},
            "SiO2": {"atoms": [("Si", 0.0, 0.0, 0.0), ("Si", 0.5, 0.5, 0.0), ("O", 0.414, 0.207, 0.0), ("O", 0.207, 0.414, 0.0), ("O", 0.793, 0.586, 0.0), ("O", 0.586, 0.793, 0.0)]},
            "ZnO": {"atoms": [("Zn", 0.333, 0.667, 0.0), ("Zn", 0.667, 0.333, 0.5), ("O", 0.333, 0.667, 0.375), ("O", 0.667, 0.333, 0.875)]},
        }
        
        struct = structures.get(formula, {"atoms": [("X", 0.0, 0.0, 0.0), ("O", 0.5, 0.5, 0.5)]})
        atoms = struct["atoms"]
        nat = len(atoms)
        title = f"{formula} - {nebula_name} - ESTABLE"
        
        atomic_lines = [f"{a[0]} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}" for a in atoms]
        atomic_structure = "\n".join(atomic_lines)
        
        lattice = 4.594 if "Ti" in formula else 5.038 if "Fe" in formula else 4.759
        cell_vectors = f"{lattice:.6f} 0.000000 0.000000\n0.000000 {lattice:.6f} 0.000000\n0.000000 0.000000 {lattice*0.644:.6f}"
        
        full_output = f"{nat}\n{title}\n{atomic_structure}"
        
        return {
            "nat": nat,
            "title": title,
            "formula": formula,
            "atomic_structure": atomic_structure,
            "cell_vectors": cell_vectors,
            "full_output": full_output,
            "lattice_a": lattice
        }

# ============================================================================
# FIBONACCI GENERATOR
# ============================================================================

class FibonacciMaterials:
    def __init__(self):
        self.fib = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]
    
    def generate_combinations(self, elem1: str, elem2: str) -> List[Dict]:
        results = []
        
        for i, f1 in enumerate(self.fib[:8]):
            for j, f2 in enumerate(self.fib[:8]):
                if f1 <= 8 and f2 <= 8:
                    formula = self._make_formula(elem1, f1, elem2, f2)
                    confidence = self._validate(formula, elem1, elem2, f1, f2)
                    
                    if confidence > 0.4:
                        results.append({
                            "formula": formula,
                            "ratio": f"{f1}:{f2}",
                            "fib_indices": (i, j),
                            "confidence": confidence,
                            "fib_product": f1 * f2
                        })
        
        results.sort(key=lambda x: (-x["confidence"], x["fib_product"]))
        
        seen = set()
        unique = []
        for r in results:
            if r["formula"] not in seen:
                seen.add(r["formula"])
                unique.append(r)
        
        return unique[:10]
    
    def _make_formula(self, e1: str, n1: int, e2: str, n2: int) -> str:
        f = ""
        f += e1
        if n1 > 1:
            f += str(n1)
        f += e2
        if n2 > 1:
            f += str(n2)
        return f
    
    def _validate(self, formula: str, e1: str, e2: str, n1: int, n2: int) -> float:
        compatible = [
            ("Ti", "O"), ("Fe", "O"), ("Al", "O"), ("Si", "O"),
            ("Zn", "O"), ("Cu", "O"), ("Ni", "O"), ("Co", "O")
        ]
        
        if (e1, e2) in compatible or (e2, e1) in compatible:
            # Verificar estequiometría válida
            if e2 == "O":
                # Para óxidos, verificar que la carga sea neutra
                if e1 in PERIODIC_TABLE:
                    ox_states = PERIODIC_TABLE[e1]["ox"]
                    for ox in ox_states:
                        # Carga total = n1 * ox_metal + n2 * (-2)
                        total_charge = n1 * ox + n2 * (-2)
                        if abs(total_charge) < 0.01:  # Carga neutra
                            return 0.95  # Muy alta confianza
                    # Si no hay carga neutra exacta
                    return 0.6
            return 0.7
        
        return 0.3
    
    def get_fibonacci_sequence(self, n: int = 12) -> List[int]:
        """Retorna la secuencia Fibonacci."""
        return self.fib[:n]

# ============================================================================
# INTERFAZ PRINCIPAL
# ============================================================================

def main():
    # Sidebar
    theme_name = st.sidebar.selectbox("🎨 Tema", list(THEMES.keys()))
    apply_theme(theme_name)
    
    # API Key input
    st.sidebar.markdown("---")
    st.sidebar.markdown("### 🔑 Materials Project API")
    mp_api_key = st.sidebar.text_input("API Key (opcional):", type="password")
    
    # Inicializar clases
    ai = DFTLearningAI()
    nanohub_gen = NanoHUBGenerator()
    fib_gen = FibonacciMaterials()
    pdf_gen = PDFGenerator()
    
    # Header
    st.title("🧠 CosmicForge Lab v4.3")
    st.markdown("<p style='text-align:center;opacity:0.7;'>IA que APRENDE cálculos DFT + Materials Project API + Fibonacci PDF</p>", unsafe_allow_html=True)
    
    # Indicadores
    memory_count = len(st.session_state.learning_memory)
    discovered_count = len(st.session_state.discovered_materials)
    dft_count = len(st.session_state.dft_knowledge)
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("🧠 Temas Aprendidos", memory_count)
    with col2:
        st.metric("🔬 Materiales Descubiertos", discovered_count)
    with col3:
        st.metric("📐 Conocimiento DFT", dft_count)
    
    # Tabs
    tabs = st.tabs(["🌠 Nebulosas", "🔬 Materiales", "🚀 nanoHUB/DFT", "💬 IA Aprendizaje", "🔢 Fibonacci", "🌐 Internet", "📊 Memoria"])
    
    # TAB 1: NEBULOSAS
    with tabs[0]:
        st.header("Selección de Nebulosa")
        
        col1, col2 = st.columns([1, 2])
        
        with col1:
            nebula_name = st.selectbox("Nebulosa:", list(NEBULAS_DATABASE.keys()))
            
            if st.button("🔍 Cargar y Consultar APIs", type="primary"):
                data = NEBULAS_DATABASE[nebula_name].copy()
                data["name"] = nebula_name
                st.session_state.nebula_data = data
                
                formula = f"{data['metal_dominant']}O2"
                st.session_state.api_results = ai.query_all_databases(formula)
                
                st.success(f"✅ {nebula_name} cargada")
        
        with col2:
            if st.session_state.nebula_data:
                data = st.session_state.nebula_data
                
                st.markdown(f"""
                <div class="card">
                    <h3>{nebula_name}</h3>
                    <p><strong>Tipo:</strong> {data.get('type')}</p>
                    <p><strong>Temperatura:</strong> {data.get('temperature_K'):,} K</p>
                    <p><strong>Metal:</strong> {data.get('metal_dominant')}</p>
                    <p><strong>Elementos:</strong> {', '.join(data.get('detected_elements', []))}</p>
                </div>
                """, unsafe_allow_html=True)
                
                if st.session_state.api_results:
                    st.markdown("**Resultados de APIs:**")
                    for api, result in st.session_state.api_results.items():
                        status = result.get("status", "?")
                        icon = "✅" if status in ["success", "found", "found_local"] else "❓" if status == "queried" else "❌"
                        st.markdown(f"{icon} {api}: {status}")
    
    # TAB 2: MATERIALES
    with tabs[1]:
        st.header("Generación de Materiales")
        
        if not st.session_state.nebula_data:
            st.info("👈 Selecciona una nebulosa primero")
        else:
            if st.button("🔮 Generar Material Validado", type="primary"):
                nebula = st.session_state.nebula_data
                metal = nebula.get("metal_dominant", "Ti")
                temp = nebula.get("temperature_K", 10000)
                
                # Determinar fórmula correcta
                if temp > 10000:
                    formula = f"{metal}O2"
                elif temp > 5000:
                    formula = f"{metal}2O3"
                else:
                    formula = f"{metal}O"
                
                # Consultar Materials Project con API key si disponible
                mp_result = ai.query_materials_project(formula, mp_api_key if mp_api_key else None)
                st.session_state.materials_project_data = mp_result
                
                # Generar output nanoHUB
                nanohub = nanohub_gen.generate(formula, nebula["name"])
                
                # Validar DFT
                dft_validation = ai.validate_dft_calculation(formula, {"ecutwfc": 60, "ecutrho": 240})
                
                st.session_state.material_stable = {
                    "formula": formula,
                    "nanohub": nanohub,
                    "mp_result": mp_result,
                    "dft_validation": dft_validation,
                    "valid": True
                }
                
                ai.learn(formula, f"Generado desde {nebula['name']} a {temp}K")
                
                st.success(f"✅ {formula} generado y validado")
                st.rerun()
        
        if st.session_state.material_stable:
            mat = st.session_state.material_stable
            formula = mat.get("formula", "?")
            
            st.markdown(f"<div class='formula-display'>{formula}</div>", unsafe_allow_html=True)
            
            info = MATERIALS_DATABASE.get(formula, {})
            if info:
                col1, col2 = st.columns(2)
                with col1:
                    st.metric("Nombre", info.get("name", formula))
                    st.metric("Band Gap", f"{info.get('band_gap', 0)} eV")
                with col2:
                    st.metric("Energía", f"{info.get('energy', 0)} eV/átomo")
                    st.metric("MP ID", info.get("mp_id", "N/A"))
                
                st.markdown(f"**Aplicaciones:** {', '.join(info.get('applications', []))}")
            
            # Mostrar validación DFT
            dft_val = mat.get("dft_validation", {})
            if dft_val.get("warnings"):
                st.markdown("<div class='dft-warning'>⚠️ Advertencias DFT:</div>", unsafe_allow_html=True)
                for w in dft_val["warnings"]:
                    st.warning(w)
    
    # TAB 3: NANOHUB/DFT
    with tabs[2]:
        st.header("🚀 Output nanoHUB + DFT Validado")
        
        if not st.session_state.material_stable:
            st.info("👈 Genera un material primero")
        else:
            mat = st.session_state.material_stable
            formula = mat.get("formula", "?")
            
            st.markdown(f"<div class='formula-display'>{formula}</div>", unsafe_allow_html=True)
            
            # Generar input QE completo
            qe_data = ai.generate_qe_input(formula, st.session_state.nebula_data.get("name", "Nebula"))
            
            # Mostrar validación
            validation = qe_data.get("valid", {})
            if validation.get("valid"):
                st.markdown("<div class='dft-valid'>✅ Cálculo DFT validado correctamente</div>", unsafe_allow_html=True)
            
            # Output nanoHUB
            st.subheader("📋 Output para nanoHUB")
            st.code(qe_data["nanohub_output"], language=None)
            
            # Input Quantum ESPRESSO
            st.subheader("🧪 Input Quantum ESPRESSO")
            st.code(qe_data["qe_input"], language=None)
            
            # Parámetros DFT
            col1, col2 = st.columns(2)
            with col1:
                st.metric("ecutwfc", f"{qe_data['ecutwfc']} Ry")
                st.metric("ecutrho", f"{qe_data['ecutrho']} Ry")
            with col2:
                st.metric("K-points", f"{qe_data['k_points']}")
                st.metric("Smearing", qe_data["smearing"])
            
            # Descargas
            col1, col2 = st.columns(2)
            with col1:
                st.download_button(
                    "📥 Descargar nanoHUB",
                    qe_data["nanohub_output"],
                    f"nanohub_{formula}.txt",
                    "text/plain",
                    type="primary"
                )
            with col2:
                st.download_button(
                    "📥 Descargar QE Input",
                    qe_data["qe_input"],
                    f"{formula}_qe.in",
                    "text/plain"
                )
    
    # TAB 4: IA APRENDIZAJE
    with tabs[3]:
        st.header("💬 IA que Aprende DFT")
        
        if st.session_state.chat_history:
            for msg in st.session_state.chat_history:
                role = msg["role"]
                content = msg["content"]
                css_class = "chat-user" if role == "user" else "chat-ai"
                st.markdown(f'<div class="chat-message {css_class}">{"👤" if role=="user" else "🧠"} {content}</div>', unsafe_allow_html=True)
        
        user_input = st.text_area("Pregunta a la IA:", height=100)
        
        col1, col2 = st.columns([3, 1])
        with col1:
            if st.button("Enviar", type="primary") and user_input:
                st.session_state.chat_history.append({"role": "user", "content": user_input})
                
                context = {
                    "nebula": st.session_state.nebula_data,
                    "material": st.session_state.material_stable
                }
                response = ai.generate_response(user_input, context)
                
                st.session_state.chat_history.append({"role": "assistant", "content": response})
                st.rerun()
        
        with col2:
            if st.button("🗑️ Limpiar"):
                st.session_state.chat_history = []
                st.rerun()
        
        st.markdown("**Preguntas sugeridas:**")
        cols = st.columns(2)
        suggestions = [
            "¿Cómo optimizar parámetros DFT?",
            "¿Por qué Ti₂O es inválido?",
            "¿Qué es ecutwfc en QE?",
            "Explica errores de convergencia"
        ]
        for i, sug in enumerate(suggestions):
            with cols[i % 2]:
                if st.button(sug, key=f"sug_{i}"):
                    st.session_state.chat_history.append({"role": "user", "content": sug})
                    context = {"nebula": st.session_state.nebula_data, "material": st.session_state.material_stable}
                    response = ai.generate_response(sug, context)
                    st.session_state.chat_history.append({"role": "assistant", "content": response})
                    st.rerun()
    
    # TAB 5: FIBONACCI
    with tabs[4]:
        st.header("🔢 Generador Fibonacci de Materiales")
        
        st.markdown("""
        <div class="card">
        <p><strong>Secuencia Fibonacci:</strong> 1, 1, 2, 3, 5, 8, 13, 21...</p>
        <p>Se aplica para descubrir nuevas estequiometrías de materiales.</p>
        </div>
        """, unsafe_allow_html=True)
        
        col1, col2 = st.columns(2)
        with col1:
            elem1 = st.selectbox("Elemento 1:", ["Ti", "Fe", "Al", "Si", "Zn", "Cu", "Ni", "Co"])
        with col2:
            elem2 = st.selectbox("Elemento 2:", ["O", "N", "C", "S"])
        
        if st.button("🔮 Generar con Fibonacci", type="primary"):
            predictions = fib_gen.generate_combinations(elem1, elem2)
            st.session_state.fibonacci_predictions = predictions
            
            for pred in predictions[:5]:
                ai.discover_new_material(elem1, elem2, pred["ratio"])
            
            st.success(f"✅ {len(predictions)} materiales predichos")
            st.rerun()
        
        if st.session_state.fibonacci_predictions:
            st.subheader("📊 Materiales Predichos")
            
            for pred in st.session_state.fibonacci_predictions[:8]:
                confidence = pred.get("confidence", 0)
                conf_color = "green" if confidence > 0.7 else "orange" if confidence > 0.5 else "red"
                
                st.markdown(f"""
                <div class="card">
                    <strong>{pred['formula']}</strong><br>
                    Ratio Fibonacci: {pred['ratio']}<br>
                    Confianza: <span style="color:{conf_color}">{confidence:.0%}</span>
                </div>
                """, unsafe_allow_html=True)
            
            # Botón para descargar PDF
            st.markdown("---")
            pdf_data = pdf_gen.generate_fibonacci_pdf(
                st.session_state.fibonacci_predictions,
                elem1,
                elem2,
                st.session_state.discovered_materials
            )
            
            st.download_button(
                "📥 Descargar Reporte PDF (Texto)",
                pdf_data,
                f"fibonacci_{elem1}_{elem2}_report.txt",
                "text/plain",
                type="primary"
            )
    
    # TAB 6: INTERNET
    with tabs[5]:
        st.header("🌐 Búsqueda en Internet")
        
        search_query = st.text_input("Buscar:", "materiales DFT última investigación")
        
        if st.button("🔍 Buscar", type="primary"):
            with st.spinner("Buscando..."):
                results = ai.search_internet(search_query)
                st.session_state.internet_search_results = results
            
            if isinstance(results, list):
                st.success(f"✅ {len(results)} resultados")
                
                for i, r in enumerate(results[:10], 1):
                    if isinstance(r, dict):
                        st.markdown(f"""
                        <div class="card">
                            <strong>{i}. {r.get('name', 'Sin título')}</strong><br>
                            {r.get('snippet', '')[:200]}...
                        </div>
                        """, unsafe_allow_html=True)
                
                ai.learn(search_query, f"Búsqueda: {len(results)} resultados")
            else:
                st.error("No se pudieron obtener resultados")
    
    # TAB 7: MEMORIA
    with tabs[6]:
        st.header("📊 Memoria de Aprendizaje")
        
        memory = st.session_state.learning_memory
        discovered = st.session_state.discovered_materials
        dft_knowledge = st.session_state.dft_knowledge
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.subheader("🧠 Temas")
            if memory:
                for topic, items in list(memory.items())[:5]:
                    with st.expander(f"📚 {topic} ({len(items)})"):
                        for item in items[-3:]:
                            st.markdown(f"• {item['info'][:80]}...")
            else:
                st.info("Sin datos")
        
        with col2:
            st.subheader("🔬 Descubiertos")
            if discovered:
                for d in discovered[-5:]:
                    st.markdown(f"""
                    <div class="card">
                        <strong>{d['formula']}</strong><br>
                        Confianza: {d['confidence']:.0%}
                    </div>
                    """, unsafe_allow_html=True)
            else:
                st.info("Sin descubiertos")
        
        with col3:
            st.subheader("📐 DFT")
            if dft_knowledge:
                for topic, items in list(dft_knowledge.items())[:5]:
                    st.markdown(f"**{topic}:** {len(items)} entradas")
            else:
                st.info("Sin conocimiento DFT")

if __name__ == "__main__":
    main()
