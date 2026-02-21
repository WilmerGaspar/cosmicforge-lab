"""
CosmicForge Lab - Main Application
Diseño de Materiales Inspirado en Firmas Astrofísicas
Versión 2.2 - Soporte para JSON de Anomaly Detector
"""

import streamlit as st
import json
import sys
import os
from pathlib import Path

# Add modules directory to path
MODULES_PATH = Path(__file__).parent / "modules"
sys.path.insert(0, str(MODULES_PATH))

from modules.physics_calculator import PhysicsCalculator
from modules.chemistry_engine import ChemistryEngine
from modules.file_generators import FileGenerator
from modules.visualizer import Visualizer
from modules.file_parser import FileParser

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
# CUSTOM CSS STYLES
# ============================================================================

st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        background: linear-gradient(90deg, #1E88E5, #7C4DFF);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        text-align: center;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.1rem;
        color: #666;
        text-align: center;
        margin-bottom: 1.5rem;
    }
    .metric-card {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1E88E5;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    .success-box {
        background-color: #d4edda;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #28a745;
        margin-top: 1rem;
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# PARSER PARA ANOMALY DETECTOR JSON
# ============================================================================

def parse_anomaly_detector_json(raw_data: dict) -> dict:
    """
    Convierte JSON de Anomaly Detector al formato que usa CosmicForge
    
    Mapeo de campos:
    - fractal_dimension → plugin_results.fractal_base.d0
    - criticality_score → plugin_results.renormalization_group.criticality_score
    - entropy → plugin_results.entropy.normalized_entropy
    - anisotropy → plugin_results.anisotropy.anisotropy_index
    - turbulence_beta → plugin_results.kolmogorov_1941.beta
    - lyapunov_max → plugin_results.lyapunov_stability.max_lyapunov
    """
    
    try:
        # Obtener plugin_results
        plugin_results = raw_data.get('plugin_results', {})
        
        # Extraer valores con valores por defecto
        fractal_base = plugin_results.get('fractal_base', {})
        kolmogorov = plugin_results.get('kolmogorov_1941', {})
        lyapunov = plugin_results.get('lyapunov_stability', {})
        renormalization = plugin_results.get('renormalization_group', {})
        anisotropy_data = plugin_results.get('anisotropy', {})
        entropy_data = plugin_results.get('entropy', {})
        
        # Construir datos en formato CosmicForge
        astro_data = {
            'object_name': raw_data.get('filename', 'Unknown').replace('.jpg', '').replace('.png', ''),
            'mode': raw_data.get('mode', 'balanced'),
            'global_score': raw_data.get('global_score', 0.5),
            
            # Campos principales para CosmicForge
            'fractal_dimension': fractal_base.get('d0', fractal_base.get('dimension_box', 1.5)),
            'criticality_score': renormalization.get('criticality_score', 0.5),
            'entropy': entropy_data.get('normalized_entropy', entropy_data.get('shannon_entropy_bits', 0.01) / 8),
            'anisotropy': anisotropy_data.get('anisotropy_index', 0.3),
            'turbulence_beta': kolmogorov.get('beta', 2.0),
            'lyapunov_max': lyapunov.get('max_lyapunov', -0.1),
            
            # Campos adicionales de Anomaly Detector
            'monte_carlo': raw_data.get('monte_carlo', {}),
            'fractal_base': fractal_base,
            'kolmogorov_1941': kolmogorov,
            'lyapunov_stability': lyapunov,
            'persistent_homology': plugin_results.get('persistent_homology', {}),
            'renormalization_group': renormalization,
            'anisotropy_details': anisotropy_data,
            'entropy_details': entropy_data,
            'report': raw_data.get('report', ''),
            
            # Indicar fuente
            'source': 'Anomaly Detector'
        }
        
        return astro_data
        
    except Exception as e:
        st.error(f"Error parseando JSON de Anomaly Detector: {e}")
        return None


def validate_anomaly_detector_json(data: dict) -> tuple:
    """Valida si el JSON tiene la estructura de Anomaly Detector"""
    
    # Verificar si tiene la estructura de Anomaly Detector
    if 'plugin_results' in data:
        return True, "Anomaly Detector JSON detectado"
    
    # Verificar si tiene campos planos (formato simple)
    required_fields = ['fractal_dimension', 'criticality_score', 'entropy', 'anisotropy', 'turbulence_beta']
    has_flat = all(field in data for field in required_fields)
    
    if has_flat:
        return True, "CosmicForge JSON simple detectado"
    
    # No tiene formato reconocido
    return False, "Formato JSON no reconocido. Se requiere formato Anomaly Detector o CosmicForge"


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def load_example_data():
    """Load example astrophysical data"""
    return {
        "object_name": "Orion Nebula M42",
        "fractal_dimension": 1.654,
        "criticality_score": 0.722,
        "entropy": 0.019,
        "anisotropy": 0.329,
        "turbulence_beta": 2.278,
        "lyapunov_max": -0.227,
        "mode": "balanced"
    }


def validate_astrophysical_data(data: dict):
    """Validate astrophysical data structure"""
    errors = []
    required_fields = ['fractal_dimension', 'criticality_score', 'entropy', 'anisotropy', 'turbulence_beta']
    
    for field in required_fields:
        if field not in data:
            errors.append(f"Campo requerido faltante: {field}")
        elif not isinstance(data[field], (int, float)):
            errors.append(f"Campo {field} debe ser numérico")
    
    # Validar rangos
    if 'fractal_dimension' in data:
        if not 0 <= data['fractal_dimension'] <= 3:
            errors.append("fractal_dimension debe estar entre 0 y 3")
    
    if 'criticality_score' in data:
        if not 0 <= data['criticality_score'] <= 1:
            errors.append("criticality_score debe estar entre 0 y 1")
    
    is_valid = len(errors) == 0
    return is_valid, errors


def display_metrics_from_data(data: dict):
    """Display astrophysical signature metrics"""
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Dimensión Fractal", f"{data.get('fractal_dimension', 0):.3f}")
        st.metric("Criticalidad", f"{data.get('criticality_score', 0):.3f}")
    
    with col2:
        st.metric("Entropía", f"{data.get('entropy', 0):.6f}")
        st.metric("Anisotropía", f"{data.get('anisotropy', 0):.3f}")
    
    with col3:
        st.metric("Turbulencia β", f"{data.get('turbulence_beta', 0):.3f}")
        st.metric("Lyapunov max", f"{data.get('lyapunov_max', 0):.4f}")
    
    with col4:
        st.metric("Objeto", str(data.get('object_name', 'Unknown'))[:20])
        st.metric("Fuente", data.get('source', 'Manual'))


def display_anomaly_detector_details(raw_data: dict):
    """Muestra detalles completos del JSON de Anomaly Detector"""
    
    if 'plugin_results' not in raw_data:
        return
    
    plugin = raw_data['plugin_results']
    
    st.subheader("🔬 Detalles de Anomaly Detector")
    
    # Crear tabs para cada plugin
    tabs = st.tabs(["🌌 Fractal", "🌀 Kolmogorov", "⚡ Lyapunov", "🔗 Topología", "📊 Entropía"])
    
    with tabs[0]:
        fractal = plugin.get('fractal_base', {})
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Dimensión D₀", f"{fractal.get('d0', 0):.4f}")
            st.metric("Dimensión D₁", f"{fractal.get('d1', 0):.4f}")
            st.metric("Dimensión D₂", f"{fractal.get('d2', 0):.4f}")
        with col2:
            st.metric("Lacunaridad", f"{fractal.get('lacunarity', 0):.4f}")
            st.metric("Índice Multifractal", f"{fractal.get('multifractality_index', 0):.4f}")
            st.metric("Complejidad", f"{fractal.get('complexity_score', 0):.4f}")
    
    with tabs[1]:
        kolmo = plugin.get('kolmogorov_1941', {})
        col1, col2 = st.columns(2)
        with col1:
            st.metric("β (turbulencia)", f"{kolmo.get('beta', 0):.4f}")
            st.metric("R²", f"{kolmo.get('r_squared', 0):.4f}")
        with col2:
            st.metric("Intermitencia", f"{kolmo.get('intermittency_factor', 0):.4f}")
            st.metric("Intensidad", f"{kolmo.get('turbulence_intensity', 0):.4f}")
    
    with tabs[2]:
        lyap = plugin.get('lyapunov_stability', {})
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Lyapunov Max", f"{lyap.get('max_lyapunov', 0):.6f}")
            st.metric("¿Caótico?", "Sí" if lyap.get('is_chaotic', False) else "No")
        with col2:
            st.metric("Entropía KS", f"{lyap.get('ks_entropy', 0):.4f}")
            st.metric("Dim. Kaplan-Yorke", f"{lyap.get('kaplan_yorke_dim', 0):.4f}")
        st.info(f"**Tipo de dinámica:** {lyap.get('dynamics_type', 'N/A')}")
    
    with tabs[3]:
        topo = plugin.get('persistent_homology', {})
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Betti₀", topo.get('betti_0', 0))
            st.metric("Betti₁", topo.get('betti_1', 0))
        with col2:
            st.metric("Característica Euler", topo.get('euler_characteristic', 0))
            st.metric("Entropía Topológica", f"{topo.get('topological_entropy', 0):.4f}")
        st.info(f"**Tipo:** {topo.get('topology_type', 'N/A')}")
    
    with tabs[4]:
        entropy = plugin.get('entropy', {})
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Shannon (bits)", f"{entropy.get('shannon_entropy_bits', 0):.6f}")
            st.metric("Entropía Normalizada", f"{entropy.get('normalized_entropy', 0):.8f}")
        with col2:
            st.metric("Entropía Máxima", entropy.get('max_possible_entropy', 8))
            st.metric("Densidad Info", f"{entropy.get('information_density', 0):.2e}")
        st.info(f"**Interpretación:** {entropy.get('interpretation', 'N/A')}")


# ============================================================================
# MAIN APPLICATION
# ============================================================================

# Header
st.markdown('<p class="main-header">🔬 CosmicForge Lab</p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Diseño de Materiales Inspirado en Firmas Astrofísicas</p>', unsafe_allow_html=True)

# =========================================================================
# SIDEBAR - DATA INPUT
# =========================================================================

st.sidebar.header("📥 Importar Datos")

input_method = st.sidebar.radio(
    "Método de entrada:",
    ["JSON de Anomaly Detector", "Datos de ejemplo", "Editor manual",
     "Archivo QE/CIF/LAMMPS"],
    help="Selecciona el tipo de entrada"
)

new_astro_data = None

if input_method == "JSON de Anomaly Detector":
    uploaded_file = st.sidebar.file_uploader(
        "Archivo JSON de Anomaly Detector",
        type=['json'],
        help="Archivo JSON generado por Anomaly Detector"
    )
    
    if uploaded_file is not None:
        try:
            # Leer archivo
            string_data = uploaded_file.read().decode('utf-8')
            raw_json = json.loads(string_data)
            
            # Guardar JSON original
            st.session_state.raw_json = raw_json
            
            # Detectar formato y convertir
            is_valid, format_msg = validate_anomaly_detector_json(raw_json)
            
            if not is_valid:
                st.sidebar.error(f"❌ {format_msg}")
            else:
                st.sidebar.success(f"✅ {format_msg}")
                
                # Si es Anomaly Detector, convertir formato
                if 'plugin_results' in raw_json:
                    new_astro_data = parse_anomaly_detector_json(raw_json)
                else:
                    new_astro_data = raw_json
                
                if new_astro_data:
                    # Validar
                    is_valid2, errors = validate_astrophysical_data(new_astro_data)
                    if is_valid2:
                        st.session_state.astro_data = new_astro_data
                        st.sidebar.success(f"✅ Cargado: {new_astro_data.get('object_name', 'Unknown')}")
                    else:
                        st.sidebar.error("❌ Error en datos:")
                        for e in errors:
                            st.sidebar.write(f"• {e}")
                            
        except json.JSONDecodeError as e:
            st.sidebar.error(f"❌ Error JSON: {e}")
        except Exception as e:
            st.sidebar.error(f"❌ Error: {e}")

elif input_method == "Datos de ejemplo":
    example_choice = st.sidebar.selectbox(
        "Ejemplo:",
        ["Orion Nebula M42", "Crab Nebula", "Ring Nebula M57", "Triangulum Galaxy (demo)"]
    )
    
    if st.sidebar.button("📋 Cargar ejemplo", type="primary"):
        if example_choice == "Orion Nebula M42":
            new_astro_data = {
                "object_name": "Orion Nebula M42",
                "fractal_dimension": 1.654,
                "criticality_score": 0.722,
                "entropy": 0.019,
                "anisotropy": 0.329,
                "turbulence_beta": 2.278,
                "lyapunov_max": -0.227,
                "mode": "balanced"
            }
        elif example_choice == "Crab Nebula":
            new_astro_data = {
                "object_name": "Crab Nebula",
                "fractal_dimension": 1.82,
                "criticality_score": 0.89,
                "entropy": 0.032,
                "anisotropy": 0.456,
                "turbulence_beta": 2.45,
                "lyapunov_max": -0.15,
                "mode": "turbulent"
            }
        elif example_choice == "Ring Nebula M57":
            new_astro_data = {
                "object_name": "Ring Nebula M57",
                "fractal_dimension": 1.45,
                "criticality_score": 0.55,
                "entropy": 0.012,
                "anisotropy": 0.28,
                "turbulence_beta": 1.95,
                "lyapunov_max": -0.35,
                "mode": "stable"
            }
        else:  # Triangulum Galaxy (demo con datos del usuario)
            new_astro_data = {
                "object_name": "Triangulum Galaxy",
                "fractal_dimension": 1.932,
                "criticality_score": 0.781,
                "entropy": 0.00000964,
                "anisotropy": 0.233,
                "turbulence_beta": 1.839,
                "lyapunov_max": 0.109,
                "mode": "balanced"
            }
        
        st.session_state.astro_data = new_astro_data
        st.sidebar.success(f"✅ {example_choice} cargado")

elif input_method == "Editor manual":
    st.sidebar.subheader("Parámetros principales")
    
    object_name = st.sidebar.text_input("Nombre del objeto", "Custom Object")
    fractal_dim = st.sidebar.slider("Dimensión Fractal", 0.0, 3.0, 1.5, 0.001)
    criticality = st.sidebar.slider("Criticalidad", 0.0, 1.0, 0.5, 0.01)
    entropy = st.sidebar.number_input("Entropía", 0.0, 1.0, 0.01, 0.0001, format="%.6f")
    anisotropy = st.sidebar.slider("Anisotropía", 0.0, 1.0, 0.3, 0.01)
    turbulence = st.sidebar.slider("Turbulencia β", 0.0, 5.0, 2.0, 0.01)
    lyapunov = st.sidebar.number_input("Lyapunov max", -1.0, 1.0, -0.1, 0.001)
    
    # Campos avanzados de Anomaly Detector
    with st.sidebar.expander("🔧 Campos avanzados (Anomaly Detector)"):
        monte_p = st.number_input("Monte Carlo p-value", 0.0, 1.0, 0.05, 0.001)
        is_significant = st.checkbox("Significativo", False)
        lyap_chaotic = st.checkbox("Es caótico", True)
        dynamics_type = st.selectbox("Tipo de dinámica", ["Caótico fuerte", "Caótico débil", "Estable"])
        betti_0 = st.number_input("Betti₀", 0, 100, 12)
        betti_1 = st.number_input("Betti₁", 0, 100, 4)
        lacunarity = st.number_input("Lacunaridad", 0.0, 2.0, 0.67, 0.01)
        intermittency = st.number_input("Intermitencia", 0.0, 1.0, 0.35, 0.01)
    
    if st.sidebar.button("✅ Aplicar parámetros", type="primary"):
        manual_data = {
            "object_name": object_name,
            "fractal_dimension": fractal_dim,
            "criticality_score": criticality,
            "entropy": entropy,
            "anisotropy": anisotropy,
            "turbulence_beta": turbulence,
            "lyapunov_max": lyapunov,
            "mode": "custom",
            "source": "Editor manual",
            # Campos avanzados
            "monte_carlo": {"p_value": monte_p, "is_significant": str(is_significant)},
            "lyapunov_stability": {"is_chaotic": lyap_chaotic, "dynamics_type": dynamics_type},
            "persistent_homology": {"betti_0": betti_0, "betti_1": betti_1},
            "fractal_base": {"lacunarity": lacunarity},
            "kolmogorov_1941": {"intermittency_factor": intermittency}
        }
        st.session_state.astro_data = manual_data
        st.sidebar.success("✅ Parámetros aplicados")

elif input_method == "Archivo QE/CIF/LAMMPS":
    file_type = st.sidebar.selectbox("Tipo de archivo", ["Quantum ESPRESSO", "CIF", "LAMMPS"])
    ext_map = {"Quantum ESPRESSO": ['in', 'pw'], "CIF": ['cif'], "LAMMPS": ['data', 'lmp']}
    
    uploaded_file = st.sidebar.file_uploader("Archivo", type=ext_map[file_type])
    
    if uploaded_file is not None:
        try:
            content = uploaded_file.read().decode('utf-8')
            parser = FileParser()
            
            if file_type == "Quantum ESPRESSO":
                parsed = parser.parse_quantum_espresso(content)
            elif file_type == "CIF":
                parsed = parser.parse_cif(content)
            else:
                parsed = parser.parse_lammps_data(content)
            
            st.session_state.astro_data = parsed
            st.sidebar.success(f"✅ {file_type} parseado")
        except Exception as e:
            st.sidebar.error(f"❌ Error: {e}")

# Botón limpiar
if st.session_state.astro_data is not None:
    if st.sidebar.button("🗑️ Limpiar datos"):
        st.session_state.astro_data = None
        st.session_state.physical_props = None
        st.session_state.recipe = None
        st.session_state.raw_json = None
        st.rerun()

# =========================================================================
# SIDEBAR - MATERIAL SELECTION
# =========================================================================

st.sidebar.header("🎯 Material")

material_type = st.sidebar.selectbox(
    "Tipo:",
    ["Óxido metálico", "Metal puro", "Cerámica", "Polímero", "Composite", "Nanopartícula"]
)

metal = st.sidebar.selectbox(
    "Metal:",
    ["Ti", "Al", "Fe", "Zn", "Cu", "Ni", "Co", "Mn", "Ag", "Au", "Pt", "Pd"]
)

with st.sidebar.expander("⚙️ Avanzado"):
    bulk_density = st.number_input("Densidad bulk (g/cm³)", value=3.0, min_value=0.1, max_value=20.0)
    bulk_conductivity = st.number_input("Conductividad (W/m·K)", value=50.0, min_value=0.1, max_value=500.0)
    bulk_elastic = st.number_input("Módulo elástico (GPa)", value=200.0, min_value=1.0, max_value=1000.0)
    supercell_size = st.slider("Supercelda", 2, 10, 5)

# =========================================================================
# MAIN CONTENT
# =========================================================================

if st.session_state.astro_data is not None:
    astro_data = st.session_state.astro_data
    
    st.header("🌌 Firma Astrofísica Detectada")
    display_metrics_from_data(astro_data)
    
    # Mostrar detalles de Anomaly Detector si existen
    if st.session_state.raw_json and 'plugin_results' in st.session_state.raw_json:
        with st.expander("📊 Ver detalles completos de Anomaly Detector", expanded=False):
            display_anomaly_detector_details(st.session_state.raw_json)
    
    if st.button("🚀 Generar Receta de Material", type="primary", use_container_width=True):
        with st.spinner("Calculando propiedades..."):
            try:
                physics_calc = PhysicsCalculator(
                    bulk_density=bulk_density,
                    bulk_conductivity=bulk_conductivity,
                    bulk_elastic=bulk_elastic
                )
                chem_engine = ChemistryEngine()
                file_gen = FileGenerator(supercell_size=supercell_size)
                
                metal_symbol = metal.split(' ')[0] if ' ' in metal else metal
                material_type_clean = material_type.split()[0].lower()
                
                physical_props = physics_calc.calculate_all_properties(astro_data)
                st.session_state.physical_props = physical_props
                
                recipe = chem_engine.generate_complete_recipe(
                    metal=metal_symbol,
                    material_type=material_type_clean,
                    astrophysical_data=astro_data
                )
                st.session_state.recipe = recipe
                
                st.success("✅ Cálculos completados!")
                st.rerun()
                
            except Exception as e:
                st.error(f"❌ Error: {e}")
                st.exception(e)

# Mostrar resultados
if st.session_state.physical_props is not None and st.session_state.recipe is not None:
    physical_props = st.session_state.physical_props
    recipe = st.session_state.recipe
    astro_data = st.session_state.astro_data
    
    viz = Visualizer()
    file_gen = FileGenerator()
    
    tab1, tab2, tab3, tab4 = st.tabs(["📊 Propiedades", "⚗️ Síntesis", "📈 Visualización", "📁 Archivos"])
    
    with tab1:
        st.subheader("Propiedades del Material")
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Porosidad", f"{physical_props.get('porosity', 0):.1%}")
            st.metric("Densidad", f"{physical_props.get('density', 0):.2f} g/cm³")
        
        with col2:
            st.metric("Cond. Térmica", f"{physical_props.get('thermal_conductivity', 0):.1f} W/m·K")
            st.metric("Módulo Elástico", f"{physical_props.get('elastic_modulus', 0):.1f} GPa")
        
        with col3:
            st.metric("Área Superficial", f"{physical_props.get('surface_area', 0):.0f} m²/g")
            st.metric("Band Gap", f"{physical_props.get('band_gap', 0):.2f} eV")
        
        with col4:
            st.metric("Calidad", f"{physical_props.get('quality_score', 0):.0%}")
            st.metric("E. Activación", f"{physical_props.get('activation_energy', 0):.3f} eV")
    
    with tab2:
        st.subheader("Receta de Síntesis")
        conditions = recipe.get('reaction_conditions', {})
        
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("### 🌡️ Condiciones")
            st.write(f"**Temperatura:** {conditions.get('temperature_C', 0):.0f}°C")
            st.write(f"**Tiempo:** {conditions.get('reaction_time_hours', 0):.1f} horas")
            st.write(f"**pH:** {conditions.get('pH', 7):.1f}")
        
        with col2:
            st.markdown("### 🧪 Precursores")
            for precursor in recipe.get('precursors', []):
                st.write(f"• {precursor}")
        
        st.markdown("### 📝 Procedimiento")
        for i, step in enumerate(recipe.get('step_by_step', []), 1):
            st.write(f"**{i}.** {step}")
    
    with tab3:
        st.subheader("Visualización")
        try:
            col1, col2 = st.columns(2)
            with col1:
                img1 = viz.visualize_crystal_lattice(
                    physical_props.get('porosity', 0.3),
                    physical_props.get('density', 2.0),
                    astro_data.get('object_name', 'Material')
                )
                st.image(img1, caption="Estructura Cristalina", use_column_width=True)
            with col2:
                img2 = viz.visualize_properties_summary(physical_props)
                st.image(img2, caption="Resumen", use_column_width=True)
        except Exception as e:
            st.warning(f"Error: {e}")
    
    with tab4:
        st.subheader("Archivos para Simuladores")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            poscar = file_gen.generate_poscar(physical_props, recipe, 'oxide')
            st.download_button("⬇️ POSCAR", poscar, "POSCAR.vasp", "text/plain")
        
        with col2:
            lammps = file_gen.generate_lammps_input(physical_props, recipe, 'oxide')
            st.download_button("⬇️ LAMMPS", lammps, "input.lmp", "text/plain")
        
        with col3:
            csv_data = file_gen.generate_properties_csv(physical_props, recipe)
            st.download_button("⬇️ CSV", csv_data, "properties.csv", "text/csv")

else:
    st.info("👈 Carga un archivo JSON de Anomaly Detector para comenzar")
    
    with st.expander("📖 Ver formato JSON esperado"):
        st.markdown("""
        **Formato 1: Anomaly Detector (recomendado)**
        - Archivo JSON con `plugin_results` que contiene:
          - `fractal_base` → dimensión fractal
          - `kolmogorov_1941` → turbulencia β
          - `lyapunov_stability` → Lyapunov max
          - `renormalization_group` → criticalidad
          - `anisotropy` → anisotropía
          - `entropy` → entropía
        
        **Formato 2: CosmicForge simple**
        """)
        st.code(json.dumps(load_example_data(), indent=2), language="json")

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666;">
    <p><strong>CosmicForge Lab v2.2</strong> | Soporte para Anomaly Detector JSON</p>
</div>
""", unsafe_allow_html=True)
