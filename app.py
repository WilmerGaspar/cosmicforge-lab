"""
CosmicForge Lab - Main Application
Diseño de Materiales Inspirado en Firmas Astrofísicas
Versión 2.1 - Bugs corregidos + mejoras
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
# SESSION STATE - Para mantener datos entre re-renderizados
# ============================================================================

if 'astro_data' not in st.session_state:
    st.session_state.astro_data = None
if 'physical_props' not in st.session_state:
    st.session_state.physical_props = None
if 'recipe' not in st.session_state:
    st.session_state.recipe = None

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
    .info-box {
        background-color: #e7f3ff;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #0d6efd;
    }
    .warning-box {
        background-color: #fff3cd;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #ffc107;
    }
    .stMetric label {
        font-size: 0.85rem;
    }
</style>
""", unsafe_allow_html=True)

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

def load_example_data_2():
    """Second example - Crab Nebula"""
    return {
        "object_name": "Crab Nebula",
        "fractal_dimension": 1.82,
        "criticality_score": 0.89,
        "entropy": 0.032,
        "anisotropy": 0.456,
        "turbulence_beta": 2.45,
        "lyapunov_max": -0.15,
        "mode": "turbulent"
    }

def load_example_data_3():
    """Third example - Ring Nebula"""
    return {
        "object_name": "Ring Nebula M57",
        "fractal_dimension": 1.45,
        "criticality_score": 0.55,
        "entropy": 0.012,
        "anisotropy": 0.28,
        "turbulence_beta": 1.95,
        "lyapunov_max": -0.35,
        "mode": "stable"
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
        st.metric("Entropía", f"{data.get('entropy', 0):.4f}")
        st.metric("Anisotropía", f"{data.get('anisotropy', 0):.3f}")
    
    with col3:
        st.metric("Turbulencia β", f"{data.get('turbulence_beta', 0):.3f}")
        st.metric("Lyapunov max", f"{data.get('lyapunov_max', 0):.3f}")
    
    with col4:
        st.metric("Objeto", str(data.get('object_name', 'Unknown'))[:20])
        st.metric("Modo", data.get('mode', 'balanced'))

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
    ["Datos de ejemplo", "JSON de Anomaly Detector", "Editor manual", 
     "Archivo Quantum ESPRESSO", "Archivo CIF", "Archivo LAMMPS"],
    help="Selecciona el tipo de archivo a importar"
)

# Variable temporal para nuevos datos
new_astro_data = None

if input_method == "Datos de ejemplo":
    st.sidebar.subheader("Selecciona un ejemplo:")
    example_choice = st.sidebar.selectbox(
        "Objeto astrofísico:",
        ["Orion Nebula M42", "Crab Nebula", "Ring Nebula M57"]
    )
    
    if st.sidebar.button("📋 Cargar ejemplo", type="primary"):
        if example_choice == "Orion Nebula M42":
            new_astro_data = load_example_data()
        elif example_choice == "Crab Nebula":
            new_astro_data = load_example_data_2()
        else:
            new_astro_data = load_example_data_3()
        
        st.session_state.astro_data = new_astro_data
        st.sidebar.success(f"✅ {example_choice} cargado")

elif input_method == "JSON de Anomaly Detector":
    uploaded_file = st.sidebar.file_uploader(
        "Archivo JSON de firma astrofísica",
        type=['json'],
        help="Archivo generado por la app Anomaly Detector"
    )
    
    if uploaded_file is not None:
        try:
            # Leer el contenido del archivo
            stringio = uploaded_file.read()
            string_data = stringio.decode('utf-8')
            
            # Parsear JSON
            astro_data_parsed = json.loads(string_data)
            
            # Validar
            is_valid, errors = validate_astrophysical_data(astro_data_parsed)
            
            if not is_valid:
                st.sidebar.error("❌ Datos inválidos:")
                for error in errors:
                    st.sidebar.write(f"• {error}")
            else:
                st.session_state.astro_data = astro_data_parsed
                st.sidebar.success("✅ Archivo JSON cargado correctamente")
                st.sidebar.json(astro_data_parsed)
                
        except json.JSONDecodeError as e:
            st.sidebar.error(f"❌ Error JSON: {e}")
            st.sidebar.info("Asegúrate que el archivo sea JSON válido")
        except Exception as e:
            st.sidebar.error(f"❌ Error: {e}")

elif input_method == "Editor manual":
    st.sidebar.subheader("Editor de parámetros")
    
    object_name = st.sidebar.text_input("Nombre del objeto", "Custom Object")
    fractal_dim = st.sidebar.slider("Dimensión Fractal", 0.0, 3.0, 1.5, 0.001)
    criticality = st.sidebar.slider("Criticalidad", 0.0, 1.0, 0.5, 0.01)
    entropy = st.sidebar.slider("Entropía", 0.0, 1.0, 0.02, 0.001)
    anisotropy = st.sidebar.slider("Anisotropía", 0.0, 1.0, 0.3, 0.01)
    turbulence = st.sidebar.slider("Turbulencia β", 0.0, 5.0, 2.0, 0.01)
    lyapunov = st.sidebar.slider("Lyapunov max", -1.0, 1.0, -0.2, 0.01)
    
    if st.sidebar.button("✅ Aplicar parámetros", type="primary"):
        manual_data = {
            "object_name": object_name,
            "fractal_dimension": fractal_dim,
            "criticality_score": criticality,
            "entropy": entropy,
            "anisotropy": anisotropy,
            "turbulence_beta": turbulence,
            "lyapunov_max": lyapunov,
            "mode": "custom"
        }
        st.session_state.astro_data = manual_data
        st.sidebar.success("✅ Parámetros aplicados")

elif input_method == "Archivo Quantum ESPRESSO":
    uploaded_file = st.sidebar.file_uploader("Archivo QE (.in, .pw)", type=['in', 'pw', 'pwi'])
    if uploaded_file is not None:
        try:
            content = uploaded_file.read().decode('utf-8')
            parser = FileParser()
            parsed_data = parser.parse_quantum_espresso(content)
            st.session_state.astro_data = parsed_data
            st.sidebar.success("✅ Archivo QE parseado")
        except Exception as e:
            st.sidebar.error(f"❌ Error: {e}")

elif input_method == "Archivo CIF":
    uploaded_file = st.sidebar.file_uploader("Archivo CIF", type=['cif'])
    if uploaded_file is not None:
        try:
            content = uploaded_file.read().decode('utf-8')
            parser = FileParser()
            parsed_data = parser.parse_cif(content)
            st.session_state.astro_data = parsed_data
            st.sidebar.success("✅ Archivo CIF parseado")
        except Exception as e:
            st.sidebar.error(f"❌ Error: {e}")

elif input_method == "Archivo LAMMPS":
    uploaded_file = st.sidebar.file_uploader("Archivo LAMMPS", type=['data', 'lmp'])
    if uploaded_file is not None:
        try:
            content = uploaded_file.read().decode('utf-8')
            parser = FileParser()
            parsed_data = parser.parse_lammps_data(content)
            st.session_state.astro_data = parsed_data
            st.sidebar.success("✅ Archivo LAMMPS parseado")
        except Exception as e:
            st.sidebar.error(f"❌ Error: {e}")

# Botón para limpiar datos
if st.session_state.astro_data is not None:
    if st.sidebar.button("🗑️ Limpiar datos"):
        st.session_state.astro_data = None
        st.session_state.physical_props = None
        st.session_state.recipe = None
        st.rerun()

# =========================================================================
# SIDEBAR - MATERIAL SELECTION
# =========================================================================

st.sidebar.header("🎯 Área de Aplicación")
area = st.sidebar.selectbox(
    "Área:",
    ["Mecánica", "Química", "Industrial", "Energía", 
     "Aeroespacial", "Biomédico", "Electrónica"]
)

st.sidebar.header("⚗️ Material")

material_type = st.sidebar.selectbox(
    "Tipo:",
    ["Óxido metálico", "Metal puro", "Cerámica", 
     "Polímero", "Composite", "Nanopartícula"]
)

metal = st.sidebar.selectbox(
    "Metal:",
    ["Ti (Titanio)", "Al (Aluminio)", "Fe (Hierro)", 
     "Zn (Zinc)", "Cu (Cobre)", "Ni (Níquel)",
     "Co (Cobalto)", "Mn (Manganeso)", "Ag (Plata)",
     "Au (Oro)", "Pt (Platino)", "Pd (Paladio)"]
)

with st.sidebar.expander("⚙️ Avanzado"):
    bulk_density = st.number_input("Densidad bulk (g/cm³)", value=3.0, min_value=0.1, max_value=20.0)
    bulk_conductivity = st.number_input("Conductividad (W/m·K)", value=50.0, min_value=0.1, max_value=500.0)
    bulk_elastic = st.number_input("Módulo elástico (GPa)", value=200.0, min_value=1.0, max_value=1000.0)
    supercell_size = st.slider("Supercelda", 2, 10, 5)

# =========================================================================
# MAIN CONTENT
# =========================================================================

# Verificar si hay datos en session_state
if st.session_state.astro_data is not None:
    astro_data = st.session_state.astro_data
    
    st.header("🌌 Firma Astrofísica Detectada")
    display_metrics_from_data(astro_data)
    
    if st.button("🚀 Generar Receta de Material", type="primary", use_container_width=True):
        with st.spinner("Calculando propiedades..."):
            try:
                # Initialize calculators
                physics_calc = PhysicsCalculator(
                    bulk_density=bulk_density,
                    bulk_conductivity=bulk_conductivity,
                    bulk_elastic=bulk_elastic
                )
                chem_engine = ChemistryEngine()
                file_gen = FileGenerator(supercell_size=supercell_size)
                viz = Visualizer()
                
                # Extract metal symbol (fix: obtener solo el símbolo)
                metal_symbol = metal.split(' ')[0].strip()
                material_type_clean = material_type.split()[0].lower()
                
                # Calculate properties
                physical_props = physics_calc.calculate_all_properties(astro_data)
                st.session_state.physical_props = physical_props
                
                # Generate recipe
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

# Mostrar resultados si existen
if st.session_state.physical_props is not None and st.session_state.recipe is not None:
    physical_props = st.session_state.physical_props
    recipe = st.session_state.recipe
    astro_data = st.session_state.astro_data
    
    # Inicializar visualizador
    viz = Visualizer()
    file_gen = FileGenerator()
    
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📊 Propiedades", "⚗️ Síntesis", "📈 Visualización", "📁 Archivos", "📋 Informe"
    ])
    
    with tab1:
        st.subheader("Propiedades del Material")
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.markdown('<div class="metric-card">', unsafe_allow_html=True)
            st.metric("Porosidad", f"{physical_props.get('porosity', 0):.1%}")
            st.metric("Densidad", f"{physical_props.get('density', 0):.2f} g/cm³")
            st.markdown('</div>', unsafe_allow_html=True)
        
        with col2:
            st.markdown('<div class="metric-card">', unsafe_allow_html=True)
            st.metric("Cond. Térmica", f"{physical_props.get('thermal_conductivity', 0):.1f} W/m·K")
            st.metric("Módulo Elástico", f"{physical_props.get('elastic_modulus', 0):.1f} GPa")
            st.markdown('</div>', unsafe_allow_html=True)
        
        with col3:
            st.markdown('<div class="metric-card">', unsafe_allow_html=True)
            st.metric("Área Superficial", f"{physical_props.get('surface_area', 0):.0f} m²/g")
            st.metric("Band Gap", f"{physical_props.get('band_gap', 0):.2f} eV")
            st.markdown('</div>', unsafe_allow_html=True)
        
        with col4:
            st.markdown('<div class="metric-card">', unsafe_allow_html=True)
            st.metric("Calidad", f"{physical_props.get('quality_score', 0):.0%}")
            st.metric("E. Activación", f"{physical_props.get('activation_energy', 0):.3f} eV")
            st.markdown('</div>', unsafe_allow_html=True)
        
        # Ecuaciones
        with st.expander("📐 Ver ecuaciones utilizadas"):
            st.markdown("""
            | Propiedad | Ecuación | Descripción |
            |-----------|----------|-------------|
            | Porosidad | Φ = 1 - (D/3)^1.5 | D = dimensión fractal |
            | Densidad | ρ = ρ_bulk × (1-Φ) | Densidad efectiva |
            | Conductividad | k = k_bulk × (1-Φ)^1.5 | Térmica |
            | Módulo elástico | E = E_bulk × (1-Φ)² | Young's modulus |
            """)
    
    with tab2:
        st.subheader("Receta de Síntesis")
        conditions = recipe.get('reaction_conditions', {})
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("### 🌡️ Condiciones")
            st.write(f"**Temperatura:** {conditions.get('temperature_C', 0):.0f}°C ({conditions.get('temperature_K', 0):.0f} K)")
            st.write(f"**Tiempo:** {conditions.get('reaction_time_hours', 0):.1f} horas")
            st.write(f"**pH:** {conditions.get('pH', 7):.1f}")
            st.write(f"**Presión:** {conditions.get('pressure_atm', 1):.2f} atm")
            st.write(f"**Agitación:** {conditions.get('stirring_speed_rpm', 0):.0f} RPM")
        
        with col2:
            st.markdown("### 🧪 Precursores")
            for precursor in recipe.get('precursors', []):
                st.write(f"• {precursor}")
        
        st.markdown("### 📝 Procedimiento")
        for i, step in enumerate(recipe.get('step_by_step', []), 1):
            st.write(f"**{i}.** {step}")
        
        with st.expander("⚠️ Seguridad"):
            for note in recipe.get('safety_notes', []):
                st.write(f"• {note}")
    
    with tab3:
        st.subheader("Visualización")
        
        try:
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("**Estructura Cristalina**")
                img1 = viz.visualize_crystal_lattice(
                    physical_props.get('porosity', 0.3),
                    physical_props.get('density', 2.0),
                    astro_data.get('object_name', 'Material')
                )
                st.image(img1, use_column_width=True)
            
            with col2:
                st.markdown("**Resumen de Propiedades**")
                img2 = viz.visualize_properties_summary(physical_props)
                st.image(img2, use_column_width=True)
                
        except Exception as e:
            st.warning(f"Error en visualización: {e}")
    
    with tab4:
        st.subheader("Archivos para Simuladores")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.markdown("**VASP / QE**")
            poscar = file_gen.generate_poscar(physical_props, recipe, 'oxide')
            st.download_button("⬇️ POSCAR", poscar, "POSCAR.vasp", "text/plain")
        
        with col2:
            st.markdown("**LAMMPS**")
            lammps = file_gen.generate_lammps_input(physical_props, recipe, 'oxide')
            st.download_button("⬇️ LAMMPS", lammps, "input.lmp", "text/plain")
        
        with col3:
            st.markdown("**Datos**")
            csv_data = file_gen.generate_properties_csv(physical_props, recipe)
            st.download_button("⬇️ CSV", csv_data, "properties.csv", "text/csv")
        
        col4, col5 = st.columns(2)
        
        with col4:
            xyz_data = file_gen.generate_xyz(physical_props, recipe)
            st.download_button("⬇️ XYZ (visualización)", xyz_data, "structure.xyz", "text/plain")
        
        with col5:
            json_export = json.dumps({
                'astrophysical_data': astro_data,
                'physical_properties': physical_props,
                'recipe': recipe
            }, indent=2, default=str)
            st.download_button("⬇️ JSON completo", json_export, "cosmicforge_data.json", "application/json")
    
    with tab5:
        st.subheader("Informe Técnico")
        report = file_gen.generate_report_metadata(physical_props, recipe)
        
        st.markdown(f"""
        ---
        # COSMICFORGE LAB - INFORME TÉCNICO
        
        **Fecha:** {report.get('timestamp', 'N/A')}  
        **Objeto Astrofísico:** {report.get('object_name', 'N/A')}  
        **Material:** {report.get('metal_name', report.get('metal', 'N/A'))} {report.get('material_type', '')}
        
        ---
        
        ## PROPIEDADES PREDICHAS
        
        | Propiedad | Valor |
        |-----------|-------|
        | Porosidad | {report.get('porosity', 'N/A')} |
        | Densidad | {report.get('density', 'N/A')} |
        | Temperatura de síntesis | {report.get('temperature', 'N/A')} |
        | Tiempo de reacción | {report.get('reaction_time', 'N/A')} |
        
        ---
        
        ## VALIDACIÓN
        
        | Nivel | Estado |
        |-------|--------|
        | Matemático | ✅ Demostrado |
        | Computacional | 🟡 Requiere DFT/MD |
        | Experimental | 🔴 Pendiente |
        ---
        """)

else:
    # Estado inicial - sin datos
    st.info("👈 Carga datos desde el panel lateral para comenzar")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("""
        ### 📄 JSON
        Archivo de firma astrofísica
        """)
    
    with col2:
        st.markdown("""
        ### 🔬 Archivos
        Quantum ESPRESSO, CIF, LAMMPS
        """)
    
    with col3:
        st.markdown("""
        ### ✏️ Editor
        Crea parámetros manualmente
        """)
    
    with st.expander("📖 Formato JSON esperado"):
        st.code(json.dumps(load_example_data(), indent=2), language="json")

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666; padding: 1rem;">
    <p><strong>CosmicForge Lab v2.1</strong> | Diseño de materiales inspirado en astrofísica</p>
    <p>🔗 <a href="https://github.com/WilmerGaspar/cosmicforge-lab" target="_blank">GitHub</a></p>
</div>
""", unsafe_allow_html=True)
