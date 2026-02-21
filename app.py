"""
CosmicForge Lab - Main Application
Diseño de Materiales Inspirado en Firmas Astrofísicas
Versión 2.0 - Código mejorado y corregido
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
# CUSTOM CSS STYLES
# ============================================================================

st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        font-weight: bold;
        background: linear-gradient(90deg, #1E88E5, #7C4DFF);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        text-align: center;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        text-align: center;
        margin-bottom: 2rem;
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
        padding: 1.5rem;
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
        st.metric("Dimensión Fractal", f"{data.get('fractal_dimension', 0):.3f}")
        st.metric("Criticalidad", f"{data.get('criticality_score', 0):.3f}")
    
    with col2:
        st.metric("Entropía", f"{data.get('entropy', 0):.3f}")
        st.metric("Anisotropía", f"{data.get('anisotropy', 0):.3f}")
    
    with col3:
        st.metric("Turbulencia β", f"{data.get('turbulence_beta', 0):.3f}")
        st.metric("Lyapunov max", f"{data.get('lyapunov_max', 0):.3f}")
    
    with col4:
        st.metric("Objeto", data.get('object_name', 'Unknown'))
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
    ["JSON de Anomaly Detector", "Datos de ejemplo", "Archivo Quantum ESPRESSO", 
     "Archivo CIF", "Archivo LAMMPS"],
    help="Selecciona el tipo de archivo a importar"
)

astro_data = None

if input_method == "JSON de Anomaly Detector":
    uploaded_file = st.sidebar.file_uploader(
        "Archivo JSON de firma astrofísica",
        type=['json'],
        help="Archivo generado por la app Anomaly Detector"
    )
    
    if uploaded_file is not None:
        try:
            astro_data = json.load(uploaded_file)
            is_valid, errors = validate_astrophysical_data(astro_data)
            if not is_valid:
                st.sidebar.error("❌ Datos inválidos:")
                for error in errors:
                    st.sidebar.write(f"• {error}")
                astro_data = None
            else:
                st.sidebar.success("✅ Archivo cargado correctamente")
        except json.JSONDecodeError as e:
            st.sidebar.error(f"❌ Error parsing JSON: {e}")

elif input_method == "Datos de ejemplo":
    if st.sidebar.button("📋 Cargar datos de ejemplo"):
        astro_data = load_example_data()
        st.sidebar.success("✅ Datos de ejemplo cargados")

elif input_method == "Archivo Quantum ESPRESSO":
    uploaded_file = st.sidebar.file_uploader("Archivo de entrada QE (.in, .pw)", type=['in', 'pw', 'pwi'])
    if uploaded_file is not None:
        try:
            content = uploaded_file.read().decode('utf-8')
            parser = FileParser()
            astro_data = parser.parse_quantum_espresso(content)
            st.sidebar.success("✅ Archivo QE parseado")
        except Exception as e:
            st.sidebar.error(f"❌ Error: {e}")

elif input_method == "Archivo CIF":
    uploaded_file = st.sidebar.file_uploader("Archivo CIF", type=['cif'])
    if uploaded_file is not None:
        try:
            content = uploaded_file.read().decode('utf-8')
            parser = FileParser()
            astro_data = parser.parse_cif(content)
            st.sidebar.success("✅ Archivo CIF parseado")
        except Exception as e:
            st.sidebar.error(f"❌ Error: {e}")

elif input_method == "Archivo LAMMPS":
    uploaded_file = st.sidebar.file_uploader("Archivo de datos LAMMPS", type=['data', 'lmp'])
    if uploaded_file is not None:
        try:
            content = uploaded_file.read().decode('utf-8')
            parser = FileParser()
            astro_data = parser.parse_lammps_data(content)
            st.sidebar.success("✅ Archivo LAMMPS parseado")
        except Exception as e:
            st.sidebar.error(f"❌ Error: {e}")

# =========================================================================
# SIDEBAR - MATERIAL SELECTION
# =========================================================================

st.sidebar.header("🎯 Área de Aplicación")
area = st.sidebar.selectbox(
    "Selecciona el área:",
    ["Mecánica", "Química", "Industrial", "Energía", 
     "Aeroespacial", "Biomédico", "Electrónica", "Custom"]
)

st.sidebar.header("⚗️ Configuración del Material")

material_type = st.sidebar.selectbox(
    "Tipo de material:",
    ["Óxido metálico (TiO₂, ZnO)", "Metal (Ti, Al, Fe)", 
     "Cerámica (SiO₂, Al₂O₃)", "Polímero", "Composite",
     "Nanopartícula", "Material poroso", "Aleación"]
)

metal = st.sidebar.selectbox(
    "Metal específico:",
    ["Titanio (Ti)", "Aluminio (Al)", "Hierro (Fe)", 
     "Zinc (Zn)", "Cobre (Cu)", "Níquel (Ni)",
     "Cobalto (Co)", "Manganeso (Mn)", "Plata (Ag)",
     "Oro (Au)", "Platino (Pt)", "Paladio (Pd)"]
)

with st.sidebar.expander("⚙️ Opciones Avanzadas"):
    bulk_density = st.number_input("Densidad bulk (g/cm³)", value=3.0, min_value=0.1, max_value=20.0, step=0.1)
    bulk_conductivity = st.number_input("Conductividad bulk (W/m·K)", value=50.0, min_value=0.1, max_value=500.0, step=1.0)
    bulk_elastic = st.number_input("Módulo elástico bulk (GPa)", value=200.0, min_value=1.0, max_value=1000.0, step=5.0)
    supercell_size = st.slider("Tamaño de supercelda", value=5, min_value=2, max_value=10)

# =========================================================================
# MAIN CONTENT
# =========================================================================

if astro_data is not None:
    st.header("🌌 Firma Astrofísica Detectada")
    display_metrics_from_data(astro_data)
    
    if st.button("🚀 Generar Receta de Material", type="primary", use_container_width=True):
        with st.spinner("Calculando propiedades desde firmas astrofísicas..."):
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
                
                # Extract metal symbol
                metal_symbol = metal.split('(')[0].strip()
                material_type_clean = material_type.split('(')[0].strip().lower()
                
                # Calculate properties
                physical_props = physics_calc.calculate_all_properties(astro_data)
                
                # Generate recipe
                recipe = chem_engine.generate_complete_recipe(
                    metal=metal_symbol,
                    material_type=material_type_clean,
                    astrophysical_data=astro_data
                )
                
                # Display results
                tab1, tab2, tab3, tab4, tab5 = st.tabs([
                    "📊 Propiedades Físicas", "⚗️ Química de Síntesis", 
                    "📈 Visualización", "📁 Archivos Generados", "📋 Informe"
                ])
                
                with tab1:
                    st.subheader("Propiedades Calculadas")
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.markdown('<div class="metric-card">', unsafe_allow_html=True)
                        st.metric("Porosidad", f"{physical_props.get('porosity', 0):.1%}")
                        st.metric("Densidad", f"{physical_props.get('density', 0):.2f} g/cm³")
                        st.markdown('</div>', unsafe_allow_html=True)
                    
                    with col2:
                        st.markdown('<div class="metric-card">', unsafe_allow_html=True)
                        st.metric("Conductividad Térmica", f"{physical_props.get('thermal_conductivity', 0):.2f} W/m·K")
                        st.metric("Módulo Elástico", f"{physical_props.get('elastic_modulus', 0):.2f} GPa")
                        st.markdown('</div>', unsafe_allow_html=True)
                    
                    with col3:
                        st.markdown('<div class="metric-card">', unsafe_allow_html=True)
                        st.metric("Energía de Activación", f"{physical_props.get('activation_energy', 0):.3f} eV")
                        st.metric("Área Superficial", f"{physical_props.get('surface_area', 0):.1f} m²/g")
                        st.markdown('</div>', unsafe_allow_html=True)
                
                with tab2:
                    st.subheader("Receta de Síntesis Química")
                    conditions = recipe.get('reaction_conditions', {})
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown("### 🌡️ Condiciones de Reacción")
                        st.write(f"**Temperatura:** {conditions.get('temperature_c', 0):.0f}°C")
                        st.write(f"**Tiempo:** {conditions.get('reaction_time_hours', 0):.1f} horas")
                        st.write(f"**pH:** {conditions.get('ph', 7):.1f}")
                    
                    with col2:
                        st.markdown("### 🧪 Precursores")
                        for precursor in recipe.get('precursors', []):
                            st.write(f"• {precursor}")
                    
                    st.markdown("### 📝 Procedimiento")
                    for i, step in enumerate(recipe.get('step_by_step', []), 1):
                        st.write(f"**{i}.** {step}")
                
                with tab3:
                    st.subheader("Visualización de Estructura")
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
                            st.image(img2, caption="Resumen de Propiedades", use_column_width=True)
                    except Exception as e:
                        st.warning(f"Visualización requiere configuración adicional: {e}")
                
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
                
                with tab5:
                    st.subheader("Informe Técnico")
                    report = file_gen.generate_report_metadata(physical_props, recipe)
                    
                    st.markdown(f"""
                    **COSMICFORGE LAB - INFORME**
                    
                    **Fecha:** {report.get('timestamp', 'N/A')}  
                    **Objeto:** {report.get('object_name', 'N/A')}  
                    **Material:** {report.get('metal_name', 'N/A')} {report.get('material_type', '')}
                    
                    | Propiedad | Valor |
                    |-----------|-------|
                    | Porosidad | {report.get('porosity', 'N/A')} |
                    | Densidad | {report.get('density', 'N/A')} |
                    | Temperatura | {report.get('temperature', 'N/A')} |
                    | Tiempo de reacción | {report.get('reaction_time', 'N/A')} |
                    """)
                
                st.markdown('<div class="success-box">✅ Receta generada exitosamente!</div>', unsafe_allow_html=True)
            
            except Exception as e:
                st.error(f"❌ Error: {e}")

else:
    st.info("👈 Importa un archivo o usa los datos de ejemplo para comenzar")
    
    with st.expander("📖 Ver formato JSON esperado"):
        st.code(json.dumps(load_example_data(), indent=2), language="json")

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666;">
    <p><strong>CosmicForge Lab v2.0</strong> | Diseño de materiales inspirado en astrofísica</p>
</div>
""", unsafe_allow_html=True)
