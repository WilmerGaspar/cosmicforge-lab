# 🔬 CosmicForge Lab

[![Version](https://img.shields.io/badge/version-2.0.0-blue.svg)](https://github.com/WilmerGaspar/cosmicforge-lab)
[![Python](https://img.shields.io/badge/python-3.8%2B-brightgreen.svg)](https://www.python.org/)

**Diseño de Materiales Inspirado en Firmas Astrofísicas**

---

## 📖 Descripción

**CosmicForge Lab** es una plataforma científica que conecta la **astrofísica** con la **ciencia de materiales**. Transforma datos de objetos cósmicos en predicciones de propiedades de materiales y recetas de síntesis.

## ✨ Características

| Característica | Descripción |
|----------------|-------------|
| 🔭 **Importación** | JSON, Quantum ESPRESSO, CIF, LAMMPS, XYZ |
| ⚗️ **Propiedades** | Porosidad, densidad, conductividad, módulo elástico |
| 🧪 **Recetas** | Síntesis química paso a paso |
| 📊 **Visualización** | Estructuras cristalinas 2D/3D |
| 📤 **Exportación** | POSCAR, LAMMPS, XYZ, CSV, CIF |

## 🚀 Instalación

```bash
# Clonar repositorio
git clone https://github.com/WilmerGaspar/cosmicforge-lab.git
cd cosmicforge-lab

# Instalar dependencias
pip install -r requirements.txt

# Ejecutar
streamlit run app.py
```

## 📚 Uso

### Ejemplo de entrada JSON

```json
{
    "object_name": "Orion Nebula M42",
    "fractal_dimension": 1.654,
    "criticality_score": 0.722,
    "entropy": 0.019,
    "anisotropy": 0.329,
    "turbulence_beta": 2.278
}
```

### Uso programático

```python
from modules.physics_calculator import PhysicsCalculator
from modules.chemistry_engine import ChemistryEngine

# Datos astrofísicos
astro_data = {
    "object_name": "Crab Nebula",
    "fractal_dimension": 1.75,
    "criticality_score": 0.85,
    "entropy": 0.025,
    "anisotropy": 0.45,
    "turbulence_beta": 2.1
}

# Calcular propiedades
physics = PhysicsCalculator()
properties = physics.calculate_all_properties(astro_data)

# Generar receta
chemistry = ChemistryEngine()
recipe = chemistry.generate_complete_recipe(metal='Ti', astrophysical_data=astro_data)
```

## 🧩 Módulos

- **PhysicsCalculator**: Cálculo de propiedades físicas
- **ChemistryEngine**: Generación de recetas químicas
- **FileGenerator**: Archivos para simuladores (VASP, LAMMPS, QE)
- **Visualizer**: Visualizaciones gráficas
- **FileParser**: Lectura de archivos externos

## 🆕 Novedades v2.0

### Correcciones
- ✅ Corregido `_file_` → `__file__`
- ✅ Corregido `METAL_PROPERTIESetal]` → `METAL_PROPERTIES[metal]`
- ✅ Parámetros corregidos en file_generators
- ✅ Añadida función `visualize_atomic_structure`

### Mejoras
- 🆕 12 metales soportados (Ti, Al, Fe, Zn, Cu, Ni, Co, Mn, Ag, Au, Pt, Pd)
- 🆕 8 precursores químicos
- 🆕 Nuevo módulo `file_parser.py`
- 🆕 Cálculo de área superficial y band gap

## 📁 Estructura

```
cosmicforge-lab/
├── app.py
├── requirements.txt
├── README.md
└── modules/
    ├── __init__.py
    ├── physics_calculator.py
    ├── chemistry_engine.py
    ├── file_generators.py
    ├── visualizer.py
    └── file_parser.py
```

## 📧 Contacto

- **GitHub**: [@WilmerGaspar](https://github.com/WilmerGaspar)
- **Repo**: [cosmicforge-lab](https://github.com/WilmerGaspar/cosmicforge-lab)

---

**Desarrollado para la ciencia de materiales** 🔬✨
