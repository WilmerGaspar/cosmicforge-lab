#!/usr/bin/env python3
"""
CosmicForge Lab - Generador nanoHUB Web Interface
==================================================
Genera el formato EXACTO que la interfaz web de nanoHUB espera.

Formato nanoHUB Web:
- Title of Run
- Atomic Structure: descripción + coordenadas fraccionales  
- Cell Vectors (Å)
- Lattice Parameter a
- Ratio b/a, c/a

Autor: Wilmer Gaspar
"""

from typing import Dict, List, Tuple
from dataclasses import dataclass

@dataclass
class NanoHUBStructure:
    """Estructura cristalina para nanoHUB"""
    id: str
    name: str
    formula: str
    description: str
    structure_type: str
    lattice_a: float  # Å
    ratio_ba: float
    ratio_ca: float
    cell_vectors: Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]
    atoms: List[Tuple[str, float, float, float]]

# Estructuras cristalinas con formato nanoHUB
NANOHUB_STRUCTURES = {
    "Si_diamond": NanoHUBStructure(
        id="Si_diamond",
        name="Si diamond",
        formula="Si",
        description="Silicon diamond structure",
        structure_type="cubic F (fcc)",
        lattice_a=5.43,
        ratio_ba=1.0,
        ratio_ca=1.0,
        cell_vectors=(
            (-2.715, 2.715, 0.0),
            (-2.715, 0.0, 2.715),
            (0.0, 2.715, 2.715)
        ),
        atoms=[
            ("Si", 0.0, 0.0, 0.0),
            ("Si", 0.25, 0.25, 0.25)
        ]
    ),
    
    "TiO2_rutile": NanoHUBStructure(
        id="TiO2_rutile",
        name="TiO2 rutile",
        formula="TiO2",
        description="TiO2 rutile structure",
        structure_type="tetragonal P",
        lattice_a=4.594,
        ratio_ba=1.0,
        ratio_ca=0.644,
        cell_vectors=(
            (4.594, 0.0, 0.0),
            (0.0, 4.594, 0.0),
            (0.0, 0.0, 2.959)
        ),
        atoms=[
            ("Ti", 0.0, 0.0, 0.0),
            ("Ti", 0.5, 0.5, 0.5),
            ("O", 0.3053, 0.3053, 0.0),
            ("O", 0.6947, 0.6947, 0.0),
            ("O", 0.8053, 0.1947, 0.5),
            ("O", 0.1947, 0.8053, 0.5)
        ]
    ),
    
    "TiO2_anatase": NanoHUBStructure(
        id="TiO2_anatase",
        name="TiO2 anatase",
        formula="TiO2",
        description="TiO2 anatase structure",
        structure_type="tetragonal I (bct)",
        lattice_a=3.785,
        ratio_ba=1.0,
        ratio_ca=2.514,
        cell_vectors=(
            (3.785, 0.0, 0.0),
            (0.0, 3.785, 0.0),
            (0.0, 0.0, 9.514)
        ),
        atoms=[
            ("Ti", 0.0, 0.0, 0.0),
            ("Ti", 0.5, 0.5, 0.5),
            ("Ti", 0.0, 0.5, 0.25),
            ("Ti", 0.5, 0.0, 0.75),
            ("O", 0.0, 0.0, 0.208),
            ("O", 0.0, 0.0, 0.792),
            ("O", 0.5, 0.5, 0.292),
            ("O", 0.5, 0.5, 0.708),
            ("O", 0.0, 0.5, 0.042),
            ("O", 0.0, 0.5, 0.458),
            ("O", 0.5, 0.0, 0.542),
            ("O", 0.5, 0.0, 0.958)
        ]
    ),
    
    "Ti2O3_corundum": NanoHUBStructure(
        id="Ti2O3_corundum",
        name="Ti2O3 corundum",
        formula="Ti2O3",
        description="Ti2O3 corundum structure",
        structure_type="trigonal R",
        lattice_a=5.148,
        ratio_ba=1.0,
        ratio_ca=2.649,
        cell_vectors=(
            (5.148, 0.0, 0.0),
            (-2.574, 4.458, 0.0),
            (0.0, 0.0, 13.642)
        ),
        atoms=[
            ("Ti", 0.0, 0.0, 0.0),
            ("Ti", 0.0, 0.0, 0.333),
            ("Ti", 0.333, 0.667, 0.167),
            ("Ti", 0.667, 0.333, 0.500),
            ("O", 0.306, 0.0, 0.083),
            ("O", 0.0, 0.306, 0.083),
            ("O", 0.694, 0.0, 0.250),
            ("O", 0.0, 0.694, 0.250),
            ("O", 0.306, 0.0, 0.417),
            ("O", 0.0, 0.306, 0.417)
        ]
    ),
    
    "TiO_rocksalt": NanoHUBStructure(
        id="TiO_rocksalt",
        name="TiO rocksalt",
        formula="TiO",
        description="TiO rocksalt structure",
        structure_type="cubic F (fcc)",
        lattice_a=4.177,
        ratio_ba=1.0,
        ratio_ca=1.0,
        cell_vectors=(
            (2.089, 2.089, 0.0),
            (2.089, 0.0, 2.089),
            (0.0, 2.089, 2.089)
        ),
        atoms=[
            ("Ti", 0.0, 0.0, 0.0),
            ("O", 0.5, 0.5, 0.5)
        ]
    ),
    
    "ZnO_wurtzite": NanoHUBStructure(
        id="ZnO_wurtzite",
        name="ZnO wurtzite",
        formula="ZnO",
        description="ZnO wurtzite structure",
        structure_type="hexagonal P",
        lattice_a=3.250,
        ratio_ba=1.0,
        ratio_ca=1.603,
        cell_vectors=(
            (3.250, 0.0, 0.0),
            (-1.625, 2.814, 0.0),
            (0.0, 0.0, 5.207)
        ),
        atoms=[
            ("Zn", 0.333, 0.667, 0.0),
            ("Zn", 0.667, 0.333, 0.5),
            ("O", 0.333, 0.667, 0.375),
            ("O", 0.667, 0.333, 0.875)
        ]
    ),
    
    "Fe2O3_hematite": NanoHUBStructure(
        id="Fe2O3_hematite",
        name="Fe2O3 hematite",
        formula="Fe2O3",
        description="Fe2O3 hematite structure",
        structure_type="trigonal R",
        lattice_a=5.038,
        ratio_ba=1.0,
        ratio_ca=2.734,
        cell_vectors=(
            (5.038, 0.0, 0.0),
            (-2.519, 4.363, 0.0),
            (0.0, 0.0, 13.772)
        ),
        atoms=[
            ("Fe", 0.0, 0.0, 0.0),
            ("Fe", 0.0, 0.0, 0.333),
            ("Fe", 0.333, 0.667, 0.167),
            ("Fe", 0.667, 0.333, 0.500),
            ("O", 0.306, 0.0, 0.083),
            ("O", 0.0, 0.306, 0.083),
            ("O", 0.694, 0.0, 0.250),
            ("O", 0.0, 0.694, 0.250),
            ("O", 0.306, 0.0, 0.417),
            ("O", 0.0, 0.306, 0.417)
        ]
    ),
    
    "Al2O3_corundum": NanoHUBStructure(
        id="Al2O3_corundum",
        name="Al2O3 corundum",
        formula="Al2O3",
        description="Al2O3 corundum structure",
        structure_type="trigonal R",
        lattice_a=4.759,
        ratio_ba=1.0,
        ratio_ca=2.730,
        cell_vectors=(
            (4.759, 0.0, 0.0),
            (-2.380, 4.122, 0.0),
            (0.0, 0.0, 12.991)
        ),
        atoms=[
            ("Al", 0.0, 0.0, 0.0),
            ("Al", 0.0, 0.0, 0.333),
            ("Al", 0.333, 0.667, 0.167),
            ("Al", 0.667, 0.333, 0.500),
            ("O", 0.306, 0.0, 0.083),
            ("O", 0.0, 0.306, 0.083),
            ("O", 0.694, 0.0, 0.250),
            ("O", 0.0, 0.694, 0.250),
            ("O", 0.306, 0.0, 0.417),
            ("O", 0.0, 0.306, 0.417)
        ]
    ),
    
    "Cu2O_cuprite": NanoHUBStructure(
        id="Cu2O_cuprite",
        name="Cu2O cuprite",
        formula="Cu2O",
        description="Cu2O cuprite structure",
        structure_type="cubic P",
        lattice_a=4.267,
        ratio_ba=1.0,
        ratio_ca=1.0,
        cell_vectors=(
            (4.267, 0.0, 0.0),
            (0.0, 4.267, 0.0),
            (0.0, 0.0, 4.267)
        ),
        atoms=[
            ("O", 0.0, 0.0, 0.0),
            ("Cu", 0.25, 0.25, 0.25),
            ("Cu", 0.75, 0.75, 0.75)
        ]
    ),
    
    "NiO_rocksalt": NanoHUBStructure(
        id="NiO_rocksalt",
        name="NiO rocksalt",
        formula="NiO",
        description="NiO rocksalt structure",
        structure_type="cubic F (fcc)",
        lattice_a=4.177,
        ratio_ba=1.0,
        ratio_ca=1.0,
        cell_vectors=(
            (2.089, 2.089, 0.0),
            (2.089, 0.0, 2.089),
            (0.0, 2.089, 2.089)
        ),
        atoms=[
            ("Ni", 0.0, 0.0, 0.0),
            ("O", 0.5, 0.5, 0.5)
        ]
    ),
    
    "Co3O4_spinel": NanoHUBStructure(
        id="Co3O4_spinel",
        name="Co3O4 spinel",
        formula="Co3O4",
        description="Co3O4 spinel structure",
        structure_type="cubic F (fcc)",
        lattice_a=8.084,
        ratio_ba=1.0,
        ratio_ca=1.0,
        cell_vectors=(
            (8.084, 0.0, 0.0),
            (0.0, 8.084, 0.0),
            (0.0, 0.0, 8.084)
        ),
        atoms=[
            ("Co", 0.0, 0.0, 0.0),
            ("Co", 0.125, 0.125, 0.125),
            ("O", 0.25, 0.25, 0.25)
        ]
    )
}

def generate_nanohub_format(structure_id: str, object_name: str = None) -> Dict:
    """
    Genera el formato para nanoHUB web interface.
    
    Args:
        structure_id: ID de la estructura (ej: "TiO2_rutile")
        object_name: Nombre del objeto astrofísico (opcional)
    
    Returns:
        Diccionario con los campos para nanoHUB
    """
    structure = NANOHUB_STRUCTURES.get(structure_id)
    if not structure:
        raise ValueError(f"Estructura no encontrada: {structure_id}")
    
    # Title
    title = f"{structure.formula} band structure"
    if object_name:
        title = f"{structure.formula} - {object_name}"
    
    # Atomic Structure (descripción + coordenadas)
    atomic_lines = [structure.description]
    for elem, x, y, z in structure.atoms:
        atomic_lines.append(f"{elem} {x:.4f} {y:.4f} {z:.4f}")
    atomic_structure = "\n".join(atomic_lines)
    
    # Cell Vectors
    vector_lines = []
    for vec in structure.cell_vectors:
        vector_lines.append(f"{vec[0]:.3f} {vec[1]:.3f} {vec[2]:.3f}")
    cell_vectors = "\n".join(vector_lines)
    
    return {
        "title": title,
        "premade_structure": structure.name,
        "atomic_coordinates": "Fractional",
        "structure_type": structure.structure_type,
        "atomic_structure": atomic_structure,
        "cell_vectors": cell_vectors,
        "lattice_a": structure.lattice_a,
        "ratio_ba": structure.ratio_ba,
        "ratio_ca": structure.ratio_ca,
        "nat": len(structure.atoms),
        "formula": structure.formula
    }

def get_available_structures() -> List[Dict]:
    """Retorna lista de estructuras disponibles."""
    return [
        {
            "id": s.id,
            "name": s.name,
            "formula": s.formula,
            "nat": len(s.atoms)
        }
        for s in NANOHUB_STRUCTURES.values()
    ]

def print_nanohub_instructions():
    """Imprime instrucciones para usar en nanoHUB."""
    return """
================================================================================
CÓMO USAR EN nanoHUB
================================================================================

1. Ve a: https://nanohub.org/tools/qe

2. En la pestaña "Input Geometry":
   - Premade atomistic structure: Selecciona "User defined"
   - Atomic Coordinates: Selecciona "Fractional"

3. Copia cada campo:

   a) Title of Run:
      [Pega el título generado]
   
   b) Atomic Structure:
      [Pega todo el contenido - descripción + coordenadas]
   
   c) Cell Vectors (Å):
      [Pega las 3 líneas de vectores]
   
   d) Lattice Parameter a (Å):
      [Ingresa el valor]
   
   e) Ratio b/a:
      [Ingresa el valor]
   
   f) Ratio c/a:
      [Ingresa el valor]

4. Configura los demás parámetros:
   - Energy cutoff: 60 Ry
   - K-points: 6x6x6
   - Smearing: mp

5. Ejecuta la simulación

================================================================================
"""

# Test
if __name__ == "__main__":
    print("=" * 60)
    print("Generador nanoHUB Web Interface - Test")
    print("=" * 60)
    
    # Generar formato para Si diamond
    result = generate_nanohub_format("Si_diamond", "Test Object")
    
    print("\n--- FORMATO GENERADO ---\n")
    print(f"Title: {result['title']}")
    print(f"\nAtomic Structure:\n{result['atomic_structure']}")
    print(f"\nCell Vectors:\n{result['cell_vectors']}")
    print(f"\nLattice a: {result['lattice_a']} Å")
    print(f"Ratio b/a: {result['ratio_ba']}")
    print(f"Ratio c/a: {result['ratio_ca']}")
    
    print(print_nanohub_instructions())
