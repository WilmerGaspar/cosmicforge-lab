"""
File Generators Module - Creates simulation input files
Versión 2.0 - Código mejorado
"""

import numpy as np
from typing import Dict
from datetime import datetime


class FileGenerator:
    """Generate various file formats for simulations"""
    
    ELEMENTS = {
        'Fe': {'mass': 55.845}, 'Ti': {'mass': 47.867}, 'Al': {'mass': 26.982},
        'Zn': {'mass': 65.38}, 'Cu': {'mass': 63.546}, 'Ni': {'mass': 58.693},
        'Co': {'mass': 58.933}, 'O': {'mass': 15.999}, 'Si': {'mass': 28.086},
        'Mn': {'mass': 54.938}, 'Ag': {'mass': 107.868}, 'Au': {'mass': 196.967},
        'Pt': {'mass': 195.084}, 'Pd': {'mass': 106.42},
    }
    
    def __init__(self, lattice_constant: float = 4.0, supercell_size: int = 5):
        self.lattice_constant = lattice_constant
        self.supercell_size = supercell_size
    
    def generate_poscar(self, material_properties: Dict, recipe: Dict, material_type: str = 'oxide') -> str:
        """Generate POSCAR file for VASP/Quantum ESPRESSO"""
        porosity = material_properties.get('porosity', 0.3)
        object_name = material_properties.get('object_name', 'Material')
        metal = recipe.get('metal', 'Fe')
        
        n_atoms = int((1 - porosity) * self.supercell_size ** 3 * 4)
        n_atoms = max(n_atoms, 8)
        
        lattice = self.lattice_constant * self.supercell_size
        
        poscar = [
            f"# POSCAR - CosmicForge Lab",
            f"# Source: {object_name}",
            f"{object_name}",
            "1.0",
            f"  {lattice:.16f}  0.0  0.0",
            f"  0.0  {lattice:.16f}  0.0",
            f"  0.0  0.0  {lattice:.16f}",
            f"  {metal}  O",
            f"  {max(1, int(n_atoms*0.4))}  {n_atoms - max(1, int(n_atoms*0.4))}",
            "Cartesian"
        ]
        
        np.random.seed(hash(object_name) % (2**32))
        for i in range(n_atoms):
            x = np.random.uniform(0, lattice)
            y = np.random.uniform(0, lattice)
            z = np.random.uniform(0, lattice)
            poscar.append(f"  {x:.16f}  {y:.16f}  {z:.16f}")
        
        return "\n".join(poscar)
    
    def generate_lammps_input(self, material_properties: Dict, recipe: Dict, material_type: str = 'oxide') -> str:
        """Generate LAMMPS input file"""
        object_name = material_properties.get('object_name', 'Material')
        conditions = recipe.get('reaction_conditions', {})
        temp_k = conditions.get('temperature_K', 300)
        porosity = material_properties.get('porosity', 0.3)
        
        lattice = self.lattice_constant * self.supercell_size
        n_atoms = int((1 - porosity) * 100)
        
        lammps = [
            "# LAMMPS input - CosmicForge Lab",
            f"# Source: {object_name}",
            "units metal",
            "atom_style atomic",
            "boundary p p p",
            f"region box block 0 {lattice:.4f} 0 {lattice:.4f} 0 {lattice:.4f}",
            "create_box 2 box",
            "pair_style buck/coul/long 10.0",
            "pair_coeff 1 1 0.0 1.0 0.0",
            "pair_coeff 1 2 1288.0 2.768 0.0",
            "pair_coeff 2 2 22764.0 0.149 0.0",
            "kspace_style ewald 1.0e-4",
            f"velocity all create {temp_k} 12345",
            "timestep 0.001",
            "minimize 1.0e-4 1.0e-6 100 1000",
            f"fix 1 all nvt temp {temp_k} {temp_k} 0.1",
            "thermo 100",
            "dump 1 all xyz 100 trajectory.xyz",
            "run 1000",
            "write_data final_structure.data"
        ]
        
        return "\n".join(lammps)
    
    def generate_xyz(self, material_properties: Dict, recipe: Dict) -> str:
        """Generate XYZ file for visualization"""
        porosity = material_properties.get('porosity', 0.3)
        object_name = material_properties.get('object_name', 'Material')
        metal = recipe.get('metal', 'Fe')
        
        n_atoms = max(50, int((1 - porosity) * 200))
        
        xyz = [str(n_atoms), f"{object_name} - Porosity: {porosity:.3f}"]
        
        np.random.seed(hash(object_name) % (2**32))
        for i in range(n_atoms):
            atom = metal if i < n_atoms * 0.4 else "O"
            x = np.random.uniform(-10, 10)
            y = np.random.uniform(-10, 10)
            z = np.random.uniform(-10, 10)
            xyz.append(f"{atom} {x:.6f} {y:.6f} {z:.6f}")
        
        return "\n".join(xyz)
    
    def generate_properties_csv(self, material_properties: Dict, recipe: Dict) -> str:
        """Generate CSV file with material properties"""
        lines = [
            "# CosmicForge Lab - Material Properties",
            f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "## MATERIAL PROPERTIES",
            f"Porosity,{material_properties.get('porosity', 0):.6f}",
            f"Density [g/cm³],{material_properties.get('density', 0):.6f}",
            f"Thermal Conductivity [W/m·K],{material_properties.get('thermal_conductivity', 0):.6f}",
            f"Elastic Modulus [GPa],{material_properties.get('elastic_modulus', 0):.6f}",
            f"Surface Area [m²/g],{material_properties.get('surface_area', 0):.2f}",
            f"Band Gap [eV],{material_properties.get('band_gap', 0):.3f}",
            "",
            "## SYNTHESIS CONDITIONS"
        ]
        
        conditions = recipe.get('reaction_conditions', {})
        lines.extend([
            f"Temperature [°C],{conditions.get('temperature_C', 0):.2f}",
            f"Reaction Time [hours],{conditions.get('reaction_time_hours', 0):.2f}",
            f"pH,{conditions.get('ph', 0):.2f}"
        ])
        
        return "\n".join(lines)
    
    def generate_report_metadata(self, material_properties: Dict, recipe: Dict) -> Dict:
        """Generate metadata for reports"""
        conditions = recipe.get('reaction_conditions', {})
        
        return {
            'object_name': material_properties.get('object_name', 'Unknown'),
            'timestamp': datetime.now().isoformat(),
            'material_type': recipe.get('material_type', 'Unknown'),
            'metal': recipe.get('metal', 'Unknown'),
            'metal_name': recipe.get('metal_name', 'Unknown'),
            'porosity': f"{material_properties.get('porosity', 0):.4f}",
            'density': f"{material_properties.get('density', 0):.4f} g/cm³",
            'temperature': f"{conditions.get('temperature_C', 0):.1f} °C",
            'reaction_time': f"{conditions.get('reaction_time_hours', 0):.2f} hours"
        }
