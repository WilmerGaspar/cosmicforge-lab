#!/usr/bin/env python3
"""
CosmicForge Lab - Tarea LAMMPS
Generación de archivos de entrada para simulación Molecular Dynamics

Archivos generados:
- input.in: Archivo de entrada principal
- data.lammps: Archivo de estructura molecular

Ejemplo de uso:
    from tasks.task_lammps import LAMMPSTask
    
    task = LAMMPSTask(metal="Ti", physical_props={"density": 4.5})
    files = task.generate_all()
    task.save_files(output_dir="./lammps_calculation")
"""

import os
from pathlib import Path
from typing import Dict, Optional

class LAMMPSTask:
    """Genera archivos de entrada para simulaciones LAMMPS."""
    
    ELEMENT_MASS = {
        "Ti": 47.867, "Al": 26.982, "Fe": 55.845, "Zn": 65.38,
        "Cu": 63.546, "Ni": 58.693, "Co": 58.933, "Mn": 54.938,
        "Cr": 52.00, "Mo": 95.95, "W": 183.84, "Ag": 107.868,
        "Au": 196.967, "Pt": 195.084, "Pd": 106.42, "O": 15.999,
    }
    
    def __init__(self, metal: str, physical_props: Optional[Dict] = None, 
                 sim_type: str = "md"):
        self.metal = metal.capitalize() if len(metal) == 2 else metal.upper()
        self.physical_props = physical_props or {}
        self.sim_type = sim_type
        self.mass_metal = self.ELEMENT_MASS.get(self.metal, 50.0)
    
    def generate_input(self, temperature: float = 300.0, 
                       n_steps: int = 50000) -> str:
        """Genera el archivo de entrada input.in."""
        
        density = self.physical_props.get("density", 4.5)
        lattice_const = (4.0 / density) ** (1/3) * 3.0
        
        return f"""# LAMMPS Input - {self.metal}O2
# Generado por CosmicForge Lab
# Simulación: {self.sim_type.upper()}

# ================== INICIALIZACIÓN ==================
units           metal
dimension       3
boundary        p p p
atom_style      full
neighbor        2.0 bin

# ================== ESTRUCTURA ==================
lattice         fcc {lattice_const:.4f}
region          box block 0 10 0 10 0 10
create_box      2 box
create_atoms    1 box

# ================== MASA ==================
mass            1 {self.mass_metal:.3f}   # {self.metal}
mass            2 15.999   # O

# ================== POTENCIAL ==================
pair_style      buck/coul/long 12.0
pair_coeff      1 1 0.001 0.1 0.0    # {self.metal}-{self.metal}
pair_coeff      2 2 0.001 0.1 0.0    # O-O
pair_coeff      1 2 0.001 0.2 0.0    # {self.metal}-O

kspace_style    pppm 1.0e-4

# ================== CONFIGURACIÓN ==================
velocity        all create {temperature} 12345 mom yes rot yes dist gaussian

# ================== MINIMIZACIÓN ==================
minimize        1.0e-4 1.0e-6 100 1000
reset_timestep  0

# ================== TERMO ==================
thermo          1000
thermo_style    custom step temp pe ke etotal press vol density

# ================== EQUILIBRACIÓN NVT ==================
fix             1 all nvt temp {temperature} {temperature} 0.1
timestep        0.001
run             10000

# ================== PRODUCCIÓN ==================
unfix           1
fix             2 all npt temp {temperature} {temperature} 0.1 iso 0 0 1.0
run             {n_steps}

# ================== OUTPUT ==================
write_data      {self.metal}O2_final.data
write_restart   {self.metal}O2.restart
"""

    def generate_data(self) -> str:
        """Genera archivo de datos LAMMPS."""
        density = self.physical_props.get("density", 4.5)
        lattice = (4.0 / density) ** (1/3) * 3.0
        
        return f"""LAMMPS data file - {self.metal}O2 - CosmicForge Lab

1000 atoms
2 atom types

0.0 {lattice * 10:.4f} xlo xhi
0.0 {lattice * 10:.4f} ylo yhi
0.0 {lattice * 10:.4f} zlo zhi

Masses

1 {self.mass_metal:.3f}  # {self.metal}
2 15.999  # O

Atoms  # full

1 1 1  0.0  0.0  0.0  0.0
2 1 1  0.0  1.0  1.0  0.0
...
"""

    def generate_minimization(self) -> str:
        """Genera input para minimización de energía."""
        return f"""# LAMMPS Minimización - {self.metal}O2
# CosmicForge Lab

units           metal
dimension       3
boundary        p p p
atom_style      full

lattice         fcc 4.0
region          box block 0 5 0 5 0 5
create_box      2 box

mass 1 {self.mass_metal:.3f}
mass 2 15.999

pair_style      buck/coul/long 12.0
pair_coeff      * * 0.001 0.1 0.0
kspace_style    pppm 1.0e-4

minimize        1.0e-6 1.0e-8 10000 100000
write_data      minimized.data
"""

    def generate_npt(self, temperature: float = 300.0, 
                     pressure: float = 0.0) -> str:
        """Genera input para NPT ensemble."""
        return self.generate_input(temperature, 100000)  # Más pasos para NPT

    def generate_all(self) -> Dict[str, str]:
        """Genera todos los archivos."""
        return {
            "input.in": self.generate_input(),
            "data.lammps": self.generate_data(),
            "minimize.in": self.generate_minimization(),
        }
    
    def save_files(self, output_dir: str = "./lammps_calc") -> None:
        """Guarda todos los archivos."""
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        for filename, content in self.generate_all().items():
            with open(f"{output_dir}/{filename}", "w") as f:
                f.write(content)
        print(f"✓ Archivos LAMMPS guardados en {output_dir}/")


# EJEMPLO DE USO
if __name__ == "__main__":
    print("=" * 50)
    print("CosmicForge Lab - LAMMPS Task Example")
    print("=" * 50)
    
    task = LAMMPSTask(metal="Ti", physical_props={"density": 4.5})
    
    print("\n--- input.in ---")
    print(task.generate_input())
    
    # task.save_files("./lammps_TiO2")
