#!/usr/bin/env python3
"""
CosmicForge Lab - Tarea Quantum ESPRESSO
Generación de archivos de entrada para cálculos DFT

Archivos generados:
- pw.in: Input para pw.x (SCF, relajación, bandas)
- bands.in: Input para bands.x
- dos.in: Input para dos.x
- ph.in: Input para ph.x (fonones)

Ejemplo de uso:
    from tasks.task_quantum_espresso import QuantumESPRESSOTask
    
    task = QuantumESPRESSOTask(metal="Ti", physical_props={"density": 4.5})
    files = task.generate_all()
    task.save_files(output_dir="./qe_calculation")
"""

import os
from pathlib import Path
from typing import Dict, Optional

class QuantumESPRESSOTask:
    """Genera archivos de entrada para Quantum ESPRESSO."""
    
    ELEMENTS = {
        "Ti": {"mw": 47.867, "pseudo": "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF"},
        "Al": {"mw": 26.982, "pseudo": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"},
        "Fe": {"mw": 55.845, "pseudo": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"},
        "Zn": {"mw": 65.38, "pseudo": "Zn.pbe-dn-kjpaw_psl.1.0.0.UPF"},
        "Cu": {"mw": 63.546, "pseudo": "Cu.pbe-dn-kjpaw_psl.1.0.0.UPF"},
        "Ni": {"mw": 58.693, "pseudo": "Ni.pbe-n-kjpaw_psl.1.0.0.UPF"},
        "Co": {"mw": 58.933, "pseudo": "Co.pbe-n-kjpaw_psl.1.0.0.UPF"},
        "Mn": {"mw": 54.938, "pseudo": "Mn.pbe-spn-kjpaw_psl.1.0.0.UPF"},
        "Cr": {"mw": 52.00, "pseudo": "Cr.pbe-spn-kjpaw_psl.1.0.0.UPF"},
        "Mo": {"mw": 95.95, "pseudo": "Mo.pbe-spn-kjpaw_psl.1.0.0.UPF"},
        "W": {"mw": 183.84, "pseudo": "W.pbe-spn-kjpaw_psl.1.0.0.UPF"},
        "Ag": {"mw": 107.868, "pseudo": "Ag.pbe-n-kjpaw_psl.1.0.0.UPF"},
        "Au": {"mw": 196.967, "pseudo": "Au.pbe-n-kjpaw_psl.1.0.0.UPF"},
        "Pt": {"mw": 195.084, "pseudo": "Pt.pbe-n-kjpaw_psl.1.0.0.UPF"},
        "Pd": {"mw": 106.42, "pseudo": "Pd.pbe-n-kjpaw_psl.1.0.0.UPF"},
    }
    
    def __init__(self, metal: str, physical_props: Optional[Dict] = None,
                 calc_type: str = "scf"):
        self.metal = metal.capitalize() if len(metal) == 2 else metal.upper()
        self.physical_props = physical_props or {}
        self.calc_type = calc_type
        self.element_data = self.ELEMENTS.get(self.metal, 
                                               {"mw": 50.0, "pseudo": f"{self.metal}.pbe-n-kjpaw_psl.1.0.0.UPF"})
    
    def generate_scf(self) -> str:
        """Genera input para cálculo SCF."""
        return f"""&CONTROL
  calculation = 'scf'
  prefix = '{self.metal}O2'
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
 {self.metal} {self.element_data['mw']:.3f} {self.element_data['pseudo']}
 O  15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
 {self.metal} 0.000 0.000 0.000
 O  0.305 0.305 0.000
 O -0.305 -0.305 0.000

K_POINTS automatic
 8 8 8 0 0 0
"""

    def generate_relax(self) -> str:
        """Genera input para relajación estructural."""
        return f"""&CONTROL
  calculation = 'vc-relax'
  prefix = '{self.metal}O2'
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

&IONS
  ion_dynamics = 'bfgs'
/

&CELL
  cell_dynamics = 'bfgs'
  press = 0.0
/

ATOMIC_SPECIES
 {self.metal} {self.element_data['mw']:.3f} {self.element_data['pseudo']}
 O  15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
 {self.metal} 0.000 0.000 0.000
 O  0.305 0.305 0.000
 O -0.305 -0.305 0.000

K_POINTS automatic
 6 6 6 0 0 0

CELL_PARAMETERS cubic
 8.0 0.0 0.0
 0.0 8.0 0.0
 0.0 0.0 8.0
"""

    def generate_bands(self) -> str:
        """Genera input para estructura de bandas."""
        return f"""&CONTROL
  calculation = 'bands'
  prefix = '{self.metal}O2'
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
  nbnd = 20
  occupations = 'smearing'
  smearing = 'mp'
  degauss = 0.02
/

&ELECTRONS
  conv_thr = 1.0d-8
/

ATOMIC_SPECIES
 {self.metal} {self.element_data['mw']:.3f} {self.element_data['pseudo']}
 O  15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
 {self.metal} 0.000 0.000 0.000
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

    def generate_dos(self) -> str:
        """Genera input para densidad de estados."""
        return f"""&DOS
  prefix = '{self.metal}O2'
  outdir = './tmp'
  fildos = '{self.metal}O2.dos'
  Emin = -10.0
  Emax = 10.0
  DeltaE = 0.01
/
"""

    def generate_phonon(self) -> str:
        """Genera input para cálculo de fonones."""
        return f"""&INPUTPH
  prefix = '{self.metal}O2'
  outdir = './tmp'
  fildyn = '{self.metal}O2.dyn'
  ldisp = .true.
  nq1 = 2
  nq2 = 2
  nq3 = 2
/
"""

    def generate_projwfc(self) -> str:
        """Genera input para proyección de densidad de estados."""
        return f"""&PROJWFC
  prefix = '{self.metal}O2'
  outdir = './tmp'
  filpdos = '{self.metal}O2.pdos'
  Emin = -10.0
  Emax = 10.0
  DeltaE = 0.01
/
"""

    def generate_run_script(self) -> str:
        """Genera script de ejecución."""
        return f"""#!/bin/bash
# Script de ejecución Quantum ESPRESSO
# Generado por CosmicForge Lab

# Crear directorios
mkdir -p tmp pseudo

# Ejecutar SCF
echo "Ejecutando SCF..."
mpirun -np 4 pw.x < {self.metal}O2_scf.in > {self.metal}O2_scf.out

# Ejecutar relajación (opcional)
# mpirun -np 4 pw.x < {self.metal}O2_relax.in > {self.metal}O2_relax.out

# Ejecutar bandas (opcional)
# mpirun -np 4 pw.x < {self.metal}O2_bands.in > {self.metal}O2_bands.out
# mpirun -np 4 bands.x < bands.in > bands.out

# Ejecutar DOS (opcional)
# mpirun -np 4 dos.x < dos.in > dos.out

echo "Cálculo completado"
"""

    def generate_all(self) -> Dict[str, str]:
        """Genera todos los archivos."""
        files = {
            f"{self.metal}O2_scf.in": self.generate_scf(),
            f"{self.metal}O2_relax.in": self.generate_relax(),
            f"{self.metal}O2_bands.in": self.generate_bands(),
            "dos.in": self.generate_dos(),
            "phonon.in": self.generate_phonon(),
            "projwfc.in": self.generate_projwfc(),
            "run.sh": self.generate_run_script(),
        }
        return files
    
    def save_files(self, output_dir: str = "./qe_calc") -> None:
        """Guarda todos los archivos."""
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        for filename, content in self.generate_all().items():
            with open(f"{output_dir}/{filename}", "w") as f:
                f.write(content)
        print(f"✓ Archivos Quantum ESPRESSO guardados en {output_dir}/")


# EJEMPLO DE USO
if __name__ == "__main__":
    print("=" * 50)
    print("CosmicForge Lab - Quantum ESPRESSO Task")
    print("=" * 50)
    
    task = QuantumESPRESSOTask(metal="Ti")
    
    print("\n--- SCF Input ---")
    print(task.generate_scf())
    
    print("\n--- Bands Input ---")
    print(task.generate_bands())
    
    # task.save_files("./qe_TiO2")
