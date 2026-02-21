"""
File Parser Module - Parse external file formats
Versión 1.0 - Nuevo módulo para CosmicForge Lab
"""

import numpy as np
from typing import Dict


class FileParser:
    """Parse various computational materials science file formats"""
    
    def parse_quantum_espresso(self, content: str) -> Dict:
        """Parse Quantum ESPRESSO input file"""
        data = {
            'source': 'Quantum ESPRESSO',
            'object_name': 'QE_Structure',
            'fractal_dimension': 1.5,
            'criticality_score': 0.5,
            'entropy': 0.01,
            'anisotropy': 0.3,
            'turbulence_beta': 2.0,
            'lyapunov_max': -0.1,
            'mode': 'computed',
            'raw_parameters': {}
        }
        
        lines = content.strip().split('\n')
        
        in_system = False
        for line in lines:
            line_lower = line.lower().strip()
            
            if '&system' in line_lower:
                in_system = True
                continue
            elif line_lower.startswith('/'):
                in_system = False
                continue
            
            if in_system:
                if 'ibrav' in line_lower:
                    try:
                        ibrav = int(line.split('=')[1].strip().strip(','))
                        data['anisotropy'] = 0.3 if ibrav == 1 else 0.5
                    except:
                        pass
                
                if 'a' in line_lower and '=' in line_lower and 'alpha' not in line_lower:
                    try:
                        a = float(line.split('=')[1].strip().strip(','))
                        data['fractal_dimension'] = np.clip(1.5 + (a - 4) * 0.1, 0, 3)
                    except:
                        pass
                
                if 'nat' in line_lower:
                    try:
                        nat = int(line.split('=')[1].strip().strip(','))
                        data['criticality_score'] = np.clip(nat / 100, 0.1, 1.0)
                        data['entropy'] = 0.01 + nat * 0.001
                    except:
                        pass
        
        data['turbulence_beta'] = 2.0 + data['anisotropy']
        return data
    
    def parse_cif(self, content: str) -> Dict:
        """Parse CIF file"""
        data = {
            'source': 'CIF',
            'object_name': 'CIF_Structure',
            'fractal_dimension': 1.5,
            'criticality_score': 0.5,
            'entropy': 0.01,
            'anisotropy': 0.3,
            'turbulence_beta': 2.0,
            'lyapunov_max': -0.1,
            'mode': 'computed',
            'raw_parameters': {}
        }
        
        lines = content.strip().split('\n')
        
        for line in lines:
            line_lower = line.lower().strip()
            
            if '_chemical_name_common' in line_lower:
                parts = line.split("'")
                if len(parts) >= 2:
                    data['object_name'] = parts[1]
            
            if '_chemical_formula_sum' in line_lower:
                parts = line.split("'")
                if len(parts) >= 2:
                    data['object_name'] = parts[1]
            
            if '_cell_length_a' in line_lower:
                try:
                    val = float(line.split()[1])
                    data['fractal_dimension'] = np.clip(1.5 + (val - 4) * 0.05, 0, 3)
                except:
                    pass
            
            if '_cell_angle_alpha' in line_lower:
                try:
                    val = float(line.split()[1])
                    data['anisotropy'] = abs(val - 90) / 90
                except:
                    pass
        
        data['turbulence_beta'] = 2.0 + data['anisotropy']
        return data
    
    def parse_lammps_data(self, content: str) -> Dict:
        """Parse LAMMPS data file"""
        data = {
            'source': 'LAMMPS',
            'object_name': 'LAMMPS_Structure',
            'fractal_dimension': 1.5,
            'criticality_score': 0.5,
            'entropy': 0.01,
            'anisotropy': 0.3,
            'turbulence_beta': 2.0,
            'lyapunov_max': -0.1,
            'mode': 'computed',
            'raw_parameters': {}
        }
        
        lines = content.strip().split('\n')
        
        for line in lines:
            line_lower = line.lower().strip()
            
            if 'atoms' in line_lower:
                try:
                    parts = line.split()
                    natoms = int(parts[0])
                    data['criticality_score'] = np.clip(natoms / 1000, 0.1, 1.0)
                    data['entropy'] = 0.01 + natoms * 0.0001
                    data['object_name'] = f"LAMMPS_{natoms}atoms"
                except:
                    pass
            
            if 'atom types' in line_lower:
                try:
                    ntypes = int(line.split()[0])
                    data['anisotropy'] = np.clip(ntypes / 10, 0.1, 0.9)
                except:
                    pass
        
        data['turbulence_beta'] = 2.0 + data['anisotropy']
        return data
    
    def parse_xyz(self, content: str) -> Dict:
        """Parse XYZ file"""
        data = {
            'source': 'XYZ',
            'object_name': 'XYZ_Structure',
            'fractal_dimension': 1.5,
            'criticality_score': 0.5,
            'entropy': 0.01,
            'anisotropy': 0.3,
            'turbulence_beta': 2.0,
            'lyapunov_max': -0.1,
            'mode': 'computed',
            'raw_parameters': {}
        }
        
        lines = content.strip().split('\n')
        
        if len(lines) >= 2:
            try:
                natoms = int(lines[0].strip())
                data['criticality_score'] = np.clip(natoms / 100, 0.1, 1.0)
                data['entropy'] = 0.01 + natoms * 0.001
            except:
                pass
            
            data['object_name'] = lines[1].strip()[:50]
        
        data['turbulence_beta'] = 2.0 + data['anisotropy']
        return data
