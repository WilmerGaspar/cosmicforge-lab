"""
Chemistry Engine Module - Generates chemical synthesis recipes
Versión 2.0 - Código corregido y mejorado
"""

import numpy as np
from typing import Dict, List, Optional


class ChemistryEngine:
    """Generate chemical synthesis recipes and conditions"""
    
    PRECURSORS = {
        'nitrate': {'formula': 'NO₃', 'mw': 62.0, 'valence': -1, 'category': 'oxidizer'},
        'hydroxide': {'formula': 'OH', 'mw': 17.0, 'valence': -1, 'category': 'base'},
        'oxide': {'formula': 'O', 'mw': 16.0, 'valence': -2, 'category': 'oxidizer'},
        'chloride': {'formula': 'Cl', 'mw': 35.5, 'valence': -1, 'category': 'halogen'},
        'sulfate': {'formula': 'SO₄', 'mw': 96.0, 'valence': -2, 'category': 'oxidizer'},
        'citrate': {'formula': 'C₆H₅O₇', 'mw': 189.0, 'valence': -3, 'category': 'chelator'},
        'acetate': {'formula': 'CH₃COO', 'mw': 59.0, 'valence': -1, 'category': 'organic'},
        'carbonate': {'formula': 'CO₃', 'mw': 60.0, 'valence': -2, 'category': 'base'}
    }
    
    # CORREGIDO: Base de datos completa de metales
    METAL_PROPERTIES = {
        'Ti': {'name': 'Titanio', 'mw': 47.867, 'oxidation_states': [2, 3, 4], 'common_oxidation': 4, 'color': 'white'},
        'Al': {'name': 'Aluminio', 'mw': 26.982, 'oxidation_states': [3], 'common_oxidation': 3, 'color': 'white'},
        'Fe': {'name': 'Hierro', 'mw': 55.845, 'oxidation_states': [2, 3], 'common_oxidation': 3, 'color': 'brown'},
        'Zn': {'name': 'Zinc', 'mw': 65.38, 'oxidation_states': [2], 'common_oxidation': 2, 'color': 'white'},
        'Cu': {'name': 'Cobre', 'mw': 63.546, 'oxidation_states': [1, 2], 'common_oxidation': 2, 'color': 'blue'},
        'Ni': {'name': 'Níquel', 'mw': 58.693, 'oxidation_states': [2, 3], 'common_oxidation': 2, 'color': 'green'},
        'Co': {'name': 'Cobalto', 'mw': 58.933, 'oxidation_states': [2, 3], 'common_oxidation': 2, 'color': 'pink'},
        'Mn': {'name': 'Manganeso', 'mw': 54.938, 'oxidation_states': [2, 3, 4], 'common_oxidation': 4, 'color': 'purple'},
        'Ag': {'name': 'Plata', 'mw': 107.868, 'oxidation_states': [1], 'common_oxidation': 1, 'color': 'white'},
        'Au': {'name': 'Oro', 'mw': 196.967, 'oxidation_states': [1, 3], 'common_oxidation': 3, 'color': 'yellow'},
        'Pt': {'name': 'Platino', 'mw': 195.084, 'oxidation_states': [2, 4], 'common_oxidation': 4, 'color': 'gray'},
        'Pd': {'name': 'Paladio', 'mw': 106.42, 'oxidation_states': [2, 4], 'common_oxidation': 2, 'color': 'white'}
    }
    
    def __init__(self, temperature_offset: float = 0.0):
        self.temperature_offset = temperature_offset
    
    def select_precursors(self, material_type: str, turbulence_beta: float) -> List[str]:
        base_precursors = {
            'óxido': ['nitrate', 'hydroxide'],
            'oxide': ['nitrate', 'hydroxide'],
            'metal': ['chloride', 'hydroxide'],
            'cerámica': ['nitrate', 'oxide', 'hydroxide'],
            'ceramic': ['nitrate', 'oxide', 'hydroxide'],
            'polímero': ['acetate', 'citrate'],
            'polymer': ['acetate', 'citrate'],
            'composite': ['nitrate', 'citrate', 'hydroxide'],
            'nanopartícula': ['nitrate', 'citrate', 'hydroxide'],
            'material poroso': ['nitrate', 'carbonate', 'citrate'],
            'aleación': ['nitrate', 'chloride', 'citrate']
        }
        
        precursors = base_precursors.get(material_type.lower(), ['nitrate', 'oxide'])
        
        if turbulence_beta > 2.5 and 'citrate' not in precursors:
            precursors.append('citrate')
        
        return precursors
    
    def generate_synthesis_equation(self, metal: str, precursors: List[str], temperature: float) -> Dict:
        # CORREGIDO: Validación correcta del metal
        if metal not in self.METAL_PROPERTIES:
            metal = 'Fe'  # Default
        
        # CORREGIDO: Acceso correcto al diccionario
        metal_props = self.METAL_PROPERTIES[metal]
        oxidation_state = metal_props['common_oxidation']
        
        precursor_formula = ' + '.join([self.PRECURSORS[p]['formula'] for p in precursors if p in self.PRECURSORS])
        
        products = f"{metal}Ox"
        equation = f"{metal}({precursor_formula})x → {products} + byproducts"
        
        stoichiometry = self._balance_equation(metal, precursors, oxidation_state)
        
        return {
            'equation': equation,
            'reactants': precursors,
            'products': products,
            'stoichiometry': stoichiometry,
            'temperature': temperature
        }
    
    def _balance_equation(self, metal: str, precursors: List[str], oxidation_state: int) -> Dict:
        if metal not in self.METAL_PROPERTIES:
            metal = 'Fe'
        
        metal_props = self.METAL_PROPERTIES[metal]
        
        stoich = {
            'metal': metal,
            'metal_coefficient': 1,
            'metal_amount_g': metal_props['mw'],
            'metal_name': metal_props['name']
        }
        
        for precursor in precursors:
            if precursor not in self.PRECURSORS:
                continue
            precursor_data = self.PRECURSORS[precursor]
            coefficient = abs(precursor_data['valence'] * oxidation_state) / max(abs(precursor_data['valence']), 1)
            
            stoich[precursor] = {
                'coefficient': coefficient,
                'amount_g': coefficient * precursor_data['mw'],
                'molar_ratio': coefficient,
                'formula': precursor_data['formula']
            }
        
        return stoich
    
    def calculate_reaction_conditions(self, criticality_score: float, turbulence_beta: float,
                                     entropy: float, anisotropy: float) -> Dict:
        T_base = 273.15 + (criticality_score * 200)
        T_turbulence_factor = 1 + (turbulence_beta - 2.0) * 0.05
        temperature_K = T_base * T_turbulence_factor + self.temperature_offset
        temperature_C = temperature_K - 273.15
        
        temperature_C = np.clip(temperature_C, 25, 500)
        temperature_K = temperature_C + 273.15
        
        pressure_atm = np.clip(1.0 + (anisotropy * 2.0), 1.0, 10.0)
        pH = np.clip(7.0 - (entropy * 3.0), 0, 14)
        reaction_time_hours = np.clip(1.0 + (1 - criticality_score) * 24, 0.5, 72)
        stirring_speed_rpm = np.clip(100 + (turbulence_beta * 200), 100, 1000)
        
        return {
            'temperature_K': round(temperature_K, 2),
            'temperature_C': round(temperature_C, 2),
            'pressure_atm': round(pressure_atm, 2),
            'pH': round(pH, 2),
            'reaction_time_hours': round(reaction_time_hours, 2),
            'stirring_speed_rpm': round(stirring_speed_rpm, 0)
        }
    
    def generate_complete_recipe(self, metal: str = 'Fe', material_type: str = 'oxide',
                                astrophysical_data: Optional[Dict] = None) -> Dict:
        if astrophysical_data is None:
            astrophysical_data = {
                'object_name': 'Unknown',
                'criticality_score': 0.5,
                'turbulence_beta': 2.0,
                'entropy': 0.01,
                'anisotropy': 0.5
            }
        
        required_fields = ['criticality_score', 'turbulence_beta', 'entropy', 'anisotropy']
        for field in required_fields:
            if field not in astrophysical_data:
                astrophysical_data[field] = 0.5
        
        precursors = self.select_precursors(material_type, astrophysical_data.get('turbulence_beta', 2.0))
        
        conditions = self.calculate_reaction_conditions(
            astrophysical_data['criticality_score'],
            astrophysical_data['turbulence_beta'],
            astrophysical_data['entropy'],
            astrophysical_data['anisotropy']
        )
        
        equation = self.generate_synthesis_equation(metal, precursors, conditions['temperature_C'])
        
        # CORREGIDO: Validación correcta del metal
        if metal in self.METAL_PROPERTIES:
            metal_color = self.METAL_PROPERTIES[metal]['color']
            metal_name = self.METAL_PROPERTIES[metal]['name']
        else:
            metal_color = 'unknown'
            metal_name = metal
        
        return {
            'object_name': astrophysical_data.get('object_name', 'Unknown'),
            'material_type': material_type,
            'metal': metal,
            'metal_name': metal_name,
            'base_metal_color': metal_color,
            'precursors': precursors,
            'synthesis_equation': equation,
            'reaction_conditions': conditions,
            'step_by_step': self._generate_procedure(precursors, conditions, metal),
            'safety_notes': self._generate_safety_notes(precursors),
            'expected_product': f"{metal} oxide/hydroxide nanoparticles"
        }
    
    def _generate_procedure(self, precursors: List[str], conditions: Dict, metal: str = 'Fe') -> List[str]:
        metal_name = self.METAL_PROPERTIES.get(metal, {}).get('name', metal)
        
        return [
            f"1. Medir todos los precursores químicos según las proporciones estequiométricas",
            f"2. Disolver la sal de {metal_name} en agua desionizada (100 mL)",
            "3. Agregar el agente precipitante/base gota a gota mientras agita",
            f"4. Calentar la solución a {conditions['temperature_C']:.0f}°C",
            f"5. Mantener la temperatura por {conditions['reaction_time_hours']:.1f} horas",
            f"6. Agitar a {conditions['stirring_speed_rpm']:.0f} RPM durante todo el proceso",
            "7. Permitir que la mezcla se enfríe a temperatura ambiente",
            "8. Filtrar el producto usando embudo Büchner",
            "9. Lavar el precipitado con agua destilada 3 veces",
            "10. Secar en estufa a 120°C por 2 horas",
            "11. Caracterizar el producto usando SEM, XRD, BET, TEM"
        ]
    
    def _generate_safety_notes(self, precursors: List[str]) -> List[str]:
        safety = [
            "Usar en campana de extracción bien ventilada",
            "Usar gafas de seguridad y guantes resistentes a químicos",
            "Tener disponible ducha de emergencia y lavaojos"
        ]
        
        for precursor in precursors:
            if precursor in self.PRECURSORS:
                safety.append(f"{precursor.capitalize()}: manejar con precaución")
        
        safety.append("Desechar residuos según regulaciones locales")
        
        return safety
    
    def get_metal_info(self, metal: str) -> Dict:
        if metal in self.METAL_PROPERTIES:
            return self.METAL_PROPERTIES[metal]
        return {'error': f"Metal {metal} not found"}
    
    def list_available_metals(self) -> List[str]:
        return list(self.METAL_PROPERTIES.keys())
