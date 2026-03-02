"""
Physics Calculator Module - Calculates material properties
Versión 2.0 - Código mejorado
"""

import numpy as np
from typing import Dict, Optional


class PhysicsCalculator:
    """Calculate material properties based on astrophysical signatures"""
    
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    k_B_eV = 8.617333262e-5  # Boltzmann constant (eV/K)
    
    def __init__(self, bulk_density: float = 3.0, bulk_conductivity: float = 50.0, 
                 bulk_elastic: float = 200.0, tau_0: float = 1.0, xi_0: float = 10.0):
        """
        Initialize physics calculator with bulk material properties
        
        Args:
            bulk_density: Bulk density in g/cm³
            bulk_conductivity: Bulk thermal conductivity in W/m·K
            bulk_elastic: Bulk elastic modulus in GPa
            tau_0: Reference process time in seconds
            xi_0: Reference correlation scale in nm
        """
        if bulk_density <= 0:
            raise ValueError("bulk_density must be positive")
        if bulk_conductivity <= 0:
            raise ValueError("bulk_conductivity must be positive")
        if bulk_elastic <= 0:
            raise ValueError("bulk_elastic must be positive")
        
        self.rho_bulk = bulk_density
        self.k_bulk = bulk_conductivity
        self.E_bulk = bulk_elastic
        self.tau_0 = tau_0
        self.xi_0 = xi_0
    
    def calculate_porosity(self, fractal_dimension: float) -> float:
        """Calculate porosity: Φ = 1 - (D/3)^1.5"""
        fractal_dimension = np.clip(fractal_dimension, 0, 3)
        
        if fractal_dimension <= 0:
            return 1.0
        if fractal_dimension >= 3:
            return 0.0
        
        porosity = 1 - (fractal_dimension / 3) ** 1.5
        return float(np.clip(porosity, 0, 1))
    
    def calculate_density(self, fractal_dimension: float) -> float:
        """Calculate effective density: ρ = ρ_bulk × (1-Φ)"""
        porosity = self.calculate_porosity(fractal_dimension)
        return float(self.rho_bulk * (1 - porosity))
    
    def calculate_thermal_conductivity(self, fractal_dimension: float) -> float:
        """Calculate thermal conductivity: k = k_bulk × (1-Φ)^1.5"""
        porosity = self.calculate_porosity(fractal_dimension)
        return float(self.k_bulk * (1 - porosity) ** 1.5)
    
    def calculate_elastic_modulus(self, fractal_dimension: float) -> float:
        """Calculate elastic modulus: E = E_bulk × (1-Φ)^2"""
        porosity = self.calculate_porosity(fractal_dimension)
        return float(self.E_bulk * (1 - porosity) ** 2)
    
    def calculate_process_time(self, entropy: float, entropy_max: float = 1.0) -> float:
        """Calculate process time: τ = τ₀ × exp(S/S_max)"""
        if entropy_max <= 0:
            entropy_max = 1.0
        entropy = np.clip(entropy, 0, entropy_max)
        return float(self.tau_0 * np.exp(entropy / entropy_max))
    
    def calculate_correlation_length(self, criticality_score: float) -> float:
        """Calculate correlation length: ξ = ξ₀ / criticality_score"""
        criticality_score = np.clip(criticality_score, 0.001, 1.0)
        return float(self.xi_0 / criticality_score)
    
    def calculate_activation_energy(self, anisotropy: float, temperature: float = 298.15) -> float:
        """Calculate activation energy: Ea = kT × ln(1/|λ|)"""
        anisotropy = np.clip(anisotropy, 0.01, 1.0)
        return float(self.k_B_eV * temperature * np.log(1 / anisotropy))
    
    def calculate_surface_area(self, porosity: float, particle_size: float = 10.0) -> float:
        """Calculate specific surface area (BET)"""
        density = self.rho_bulk * (1 - porosity)
        if density <= 0:
            density = 0.1
        
        d_m = particle_size * 1e-9
        rho_kg_m3 = density * 1000
        
        ssa = 6 / (rho_kg_m3 * d_m)
        ssa = ssa / 1000
        
        return float(np.clip(ssa, 0, 2000))
    
    def calculate_band_gap(self, porosity: float, base_band_gap: float = 3.2) -> float:
        """Estimate band gap modification due to quantum confinement"""
        confinement_factor = 1 + 0.5 * porosity
        return float(np.clip(base_band_gap * confinement_factor, 0, 10))
    
    def calculate_all_properties(self, astrophysical_data: Dict) -> Dict:
        """Calculate all material properties from astrophysical signatures"""
        required_fields = ['fractal_dimension', 'criticality_score', 'entropy', 'anisotropy', 'turbulence_beta']
        
        for field in required_fields:
            if field not in astrophysical_data:
                raise ValueError(f"Missing required field: {field}")
        
        fractal_dim = float(astrophysical_data['fractal_dimension'])
        criticality = float(astrophysical_data['criticality_score'])
        entropy = float(astrophysical_data['entropy'])
        anisotropy = float(astrophysical_data['anisotropy'])
        turbulence = float(astrophysical_data['turbulence_beta'])
        
        properties = {
            'object_name': astrophysical_data.get('object_name', 'Unknown'),
            'input_parameters': {
                'fractal_dimension': fractal_dim,
                'criticality_score': criticality,
                'entropy': entropy,
                'anisotropy': anisotropy,
                'turbulence_beta': turbulence
            }
        }
        
        properties['porosity'] = self.calculate_porosity(fractal_dim)
        properties['density'] = self.calculate_density(fractal_dim)
        properties['thermal_conductivity'] = self.calculate_thermal_conductivity(fractal_dim)
        properties['elastic_modulus'] = self.calculate_elastic_modulus(fractal_dim)
        properties['process_time'] = self.calculate_process_time(entropy)
        properties['correlation_length'] = self.calculate_correlation_length(criticality)
        properties['activation_energy'] = self.calculate_activation_energy(anisotropy)
        properties['surface_area'] = self.calculate_surface_area(properties['porosity'])
        properties['band_gap'] = self.calculate_band_gap(properties['porosity'])
        
        # Quality score
        properties['quality_score'] = self._calculate_quality_score(properties)
        
        return properties
    
    def _calculate_quality_score(self, properties: Dict) -> float:
        """Calculate overall quality score"""
        density = properties.get('density', 0)
        porosity = properties.get('porosity', 0)
        
        density_score = max(0, 1 - abs(density - 3.5) / 5)
        porosity_score = max(0, 1 - abs(porosity - 0.35) / 0.5)
        
        score = 0.5 * density_score + 0.5 * porosity_score
        return float(np.clip(score, 0, 1))
