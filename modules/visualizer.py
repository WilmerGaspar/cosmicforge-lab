"""
Visualizer Module - Create visualizations
Versión 2.0 - Código mejorado
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from typing import Dict, Tuple
import io
from PIL import Image
import warnings

warnings.filterwarnings('ignore', category=UserWarning)


class Visualizer:
    """Generate visualizations of material structures and properties"""
    
    def __init__(self, dpi: int = 100, figsize: Tuple[int, int] = (10, 8)):
        self.dpi = dpi
        self.figsize = figsize
        try:
            plt.style.use('seaborn-v0_8-whitegrid')
        except:
            pass
    
    def visualize_crystal_lattice(self, porosity: float, density: float, 
                                  object_name: str = "Material") -> Image.Image:
        """Create 2D visualization of crystal lattice"""
        fig, ax = plt.subplots(figsize=self.figsize, dpi=self.dpi)
        
        lattice_size = int(10 * (1 - porosity) + 5)
        positions = []
        np.random.seed(hash(object_name) % (2**32))
        
        for i in range(lattice_size):
            for j in range(lattice_size):
                if np.random.random() > porosity:
                    x = i + np.random.uniform(-0.2, 0.2)
                    y = j + np.random.uniform(-0.2, 0.2)
                    positions.append((x, y))
        
        if positions:
            x_coords, y_coords = zip(*positions)
            ax.scatter(x_coords, y_coords, s=200, c='#FF6B6B', alpha=0.8, 
                      edgecolors='black', linewidth=1.5, label='Átomos')
        
        n_pores = max(1, int(lattice_size ** 2 * porosity * 0.3))
        for _ in range(n_pores):
            pore_x = np.random.uniform(0, lattice_size)
            pore_y = np.random.uniform(0, lattice_size)
            pore_r = 0.4 + 0.3 * porosity
            circle = Circle((pore_x, pore_y), pore_r, fill=True, 
                           color='lightblue', alpha=0.6, edgecolor='blue')
            ax.add_patch(circle)
        
        ax.set_xlim(-1, lattice_size + 1)
        ax.set_ylim(-1, lattice_size + 1)
        ax.set_aspect('equal')
        ax.set_xlabel('Posición X (Å)', fontsize=12)
        ax.set_ylabel('Posición Y (Å)', fontsize=12)
        ax.set_title(f'Estructura Cristalina - {object_name}\n' + 
                    f'Porosidad: {porosity:.1%} | Densidad: {density:.3f} g/cm³', 
                    fontsize=14, fontweight='bold')
        ax.legend(loc='upper right')
        
        img_buf = io.BytesIO()
        plt.savefig(img_buf, format='png', dpi=self.dpi, bbox_inches='tight', facecolor='white')
        img_buf.seek(0)
        img = Image.open(img_buf)
        plt.close(fig)
        
        return img
    
    def visualize_density_map(self, material_properties: Dict, size: int = 30) -> Image.Image:
        """Create 2D density map"""
        fig, ax = plt.subplots(figsize=self.figsize, dpi=self.dpi)
        
        density = material_properties.get('density', 2.0)
        porosity = material_properties.get('porosity', 0.3)
        object_name = material_properties.get('object_name', 'Material')
        
        np.random.seed(hash(object_name) % (2**32))
        
        x = np.linspace(0, 10, size)
        y = np.linspace(0, 10, size)
        X, Y = np.meshgrid(x, y)
        
        density_field = density * (1 + 0.3 * np.sin(X * 2) * np.cos(Y * 2) + 
                                   0.1 * np.random.randn(size, size))
        density_field = np.maximum(density_field, 0.1)
        
        im = ax.imshow(density_field, cmap='YlOrRd', origin='lower', interpolation='bilinear')
        plt.colorbar(im, ax=ax, label='Densidad local (g/cm³)', shrink=0.8)
        
        ax.set_title(f'Mapa de Densidad - {object_name}', fontsize=14, fontweight='bold')
        
        img_buf = io.BytesIO()
        plt.savefig(img_buf, format='png', dpi=self.dpi, bbox_inches='tight', facecolor='white')
        img_buf.seek(0)
        img = Image.open(img_buf)
        plt.close(fig)
        
        return img
    
    def visualize_properties_summary(self, material_properties: Dict) -> Image.Image:
        """Create summary visualization"""
        fig, axes = plt.subplots(2, 3, figsize=(16, 10), dpi=self.dpi)
        
        object_name = material_properties.get('object_name', 'Material')
        fig.suptitle(f"Resumen de Propiedades - {object_name}", fontsize=16, fontweight='bold')
        
        # Porosity
        ax = axes[0, 0]
        porosity = material_properties.get('porosity', 0.3)
        ax.bar(['Porosidad'], [porosity * 100], color='#FF6B6B')
        ax.set_ylabel('%')
        ax.set_title('Porosidad', fontweight='bold')
        
        # Density
        ax = axes[0, 1]
        density = material_properties.get('density', 2.0)
        ax.bar(['Densidad'], [density], color='#4ECDC4')
        ax.set_ylabel('g/cm³')
        ax.set_title('Densidad', fontweight='bold')
        
        # Thermal conductivity
        ax = axes[0, 2]
        k = material_properties.get('thermal_conductivity', 30)
        ax.bar(['Conductividad'], [k], color='#45B7D1')
        ax.set_ylabel('W/m·K')
        ax.set_title('Conductividad Térmica', fontweight='bold')
        
        # Elastic modulus
        ax = axes[1, 0]
        elastic = material_properties.get('elastic_modulus', 150)
        ax.bar(['Módulo Elástico'], [elastic], color='#96CEB4')
        ax.set_ylabel('GPa')
        ax.set_title('Módulo Elástico', fontweight='bold')
        
        # Surface area
        ax = axes[1, 1]
        ssa = material_properties.get('surface_area', 100)
        ax.bar(['Área Superficial'], [ssa], color='#DDA0DD')
        ax.set_ylabel('m²/g')
        ax.set_title('Área Superficial (BET)', fontweight='bold')
        
        # Band gap
        ax = axes[1, 2]
        bg = material_properties.get('band_gap', 3.0)
        ax.bar(['Band Gap'], [bg], color='#FFD93D')
        ax.set_ylabel('eV')
        ax.set_title('Band Gap', fontweight='bold')
        
        plt.tight_layout()
        
        img_buf = io.BytesIO()
        plt.savefig(img_buf, format='png', dpi=self.dpi, bbox_inches='tight', facecolor='white')
        img_buf.seek(0)
        img = Image.open(img_buf)
        plt.close(fig)
        
        return img
    
    def save_image(self, img: Image.Image, filename: str) -> str:
        """Save image to file"""
        img.save(filename, quality=95)
        return filename
