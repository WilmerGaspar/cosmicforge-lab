"""
CosmicForge Lab Modules
=======================
Versión 2.0

Modules:
    - physics_calculator: Cálculo de propiedades físicas
    - chemistry_engine: Generación de recetas químicas
    - file_generators: Generación de archivos para simuladores
    - visualizer: Visualizaciones de estructuras y propiedades
    - file_parser: Lectura de archivos externos (QE, CIF, LAMMPS)
"""

from .physics_calculator import PhysicsCalculator
from .chemistry_engine import ChemistryEngine
from .file_generators import FileGenerator
from .visualizer import Visualizer
from .file_parser import FileParser

__version__ = "2.0.0"
__author__ = "WilmerGaspar"

__all__ = ['PhysicsCalculator', 'ChemistryEngine', 'FileGenerator', 'Visualizer', 'FileParser']
