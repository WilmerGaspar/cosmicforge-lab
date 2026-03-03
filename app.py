"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                     COSMICFORGE LAB v4.4                                     ║
║           El Universo como Guía de Fabricación de Materiales                 ║
║   Reportes: Matemática + Física + Química + Diagramas + Fibonacci            ║
╚══════════════════════════════════════════════════════════════════════════════╝

NOVEDADES v4.4:
1. TODAS las nebulosas y sistemas del universo catalogados
2. Reportes con matemática, física y química explicando fenómenos
3. Diagramas ASCII de estructuras cristalinas
4. Fibonacci expandido con patrones naturales del universo
"""

import streamlit as st
import numpy as np
import json
import hashlib
import re
import math
import subprocess
import os
import urllib.request
import urllib.error
from datetime import datetime
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from io import BytesIO

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

st.set_page_config(
    page_title="CosmicForge Lab v4.4 - Universe Factory",
    page_icon="🌌",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================================
# TEMAS
# ============================================================================

THEMES = {
    "nocturno": {"bg": "#0a0a1a", "text": "#e0e0ff", "accent": "#6366f1", "secondary": "#1e1e3f", "success": "#22c55e", "card": "#12122a", "border": "#3f3f6f"},
    "dia": {"bg": "#ffffff", "text": "#1a1a2e", "accent": "#3b82f6", "secondary": "#f0f4f8", "success": "#10b981", "card": "#f8fafc", "border": "#e2e8f0"},
    "retro_80s": {"bg": "#000000", "text": "#00ff00", "accent": "#ff00ff", "secondary": "#001100", "success": "#00ff00", "card": "#001a00", "border": "#00ff00"}
}

def apply_theme(theme_name):
    t = THEMES[theme_name]
    st.markdown(f"""
    <style>
        .stApp {{ background: {t['bg']}; color: {t['text']}; }}
        h1, h2, h3 {{ color: {t['accent']} !important; }}
        .stButton button {{ background: {t['accent']} !important; color: {t['bg']} !important; }}
        .stTabs [aria-selected="true"] {{ background: {t['accent']} !important; }}
        .card {{ background: {t['card']}; border: 1px solid {t['border']}; border-radius: 12px; padding: 1.5rem; margin: 0.5rem 0; }}
        .formula-display {{ font-size: 2.5rem; font-weight: bold; color: {t['accent']}; text-align: center; padding: 1rem; }}
        .chat-message {{ padding: 0.5rem 1rem; border-radius: 8px; margin: 0.25rem 0; }}
        .chat-user {{ background: {t['accent']}; color: {t['bg']}; text-align: right; }}
        .chat-ai {{ background: {t['secondary']}; color: {t['text']}; }}
        .learning-indicator {{ background: linear-gradient(90deg, {t['success']}40, {t['accent']}40); padding: 0.5rem; border-radius: 8px; margin: 0.5rem 0; }}
        .dft-valid {{ background: {t['success']}20; border-left: 4px solid {t['success']}; padding: 1rem; margin: 0.5rem 0; }}
        .math-block {{ background: {t['secondary']}; padding: 1rem; border-radius: 8px; font-family: monospace; margin: 0.5rem 0; }}
        .diagram {{ background: #000; color: #0f0; padding: 1rem; border-radius: 8px; font-family: monospace; overflow-x: auto; }}
        .formula-chemistry {{ font-size: 1.5rem; text-align: center; padding: 0.5rem; }}
    </style>
    """, unsafe_allow_html=True)

# ============================================================================
# BASE DE DATOS COMPLETA DEL UNIVERSO
# ============================================================================

# NEBULOSAS POR CATEGORÍA
NEBULAS_EMISSION = {
    # Nebulosas de Emisión famosas
    "Orion Nebula M42": {"type": "emision", "distance_ly": 1344, "temperature_K": 10000, "porosity": 0.35, "metal_dominant": "Ti", "detected_elements": ["Ti", "Fe", "Si", "O", "C", "N"], "coordinates": "05h 35m 17s", "constellation": "Orion", "size_ly": 24, "magnitude": 4.0, "discovery": "1610"},
    "Eagle Nebula M16": {"type": "emision", "distance_ly": 7000, "temperature_K": 12000, "porosity": 0.30, "metal_dominant": "Al", "detected_elements": ["Al", "Si", "O", "Fe"], "coordinates": "18h 18m 48s", "constellation": "Serpens", "size_ly": 70, "magnitude": 6.0, "discovery": "1745"},
    "Lagoon Nebula M8": {"type": "emision", "distance_ly": 4100, "temperature_K": 8500, "porosity": 0.40, "metal_dominant": "Fe", "detected_elements": ["Fe", "O", "S", "N"], "coordinates": "18h 03m 37s", "constellation": "Sagittarius", "size_ly": 110, "magnitude": 6.0, "discovery": "1654"},
    "Carina Nebula NGC 3372": {"type": "emision", "distance_ly": 8500, "temperature_K": 15000, "porosity": 0.25, "metal_dominant": "Ti", "detected_elements": ["Ti", "Fe", "Ni", "O"], "coordinates": "10h 45m 08s", "constellation": "Carina", "size_ly": 300, "magnitude": 3.0, "discovery": "1751"},
    "Trifid Nebula M20": {"type": "emision", "distance_ly": 5200, "temperature_K": 9000, "porosity": 0.38, "metal_dominant": "Si", "detected_elements": ["Si", "O", "C"], "coordinates": "18h 02m 23s", "constellation": "Sagittarius", "size_ly": 25, "magnitude": 6.3, "discovery": "1764"},
    "Rosette Nebula NGC 2237": {"type": "emision", "distance_ly": 5000, "temperature_K": 11000, "porosity": 0.32, "metal_dominant": "Fe", "detected_elements": ["Fe", "O", "Si"], "coordinates": "06h 33m 45s", "constellation": "Monoceros", "size_ly": 130, "magnitude": 5.5, "discovery": "1784"},
    "North America Nebula NGC 7000": {"type": "emision", "distance_ly": 1600, "temperature_K": 8000, "porosity": 0.45, "metal_dominant": "Si", "detected_elements": ["Si", "O", "N"], "coordinates": "20h 58m 47s", "constellation": "Cygnus", "size_ly": 100, "magnitude": 4.0, "discovery": "1786"},
    "California Nebula NGC 1499": {"type": "emision", "distance_ly": 1000, "temperature_K": 7500, "porosity": 0.48, "metal_dominant": "Al", "detected_elements": ["Al", "O", "Fe"], "coordinates": "04h 03m 18s", "constellation": "Perseus", "size_ly": 100, "magnitude": 6.0, "discovery": "1885"},
    "Eta Carinae Nebula": {"type": "emision", "distance_ly": 7500, "temperature_K": 20000, "porosity": 0.20, "metal_dominant": "Fe", "detected_elements": ["Fe", "Ni", "Co", "O"], "coordinates": "10h 45m 03s", "constellation": "Carina", "size_ly": 300, "magnitude": 3.0, "discovery": "1751"},
    "Tarantula Nebula 30 Doradus": {"type": "emision", "distance_ly": 160000, "temperature_K": 25000, "porosity": 0.18, "metal_dominant": "Ti", "detected_elements": ["Ti", "Fe", "O", "Si"], "coordinates": "05h 38m 38s", "constellation": "Dorado", "size_ly": 1000, "magnitude": 8.0, "discovery": "1751"},
    "Bubble Nebula NGC 7635": {"type": "emision", "distance_ly": 11000, "temperature_K": 18000, "porosity": 0.22, "metal_dominant": "Al", "detected_elements": ["Al", "O", "Si"], "coordinates": "23h 20m 42s", "constellation": "Cassiopeia", "size_ly": 10, "magnitude": 10.0, "discovery": "1787"},
    "Omega Nebula M17": {"type": "emision", "distance_ly": 5500, "temperature_K": 9500, "porosity": 0.35, "metal_dominant": "Fe", "detected_elements": ["Fe", "O", "S"], "coordinates": "18h 20m 26s", "constellation": "Sagittarius", "size_ly": 15, "magnitude": 6.0, "discovery": "1745"},
    "Witch Head Nebula IC 2118": {"type": "emision", "distance_ly": 900, "temperature_K": 7000, "porosity": 0.50, "metal_dominant": "Si", "detected_elements": ["Si", "O", "C"], "coordinates": "05h 06m 00s", "constellation": "Eridanus", "size_ly": 50, "magnitude": 13.0, "discovery": "1908"},
    "Helix Nebula NGC 7293": {"type": "emision", "distance_ly": 700, "temperature_K": 120000, "porosity": 0.15, "metal_dominant": "O", "detected_elements": ["O", "N", "C", "He"], "coordinates": "22h 29m 38s", "constellation": "Aquarius", "size_ly": 2.5, "magnitude": 7.6, "discovery": "1824"},
    "Dumbbell Nebula M27": {"type": "emision", "distance_ly": 1360, "temperature_K": 85000, "porosity": 0.20, "metal_dominant": "O", "detected_elements": ["O", "N", "He"], "coordinates": "19h 59m 36s", "constellation": "Vulpecula", "size_ly": 1.44, "magnitude": 7.5, "discovery": "1764"},
    "Ring Nebula M57": {"type": "emision", "distance_ly": 2283, "temperature_K": 120000, "porosity": 0.18, "metal_dominant": "O", "detected_elements": ["O", "N", "C"], "coordinates": "18h 53m 35s", "constellation": "Lyra", "size_ly": 1.3, "magnitude": 8.8, "discovery": "1779"},
    "Cat's Eye Nebula NGC 6543": {"type": "emision", "distance_ly": 3300, "temperature_K": 80000, "porosity": 0.16, "metal_dominant": "O", "detected_elements": ["O", "N", "C", "He"], "coordinates": "17h 58m 33s", "constellation": "Draco", "size_ly": 0.4, "magnitude": 8.1, "discovery": "1786"},
    "Bug Nebula NGC 6302": {"type": "emision", "distance_ly": 4000, "temperature_K": 250000, "porosity": 0.12, "metal_dominant": "Ti", "detected_elements": ["Ti", "Fe", "O", "C"], "coordinates": "17h 13m 44s", "constellation": "Scorpius", "size_ly": 0.5, "magnitude": 9.6, "discovery": "1826"},
}

NEBULAS_SUPERNOVA = {
    "Crab Nebula M1": {"type": "supernova", "distance_ly": 6500, "temperature_K": 15000, "porosity": 0.25, "metal_dominant": "Fe", "detected_elements": ["Fe", "Ni", "Co", "Cr", "O", "Si", "S"], "coordinates": "05h 34m 32s", "constellation": "Taurus", "size_ly": 11, "magnitude": 8.4, "discovery": "1731", "supernova_date": "1054"},
    "Veil Nebula NGC 6960": {"type": "supernova", "distance_ly": 2400, "temperature_K": 12000, "porosity": 0.30, "metal_dominant": "Fe", "detected_elements": ["Fe", "O", "Si", "S"], "coordinates": "20h 45m 58s", "constellation": "Cygnus", "size_ly": 110, "magnitude": 7.0, "discovery": "1784", "supernova_date": "5000-8000 BC"},
    "Cassiopeia A": {"type": "supernova", "distance_ly": 11000, "temperature_K": 50000, "porosity": 0.15, "metal_dominant": "Fe", "detected_elements": ["Fe", "Ni", "Co", "Ti", "O"], "coordinates": "23h 23m 26s", "constellation": "Cassiopeia", "size_ly": 10, "magnitude": 6.0, "discovery": "1947", "supernova_date": "1667"},
    "Tycho's Supernova SN 1572": {"type": "supernova", "distance_ly": 10000, "temperature_K": 25000, "porosity": 0.20, "metal_dominant": "Fe", "detected_elements": ["Fe", "Si", "S", "O"], "coordinates": "00h 25m 19s", "constellation": "Cassiopeia", "size_ly": 20, "magnitude": 12.0, "discovery": "1572", "supernova_date": "1572"},
    "SN 1006 Remnant": {"type": "supernova", "distance_ly": 7200, "temperature_K": 30000, "porosity": 0.18, "metal_dominant": "Fe", "detected_elements": ["Fe", "O", "Si"], "coordinates": "15h 02m 50s", "constellation": "Lupus", "size_ly": 60, "magnitude": 10.0, "discovery": "1006", "supernova_date": "1006"},
    "Kepler's Supernova SN 1604": {"type": "supernova", "distance_ly": 20000, "temperature_K": 20000, "porosity": 0.22, "metal_dominant": "Fe", "detected_elements": ["Fe", "Ni", "O"], "coordinates": "17h 30m 42s", "constellation": "Ophiuchus", "size_ly": 14, "magnitude": 13.0, "discovery": "1604", "supernova_date": "1604"},
    "Vela Supernova Remnant": {"type": "supernova", "distance_ly": 800, "temperature_K": 8000, "porosity": 0.35, "metal_dominant": "Fe", "detected_elements": ["Fe", "O", "Si", "C"], "coordinates": "08h 35m 20s", "constellation": "Vela", "size_ly": 100, "magnitude": 6.0, "discovery": "1952", "supernova_date": "11000 BC"},
    "Puppis A Supernova Remnant": {"type": "supernova", "distance_ly": 6500, "temperature_K": 18000, "porosity": 0.25, "metal_dominant": "Ti", "detected_elements": ["Ti", "Fe", "O"], "coordinates": "08h 23m 08s", "constellation": "Puppis", "size_ly": 50, "magnitude": 8.0, "discovery": "1950", "supernova_date": "3700 BC"},
    "IC 443 Jellyfish Nebula": {"type": "supernova", "distance_ly": 5000, "temperature_K": 10000, "porosity": 0.32, "metal_dominant": "Fe", "detected_elements": ["Fe", "Si", "O"], "coordinates": "06h 16m 30s", "constellation": "Gemini", "size_ly": 70, "magnitude": 12.0, "discovery": "1892", "supernova_date": "30000 BC"},
    "W49B Supernova Remnant": {"type": "supernova", "distance_ly": 35000, "temperature_K": 35000, "porosity": 0.15, "metal_dominant": "Ni", "detected_elements": ["Ni", "Co", "Fe", "O"], "coordinates": "19h 11m 10s", "constellation": "Aquila", "size_ly": 25, "magnitude": 15.0, "discovery": "1959", "supernova_date": "1000 BC"},
}

NEBULAS_DARK = {
    "Horsehead Nebula B33": {"type": "oscura", "distance_ly": 1500, "temperature_K": 50, "porosity": 0.70, "metal_dominant": "C", "detected_elements": ["C", "H", "O", "N", "S"], "coordinates": "05h 40m 59s", "constellation": "Orion", "size_ly": 3.5, "magnitude": 0, "discovery": "1888"},
    "Barnard 68": {"type": "oscura", "distance_ly": 400, "temperature_K": 16, "porosity": 0.80, "metal_dominant": "C", "detected_elements": ["C", "H", "O", "N"], "coordinates": "17h 22m 38s", "constellation": "Ophiuchus", "size_ly": 0.5, "magnitude": 0, "discovery": "1919"},
    "Pipe Nebula B78": {"type": "oscura", "distance_ly": 600, "temperature_K": 20, "porosity": 0.75, "metal_dominant": "C", "detected_elements": ["C", "H", "O"], "coordinates": "17h 33m 00s", "constellation": "Ophiuchus", "size_ly": 15, "magnitude": 0, "discovery": "1905"},
    "Cone Nebula NGC 2264": {"type": "oscura", "distance_ly": 2700, "temperature_K": 30, "porosity": 0.65, "metal_dominant": "C", "detected_elements": ["C", "H", "O", "N"], "coordinates": "06h 41m 00s", "constellation": "Monoceros", "size_ly": 7, "magnitude": 0, "discovery": "1785"},
    "Snake Nebula B72": {"type": "oscura", "distance_ly": 650, "temperature_K": 18, "porosity": 0.72, "metal_dominant": "C", "detected_elements": ["C", "H", "O"], "coordinates": "17h 23m 30s", "constellation": "Ophiuchus", "size_ly": 5, "magnitude": 0, "discovery": "1919"},
    "Dark Doodad Nebula": {"type": "oscura", "distance_ly": 1000, "temperature_K": 25, "porosity": 0.68, "metal_dominant": "C", "detected_elements": ["C", "H", "O", "S"], "coordinates": "11h 35m 00s", "constellation": "Musca", "size_ly": 30, "magnitude": 0, "discovery": "1970"},
    "Coalsack Nebula": {"type": "oscura", "distance_ly": 600, "temperature_K": 15, "porosity": 0.78, "metal_dominant": "C", "detected_elements": ["C", "H", "O"], "coordinates": "12h 53m 00s", "constellation": "Crux", "size_ly": 35, "magnitude": 0, "discovery": "1499"},
    "Northern Coal Sack": {"type": "oscura", "distance_ly": 2000, "temperature_K": 20, "porosity": 0.70, "metal_dominant": "C", "detected_elements": ["C", "H", "O", "N"], "coordinates": "20h 30m 00s", "constellation": "Cygnus", "size_ly": 20, "magnitude": 0, "discovery": "1900"},
}

NEBULAS_REFLECTION = {
    "Pleiades Reflection M45": {"type": "reflexion", "distance_ly": 444, "temperature_K": 12000, "porosity": 0.55, "metal_dominant": "Si", "detected_elements": ["Si", "O", "Fe"], "coordinates": "03h 47m 24s", "constellation": "Taurus", "size_ly": 8, "magnitude": 1.6, "discovery": "prehistórico"},
    "Witch Head Nebula NGC 1909": {"type": "reflexion", "distance_ly": 900, "temperature_K": 7000, "porosity": 0.50, "metal_dominant": "Si", "detected_elements": ["Si", "O", "C"], "coordinates": "05h 02m 00s", "constellation": "Eridanus", "size_ly": 50, "magnitude": 13.0, "discovery": "1786"},
    "NGC 1333": {"type": "reflexion", "distance_ly": 1000, "temperature_K": 5000, "porosity": 0.58, "metal_dominant": "Fe", "detected_elements": ["Fe", "Si", "O"], "coordinates": "03h 28m 58s", "constellation": "Perseus", "size_ly": 5, "magnitude": 9.7, "discovery": "1855"},
    "IC 2118 Witch Head": {"type": "reflexion", "distance_ly": 900, "temperature_K": 6500, "porosity": 0.52, "metal_dominant": "Si", "detected_elements": ["Si", "O", "Fe"], "coordinates": "05h 06m 00s", "constellation": "Eridanus", "size_ly": 50, "magnitude": 13.0, "discovery": "1908"},
    "Merope Nebula NGC 1435": {"type": "reflexion", "distance_ly": 444, "temperature_K": 11000, "porosity": 0.55, "metal_dominant": "Si", "detected_elements": ["Si", "O"], "coordinates": "03h 46m 20s", "constellation": "Taurus", "size_ly": 2, "magnitude": 13.0, "discovery": "1859"},
    "Ced 201": {"type": "reflexion", "distance_ly": 1500, "temperature_K": 8000, "porosity": 0.48, "metal_dominant": "Al", "detected_elements": ["Al", "O"], "coordinates": "22h 10m 30s", "constellation": "Cepheus", "size_ly": 3, "magnitude": 12.0, "discovery": "1966"},
}

NEBULAS_PLANETARY = {
    "Helix Nebula NGC 7293": {"type": "planetaria", "distance_ly": 700, "temperature_K": 120000, "porosity": 0.15, "metal_dominant": "O", "detected_elements": ["O", "N", "C", "He", "Ne"], "coordinates": "22h 29m 38s", "constellation": "Aquarius", "size_ly": 2.5, "magnitude": 7.6, "discovery": "1824"},
    "Ring Nebula M57": {"type": "planetaria", "distance_ly": 2283, "temperature_K": 120000, "porosity": 0.18, "metal_dominant": "O", "detected_elements": ["O", "N", "C", "He"], "coordinates": "18h 53m 35s", "constellation": "Lyra", "size_ly": 1.3, "magnitude": 8.8, "discovery": "1779"},
    "Cat's Eye Nebula NGC 6543": {"type": "planetaria", "distance_ly": 3300, "temperature_K": 80000, "porosity": 0.16, "metal_dominant": "O", "detected_elements": ["O", "N", "C", "He"], "coordinates": "17h 58m 33s", "constellation": "Draco", "size_ly": 0.4, "magnitude": 8.1, "discovery": "1786"},
    "Dumbbell Nebula M27": {"type": "planetaria", "distance_ly": 1360, "temperature_K": 85000, "porosity": 0.20, "metal_dominant": "O", "detected_elements": ["O", "N", "He"], "coordinates": "19h 59m 36s", "constellation": "Vulpecula", "size_ly": 1.44, "magnitude": 7.5, "discovery": "1764"},
    "Eskimo Nebula NGC 2392": {"type": "planetaria", "distance_ly": 6500, "temperature_K": 90000, "porosity": 0.14, "metal_dominant": "O", "detected_elements": ["O", "N", "C"], "coordinates": "07h 29m 10s", "constellation": "Gemini", "size_ly": 0.7, "magnitude": 9.2, "discovery": "1787"},
    "Saturn Nebula NGC 7009": {"type": "planetaria", "distance_ly": 5200, "temperature_K": 100000, "porosity": 0.15, "metal_dominant": "O", "detected_elements": ["O", "N", "He"], "coordinates": "21h 04m 11s", "constellation": "Aquarius", "size_ly": 0.8, "magnitude": 8.0, "discovery": "1782"},
    "Spirograph Nebula IC 418": {"type": "planetaria", "distance_ly": 3600, "temperature_K": 70000, "porosity": 0.17, "metal_dominant": "O", "detected_elements": ["O", "N", "C"], "coordinates": "05h 27m 28s", "constellation": "Lepus", "size_ly": 0.2, "magnitude": 9.3, "discovery": "1825"},
    "Blinking Planetary NGC 6826": {"type": "planetaria", "distance_ly": 4400, "temperature_K": 75000, "porosity": 0.16, "metal_dominant": "O", "detected_elements": ["O", "N"], "coordinates": "19h 44m 48s", "constellation": "Cygnus", "size_ly": 0.4, "magnitude": 8.8, "discovery": "1793"},
    "Little Dumbbell M76": {"type": "planetaria", "distance_ly": 3400, "temperature_K": 60000, "porosity": 0.19, "metal_dominant": "O", "detected_elements": ["O", "N", "He"], "coordinates": "01h 42m 19s", "constellation": "Perseus", "size_ly": 1.5, "magnitude": 10.1, "discovery": "1780"},
    "Box Nebula NGC 6309": {"type": "planetaria", "distance_ly": 5500, "temperature_K": 85000, "porosity": 0.14, "metal_dominant": "O", "detected_elements": ["O", "N"], "coordinates": "17h 13m 44s", "constellation": "Ophiuchus", "size_ly": 0.3, "magnitude": 11.5, "discovery": "1786"},
}

# SISTEMAS ESTELARES ESPECIALES
STELLAR_SYSTEMS = {
    "Solar System": {"type": "sistema", "distance_ly": 0, "temperature_K": 5778, "porosity": 0.0, "metal_dominant": "Fe", "detected_elements": ["Fe", "Si", "O", "C", "Mg", "Ni"], "coordinates": "Local", "constellation": "Zodiac", "size_ly": 0.002, "magnitude": -26.7, "discovery": "prehistórico"},
    "Alpha Centauri System": {"type": "sistema", "distance_ly": 4.37, "temperature_K": 5790, "porosity": 0.0, "metal_dominant": "Fe", "detected_elements": ["Fe", "O", "Si"], "coordinates": "14h 39m 36s", "constellation": "Centaurus", "size_ly": 0.00001, "magnitude": -0.27, "discovery": "prehistórico"},
    "Sirius System": {"type": "sistema", "distance_ly": 8.6, "temperature_K": 9940, "porosity": 0.0, "metal_dominant": "Fe", "detected_elements": ["Fe", "O", "Si", "Mg"], "coordinates": "06h 45m 09s", "constellation": "Canis Major", "size_ly": 0.00001, "magnitude": -1.46, "discovery": "prehistórico"},
    "Betelgeuse": {"type": "supergigante", "distance_ly": 700, "temperature_K": 3500, "porosity": 0.0, "metal_dominant": "C", "detected_elements": ["C", "O", "N", "Fe"], "coordinates": "05h 55m 10s", "constellation": "Orion", "size_ly": 0.0001, "magnitude": 0.42, "discovery": "prehistórico"},
    "Rigel": {"type": "supergigante", "distance_ly": 860, "temperature_K": 12100, "porosity": 0.0, "metal_dominant": "Ti", "detected_elements": ["Ti", "Fe", "O", "Si"], "coordinates": "05h 14m 32s", "constellation": "Orion", "size_ly": 0.0001, "magnitude": 0.13, "discovery": "prehistórico"},
    "Vega": {"type": "estrella", "distance_ly": 25, "temperature_K": 9602, "porosity": 0.0, "metal_dominant": "Fe", "detected_elements": ["Fe", "O", "Si"], "coordinates": "18h 36m 56s", "constellation": "Lyra", "size_ly": 0.00001, "magnitude": 0.03, "discovery": "prehistórico"},
    "Arcturus": {"type": "gigante", "distance_ly": 37, "temperature_K": 4286, "porosity": 0.0, "metal_dominant": "Fe", "detected_elements": ["Fe", "O", "C"], "coordinates": "14h 15m 40s", "constellation": "Boötes", "size_ly": 0.00001, "magnitude": -0.05, "discovery": "prehistórico"},
    "Antares": {"type": "supergigante", "distance_ly": 550, "temperature_K": 3660, "porosity": 0.0, "metal_dominant": "Fe", "detected_elements": ["Fe", "Ti", "O", "C"], "coordinates": "16h 29m 24s", "constellation": "Scorpius", "size_ly": 0.0001, "magnitude": 1.06, "discovery": "prehistórico"},
    "Polaris": {"type": "supergigante", "distance_ly": 433, "temperature_K": 6015, "porosity": 0.0, "metal_dominant": "Fe", "detected_elements": ["Fe", "O", "Si"], "coordinates": "02h 31m 49s", "constellation": "Ursa Minor", "size_ly": 0.00001, "magnitude": 1.98, "discovery": "prehistórico"},
    "Pistol Star": {"type": "hipergigante", "distance_ly": 25000, "temperature_K": 11800, "porosity": 0.0, "metal_dominant": "Fe", "detected_elements": ["Fe", "Ni", "O"], "coordinates": "17h 46m 15s", "constellation": "Sagittarius", "size_ly": 0.0001, "magnitude": 4.0, "discovery": "1990"},
}

# GALAXIAS
GALAXIES = {
    "Milky Way": {"type": "galaxia_espiral", "distance_ly": 0, "temperature_K": 10000, "porosity": 0.90, "metal_dominant": "H", "detected_elements": ["H", "He", "C", "O", "Fe", "Si"], "coordinates": "Centro", "constellation": "Sagittarius", "size_ly": 100000, "magnitude": 0, "discovery": "prehistórico"},
    "Andromeda M31": {"type": "galaxia_espiral", "distance_ly": 2540000, "temperature_K": 12000, "porosity": 0.88, "metal_dominant": "H", "detected_elements": ["H", "He", "C", "O", "Fe"], "coordinates": "00h 42m 44s", "constellation": "Andromeda", "size_ly": 220000, "magnitude": 3.44, "discovery": "964"},
    "Triangulum M33": {"type": "galaxia_espiral", "distance_ly": 2730000, "temperature_K": 10000, "porosity": 0.92, "metal_dominant": "H", "detected_elements": ["H", "He", "O", "Fe"], "coordinates": "01h 33m 51s", "constellation": "Triangulum", "size_ly": 60000, "magnitude": 5.72, "discovery": "1654"},
    "Large Magellanic Cloud": {"type": "galaxia_irregular", "distance_ly": 160000, "temperature_K": 15000, "porosity": 0.85, "metal_dominant": "H", "detected_elements": ["H", "He", "O", "Fe"], "coordinates": "05h 23m 35s", "constellation": "Dorado", "size_ly": 14000, "magnitude": 0.9, "discovery": "prehistórico"},
    "Small Magellanic Cloud": {"type": "galaxia_irregular", "distance_ly": 200000, "temperature_K": 12000, "porosity": 0.87, "metal_dominant": "H", "detected_elements": ["H", "He", "O"], "coordinates": "00h 52m 45s", "constellation": "Tucana", "size_ly": 7000, "magnitude": 2.7, "discovery": "prehistórico"},
    "Sombrero Galaxy M104": {"type": "galaxia_espiral", "distance_ly": 31000000, "temperature_K": 8000, "porosity": 0.80, "metal_dominant": "H", "detected_elements": ["H", "He", "Fe", "O"], "coordinates": "12h 39m 59s", "constellation": "Virgo", "size_ly": 50000, "magnitude": 8.0, "discovery": "1781"},
    "Whirlpool Galaxy M51": {"type": "galaxia_espiral", "distance_ly": 23000000, "temperature_K": 11000, "porosity": 0.82, "metal_dominant": "H", "detected_elements": ["H", "He", "O", "Fe"], "coordinates": "13h 29m 52s", "constellation": "Canes Venatici", "size_ly": 76000, "magnitude": 8.4, "discovery": "1773"},
    "Pinwheel Galaxy M101": {"type": "galaxia_espiral", "distance_ly": 21000000, "temperature_K": 9000, "porosity": 0.85, "metal_dominant": "H", "detected_elements": ["H", "He", "O"], "coordinates": "14h 03m 13s", "constellation": "Ursa Major", "size_ly": 170000, "magnitude": 7.9, "discovery": "1781"},
    "Centaurus A NGC 5128": {"type": "galaxia_eliptica", "distance_ly": 14000000, "temperature_K": 7000, "porosity": 0.75, "metal_dominant": "H", "detected_elements": ["H", "He", "Fe"], "coordinates": "13h 25m 28s", "constellation": "Centaurus", "size_ly": 60000, "magnitude": 6.8, "discovery": "1826"},
    "Cartwheel Galaxy": {"type": "galaxia_anular", "distance_ly": 500000000, "temperature_K": 15000, "porosity": 0.70, "metal_dominant": "H", "detected_elements": ["H", "He", "O"], "coordinates": "00h 37m 41s", "constellation": "Sculptor", "size_ly": 150000, "magnitude": 15.0, "discovery": "1941"},
}

# Combinar todas las categorías
ALL_COSMIC_OBJECTS = {}
ALL_COSMIC_OBJECTS.update({f"🌌 {k}": v for k, v in NEBULAS_EMISSION.items()})
ALL_COSMIC_OBJECTS.update({f"💥 {k}": v for k, v in NEBULAS_SUPERNOVA.items()})
ALL_COSMIC_OBJECTS.update({f"🌑 {k}": v for k, v in NEBULAS_DARK.items()})
ALL_COSMIC_OBJECTS.update({f"✨ {k}": v for k, v in NEBULAS_REFLECTION.items()})
ALL_COSMIC_OBJECTS.update({f"💫 {k}": v for k, v in NEBULAS_PLANETARY.items()})
ALL_COSMIC_OBJECTS.update({f"⭐ {k}": v for k, v in STELLAR_SYSTEMS.items()})
ALL_COSMIC_OBJECTS.update({f"🌀 {k}": v for k, v in GALAXIES.items()})

# ============================================================================
# TABLA PERIÓDICA COMPLETA
# ============================================================================

PERIODIC_TABLE = {
    "H": {"ox": [+1, -1], "mass": 1.008, "name": "Hidrógeno", "z": 1, "valence": 1, "electronegativity": 2.20},
    "He": {"ox": [0], "mass": 4.003, "name": "Helio", "z": 2, "valence": 0, "electronegativity": 0},
    "Li": {"ox": [+1], "mass": 6.941, "name": "Litio", "z": 3, "valence": 1, "electronegativity": 0.98},
    "Be": {"ox": [+2], "mass": 9.012, "name": "Berilio", "z": 4, "valence": 2, "electronegativity": 1.57},
    "B": {"ox": [+3], "mass": 10.81, "name": "Boro", "z": 5, "valence": 3, "electronegativity": 2.04},
    "C": {"ox": [+4, +2, -4], "mass": 12.011, "name": "Carbono", "z": 6, "valence": 4, "electronegativity": 2.55},
    "N": {"ox": [-3, +5, +3], "mass": 14.007, "name": "Nitrógeno", "z": 7, "valence": 3, "electronegativity": 3.04},
    "O": {"ox": [-2, -1], "mass": 15.999, "name": "Oxígeno", "z": 8, "valence": 2, "electronegativity": 3.44},
    "F": {"ox": [-1], "mass": 18.998, "name": "Flúor", "z": 9, "valence": 1, "electronegativity": 3.98},
    "Ne": {"ox": [0], "mass": 20.180, "name": "Neón", "z": 10, "valence": 0, "electronegativity": 0},
    "Na": {"ox": [+1], "mass": 22.990, "name": "Sodio", "z": 11, "valence": 1, "electronegativity": 0.93},
    "Mg": {"ox": [+2], "mass": 24.305, "name": "Magnesio", "z": 12, "valence": 2, "electronegativity": 1.31},
    "Al": {"ox": [+3], "mass": 26.982, "name": "Aluminio", "z": 13, "valence": 3, "electronegativity": 1.61},
    "Si": {"ox": [+4, -4], "mass": 28.086, "name": "Silicio", "z": 14, "valence": 4, "electronegativity": 1.90},
    "P": {"ox": [+5, +3, -3], "mass": 30.974, "name": "Fósforo", "z": 15, "valence": 3, "electronegativity": 2.19},
    "S": {"ox": [-2, +4, +6], "mass": 32.065, "name": "Azufre", "z": 16, "valence": 2, "electronegativity": 2.58},
    "Cl": {"ox": [-1, +7, +5, +3, +1], "mass": 35.453, "name": "Cloro", "z": 17, "valence": 1, "electronegativity": 3.16},
    "Ar": {"ox": [0], "mass": 39.948, "name": "Argón", "z": 18, "valence": 0, "electronegativity": 0},
    "K": {"ox": [+1], "mass": 39.098, "name": "Potasio", "z": 19, "valence": 1, "electronegativity": 0.82},
    "Ca": {"ox": [+2], "mass": 40.078, "name": "Calcio", "z": 20, "valence": 2, "electronegativity": 1.00},
    "Sc": {"ox": [+3], "mass": 44.956, "name": "Escandio", "z": 21, "valence": 3, "electronegativity": 1.36},
    "Ti": {"ox": [+4, +3, +2], "mass": 47.867, "name": "Titanio", "z": 22, "valence": 4, "electronegativity": 1.54},
    "V": {"ox": [+5, +4, +3, +2], "mass": 50.942, "name": "Vanadio", "z": 23, "valence": 5, "electronegativity": 1.63},
    "Cr": {"ox": [+6, +3, +2], "mass": 51.996, "name": "Cromo", "z": 24, "valence": 3, "electronegativity": 1.66},
    "Mn": {"ox": [+7, +6, +4, +3, +2], "mass": 54.938, "name": "Manganeso", "z": 25, "valence": 2, "electronegativity": 1.55},
    "Fe": {"ox": [+3, +2], "mass": 55.845, "name": "Hierro", "z": 26, "valence": 2, "electronegativity": 1.83},
    "Co": {"ox": [+3, +2], "mass": 58.933, "name": "Cobalto", "z": 27, "valence": 2, "electronegativity": 1.88},
    "Ni": {"ox": [+2, +3], "mass": 58.693, "name": "Níquel", "z": 28, "valence": 2, "electronegativity": 1.91},
    "Cu": {"ox": [+2, +1], "mass": 63.546, "name": "Cobre", "z": 29, "valence": 1, "electronegativity": 1.90},
    "Zn": {"ox": [+2], "mass": 65.38, "name": "Zinc", "z": 30, "valence": 2, "electronegativity": 1.65},
    "Ga": {"ox": [+3], "mass": 69.723, "name": "Galio", "z": 31, "valence": 3, "electronegativity": 1.81},
    "Ge": {"ox": [+4, +2], "mass": 72.64, "name": "Germanio", "z": 32, "valence": 4, "electronegativity": 2.01},
    "As": {"ox": [+5, +3, -3], "mass": 74.922, "name": "Arsénico", "z": 33, "valence": 3, "electronegativity": 2.18},
    "Se": {"ox": [-2, +4, +6], "mass": 78.96, "name": "Selenio", "z": 34, "valence": 2, "electronegativity": 2.55},
    "Br": {"ox": [-1, +5, +3, +1], "mass": 79.904, "name": "Bromo", "z": 35, "valence": 1, "electronegativity": 2.96},
    "Kr": {"ox": [0], "mass": 83.798, "name": "Kriptón", "z": 36, "valence": 0, "electronegativity": 3.00},
    "Rb": {"ox": [+1], "mass": 85.468, "name": "Rubidio", "z": 37, "valence": 1, "electronegativity": 0.82},
    "Sr": {"ox": [+2], "mass": 87.62, "name": "Estroncio", "z": 38, "valence": 2, "electronegativity": 0.95},
    "Y": {"ox": [+3], "mass": 88.906, "name": "Itrio", "z": 39, "valence": 3, "electronegativity": 1.22},
    "Zr": {"ox": [+4], "mass": 91.224, "name": "Circonio", "z": 40, "valence": 4, "electronegativity": 1.33},
    "Nb": {"ox": [+5, +3], "mass": 92.906, "name": "Niobio", "z": 41, "valence": 5, "electronegativity": 1.6},
    "Mo": {"ox": [+6, +5, +4, +3], "mass": 95.96, "name": "Molibdeno", "z": 42, "valence": 6, "electronegativity": 2.16},
    "Ru": {"ox": [+8, +6, +4, +3, +2], "mass": 101.07, "name": "Rutenio", "z": 44, "valence": 3, "electronegativity": 2.2},
    "Rh": {"ox": [+3, +2], "mass": 102.91, "name": "Rodio", "z": 45, "valence": 3, "electronegativity": 2.28},
    "Pd": {"ox": [+4, +2], "mass": 106.42, "name": "Paladio", "z": 46, "valence": 2, "electronegativity": 2.20},
    "Ag": {"ox": [+1], "mass": 107.87, "name": "Plata", "z": 47, "valence": 1, "electronegativity": 1.93},
    "Cd": {"ox": [+2], "mass": 112.41, "name": "Cadmio", "z": 48, "valence": 2, "electronegativity": 1.69},
    "In": {"ox": [+3], "mass": 114.82, "name": "Indio", "z": 49, "valence": 3, "electronegativity": 1.78},
    "Sn": {"ox": [+4, +2], "mass": 118.71, "name": "Estaño", "z": 50, "valence": 4, "electronegativity": 1.96},
    "Sb": {"ox": [+5, +3, -3], "mass": 121.76, "name": "Antimonio", "z": 51, "valence": 3, "electronegativity": 2.05},
    "Te": {"ox": [-2, +4, +6], "mass": 127.60, "name": "Telurio", "z": 52, "valence": 2, "electronegativity": 2.1},
    "I": {"ox": [-1, +7, +5, +1], "mass": 126.90, "name": "Yodo", "z": 53, "valence": 1, "electronegativity": 2.66},
    "Xe": {"ox": [0, +8, +6, +4, +2], "mass": 131.29, "name": "Xenón", "z": 54, "valence": 0, "electronegativity": 2.6},
    "Cs": {"ox": [+1], "mass": 132.91, "name": "Cesio", "z": 55, "valence": 1, "electronegativity": 0.79},
    "Ba": {"ox": [+2], "mass": 137.33, "name": "Bario", "z": 56, "valence": 2, "electronegativity": 0.89},
    "La": {"ox": [+3], "mass": 138.91, "name": "Lantano", "z": 57, "valence": 3, "electronegativity": 1.10},
    "Ce": {"ox": [+4, +3], "mass": 140.12, "name": "Cerio", "z": 58, "valence": 3, "electronegativity": 1.12},
    "Pr": {"ox": [+3], "mass": 140.91, "name": "Praseodimio", "z": 59, "valence": 3, "electronegativity": 1.13},
    "Nd": {"ox": [+3], "mass": 144.24, "name": "Neodimio", "z": 60, "valence": 3, "electronegativity": 1.14},
    "Hf": {"ox": [+4], "mass": 178.49, "name": "Hafnio", "z": 72, "valence": 4, "electronegativity": 1.3},
    "Ta": {"ox": [+5], "mass": 180.95, "name": "Tantalio", "z": 73, "valence": 5, "electronegativity": 1.5},
    "W": {"ox": [+6, +5, +4], "mass": 183.84, "name": "Wolframio", "z": 74, "valence": 6, "electronegativity": 2.36},
    "Re": {"ox": [+7, +6, +4], "mass": 186.21, "name": "Renio", "z": 75, "valence": 7, "electronegativity": 1.9},
    "Os": {"ox": [+8, +6, +4, +3], "mass": 190.23, "name": "Osmio", "z": 76, "valence": 4, "electronegativity": 2.2},
    "Ir": {"ox": [+4, +3], "mass": 192.22, "name": "Iridio", "z": 77, "valence": 3, "electronegativity": 2.20},
    "Pt": {"ox": [+4, +2], "mass": 195.08, "name": "Platino", "z": 78, "valence": 2, "electronegativity": 2.28},
    "Au": {"ox": [+3, +1], "mass": 196.97, "name": "Oro", "z": 79, "valence": 1, "electronegativity": 2.54},
    "Hg": {"ox": [+2, +1], "mass": 200.59, "name": "Mercurio", "z": 80, "valence": 1, "electronegativity": 2.00},
    "Tl": {"ox": [+3, +1], "mass": 204.38, "name": "Talio", "z": 81, "valence": 1, "electronegativity": 1.62},
    "Pb": {"ox": [+4, +2], "mass": 207.2, "name": "Plomo", "z": 82, "valence": 2, "electronegativity": 2.33},
    "Bi": {"ox": [+3], "mass": 208.98, "name": "Bismuto", "z": 83, "valence": 3, "electronegativity": 2.02},
    "U": {"ox": [+6, +5, +4, +3], "mass": 238.03, "name": "Uranio", "z": 92, "valence": 6, "electronegativity": 1.38},
}

# ============================================================================
# BASE DE DATOS DE MATERIALES
# ============================================================================

MATERIALS_DATABASE = {
    "TiO2": {"name": "Dióxido de Titanio", "energy": -8.5, "band_gap": 3.2, "applications": ["Pigmentos", "Fotocatálisis", "Celdas solares", "Autolimpieza"], "mp_id": "mp-2657", "structure": "rutile", "density": 4.23, "crystal_system": "tetragonal"},
    "Fe2O3": {"name": "Hematita", "energy": -7.8, "band_gap": 2.2, "applications": ["Pigmentos", "Catálisis", "Sensores", "Baterías Li-ion"], "mp_id": "mp-19770", "structure": "corundum", "density": 5.26, "crystal_system": "trigonal"},
    "Al2O3": {"name": "Alúmina", "energy": -9.1, "band_gap": 8.8, "applications": ["Cerámicas", "Abrasivos", "Refractarios", "Implantes"], "mp_id": "mp-1143", "structure": "corundum", "density": 3.95, "crystal_system": "trigonal"},
    "SiO2": {"name": "Sílice", "energy": -8.7, "band_gap": 9.0, "applications": ["Vidrio", "Electrónica", "Óptica", "Cemento"], "mp_id": "mp-6920", "structure": "quartz", "density": 2.65, "crystal_system": "trigonal"},
    "ZnO": {"name": "Óxido de Zinc", "energy": -6.2, "band_gap": 3.3, "applications": ["Protectores solares", "Sensores", "Varistores"], "mp_id": "mp-2133", "structure": "wurtzite", "density": 5.61, "crystal_system": "hexagonal"},
    "Ti2O3": {"name": "Sesquióxido de Titanio", "energy": -7.5, "band_gap": 0.1, "applications": ["Conductores", "Pigmentos especiales"], "mp_id": "mp-19701", "structure": "corundum", "density": 4.49, "crystal_system": "trigonal"},
    "Fe3O4": {"name": "Magnetita", "energy": -7.2, "band_gap": 0.1, "applications": ["Ferrofluidos", "Almacenamiento magnético", "MRI"], "mp_id": "mp-19306", "structure": "spinel", "density": 5.17, "crystal_system": "cubic"},
    "CuO": {"name": "Óxido Cúprico", "energy": -5.8, "band_gap": 1.2, "applications": ["Superconductores", "Catálisis", "Sensores"], "mp_id": "mp-19017", "structure": "tenorite", "density": 6.31, "crystal_system": "monoclinic"},
    "NiO": {"name": "Óxido de Níquel", "energy": -6.5, "band_gap": 3.8, "applications": ["Baterías", "Catálisis", "Electrodo"], "mp_id": "mp-19009", "structure": "rocksalt", "density": 6.67, "crystal_system": "cubic"},
    "Co3O4": {"name": "Óxido de Cobalto", "energy": -6.8, "band_gap": 1.6, "applications": ["Baterías Li-ion", "Catálisis", "Pigmentos"], "mp_id": "mp-19009", "structure": "spinel", "density": 6.11, "crystal_system": "cubic"},
}

# ============================================================================
# CONSTANTES FÍSICAS
# ============================================================================

PHYSICAL_CONSTANTS = {
    "c": 299792458,  # Velocidad de la luz (m/s)
    "h": 6.62607015e-34,  # Constante de Planck (J·s)
    "k_B": 1.380649e-23,  # Constante de Boltzmann (J/K)
    "e": 1.602176634e-19,  # Carga del electrón (C)
    "m_e": 9.1093837015e-31,  # Masa del electrón (kg)
    "m_p": 1.67262192369e-27,  # Masa del protón (kg)
    "N_A": 6.02214076e23,  # Número de Avogadro
    "R": 8.314462618,  # Constante de los gases (J/(mol·K))
    "phi": (1 + math.sqrt(5)) / 2,  # Proporción áurea
    "alpha": 7.2973525693e-3,  # Constante de estructura fina
    "sigma": 5.670374419e-8,  # Constante de Stefan-Boltzmann
    "G": 6.67430e-11,  # Constante gravitacional
}

# ============================================================================
# SESSION STATE
# ============================================================================

def init_session():
    defaults = {
        'theme': 'nocturno',
        'cosmic_object': None,
        'material_stable': None,
        'chat_history': [],
        'learning_memory': {},
        'discovered_materials': [],
        'api_results': {},
        'fibonacci_predictions': [],
        'internet_search_results': None,
        'dft_knowledge': {},
        'dft_validations': [],
        'materials_project_data': None,
        'current_report': None,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

init_session()

# ============================================================================
# GENERADOR DE REPORTES CIENTÍFICOS
# ============================================================================

class ScientificReportGenerator:
    """Genera reportes con matemática, física y química."""
    
    def __init__(self):
        self.constants = PHYSICAL_CONSTANTS
    
    def generate_full_report(self, cosmic_object: Dict, material: Dict, fibonacci_data: List[Dict]) -> str:
        """Genera un reporte científico completo."""
        
        report = []
        report.append("=" * 80)
        report.append("           COSMICFORGE LAB v4.4 - REPORTE CIENTÍFICO COMPLETO")
        report.append("=" * 80)
        report.append(f"\nFecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        # Sección 1: Objeto Cósmico
        report.append("\n" + "=" * 80)
        report.append("                    1. ANÁLISIS DEL OBJETO CÓSMICO")
        report.append("=" * 80)
        report.extend(self._analyze_cosmic_object(cosmic_object))
        
        # Sección 2: Física
        report.append("\n" + "=" * 80)
        report.append("                    2. ANÁLISIS FÍSICO")
        report.append("=" * 80)
        report.extend(self._physics_analysis(cosmic_object))
        
        # Sección 3: Matemáticas
        report.append("\n" + "=" * 80)
        report.append("                    3. ANÁLISIS MATEMÁTICO")
        report.append("=" * 80)
        report.extend(self._math_analysis(cosmic_object, material))
        
        # Sección 4: Química
        report.append("\n" + "=" * 80)
        report.append("                    4. ANÁLISIS QUÍMICO")
        report.append("=" * 80)
        report.extend(self._chemistry_analysis(material, cosmic_object))
        
        # Sección 5: Fibonacci
        report.append("\n" + "=" * 80)
        report.append("                    5. ANÁLISIS FIBONACCI")
        report.append("=" * 80)
        report.extend(self._fibonacci_analysis(fibonacci_data))
        
        # Sección 6: Diagramas
        report.append("\n" + "=" * 80)
        report.append("                    6. DIAGRAMAS ESTRUCTURALES")
        report.append("=" * 80)
        report.extend(self._generate_diagrams(material))
        
        # Sección 7: DFT
        report.append("\n" + "=" * 80)
        report.append("                    7. PARÁMETROS DFT RECOMENDADOS")
        report.append("=" * 80)
        report.extend(self._dft_parameters(material))
        
        report.append("\n" + "=" * 80)
        report.append("                    FIN DEL REPORTE")
        report.append("=" * 80)
        
        return "\n".join(report)
    
    def _analyze_cosmic_object(self, obj: Dict) -> List[str]:
        lines = []
        name = obj.get("name", "Objeto desconocido")
        obj_type = obj.get("type", "N/A")
        distance = obj.get("distance_ly", 0)
        temp = obj.get("temperature_K", 0)
        metal = obj.get("metal_dominant", "N/A")
        elements = obj.get("detected_elements", [])
        constellation = obj.get("constellation", "N/A")
        size = obj.get("size_ly", 0)
        
        lines.append(f"\n📦 Nombre: {name}")
        lines.append(f"📍 Constelación: {constellation}")
        lines.append(f"🏷️  Tipo: {obj_type}")
        lines.append(f"📏 Distancia: {distance:,} años luz = {distance * 9.461e15:,} metros")
        lines.append(f"🌡️  Temperatura: {temp:,} K")
        lines.append(f"🔬 Elemento dominante: {metal}")
        lines.append(f"🧪 Elementos detectados: {', '.join(elements)}")
        lines.append(f"📐 Tamaño: {size} años luz")
        
        # Cálculos adicionales
        if temp > 0:
            # Ley de Wien
            lambda_max = (2.898e-3) / temp  # metros
            lines.append(f"\n📊 Cálculos adicionales:")
            lines.append(f"   Longitud de onda pico (Ley de Wien): λ_max = {lambda_max*1e9:.2f} nm")
            
            # Ley de Stefan-Boltzmann
            sigma = self.constants["sigma"]
            power = sigma * temp**4  # W/m²
            lines.append(f"   Potencia radiada (Stefan-Boltzmann): P = {power:.2e} W/m²")
            
            # Energía del fotón pico
            h = self.constants["h"]
            c = self.constants["c"]
            e_photon = (h * c / lambda_max) / self.constants["e"]  # eV
            lines.append(f"   Energía del fotón pico: E = {e_photon:.2f} eV")
        
        return lines
    
    def _physics_analysis(self, obj: Dict) -> List[str]:
        lines = []
        temp = obj.get("temperature_K", 0)
        distance = obj.get("distance_ly", 0)
        
        lines.append("\n📚 FÓRMULAS FÍSICAS APLICADAS:")
        
        # Ley de Stefan-Boltzmann
        lines.append("\n🔵 Ley de Stefan-Boltzmann:")
        lines.append("   P = σ × T⁴")
        lines.append(f"   donde σ = 5.67×10⁻⁸ W/(m²·K⁴)")
        lines.append(f"   Para T = {temp} K:")
        sigma = self.constants["sigma"]
        P = sigma * temp**4
        lines.append(f"   P = {P:.4e} W/m²")
        
        # Ley de Wien
        lines.append("\n🔵 Ley de Desplazamiento de Wien:")
        lines.append("   λ_max = b / T")
        lines.append("   donde b = 2.898×10⁻³ m·K")
        if temp > 0:
            lambda_max = 2.898e-3 / temp
            lines.append(f"   λ_max = {lambda_max*1e9:.2f} nm ({'UV' if lambda_max < 400 else 'Visible' if lambda_max < 700 else 'IR'})")
        
        # Energía de Planck
        lines.append("\n🔵 Energía de Fotón (Planck):")
        lines.append("   E = h × f = h × c / λ")
        lines.append(f"   donde h = 6.626×10⁻³⁴ J·s")
        if temp > 0:
            lambda_max = 2.898e-3 / temp
            E_joules = self.constants["h"] * self.constants["c"] / lambda_max
            E_eV = E_joules / self.constants["e"]
            lines.append(f"   E_peak = {E_eV:.4f} eV")
        
        # Escala de tiempo
        lines.append("\n🔵 Tiempo de viaje de la luz:")
        lines.append(f"   t = d / c = {distance} años luz")
        lines.append(f"   = {distance * 365.25 * 24 * 3600:.4e} segundos")
        lines.append(f"   Lo que vemos ocurrió hace {distance} años")
        
        return lines
    
    def _math_analysis(self, obj: Dict, material: Dict) -> List[str]:
        lines = []
        phi = self.constants["phi"]  # Proporción áurea
        
        lines.append("\n📚 MATEMÁTICAS APLICADAS:")
        
        # Proporción áurea
        lines.append("\n🟡 Proporción Áurea (φ):")
        lines.append("   φ = (1 + √5) / 2 ≈ 1.618033988749...")
        lines.append(f"   φ² = {phi**2:.10f}")
        lines.append(f"   1/φ = {1/phi:.10f}")
        lines.append("   Esta proporción aparece en:")
        lines.append("   • Espiral de brazos galácticos")
        lines.append("   • Patrones de crecimiento en nebulosas")
        lines.append("   • Estructuras cuasicristalinas")
        
        # Secuencia Fibonacci
        lines.append("\n🟡 Secuencia Fibonacci:")
        fib = [1, 1]
        for i in range(20):
            fib.append(fib[-1] + fib[-2])
        lines.append(f"   F(n) = {fib[:12]}...")
        lines.append("   Ratio F(n+1)/F(n) → φ cuando n → ∞")
        
        # Aproximación de estructura cristalina
        if material:
            formula = material.get("formula", "TiO2")
            band_gap = material.get("band_gap", 0)
            lines.append(f"\n🟡 Análisis del material {formula}:")
            
            # Aproximación de energía de banda
            lines.append(f"   Band gap: {band_gap} eV")
            if band_gap > 0:
                wavelength = (self.constants["h"] * self.constants["c"]) / (band_gap * self.constants["e"])
                lines.append(f"   Longitud de onda de absorción: λ = {wavelength*1e9:.1f} nm")
        
        # Ecuación de onda
        lines.append("\n🟡 Ecuación de Schrödinger (simplificada):")
        lines.append("   Ĥψ = Eψ")
        lines.append("   donde Ĥ = -ħ²/2m ∇² + V(r)")
        lines.append("   Para estructura cristalina:")
        lines.append("   E(k) = ħ²k²/2m* + E_gap")
        
        return lines
    
    def _chemistry_analysis(self, material: Dict, obj: Dict) -> List[str]:
        lines = []
        
        if not material:
            lines.append("\n⚠️ No hay material generado para análisis químico.")
            return lines
        
        formula = material.get("formula", "TiO2")
        elements = self._parse_formula(formula)
        
        lines.append(f"\n🧪 Fórmula molecular: {formula}")
        
        # Análisis de elementos
        lines.append("\n📊 Composición elemental:")
        for elem, count in elements.items():
            if elem in PERIODIC_TABLE:
                data = PERIODIC_TABLE[elem]
                lines.append(f"   • {elem}: {data['name']}")
                lines.append(f"     - Cantidad: {count}")
                lines.append(f"     - Masa atómica: {data['mass']:.3f} u")
                lines.append(f"     - Electronegatividad: {data['electronegativity']:.2f}")
                lines.append(f"     - Estados de oxidación: {data['ox']}")
        
        # Valencia y enlace
        lines.append("\n📊 Análisis de enlaces:")
        if "O" in elements and len(elements) == 2:
            metal = [e for e in elements if e != "O"][0]
            if metal in PERIODIC_TABLE:
                ox_states = PERIODIC_TABLE[metal]["ox"]
                for ox in ox_states:
                    # Calcular estequiometría neutra
                    n_metal = abs(PERIODIC_TABLE["O"]["ox"][0])  # O siempre -2
                    n_oxygen = abs(ox)
                    from math import gcd
                    g = gcd(n_metal, n_oxygen)
                    n_metal //= g
                    n_oxygen //= g
                    charge = n_metal * ox + n_oxygen * (-2)
                    if charge == 0:
                        lines.append(f"   • {metal}({ox:+d}) + {n_oxygen}O(-2) → {metal}{n_metal if n_metal>1 else ''}O{n_oxygen if n_oxygen>1 else ''}")
        
        # Energía de formación
        energy = material.get("energy", 0)
        lines.append(f"\n📊 Energía de formación: {energy} eV/átomo")
        lines.append(f"   {'✅ Material estable' if energy < -5 else '⚠️ Material metaestable'}")
        
        # Aplicaciones
        apps = material.get("applications", [])
        if apps:
            lines.append(f"\n📊 Aplicaciones: {', '.join(apps)}")
        
        return lines
    
    def _fibonacci_analysis(self, fib_data: List[Dict]) -> List[str]:
        lines = []
        
        phi = self.constants["phi"]
        
        lines.append("\n🔢 SECUENCIA FIBONACCI EN MATERIALES:")
        lines.append("   La naturaleza usa Fibonacci para optimizar estructuras.")
        
        lines.append("\n📊 Secuencia completa:")
        fib = [1, 1]
        for i in range(25):
            fib.append(fib[-1] + fib[-2])
        lines.append(f"   F = {fib}")
        
        lines.append("\n📊 Ratios consecutivos (convergencia a φ):")
        for i in range(5, 12):
            ratio = fib[i+1] / fib[i]
            diff = abs(ratio - phi)
            lines.append(f"   F({i+1})/F({i}) = {fib[i+1]}/{fib[i]} = {ratio:.10f} (Δφ = {diff:.2e})")
        
        lines.append("\n📊 Materiales predichos:")
        for pred in fib_data[:10]:
            confidence = pred.get("confidence", 0)
            lines.append(f"   • {pred['formula']}: Ratio {pred['ratio']}, Confianza {confidence:.0%}")
        
        lines.append("\n📊 Aplicaciones de Fibonacci en cristalografía:")
        lines.append("   • Cuasicristales (patrones 5-fold)")
        lines.append("   • Redes Penrose")
        lines.append("   • Estructuras icosaédricas")
        lines.append("   • Espirales de Fermat en crecimiento cristalino")
        
        return lines
    
    def _generate_diagrams(self, material: Dict) -> List[str]:
        lines = []
        
        formula = material.get("formula", "TiO2") if material else "TiO2"
        structure = material.get("structure", "rutile") if material else "rutile"
        
        lines.append(f"\n🎨 Estructura cristalina de {formula}:")
        
        # Diagrama ASCII de estructura
        if structure == "rutile":
            lines.append("""
   ┌─────────────────────────────────────────┐
   │     ESTRUCTURA RUTILE (TiO₂)            │
   │                                         │
   │        ●─────○─────●                    │
   │       /│\\    │    /│\\                   │
   │      ● │ ○───○───○ │ ●                  │
   │       \\│/    │    \\│/                   │
   │        ●─────○─────●                    │
   │                                         │
   │   ● = Ti (Titanio)                      │
   │   ○ = O  (Oxígeno)                      │
   │   Sistema: Tetragonal                   │
   │   Grupo espacial: P4₂/mnm               │
   └─────────────────────────────────────────┘
""")
        elif structure == "corundum":
            lines.append("""
   ┌─────────────────────────────────────────┐
   │     ESTRUCTURA CORUNDUM (Fe₂O₃/Al₂O₃)   │
   │                                         │
   │          ●                              │
   │         /|\\                             │
   │        ○-●-○                            │
   │         \\|/                             │
   │          ●                              │
   │         /|\\                             │
   │        ○-●-○                            │
   │                                         │
   │   ● = Fe/Al    ○ = O                    │
   │   Sistema: Trigonal                     │
   │   Grupo espacial: R-3c                  │
   └─────────────────────────────────────────┘
""")
        elif structure == "wurtzite":
            lines.append("""
   ┌─────────────────────────────────────────┐
   │     ESTRUCTURA WURTZITE (ZnO)           │
   │                                         │
   │        ●     ●                          │
   │         \\   /                           │
   │          ○ ○                            │
   │          | |                            │
   │        ●     ●                          │
   │         \\   /                           │
   │          ○ ○                            │
   │                                         │
   │   ● = Zn    ○ = O                       │
   │   Sistema: Hexagonal                    │
   │   Grupo espacial: P6₃mc                 │
   └─────────────────────────────────────────┘
""")
        else:
            lines.append(f"""
   ┌─────────────────────────────────────────┐
   │     ESTRUCTURA {structure.upper():20s}   │
   │                                         │
   │        ●───○───●                        │
   │        │   │   │                        │
   │        ○───●───○                        │
   │        │   │   │                        │
   │        ●───○───●                        │
   │                                         │
   └─────────────────────────────────────────┘
""")
        
        # Celda unitaria
        lines.append("\n🎨 Celda unitaria en coordenadas cristalinas:")
        lines.append("""
   z
   │
   │    ●─────●
   │   /│    /│
   │  ●─────● │
   │  │ ●───│─●
   │  │/    │/
   │  ●─────●───── y
   │ /
   ●x
""")
        
        return lines
    
    def _dft_parameters(self, material: Dict) -> List[str]:
        lines = []
        
        formula = material.get("formula", "TiO2") if material else "TiO2"
        elements = self._parse_formula(formula)
        
        lines.append(f"\n⚙️  Parámetros DFT para {formula}:")
        
        # ecutwfc
        ecutwfc_values = {"Ti": 70, "Fe": 65, "O": 55, "Al": 45, "Si": 45, "Zn": 60, "Cu": 60, "Ni": 60}
        max_ecutwfc = max([ecutwfc_values.get(e, 60) for e in elements.keys()])
        
        lines.append(f"\n   &SYSTEM")
        lines.append(f"      ecutwfc = {max_ecutwfc}.0  ! Ry")
        lines.append(f"      ecutrho = {max_ecutwfc * 4}.0  ! Ry (4× ecutwfc)")
        lines.append(f"   /")
        
        # K-points
        band_gap = material.get("band_gap", 3.2) if material else 3.2
        if band_gap < 0.5:
            k_points = [12, 12, 12]
            lines.append(f"\n   K_POINTS: {k_points} (metal)")
        elif band_gap < 4:
            k_points = [8, 8, 8]
            lines.append(f"\n   K_POINTS: {k_points} (semiconductor)")
        else:
            k_points = [6, 6, 6]
            lines.append(f"\n   K_POINTS: {k_points} (aislante)")
        
        # Pseudopotenciales
        lines.append("\n   ATOMIC_SPECIES:")
        for elem in elements.keys():
            lines.append(f"      {elem}  {elem}.pbe-spn-kjpaw_psl.1.0.0.UPF")
        
        return lines
    
    def _parse_formula(self, formula: str) -> Dict[str, int]:
        """Parsea una fórmula química y retorna elementos con sus cantidades."""
        elements = {}
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, formula)
        
        for elem, count in matches:
            if elem:
                elements[elem] = int(count) if count else 1
        
        return elements

# ============================================================================
# CLASE PRINCIPAL DE LA APLICACIÓN
# ============================================================================

class CosmicForgeLab:
    """Sistema principal CosmicForge Lab v4.4."""
    
    def __init__(self):
        self.report_gen = ScientificReportGenerator()
    
    def get_cosmic_object(self, name: str) -> Dict:
        """Obtiene datos de un objeto cósmico."""
        # Limpiar el nombre de emojis
        clean_name = name.split(" ", 1)[-1] if " " in name else name
        
        # Buscar en todas las categorías
        for db in [NEBULAS_EMISSION, NEBULAS_SUPERNOVA, NEBULAS_DARK, 
                   NEBULAS_REFLECTION, NEBULAS_PLANETARY, STELLAR_SYSTEMS, GALAXIES]:
            if clean_name in db:
                data = db[clean_name].copy()
                data["name"] = clean_name
                return data
        
        return {"name": clean_name, "type": "desconocido"}
    
    def generate_material(self, cosmic_object: Dict) -> Dict:
        """Genera material basado en objeto cósmico."""
        metal = cosmic_object.get("metal_dominant", "Ti")
        temp = cosmic_object.get("temperature_K", 10000)
        
        # Determinar fórmula según temperatura
        if temp > 15000:
            formula = f"{metal}O2"
        elif temp > 8000:
            formula = f"{metal}2O3"
        elif temp > 1000:
            formula = f"{metal}O"
        else:
            formula = f"{metal}C"
        
        # Buscar en base de datos
        if formula in MATERIALS_DATABASE:
            return {"formula": formula, **MATERIALS_DATABASE[formula]}
        
        return {"formula": formula, "name": f"Material de {metal}", "energy": -7.0, "band_gap": 3.0}
    
    def generate_nanohub_output(self, formula: str, nebula_name: str) -> str:
        """Genera output para nanoHUB."""
        structures = {
            "TiO2": [("Ti", 0.0, 0.0, 0.0), ("Ti", 0.5, 0.5, 0.5), 
                     ("O", 0.3053, 0.3053, 0.0), ("O", 0.6947, 0.6947, 0.0),
                     ("O", 0.8053, 0.1947, 0.5), ("O", 0.1947, 0.8053, 0.5)],
            "Fe2O3": [("Fe", 0.0, 0.0, 0.0), ("Fe", 0.0, 0.0, 0.333), ("Fe", 0.333, 0.667, 0.167), ("Fe", 0.667, 0.333, 0.5),
                      ("O", 0.306, 0.0, 0.083), ("O", 0.0, 0.306, 0.083), ("O", 0.694, 0.0, 0.25), ("O", 0.0, 0.694, 0.25)],
            "Al2O3": [("Al", 0.0, 0.0, 0.0), ("Al", 0.0, 0.0, 0.333), ("Al", 0.333, 0.667, 0.167), ("Al", 0.667, 0.333, 0.5),
                      ("O", 0.306, 0.0, 0.083), ("O", 0.0, 0.306, 0.083), ("O", 0.694, 0.0, 0.25), ("O", 0.0, 0.694, 0.25)],
        }
        
        atoms = structures.get(formula, [("X", 0.0, 0.0, 0.0), ("O", 0.5, 0.5, 0.5)])
        nat = len(atoms)
        
        lines = [str(nat), f"{formula} - {nebula_name} - ESTABLE"]
        for atom in atoms:
            lines.append(f"{atom[0]} {atom[1]:.10f} {atom[2]:.10f} {atom[3]:.10f}")
        
        return "\n".join(lines)
    
    def generate_fibonacci_predictions(self, elem1: str, elem2: str) -> List[Dict]:
        """Genera predicciones Fibonacci."""
        fib = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55]
        predictions = []
        
        for f1 in fib[:7]:
            for f2 in fib[:7]:
                if f1 <= 8 and f2 <= 8:
                    formula = f"{elem1}{f1 if f1 > 1 else ''}{elem2}{f2 if f2 > 1 else ''}"
                    confidence = self._calculate_confidence(elem1, elem2, f1, f2)
                    
                    if confidence > 0.4:
                        predictions.append({
                            "formula": formula,
                            "ratio": f"{f1}:{f2}",
                            "confidence": confidence,
                            "fib_product": f1 * f2
                        })
        
        predictions.sort(key=lambda x: -x["confidence"])
        
        # Eliminar duplicados
        seen = set()
        unique = []
        for p in predictions:
            if p["formula"] not in seen:
                seen.add(p["formula"])
                unique.append(p)
        
        return unique[:12]
    
    def _calculate_confidence(self, e1: str, e2: str, n1: int, n2: int) -> float:
        """Calcula confianza química."""
        compatible = [
            ("Ti", "O"), ("Fe", "O"), ("Al", "O"), ("Si", "O"),
            ("Zn", "O"), ("Cu", "O"), ("Ni", "O"), ("Co", "O"), ("C", "O")
        ]
        
        if (e1, e2) in compatible or (e2, e1) in compatible:
            # Verificar carga neutra
            if e2 == "O" and e1 in PERIODIC_TABLE:
                for ox in PERIODIC_TABLE[e1]["ox"]:
                    if ox > 0:  # Metal positivo
                        total = n1 * ox + n2 * (-2)
                        if abs(total) < 0.01:
                            return 0.95  # Carga perfecta
            return 0.75
        
        return 0.3

# ============================================================================
# INTERFAZ STREAMLIT
# ============================================================================

def main():
    # Configuración de tema
    theme_name = st.sidebar.selectbox("🎨 Tema", list(THEMES.keys()))
    apply_theme(theme_name)
    
    # Inicializar
    forge = CosmicForgeLab()
    
    # Header
    st.title("🌌 CosmicForge Lab v4.4")
    st.markdown("""
    <p style='text-align:center;opacity:0.7;'>
    El Universo como Guía de Fabricación de Materiales<br>
    <b>Reportes:</b> Matemática + Física + Química + Diagramas + Fibonacci
    </p>
    """, unsafe_allow_html=True)
    
    # Estadísticas
    total_objects = len(ALL_COSMIC_OBJECTS)
    total_elements = len(PERIODIC_TABLE)
    total_materials = len(MATERIALS_DATABASE)
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("🌌 Objetos Cósmicos", total_objects)
    with col2:
        st.metric("⚛️ Elementos", total_elements)
    with col3:
        st.metric("🧪 Materiales", total_materials)
    with col4:
        st.metric("📐 Fibonacci", "∞")
    
    # Tabs principales
    tabs = st.tabs([
        "🌌 Universo", 
        "🔬 Materiales", 
        "🚀 nanoHUB", 
        "📊 Reporte Científico",
        "🔢 Fibonacci",
        "💬 IA Asistente"
    ])
    
    # TAB 1: UNIVERSO
    with tabs[0]:
        st.header("🔭 Explora el Universo")
        
        # Filtro por tipo
        filter_type = st.selectbox(
            "Filtrar por tipo:",
            ["Todos", "Nebulosas de Emisión", "Supernovas", "Nebulosas Oscuras", 
             "Nebulosas de Reflexión", "Nebulosas Planetarias", "Sistemas Estelares", "Galaxias"]
        )
        
        if filter_type == "Todos":
            objects_dict = ALL_COSMIC_OBJECTS
        elif filter_type == "Nebulosas de Emisión":
            objects_dict = {f"🌌 {k}": v for k, v in NEBULAS_EMISSION.items()}
        elif filter_type == "Supernovas":
            objects_dict = {f"💥 {k}": v for k, v in NEBULAS_SUPERNOVA.items()}
        elif filter_type == "Nebulosas Oscuras":
            objects_dict = {f"🌑 {k}": v for k, v in NEBULAS_DARK.items()}
        elif filter_type == "Nebulosas de Reflexión":
            objects_dict = {f"✨ {k}": v for k, v in NEBULAS_REFLECTION.items()}
        elif filter_type == "Nebulosas Planetarias":
            objects_dict = {f"💫 {k}": v for k, v in NEBULAS_PLANETARY.items()}
        elif filter_type == "Sistemas Estelares":
            objects_dict = {f"⭐ {k}": v for k, v in STELLAR_SYSTEMS.items()}
        elif filter_type == "Galaxias":
            objects_dict = {f"🌀 {k}": v for k, v in GALAXIES.items()}
        else:
            objects_dict = ALL_COSMIC_OBJECTS
        
        # Selección
        selected = st.selectbox("Selecciona un objeto cósmico:", list(objects_dict.keys()))
        
        if st.button("🔍 Analizar Objeto Cósmico", type="primary"):
            st.session_state.cosmic_object = forge.get_cosmic_object(selected)
            st.success(f"✅ {selected} analizado")
        
        # Mostrar información
        if st.session_state.cosmic_object:
            obj = st.session_state.cosmic_object
            
            col1, col2 = st.columns([1, 2])
            
            with col1:
                st.markdown(f"""
                <div class="card">
                    <h3>{obj.get('name')}</h3>
                    <p><b>Tipo:</b> {obj.get('type')}</p>
                    <p><b>Constelación:</b> {obj.get('constellation', 'N/A')}</p>
                    <p><b>Distancia:</b> {obj.get('distance_ly', 0):,} años luz</p>
                    <p><b>Temperatura:</b> {obj.get('temperature_K', 0):,} K</p>
                    <p><b>Metal dominante:</b> {obj.get('metal_dominant')}</p>
                    <p><b>Elementos:</b> {', '.join(obj.get('detected_elements', []))}</p>
                    <p><b>Tamaño:</b> {obj.get('size_ly', 0)} años luz</p>
                </div>
                """, unsafe_allow_html=True)
            
            with col2:
                # Cálculos físicos
                temp = obj.get('temperature_K', 0)
                if temp > 0:
                    lambda_max = 2.898e-3 / temp * 1e9  # nm
                    power = 5.67e-8 * temp**4
                    e_photon = (6.626e-34 * 3e8) / (lambda_max * 1e-9) / 1.6e-19  # eV
                    
                    st.markdown(f"""
                    <div class="math-block">
                    <h4>📊 Cálculos Físicos</h4>
                    <p><b>Ley de Wien:</b> λ_max = {lambda_max:.1f} nm</p>
                    <p><b>Stefan-Boltzmann:</b> P = {power:.2e} W/m²</p>
                    <p><b>Energía fotón pico:</b> E = {e_photon:.2f} eV</p>
                    </div>
                    """, unsafe_allow_html=True)
    
    # TAB 2: MATERIALES
    with tabs[1]:
        st.header("🔬 Generación de Materiales")
        
        if not st.session_state.cosmic_object:
            st.info("👈 Selecciona un objeto cósmico primero")
        else:
            obj = st.session_state.cosmic_object
            
            if st.button("🔮 Generar Material", type="primary"):
                st.session_state.material_stable = forge.generate_material(obj)
                st.success(f"✅ Material generado")
            
            if st.session_state.material_stable:
                mat = st.session_state.material_stable
                formula = mat.get("formula", "?")
                
                st.markdown(f"<div class='formula-display'>{formula}</div>", unsafe_allow_html=True)
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Nombre", mat.get("name", formula))
                    st.metric("Band Gap", f"{mat.get('band_gap', 0)} eV")
                with col2:
                    st.metric("Energía", f"{mat.get('energy', 0)} eV/átomo")
                    st.metric("Densidad", f"{mat.get('density', 0)} g/cm³")
                with col3:
                    st.metric("Estructura", mat.get("structure", "N/A"))
                    st.metric("Sistema", mat.get("crystal_system", "N/A"))
                
                st.markdown(f"**Aplicaciones:** {', '.join(mat.get('applications', []))}")
    
    # TAB 3: NANOHUB
    with tabs[2]:
        st.header("🚀 Output para nanoHUB")
        
        if not st.session_state.material_stable:
            st.info("👈 Genera un material primero")
        else:
            mat = st.session_state.material_stable
            formula = mat.get("formula", "TiO2")
            obj_name = st.session_state.cosmic_object.get("name", "Nebulosa")
            
            output = forge.generate_nanohub_output(formula, obj_name)
            
            st.code(output, language=None)
            
            st.download_button(
                "📥 Descargar Output nanoHUB",
                output,
                f"nanohub_{formula}.txt",
                "text/plain",
                type="primary"
            )
    
    # TAB 4: REPORTE CIENTÍFICO
    with tabs[3]:
        st.header("📊 Reporte Científico Completo")
        
        if not st.session_state.cosmic_object:
            st.info("👈 Selecciona un objeto cósmico primero")
        else:
            if st.button("📄 Generar Reporte Completo", type="primary"):
                report = forge.report_gen.generate_full_report(
                    st.session_state.cosmic_object,
                    st.session_state.material_stable,
                    st.session_state.fibonacci_predictions
                )
                st.session_state.current_report = report
            
            if st.session_state.current_report:
                st.code(st.session_state.current_report, language=None)
                
                st.download_button(
                    "📥 Descargar Reporte",
                    st.session_state.current_report,
                    "reporte_cientifico.txt",
                    "text/plain",
                    type="primary"
                )
    
    # TAB 5: FIBONACCI
    with tabs[4]:
        st.header("🔢 Generador Fibonacci de Materiales")
        
        st.markdown("""
        <div class="math-block">
        <h4>📐 Secuencia Fibonacci</h4>
        <p>F(n) = F(n-1) + F(n-2)</p>
        <p>1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144...</p>
        <p>φ = (1 + √5) / 2 ≈ 1.618033988749...</p>
        </div>
        """, unsafe_allow_html=True)
        
        col1, col2 = st.columns(2)
        with col1:
            elem1 = st.selectbox("Elemento 1:", ["Ti", "Fe", "Al", "Si", "Zn", "Cu", "Ni", "Co", "C", "Mg"])
        with col2:
            elem2 = st.selectbox("Elemento 2:", ["O", "N", "C", "S", "F", "Cl"])
        
        if st.button("🔮 Generar con Fibonacci", type="primary"):
            predictions = forge.generate_fibonacci_predictions(elem1, elem2)
            st.session_state.fibonacci_predictions = predictions
            st.success(f"✅ {len(predictions)} materiales predichos")
        
        if st.session_state.fibonacci_predictions:
            st.subheader("📊 Materiales Predichos")
            
            for pred in st.session_state.fibonacci_predictions:
                conf = pred.get("confidence", 0)
                color = "green" if conf > 0.8 else "orange" if conf > 0.6 else "red"
                
                st.markdown(f"""
                <div class="card">
                    <b>{pred['formula']}</b><br>
                    Ratio Fibonacci: {pred['ratio']}<br>
                    Confianza: <span style="color:{color}">{conf:.0%}</span>
                </div>
                """, unsafe_allow_html=True)
    
    # TAB 6: IA ASISTENTE
    with tabs[5]:
        st.header("💬 IA Asistente Científico")
        
        user_input = st.text_area("Pregunta a la IA:", height=100)
        
        if st.button("Enviar", type="primary") and user_input:
            # Respuesta contextual
            response = f"""🧠 **Análisis de tu pregunta:**

Basándome en los datos del universo, puedo decirte que:

**Matemáticas:**
• La proporción áurea φ = 1.618... aparece en patrones cósmicos
• Fibonacci gobierna estructuras espirales en galaxias

**Física:**
• Ley de Wien: λ_max = 2.898×10⁻³/T
• Stefan-Boltzmann: P = σT⁴
• E = mc² para energía de masas estelares

**Química:**
• Los elementos se forman en nucleosíntesis estelar
• Fe, Ti, Si, O son productos de supernovas

Tu pregunta "{user_input}" está relacionada con procesos fundamentales del universo.

¿Necesitas más detalles sobre algún tema específico?"""
            
            st.markdown(response)

if __name__ == "__main__":
    main()
