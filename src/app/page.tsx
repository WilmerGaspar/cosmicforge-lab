'use client'

import { useState, useCallback, useEffect } from 'react'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { Button } from '@/components/ui/button'
import { Label } from '@/components/ui/label'
import { Input } from '@/components/ui/input'
import { Slider } from '@/components/ui/slider'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'
import { Badge } from '@/components/ui/badge'
import { Separator } from '@/components/ui/separator'
import { Accordion, AccordionContent, AccordionItem, AccordionTrigger } from '@/components/ui/accordion'
import { Progress } from '@/components/ui/progress'
import { Switch } from '@/components/ui/switch'
import { 
  Sparkles, FlaskConical, Atom, FileText, Download, Database, 
  Thermometer, Gauge, Beaker, Shield, AlertTriangle, CheckCircle2,
  ExternalLink, Loader2, Settings, BookOpen, Satellite, Radiation,
  Sun, Moon, Zap, Activity, BarChart3, Radar, RefreshCw, Play,
  Cloud, Cpu, LineChart, Target, Trophy, Star, TrendingUp
} from 'lucide-react'

// Astrophysical Examples Database (Extended)
const ASTROPHYSICAL_EXAMPLES = {
  "Orion Nebula M42": {
    fractal_dimension: 1.654, criticality_score: 0.722, entropy: 0.019,
    anisotropy: 0.329, turbulence_beta: 2.278, lyapunov_max: -0.227,
    mode: "balanced", type: "Nebulosa de Emisión", distance: "1,344 ly",
    temperature: "10,000 K", luminosity: "10^5 L☉"
  },
  "Crab Nebula M1": {
    fractal_dimension: 1.82, criticality_score: 0.89, entropy: 0.032,
    anisotropy: 0.456, turbulence_beta: 2.45, lyapunov_max: -0.15,
    mode: "turbulent", type: "Remanente Supernova", distance: "6,500 ly",
    temperature: "1,600 K", luminosity: "10^5 L☉"
  },
  "Ring Nebula M57": {
    fractal_dimension: 1.45, criticality_score: 0.55, entropy: 0.012,
    anisotropy: 0.28, turbulence_beta: 1.95, lyapunov_max: -0.35,
    mode: "stable", type: "Nebulosa Planetaria", distance: "2,283 ly",
    temperature: "125,000 K", luminosity: "200 L☉"
  },
  "Triangulum Galaxy M33": {
    fractal_dimension: 1.93, criticality_score: 0.78, entropy: 0.00000964,
    anisotropy: 0.23, turbulence_beta: 1.84, lyapunov_max: 0.11,
    mode: "balanced", type: "Galaxia Espiral", distance: "2.73M ly",
    temperature: "10 K", luminosity: "3×10^9 L☉"
  },
  "Andromeda Galaxy M31": {
    fractal_dimension: 1.88, criticality_score: 0.85, entropy: 0.000012,
    anisotropy: 0.19, turbulence_beta: 1.92, lyapunov_max: 0.08,
    mode: "balanced", type: "Galaxia Espiral", distance: "2.537M ly",
    temperature: "10 K", luminosity: "2.6×10^10 L☉"
  },
  "Eagle Nebula M16": {
    fractal_dimension: 1.78, criticality_score: 0.81, entropy: 0.028,
    anisotropy: 0.42, turbulence_beta: 2.35, lyapunov_max: -0.18,
    mode: "turbulent", type: "Nebulosa de Emisión", distance: "7,000 ly",
    temperature: "8,000 K", luminosity: "10^6 L☉"
  },
  "Tarantula Nebula": {
    fractal_dimension: 1.91, criticality_score: 0.93, entropy: 0.045,
    anisotropy: 0.58, turbulence_beta: 2.68, lyapunov_max: 0.15,
    mode: "turbulent", type: "Región HII", distance: "160,000 ly",
    temperature: "15,000 K", luminosity: "10^7 L☉"
  },
  "Pillars of Creation": {
    fractal_dimension: 1.85, criticality_score: 0.88, entropy: 0.035,
    anisotropy: 0.52, turbulence_beta: 2.55, lyapunov_max: -0.12,
    mode: "turbulent", type: "Nube Molecular", distance: "6,500 ly",
    temperature: "10 K", luminosity: "N/A"
  },
  "Whirlpool Galaxy M51": {
    fractal_dimension: 1.95, criticality_score: 0.91, entropy: 0.000015,
    anisotropy: 0.22, turbulence_beta: 1.98, lyapunov_max: 0.12,
    mode: "turbulent", type: "Galaxia Interactuante", distance: "23M ly",
    temperature: "10 K", luminosity: "10^10 L☉"
  },
  "Helix Nebula NGC 7293": {
    fractal_dimension: 1.48, criticality_score: 0.58, entropy: 0.011,
    anisotropy: 0.31, turbulence_beta: 2.02, lyapunov_max: -0.32,
    mode: "stable", type: "Nebulosa Planetaria", distance: "655 ly",
    temperature: "100,000 K", luminosity: "100 L☉"
  },
  "Betelgeuse": {
    fractal_dimension: 1.72, criticality_score: 0.95, entropy: 0.055,
    anisotropy: 0.65, turbulence_beta: 2.85, lyapunov_max: 0.25,
    mode: "turbulent", type: "Supergigante Roja", distance: "700 ly",
    temperature: "3,500 K", luminosity: "1.2×10^5 L☉"
  },
  "Sirius Binary System": {
    fractal_dimension: 1.55, criticality_score: 0.62, entropy: 0.015,
    anisotropy: 0.38, turbulence_beta: 2.15, lyapunov_max: -0.08,
    mode: "balanced", type: "Sistema Binario", distance: "8.6 ly",
    temperature: "9,940 K", luminosity: "25.4 L☉"
  },
  "Jupiter Great Red Spot": {
    fractal_dimension: 1.68, criticality_score: 0.75, entropy: 0.022,
    anisotropy: 0.42, turbulence_beta: 2.32, lyapunov_max: -0.05,
    mode: "balanced", type: "Tormenta Atmosférica", distance: "4.2 AU",
    temperature: "110 K", luminosity: "N/A"
  },
  "Saturn Rings": {
    fractal_dimension: 1.92, criticality_score: 0.68, entropy: 0.008,
    anisotropy: 0.15, turbulence_beta: 1.45, lyapunov_max: -0.42,
    mode: "stable", type: "Sistema de Anillos", distance: "9.5 AU",
    temperature: "90 K", luminosity: "N/A"
  },
  "Vela Supernova Remnant": {
    fractal_dimension: 1.79, criticality_score: 0.87, entropy: 0.038,
    anisotropy: 0.48, turbulence_beta: 2.52, lyapunov_max: 0.08,
    mode: "turbulent", type: "Remanente Supernova", distance: "936 ly",
    temperature: "800 K", luminosity: "10^4 L☉"
  }
}

// Elements Database (Extended)
const ELEMENTS_DATABASE: Record<string, { name: string; atomic_num: number; mw: number; mp: number; bp: number; density: number; compatible: string[]; category: string }> = {
  "Ti": { name: "Titanio", atomic_num: 22, mw: 47.867, mp: 1668, bp: 3287, density: 4.506, compatible: ["Al", "V", "Fe", "Ni", "O", "N", "C"], category: "transition" },
  "Al": { name: "Aluminio", atomic_num: 13, mw: 26.982, mp: 660, bp: 2519, density: 2.70, compatible: ["Ti", "Cu", "Mg", "Si", "Zn", "O"], category: "post-transition" },
  "Fe": { name: "Hierro", atomic_num: 26, mw: 55.845, mp: 1538, bp: 2862, density: 7.874, compatible: ["Ni", "Cr", "Co", "Mn", "C", "O"], category: "transition" },
  "Zn": { name: "Zinc", atomic_num: 30, mw: 65.38, mp: 420, bp: 907, density: 7.14, compatible: ["Cu", "Al", "O", "S"], category: "transition" },
  "Cu": { name: "Cobre", atomic_num: 29, mw: 63.546, mp: 1085, bp: 2562, density: 8.96, compatible: ["Ni", "Zn", "Al", "Sn", "O"], category: "transition" },
  "Ni": { name: "Níquel", atomic_num: 28, mw: 58.693, mp: 1455, bp: 2913, density: 8.91, compatible: ["Fe", "Co", "Cu", "Cr", "Ti", "O"], category: "transition" },
  "Co": { name: "Cobalto", atomic_num: 27, mw: 58.933, mp: 1495, bp: 2927, density: 8.86, compatible: ["Ni", "Fe", "Mn", "O"], category: "transition" },
  "Mn": { name: "Manganeso", atomic_num: 25, mw: 54.938, mp: 1246, bp: 2061, density: 7.47, compatible: ["Fe", "Co", "O"], category: "transition" },
  "Ag": { name: "Plata", atomic_num: 47, mw: 107.868, mp: 962, bp: 2162, density: 10.49, compatible: ["Cu", "Au", "Pd", "O"], category: "transition" },
  "Au": { name: "Oro", atomic_num: 79, mw: 196.967, mp: 1064, bp: 2856, density: 19.30, compatible: ["Ag", "Pt", "Pd", "Cu"], category: "transition" },
  "Pt": { name: "Platino", atomic_num: 78, mw: 195.084, mp: 1768, bp: 3825, density: 21.45, compatible: ["Pd", "Au", "Rh", "Ir", "O"], category: "transition" },
  "Pd": { name: "Paladio", atomic_num: 46, mw: 106.42, mp: 1555, bp: 2963, density: 12.02, compatible: ["Pt", "Au", "Ag", "Ni", "O"], category: "transition" },
  "Cr": { name: "Cromo", atomic_num: 24, mw: 52.00, mp: 1907, bp: 2671, density: 7.19, compatible: ["Fe", "Ni", "Co", "O"], category: "transition" },
  "Mo": { name: "Molibdeno", atomic_num: 42, mw: 95.95, mp: 2623, bp: 4639, density: 10.28, compatible: ["Ti", "Fe", "Ni", "Cr", "O"], category: "transition" },
  "W": { name: "Wolframio", atomic_num: 74, mw: 183.84, mp: 3422, bp: 5930, density: 19.25, compatible: ["Ti", "Mo", "Cr", "O"], category: "transition" }
}

// Alloy Systems Database (Extended)
const ALLOY_SYSTEMS: Record<string, {
  name: string;
  ratio: string;
  melting_point: string;
  density: string;
  applications: string[];
  synthesis: string[];
  equipment: string[];
  precursors: string[];
  safety: string[];
  commercial_name?: string;
  category: string;
}> = {
  "Ti-Al": {
    name: "Titanio-Aluminio (Gamma-TiAl)",
    ratio: "Ti:Al = 1:1",
    melting_point: "1460°C",
    density: "3.76 g/cm³",
    applications: ["Turbinas aeronáuticas", "Componentes automotrices", "Implantes médicos"],
    synthesis: [
      "1. Fundir Ti puro (99.9%) en horno de arco eléctrico bajo argón",
      "2. Añadir Al puro (99.99%) en proporción 1:1 atómica",
      "3. Homogeneizar a 1400°C por 4 horas",
      "4. Enfriar lentamente (5°C/min) hasta 1000°C",
      "5. Tratamiento térmico a 900°C por 2 horas",
      "6. Enfriar al aire"
    ],
    equipment: ["Horno de arco eléctrico", "Atmósfera de argón", "Crisol de grafito", "Horno de tratamiento térmico"],
    precursors: ["Ti esponja (99.9%)", "Al lingote (99.99%)"],
    safety: ["Usar guantes refractarios", "Atmósfera inerte obligatoria", "Ventilación adecuada"],
    commercial_name: "Gamma-Met",
    category: "intermetallic"
  },
  "Ti-Al-V (TA6V)": {
    name: "Titanio TA6V (Ti-6Al-4V)",
    ratio: "Ti:Al:V = 90:6:4",
    melting_point: "1650°C",
    density: "4.43 g/cm³",
    applications: ["Aeroespacial", "Implantes ortopédicos", "Industria química"],
    synthesis: [
      "1. Preparar aleación maestra Ti-Al (60:40) previamente",
      "2. Fundir Ti esponja en horno de arco bajo vacío",
      "3. Añadir aleación maestra Ti-Al y V metálico",
      "4. Fundir y voltear 3-5 veces para homogeneizar",
      "5. Colar en molde de grafito precalentado",
      "6. Forjar a 950°C para eliminar segregación",
      "7. Tratamiento térmico: solución a 950°C + envejecimiento a 540°C"
    ],
    equipment: ["Horno de arco bajo vacío", "Crisol de cobre refrigerado", "Prensa de forja", "Horno de tratamiento"],
    precursors: ["Ti esponja grado aeronáutico", "Al lingote", "V metálico"],
    safety: ["Vacío alto (<10⁻³ mbar)", "Protección radiación UV", "Evitar contaminación con O₂, N₂, H₂"],
    commercial_name: "Ti-6Al-4V / Grade 5",
    category: "aerospace"
  },
  "Fe-Ni (Invar)": {
    name: "Invar (Fe-36Ni)",
    ratio: "Fe:Ni = 64:36",
    melting_point: "1425°C",
    density: "8.12 g/cm³",
    applications: ["Instrumentos de precisión", "Relojes", "Sellos térmicos"],
    synthesis: [
      "1. Fundir Fe electrolítico en horno de inducción",
      "2. Añadir Ni electrolítico gradualmente",
      "3. Mantener a 1500°C por 30 min bajo argón",
      "4. Colar en lingotera precalentada",
      "5. Laminar en caliente a 1000°C",
      "6. Recocer a 850°C por 1 hora",
      "7. Enfriar lentamente en horno"
    ],
    equipment: ["Horno de inducción", "Atmósfera inerte", "Laminador", "Horno de recocido"],
    precursors: ["Fe electrolítico", "Ni electrolítico (99.9%)"],
    safety: ["Evitar oxidación", "Controlar temperatura exacta"],
    commercial_name: "Invar 36",
    category: "precision"
  },
  "Cu-Ni (Cuproníquel)": {
    name: "Cuproníquel 70/30",
    ratio: "Cu:Ni = 70:30",
    melting_point: "1170°C",
    density: "8.94 g/cm³",
    applications: ["Componentes marinos", "Monedas", "Intercambiadores de calor"],
    synthesis: [
      "1. Fundir Cu electrolítico en horno de inducción",
      "2. Añadir Ni gradualmente bajo fundente bórax",
      "3. Desoxidar con Cu-P (0.02%)",
      "4. Colar a 1250°C en molde metálico",
      "5. Laminar en caliente a 800°C",
      "6. Recocer a 650°C por 2 horas"
    ],
    equipment: ["Horno de inducción", "Fundente bórax", "Laminador"],
    precursors: ["Cu cátodo (99.99%)", "Ni electrolítico"],
    safety: ["Fundente para evitar oxidación", "Ventilación"],
    commercial_name: "C70600 / CuNi 70/30",
    category: "marine"
  },
  "Al-Cu (Duraluminio)": {
    name: "Duraluminio 2024",
    ratio: "Al:Cu:Mg:Mn = 93.5:4.4:1.5:0.6",
    melting_point: "660°C",
    density: "2.78 g/cm³",
    applications: ["Estructuras aeronáuticas", "Componentes automotrices", "Herramientas"],
    synthesis: [
      "1. Fundir Al primario en crisol de grafito",
      "2. Añadir Cu, Mg y Mn como aleaciones maestras",
      "3. Degasear con cloro o nitrógeno",
      "4. Refinar con fundente salino",
      "5. Colar a 700°C en molde metálico",
      "6. Tratamiento T6: solución 495°C + temple + envejecimiento 190°C"
    ],
    equipment: ["Horno de resistencia", "Crisol grafito", "Molde metálico", "Horno de tratamiento"],
    precursors: ["Al primario", "Cu electrolítico", "Mg lingote", "Mn metálico"],
    safety: ["Evitar humedad (explosión)", "Degaseado obligatorio"],
    commercial_name: "AA2024-T6",
    category: "aerospace"
  },
  "TiO2 (Titania)": {
    name: "Óxido de Titanio (TiO₂)",
    ratio: "Ti:O = 1:2",
    melting_point: "1843°C",
    density: "4.23 g/cm³",
    applications: ["Pigmentos", "Fotocatálisis", "Celdas solares", "Recubrimientos"],
    synthesis: [
      "MÉTODO SOL-GEL:",
      "1. Disolver precursor Ti (isopropóxido de Ti) en etanol",
      "2. Añadir agua destilada lentamente con agitación",
      "3. Ajustar pH a 2 con HNO₃",
      "4. Envejecer gel 24 horas",
      "5. Secar a 100°C por 12 horas",
      "6. Calcinación: 400°C (anatasa) o 800°C (rutilo)"
    ],
    equipment: ["Matraz reacción", "Agitador magnético", "Horno mufla", "Autoclave (opcional)"],
    precursors: ["Isopropóxido de Ti", "TiCl₄", "Etanol", "HNO₃"],
    safety: ["Guantes y gafas", "Campana extractora", "TiCl₄ es corrosivo"],
    category: "ceramic"
  },
  "Ni-Ti (Nitinol)": {
    name: "Nitinol (NiTi - Memoria de Forma)",
    ratio: "Ni:Ti = 50:50",
    melting_point: "1310°C",
    density: "6.45 g/cm³",
    applications: ["Stents cardiovasculares", "Actuadores", "Implantes médicos", "Aparatos ortodónticos"],
    synthesis: [
      "1. Fundir Ni y Ti en horno de arco bajo argón",
      "2. Voltear lingote 5-6 veces para homogeneizar",
      "3. Forjar a 850°C",
      "4. Recocer a 500°C para ajustar temperatura de transición",
      "5. Enfriar en agua para fijar estructura",
      "6. Tratamiento térmico final para ajustar Af"
    ],
    equipment: ["Horno de arco", "Atmósfera inerte", "Horno de tratamiento", "Baño de sal"],
    precursors: ["Ni electrolítico (99.99%)", "Ti esponja (99.9%)"],
    safety: ["Control preciso de composición", "Evitar contaminación", "Atmósfera inerte"],
    commercial_name: "Nitinol SE508",
    category: "biomedical"
  },
  "Co-Cr-Mo (Vitallium)": {
    name: "Vitallium (Co-Cr-Mo)",
    ratio: "Co:Cr:Mo = 60:30:5",
    melting_point: "1490°C",
    density: "8.29 g/cm³",
    applications: ["Implantes dentales", "Prótesis articulares", "Componentes aeroespaciales"],
    synthesis: [
      "1. Fundir Co en horno de inducción bajo vacío",
      "2. Añadir Cr y Mo gradualmente",
      "3. Homogeneizar a 1500°C por 1 hora",
      "4. Colar en molde de investment",
      "5. Tratamiento térmico: solución 1200°C + envejecimiento",
      "6. Pulido y pasivación"
    ],
    equipment: ["Horno de inducción al vacío", "Moldes cerámicos", "Horno de tratamiento"],
    precursors: ["Co electrolítico", "Cr electrolítico", "Mo metálico"],
    safety: ["Vacío alto", "Protección contra humos metálicos"],
    commercial_name: "Vitallium / F75",
    category: "biomedical"
  }
}

// Commercial materials for comparison
const COMMERCIAL_MATERIALS = [
  { name: "Ti-6Al-4V (Grade 5)", strength: 950, weight: 4.43, cost: 85, thermal: 7.3, category: "Aerospace" },
  { name: "Al 7075-T6", strength: 570, weight: 2.81, cost: 35, thermal: 130, category: "Aerospace" },
  { name: "Steel 4340", strength: 1080, weight: 7.85, cost: 25, thermal: 44, category: "Structural" },
  { name: "Inconel 718", strength: 1200, weight: 8.19, cost: 120, thermal: 11, category: "High-Temp" },
  { name: "Al 2024-T3", strength: 483, weight: 2.78, cost: 30, thermal: 120, category: "Aerospace" },
  { name: "Mg AZ31", strength: 240, weight: 1.77, cost: 45, thermal: 96, category: "Lightweight" }
]

// Manin philosophical quotes
const MANIN_QUOTES = [
  "La matemática es la parte más pura de la cultura humana, y sin embargo está profundamente conectada con la realidad física.",
  "La demostrabilidad es un concepto que trasciende la lógica formal: es la búsqueda de certeza en un universo de incertidumbre.",
  "Entre lo demostrable y lo indemostrable existe una zona fronteriza donde nace la creatividad científica.",
  "La física teórica y las matemáticas comparten un misterio común: ¿por qué sus abstracciones describen tan bien la realidad?",
  "La verdad matemática tiene un estatus ontológico único: existe independientemente de nuestra capacidad para demostrarla."
]

interface AstrophysicalData {
  object_name: string;
  fractal_dimension: number;
  criticality_score: number;
  entropy: number;
  anisotropy: number;
  turbulence_beta: number;
  lyapunov_max: number;
  mode: string;
  type?: string;
  distance?: string;
  temperature?: string;
  luminosity?: string;
  source: string;
}

interface PhysicalProperties {
  porosity: number;
  density: number;
  thermal_conductivity: number;
  elastic_modulus: number;
  surface_area: number;
  band_gap: number;
  quality_score: number;
  activation_energy: number;
}

interface Recipe {
  metal: string;
  material_type: string;
  reaction_conditions: {
    temperature_C: number;
    reaction_time_hours: number;
    pH: number;
  };
  precursors: string[];
  step_by_step: string[];
}

interface SimulationResult {
  vacuum_stability: number;
  radiation_degradation: number;
  microgravity_effect: number;
  thermal_cycle_life: number;
  estimated_lifetime: string;
  recommendations: string[];
}

interface OptimizationResult {
  best_composition: string;
  score: number;
  parameters: {
    temperature: number;
    time: number;
    ph: number;
  };
  generations: number;
}

interface DFTResult {
  total_energy: number;
  band_gap: number;
  optimized_lattice: number;
  fermi_energy: number;
  is_stable: boolean;
  calculation_time: string;
}

export default function CosmicForgeLab() {
  // State
  const [selectedExample, setSelectedExample] = useState<string>("Orion Nebula M42")
  const [selectedMetal, setSelectedMetal] = useState<string>("Ti")
  const [materialType, setMaterialType] = useState<string>("Oxido")
  const [selectedAlloy, setSelectedAlloy] = useState<string>("")
  const [useAlloy, setUseAlloy] = useState(false)
  
  const [astroData, setAstroData] = useState<AstrophysicalData | null>(null)
  const [physicalProps, setPhysicalProps] = useState<PhysicalProperties | null>(null)
  const [recipe, setRecipe] = useState<Recipe | null>(null)
  const [isLoading, setIsLoading] = useState(false)
  
  // New states
  const [simulationResult, setSimulationResult] = useState<SimulationResult | null>(null)
  const [optimizationResult, setOptimizationResult] = useState<OptimizationResult | null>(null)
  const [dftResult, setDftResult] = useState<DFTResult | null>(null)
  const [materialsProjectData, setMaterialsProjectData] = useState<any>(null)
  const [activeTab, setActiveTab] = useState("properties")
  const [simulationLoading, setSimulationLoading] = useState(false)
  const [optimizationLoading, setOptimizationLoading] = useState(false)
  const [dftLoading, setDftLoading] = useState(false)
  const [dftProgress, setDftProgress] = useState(0)
  
  // NASA-style dark mode
  const [darkMode, setDarkMode] = useState(true)
  
  // Manual editor state
  const [manualName, setManualName] = useState("Objeto Personalizado")
  const [manualFd, setManualFd] = useState([1.5])
  const [manualCs, setManualCs] = useState([0.5])
  const [manualEntropy, setManualEntropy] = useState("0.01")
  const [manualAnisotropy, setManualAnisotropy] = useState([0.3])
  const [manualTurbulence, setManualTurbulence] = useState([2.0])
  const [manualLyapunov, setManualLyapunov] = useState("-0.1")
  
  const [inputMethod, setInputMethod] = useState<"examples" | "manual">("examples")

  // Load example
  const loadExample = useCallback(() => {
    const example = ASTROPHYSICAL_EXAMPLES[selectedExample as keyof typeof ASTROPHYSICAL_EXAMPLES]
    if (example) {
      setAstroData({
        ...example,
        object_name: selectedExample,
        source: 'Base de Datos'
      })
    }
  }, [selectedExample])

  // Apply manual data
  const applyManualData = useCallback(() => {
    setAstroData({
      object_name: manualName,
      fractal_dimension: manualFd[0],
      criticality_score: manualCs[0],
      entropy: parseFloat(manualEntropy) || 0.01,
      anisotropy: manualAnisotropy[0],
      turbulence_beta: manualTurbulence[0],
      lyapunov_max: parseFloat(manualLyapunov) || -0.1,
      mode: 'custom',
      source: 'Editor'
    })
  }, [manualName, manualFd, manualCs, manualEntropy, manualAnisotropy, manualTurbulence, manualLyapunov])

  // Generate recipe
  const generateRecipe = useCallback(async () => {
    if (!astroData) return
    
    setIsLoading(true)
    try {
      const response = await fetch('/api/cosmicforge', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          action: 'generate',
          astroData,
          metal: selectedMetal,
          materialType,
          alloyKey: useAlloy ? selectedAlloy : null
        })
      })
      
      const data = await response.json()
      if (data.success) {
        setPhysicalProps(data.physicalProps)
        setRecipe(data.recipe)
      }
    } catch (error) {
      console.error('Error generating recipe:', error)
    } finally {
      setIsLoading(false)
    }
  }, [astroData, selectedMetal, materialType, useAlloy, selectedAlloy])

  // Run extreme conditions simulation
  const runSimulation = useCallback(async () => {
    if (!physicalProps) return
    
    setSimulationLoading(true)
    try {
      const response = await fetch('/api/cosmicforge', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          action: 'simulate',
          physicalProps,
          metal: selectedMetal,
          materialType
        })
      })
      
      const data = await response.json()
      if (data.success) {
        setSimulationResult(data.simulation)
      }
    } catch (error) {
      console.error('Error running simulation:', error)
    } finally {
      setSimulationLoading(false)
    }
  }, [physicalProps, selectedMetal, materialType])

  // Run optimization
  const runOptimization = useCallback(async () => {
    if (!astroData) return
    
    setOptimizationLoading(true)
    try {
      const response = await fetch('/api/cosmicforge', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          action: 'optimize',
          astroData,
          metal: selectedMetal,
          materialType
        })
      })
      
      const data = await response.json()
      if (data.success) {
        setOptimizationResult(data.optimization)
      }
    } catch (error) {
      console.error('Error running optimization:', error)
    } finally {
      setOptimizationLoading(false)
    }
  }, [astroData, selectedMetal, materialType])

  // Run DFT calculation
  const runDFT = useCallback(async () => {
    if (!physicalProps || !recipe) return
    
    setDftLoading(true)
    setDftProgress(0)
    
    // Simulate progress
    const progressInterval = setInterval(() => {
      setDftProgress(prev => {
        if (prev >= 95) {
          clearInterval(progressInterval)
          return prev
        }
        return prev + Math.random() * 10
      })
    }, 500)
    
    try {
      const response = await fetch('/api/cosmicforge', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          action: 'dft',
          physicalProps,
          recipe,
          metal: selectedMetal
        })
      })
      
      const data = await response.json()
      if (data.success) {
        setDftResult(data.dft)
      }
    } catch (error) {
      console.error('Error running DFT:', error)
    } finally {
      clearInterval(progressInterval)
      setDftProgress(100)
      setTimeout(() => setDftLoading(false), 500)
    }
  }, [physicalProps, recipe, selectedMetal])

  // Query Materials Project
  const queryMaterialsProject = useCallback(async () => {
    if (!selectedMetal) return
    
    try {
      const response = await fetch('/api/cosmicforge', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          action: 'materials_project',
          metal: selectedMetal,
          formula: `${selectedMetal}O2`
        })
      })
      
      const data = await response.json()
      if (data.success) {
        setMaterialsProjectData(data.results)
      }
    } catch (error) {
      console.error('Error querying Materials Project:', error)
    }
  }, [selectedMetal])

  // Download PDF
  const downloadPDF = useCallback(async () => {
    if (!astroData || !physicalProps || !recipe) return
    
    try {
      const response = await fetch('/api/cosmicforge', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          action: 'pdf',
          astroData,
          physicalProps,
          recipe,
          alloyKey: selectedAlloy,
          simulation: simulationResult,
          optimization: optimizationResult,
          dft: dftResult
        })
      })
      
      const blob = await response.blob()
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = `CosmicForge_Report_${astroData.object_name.replace(/\s+/g, '_')}.txt`
      document.body.appendChild(a)
      a.click()
      document.body.removeChild(a)
      URL.revokeObjectURL(url)
    } catch (error) {
      console.error('Error downloading PDF:', error)
    }
  }, [astroData, physicalProps, recipe, selectedAlloy, simulationResult, optimizationResult, dftResult])

  // Clear data
  const clearData = useCallback(() => {
    setAstroData(null)
    setPhysicalProps(null)
    setRecipe(null)
    setSimulationResult(null)
    setOptimizationResult(null)
    setDftResult(null)
    setMaterialsProjectData(null)
  }, [])

  const getAlloyOptions = () => {
    return Object.keys(ALLOY_SYSTEMS).filter(key => key.includes(selectedMetal))
  }

  const currentAlloy = selectedAlloy ? ALLOY_SYSTEMS[selectedAlloy] : null
  const randomQuote = MANIN_QUOTES[Math.floor(Math.random() * MANIN_QUOTES.length)]

  // Theme classes
  const themeClasses = darkMode 
    ? 'min-h-screen bg-[#0a0a0a] text-white'
    : 'min-h-screen bg-gradient-to-br from-slate-50 via-blue-50 to-indigo-50'
    
  const cardClasses = darkMode
    ? 'shadow-sm border-[#333] bg-[#111]'
    : 'shadow-sm border-slate-200'

  return (
    <div className={themeClasses}>
      {/* Header */}
      <header className={`sticky top-0 z-50 ${darkMode ? 'bg-[#0a0a0a]/90 border-[#333]' : 'bg-white/80 border-slate-200'} backdrop-blur-sm border-b`}>
        <div className="max-w-7xl mx-auto px-4 py-4">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-3">
              <div className={`p-2 rounded-xl ${darkMode ? 'bg-gradient-to-br from-[#ff6b00] to-[#ff8c00]' : 'bg-gradient-to-br from-blue-500 to-indigo-600'}`}>
                <Atom className="w-8 h-8 text-white" />
              </div>
              <div>
                <h1 className={`text-2xl font-bold ${darkMode ? 'text-white' : 'text-slate-800'}`}>CosmicForge Lab</h1>
                <p className={`text-sm ${darkMode ? 'text-gray-400' : 'text-slate-500'}`}>Diseño de Materiales Inspirado en Firmas Astrofísicas</p>
              </div>
            </div>
            <div className="flex items-center gap-3">
              {/* Dark mode toggle */}
              <div className="flex items-center gap-2">
                <Sun className={`w-4 h-4 ${darkMode ? 'text-gray-500' : 'text-amber-500'}`} />
                <Switch checked={darkMode} onCheckedChange={setDarkMode} />
                <Moon className={`w-4 h-4 ${darkMode ? 'text-[#ff6b00]' : 'text-gray-500'}`} />
              </div>
              <Badge variant="outline" className={`text-sm px-3 py-1 ${darkMode ? 'border-[#ff6b00] text-[#ff6b00]' : ''}`}>
                v3.5 NASA Edition
              </Badge>
            </div>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-4 py-6">
        <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
          {/* Sidebar */}
          <aside className="lg:col-span-1 space-y-4">
            {/* Input Method Selection */}
            <Card className={cardClasses}>
              <CardHeader className="pb-3">
                <CardTitle className={`text-lg flex items-center gap-2 ${darkMode ? 'text-white' : ''}`}>
                  <Database className={`w-5 h-5 ${darkMode ? 'text-[#ff6b00]' : 'text-blue-500'}`} />
                  Importar Datos
                </CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <Tabs value={inputMethod} onValueChange={(v) => setInputMethod(v as "examples" | "manual")}>
                  <TabsList className="grid w-full grid-cols-2">
                    <TabsTrigger value="examples">Ejemplos</TabsTrigger>
                    <TabsTrigger value="manual">Manual</TabsTrigger>
                  </TabsList>
                  
                  <TabsContent value="examples" className="space-y-3 mt-3">
                    <Label className={darkMode ? 'text-gray-300' : ''}>Objeto Astrofísico</Label>
                    <Select value={selectedExample} onValueChange={setSelectedExample}>
                      <SelectTrigger className={darkMode ? 'bg-[#1a1a1a] border-[#333]' : ''}>
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        {Object.keys(ASTROPHYSICAL_EXAMPLES).map(name => (
                          <SelectItem key={name} value={name}>{name}</SelectItem>
                        ))}
                      </SelectContent>
                    </Select>
                    
                    {ASTROPHYSICAL_EXAMPLES[selectedExample as keyof typeof ASTROPHYSICAL_EXAMPLES] && (
                      <div className={`p-3 rounded-lg text-sm ${darkMode ? 'bg-[#1a1a1a]' : 'bg-blue-50'}`}>
                        <p><strong>Tipo:</strong> {ASTROPHYSICAL_EXAMPLES[selectedExample as keyof typeof ASTROPHYSICAL_EXAMPLES].type}</p>
                        <p><strong>Distancia:</strong> {ASTROPHYSICAL_EXAMPLES[selectedExample as keyof typeof ASTROPHYSICAL_EXAMPLES].distance}</p>
                        <p><strong>Temp:</strong> {ASTROPHYSICAL_EXAMPLES[selectedExample as keyof typeof ASTROPHYSICAL_EXAMPLES].temperature}</p>
                      </div>
                    )}
                    
                    <Button onClick={loadExample} className="w-full" style={darkMode ? { backgroundColor: '#ff6b00' } : {}}>
                      Cargar Ejemplo
                    </Button>
                  </TabsContent>
                  
                  <TabsContent value="manual" className="space-y-3 mt-3">
                    <div>
                      <Label className={darkMode ? 'text-gray-300' : ''}>Nombre del Objeto</Label>
                      <Input value={manualName} onChange={(e) => setManualName(e.target.value)} className={darkMode ? 'bg-[#1a1a1a] border-[#333]' : ''} />
                    </div>
                    
                    <div>
                      <Label className={darkMode ? 'text-gray-300' : ''}>Dim. Fractal: {manualFd[0].toFixed(3)}</Label>
                      <Slider value={manualFd} onValueChange={setManualFd} min={0} max={3} step={0.001} />
                    </div>
                    
                    <div>
                      <Label className={darkMode ? 'text-gray-300' : ''}>Criticalidad: {manualCs[0].toFixed(3)}</Label>
                      <Slider value={manualCs} onValueChange={setManualCs} min={0} max={1} step={0.001} />
                    </div>
                    
                    <div>
                      <Label className={darkMode ? 'text-gray-300' : ''}>Entropía</Label>
                      <Input type="number" value={manualEntropy} onChange={(e) => setManualEntropy(e.target.value)} step="0.000001" className={darkMode ? 'bg-[#1a1a1a] border-[#333]' : ''} />
                    </div>
                    
                    <Button onClick={applyManualData} className="w-full" style={darkMode ? { backgroundColor: '#ff6b00' } : {}}>
                      Aplicar Parámetros
                    </Button>
                  </TabsContent>
                </Tabs>
              </CardContent>
            </Card>

            {/* Material Configuration */}
            <Card className={cardClasses}>
              <CardHeader className="pb-3">
                <CardTitle className={`text-lg flex items-center gap-2 ${darkMode ? 'text-white' : ''}`}>
                  <Settings className={`w-5 h-5 ${darkMode ? 'text-[#ff6b00]' : 'text-green-500'}`} />
                  Configuración
                </CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <div>
                  <Label className={darkMode ? 'text-gray-300' : ''}>Tipo de Material</Label>
                  <Select value={materialType} onValueChange={setMaterialType}>
                    <SelectTrigger className={darkMode ? 'bg-[#1a1a1a] border-[#333]' : ''}>
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="Oxido">Óxido</SelectItem>
                      <SelectItem value="Metal">Metal</SelectItem>
                      <SelectItem value="Ceramica">Cerámica</SelectItem>
                      <SelectItem value="Nanoparticula">Nanopartícula</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
                
                <div>
                  <Label className={darkMode ? 'text-gray-300' : ''}>Metal Base</Label>
                  <Select value={selectedMetal} onValueChange={setSelectedMetal}>
                    <SelectTrigger className={darkMode ? 'bg-[#1a1a1a] border-[#333]' : ''}>
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      {Object.entries(ELEMENTS_DATABASE).map(([symbol, data]) => (
                        <SelectItem key={symbol} value={symbol}>
                          {symbol} - {data.name}
                        </SelectItem>
                      ))}
                    </SelectContent>
                  </Select>
                </div>
                
                {ELEMENTS_DATABASE[selectedMetal] && (
                  <div className={`p-3 rounded-lg text-sm ${darkMode ? 'bg-[#1a1a1a]' : 'bg-green-50'}`}>
                    <p><strong>Compatible con:</strong></p>
                    <div className="flex flex-wrap gap-1 mt-1">
                      {ELEMENTS_DATABASE[selectedMetal].compatible.map(el => (
                        <Badge key={el} variant="secondary" className="text-xs">{el}</Badge>
                      ))}
                    </div>
                  </div>
                )}
                
                <Separator className={darkMode ? 'bg-[#333]' : ''} />
                
                <div className="flex items-center space-x-2">
                  <input
                    type="checkbox"
                    id="useAlloy"
                    checked={useAlloy}
                    onChange={(e) => setUseAlloy(e.target.checked)}
                    className="h-4 w-4 rounded"
                  />
                  <Label htmlFor="useAlloy" className={darkMode ? 'text-gray-300' : ''}>Usar aleación predefinida</Label>
                </div>
                
                {useAlloy && (
                  <div>
                    <Label className={darkMode ? 'text-gray-300' : ''}>Sistema de Aleación</Label>
                    <Select value={selectedAlloy} onValueChange={setSelectedAlloy}>
                      <SelectTrigger className={darkMode ? 'bg-[#1a1a1a] border-[#333]' : ''}>
                        <SelectValue placeholder="Seleccionar..." />
                      </SelectTrigger>
                      <SelectContent>
                        {Object.keys(ALLOY_SYSTEMS).map(key => (
                          <SelectItem key={key} value={key}>
                            {ALLOY_SYSTEMS[key].name}
                          </SelectItem>
                        ))}
                      </SelectContent>
                    </Select>
                  </div>
                )}
              </CardContent>
            </Card>

            {/* Clear Button */}
            {astroData && (
              <Button onClick={clearData} variant="outline" className={`w-full ${darkMode ? 'border-[#333] text-gray-300' : ''}`}>
                Limpiar Datos
              </Button>
            )}
          </aside>

          {/* Main Content */}
          <div className="lg:col-span-3 space-y-6">
            {/* Astro Data Display */}
            {astroData ? (
              <>
                <Card className={cardClasses}>
                  <CardHeader>
                    <CardTitle className={`flex items-center gap-2 ${darkMode ? 'text-white' : ''}`}>
                      <Sparkles className={`w-5 h-5 ${darkMode ? 'text-[#ff6b00]' : 'text-amber-500'}`} />
                      Firma Astrofísica Detectada
                    </CardTitle>
                    <CardDescription className={darkMode ? 'text-gray-400' : ''}>
                      Objeto: {astroData.object_name} | Fuente: {astroData.source}
                    </CardDescription>
                  </CardHeader>
                  <CardContent>
                    <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                      {[
                        { label: 'Dimensión Fractal', value: astroData.fractal_dimension.toFixed(4) },
                        { label: 'Criticalidad', value: astroData.criticality_score.toFixed(4) },
                        { label: 'Entropía', value: astroData.entropy.toFixed(6) },
                        { label: 'Anisotropía', value: astroData.anisotropy.toFixed(4) },
                        { label: 'Turbulencia β', value: astroData.turbulence_beta.toFixed(4) },
                        { label: 'Lyapunov máx', value: astroData.lyapunov_max.toFixed(6) },
                        { label: 'Modo', value: astroData.mode, badge: true },
                        { label: 'Tipo', value: astroData.type || 'N/A' }
                      ].map((item, i) => (
                        <div key={i} className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                          <p className={`text-sm ${darkMode ? 'text-gray-400' : 'text-slate-500'}`}>{item.label}</p>
                          {item.badge ? (
                            <Badge variant={astroData.mode === 'turbulent' ? 'destructive' : astroData.mode === 'stable' ? 'default' : 'secondary'}>
                              {item.value}
                            </Badge>
                          ) : (
                            <p className={`text-xl font-bold ${darkMode ? 'text-white' : 'text-slate-800'}`}>{item.value}</p>
                          )}
                        </div>
                      ))}
                    </div>
                    
                    <Button 
                      onClick={generateRecipe} 
                      className="w-full mt-4" 
                      size="lg"
                      disabled={isLoading}
                      style={darkMode ? { backgroundColor: '#ff6b00' } : {}}
                    >
                      {isLoading ? (
                        <>
                          <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                          Calculando...
                        </>
                      ) : (
                        <>
                          <FlaskConical className="w-4 h-4 mr-2" />
                          Generar Receta de Material
                        </>
                      )}
                    </Button>
                  </CardContent>
                </Card>

                {/* Results Tabs */}
                {physicalProps && recipe && (
                  <Tabs value={activeTab} onValueChange={setActiveTab} className="space-y-4">
                    <TabsList className={`grid w-full grid-cols-7 ${darkMode ? 'bg-[#1a1a1a]' : ''}`}>
                      <TabsTrigger value="properties">Propiedades</TabsTrigger>
                      <TabsTrigger value="simulation">Simulación</TabsTrigger>
                      <TabsTrigger value="optimization">Optimización</TabsTrigger>
                      <TabsTrigger value="dft">DFT</TabsTrigger>
                      <TabsTrigger value="comparison">Comparación</TabsTrigger>
                      <TabsTrigger value="files">Archivos</TabsTrigger>
                      <TabsTrigger value="pdf">PDF</TabsTrigger>
                    </TabsList>

                    {/* Properties Tab */}
                    <TabsContent value="properties">
                      <Card className={cardClasses}>
                        <CardHeader>
                          <CardTitle className={darkMode ? 'text-white' : ''}>Propiedades del Material</CardTitle>
                        </CardHeader>
                        <CardContent>
                          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                            {[
                              { label: 'Porosidad', value: `${(physicalProps.porosity * 100).toFixed(1)}%`, color: 'blue' },
                              { label: 'Densidad', value: `${physicalProps.density.toFixed(3)} g/cm³`, color: 'green' },
                              { label: 'Cond. Térmica', value: `${physicalProps.thermal_conductivity.toFixed(2)} W/mK`, color: 'amber' },
                              { label: 'Mód. Elástico', value: `${physicalProps.elastic_modulus.toFixed(2)} GPa`, color: 'purple' },
                              { label: 'Área Superficial', value: `${physicalProps.surface_area.toFixed(1)} m²/g`, color: 'rose' },
                              { label: 'Band Gap', value: `${physicalProps.band_gap.toFixed(3)} eV`, color: 'cyan' },
                              { label: 'Calidad', value: `${(physicalProps.quality_score * 100).toFixed(1)}%`, color: 'indigo' },
                              { label: 'E. Activación', value: `${physicalProps.activation_energy.toFixed(4)} eV`, color: 'teal' }
                            ].map((item, i) => (
                              <div key={i} className={`p-4 rounded-lg ${darkMode ? `bg-[#1a1a1a] border-l-4 border-${item.color}-500` : `bg-${item.color}-50`}`}>
                                <p className={`text-sm ${darkMode ? 'text-gray-400' : ''}`}>{item.label}</p>
                                <p className={`text-2xl font-bold ${darkMode ? 'text-white' : ''}`}>{item.value}</p>
                              </div>
                            ))}
                          </div>

                          <Separator className={`my-4 ${darkMode ? 'bg-[#333]' : ''}`} />

                          <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                            <h4 className="font-semibold mb-2 flex items-center gap-2">
                              <ExternalLink className="w-4 h-4" />
                              Validación con Bases de Datos
                            </h4>
                            <div className="flex flex-wrap gap-4">
                              <a 
                                href={`https://materialsproject.org/materials?search=${selectedMetal}O2`} 
                                target="_blank" 
                                rel="noopener noreferrer"
                                className={`text-sm ${darkMode ? 'text-[#ff6b00]' : 'text-blue-600'} hover:underline`}
                              >
                                Materials Project
                              </a>
                              <a 
                                href={`https://www.aflowlib.org/?search=${selectedMetal}O2`} 
                                target="_blank" 
                                rel="noopener noreferrer"
                                className={`text-sm ${darkMode ? 'text-[#ff6b00]' : 'text-blue-600'} hover:underline`}
                              >
                                AFLOW Library
                              </a>
                              <a 
                                href={`http://oqmd.org/materials?search=${selectedMetal}O2`} 
                                target="_blank" 
                                rel="noopener noreferrer"
                                className={`text-sm ${darkMode ? 'text-[#ff6b00]' : 'text-blue-600'} hover:underline`}
                              >
                                OQMD Database
                              </a>
                            </div>
                            <Button 
                              onClick={queryMaterialsProject} 
                              variant="outline" 
                              size="sm" 
                              className={`mt-2 ${darkMode ? 'border-[#333]' : ''}`}
                            >
                              <Database className="w-3 h-3 mr-1" />
                              Consultar API Materials Project
                            </Button>
                            
                            {materialsProjectData && (
                              <div className={`mt-3 p-3 rounded-lg ${darkMode ? 'bg-[#222]' : 'bg-blue-50'}`}>
                                <p className="text-sm font-semibold">Materiales encontrados: {materialsProjectData.count}</p>
                                <p className="text-xs">Mejor match: {materialsProjectData.bestMatch}</p>
                              </div>
                            )}
                          </div>
                        </CardContent>
                      </Card>
                    </TabsContent>

                    {/* Simulation Tab */}
                    <TabsContent value="simulation">
                      <Card className={cardClasses}>
                        <CardHeader>
                          <CardTitle className={`flex items-center gap-2 ${darkMode ? 'text-white' : ''}`}>
                            <Satellite className={`w-5 h-5 ${darkMode ? 'text-[#ff6b00]' : 'text-blue-500'}`} />
                            Simulación de Condiciones Extremas
                          </CardTitle>
                          <CardDescription className={darkMode ? 'text-gray-400' : ''}>
                            Evalúa el material en ambientes espaciales
                          </CardDescription>
                        </CardHeader>
                        <CardContent>
                          {!simulationResult ? (
                            <div className="space-y-4">
                              <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-center">
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <Satellite className={`w-8 h-8 mx-auto mb-2 ${darkMode ? 'text-[#ff6b00]' : 'text-blue-500'}`} />
                                  <p className="text-sm font-medium">Vacío Espacial</p>
                                </div>
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <Radiation className="w-8 h-8 mx-auto mb-2 text-red-500" />
                                  <p className="text-sm font-medium">Radiación</p>
                                </div>
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <Activity className="w-8 h-8 mx-auto mb-2 text-green-500" />
                                  <p className="text-sm font-medium">Microgravedad</p>
                                </div>
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <Thermometer className="w-8 h-8 mx-auto mb-2 text-amber-500" />
                                  <p className="text-sm font-medium">Ciclo Térmico</p>
                                </div>
                              </div>
                              
                              <Button 
                                onClick={runSimulation} 
                                className="w-full" 
                                disabled={simulationLoading}
                                style={darkMode ? { backgroundColor: '#ff6b00' } : {}}
                              >
                                {simulationLoading ? (
                                  <>
                                    <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                    Simulando...
                                  </>
                                ) : (
                                  <>
                                    <Play className="w-4 h-4 mr-2" />
                                    Ejecutar Simulación
                                  </>
                                )}
                              </Button>
                            </div>
                          ) : (
                            <div className="space-y-4">
                              <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-blue-50'}`}>
                                  <p className="text-sm text-gray-400">Estabilidad en Vacío</p>
                                  <p className={`text-2xl font-bold ${simulationResult.vacuum_stability > 0.7 ? 'text-green-500' : 'text-yellow-500'}`}>
                                    {(simulationResult.vacuum_stability * 100).toFixed(1)}%
                                  </p>
                                </div>
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-red-50'}`}>
                                  <p className="text-sm text-gray-400">Degradación Radiación</p>
                                  <p className={`text-2xl font-bold ${simulationResult.radiation_degradation < 0.3 ? 'text-green-500' : 'text-red-500'}`}>
                                    {(simulationResult.radiation_degradation * 100).toFixed(1)}%
                                  </p>
                                </div>
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-green-50'}`}>
                                  <p className="text-sm text-gray-400">Efecto Microgravedad</p>
                                  <p className={`text-2xl font-bold ${simulationResult.microgravity_effect < 0.2 ? 'text-green-500' : 'text-yellow-500'}`}>
                                    {(simulationResult.microgravity_effect * 100).toFixed(1)}%
                                  </p>
                                </div>
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-amber-50'}`}>
                                  <p className="text-sm text-gray-400">Vida Ciclos Térmicos</p>
                                  <p className="text-2xl font-bold text-blue-500">
                                    {simulationResult.thermal_cycle_life}
                                  </p>
                                </div>
                              </div>
                              
                              <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                <h4 className="font-semibold mb-2">Vida Útil Estimada</h4>
                                <p className={`text-2xl font-bold ${darkMode ? 'text-[#ff6b00]' : 'text-blue-600'}`}>
                                  {simulationResult.estimated_lifetime}
                                </p>
                              </div>
                              
                              <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                <h4 className="font-semibold mb-2">Recomendaciones</h4>
                                <ul className="space-y-1 text-sm">
                                  {simulationResult.recommendations.map((rec, i) => (
                                    <li key={i} className="flex items-start gap-2">
                                      <CheckCircle2 className="w-4 h-4 text-green-500 mt-0.5" />
                                      {rec}
                                    </li>
                                  ))}
                                </ul>
                              </div>
                              
                              <Button onClick={() => setSimulationResult(null)} variant="outline" className={`w-full ${darkMode ? 'border-[#333]' : ''}`}>
                                <RefreshCw className="w-4 h-4 mr-2" />
                                Nueva Simulación
                              </Button>
                            </div>
                          )}
                        </CardContent>
                      </Card>
                    </TabsContent>

                    {/* Optimization Tab */}
                    <TabsContent value="optimization">
                      <Card className={cardClasses}>
                        <CardHeader>
                          <CardTitle className={`flex items-center gap-2 ${darkMode ? 'text-white' : ''}`}>
                            <Target className={`w-5 h-5 ${darkMode ? 'text-[#ff6b00]' : 'text-purple-500'}`} />
                            Optimización con Algoritmo Genético
                          </CardTitle>
                          <CardDescription className={darkMode ? 'text-gray-400' : ''}>
                            Encuentra la mejor composición y parámetros de síntesis
                          </CardDescription>
                        </CardHeader>
                        <CardContent>
                          {!optimizationResult ? (
                            <div className="space-y-4">
                              <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                <p className="text-sm">
                                  El algoritmo genético explorará múltiples combinaciones de:
                                </p>
                                <ul className="mt-2 space-y-1 text-sm">
                                  <li>• Proporciones de elementos</li>
                                  <li>• Temperatura de síntesis</li>
                                  <li>• Tiempo de reacción</li>
                                  <li>• pH óptimo</li>
                                </ul>
                              </div>
                              
                              <Button 
                                onClick={runOptimization} 
                                className="w-full" 
                                disabled={optimizationLoading}
                                style={darkMode ? { backgroundColor: '#ff6b00' } : {}}
                              >
                                {optimizationLoading ? (
                                  <>
                                    <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                    Optimizando (50 generaciones)...
                                  </>
                                ) : (
                                  <>
                                    <TrendingUp className="w-4 h-4 mr-2" />
                                    Ejecutar Optimización
                                  </>
                                )}
                              </Button>
                            </div>
                          ) : (
                            <div className="space-y-4">
                              <div className={`p-6 rounded-lg text-center ${darkMode ? 'bg-gradient-to-r from-[#1a1a1a] to-[#222]' : 'bg-gradient-to-r from-purple-50 to-indigo-50'}`}>
                                <Trophy className="w-12 h-12 mx-auto mb-2 text-yellow-500" />
                                <h3 className="text-xl font-bold">Mejor Composición</h3>
                                <p className={`text-2xl font-bold mt-2 ${darkMode ? 'text-[#ff6b00]' : 'text-purple-600'}`}>
                                  {optimizationResult.best_composition}
                                </p>
                                <p className="text-sm text-gray-400 mt-1">
                                  Score: {(optimizationResult.score * 100).toFixed(2)}%
                                </p>
                              </div>
                              
                              <div className="grid grid-cols-3 gap-4">
                                <div className={`p-4 rounded-lg text-center ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <p className="text-sm text-gray-400">Temperatura</p>
                                  <p className="text-xl font-bold">{optimizationResult.parameters.temperature}°C</p>
                                </div>
                                <div className={`p-4 rounded-lg text-center ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <p className="text-sm text-gray-400">Tiempo</p>
                                  <p className="text-xl font-bold">{optimizationResult.parameters.time}h</p>
                                </div>
                                <div className={`p-4 rounded-lg text-center ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <p className="text-sm text-gray-400">pH</p>
                                  <p className="text-xl font-bold">{optimizationResult.parameters.ph}</p>
                                </div>
                              </div>
                              
                              <p className="text-sm text-center text-gray-400">
                                Generaciones evaluadas: {optimizationResult.generations}
                              </p>
                              
                              <Button onClick={() => setOptimizationResult(null)} variant="outline" className={`w-full ${darkMode ? 'border-[#333]' : ''}`}>
                                <RefreshCw className="w-4 h-4 mr-2" />
                                Nueva Optimización
                              </Button>
                            </div>
                          )}
                        </CardContent>
                      </Card>
                    </TabsContent>

                    {/* DFT Tab */}
                    <TabsContent value="dft">
                      <Card className={cardClasses}>
                        <CardHeader>
                          <CardTitle className={`flex items-center gap-2 ${darkMode ? 'text-white' : ''}`}>
                            <Cpu className={`w-5 h-5 ${darkMode ? 'text-[#ff6b00]' : 'text-cyan-500'}`} />
                            Cálculo DFT (Quantum ESPRESSO)
                          </CardTitle>
                          <CardDescription className={darkMode ? 'text-gray-400' : ''}>
                            Simulación de estructura electrónica desde primeros principios
                          </CardDescription>
                        </CardHeader>
                        <CardContent>
                          {!dftResult ? (
                            <div className="space-y-4">
                              <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-center">
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <Zap className={`w-8 h-8 mx-auto mb-2 ${darkMode ? 'text-[#ff6b00]' : 'text-cyan-500'}`} />
                                  <p className="text-sm font-medium">SCF</p>
                                  <p className="text-xs text-gray-400">Campo autoconsistente</p>
                                </div>
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <LineChart className="w-8 h-8 mx-auto mb-2 text-blue-500" />
                                  <p className="text-sm font-medium">Bandas</p>
                                  <p className="text-xs text-gray-400">Estructura de bandas</p>
                                </div>
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <BarChart3 className="w-8 h-8 mx-auto mb-2 text-green-500" />
                                  <p className="text-sm font-medium">DOS</p>
                                  <p className="text-xs text-gray-400">Densidad de estados</p>
                                </div>
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <Atom className="w-8 h-8 mx-auto mb-2 text-purple-500" />
                                  <p className="text-sm font-medium">Relajación</p>
                                  <p className="text-xs text-gray-400">Optimización celda</p>
                                </div>
                              </div>
                              
                              {dftLoading && (
                                <div className="space-y-2">
                                  <div className="flex items-center justify-between">
                                    <span className="text-sm">Progreso del cálculo</span>
                                    <span className="text-sm">{Math.round(dftProgress)}%</span>
                                  </div>
                                  <Progress value={dftProgress} className="h-2" />
                                  <p className="text-xs text-center text-gray-400">
                                    Ejecutando cálculo SCF en la nube...
                                  </p>
                                </div>
                              )}
                              
                              <Button 
                                onClick={runDFT} 
                                className="w-full" 
                                disabled={dftLoading}
                                style={darkMode ? { backgroundColor: '#ff6b00' } : {}}
                              >
                                {dftLoading ? (
                                  <>
                                    <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                    Calculando DFT...
                                  </>
                                ) : (
                                  <>
                                    <Cloud className="w-4 h-4 mr-2" />
                                    Ejecutar Cálculo DFT
                                  </>
                                )}
                              </Button>
                              
                              <div className="flex gap-2">
                                <Button variant="outline" size="sm" className={`flex-1 ${darkMode ? 'border-[#333]' : ''}`}>
                                  <ExternalLink className="w-3 h-3 mr-1" />
                                  Google Colab
                                </Button>
                                <Button variant="outline" size="sm" className={`flex-1 ${darkMode ? 'border-[#333]' : ''}`}>
                                  SSH Cluster
                                </Button>
                              </div>
                            </div>
                          ) : (
                            <div className="space-y-4">
                              <div className={`p-4 rounded-lg ${dftResult.is_stable ? 'bg-green-900/20 border border-green-500' : 'bg-red-900/20 border border-red-500'}`}>
                                <div className="flex items-center gap-2">
                                  {dftResult.is_stable ? (
                                    <CheckCircle2 className="w-6 h-6 text-green-500" />
                                  ) : (
                                    <AlertTriangle className="w-6 h-6 text-red-500" />
                                  )}
                                  <span className="font-semibold">
                                    {dftResult.is_stable ? 'Material Estable' : 'Material Inestable'}
                                  </span>
                                </div>
                              </div>
                              
                              <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <p className="text-sm text-gray-400">Energía Total</p>
                                  <p className="text-xl font-bold">{dftResult.total_energy.toFixed(2)} Ry</p>
                                </div>
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <p className="text-sm text-gray-400">Band Gap</p>
                                  <p className="text-xl font-bold">{dftResult.band_gap.toFixed(3)} eV</p>
                                </div>
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <p className="text-sm text-gray-400">Celda Óptima</p>
                                  <p className="text-xl font-bold">{dftResult.optimized_lattice.toFixed(3)} Å</p>
                                </div>
                                <div className={`p-4 rounded-lg ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                                  <p className="text-sm text-gray-400">E. Fermi</p>
                                  <p className="text-xl font-bold">{dftResult.fermi_energy.toFixed(3)} eV</p>
                                </div>
                              </div>
                              
                              <p className="text-sm text-center text-gray-400">
                                Tiempo de cálculo: {dftResult.calculation_time}
                              </p>
                              
                              <Button onClick={() => setDftResult(null)} variant="outline" className={`w-full ${darkMode ? 'border-[#333]' : ''}`}>
                                <RefreshCw className="w-4 h-4 mr-2" />
                                Nuevo Cálculo
                              </Button>
                            </div>
                          )}
                        </CardContent>
                      </Card>
                    </TabsContent>

                    {/* Comparison Tab */}
                    <TabsContent value="comparison">
                      <Card className={cardClasses}>
                        <CardHeader>
                          <CardTitle className={`flex items-center gap-2 ${darkMode ? 'text-white' : ''}`}>
                            <Radar className={`w-5 h-5 ${darkMode ? 'text-[#ff6b00]' : 'text-blue-500'}`} />
                            Comparación con Materiales Comerciales
                          </CardTitle>
                        </CardHeader>
                        <CardContent>
                          {/* Radar Chart Placeholder */}
                          <div className={`p-6 rounded-lg mb-4 ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                            <div className="aspect-square max-w-md mx-auto relative">
                              {/* Simple radar visualization */}
                              <svg viewBox="0 0 200 200" className="w-full h-full">
                                {/* Background hexagon */}
                                <polygon points="100,20 180,60 180,140 100,180 20,140 20,60" fill="none" stroke={darkMode ? '#333' : '#ccc'} strokeWidth="1"/>
                                <polygon points="100,40 160,70 160,130 100,160 40,130 40,70" fill="none" stroke={darkMode ? '#333' : '#ccc'} strokeWidth="1"/>
                                <polygon points="100,60 140,80 140,120 100,140 60,120 60,80" fill="none" stroke={darkMode ? '#333' : '#ccc'} strokeWidth="1"/>
                                
                                {/* Data polygon - your material */}
                                <polygon 
                                  points={`${100 + 40 * physicalProps.quality_score},60 160,${70 + 30 * (1-physicalProps.porosity)} ${160 - 20 * physicalProps.density/8},130 100,${160 - 30 * physicalProps.thermal_conductivity/10} 40,${130 - 20 * physicalProps.elastic_modulus/300} ${60 + 30 * physicalProps.surface_area/100},80`}
                                  fill={darkMode ? 'rgba(255, 107, 0, 0.3)' : 'rgba(59, 130, 246, 0.3)'}
                                  stroke={darkMode ? '#ff6b00' : '#3b82f6'}
                                  strokeWidth="2"
                                />
                                
                                {/* Labels */}
                                <text x="100" y="10" textAnchor="middle" fill={darkMode ? '#fff' : '#333'} fontSize="10">Resistencia</text>
                                <text x="190" y="65" textAnchor="start" fill={darkMode ? '#fff' : '#333'} fontSize="10">Peso</text>
                                <text x="190" y="145" textAnchor="start" fill={darkMode ? '#fff' : '#333'} fontSize="10">Costo</text>
                                <text x="100" y="195" textAnchor="middle" fill={darkMode ? '#fff' : '#333'} fontSize="10">Estabilidad</text>
                                <text x="5" y="145" textAnchor="start" fill={darkMode ? '#fff' : '#333'} fontSize="10">Térmica</text>
                                <text x="5" y="65" textAnchor="start" fill={darkMode ? '#fff' : '#333'} fontSize="10">Dureza</text>
                              </svg>
                            </div>
                          </div>

                          {/* Comparison Table */}
                          <div className="overflow-x-auto">
                            <table className="w-full text-sm">
                              <thead>
                                <tr className={darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-100'}>
                                  <th className="p-3 text-left">Material</th>
                                  <th className="p-3 text-center">Resistencia</th>
                                  <th className="p-3 text-center">Peso</th>
                                  <th className="p-3 text-center">Costo</th>
                                  <th className="p-3 text-center">Cond. Térmica</th>
                                  <th className="p-3 text-center">Ranking</th>
                                </tr>
                              </thead>
                              <tbody>
                                {/* Your material */}
                                <tr className={`${darkMode ? 'bg-[#ff6b00]/20' : 'bg-blue-50'} border-l-4 ${darkMode ? 'border-[#ff6b00]' : 'border-blue-500'}`}>
                                  <td className="p-3 font-semibold">
                                    <Star className="w-4 h-4 inline mr-1 text-yellow-500" />
                                    Tu Material ({selectedMetal}-O)
                                  </td>
                                  <td className="p-3 text-center">{(physicalProps.elastic_modulus * 0.8).toFixed(0)}</td>
                                  <td className="p-3 text-center">{physicalProps.density.toFixed(2)}</td>
                                  <td className="p-3 text-center">~50</td>
                                  <td className="p-3 text-center">{physicalProps.thermal_conductivity.toFixed(1)}</td>
                                  <td className="p-3 text-center">
                                    <Badge className="bg-green-500">NUEVO</Badge>
                                  </td>
                                </tr>
                                {/* Commercial materials */}
                                {COMMERCIAL_MATERIALS.map((mat, i) => (
                                  <tr key={i} className={`border-t ${darkMode ? 'border-[#333]' : ''}`}>
                                    <td className="p-3">{mat.name}</td>
                                    <td className="p-3 text-center">{mat.strength}</td>
                                    <td className="p-3 text-center">{mat.weight}</td>
                                    <td className="p-3 text-center">{mat.cost}</td>
                                    <td className="p-3 text-center">{mat.thermal}</td>
                                    <td className="p-3 text-center">
                                      <Badge variant="outline">{mat.category}</Badge>
                                    </td>
                                  </tr>
                                ))}
                              </tbody>
                            </table>
                          </div>

                          {/* Rankings */}
                          <div className="mt-4 grid grid-cols-3 gap-4">
                            <div className={`p-4 rounded-lg text-center ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                              <Trophy className="w-6 h-6 mx-auto mb-1 text-yellow-500" />
                              <p className="text-sm font-medium">Aeroespacial</p>
                              <p className={`text-lg font-bold ${darkMode ? 'text-[#ff6b00]' : 'text-blue-600'}`}>#2</p>
                            </div>
                            <div className={`p-4 rounded-lg text-center ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                              <Zap className="w-6 h-6 mx-auto mb-1 text-amber-500" />
                              <p className="text-sm font-medium">Energía</p>
                              <p className={`text-lg font-bold ${darkMode ? 'text-[#ff6b00]' : 'text-blue-600'}`}>#1</p>
                            </div>
                            <div className={`p-4 rounded-lg text-center ${darkMode ? 'bg-[#1a1a1a]' : 'bg-slate-50'}`}>
                              <Shield className="w-6 h-6 mx-auto mb-1 text-green-500" />
                              <p className="text-sm font-medium">Médico</p>
                              <p className={`text-lg font-bold ${darkMode ? 'text-[#ff6b00]' : 'text-blue-600'}`}>#3</p>
                            </div>
                          </div>
                        </CardContent>
                      </Card>
                    </TabsContent>

                    {/* Files Tab */}
                    <TabsContent value="files">
                      <Card className={cardClasses}>
                        <CardHeader>
                          <CardTitle className={darkMode ? 'text-white' : ''}>Archivos para Supercomputadora</CardTitle>
                          <CardDescription className={darkMode ? 'text-gray-400' : ''}>
                            Quantum ESPRESSO, VASP, LAMMPS - Listos para ejecutar
                          </CardDescription>
                        </CardHeader>
                        <CardContent>
                          <Accordion type="single" collapsible className="space-y-2">
                            {/* Quantum ESPRESSO */}
                            <AccordionItem value="qe" className={`border rounded-lg ${darkMode ? 'border-[#333]' : ''}`}>
                              <AccordionTrigger className="px-4">
                                <div className="flex items-center gap-2">
                                  <Atom className="w-5 h-5 text-blue-500" />
                                  <span>Quantum ESPRESSO Input</span>
                                  <Badge variant="secondary" className="ml-2">SCF</Badge>
                                </div>
                              </AccordionTrigger>
                              <AccordionContent className="px-4">
                                <pre className={`p-4 rounded-lg overflow-x-auto text-sm ${darkMode ? 'bg-[#0a0a0a]' : 'bg-slate-900'} text-green-400`}>
{`&CONTROL
  calculation = 'scf'
  prefix = '${selectedMetal}O2'
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
 ${selectedMetal} ${ELEMENTS_DATABASE[selectedMetal]?.mw || 50}.00 ${selectedMetal}.pbe-spn-kjpaw_psl.1.0.0.UPF
 O  15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
 ${selectedMetal} 0.000 0.000 0.000
 O  0.305 0.305 0.000
 O -0.305 -0.305 0.000

K_POINTS automatic
 6 6 6 0 0 0`}
                                </pre>
                                <Button size="sm" className="mt-2" style={darkMode ? { backgroundColor: '#ff6b00' } : {}}>
                                  <Download className="w-3 h-3 mr-1" />
                                  Descargar .in
                                </Button>
                              </AccordionContent>
                            </AccordionItem>

                            {/* VASP */}
                            <AccordionItem value="vasp" className={`border rounded-lg ${darkMode ? 'border-[#333]' : ''}`}>
                              <AccordionTrigger className="px-4">
                                <div className="flex items-center gap-2">
                                  <Atom className="w-5 h-5 text-purple-500" />
                                  <span>VASP POSCAR</span>
                                  <Badge variant="secondary" className="ml-2">Struct</Badge>
                                </div>
                              </AccordionTrigger>
                              <AccordionContent className="px-4">
                                <pre className={`p-4 rounded-lg overflow-x-auto text-sm ${darkMode ? 'bg-[#0a0a0a]' : 'bg-slate-900'} text-green-400`}>
{`${selectedMetal}O2 - Generated by CosmicForge Lab
1.0
   4.593700 0.000000 0.000000
   0.000000 4.593700 0.000000
   0.000000 0.000000 2.958700
${selectedMetal} O
1 2
Direct
 0.000000 0.000000 0.000000 ${selectedMetal}
 0.305000 0.305000 0.000000 O
-0.305000 -0.305000 0.000000 O`}
                                </pre>
                                <Button size="sm" className="mt-2" style={darkMode ? { backgroundColor: '#ff6b00' } : {}}>
                                  <Download className="w-3 h-3 mr-1" />
                                  Descargar POSCAR
                                </Button>
                              </AccordionContent>
                            </AccordionItem>

                            {/* LAMMPS */}
                            <AccordionItem value="lammps" className={`border rounded-lg ${darkMode ? 'border-[#333]' : ''}`}>
                              <AccordionTrigger className="px-4">
                                <div className="flex items-center gap-2">
                                  <Activity className="w-5 h-5 text-green-500" />
                                  <span>LAMMPS Input</span>
                                  <Badge variant="secondary" className="ml-2">MD</Badge>
                                </div>
                              </AccordionTrigger>
                              <AccordionContent className="px-4">
                                <pre className={`p-4 rounded-lg overflow-x-auto text-sm ${darkMode ? 'bg-[#0a0a0a]' : 'bg-slate-900'} text-green-400`}>
{`# LAMMPS input for ${selectedMetal}O2 - CosmicForge Lab
# Molecular Dynamics Simulation

units metal
dimension 3
boundary p p p
atom_style full

lattice fcc 4.0
region box block 0 10 0 10 0 10
create_box 2 box

mass 1 ${ELEMENTS_DATABASE[selectedMetal]?.mw || 50}
mass 2 15.999

pair_style buck/coul/long 12.0
pair_coeff 1 1 0.001 0.1 0.0
pair_coeff 2 2 0.001 0.1 0.0
pair_coeff 1 2 0.001 0.2 0.0

kspace_style pppm 1.0e-4

timestep 0.001
thermo 100

# Minimization
minimize 1.0e-4 1.0e-6 100 1000

# NVT equilibration
velocity all create 300 12345
fix 1 all nvt temp 300 300 0.1
run 10000

# Production run
run 50000`}
                                </pre>
                                <Button size="sm" className="mt-2" style={darkMode ? { backgroundColor: '#ff6b00' } : {}}>
                                  <Download className="w-3 h-3 mr-1" />
                                  Descargar .in
                                </Button>
                              </AccordionContent>
                            </AccordionItem>

                            {/* Band structure */}
                            <AccordionItem value="bands" className={`border rounded-lg ${darkMode ? 'border-[#333]' : ''}`}>
                              <AccordionTrigger className="px-4">
                                <div className="flex items-center gap-2">
                                  <LineChart className="w-5 h-5 text-cyan-500" />
                                  <span>Band Structure Input</span>
                                  <Badge variant="secondary" className="ml-2">Bands</Badge>
                                </div>
                              </AccordionTrigger>
                              <AccordionContent className="px-4">
                                <pre className={`p-4 rounded-lg overflow-x-auto text-sm ${darkMode ? 'bg-[#0a0a0a]' : 'bg-slate-900'} text-green-400`}>
{`&CONTROL
  calculation = 'bands'
  prefix = '${selectedMetal}O2'
  outdir = './tmp'
/

&SYSTEM
  ibrav = 1
  celldm(1) = 8.0
  nat = 3
  ntyp = 2
  ecutwfc = 60.0
  ecutrho = 480.0
  nbnd = 20
/

&ELECTRONS
  conv_thr = 1.0d-8
/

ATOMIC_SPECIES
 ${selectedMetal} ${ELEMENTS_DATABASE[selectedMetal]?.mw || 50}.00 ${selectedMetal}.pbe-spn-kjpaw_psl.1.0.0.UPF
 O  15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
 ${selectedMetal} 0.000 0.000 0.000
 O  0.305 0.305 0.000
 O -0.305 -0.305 0.000

K_POINTS crystal
5
  0.0  0.0  0.0  20  ! Gamma
  0.5  0.0  0.0  20  ! X
  0.5  0.5  0.0  20  ! M
  0.0  0.0  0.5  20  ! Z
  0.0  0.0  0.0  1   ! Gamma`}
                                </pre>
                              </AccordionContent>
                            </AccordionItem>
                          </Accordion>
                        </CardContent>
                      </Card>
                    </TabsContent>

                    {/* PDF Tab */}
                    <TabsContent value="pdf">
                      <Card className={cardClasses}>
                        <CardHeader>
                          <CardTitle className={darkMode ? 'text-white' : ''}>Generar Reporte Completo</CardTitle>
                        </CardHeader>
                        <CardContent className="space-y-4">
                          <p className={darkMode ? 'text-gray-400' : 'text-slate-600'}>
                            Genera un reporte técnico completo que incluye:
                          </p>
                          <ul className="space-y-2 text-sm">
                            {[
                              'Información del objeto astrofísico',
                              'Propiedades calculadas del material',
                              'Resultados de simulación de condiciones extremas',
                              'Resultados de optimización',
                              'Resultados de cálculo DFT',
                              'Comparación con materiales comerciales',
                              'Archivos de simulación adjuntos',
                              'Reflexión filosófica (Yu. I. Manin)'
                            ].map((item, i) => (
                              <li key={i} className="flex items-center gap-2">
                                <CheckCircle2 className={`w-4 h-4 ${darkMode ? 'text-[#ff6b00]' : 'text-green-500'}`} />
                                <span className={darkMode ? 'text-gray-300' : ''}>{item}</span>
                              </li>
                            ))}
                          </ul>

                          <Button onClick={downloadPDF} size="lg" className="w-full" style={darkMode ? { backgroundColor: '#ff6b00' } : {}}>
                            <Download className="w-4 h-4 mr-2" />
                            Descargar Reporte Completo
                          </Button>
                        </CardContent>
                      </Card>
                    </TabsContent>
                  </Tabs>
                )}
              </>
            ) : (
              <Card className={cardClasses}>
                <CardContent className="py-16 text-center">
                  <div className="max-w-md mx-auto">
                    <Atom className={`w-16 h-16 mx-auto mb-4 ${darkMode ? 'text-[#ff6b00]' : 'text-blue-500'}`} />
                    <h2 className={`text-xl font-semibold mb-2 ${darkMode ? 'text-white' : ''}`}>
                      Bienvenido a CosmicForge Lab
                    </h2>
                    <p className={`${darkMode ? 'text-gray-400' : 'text-slate-500'} mb-6`}>
                      Selecciona un objeto astrofísico para comenzar a diseñar materiales.
                    </p>
                    <div className={`flex items-center justify-center gap-4 text-sm ${darkMode ? 'text-gray-500' : 'text-slate-400'}`}>
                      <span className="flex items-center gap-1">
                        <Database className="w-4 h-4" />
                        15 ejemplos
                      </span>
                      <span className="flex items-center gap-1">
                        <Beaker className="w-4 h-4" />
                        15 metales
                      </span>
                      <span className="flex items-center gap-1">
                        <FlaskConical className="w-4 h-4" />
                        8 aleaciones
                      </span>
                    </div>
                  </div>
                </CardContent>
              </Card>
            )}
          </div>
        </div>

        {/* Footer Quote */}
        <Card className={`mt-8 ${darkMode ? 'bg-gradient-to-r from-[#1a1a1a] to-[#222] border-[#333]' : 'bg-gradient-to-r from-rose-50 to-pink-50'}`}>
          <CardContent className="py-6">
            <div className="flex items-start gap-4">
              <BookOpen className={`w-8 h-8 flex-shrink-0 ${darkMode ? 'text-[#ff6b00]' : 'text-rose-500'}`} />
              <div>
                <p className={`italic ${darkMode ? 'text-gray-300' : 'text-slate-700'}`}>"{randomQuote}"</p>
                <p className={`text-sm mt-2 ${darkMode ? 'text-gray-500' : 'text-slate-500'}`}>
                  — Yuri I. Manin, "Lo demostrable e indemostrable"
                </p>
              </div>
            </div>
          </CardContent>
        </Card>
      </main>

      {/* Footer */}
      <footer className={`${darkMode ? 'bg-[#0a0a0a] border-[#333]' : 'bg-slate-800'} text-slate-300 py-6 mt-8 border-t`}>
        <div className="max-w-7xl mx-auto px-4 text-center">
          <p className="text-sm">
            CosmicForge Lab v3.5 NASA Edition — Diseño de Materiales Inspirado en Firmas Astrofísicas
          </p>
          <p className="text-xs text-slate-500 mt-1">
            Edición Personal — Simulación | DFT | Optimización | Comparación
          </p>
        </div>
      </footer>
    </div>
  )
}
