'use client'

import { useState, useCallback } from 'react'
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
import { 
  Sparkles, FlaskConical, Atom, FileText, Download, Database, 
  Thermometer, Gauge, Beaker, Shield, AlertTriangle, CheckCircle2,
  ExternalLink, Loader2, Settings, BookOpen
} from 'lucide-react'

// Astrophysical Examples Database
const ASTROPHYSICAL_EXAMPLES = {
  "Orion Nebula M42": {
    fractal_dimension: 1.654, criticality_score: 0.722, entropy: 0.019,
    anisotropy: 0.329, turbulence_beta: 2.278, lyapunov_max: -0.227,
    mode: "balanced", type: "Nebulosa de Emisión", distance: "1,344 ly"
  },
  "Crab Nebula M1": {
    fractal_dimension: 1.82, criticality_score: 0.89, entropy: 0.032,
    anisotropy: 0.456, turbulence_beta: 2.45, lyapunov_max: -0.15,
    mode: "turbulent", type: "Remanente Supernova", distance: "6,500 ly"
  },
  "Ring Nebula M57": {
    fractal_dimension: 1.45, criticality_score: 0.55, entropy: 0.012,
    anisotropy: 0.28, turbulence_beta: 1.95, lyapunov_max: -0.35,
    mode: "stable", type: "Nebulosa Planetaria", distance: "2,283 ly"
  },
  "Triangulum Galaxy M33": {
    fractal_dimension: 1.93, criticality_score: 0.78, entropy: 0.00000964,
    anisotropy: 0.23, turbulence_beta: 1.84, lyapunov_max: 0.11,
    mode: "balanced", type: "Galaxia Espiral", distance: "2.73M ly"
  },
  "Andromeda Galaxy M31": {
    fractal_dimension: 1.88, criticality_score: 0.85, entropy: 0.000012,
    anisotropy: 0.19, turbulence_beta: 1.92, lyapunov_max: 0.08,
    mode: "balanced", type: "Galaxia Espiral", distance: "2.537M ly"
  },
  "Eagle Nebula M16": {
    fractal_dimension: 1.78, criticality_score: 0.81, entropy: 0.028,
    anisotropy: 0.42, turbulence_beta: 2.35, lyapunov_max: -0.18,
    mode: "turbulent", type: "Nebulosa de Emisión", distance: "7,000 ly"
  },
  "Tarantula Nebula": {
    fractal_dimension: 1.91, criticality_score: 0.93, entropy: 0.045,
    anisotropy: 0.58, turbulence_beta: 2.68, lyapunov_max: 0.15,
    mode: "turbulent", type: "Región HII", distance: "160,000 ly"
  },
  "Pillars of Creation": {
    fractal_dimension: 1.85, criticality_score: 0.88, entropy: 0.035,
    anisotropy: 0.52, turbulence_beta: 2.55, lyapunov_max: -0.12,
    mode: "turbulent", type: "Nube Molecular", distance: "6,500 ly"
  },
  "Whirlpool Galaxy M51": {
    fractal_dimension: 1.95, criticality_score: 0.91, entropy: 0.000015,
    anisotropy: 0.22, turbulence_beta: 1.98, lyapunov_max: 0.12,
    mode: "turbulent", type: "Galaxia Interactuante", distance: "23M ly"
  },
  "Helix Nebula NGC 7293": {
    fractal_dimension: 1.48, criticality_score: 0.58, entropy: 0.011,
    anisotropy: 0.31, turbulence_beta: 2.02, lyapunov_max: -0.32,
    mode: "stable", type: "Nebulosa Planetaria", distance: "655 ly"
  }
}

// Elements Database
const ELEMENTS_DATABASE: Record<string, { name: string; atomic_num: number; mw: number; mp: number; compatible: string[] }> = {
  "Ti": { name: "Titanio", atomic_num: 22, mw: 47.867, mp: 1668, compatible: ["Al", "V", "Fe", "Ni", "O", "N", "C"] },
  "Al": { name: "Aluminio", atomic_num: 13, mw: 26.982, mp: 660, compatible: ["Ti", "Cu", "Mg", "Si", "Zn", "O"] },
  "Fe": { name: "Hierro", atomic_num: 26, mw: 55.845, mp: 1538, compatible: ["Ni", "Cr", "Co", "Mn", "C", "O"] },
  "Zn": { name: "Zinc", atomic_num: 30, mw: 65.38, mp: 420, compatible: ["Cu", "Al", "O", "S"] },
  "Cu": { name: "Cobre", atomic_num: 29, mw: 63.546, mp: 1085, compatible: ["Ni", "Zn", "Al", "Sn", "O"] },
  "Ni": { name: "Níquel", atomic_num: 28, mw: 58.693, mp: 1455, compatible: ["Fe", "Co", "Cu", "Cr", "Ti", "O"] },
  "Co": { name: "Cobalto", atomic_num: 27, mw: 58.933, mp: 1495, compatible: ["Ni", "Fe", "Mn", "O"] },
  "Mn": { name: "Manganeso", atomic_num: 25, mw: 54.938, mp: 1246, compatible: ["Fe", "Co", "O"] },
  "Ag": { name: "Plata", atomic_num: 47, mw: 107.868, mp: 962, compatible: ["Cu", "Au", "Pd", "O"] },
  "Au": { name: "Oro", atomic_num: 79, mw: 196.967, mp: 1064, compatible: ["Ag", "Pt", "Pd", "Cu"] },
  "Pt": { name: "Platino", atomic_num: 78, mw: 195.084, mp: 1768, compatible: ["Pd", "Au", "Rh", "Ir", "O"] },
  "Pd": { name: "Paladio", atomic_num: 46, mw: 106.42, mp: 1555, compatible: ["Pt", "Au", "Ag", "Ni", "O"] }
}

// Alloy Systems Database
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
    safety: ["Usar guantes refractarios", "Atmósfera inerte obligatoria", "Ventilación adecuada"]
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
    safety: ["Vacío alto (<10⁻³ mbar)", "Protección radiación UV", "Evitar contaminación con O₂, N₂, H₂"]
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
    safety: ["Evitar oxidación", "Controlar temperatura exacta"]
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
    safety: ["Fundente para evitar oxidación", "Ventilación"]
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
    safety: ["Evitar humedad (explosión)", "Degaseado obligatorio"]
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
      "6. Calcinación: 400°C (anatasa) o 800°C (rutilo)",
      "",
      "MÉTODO HIDROTERMAL:",
      "1. Mezclar TiCl₄ con NaOH 10M",
      "2. Transferir a autoclave",
      "3. Calentar a 180°C por 12 horas",
      "4. Lavar con agua hasta pH neutro",
      "5. Secar a 80°C"
    ],
    equipment: ["Matraz reacción", "Agitador magnético", "Horno mufla", "Autoclave (opcional)"],
    precursors: ["Isopropóxido de Ti", "TiCl₄", "Etanol", "HNO₃"],
    safety: ["Guantes y gafas", "Campana extractora", "TiCl₄ es corrosivo"]
  },
  "ZnO (Óxido de Zinc)": {
    name: "Óxido de Zinc (ZnO)",
    ratio: "Zn:O = 1:1",
    melting_point: "1975°C",
    density: "5.61 g/cm³",
    applications: ["Protectores solares", "Sensores de gas", "Varistores", "Catalizadores"],
    synthesis: [
      "MÉTODO DE PRECIPITACIÓN:",
      "1. Disolver Zn(NO₃)₂ en agua desionizada (0.5M)",
      "2. Preparar solución NaOH 1M",
      "3. Añadir NaOH gota a gota con agitación vigorosa",
      "4. Mantener pH 10-11",
      "5. Envejecer precipitado 2 horas",
      "6. Filtrar y lavar con agua/etanol",
      "7. Secar a 80°C por 12 horas",
      "8. Calcinación: 400-600°C por 4 horas"
    ],
    equipment: ["Matraz", "Agitador", "Horno mufla", "Centrífuga"],
    precursors: ["Zn(NO₃)₂", "NaOH", "Acetato de Zn", "Metanol"],
    safety: ["Guantes", "Evitar inhalación de polvo"]
  },
  "Fe₂O₃ (Hematita)": {
    name: "Óxido de Hierro (Hematita α-Fe₂O₃)",
    ratio: "Fe:O = 2:3",
    melting_point: "1565°C",
    density: "5.26 g/cm³",
    applications: ["Pigmentos", "Catálisis", "Sensores magnéticos", "Electrodos"],
    synthesis: [
      "MÉTODO DE PRECIPITACIÓN:",
      "1. Disolver FeCl₃·6H₂O en agua (0.2M)",
      "2. Añadir NH₄OH hasta pH 8-9",
      "3. Envejecer precipitado 24 horas",
      "4. Filtrar y lavar con agua caliente",
      "5. Secar a 100°C",
      "6. Calcinación: 500-700°C por 3 horas"
    ],
    equipment: ["Matraz", "Agitador", "Horno mufla", "Autoclave"],
    precursors: ["FeCl₃·6H₂O", "Fe(NO₃)₃", "NH₄OH", "NaOH"],
    safety: ["Manchas la piel", "Usar guantes"]
  }
}

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

export default function CosmicForgeLab() {
  const [selectedExample, setSelectedExample] = useState<string>("Orion Nebula M42")
  const [selectedMetal, setSelectedMetal] = useState<string>("Ti")
  const [materialType, setMaterialType] = useState<string>("Oxido")
  const [selectedAlloy, setSelectedAlloy] = useState<string>("")
  const [useAlloy, setUseAlloy] = useState(false)
  
  const [astroData, setAstroData] = useState<AstrophysicalData | null>(null)
  const [physicalProps, setPhysicalProps] = useState<PhysicalProperties | null>(null)
  const [recipe, setRecipe] = useState<Recipe | null>(null)
  const [isLoading, setIsLoading] = useState(false)
  
  // Manual editor state
  const [manualName, setManualName] = useState("Objeto Personalizado")
  const [manualFd, setManualFd] = useState([1.5])
  const [manualCs, setManualCs] = useState([0.5])
  const [manualEntropy, setManualEntropy] = useState("0.01")
  const [manualAnisotropy, setManualAnisotropy] = useState([0.3])
  const [manualTurbulence, setManualTurbulence] = useState([2.0])
  const [manualLyapunov, setManualLyapunov] = useState("-0.1")
  
  const [inputMethod, setInputMethod] = useState<"examples" | "manual">("examples")

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
          alloyKey: selectedAlloy
        })
      })
      
      const blob = await response.blob()
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = `CosmicForge_Report_${astroData.object_name.replace(/\s+/g, '_')}.pdf`
      document.body.appendChild(a)
      a.click()
      document.body.removeChild(a)
      URL.revokeObjectURL(url)
    } catch (error) {
      console.error('Error downloading PDF:', error)
    }
  }, [astroData, physicalProps, recipe, selectedAlloy])

  const clearData = useCallback(() => {
    setAstroData(null)
    setPhysicalProps(null)
    setRecipe(null)
  }, [])

  const getAlloyOptions = () => {
    return Object.keys(ALLOY_SYSTEMS).filter(key => key.includes(selectedMetal))
  }

  const currentAlloy = selectedAlloy ? ALLOY_SYSTEMS[selectedAlloy] : null
  const randomQuote = MANIN_QUOTES[Math.floor(Math.random() * MANIN_QUOTES.length)]

  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-50 via-blue-50 to-indigo-50">
      {/* Header */}
      <header className="bg-white/80 backdrop-blur-sm border-b border-slate-200 sticky top-0 z-50">
        <div className="max-w-7xl mx-auto px-4 py-4">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-3">
              <div className="p-2 bg-gradient-to-br from-blue-500 to-indigo-600 rounded-xl">
                <Atom className="w-8 h-8 text-white" />
              </div>
              <div>
                <h1 className="text-2xl font-bold text-slate-800">CosmicForge Lab</h1>
                <p className="text-sm text-slate-500">Diseño de Materiales Inspirado en Firmas Astrofísicas</p>
              </div>
            </div>
            <Badge variant="outline" className="text-sm px-3 py-1">
              v3.1 Personal Edition
            </Badge>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-4 py-6">
        <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
          {/* Sidebar */}
          <aside className="lg:col-span-1 space-y-4">
            {/* Input Method Selection */}
            <Card className="shadow-sm border-slate-200">
              <CardHeader className="pb-3">
                <CardTitle className="text-lg flex items-center gap-2">
                  <Database className="w-5 h-5 text-blue-500" />
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
                    <Label>Objeto Astrofísico</Label>
                    <Select value={selectedExample} onValueChange={setSelectedExample}>
                      <SelectTrigger>
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        {Object.keys(ASTROPHYSICAL_EXAMPLES).map(name => (
                          <SelectItem key={name} value={name}>{name}</SelectItem>
                        ))}
                      </SelectContent>
                    </Select>
                    
                    {ASTROPHYSICAL_EXAMPLES[selectedExample as keyof typeof ASTROPHYSICAL_EXAMPLES] && (
                      <div className="p-3 bg-blue-50 rounded-lg text-sm">
                        <p><strong>Tipo:</strong> {ASTROPHYSICAL_EXAMPLES[selectedExample as keyof typeof ASTROPHYSICAL_EXAMPLES].type}</p>
                        <p><strong>Distancia:</strong> {ASTROPHYSICAL_EXAMPLES[selectedExample as keyof typeof ASTROPHYSICAL_EXAMPLES].distance}</p>
                      </div>
                    )}
                    
                    <Button onClick={loadExample} className="w-full" variant="default">
                      Cargar Ejemplo
                    </Button>
                  </TabsContent>
                  
                  <TabsContent value="manual" className="space-y-3 mt-3">
                    <div>
                      <Label>Nombre del Objeto</Label>
                      <Input value={manualName} onChange={(e) => setManualName(e.target.value)} />
                    </div>
                    
                    <div>
                      <Label>Dimensión Fractal: {manualFd[0].toFixed(3)}</Label>
                      <Slider value={manualFd} onValueChange={setManualFd} min={0} max={3} step={0.001} />
                    </div>
                    
                    <div>
                      <Label>Criticalidad: {manualCs[0].toFixed(3)}</Label>
                      <Slider value={manualCs} onValueChange={setManualCs} min={0} max={1} step={0.001} />
                    </div>
                    
                    <div>
                      <Label>Entropía</Label>
                      <Input type="number" value={manualEntropy} onChange={(e) => setManualEntropy(e.target.value)} step="0.000001" />
                    </div>
                    
                    <div>
                      <Label>Anisotropía: {manualAnisotropy[0].toFixed(3)}</Label>
                      <Slider value={manualAnisotropy} onValueChange={setManualAnisotropy} min={0} max={1} step={0.001} />
                    </div>
                    
                    <div>
                      <Label>Turbulencia β: {manualTurbulence[0].toFixed(3)}</Label>
                      <Slider value={manualTurbulence} onValueChange={setManualTurbulence} min={0} max={5} step={0.001} />
                    </div>
                    
                    <div>
                      <Label>Lyapunov máx</Label>
                      <Input type="number" value={manualLyapunov} onChange={(e) => setManualLyapunov(e.target.value)} step="0.001" />
                    </div>
                    
                    <Button onClick={applyManualData} className="w-full" variant="default">
                      Aplicar Parámetros
                    </Button>
                  </TabsContent>
                </Tabs>
              </CardContent>
            </Card>

            {/* Material Configuration */}
            <Card className="shadow-sm border-slate-200">
              <CardHeader className="pb-3">
                <CardTitle className="text-lg flex items-center gap-2">
                  <Settings className="w-5 h-5 text-green-500" />
                  Configuración del Material
                </CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <div>
                  <Label>Tipo de Material</Label>
                  <Select value={materialType} onValueChange={setMaterialType}>
                    <SelectTrigger>
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
                  <Label>Metal Base</Label>
                  <Select value={selectedMetal} onValueChange={setSelectedMetal}>
                    <SelectTrigger>
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
                  <div className="p-3 bg-green-50 rounded-lg text-sm">
                    <p><strong>Compatible con:</strong></p>
                    <div className="flex flex-wrap gap-1 mt-1">
                      {ELEMENTS_DATABASE[selectedMetal].compatible.map(el => (
                        <Badge key={el} variant="secondary" className="text-xs">{el}</Badge>
                      ))}
                    </div>
                  </div>
                )}
                
                <Separator />
                
                <div className="flex items-center space-x-2">
                  <input
                    type="checkbox"
                    id="useAlloy"
                    checked={useAlloy}
                    onChange={(e) => setUseAlloy(e.target.checked)}
                    className="h-4 w-4 rounded border-gray-300"
                  />
                  <Label htmlFor="useAlloy">Usar aleación predefinida</Label>
                </div>
                
                {useAlloy && (
                  <div>
                    <Label>Sistema de Aleación</Label>
                    <Select value={selectedAlloy} onValueChange={setSelectedAlloy}>
                      <SelectTrigger>
                        <SelectValue placeholder="Seleccionar aleación..." />
                      </SelectTrigger>
                      <SelectContent>
                        {getAlloyOptions().map(key => (
                          <SelectItem key={key} value={key}>
                            {ALLOY_SYSTEMS[key].name}
                          </SelectItem>
                        ))}
                      </SelectContent>
                    </Select>
                    
                    {currentAlloy && (
                      <div className="p-3 bg-purple-50 rounded-lg text-sm mt-2">
                        <p><strong>Ratio:</strong> {currentAlloy.ratio}</p>
                        <p><strong>P.F.:</strong> {currentAlloy.melting_point}</p>
                        <p><strong>Densidad:</strong> {currentAlloy.density}</p>
                      </div>
                    )}
                  </div>
                )}
              </CardContent>
            </Card>

            {/* Clear Button */}
            {astroData && (
              <Button onClick={clearData} variant="outline" className="w-full">
                Limpiar Datos
              </Button>
            )}
          </aside>

          {/* Main Content */}
          <div className="lg:col-span-3 space-y-6">
            {/* Astro Data Display */}
            {astroData ? (
              <>
                <Card className="shadow-sm border-slate-200">
                  <CardHeader>
                    <CardTitle className="flex items-center gap-2">
                      <Sparkles className="w-5 h-5 text-amber-500" />
                      Firma Astrofísica Detectada
                    </CardTitle>
                    <CardDescription>
                      Objeto: {astroData.object_name} | Fuente: {astroData.source}
                    </CardDescription>
                  </CardHeader>
                  <CardContent>
                    <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                      <div className="p-4 bg-slate-50 rounded-lg">
                        <p className="text-sm text-slate-500">Dimensión Fractal</p>
                        <p className="text-2xl font-bold text-slate-800">{astroData.fractal_dimension.toFixed(4)}</p>
                      </div>
                      <div className="p-4 bg-slate-50 rounded-lg">
                        <p className="text-sm text-slate-500">Criticalidad</p>
                        <p className="text-2xl font-bold text-slate-800">{astroData.criticality_score.toFixed(4)}</p>
                      </div>
                      <div className="p-4 bg-slate-50 rounded-lg">
                        <p className="text-sm text-slate-500">Entropía</p>
                        <p className="text-2xl font-bold text-slate-800">{astroData.entropy.toFixed(6)}</p>
                      </div>
                      <div className="p-4 bg-slate-50 rounded-lg">
                        <p className="text-sm text-slate-500">Anisotropía</p>
                        <p className="text-2xl font-bold text-slate-800">{astroData.anisotropy.toFixed(4)}</p>
                      </div>
                      <div className="p-4 bg-slate-50 rounded-lg">
                        <p className="text-sm text-slate-500">Turbulencia β</p>
                        <p className="text-2xl font-bold text-slate-800">{astroData.turbulence_beta.toFixed(4)}</p>
                      </div>
                      <div className="p-4 bg-slate-50 rounded-lg">
                        <p className="text-sm text-slate-500">Lyapunov máx</p>
                        <p className="text-2xl font-bold text-slate-800">{astroData.lyapunov_max.toFixed(6)}</p>
                      </div>
                      <div className="p-4 bg-slate-50 rounded-lg">
                        <p className="text-sm text-slate-500">Modo</p>
                        <Badge variant={astroData.mode === 'turbulent' ? 'destructive' : astroData.mode === 'stable' ? 'default' : 'secondary'}>
                          {astroData.mode}
                        </Badge>
                      </div>
                      <div className="p-4 bg-slate-50 rounded-lg">
                        <p className="text-sm text-slate-500">Tipo</p>
                        <p className="text-sm font-medium text-slate-700">{astroData.type || 'N/A'}</p>
                      </div>
                    </div>
                    
                    <Button 
                      onClick={generateRecipe} 
                      className="w-full mt-4" 
                      size="lg"
                      disabled={isLoading}
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

                {/* Results */}
                {physicalProps && recipe && (
                  <Tabs defaultValue="properties" className="space-y-4">
                    <TabsList className="grid w-full grid-cols-5">
                      <TabsTrigger value="properties">Propiedades</TabsTrigger>
                      <TabsTrigger value="synthesis">Síntesis</TabsTrigger>
                      <TabsTrigger value="production">Producción</TabsTrigger>
                      <TabsTrigger value="files">Archivos</TabsTrigger>
                      <TabsTrigger value="pdf">PDF</TabsTrigger>
                    </TabsList>

                    {/* Properties Tab */}
                    <TabsContent value="properties">
                      <Card className="shadow-sm border-slate-200">
                        <CardHeader>
                          <CardTitle>Propiedades del Material</CardTitle>
                        </CardHeader>
                        <CardContent>
                          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                            <div className="p-4 bg-blue-50 rounded-lg">
                              <p className="text-sm text-blue-600">Porosidad</p>
                              <p className="text-2xl font-bold text-blue-700">{(physicalProps.porosity * 100).toFixed(1)}%</p>
                            </div>
                            <div className="p-4 bg-green-50 rounded-lg">
                              <p className="text-sm text-green-600">Densidad</p>
                              <p className="text-2xl font-bold text-green-700">{physicalProps.density.toFixed(3)} g/cm³</p>
                            </div>
                            <div className="p-4 bg-amber-50 rounded-lg">
                              <p className="text-sm text-amber-600">Conductividad Térmica</p>
                              <p className="text-2xl font-bold text-amber-700">{physicalProps.thermal_conductivity.toFixed(2)} W/mK</p>
                            </div>
                            <div className="p-4 bg-purple-50 rounded-lg">
                              <p className="text-sm text-purple-600">Módulo Elástico</p>
                              <p className="text-2xl font-bold text-purple-700">{physicalProps.elastic_modulus.toFixed(2)} GPa</p>
                            </div>
                            <div className="p-4 bg-rose-50 rounded-lg">
                              <p className="text-sm text-rose-600">Área Superficial</p>
                              <p className="text-2xl font-bold text-rose-700">{physicalProps.surface_area.toFixed(1)} m²/g</p>
                            </div>
                            <div className="p-4 bg-cyan-50 rounded-lg">
                              <p className="text-sm text-cyan-600">Band Gap</p>
                              <p className="text-2xl font-bold text-cyan-700">{physicalProps.band_gap.toFixed(3)} eV</p>
                            </div>
                            <div className="p-4 bg-indigo-50 rounded-lg">
                              <p className="text-sm text-indigo-600">Calidad</p>
                              <p className="text-2xl font-bold text-indigo-700">{(physicalProps.quality_score * 100).toFixed(1)}%</p>
                            </div>
                            <div className="p-4 bg-teal-50 rounded-lg">
                              <p className="text-sm text-teal-600">Energía Activación</p>
                              <p className="text-2xl font-bold text-teal-700">{physicalProps.activation_energy.toFixed(4)} eV</p>
                            </div>
                          </div>

                          <Separator className="my-4" />

                          <div className="p-4 bg-slate-50 rounded-lg">
                            <h4 className="font-semibold mb-2 flex items-center gap-2">
                              <ExternalLink className="w-4 h-4" />
                              Validación con Bases de Datos
                            </h4>
                            <div className="flex flex-wrap gap-2">
                              <a 
                                href={`https://materialsproject.org/materials?search=${selectedMetal}O2`} 
                                target="_blank" 
                                rel="noopener noreferrer"
                                className="text-blue-600 hover:underline text-sm"
                              >
                                Materials Project - {selectedMetal}O₂
                              </a>
                              <span className="text-slate-400">|</span>
                              <a 
                                href={`https://www.aflowlib.org/?search=${selectedMetal}O2`} 
                                target="_blank" 
                                rel="noopener noreferrer"
                                className="text-blue-600 hover:underline text-sm"
                              >
                                AFLOW Library
                              </a>
                              <span className="text-slate-400">|</span>
                              <a 
                                href={`http://oqmd.org/materials?search=${selectedMetal}O2`} 
                                target="_blank" 
                                rel="noopener noreferrer"
                                className="text-blue-600 hover:underline text-sm"
                              >
                                OQMD Database
                              </a>
                            </div>
                          </div>
                        </CardContent>
                      </Card>
                    </TabsContent>

                    {/* Synthesis Tab */}
                    <TabsContent value="synthesis">
                      <Card className="shadow-sm border-slate-200">
                        <CardHeader>
                          <CardTitle>Receta de Síntesis</CardTitle>
                        </CardHeader>
                        <CardContent>
                          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                            <div>
                              <h4 className="font-semibold mb-3 flex items-center gap-2">
                                <Thermometer className="w-4 h-4 text-orange-500" />
                                Condiciones de Reacción
                              </h4>
                              <div className="space-y-2">
                                <p><strong>Temperatura:</strong> {recipe.reaction_conditions.temperature_C.toFixed(0)}°C</p>
                                <p><strong>Tiempo:</strong> {recipe.reaction_conditions.reaction_time_hours.toFixed(2)} horas</p>
                                <p><strong>pH:</strong> {recipe.reaction_conditions.pH.toFixed(1)}</p>
                              </div>
                            </div>
                            <div>
                              <h4 className="font-semibold mb-3 flex items-center gap-2">
                                <Beaker className="w-4 h-4 text-blue-500" />
                                Precursores
                              </h4>
                              <ul className="space-y-1">
                                {recipe.precursors.map((p, i) => (
                                  <li key={i} className="flex items-center gap-2">
                                    <CheckCircle2 className="w-4 h-4 text-green-500" />
                                    {p}
                                  </li>
                                ))}
                              </ul>
                            </div>
                          </div>

                          <Separator className="my-4" />

                          <h4 className="font-semibold mb-3">Procedimiento General</h4>
                          <ol className="space-y-2">
                            {recipe.step_by_step.map((step, i) => (
                              <li key={i} className="flex gap-3 p-3 bg-slate-50 rounded-lg">
                                <span className="flex-shrink-0 w-6 h-6 bg-blue-500 text-white rounded-full flex items-center justify-center text-sm">
                                  {i + 1}
                                </span>
                                <span>{step}</span>
                              </li>
                            ))}
                          </ol>
                        </CardContent>
                      </Card>
                    </TabsContent>

                    {/* Production Tab */}
                    <TabsContent value="production">
                      <Card className="shadow-sm border-slate-200">
                        <CardHeader>
                          <CardTitle>Guía de Producción Detallada</CardTitle>
                        </CardHeader>
                        <CardContent>
                          {currentAlloy ? (
                            <>
                              <div className="p-4 bg-gradient-to-r from-blue-50 to-indigo-50 rounded-lg mb-4">
                                <h3 className="text-xl font-bold text-slate-800">{currentAlloy.name}</h3>
                                <div className="grid grid-cols-3 gap-4 mt-2 text-sm">
                                  <div><strong>Ratio:</strong> {currentAlloy.ratio}</div>
                                  <div><strong>P.F.:</strong> {currentAlloy.melting_point}</div>
                                  <div><strong>Densidad:</strong> {currentAlloy.density}</div>
                                </div>
                              </div>

                              <div className="mb-4">
                                <h4 className="font-semibold mb-2">Aplicaciones</h4>
                                <div className="flex flex-wrap gap-2">
                                  {currentAlloy.applications.map((app, i) => (
                                    <Badge key={i} variant="secondary">{app}</Badge>
                                  ))}
                                </div>
                              </div>

                              <h4 className="font-semibold mb-2">Proceso de Síntesis Paso a Paso</h4>
                              <div className="space-y-2 mb-4">
                                {currentAlloy.synthesis.map((step, i) => (
                                  step.trim() ? (
                                    <div key={i} className="p-3 bg-slate-50 rounded-lg border-l-4 border-blue-400">
                                      {step}
                                    </div>
                                  ) : <div key={i} className="h-4" />
                                ))}
                              </div>

                              <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                                <div className="p-4 bg-slate-50 rounded-lg">
                                  <h5 className="font-semibold mb-2 flex items-center gap-2">
                                    <Gauge className="w-4 h-4 text-blue-500" />
                                    Equipamiento
                                  </h5>
                                  <ul className="space-y-1 text-sm">
                                    {currentAlloy.equipment.map((eq, i) => (
                                      <li key={i}>• {eq}</li>
                                    ))}
                                  </ul>
                                </div>
                                <div className="p-4 bg-slate-50 rounded-lg">
                                  <h5 className="font-semibold mb-2 flex items-center gap-2">
                                    <Beaker className="w-4 h-4 text-green-500" />
                                    Precursores
                                  </h5>
                                  <ul className="space-y-1 text-sm">
                                    {currentAlloy.precursors.map((p, i) => (
                                      <li key={i}>• {p}</li>
                                    ))}
                                  </ul>
                                </div>
                                <div className="p-4 bg-slate-50 rounded-lg">
                                  <h5 className="font-semibold mb-2 flex items-center gap-2">
                                    <Shield className="w-4 h-4 text-red-500" />
                                    Seguridad
                                  </h5>
                                  <ul className="space-y-1 text-sm">
                                    {currentAlloy.safety.map((s, i) => (
                                      <li key={i}>• {s}</li>
                                    ))}
                                  </ul>
                                </div>
                              </div>
                            </>
                          ) : (
                            <div className="text-center py-8 text-slate-500">
                              <p>Selecciona una aleación predefinida para ver la guía de producción detallada.</p>
                              <p className="text-sm mt-2">O usa la pestaña "Síntesis" para ver el procedimiento general.</p>
                            </div>
                          )}
                        </CardContent>
                      </Card>
                    </TabsContent>

                    {/* Files Tab */}
                    <TabsContent value="files">
                      <Card className="shadow-sm border-slate-200">
                        <CardHeader>
                          <CardTitle>Archivos de Simulación</CardTitle>
                          <CardDescription>Archivos para Quantum ESPRESSO, VASP y LAMMPS</CardDescription>
                        </CardHeader>
                        <CardContent>
                          <Accordion type="single" collapsible className="space-y-2">
                            <AccordionItem value="qe" className="border rounded-lg">
                              <AccordionTrigger className="px-4">Quantum ESPRESSO Input</AccordionTrigger>
                              <AccordionContent className="px-4">
                                <pre className="bg-slate-900 text-green-400 p-4 rounded-lg overflow-x-auto text-sm">
{`&CONTROL
  calculation = 'scf'
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
/

&ELECTRONS
  conv_thr = 1.0d-8
/

ATOMIC_SPECIES
 ${selectedMetal} ${ELEMENTS_DATABASE[selectedMetal]?.mw || 50}.00 Ti.pbe-spn-kjpaw_psl.1.0.0.UPF
 O  15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
 ${selectedMetal} 0.000 0.000 0.000
 O  0.305 0.305 0.000
 O -0.305 -0.305 0.000

K_POINTS automatic
 4 4 4 0 0 0`}
                                </pre>
                              </AccordionContent>
                            </AccordionItem>

                            <AccordionItem value="vasp" className="border rounded-lg">
                              <AccordionTrigger className="px-4">VASP POSCAR</AccordionTrigger>
                              <AccordionContent className="px-4">
                                <pre className="bg-slate-900 text-green-400 p-4 rounded-lg overflow-x-auto text-sm">
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
                              </AccordionContent>
                            </AccordionItem>

                            <AccordionItem value="lammps" className="border rounded-lg">
                              <AccordionTrigger className="px-4">LAMMPS Input</AccordionTrigger>
                              <AccordionContent className="px-4">
                                <pre className="bg-slate-900 text-green-400 p-4 rounded-lg overflow-x-auto text-sm">
{`# LAMMPS input for ${selectedMetal}O2 - CosmicForge Lab
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

timestep 0.001
thermo 100

minimize 1.0e-4 1.0e-6 100 1000
run 1000`}
                                </pre>
                              </AccordionContent>
                            </AccordionItem>
                          </Accordion>
                        </CardContent>
                      </Card>
                    </TabsContent>

                    {/* PDF Tab */}
                    <TabsContent value="pdf">
                      <Card className="shadow-sm border-slate-200">
                        <CardHeader>
                          <CardTitle>Generar Reporte PDF</CardTitle>
                        </CardHeader>
                        <CardContent className="space-y-4">
                          <p className="text-slate-600">
                            Genera un reporte técnico completo en formato PDF que incluye:
                          </p>
                          <ul className="space-y-2 text-sm">
                            <li className="flex items-center gap-2">
                              <CheckCircle2 className="w-4 h-4 text-green-500" />
                              Información del objeto astrofísico
                            </li>
                            <li className="flex items-center gap-2">
                              <CheckCircle2 className="w-4 h-4 text-green-500" />
                              Propiedades calculadas del material
                            </li>
                            <li className="flex items-center gap-2">
                              <CheckCircle2 className="w-4 h-4 text-green-500" />
                              Guía de producción detallada
                            </li>
                            <li className="flex items-center gap-2">
                              <CheckCircle2 className="w-4 h-4 text-green-500" />
                              Lista de equipamiento y precursores
                            </li>
                            <li className="flex items-center gap-2">
                              <CheckCircle2 className="w-4 h-4 text-green-500" />
                              Reflexión filosófica (Yu. I. Manin)
                            </li>
                          </ul>

                          <Button onClick={downloadPDF} size="lg" className="w-full">
                            <Download className="w-4 h-4 mr-2" />
                            Descargar PDF
                          </Button>
                        </CardContent>
                      </Card>
                    </TabsContent>
                  </Tabs>
                )}
              </>
            ) : (
              <Card className="shadow-sm border-slate-200">
                <CardContent className="py-16 text-center">
                  <div className="max-w-md mx-auto">
                    <Atom className="w-16 h-16 mx-auto text-blue-500 mb-4" />
                    <h2 className="text-xl font-semibold mb-2">Bienvenido a CosmicForge Lab</h2>
                    <p className="text-slate-500 mb-6">
                      Selecciona un objeto astrofísico de la base de datos o ingresa parámetros 
                      manuales para comenzar a diseñar materiales inspirados en firmas cósmicas.
                    </p>
                    <div className="flex items-center justify-center gap-4 text-sm text-slate-400">
                      <span className="flex items-center gap-1">
                        <Database className="w-4 h-4" />
                        10 ejemplos
                      </span>
                      <span className="flex items-center gap-1">
                        <Beaker className="w-4 h-4" />
                        12 metales
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
        <Card className="mt-8 shadow-sm border-slate-200 bg-gradient-to-r from-rose-50 to-pink-50">
          <CardContent className="py-6">
            <div className="flex items-start gap-4">
              <BookOpen className="w-8 h-8 text-rose-500 flex-shrink-0" />
              <div>
                <p className="italic text-slate-700">"{randomQuote}"</p>
                <p className="text-sm text-slate-500 mt-2">— Yuri I. Manin, "Lo demostrable e indemostrable"</p>
              </div>
            </div>
          </CardContent>
        </Card>
      </main>

      {/* Footer */}
      <footer className="bg-slate-800 text-slate-300 py-6 mt-8">
        <div className="max-w-7xl mx-auto px-4 text-center">
          <p className="text-sm">
            CosmicForge Lab v3.1 — Diseño de Materiales Inspirado en Firmas Astrofísicas
          </p>
          <p className="text-xs text-slate-500 mt-1">
            Edición Personal — Uso Exclusivo
          </p>
        </div>
      </footer>
    </div>
  )
}
