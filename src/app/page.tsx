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
import { Textarea } from '@/components/ui/textarea'
import { 
  Sparkles, FlaskConical, Atom, FileText, Download, Database, 
  Thermometer, Gauge, Beaker, Shield, AlertTriangle, CheckCircle2,
  ExternalLink, Loader2, Settings, BookOpen, Play, Code, Copy, 
  CheckCircle, XCircle, AlertCircle
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
  "Si": { name: "Silicio", atomic_num: 14, mw: 28.086, mp: 1414, compatible: ["O", "C", "Ge"] }
}

// Crystal Structures for nanoHUB
const CRYSTAL_STRUCTURES_LIST = [
  { id: "Si_diamond", name: "Silicio Diamond", formula: "Si", system: "cubic", nat: 2 },
  { id: "TiO2_rutile", name: "TiO2 Rutilo", formula: "TiO2", system: "tetragonal", nat: 6 },
  { id: "TiO2_anatase", name: "TiO2 Anatasa", formula: "TiO2", system: "tetragonal", nat: 12 },
  { id: "Ti2O3_corundum", name: "Ti2O3 Corundum", formula: "Ti2O3", system: "trigonal", nat: 10 },
  { id: "TiO_rocksalt", name: "TiO Rocksalt", formula: "TiO", system: "cubic", nat: 2 },
  { id: "ZnO_wurtzite", name: "ZnO Wurtzita", formula: "ZnO", system: "hexagonal", nat: 4 },
  { id: "Fe2O3_hematite", name: "Fe2O3 Hematita", formula: "Fe2O3", system: "trigonal", nat: 10 },
  { id: "Al2O3_corundum", name: "Al2O3 Corundum", formula: "Al2O3", system: "trigonal", nat: 10 }
]

// Manin philosophical quotes
const MANIN_QUOTES = [
  "La matemática es la parte más pura de la cultura humana, y sin embargo está profundamente conectada con la realidad física.",
  "La demostrabilidad es un concepto que trasciende la lógica formal: es la búsqueda de certeza en un universo de incertidumbre.",
  "Entre lo demostrable y lo indemostrable existe una zona fronteriza donde nace la creatividad científica."
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

interface QEFiles {
  pwInput: string;
  runScript: string;
  readme: string;
}

interface MaterialStatus {
  status: 'EXISTENTE' | 'FUTURO' | 'DESCONOCIDO';
  color: string;
  icon: string;
  description: string;
  sources: string[];
}

export default function CosmicForgeLab() {
  const [selectedExample, setSelectedExample] = useState<string>("Orion Nebula M42")
  const [selectedMetal, setSelectedMetal] = useState<string>("Ti")
  const [materialType, setMaterialType] = useState<string>("Oxido")
  
  const [astroData, setAstroData] = useState<AstrophysicalData | null>(null)
  const [physicalProps, setPhysicalProps] = useState<PhysicalProperties | null>(null)
  const [recipe, setRecipe] = useState<Recipe | null>(null)
  const [isLoading, setIsLoading] = useState(false)
  
  // QE Generator state
  const [selectedStructure, setSelectedStructure] = useState<string>("TiO2_rutile")
  const [qeFiles, setQeFiles] = useState<QEFiles | null>(null)
  const [qeStructure, setQeStructure] = useState<any>(null)
  const [materialStatus, setMaterialStatus] = useState<MaterialStatus | null>(null)
  const [qeParams, setQeParams] = useState({
    ecutwfc: 60,
    kpoints: 6,
    calculationType: 'scf'
  })
  const [isQELoading, setIsQELoading] = useState(false)
  const [copied, setCopied] = useState<string | null>(null)
  
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
          materialType
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
  }, [astroData, selectedMetal, materialType])

  const generateQEFiles = useCallback(async () => {
    setIsQELoading(true)
    try {
      const response = await fetch('/api/cosmicforge', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          action: 'qe_generate',
          structureId: selectedStructure,
          astroData,
          metal: selectedMetal,
          params: {
            ecutwfc: qeParams.ecutwfc,
            kpoints: [qeParams.kpoints, qeParams.kpoints, qeParams.kpoints],
            calculationType: qeParams.calculationType
          }
        })
      })
      
      const data = await response.json()
      if (data.success) {
        setQeFiles(data.qeFiles)
        setQeStructure(data.structure)
        setMaterialStatus(data.materialStatus)
      }
    } catch (error) {
      console.error('Error generating QE files:', error)
    } finally {
      setIsQELoading(false)
    }
  }, [selectedStructure, astroData, selectedMetal, qeParams])

  const downloadQEInput = useCallback(() => {
    if (!qeFiles || !qeStructure) return
    
    const blob = new Blob([qeFiles.pwInput], { type: 'text/plain' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `${qeStructure.formula}_${qeStructure.id}_qe.in`
    document.body.appendChild(a)
    a.click()
    document.body.removeChild(a)
    URL.revokeObjectURL(url)
  }, [qeFiles, qeStructure])

  const downloadAllQEFiles = useCallback(() => {
    if (!qeFiles || !qeStructure) return
    
    // Download input file
    const inputBlob = new Blob([qeFiles.pwInput], { type: 'text/plain' })
    const inputUrl = URL.createObjectURL(inputBlob)
    const a1 = document.createElement('a')
    a1.href = inputUrl
    a1.download = `pw_scf.in`
    document.body.appendChild(a1)
    a1.click()
    document.body.removeChild(a1)
    URL.revokeObjectURL(inputUrl)
    
    // Download run script
    setTimeout(() => {
      const scriptBlob = new Blob([qeFiles.runScript], { type: 'text/plain' })
      const scriptUrl = URL.createObjectURL(scriptBlob)
      const a2 = document.createElement('a')
      a2.href = scriptUrl
      a2.download = `run_qe.sh`
      document.body.appendChild(a2)
      a2.click()
      document.body.removeChild(a2)
      URL.revokeObjectURL(scriptUrl)
    }, 100)
  }, [qeFiles, qeStructure])

  const copyToClipboard = useCallback((text: string, key: string) => {
    navigator.clipboard.writeText(text)
    setCopied(key)
    setTimeout(() => setCopied(null), 2000)
  }, [])

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
          recipe
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
  }, [astroData, physicalProps, recipe])

  const clearData = useCallback(() => {
    setAstroData(null)
    setPhysicalProps(null)
    setRecipe(null)
    setQeFiles(null)
    setQeStructure(null)
    setMaterialStatus(null)
  }, [])

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
                <p className="text-sm text-slate-500">Diseño de Materiales + Quantum ESPRESSO para nanoHUB</p>
              </div>
            </div>
            <Badge variant="outline" className="text-sm px-3 py-1">
              v3.2 nanoHUB Edition
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
                  Configuración
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
                  <Tabs defaultValue="nanohub" className="space-y-4">
                    <TabsList className="grid w-full grid-cols-5">
                      <TabsTrigger value="nanohub" className="flex items-center gap-1">
                        <Play className="w-4 h-4" />
                        nanoHUB/QE
                      </TabsTrigger>
                      <TabsTrigger value="properties">Propiedades</TabsTrigger>
                      <TabsTrigger value="synthesis">Síntesis</TabsTrigger>
                      <TabsTrigger value="files">Archivos</TabsTrigger>
                      <TabsTrigger value="pdf">PDF</TabsTrigger>
                    </TabsList>

                    {/* nanoHUB/QE Tab - NEW */}
                    <TabsContent value="nanohub">
                      <div className="space-y-4">
                        {/* Structure Selection */}
                        <Card className="shadow-sm border-slate-200">
                          <CardHeader>
                            <CardTitle className="flex items-center gap-2">
                              <Atom className="w-5 h-5 text-purple-500" />
                              Estructura Cristalina
                            </CardTitle>
                            <CardDescription>
                              Selecciona la estructura para Quantum ESPRESSO
                            </CardDescription>
                          </CardHeader>
                          <CardContent className="space-y-4">
                            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                              <div>
                                <Label>Estructura Predefinida</Label>
                                <Select value={selectedStructure} onValueChange={setSelectedStructure}>
                                  <SelectTrigger>
                                    <SelectValue />
                                  </SelectTrigger>
                                  <SelectContent>
                                    {CRYSTAL_STRUCTURES_LIST.map(s => (
                                      <SelectItem key={s.id} value={s.id}>
                                        {s.name} ({s.formula}) - {s.nat} átomos
                                      </SelectItem>
                                    ))}
                                  </SelectContent>
                                </Select>
                              </div>
                              <div>
                                <Label>Tipo de Cálculo</Label>
                                <Select 
                                  value={qeParams.calculationType} 
                                  onValueChange={(v) => setQeParams({...qeParams, calculationType: v})}
                                >
                                  <SelectTrigger>
                                    <SelectValue />
                                  </SelectTrigger>
                                  <SelectContent>
                                    <SelectItem value="scf">SCF (Energía)</SelectItem>
                                    <SelectItem value="relax">Relax (Geometría)</SelectItem>
                                    <SelectItem value="vc-relax">VC-Relax (Celda)</SelectItem>
                                  </SelectContent>
                                </Select>
                              </div>
                            </div>
                            
                            <div className="grid grid-cols-2 gap-4">
                              <div>
                                <Label>Energía de corte (Ry): {qeParams.ecutwfc}</Label>
                                <Slider 
                                  value={[qeParams.ecutwfc]} 
                                  onValueChange={(v) => setQeParams({...qeParams, ecutwfc: v[0]})} 
                                  min={30} max={100} step={5} 
                                />
                              </div>
                              <div>
                                <Label>Puntos K: {qeParams.kpoints}x{qeParams.kpoints}x{qeParams.kpoints}</Label>
                                <Slider 
                                  value={[qeParams.kpoints]} 
                                  onValueChange={(v) => setQeParams({...qeParams, kpoints: v[0]})} 
                                  min={2} max={12} step={1} 
                                />
                              </div>
                            </div>
                            
                            <Button 
                              onClick={generateQEFiles} 
                              className="w-full" 
                              size="lg"
                              disabled={isQELoading}
                            >
                              {isQELoading ? (
                                <>
                                  <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                  Generando...
                                </>
                              ) : (
                                <>
                                  <Code className="w-4 h-4 mr-2" />
                                  Generar Archivos QE para nanoHUB
                                </>
                              )}
                            </Button>
                          </CardContent>
                        </Card>

                        {/* QE Output */}
                        {qeFiles && qeStructure && (
                          <>
                            {/* Material Status */}
                            {materialStatus && (
                              <Card className={`shadow-sm border-2 ${
                                materialStatus.status === 'EXISTENTE' ? 'border-green-400 bg-green-50' :
                                materialStatus.status === 'FUTURO' ? 'border-yellow-400 bg-yellow-50' :
                                'border-red-400 bg-red-50'
                              }`}>
                                <CardContent className="p-4">
                                  <div className="flex items-center justify-between">
                                    <div className="flex items-center gap-3">
                                      <span className="text-3xl">{materialStatus.color}</span>
                                      <div>
                                        <h4 className="font-bold text-lg">{materialStatus.status}</h4>
                                        <p className="text-sm text-slate-600">{materialStatus.description}</p>
                                      </div>
                                    </div>
                                    <div className="text-right">
                                      <p className="text-sm font-medium">{qeStructure.formula}</p>
                                      <p className="text-xs text-slate-500">{qeStructure.system} | {qeStructure.nat} átomos</p>
                                    </div>
                                  </div>
                                  {materialStatus.sources.length > 0 && (
                                    <div className="mt-3 flex flex-wrap gap-2">
                                      {materialStatus.sources.map((s, i) => (
                                        <Badge key={i} variant="outline">{s}</Badge>
                                      ))}
                                    </div>
                                  )}
                                </CardContent>
                              </Card>
                            )}

                            {/* Input File */}
                            <Card className="shadow-sm border-slate-200">
                              <CardHeader className="pb-2">
                                <div className="flex items-center justify-between">
                                  <CardTitle className="text-lg flex items-center gap-2">
                                    <FileText className="w-5 h-5 text-blue-500" />
                                    pw_scf.in (Input QE)
                                  </CardTitle>
                                  <div className="flex gap-2">
                                    <Button 
                                      variant="outline" 
                                      size="sm"
                                      onClick={() => copyToClipboard(qeFiles.pwInput, 'input')}
                                    >
                                      {copied === 'input' ? (
                                        <CheckCircle className="w-4 h-4 mr-1" />
                                      ) : (
                                        <Copy className="w-4 h-4 mr-1" />
                                      )}
                                      Copiar
                                    </Button>
                                    <Button 
                                      variant="outline" 
                                      size="sm"
                                      onClick={downloadQEInput}
                                    >
                                      <Download className="w-4 h-4 mr-1" />
                                      Descargar
                                    </Button>
                                  </div>
                                </div>
                              </CardHeader>
                              <CardContent>
                                <pre className="bg-slate-900 text-green-400 p-4 rounded-lg overflow-x-auto text-xs max-h-96 overflow-y-auto">
                                  {qeFiles.pwInput}
                                </pre>
                              </CardContent>
                            </Card>

                            {/* Run Script */}
                            <Card className="shadow-sm border-slate-200">
                              <CardHeader className="pb-2">
                                <div className="flex items-center justify-between">
                                  <CardTitle className="text-lg flex items-center gap-2">
                                    <Play className="w-5 h-5 text-green-500" />
                                    run_qe.sh (Script)
                                  </CardTitle>
                                  <Button 
                                    variant="outline" 
                                    size="sm"
                                    onClick={() => copyToClipboard(qeFiles.runScript, 'script')}
                                  >
                                    {copied === 'script' ? (
                                      <CheckCircle className="w-4 h-4 mr-1" />
                                    ) : (
                                      <Copy className="w-4 h-4 mr-1" />
                                    )}
                                    Copiar
                                  </Button>
                                </div>
                              </CardHeader>
                              <CardContent>
                                <pre className="bg-slate-900 text-yellow-400 p-4 rounded-lg overflow-x-auto text-xs max-h-64 overflow-y-auto">
                                  {qeFiles.runScript}
                                </pre>
                              </CardContent>
                            </Card>

                            {/* Download All */}
                            <div className="flex gap-4">
                              <Button onClick={downloadAllQEFiles} className="flex-1" size="lg">
                                <Download className="w-4 h-4 mr-2" />
                                Descargar Todos los Archivos
                              </Button>
                            </div>

                            {/* Instructions */}
                            <Card className="shadow-sm border-slate-200 bg-blue-50">
                              <CardContent className="p-4">
                                <h4 className="font-semibold mb-2 flex items-center gap-2">
                                  <BookOpen className="w-4 h-4" />
                                  Cómo usar en nanoHUB
                                </h4>
                                <ol className="space-y-2 text-sm">
                                  <li>1. Ve a <a href="https://nanohub.org/tools/qe" target="_blank" rel="noopener noreferrer" className="text-blue-600 underline">nanoHUB Quantum ESPRESSO</a></li>
                                  <li>2. Crea una nueva simulación</li>
                                  <li>3. Copia el contenido de <code className="bg-slate-200 px-1 rounded">pw_scf.in</code> al editor</li>
                                  <li>4. O sube el archivo descargado directamente</li>
                                  <li>5. Ejecuta la simulación</li>
                                </ol>
                              </CardContent>
                            </Card>
                          </>
                        )}
                      </div>
                    </TabsContent>

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
                                Materials Project
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

                    {/* Files Tab */}
                    <TabsContent value="files">
                      <Card className="shadow-sm border-slate-200">
                        <CardHeader>
                          <CardTitle>Archivos de Simulación</CardTitle>
                          <CardDescription>Plantillas para Quantum ESPRESSO, VASP y LAMMPS</CardDescription>
                        </CardHeader>
                        <CardContent>
                          <Accordion type="single" collapsible className="space-y-2">
                            <AccordionItem value="qe" className="border rounded-lg">
                              <AccordionTrigger className="px-4">Quantum ESPRESSO Input (Plantilla)</AccordionTrigger>
                              <AccordionContent className="px-4">
                                <pre className="bg-slate-900 text-green-400 p-4 rounded-lg overflow-x-auto text-sm">
{`&CONTROL
  calculation = 'scf'
  prefix = '${selectedMetal}O2'
  outdir = './tmp'
/

&SYSTEM
  ibrav = 0
  nat = 6
  ntyp = 2
  ecutwfc = 60.0
  ecutrho = 480.0
/

&ELECTRONS
  conv_thr = 1.0d-8
/

ATOMIC_SPECIES
 Ti  47.867  Ti.pbe-spn-kjpaw_psl.1.0.0.UPF
 O  15.999  O.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
 Ti  0.000  0.000  0.000
 Ti  0.500  0.500  0.500
 O  0.305  0.305  0.000
 O  0.695  0.695  0.000
 O  0.805  0.195  0.500
 O  0.195  0.805  0.500

CELL_PARAMETERS angstrom
  4.594  0.000  0.000
  0.000  4.594  0.000
  0.000  0.000  2.959

K_POINTS automatic
  6 6 6 0 0 0`}
                                </pre>
                              </AccordionContent>
                            </AccordionItem>

                            <AccordionItem value="vasp" className="border rounded-lg">
                              <AccordionTrigger className="px-4">VASP POSCAR</AccordionTrigger>
                              <AccordionContent className="px-4">
                                <pre className="bg-slate-900 text-green-400 p-4 rounded-lg overflow-x-auto text-sm">
{`${selectedMetal}O2 - CosmicForge Lab
1.0
   4.594  0.000  0.000
   0.000  4.594  0.000
   0.000  0.000  2.959
${selectedMetal} O
2 4
Direct
  0.000  0.000  0.000 ${selectedMetal}
  0.500  0.500  0.500 ${selectedMetal}
  0.305  0.305  0.000 O
  0.695  0.695  0.000 O
  0.805  0.195  0.500 O
  0.195  0.805  0.500 O`}
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
                          <CardTitle>Generar Reporte</CardTitle>
                        </CardHeader>
                        <CardContent className="space-y-4">
                          <Button onClick={downloadPDF} size="lg" className="w-full">
                            <Download className="w-4 h-4 mr-2" />
                            Descargar Reporte TXT
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
                      manuales para comenzar a diseñar materiales y generar archivos para nanoHUB.
                    </p>
                    <div className="flex items-center justify-center gap-4 text-sm text-slate-400">
                      <span className="flex items-center gap-1">
                        <Database className="w-4 h-4" />
                        {Object.keys(ASTROPHYSICAL_EXAMPLES).length} ejemplos
                      </span>
                      <span className="flex items-center gap-1">
                        <Play className="w-4 h-4" />
                        nanoHUB
                      </span>
                    </div>
                  </div>
                </CardContent>
              </Card>
            )}
          </div>
        </div>

        {/* Footer */}
        <footer className="mt-8 py-6 border-t border-slate-200 text-center text-sm text-slate-500">
          <p className="italic mb-2">"{randomQuote}"</p>
          <p>— Yuri I. Manin, "Lo demostrable e indemostrable"</p>
          <p className="mt-4">CosmicForge Lab v3.2 | Diseño de Materiales + Quantum ESPRESSO para nanoHUB</p>
        </footer>
      </main>
    </div>
  )
}
