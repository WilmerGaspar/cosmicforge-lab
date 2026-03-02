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
import { 
  Sparkles, FlaskConical, Atom, Download, Database, 
  Thermometer, Gauge, Beaker, CheckCircle2,
  ExternalLink, Loader2, Settings, Play, Copy, 
  CheckCircle
} from 'lucide-react'
import {
  NANOHUB_STRUCTURES,
  generateNanoHUBFormat,
  getNanoHUBStructureList,
  type NanoHUBStructure
} from '@/lib/nanohub-generator'

// Astrophysical Examples
const ASTROPHYSICAL_EXAMPLES = {
  "Orion Nebula M42": {
    fractal_dimension: 1.654, criticality_score: 0.722, entropy: 0.019,
    anisotropy: 0.329, turbulence_beta: 2.278, lyapunov_max: -0.227,
    mode: "balanced", type: "emission", distance: "1,344 ly"
  },
  "Crab Nebula M1": {
    fractal_dimension: 1.82, criticality_score: 0.89, entropy: 0.032,
    anisotropy: 0.456, turbulence_beta: 2.45, lyapunov_max: -0.15,
    mode: "turbulent", type: "supernova", distance: "6,500 ly"
  },
  "Ring Nebula M57": {
    fractal_dimension: 1.45, criticality_score: 0.55, entropy: 0.012,
    anisotropy: 0.28, turbulence_beta: 1.95, lyapunov_max: -0.35,
    mode: "stable", type: "planetary", distance: "2,283 ly"
  }
}

// Elements Database
const ELEMENTS_DATABASE: Record<string, { name: string; mw: number }> = {
  "Ti": { name: "Titanio", mw: 47.867 },
  "Al": { name: "Aluminio", mw: 26.982 },
  "Fe": { name: "Hierro", mw: 55.845 },
  "Zn": { name: "Zinc", mw: 65.38 },
  "Si": { name: "Silicio", mw: 28.086 }
}

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

export default function CosmicForgeLab() {
  const [selectedExample, setSelectedExample] = useState<string>("Orion Nebula M42")
  const [selectedMetal, setSelectedMetal] = useState<string>("Ti")
  
  const [astroData, setAstroData] = useState<AstrophysicalData | null>(null)
  const [isLoading, setIsLoading] = useState(false)
  
  // nanoHUB state
  const [selectedStructure, setSelectedStructure] = useState<string>("TiO2_rutile")
  const [nanohubFormat, setNanohubFormat] = useState<ReturnType<typeof generateNanoHUBFormat> | null>(null)
  const [copied, setCopied] = useState<string | null>(null)
  
  const [inputMethod, setInputMethod] = useState<"examples" | "manual">("examples")
  const [manualName, setManualName] = useState("Objeto Personalizado")
  const [manualFd, setManualFd] = useState([1.5])
  const [manualCs, setManualCs] = useState([0.5])

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

  const generateForNanoHUB = useCallback(() => {
    const structure = NANOHUB_STRUCTURES[selectedStructure]
    if (structure) {
      const format = generateNanoHUBFormat(structure, astroData?.object_name)
      setNanohubFormat(format)
    }
  }, [selectedStructure, astroData])

  const copyToClipboard = useCallback((text: string, key: string) => {
    navigator.clipboard.writeText(text)
    setCopied(key)
    setTimeout(() => setCopied(null), 2000)
  }, [])

  const generateRecipe = useCallback(async () => {
    if (!astroData) return
    setIsLoading(true)
    // Simular cálculo
    setTimeout(() => {
      setIsLoading(false)
      generateForNanoHUB()
    }, 1000)
  }, [astroData, generateForNanoHUB])

  const clearData = useCallback(() => {
    setAstroData(null)
    setNanohubFormat(null)
  }, [])

  const structureList = getNanoHUBStructureList()
  const currentStructure = NANOHUB_STRUCTURES[selectedStructure]

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
                <p className="text-sm text-slate-500">Generador para nanoHUB Quantum ESPRESSO</p>
              </div>
            </div>
            <Badge variant="outline" className="text-sm px-3 py-1">
              v3.3 nanoHUB Edition
            </Badge>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-4 py-6">
        <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
          {/* Sidebar */}
          <aside className="lg:col-span-1 space-y-4">
            {/* Input Selection */}
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
                      <Label>Dim. Fractal: {manualFd[0].toFixed(3)}</Label>
                      <Slider value={manualFd} onValueChange={setManualFd} min={0} max={3} step={0.001} />
                    </div>
                    
                    <Button onClick={() => {
                      setAstroData({
                        object_name: manualName,
                        fractal_dimension: manualFd[0],
                        criticality_score: manualCs[0],
                        entropy: 0.01,
                        anisotropy: 0.3,
                        turbulence_beta: 2.0,
                        lyapunov_max: -0.1,
                        mode: 'custom',
                        source: 'Editor'
                      })
                    }} className="w-full" variant="default">
                      Aplicar
                    </Button>
                  </TabsContent>
                </Tabs>
              </CardContent>
            </Card>

            {/* Structure Selection */}
            <Card className="shadow-sm border-slate-200">
              <CardHeader className="pb-3">
                <CardTitle className="text-lg flex items-center gap-2">
                  <Atom className="w-5 h-5 text-purple-500" />
                  Estructura Cristalina
                </CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <div>
                  <Label>Seleccionar Estructura</Label>
                  <Select value={selectedStructure} onValueChange={setSelectedStructure}>
                    <SelectTrigger>
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      {structureList.map(s => (
                        <SelectItem key={s.id} value={s.id}>
                          {s.name} ({s.formula}) - {s.nat} átomos
                        </SelectItem>
                      ))}
                    </SelectContent>
                  </Select>
                </div>
                
                {currentStructure && (
                  <div className="p-3 bg-purple-50 rounded-lg text-sm">
                    <p><strong>Fórmula:</strong> {currentStructure.formula}</p>
                    <p><strong>Tipo:</strong> {currentStructure.structureType}</p>
                    <p><strong>Átomos:</strong> {currentStructure.atoms.length}</p>
                    <p><strong>a:</strong> {currentStructure.latticeParameter} Å</p>
                  </div>
                )}
              </CardContent>
            </Card>

            {astroData && (
              <Button onClick={clearData} variant="outline" className="w-full">
                Limpiar
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
                      Firma Astrofísica: {astroData.object_name}
                    </CardTitle>
                  </CardHeader>
                  <CardContent>
                    <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mb-4">
                      <div className="p-3 bg-slate-50 rounded-lg">
                        <p className="text-xs text-slate-500">Dim. Fractal</p>
                        <p className="text-xl font-bold">{astroData.fractal_dimension.toFixed(4)}</p>
                      </div>
                      <div className="p-3 bg-slate-50 rounded-lg">
                        <p className="text-xs text-slate-500">Criticalidad</p>
                        <p className="text-xl font-bold">{astroData.criticality_score.toFixed(4)}</p>
                      </div>
                      <div className="p-3 bg-slate-50 rounded-lg">
                        <p className="text-xs text-slate-500">Tipo</p>
                        <p className="text-xl font-bold">{astroData.type || 'N/A'}</p>
                      </div>
                      <div className="p-3 bg-slate-50 rounded-lg">
                        <p className="text-xs text-slate-500">Distancia</p>
                        <p className="text-xl font-bold">{astroData.distance || 'N/A'}</p>
                      </div>
                    </div>
                    
                    <Button onClick={generateForNanoHUB} className="w-full" size="lg">
                      <Play className="w-4 h-4 mr-2" />
                      Generar Formato para nanoHUB
                    </Button>
                  </CardContent>
                </Card>

                {/* nanoHUB Format Output */}
                {nanohubFormat && (
                  <Card className="shadow-sm border-green-200 bg-green-50">
                    <CardHeader>
                      <CardTitle className="flex items-center gap-2 text-green-700">
                        <CheckCircle2 className="w-5 h-5" />
                        Formato para nanoHUB
                      </CardTitle>
                      <CardDescription>
                        Copia cada campo a la interfaz de nanoHUB
                      </CardDescription>
                    </CardHeader>
                    <CardContent className="space-y-4">
                      {/* Title */}
                      <div className="space-y-2">
                        <div className="flex items-center justify-between">
                          <Label className="font-semibold">Title of Run</Label>
                          <Button 
                            variant="ghost" 
                            size="sm"
                            onClick={() => copyToClipboard(nanohubFormat.title, 'title')}
                          >
                            {copied === 'title' ? <CheckCircle className="w-4 h-4 text-green-500" /> : <Copy className="w-4 h-4" />}
                          </Button>
                        </div>
                        <Input readOnly value={nanohubFormat.title} className="bg-white" />
                      </div>

                      {/* Atomic Structure */}
                      <div className="space-y-2">
                        <div className="flex items-center justify-between">
                          <Label className="font-semibold">Atomic Structure (copiar todo)</Label>
                          <Button 
                            variant="ghost" 
                            size="sm"
                            onClick={() => copyToClipboard(nanohubFormat.atomicStructure, 'atoms')}
                          >
                            {copied === 'atoms' ? <CheckCircle className="w-4 h-4 text-green-500" /> : <Copy className="w-4 h-4" />}
                          </Button>
                        </div>
                        <pre className="bg-white p-3 rounded-lg border text-sm overflow-x-auto">
{nanohubFormat.atomicStructure}
                        </pre>
                      </div>

                      {/* Cell Vectors */}
                      <div className="space-y-2">
                        <div className="flex items-center justify-between">
                          <Label className="font-semibold">Cell Vectors (Å)</Label>
                          <Button 
                            variant="ghost" 
                            size="sm"
                            onClick={() => copyToClipboard(nanohubFormat.cellVectors, 'vectors')}
                          >
                            {copied === 'vectors' ? <CheckCircle className="w-4 h-4 text-green-500" /> : <Copy className="w-4 h-4" />}
                          </Button>
                        </div>
                        <pre className="bg-white p-3 rounded-lg border text-sm">
{nanohubFormat.cellVectors}
                        </pre>
                      </div>

                      {/* Lattice Parameters */}
                      <div className="grid grid-cols-3 gap-4">
                        <div className="space-y-2">
                          <Label className="font-semibold">Lattice Parameter a (Å)</Label>
                          <Input readOnly value={nanohubFormat.latticeParameterA} className="bg-white" />
                        </div>
                        <div className="space-y-2">
                          <Label className="font-semibold">Ratio b/a</Label>
                          <Input readOnly value={nanohubFormat.ratioBA} className="bg-white" />
                        </div>
                        <div className="space-y-2">
                          <Label className="font-semibold">Ratio c/a</Label>
                          <Input readOnly value={nanohubFormat.ratioCA} className="bg-white" />
                        </div>
                      </div>

                      {/* Structure Type */}
                      <div className="grid grid-cols-2 gap-4">
                        <div className="space-y-2">
                          <Label className="font-semibold">Structure Type</Label>
                          <Input readOnly value={nanohubFormat.structureType} className="bg-white" />
                        </div>
                        <div className="space-y-2">
                          <Label className="font-semibold">Atomic Coordinates</Label>
                          <Input readOnly value={nanohubFormat.atomicCoordinates} className="bg-white" />
                        </div>
                      </div>

                      <Separator />

                      {/* Instructions */}
                      <div className="p-4 bg-blue-100 rounded-lg">
                        <h4 className="font-semibold mb-2 flex items-center gap-2">
                          <ExternalLink className="w-4 h-4" />
                          Cómo usar en nanoHUB:
                        </h4>
                        <ol className="space-y-1 text-sm list-decimal list-inside">
                          <li>Ve a <a href="https://nanohub.org/tools/qe" target="_blank" rel="noopener noreferrer" className="text-blue-600 underline">nanoHUB QE</a></li>
                          <li>En "Input Geometry", selecciona <strong>"User defined"</strong></li>
                          <li>En "Atomic Coordinates" selecciona <strong>"Fractional"</strong></li>
                          <li>Copia el <strong>Title of Run</strong></li>
                          <li>Copia todo el campo <strong>Atomic Structure</strong></li>
                          <li>Copia los <strong>Cell Vectors</strong></li>
                          <li>Ingresa el <strong>Lattice Parameter a</strong></li>
                          <li>Ingresa los <strong>ratios b/a</strong> y <strong>c/a</strong></li>
                          <li>Ejecuta la simulación</li>
                        </ol>
                      </div>
                    </CardContent>
                  </Card>
                )}
              </>
            ) : (
              <Card className="shadow-sm border-slate-200">
                <CardContent className="py-16 text-center">
                  <div className="max-w-md mx-auto">
                    <Atom className="w-16 h-16 mx-auto text-blue-500 mb-4" />
                    <h2 className="text-xl font-semibold mb-2">CosmicForge Lab</h2>
                    <p className="text-slate-500 mb-6">
                      Carga un objeto astrofísico para generar el formato de nanoHUB.
                    </p>
                    <p className="text-sm text-slate-400">
                      O selecciona una estructura cristalina y haz clic en "Generar Formato para nanoHUB"
                    </p>
                    
                    <Button onClick={generateForNanoHUB} className="mt-4" variant="outline">
                      <Play className="w-4 h-4 mr-2" />
                      Generar sin datos astrofísicos
                    </Button>
                  </div>
                </CardContent>
              </Card>
            )}
          </div>
        </div>

        {/* Footer */}
        <footer className="mt-8 py-4 text-center text-sm text-slate-500">
          CosmicForge Lab v3.3 | Generador de formato para nanoHUB Quantum ESPRESSO
        </footer>
      </main>
    </div>
  )
}
