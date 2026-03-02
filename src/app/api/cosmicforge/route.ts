import { NextRequest, NextResponse } from 'next/server'
import {
  CRYSTAL_STRUCTURES,
  PSEUDOPOTENTIALS,
  generateQEInput,
  generateRunScript,
  generateReadme,
  nebulaToStructure,
  astroToDFTParams,
  getMaterialStatus,
  DEFAULT_PARAMS,
  type CrystalStructure,
  type SimulationParams
} from '@/lib/qe-generator'

// Physics calculation functions
function calculatePhysicalProperties(astroData: {
  fractal_dimension: number;
  criticality_score: number;
  entropy: number;
  anisotropy: number;
  turbulence_beta: number;
  lyapunov_max: number;
}): {
  porosity: number;
  density: number;
  thermal_conductivity: number;
  elastic_modulus: number;
  surface_area: number;
  band_gap: number;
  quality_score: number;
  activation_energy: number;
} {
  // Porosity: inversely related to fractal dimension
  const porosity = Math.max(0.05, Math.min(0.95, 1.5 - astroData.fractal_dimension * 0.7))
  
  // Density: related to criticality and inversely to porosity
  const base_density = 4.5 // g/cm³ typical for TiO2
  const density = base_density * (1 - porosity * 0.4) * (0.8 + astroData.criticality_score * 0.4)
  
  // Thermal conductivity: related to fractal dimension and anisotropy
  const thermal_conductivity = 2.5 * Math.pow(astroData.fractal_dimension / 1.5, 0.8) * (1 - astroData.anisotropy * 0.3)
  
  // Elastic modulus: related to criticality and fractal dimension
  const elastic_modulus = 230 * astroData.criticality_score * (astroData.fractal_dimension / 1.5)
  
  // Surface area: inversely related to entropy
  const surface_area = Math.max(10, 200 * (1 - astroData.entropy * 20) * (1 + astroData.anisotropy * 0.5))
  
  // Band gap: related to criticality and anisotropy
  const band_gap = 3.2 * astroData.criticality_score * (1 + astroData.anisotropy * 0.2)
  
  // Quality score: overall assessment
  const quality_score = (
    astroData.criticality_score * 0.3 +
    (1 - porosity) * 0.2 +
    Math.min(1, surface_area / 100) * 0.2 +
    Math.min(1, elastic_modulus / 200) * 0.15 +
    Math.min(1, band_gap / 3) * 0.15
  )
  
  // Activation energy
  const activation_energy = 0.5 + Math.abs(astroData.lyapunov_max) * 0.3 + astroData.turbulence_beta * 0.02
  
  return {
    porosity,
    density: Math.max(1, density),
    thermal_conductivity: Math.max(0.1, thermal_conductivity),
    elastic_modulus: Math.max(10, elastic_modulus),
    surface_area: Math.max(5, surface_area),
    band_gap: Math.max(0.1, band_gap),
    quality_score: Math.min(1, quality_score),
    activation_energy: Math.max(0.1, activation_energy)
  }
}

// Recipe generation
function generateRecipe(
  metal: string,
  materialType: string,
  astroData: {
    fractal_dimension: number;
    criticality_score: number;
    entropy: number;
    anisotropy: number;
    turbulence_beta: number;
  }
): {
  metal: string;
  material_type: string;
  reaction_conditions: {
    temperature_C: number;
    reaction_time_hours: number;
    pH: number;
  };
  precursors: string[];
  step_by_step: string[];
} {
  // Calculate reaction conditions based on astrophysical parameters
  const temperature_C = 300 + astroData.fractal_dimension * 150 + astroData.turbulence_beta * 50
  const reaction_time_hours = 2 + astroData.criticality_score * 4
  const pH = 7 + astroData.anisotropy * 3 - astroData.entropy * 100
  
  // Precursors based on metal
  const precursorMap: Record<string, string[]> = {
    'Ti': ['Isopropóxido de Titanio (IV)', 'TiCl₄', 'Nitrato de Titanio'],
    'Al': ['Isopropóxido de Aluminio', 'Al(NO₃)₃·9H₂O', 'AlCl₃'],
    'Fe': ['FeCl₃·6H₂O', 'Fe(NO₃)₃·9H₂O', 'Sulfato de Hierro (III)'],
    'Zn': ['Zn(NO₃)₂·6H₂O', 'Acetato de Zinc', 'ZnCl₂'],
    'Cu': ['Cu(NO₃)₂·3H₂O', 'Sulfato de Cobre (II)', 'Acetato de Cobre'],
    'Ni': ['Ni(NO₃)₂·6H₂O', 'NiCl₂·6H₂O', 'Acetato de Níquel'],
    'Co': ['Co(NO₃)₂·6H₂O', 'CoCl₂·6H₂O', 'Acetato de Cobalto'],
    'Mn': ['Mn(NO₃)₂·4H₂O', 'MnCl₂·4H₂O', 'Acetato de Manganeso'],
    'Ag': ['AgNO₃', 'Acetato de Plata', 'Nitrato de Plata'],
    'Au': ['HAuCl₄·3H₂O', 'Cloruro de Oro (III)', 'Acetato de Oro'],
    'Pt': ['H₂PtCl₆·6H₂O', 'Cloroplatinato', 'Acetilacetonato de Pt'],
    'Pd': ['PdCl₂', 'Pd(NO₃)₂', 'Acetato de Paladio']
  }
  
  const precursors = precursorMap[metal] || ['Precursor metálico', 'Agua desionizada', 'Base (NaOH/NH₄OH)']
  
  // Step by step procedure
  const step_by_step = [
    `Preparar solución 0.5M del precursor de ${metal} en agua desionizada`,
    `Ajustar pH a ${pH.toFixed(1)} con NaOH o NH₄OH según el caso`,
    `Calentar la solución a ${temperature_C.toFixed(0)}°C con agitación constante`,
    `Mantener la temperatura por ${reaction_time_hours.toFixed(1)} horas`,
    'Enfriar lentamente a temperatura ambiente',
    'Filtrar el precipitado y lavar con agua/etanol (3x)',
    'Secar en estufa a 100°C por 12 horas',
    `Calcinación en horno mufla a ${(temperature_C + 100).toFixed(0)}°C por 4 horas`,
    'Caracterizar por XRD, SEM y BET'
  ]
  
  return {
    metal,
    material_type: materialType,
    reaction_conditions: {
      temperature_C: Math.round(temperature_C),
      reaction_time_hours: Math.round(reaction_time_hours * 10) / 10,
      pH: Math.round(pH * 10) / 10
    },
    precursors,
    step_by_step
  }
}

// PDF Report Generation
function generatePDFReport(
  astroData: any,
  physicalProps: any,
  recipe: any,
  alloyKey: string | null
): string {
  const date = new Date().toLocaleDateString('es-ES', { 
    year: 'numeric', 
    month: 'long', 
    day: 'numeric',
    hour: '2-digit',
    minute: '2-digit'
  })
  
  const maninQuotes = [
    "La matemática es la parte más pura de la cultura humana, y sin embargo está profundamente conectada con la realidad física.",
    "La demostrabilidad es un concepto que trasciende la lógica formal: es la búsqueda de certeza en un universo de incertidumbre.",
    "Entre lo demostrable y lo indemostrable existe una zona fronteriza donde nace la creatividad científica."
  ]
  const quote = maninQuotes[Math.floor(Math.random() * maninQuotes.length)]
  
  const report = `
================================================================================
                        COSMICFORGE LAB - INFORME TÉCNICO
                              v3.1 Personal Edition
================================================================================

FECHA: ${date}
OBJETO ASTROFÍSICO: ${astroData.object_name}
FUENTE: ${astroData.source}

================================================================================
                           FIRMA ASTROFÍSICA
================================================================================

  Dimensión Fractal:     ${astroData.fractal_dimension.toFixed(6)}
  Criticalidad:          ${astroData.criticality_score.toFixed(6)}
  Entropía:              ${astroData.entropy.toFixed(9)}
  Anisotropía:           ${astroData.anisotropy.toFixed(6)}
  Turbulencia β:         ${astroData.turbulence_beta.toFixed(6)}
  Lyapunov máximo:       ${astroData.lyapunov_max.toFixed(6)}
  Modo:                  ${astroData.mode}

================================================================================
                        PROPIEDADES DEL MATERIAL
================================================================================

  Porosidad:             ${(physicalProps.porosity * 100).toFixed(2)}%
  Densidad:              ${physicalProps.density.toFixed(4)} g/cm³
  Conductividad Térmica: ${physicalProps.thermal_conductivity.toFixed(4)} W/mK
  Módulo Elástico:       ${physicalProps.elastic_modulus.toFixed(4)} GPa
  Área Superficial:      ${physicalProps.surface_area.toFixed(2)} m²/g
  Band Gap:              ${physicalProps.band_gap.toFixed(4)} eV
  Calidad:               ${(physicalProps.quality_score * 100).toFixed(2)}%
  Energía de Activación: ${physicalProps.activation_energy.toFixed(6)} eV

================================================================================
                        RECETA DE SÍNTESIS
================================================================================

  Material: ${recipe.metal} - ${recipe.material_type}

  CONDICIONES DE REACCIÓN:
    Temperatura:    ${recipe.reaction_conditions.temperature_C}°C
    Tiempo:         ${recipe.reaction_conditions.reaction_time_hours} horas
    pH:             ${recipe.reaction_conditions.pH}

  PRECURSORES:
${recipe.precursors.map((p: string) => `    • ${p}`).join('\n')}

  PROCEDIMIENTO:
${recipe.step_by_step.map((s: string, i: number) => `    ${i + 1}. ${s}`).join('\n')}

================================================================================
                        REFLEXIÓN CIENTÍFICA
================================================================================

  "${quote}"

                                        — Yuri I. Manin
                                          "Lo demostrable e indemostrable"

================================================================================
                        ANÁLISIS EPISTEMOLÓGICO
================================================================================

  El material predicho a partir de ${astroData.object_name} representa una
  hipótesis científica derivada de analogías matemáticas entre estructuras
  cósmicas y propiedades materiales.

  Nivel de certeza matemática: ${physicalProps.quality_score > 0.6 ? 'ALTO' : 'MODERADO'}
  - Las ecuaciones están formalmente derivadas.

  Nivel de certeza física: REQUIERE VALIDACIÓN
  - La traducción astrofísica a material es una hipótesis por verificar.

  PRÓXIMOS PASOS RECOMENDADOS:
    1. Simulación DFT con VASP o Quantum ESPRESSO
    2. Validación experimental de propiedades
    3. Comparación con bases de datos (Materials Project, AFLOW)

================================================================================
                        ENLACES DE VALIDACIÓN
================================================================================

  Materials Project:  https://materialsproject.org/materials?search=${recipe.metal}O2
  AFLOW Library:      https://www.aflowlib.org/?search=${recipe.metal}O2
  OQMD Database:      http://oqmd.org/materials?search=${recipe.metal}O2

================================================================================

                          Generado por CosmicForge Lab
                               Edición Personal

================================================================================
`
  return report
}

export async function POST(request: NextRequest) {
  try {
    const body = await request.json()
    const { action, astroData, metal, materialType, alloyKey, physicalProps, recipe, structureId, params } = body

    if (action === 'generate') {
      // Calculate physical properties
      const physicalPropsResult = calculatePhysicalProperties(astroData)
      
      // Generate recipe
      const recipeResult = generateRecipe(metal, materialType, astroData)
      
      return NextResponse.json({
        success: true,
        physicalProps: physicalPropsResult,
        recipe: recipeResult
      })
    }
    
    if (action === 'pdf') {
      // Generate PDF report
      const reportContent = generatePDFReport(astroData, physicalProps, recipe, alloyKey)
      
      // Return as downloadable text file
      return new NextResponse(reportContent, {
        headers: {
          'Content-Type': 'text/plain; charset=utf-8',
          'Content-Disposition': `attachment; filename="CosmicForge_Report_${astroData.object_name.replace(/\s+/g, '_')}.txt"`
        }
      })
    }

    // New: Generate QE input files for nanoHUB
    if (action === 'qe_generate') {
      const targetStructureId = structureId || (astroData ? nebulaToStructure(astroData.type || 'emission', metal) : 'TiO2_rutile')
      const structure = CRYSTAL_STRUCTURES[targetStructureId]
      
      if (!structure) {
        return NextResponse.json({ 
          success: false, 
          error: `Unknown structure: ${targetStructureId}` 
        }, { status: 400 })
      }

      // Calculate DFT params from astro data if available
      let simParams: SimulationParams = { ...DEFAULT_PARAMS }
      if (astroData) {
        const astroParams = astroToDFTParams(astroData)
        simParams = { ...DEFAULT_PARAMS, ...astroParams }
      }
      
      // Override with user params if provided
      if (params) {
        simParams = { ...simParams, ...params }
      }
      
      // Set title
      simParams.title = astroData?.object_name 
        ? `${structure.formula} - ${astroData.object_name}`
        : `${structure.formula} DFT Calculation`
      
      // Generate files
      const qeInput = generateQEInput(structure, simParams)
      const runScript = generateRunScript(structure, simParams)
      const readme = generateReadme(structure, simParams)
      const status = getMaterialStatus(structure.formula)
      
      return NextResponse.json({
        success: true,
        qeFiles: {
          pwInput: qeInput,
          runScript: runScript,
          readme: readme
        },
        structure: {
          id: structure.id,
          name: structure.name,
          formula: structure.formula,
          system: structure.system,
          spaceGroup: structure.spaceGroup,
          nat: structure.atoms.length,
          latticeParameter: structure.latticeParameter
        },
        params: simParams,
        materialStatus: status
      })
    }

    // New: Get available structures
    if (action === 'get_structures') {
      const structures = Object.values(CRYSTAL_STRUCTURES).map(s => ({
        id: s.id,
        name: s.name,
        formula: s.formula,
        system: s.system,
        nat: s.atoms.length
      }))
      
      return NextResponse.json({
        success: true,
        structures
      })
    }

    // New: Download QE files as ZIP
    if (action === 'qe_download') {
      const targetStructureId = structureId || 'TiO2_rutile'
      const structure = CRYSTAL_STRUCTURES[targetStructureId]
      
      if (!structure) {
        return NextResponse.json({ 
          success: false, 
          error: `Unknown structure: ${targetStructureId}` 
        }, { status: 400 })
      }

      let simParams: SimulationParams = { ...DEFAULT_PARAMS }
      if (astroData) {
        const astroParams = astroToDFTParams(astroData)
        simParams = { ...DEFAULT_PARAMS, ...astroParams }
      }
      if (params) {
        simParams = { ...simParams, ...params }
      }
      simParams.title = astroData?.object_name 
        ? `${structure.formula} - ${astroData.object_name}`
        : `${structure.formula} DFT Calculation`
      
      const qeInput = generateQEInput(structure, simParams)
      const runScript = generateRunScript(structure, simParams)
      const readme = generateReadme(structure, simParams)
      
      // Return individual file for direct download
      const filename = `${structure.formula}_${targetStructureId}_qe.in`
      
      return new NextResponse(qeInput, {
        headers: {
          'Content-Type': 'text/plain; charset=utf-8',
          'Content-Disposition': `attachment; filename="${filename}"`
        }
      })
    }

    return NextResponse.json({ success: false, error: 'Invalid action' }, { status: 400 })
  } catch (error) {
    console.error('API Error:', error)
    return NextResponse.json({ success: false, error: 'Server error' }, { status: 500 })
  }
}
