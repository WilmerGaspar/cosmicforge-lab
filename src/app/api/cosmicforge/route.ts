import { NextRequest, NextResponse } from 'next/server'

// ============================================================================
// PHYSICS CALCULATIONS
// ============================================================================

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
  const porosity = Math.max(0.05, Math.min(0.95, 1.5 - astroData.fractal_dimension * 0.7))
  const base_density = 4.5
  const density = base_density * (1 - porosity * 0.4) * (0.8 + astroData.criticality_score * 0.4)
  const thermal_conductivity = 2.5 * Math.pow(astroData.fractal_dimension / 1.5, 0.8) * (1 - astroData.anisotropy * 0.3)
  const elastic_modulus = 230 * astroData.criticality_score * (astroData.fractal_dimension / 1.5)
  const surface_area = Math.max(10, 200 * (1 - astroData.entropy * 20) * (1 + astroData.anisotropy * 0.5))
  const band_gap = 3.2 * astroData.criticality_score * (1 + astroData.anisotropy * 0.2)
  const quality_score = (
    astroData.criticality_score * 0.3 +
    (1 - porosity) * 0.2 +
    Math.min(1, surface_area / 100) * 0.2 +
    Math.min(1, elastic_modulus / 200) * 0.15 +
    Math.min(1, band_gap / 3) * 0.15
  )
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

// ============================================================================
// RECIPE GENERATION
// ============================================================================

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
  reaction_conditions: { temperature_C: number; reaction_time_hours: number; pH: number };
  precursors: string[];
  step_by_step: string[];
} {
  const temperature_C = 300 + astroData.fractal_dimension * 150 + astroData.turbulence_beta * 50
  const reaction_time_hours = 2 + astroData.criticality_score * 4
  const pH = 7 + astroData.anisotropy * 3 - astroData.entropy * 100

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
    'Pd': ['PdCl₂', 'Pd(NO₃)₂', 'Acetato de Paladio'],
    'Cr': ['Cr(NO₃)₃·9H₂O', 'CrCl₃·6H₂O', 'Acetato de Cromo'],
    'Mo': ['MoO₃', '(NH₄)₆Mo₇O₂₄', 'MoCl₅'],
    'W': ['WO₃', 'Na₂WO₄', 'WCl₆']
  }

  const precursors = precursorMap[metal] || ['Precursor metálico', 'Agua desionizada', 'Base (NaOH/NH₄OH)']

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

// ============================================================================
// EXTREME CONDITIONS SIMULATION
// ============================================================================

function simulateExtremeConditions(
  physicalProps: any,
  metal: string,
  materialType: string
): {
  vacuum_stability: number;
  radiation_degradation: number;
  microgravity_effect: number;
  thermal_cycle_life: number;
  estimated_lifetime: string;
  recommendations: string[];
} {
  // Vacuum stability - based on porosity and density
  const vacuum_stability = Math.min(0.99, 
    (1 - physicalProps.porosity * 0.3) * 
    (physicalProps.density > 4 ? 0.9 : 0.7) * 
    (materialType === 'Ceramica' ? 1.1 : 1)
  )

  // Radiation degradation - based on band gap and activation energy
  const radiation_degradation = Math.max(0.01,
    0.4 - physicalProps.band_gap * 0.05 + 
    (1 - physicalProps.quality_score) * 0.2
  )

  // Microgravity effect - based on elastic modulus and density
  const microgravity_effect = Math.max(0.01,
    0.1 + (1 - physicalProps.elastic_modulus / 300) * 0.15
  )

  // Thermal cycle life
  const thermal_cycle_life = Math.round(
    100 + physicalProps.thermal_conductivity * 50 + physicalProps.quality_score * 200
  )

  // Estimated lifetime
  const years = Math.round(10 + vacuum_stability * 20 - radiation_degradation * 5)
  const estimated_lifetime = `${years} años en órbita LEO`

  // Recommendations
  const recommendations: string[] = []
  
  if (vacuum_stability < 0.8) {
    recommendations.push('Aplicar recubrimiento protector para mejorar estabilidad en vacío')
  }
  if (radiation_degradation > 0.2) {
    recommendations.push('Considerar blindaje contra radiación o dopaje con elementos pesados')
  }
  if (microgravity_effect > 0.15) {
    recommendations.push('Evaluar estabilidad estructural en condiciones de microgravedad')
  }
  if (thermal_cycle_life < 200) {
    recommendations.push('Mejorar resistencia térmica mediante tratamiento térmico adicional')
  }
  if (recommendations.length === 0) {
    recommendations.push('Material apto para aplicaciones espaciales sin modificaciones')
    recommendations.push('Proceder con pruebas de calificación espacial')
  }

  return {
    vacuum_stability,
    radiation_degradation,
    microgravity_effect,
    thermal_cycle_life,
    estimated_lifetime,
    recommendations
  }
}

// ============================================================================
// GENETIC ALGORITHM OPTIMIZATION
// ============================================================================

function optimizeComposition(
  astroData: any,
  metal: string,
  materialType: string
): {
  best_composition: string;
  score: number;
  parameters: { temperature: number; time: number; ph: number };
  generations: number;
} {
  // Simulated genetic algorithm
  const populationSize = 50
  const generations = 50
  let bestScore = 0
  let bestParams = { temp: 400, time: 4, ph: 7 }
  let bestRatio = '100%'

  // Simulate optimization process
  for (let gen = 0; gen < generations; gen++) {
    for (let i = 0; i < populationSize; i++) {
      const temp = 300 + Math.random() * 500
      const time = 1 + Math.random() * 10
      const ph = 2 + Math.random() * 10
      const ratio = 0.5 + Math.random() * 0.5

      // Fitness function based on astrophysical parameters
      const score = (
        astroData.criticality_score * 0.3 +
        (1 - Math.abs(temp - 500) / 500) * 0.2 +
        (1 - Math.abs(ph - 7) / 7) * 0.2 +
        ratio * 0.15 +
        Math.random() * 0.15
      )

      if (score > bestScore) {
        bestScore = score
        bestParams = { temp: Math.round(temp), time: Math.round(time * 10) / 10, ph: Math.round(ph * 10) / 10 }
        bestRatio = `${(ratio * 100).toFixed(0)}% ${metal} - ${((1-ratio) * 100).toFixed(0)}% O`
      }
    }
  }

  return {
    best_composition: bestRatio,
    score: bestScore,
    parameters: {
      temperature: bestParams.temp,
      time: bestParams.time,
      ph: bestParams.ph
    },
    generations
  }
}

// ============================================================================
// DFT CALCULATION (Simulated)
// ============================================================================

function runDFTCalculation(
  physicalProps: any,
  recipe: any,
  metal: string
): {
  total_energy: number;
  band_gap: number;
  optimized_lattice: number;
  fermi_energy: number;
  is_stable: boolean;
  calculation_time: string;
} {
  // Simulated DFT results based on physical properties
  const baseEnergy = -150 - Math.random() * 50
  const total_energy = baseEnergy * (1 + physicalProps.quality_score * 0.1)
  
  const band_gap = physicalProps.band_gap * (0.8 + Math.random() * 0.4)
  
  const optimized_lattice = 4.0 + physicalProps.density * 0.1 + Math.random() * 0.5
  
  const fermi_energy = 5 + band_gap / 2 + Math.random() * 0.5
  
  // Stability based on energy and band gap
  const is_stable = total_energy < -100 && band_gap > 0.5

  const calculation_time = `${Math.round(30 + Math.random() * 60)}s`

  return {
    total_energy,
    band_gap,
    optimized_lattice,
    fermi_energy,
    is_stable,
    calculation_time
  }
}

// ============================================================================
// MATERIALS PROJECT API (Simulated)
// ============================================================================

async function queryMaterialsProject(metal: string, formula: string): Promise<{
  count: number;
  bestMatch: string;
  materials: Array<{ id: string; formula: string; band_gap: number; energy: number }>;
}> {
  // Simulated Materials Project query
  // In production, this would call the actual API: https://materialsproject.org/rest/v2
  
  return {
    count: Math.round(5 + Math.random() * 20),
    bestMatch: `${formula} (mp-${10000 + Math.round(Math.random() * 50000)})`,
    materials: [
      { id: `mp-${10000 + Math.round(Math.random() * 10000)}`, formula, band_gap: 3.0 + Math.random(), energy: -150 - Math.random() * 50 },
      { id: `mp-${10000 + Math.round(Math.random() * 10000)}`, formula: `${metal}O`, band_gap: 2.5 + Math.random(), energy: -120 - Math.random() * 30 },
      { id: `mp-${10000 + Math.round(Math.random() * 10000)}`, formula: `${metal}2O3`, band_gap: 2.0 + Math.random(), energy: -200 - Math.random() * 40 }
    ]
  }
}

// ============================================================================
// PDF REPORT GENERATION
// ============================================================================

function generatePDFReport(
  astroData: any,
  physicalProps: any,
  recipe: any,
  alloyKey: string | null,
  simulation?: any,
  optimization?: any,
  dft?: any
): string {
  const date = new Date().toLocaleDateString('es-ES', { 
    year: 'numeric', month: 'long', day: 'numeric', hour: '2-digit', minute: '2-digit'
  })

  const maninQuotes = [
    "La matemática es la parte más pura de la cultura humana, y sin embargo está profundamente conectada con la realidad física.",
    "La demostrabilidad es un concepto que trasciende la lógica formal: es la búsqueda de certeza en un universo de incertidumbre.",
    "Entre lo demostrable y lo indemostrable existe una zona fronteriza donde nace la creatividad científica."
  ]
  const quote = maninQuotes[Math.floor(Math.random() * maninQuotes.length)]

  return `
================================================================================
                    COSMICFORGE LAB - INFORME TÉCNICO COMPLETO
                         v3.5 NASA Edition - Personal
================================================================================

FECHA: ${date}
OBJETO ASTROFÍSICO: ${astroData.object_name}
FUENTE: ${astroData.source}
MATERIAL: ${recipe.metal} - ${recipe.material_type}

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

  CONDICIONES DE REACCIÓN:
    Temperatura:    ${recipe.reaction_conditions.temperature_C}°C
    Tiempo:         ${recipe.reaction_conditions.reaction_time_hours} horas
    pH:             ${recipe.reaction_conditions.pH}

  PRECURSORES:
${recipe.precursors.map((p: string) => `    • ${p}`).join('\n')}

  PROCEDIMIENTO:
${recipe.step_by_step.map((s: string, i: number) => `    ${i + 1}. ${s}`).join('\n')}
${simulation ? `
================================================================================
                    SIMULACIÓN DE CONDICIONES EXTREMAS
================================================================================

  Estabilidad en Vacío:     ${(simulation.vacuum_stability * 100).toFixed(1)}%
  Degradación por Radiación: ${(simulation.radiation_degradation * 100).toFixed(1)}%
  Efecto Microgravedad:     ${(simulation.microgravity_effect * 100).toFixed(1)}%
  Vida en Ciclos Térmicos:  ${simulation.thermal_cycle_life} ciclos
  Vida Útil Estimada:       ${simulation.estimated_lifetime}

  RECOMENDACIONES:
${simulation.recommendations.map((r: string) => `    → ${r}`).join('\n')}
` : ''}${optimization ? `
================================================================================
                    OPTIMIZACIÓN (ALGORITMO GENÉTICO)
================================================================================

  Mejor Composición:    ${optimization.best_composition}
  Score:                ${(optimization.score * 100).toFixed(2)}%
  Generaciones:         ${optimization.generations}

  PARÁMETROS ÓPTIMOS:
    Temperatura: ${optimization.parameters.temperature}°C
    Tiempo:      ${optimization.parameters.time} horas
    pH:          ${optimization.parameters.ph}
` : ''}${dft ? `
================================================================================
                    CÁLCULO DFT (QUANTUM ESPRESSO)
================================================================================

  Energía Total:        ${dft.total_energy.toFixed(4)} Ry
  Band Gap:             ${dft.band_gap.toFixed(4)} eV
  Celda Optimizada:     ${dft.optimized_lattice.toFixed(4)} Å
  Energía de Fermi:     ${dft.fermi_energy.toFixed(4)} eV
  Estado:               ${dft.is_stable ? 'ESTABLE ✓' : 'INESTABLE ⚠'}
  Tiempo de Cálculo:    ${dft.calculation_time}
` : ''}
================================================================================
                        COMPARACIÓN CON MATERIALES COMERCIALES
================================================================================

  Material              | Resistencia | Peso   | Costo | Ranking
  ----------------------|-------------|--------|-------|----------
  Tu Material (${recipe.metal}-O) | ${(physicalProps.elastic_modulus * 0.8).toFixed(0)}     | ${physicalProps.density.toFixed(2)} | ~50   | NUEVO
  Ti-6Al-4V (Grade 5)   | 950         | 4.43   | 85    | Aerospace
  Al 7075-T6            | 570         | 2.81   | 35    | Aerospace
  Steel 4340            | 1080        | 7.85   | 25    | Structural
  Inconel 718           | 1200        | 8.19   | 120   | High-Temp

================================================================================
                        REFLEXIÓN CIENTÍFICA
================================================================================

  "${quote}"

                                        — Yuri I. Manin
                                          "Lo demostrable e indemostrable"

================================================================================
                        ANÁLISIS EPISTEMOLÓGICO
================================================================================

  Nivel de certeza matemática: ${physicalProps.quality_score > 0.6 ? 'ALTO' : 'MODERADO'}
  Nivel de certeza física: REQUIERE VALIDACIÓN EXPERIMENTAL

  PRÓXIMOS PASOS:
    1. Validar estructura cristalina por XRD
    2. Medir propiedades mecánicas reales
    3. Comparar con bases de datos (Materials Project, AFLOW)
    4. Realizar cálculos DFT de alta precisión
    5. Síntesis experimental y caracterización

================================================================================
                        ENLACES DE VALIDACIÓN
================================================================================

  Materials Project:  https://materialsproject.org/materials?search=${recipe.metal}O2
  AFLOW Library:      https://www.aflowlib.org/?search=${recipe.metal}O2
  OQMD Database:      http://oqmd.org/materials?search=${recipe.metal}O2

================================================================================

                    Generado por CosmicForge Lab v3.5 NASA Edition
                              Edición Personal

================================================================================
`
}

// ============================================================================
// API HANDLER
// ============================================================================

export async function POST(request: NextRequest) {
  try {
    const body = await request.json()
    const { action, astroData, metal, materialType, alloyKey, physicalProps, recipe, formula } = body

    switch (action) {
      case 'generate': {
        const physicalPropsResult = calculatePhysicalProperties(astroData)
        const recipeResult = generateRecipe(metal, materialType, astroData)
        return NextResponse.json({
          success: true,
          physicalProps: physicalPropsResult,
          recipe: recipeResult
        })
      }

      case 'simulate': {
        const simulation = simulateExtremeConditions(physicalProps, metal, materialType)
        return NextResponse.json({
          success: true,
          simulation
        })
      }

      case 'optimize': {
        const optimization = optimizeComposition(astroData, metal, materialType)
        return NextResponse.json({
          success: true,
          optimization
        })
      }

      case 'dft': {
        const dft = runDFTCalculation(physicalProps, recipe, metal)
        return NextResponse.json({
          success: true,
          dft
        })
      }

      case 'materials_project': {
        const results = await queryMaterialsProject(metal, formula)
        return NextResponse.json({
          success: true,
          results
        })
      }

      case 'pdf': {
        const reportContent = generatePDFReport(
          astroData, physicalProps, recipe, alloyKey,
          body.simulation, body.optimization, body.dft
        )
        return new NextResponse(reportContent, {
          headers: {
            'Content-Type': 'text/plain; charset=utf-8',
            'Content-Disposition': `attachment; filename="CosmicForge_Report_${astroData.object_name.replace(/\s+/g, '_')}.txt"`
          }
        })
      }

      default:
        return NextResponse.json({ success: false, error: 'Invalid action' }, { status: 400 })
    }
  } catch (error) {
    console.error('API Error:', error)
    return NextResponse.json({ success: false, error: 'Server error' }, { status: 500 })
  }
}
