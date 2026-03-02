/**
 * CosmicForge Lab - Generador Quantum ESPRESSO v3.0
 * ==================================================
 * Genera inputs QE 100% compatibles con nanoHUB.
 * 
 * FORMATO nanoHUB:
 * - Título descriptivo
 * - Estructura atómica con coordenadas fraccionales
 * - Vectores de celda en Ångstroms
 * - Parámetros de red correctos
 */

// Estructuras cristalinas predefinidas
export interface CrystalStructure {
  id: string;
  name: string;
  formula: string;
  system: 'cubic' | 'tetragonal' | 'hexagonal' | 'trigonal' | 'orthorhombic' | 'monoclinic' | 'triclinic';
  spaceGroup: string;
  latticeParameter: number; // en Å
  ratioBA?: number;
  ratioCA?: number;
  cellVectors: [[number, number, number], [number, number, number], [number, number, number]];
  atoms: Array<{
    element: string;
    x: number;
    y: number;
    z: number;
  }>;
  description: string;
}

// Estructuras cristalinas válidas para nanoHUB
export const CRYSTAL_STRUCTURES: Record<string, CrystalStructure> = {
  // Silicio diamante - Referencia nanoHUB
  "Si_diamond": {
    id: "Si_diamond",
    name: "Silicio Diamond",
    formula: "Si",
    system: "cubic",
    spaceGroup: "Fd-3m",
    latticeParameter: 5.43,
    ratioBA: 1,
    ratioCA: 1,
    cellVectors: [
      [-2.715, 2.715, 0.0],
      [-2.715, 0.0, 2.715],
      [0.0, 2.715, 2.715]
    ],
    atoms: [
      { element: "Si", x: 0.0, y: 0.0, z: 0.0 },
      { element: "Si", x: 0.25, y: 0.25, z: 0.25 }
    ],
    description: "Silicon diamond structure"
  },
  
  // TiO2 Rutilo - CORREGIDO
  "TiO2_rutile": {
    id: "TiO2_rutile",
    name: "TiO2 Rutilo",
    formula: "TiO2",
    system: "tetragonal",
    spaceGroup: "P4_2/mnm",
    latticeParameter: 4.594, // a en Å
    ratioBA: 1,
    ratioCA: 0.644, // c/a
    cellVectors: [
      [4.594, 0.0, 0.0],
      [0.0, 4.594, 0.0],
      [0.0, 0.0, 2.959]
    ],
    atoms: [
      { element: "Ti", x: 0.0, y: 0.0, z: 0.0 },
      { element: "Ti", x: 0.5, y: 0.5, z: 0.5 },
      { element: "O", x: 0.3053, y: 0.3053, z: 0.0 },
      { element: "O", x: 0.6947, y: 0.6947, z: 0.0 },
      { element: "O", x: 0.8053, y: 0.1947, z: 0.5 },
      { element: "O", x: 0.1947, y: 0.8053, z: 0.5 }
    ],
    description: "TiO2 rutile structure"
  },
  
  // TiO2 Anatasa
  "TiO2_anatase": {
    id: "TiO2_anatase",
    name: "TiO2 Anatasa",
    formula: "TiO2",
    system: "tetragonal",
    spaceGroup: "I4_1/amd",
    latticeParameter: 3.785, // a en Å
    ratioBA: 1,
    ratioCA: 2.514, // c/a
    cellVectors: [
      [3.785, 0.0, 0.0],
      [0.0, 3.785, 0.0],
      [0.0, 0.0, 9.514]
    ],
    atoms: [
      { element: "Ti", x: 0.0, y: 0.0, z: 0.0 },
      { element: "Ti", x: 0.5, y: 0.5, z: 0.5 },
      { element: "Ti", x: 0.0, y: 0.5, z: 0.25 },
      { element: "Ti", x: 0.5, y: 0.0, z: 0.75 },
      { element: "O", x: 0.0, y: 0.0, z: 0.208 },
      { element: "O", x: 0.0, y: 0.0, z: 0.792 },
      { element: "O", x: 0.5, y: 0.5, z: 0.292 },
      { element: "O", x: 0.5, y: 0.5, z: 0.708 },
      { element: "O", x: 0.0, y: 0.5, z: 0.042 },
      { element: "O", x: 0.0, y: 0.5, z: 0.458 },
      { element: "O", x: 0.5, y: 0.0, z: 0.542 },
      { element: "O", x: 0.5, y: 0.0, z: 0.958 }
    ],
    description: "TiO2 anatase structure"
  },
  
  // Ti2O3 Corundum
  "Ti2O3_corundum": {
    id: "Ti2O3_corundum",
    name: "Ti2O3 Corundum",
    formula: "Ti2O3",
    system: "trigonal",
    spaceGroup: "R-3c",
    latticeParameter: 5.148, // a en Å
    ratioBA: 1,
    ratioCA: 1.476, // c/a ratio para celda hexagonal
    cellVectors: [
      [5.148, 0.0, 0.0],
      [-2.574, 4.458, 0.0],
      [0.0, 0.0, 13.642]
    ],
    atoms: [
      { element: "Ti", x: 0.0, y: 0.0, z: 0.0 },
      { element: "Ti", x: 0.0, y: 0.0, z: 0.333 },
      { element: "Ti", x: 0.333, y: 0.667, z: 0.167 },
      { element: "Ti", x: 0.667, y: 0.333, z: 0.500 },
      { element: "O", x: 0.306, y: 0.0, z: 0.083 },
      { element: "O", x: 0.0, y: 0.306, z: 0.083 },
      { element: "O", x: 0.694, y: 0.0, z: 0.250 },
      { element: "O", x: 0.0, y: 0.694, z: 0.250 },
      { element: "O", x: 0.306, y: 0.0, z: 0.417 },
      { element: "O", x: 0.0, y: 0.306, z: 0.417 }
    ],
    description: "Ti2O3 corundum structure"
  },
  
  // TiO Rocksalt
  "TiO_rocksalt": {
    id: "TiO_rocksalt",
    name: "TiO Rocksalt",
    formula: "TiO",
    system: "cubic",
    spaceGroup: "Fm-3m",
    latticeParameter: 4.177,
    ratioBA: 1,
    ratioCA: 1,
    cellVectors: [
      [4.177, 0.0, 0.0],
      [0.0, 4.177, 0.0],
      [0.0, 0.0, 4.177]
    ],
    atoms: [
      { element: "Ti", x: 0.0, y: 0.0, z: 0.0 },
      { element: "O", x: 0.5, y: 0.5, z: 0.5 }
    ],
    description: "TiO rocksalt structure"
  },
  
  // ZnO Wurtzita
  "ZnO_wurtzite": {
    id: "ZnO_wurtzite",
    name: "ZnO Wurtzita",
    formula: "ZnO",
    system: "hexagonal",
    spaceGroup: "P6_3mc",
    latticeParameter: 3.250,
    ratioBA: 1,
    ratioCA: 1.603,
    cellVectors: [
      [3.250, 0.0, 0.0],
      [-1.625, 2.814, 0.0],
      [0.0, 0.0, 5.207]
    ],
    atoms: [
      { element: "Zn", x: 0.333, y: 0.667, z: 0.0 },
      { element: "Zn", x: 0.667, y: 0.333, z: 0.5 },
      { element: "O", x: 0.333, y: 0.667, z: 0.375 },
      { element: "O", x: 0.667, y: 0.333, z: 0.875 }
    ],
    description: "ZnO wurtzite structure"
  },
  
  // Fe2O3 Hematita
  "Fe2O3_hematite": {
    id: "Fe2O3_hematite",
    name: "Fe2O3 Hematita",
    formula: "Fe2O3",
    system: "trigonal",
    spaceGroup: "R-3c",
    latticeParameter: 5.038,
    ratioBA: 1,
    ratioCA: 1.357,
    cellVectors: [
      [5.038, 0.0, 0.0],
      [-2.519, 4.363, 0.0],
      [0.0, 0.0, 13.772]
    ],
    atoms: [
      { element: "Fe", x: 0.0, y: 0.0, z: 0.0 },
      { element: "Fe", x: 0.0, y: 0.0, z: 0.333 },
      { element: "Fe", x: 0.333, y: 0.667, z: 0.167 },
      { element: "Fe", x: 0.667, y: 0.333, z: 0.500 },
      { element: "O", x: 0.306, y: 0.0, z: 0.083 },
      { element: "O", x: 0.0, y: 0.306, z: 0.083 },
      { element: "O", x: 0.694, y: 0.0, z: 0.250 },
      { element: "O", x: 0.0, y: 0.694, z: 0.250 },
      { element: "O", x: 0.306, y: 0.0, z: 0.417 },
      { element: "O", x: 0.0, y: 0.306, z: 0.417 }
    ],
    description: "Fe2O3 hematite structure"
  },
  
  // Al2O3 Corundum (Zafiro)
  "Al2O3_corundum": {
    id: "Al2O3_corundum",
    name: "Al2O3 Corundum",
    formula: "Al2O3",
    system: "trigonal",
    spaceGroup: "R-3c",
    latticeParameter: 4.759,
    ratioBA: 1,
    ratioCA: 1.306,
    cellVectors: [
      [4.759, 0.0, 0.0],
      [-2.380, 4.122, 0.0],
      [0.0, 0.0, 12.991]
    ],
    atoms: [
      { element: "Al", x: 0.0, y: 0.0, z: 0.0 },
      { element: "Al", x: 0.0, y: 0.0, z: 0.333 },
      { element: "Al", x: 0.333, y: 0.667, z: 0.167 },
      { element: "Al", x: 0.667, y: 0.333, z: 0.500 },
      { element: "O", x: 0.306, y: 0.0, z: 0.083 },
      { element: "O", x: 0.0, y: 0.306, z: 0.083 },
      { element: "O", x: 0.694, y: 0.0, z: 0.250 },
      { element: "O", x: 0.0, y: 0.694, z: 0.250 },
      { element: "O", x: 0.306, y: 0.0, z: 0.417 },
      { element: "O", x: 0.0, y: 0.306, z: 0.417 }
    ],
    description: "Al2O3 corundum structure"
  }
};

// Pseudopotenciales PSL
export const PSEUDOPOTENTIALS: Record<string, { mass: number; filename: string; url: string }> = {
  "Si": { mass: 28.086, filename: "Si.pbe-n-kjpaw_psl.1.0.0.UPF", url: "https://www.quantum-espresso.org/upf_files/Si.pbe-n-kjpaw_psl.1.0.0.UPF" },
  "Ti": { mass: 47.867, filename: "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF", url: "https://www.quantum-espresso.org/upf_files/Ti.pbe-spn-kjpaw_psl.1.0.0.UPF" },
  "O": { mass: 15.999, filename: "O.pbe-n-kjpaw_psl.1.0.0.UPF", url: "https://www.quantum-espresso.org/upf_files/O.pbe-n-kjpaw_psl.1.0.0.UPF" },
  "Zn": { mass: 65.38, filename: "Zn.pbe-dn-kjpaw_psl.1.0.0.UPF", url: "https://www.quantum-espresso.org/upf_files/Zn.pbe-dn-kjpaw_psl.1.0.0.UPF" },
  "Fe": { mass: 55.845, filename: "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF", url: "https://www.quantum-espresso.org/upf_files/Fe.pbe-spn-kjpaw_psl.1.0.0.UPF" },
  "Al": { mass: 26.982, filename: "Al.pbe-n-kjpaw_psl.1.0.0.UPF", url: "https://www.quantum-espresso.org/upf_files/Al.pbe-n-kjpaw_psl.1.0.0.UPF" },
  "Cu": { mass: 63.546, filename: "Cu.pbe-dn-kjpaw_psl.1.0.0.UPF", url: "https://www.quantum-espresso.org/upf_files/Cu.pbe-dn-kjpaw_psl.1.0.0.UPF" },
  "Ni": { mass: 58.693, filename: "Ni.pbe-dn-kjpaw_psl.1.0.0.UPF", url: "https://www.quantum-espresso.org/upf_files/Ni.pbe-dn-kjpaw_psl.1.0.0.UPF" },
  "Co": { mass: 58.933, filename: "Co.pbe-dn-kjpaw_psl.1.0.0.UPF", url: "https://www.quantum-espresso.org/upf_files/Co.pbe-dn-kjpaw_psl.1.0.0.UPF" },
  "Mn": { mass: 54.938, filename: "Mn.pbe-dn-kjpaw_psl.1.0.0.UPF", url: "https://www.quantum-espresso.org/upf_files/Mn.pbe-dn-kjpaw_psl.1.0.0.UPF" }
};

// Parámetros de simulación
export interface SimulationParams {
  title: string;
  structureId: string;
  calculationType: 'scf' | 'relax' | 'vc-relax' | 'bands' | 'dos';
  ecutwfc: number; // Energía de corte en Ry
  ecutrho: number; // Densidad de corte en Ry
  kpoints: [number, number, number]; // Puntos k
  degauss: number; // Smearing
  convThr: number; // Convergencia
  mixingBeta: number;
  electronMaxstep: number;
}

// Parámetros por defecto
export const DEFAULT_PARAMS: SimulationParams = {
  title: "DFT Calculation",
  structureId: "TiO2_rutile",
  calculationType: 'scf',
  ecutwfc: 60,
  ecutrho: 480,
  kpoints: [6, 6, 6],
  degauss: 0.02,
  convThr: 1.0e-8,
  mixingBeta: 0.7,
  electronMaxstep: 100
};

/**
 * Genera el archivo pw_scf.in para Quantum ESPRESSO
 * Formato compatible con nanoHUB
 */
export function generateQEInput(structure: CrystalStructure, params: SimulationParams): string {
  const lines: string[] = [];
  
  // &CONTROL
  lines.push("&CONTROL");
  lines.push(`  calculation = '${params.calculationType}'`);
  lines.push(`  title = '${params.title}'`);
  lines.push("  prefix = 'material'");
  lines.push("  outdir = './tmp'");
  lines.push("  pseudo_dir = './pseudo'");
  lines.push("  tprnfor = .true.");
  lines.push("  tstress = .true.");
  lines.push("/");
  lines.push("");
  
  // &SYSTEM
  lines.push("&SYSTEM");
  lines.push("  ibrav = 0"); // Usamos CELL_PARAMETERS
  lines.push(`  nat = ${structure.atoms.length}`);
  
  // Contar tipos de átomos
  const uniqueElements = [...new Set(structure.atoms.map(a => a.element))];
  lines.push(`  ntyp = ${uniqueElements.length}`);
  lines.push(`  ecutwfc = ${params.ecutwfc}`);
  lines.push(`  ecutrho = ${params.ecutrho}`);
  lines.push("  occupations = 'smearing'");
  lines.push("  smearing = 'mp'");
  lines.push(`  degauss = ${params.degauss}`);
  lines.push("/");
  lines.push("");
  
  // &ELECTRONS
  lines.push("&ELECTRONS");
  lines.push(`  conv_thr = ${params.convThr.toExponential(2)}`);
  lines.push(`  mixing_beta = ${params.mixingBeta}`);
  lines.push(`  electron_maxstep = ${params.electronMaxstep}`);
  lines.push("/");
  lines.push("");
  
  // &IONS (para relax)
  if (params.calculationType === 'relax' || params.calculationType === 'vc-relax') {
    lines.push("&IONS");
    lines.push("  ion_dynamics = 'bfgs'");
    lines.push("/");
    lines.push("");
    
    if (params.calculationType === 'vc-relax') {
      lines.push("&CELL");
      lines.push("  cell_dynamics = 'bfgs'");
      lines.push("  press = 0.0");
      lines.push("/");
      lines.push("");
    }
  }
  
  // ATOMIC_SPECIES
  lines.push("ATOMIC_SPECIES");
  for (const elem of uniqueElements.sort()) {
    const pseudo = PSEUDOPOTENTIALS[elem];
    if (pseudo) {
      lines.push(`  ${elem}  ${pseudo.mass.toFixed(3)}  ${pseudo.filename}`);
    } else {
      lines.push(`  ${elem}  50.000  ${elem}.UPF`);
    }
  }
  lines.push("");
  
  // ATOMIC_POSITIONS (crystal = coordenadas fraccionales)
  lines.push("ATOMIC_POSITIONS crystal");
  for (const atom of structure.atoms) {
    lines.push(`  ${atom.element}  ${atom.x.toFixed(6)}  ${atom.y.toFixed(6)}  ${atom.z.toFixed(6)}`);
  }
  lines.push("");
  
  // CELL_PARAMETERS (en Ångstroms)
  lines.push("CELL_PARAMETERS angstrom");
  for (const vec of structure.cellVectors) {
    lines.push(`  ${vec[0].toFixed(6)}  ${vec[1].toFixed(6)}  ${vec[2].toFixed(6)}`);
  }
  lines.push("");
  
  // K_POINTS
  lines.push("K_POINTS automatic");
  lines.push(`  ${params.kpoints[0]} ${params.kpoints[1]} ${params.kpoints[2]}  0 0 0`);
  
  return lines.join("\n");
}

/**
 * Genera script de ejecución para nanoHUB
 */
export function generateRunScript(structure: CrystalStructure, params: SimulationParams): string {
  const uniqueElements = [...new Set(structure.atoms.map(a => a.element))];
  const pseudoDownloads = uniqueElements
    .map(elem => {
      const pseudo = PSEUDOPOTENTIALS[elem];
      if (pseudo) {
        return `wget -q "${pseudo.url}" -O ./pseudo/${pseudo.filename}`;
      }
      return null;
    })
    .filter(Boolean)
    .join("\n");

  return `#!/bin/bash
# ============================================================================
# CosmicForge Lab - Quantum ESPRESSO Script
# Material: ${structure.formula} (${structure.name})
# Calculation: ${params.calculationType}
# ============================================================================

echo "=========================================="
echo "CosmicForge Lab - Quantum ESPRESSO"
echo "=========================================="
echo "Material: ${structure.formula} (${structure.name})"
echo "Structure: ${structure.system} (${structure.spaceGroup})"
echo "Calculation: ${params.calculationType}"
echo "=========================================="

# Create directories
mkdir -p ./tmp ./pseudo

# Download pseudopotentials
echo "Downloading pseudopotentials..."
${pseudoDownloads}

# Run calculation
echo ""
echo "Running ${params.calculationType} calculation..."
pw.x < pw_scf.in > pw_scf.out

# Check convergence
if grep -q "convergence has been achieved" pw_scf.out; then
    echo "OK: Calculation converged successfully!"
    
    # Extract total energy
    ENERGY=$(grep "!.*total energy" pw_scf.out | tail -1 | awk '{print $5}')
    echo "Total energy: $ENERGY Ry"
    
    # Extract optimized cell if relax
    if [ "${params.calculationType}" = "relax" ] || [ "${params.calculationType}" = "vc-relax" ]; then
        echo ""
        echo "Optimized structure:"
        grep "CELL_PARAMETERS" -A 3 pw_scf.out | tail -3
        grep "ATOMIC_POSITIONS" -A ${structure.atoms.length} pw_scf.out | tail -${structure.atoms.length}
    fi
    
    echo ""
    echo "=========================================="
    echo "Calculation completed successfully!"
    echo "=========================================="
else
    echo "ERROR: Calculation did not converge"
    echo ""
    echo "Suggestions:"
    echo "  1. Reduce mixing_beta to 0.3-0.5"
    echo "  2. Increase electron_maxstep to 200"
    echo "  3. Reduce conv_thr to 1.0e-6"
    echo "  4. Increase ecutwfc by 10-20 Ry"
    echo "  5. Use mixing_mode = 'local-TF'"
    exit 1
fi
`;
}

/**
 * Genera README con información del material
 */
export function generateReadme(structure: CrystalStructure, params: SimulationParams): string {
  return `# CosmicForge Lab - Material Generated

## Crystal Structure

| Property | Value |
|----------|-------|
| Formula | ${structure.formula} |
| Name | ${structure.name} |
| System | ${structure.system} |
| Space Group | ${structure.spaceGroup} |
| Lattice Parameter (a) | ${structure.latticeParameter} Å |
${structure.ratioBA ? `| Ratio b/a | ${structure.ratioBA} |` : ''}
${structure.ratioCA ? `| Ratio c/a | ${structure.ratioCA} |` : ''}

## Atomic Positions (Fractional)

| Element | x | y | z |
|---------|---|---|---|
${structure.atoms.map(a => `| ${a.element} | ${a.x.toFixed(6)} | ${a.y.toFixed(6)} | ${a.z.toFixed(6)} |`).join('\n')}

## Cell Vectors (Å)

| Vector | x | y | z |
|--------|---|---|---|
${structure.cellVectors.map((v, i) => `| ${['a1', 'a2', 'a3'][i]} | ${v[0].toFixed(6)} | ${v[1].toFixed(6)} | ${v[2].toFixed(6)} |`).join('\n')}

## DFT Parameters

| Parameter | Value |
|-----------|-------|
| Calculation | ${params.calculationType} |
| Cutoff (ecutwfc) | ${params.ecutwfc} Ry |
| Density (ecutrho) | ${params.ecutrho} Ry |
| k-points | ${params.kpoints.join(' × ')} |
| XC Functional | PBE (GGA) |
| Convergence | ${params.convThr.toExponential(2)} Ry |

## Files Included

1. \`pw_scf.in\` - Input file for Quantum ESPRESSO
2. \`run_qe.sh\` - Execution script
3. \`README.md\` - This file

## Running on nanoHUB

\`\`\`bash
chmod +x run_qe.sh
./run_qe.sh
\`\`\`

---
*Generated by CosmicForge Lab v3.0*
`;
}

/**
 * Mapea tipo de nebulosa a estructura cristalina
 */
export function nebulaToStructure(nebulaType: string, metal: string = "Ti"): string {
  const mapping: Record<string, string> = {
    "emission": "TiO2_rutile",
    "reflection": "TiO2_anatase",
    "dark": "Ti2O3_corundum",
    "planetary": "TiO2_rutile",
    "supernova": "TiO_rocksalt",
    "spiral_galaxy": "TiO2_rutile"
  };
  
  // Si es Zn, usar ZnO
  if (metal === "Zn") return "ZnO_wurtzite";
  // Si es Fe, usar Fe2O3
  if (metal === "Fe") return "Fe2O3_hematite";
  // Si es Al, usar Al2O3
  if (metal === "Al") return "Al2O3_corundum";
  // Si es Si, usar Si_diamond
  if (metal === "Si") return "Si_diamond";
  
  return mapping[nebulaType] || "TiO2_rutile";
}

/**
 * Calcula parámetros DFT desde datos astrofísicos
 */
export function astroToDFTParams(astroData: {
  fractal_dimension: number;
  criticality_score: number;
  entropy: number;
}): Partial<SimulationParams> {
  const { fractal_dimension, criticality_score, entropy } = astroData;
  
  // Energía de corte: 45-80 Ry basado en criticalidad
  const ecutwfc = Math.round(45 + criticality_score * 35);
  
  // Puntos k: 4-9 basado en dimensión fractal
  const k = Math.round(4 + (fractal_dimension - 1) * 2.5);
  
  // Convergencia basada en entropía
  const convThr = Math.pow(10, -8 + entropy * 2);
  
  return {
    ecutwfc,
    ecutrho: ecutwfc * 4,
    kpoints: [k, k, k] as [number, number, number],
    convThr
  };
}

/**
 * Valida la fórmula química
 */
export function validateFormula(formula: string): { valid: boolean; elements: string[]; error?: string } {
  const elements: string[] = [];
  
  // Parsear fórmula (simplificado)
  const regex = /([A-Z][a-z]?)(\d*)/g;
  let match;
  
  while ((match = regex.exec(formula)) !== null) {
    if (match[1]) {
      elements.push(match[1]);
    }
  }
  
  if (elements.length === 0) {
    return { valid: false, elements: [], error: "Invalid formula format" };
  }
  
  // Verificar que los elementos existan
  for (const elem of elements) {
    if (!PSEUDOPOTENTIALS[elem]) {
      return { 
        valid: false, 
        elements, 
        error: `Unknown element: ${elem}. No pseudopotential available.` 
      };
    }
  }
  
  return { valid: true, elements };
}

/**
 * Estado del material
 */
export function getMaterialStatus(formula: string): {
  status: 'EXISTENTE' | 'FUTURO' | 'DESCONOCIDO';
  color: string;
  icon: string;
  description: string;
  sources: string[];
} {
  const existing = ["Si", "TiO2", "ZnO", "Fe2O3", "Al2O3", "Ti2O3"];
  const predicted = ["TiO", "Ti3O5"];
  
  if (existing.includes(formula)) {
    return {
      status: "EXISTENTE",
      color: "🟢",
      icon: "📚",
      description: "En bases de datos, sintetizable",
      sources: ["Materials Project", "AFLOW", "OQMD", "COD"]
    };
  } else if (predicted.includes(formula)) {
    return {
      status: "FUTURO",
      color: "🟡",
      icon: "🔮",
      description: "Predicho, requiere validación",
      sources: ["GNoME", "NOMAD"]
    };
  } else {
    return {
      status: "DESCONOCIDO",
      color: "🔴",
      icon: "❓",
      description: "Nuevo, requiere validación completa",
      sources: []
    };
  }
}
