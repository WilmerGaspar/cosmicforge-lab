/**
 * CosmicForge Lab - Generador para nanoHUB Web Interface
 * =======================================================
 * Genera el formato EXACTO que la interfaz web de nanoHUB espera.
 * 
 * CAMPOS DE nanoHUB:
 * - Title of Run
 * - Atomic Structure (descripción + coordenadas fraccionales)
 * - Cell Vectors (Å)
 * - Lattice Parameter a (Å)
 * - Ratio b/a, c/a
 */

export interface NanoHUBStructure {
  id: string;
  name: string;
  formula: string;
  description: string;
  structureType: string;
  latticeParameter: number; // a en Å
  ratioBA: number;
  ratioCA: number;
  cellVectors: [[number, number, number], [number, number, number], [number, number, number]];
  atoms: Array<{
    element: string;
    x: number;
    y: number;
    z: number;
  }>;
}

// Estructuras cristalinas con formato nanoHUB
export const NANOHUB_STRUCTURES: Record<string, NanoHUBStructure> = {
  "Si_diamond": {
    id: "Si_diamond",
    name: "Si diamond",
    formula: "Si",
    description: "Silicon diamond structure",
    structureType: "cubic F (fcc)",
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
    ]
  },
  
  "TiO2_rutile": {
    id: "TiO2_rutile",
    name: "TiO2 rutile",
    formula: "TiO2",
    description: "TiO2 rutile structure",
    structureType: "tetragonal P",
    latticeParameter: 4.594,
    ratioBA: 1,
    ratioCA: 0.644,
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
    ]
  },
  
  "TiO2_anatase": {
    id: "TiO2_anatase",
    name: "TiO2 anatase",
    formula: "TiO2",
    description: "TiO2 anatase structure",
    structureType: "tetragonal I (bct)",
    latticeParameter: 3.785,
    ratioBA: 1,
    ratioCA: 2.514,
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
    ]
  },
  
  "Ti2O3_corundum": {
    id: "Ti2O3_corundum",
    name: "Ti2O3 corundum",
    formula: "Ti2O3",
    description: "Ti2O3 corundum structure",
    structureType: "trigonal R",
    latticeParameter: 5.148,
    ratioBA: 1,
    ratioCA: 2.649,
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
    ]
  },
  
  "TiO_rocksalt": {
    id: "TiO_rocksalt",
    name: "TiO rocksalt",
    formula: "TiO",
    description: "TiO rocksalt structure",
    structureType: "cubic F (fcc)",
    latticeParameter: 4.177,
    ratioBA: 1,
    ratioCA: 1,
    cellVectors: [
      [2.089, 2.089, 0.0],
      [2.089, 0.0, 2.089],
      [0.0, 2.089, 2.089]
    ],
    atoms: [
      { element: "Ti", x: 0.0, y: 0.0, z: 0.0 },
      { element: "O", x: 0.5, y: 0.5, z: 0.5 }
    ]
  },
  
  "ZnO_wurtzite": {
    id: "ZnO_wurtzite",
    name: "ZnO wurtzite",
    formula: "ZnO",
    description: "ZnO wurtzite structure",
    structureType: "hexagonal P",
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
    ]
  },
  
  "Fe2O3_hematite": {
    id: "Fe2O3_hematite",
    name: "Fe2O3 hematite",
    formula: "Fe2O3",
    description: "Fe2O3 hematite structure",
    structureType: "trigonal R",
    latticeParameter: 5.038,
    ratioBA: 1,
    ratioCA: 2.734,
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
    ]
  },
  
  "Al2O3_corundum": {
    id: "Al2O3_corundum",
    name: "Al2O3 corundum",
    formula: "Al2O3",
    description: "Al2O3 corundum structure",
    structureType: "trigonal R",
    latticeParameter: 4.759,
    ratioBA: 1,
    ratioCA: 2.730,
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
    ]
  }
};

/**
 * Genera el contenido del campo "Atomic Structure" para nanoHUB
 * Formato: descripción + coordenadas
 */
export function generateAtomicStructure(structure: NanoHUBStructure): string {
  const lines: string[] = [];
  
  // Primera línea: descripción
  lines.push(structure.description);
  
  // Siguientes líneas: átomos con coordenadas fraccionales
  for (const atom of structure.atoms) {
    lines.push(`${atom.element} ${atom.x.toFixed(4)} ${atom.y.toFixed(4)} ${atom.z.toFixed(4)}`);
  }
  
  return lines.join('\n');
}

/**
 * Genera el contenido del campo "Cell Vectors" para nanoHUB
 */
export function generateCellVectors(structure: NanoHUBStructure): string {
  const lines: string[] = [];
  
  for (const vec of structure.cellVectors) {
    lines.push(`${vec[0].toFixed(3)} ${vec[1].toFixed(3)} ${vec[2].toFixed(3)}`);
  }
  
  return lines.join('\n');
}

/**
 * Genera el título para nanoHUB
 */
export function generateTitle(structure: NanoHUBStructure, objectName?: string): string {
  if (objectName) {
    return `${structure.formula} - ${objectName}`;
  }
  return `${structure.formula} band structure`;
}

/**
 * Estructura completa formateada para mostrar al usuario
 */
export function generateNanoHUBFormat(structure: NanoHUBStructure, objectName?: string): {
  title: string;
  premadeStructure: string;
  atomicCoordinates: string;
  structureType: string;
  atomicStructure: string;
  cellVectors: string;
  latticeParameterA: number;
  ratioBA: number;
  ratioCA: number;
} {
  return {
    title: generateTitle(structure, objectName),
    premadeStructure: structure.name,
    atomicCoordinates: "Fractional",
    structureType: structure.structureType,
    atomicStructure: generateAtomicStructure(structure),
    cellVectors: generateCellVectors(structure),
    latticeParameterA: structure.latticeParameter,
    ratioBA: structure.ratioBA,
    ratioCA: structure.ratioCA
  };
}

/**
 * Lista de estructuras disponibles
 */
export function getNanoHUBStructureList(): Array<{ id: string; name: string; formula: string; nat: number }> {
  return Object.values(NANOHUB_STRUCTURES).map(s => ({
    id: s.id,
    name: s.name,
    formula: s.formula,
    nat: s.atoms.length
  }));
}
