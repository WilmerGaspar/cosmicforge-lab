import { timeTool } from './time';

// Interfaz que deben cumplir todas las herramientas implementadas
export interface Tool {
  definition: {
    type: "function";
    function: {
      name: string;
      description: string;
      parameters: Record<string, any>;
    };
  };
  execute: (args: any) => Promise<any> | any;
}

// Registro de todas las herramientas disponibles
export const toolsRegistry: Record<string, Tool> = {
  [timeTool.definition.function.name]: timeTool,
};

// Devuelve las definiciones de las herramientas en el formato que espera Groq / OpenAI
export function getToolsDefinitions() {
  return Object.values(toolsRegistry).map(tool => tool.definition);
}

// Ejecuta una herramienta dado su nombre y argumentos
export async function executeTool(name: string, argsStr: string): Promise<string> {
  const tool = toolsRegistry[name];
  if (!tool) {
    return JSON.stringify({ error: `La herramienta ${name} no existe.` });
  }

  try {
    const args = JSON.parse(argsStr);
    const result = await tool.execute(args);
    return typeof result === 'string' ? result : JSON.stringify(result);
  } catch (error: any) {
    console.error(`❌ Error ejecutando herramienta ${name}:`, error);
    return JSON.stringify({ error: `Error al ejecutar ${name}: ${error.message}` });
  }
}
