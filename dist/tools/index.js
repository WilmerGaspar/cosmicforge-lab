"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.toolsRegistry = void 0;
exports.getToolsDefinitions = getToolsDefinitions;
exports.executeTool = executeTool;
const time_1 = require("./time");
// Registro de todas las herramientas disponibles
exports.toolsRegistry = {
    [time_1.timeTool.definition.function.name]: time_1.timeTool,
};
// Devuelve las definiciones de las herramientas en el formato que espera Groq / OpenAI
function getToolsDefinitions() {
    return Object.values(exports.toolsRegistry).map(tool => tool.definition);
}
// Ejecuta una herramienta dado su nombre y argumentos
async function executeTool(name, argsStr) {
    const tool = exports.toolsRegistry[name];
    if (!tool) {
        return JSON.stringify({ error: `La herramienta ${name} no existe.` });
    }
    try {
        const args = JSON.parse(argsStr);
        const result = await tool.execute(args);
        return typeof result === 'string' ? result : JSON.stringify(result);
    }
    catch (error) {
        console.error(`❌ Error ejecutando herramienta ${name}:`, error);
        return JSON.stringify({ error: `Error al ejecutar ${name}: ${error.message}` });
    }
}
