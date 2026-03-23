"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.processUserMessage = processUserMessage;
const client_1 = require("../llm/client");
const tools_1 = require("../tools");
const database_1 = require("../db/database");
const MAX_ITERATIONS = 5;
// Instrucción del sistema
const SYSTEM_PROMPT = `Eres OpenGravity, un asistente personal de inteligencia artificial privado que opera a través de Telegram y se ejecuta de manera local. 
Eres servicial, conciso, y seguro. Sólo hablas español.
Puedes acceder a herramientas para ayudar al usuario de forma precisa. 
Responde de manera natural como un asistente en un chat. Usa formato markdown para destacar información (negritas, código, listas, etc.).`;
async function processUserMessage(userId, text) {
    // 1. Guardar mensaje del usuario en persistencia
    await (0, database_1.addMessage)(userId, 'user', text);
    let iterations = 0;
    while (iterations < MAX_ITERATIONS) {
        iterations++;
        // 2. Obtener historial actualizado
        const history = await (0, database_1.getHistory)(userId, 30);
        // Mapear historial al formato del SDK de Groq / OpenAI
        const messages = [
            { role: 'system', content: SYSTEM_PROMPT },
            ...history.map((msg) => {
                const payload = { role: msg.role };
                if (msg.content)
                    payload.content = msg.content;
                if (msg.role === 'assistant' && msg.tool_calls) {
                    payload.tool_calls = JSON.parse(msg.tool_calls);
                }
                else if (msg.role === 'tool' && msg.tool_call_id) {
                    payload.tool_call_id = msg.tool_call_id;
                }
                return payload;
            })
        ];
        try {
            // 3. Obtener respuesta del LLM
            const response = await (0, client_1.chatCompletion)(messages);
            // 4. Evaluar si hay llamadas a herramientas
            if (response.tool_calls && response.tool_calls.length > 0) {
                // Guardar la intención de llamar a herramientas
                await (0, database_1.addMessage)(userId, 'assistant', response.content, JSON.stringify(response.tool_calls));
                // Por cada herramienta llamada, ejecutarla
                for (const toolCall of response.tool_calls) {
                    const fnName = toolCall.function.name;
                    const fnArgs = toolCall.function.arguments;
                    console.log(`[Agent] Ejecutando herramienta: ${fnName}(${fnArgs})`);
                    const resultStr = await (0, tools_1.executeTool)(fnName, fnArgs);
                    // Guardar resultado de herramienta
                    await (0, database_1.addMessage)(userId, 'tool', resultStr, null, toolCall.id);
                }
                // Continuar el bucle para que el LLM genere una respuesta con la nueva información
                continue;
            }
            // 5. Si no hay herramientas, es un mensaje final del asistente
            if (response.content) {
                await (0, database_1.addMessage)(userId, 'assistant', response.content);
                return response.content;
            }
            // Si por alguna razón la respuesta está vacía pero no es tool call
            return "Hubo un error al procesar tu respuesta (respuesta vacía del LLM).";
        }
        catch (error) {
            console.error("[Agent Loop] Error:", error);
            return `❌ Ocurrió un error al procesar tu solicitud: ${error.message || 'Error desconocido'}`;
        }
    }
    // 6. Si supera las iteraciones (Loop infinito prevención)
    return "Se alcanzó el límite máximo de operaciones internas. Intenta reformular tu pregunta.";
}
