"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.timeTool = void 0;
exports.timeTool = {
    definition: {
        type: "function",
        function: {
            name: "get_current_time",
            description: "Obtiene la fecha y hora actual del sistema donde corre el agente.",
            parameters: {
                type: "object",
                properties: {},
                required: [],
            },
        },
    },
    execute: () => {
        const now = new Date();
        return {
            success: true,
            time: now.toISOString(),
            local_time: now.toLocaleString('es-ES'),
        };
    },
};
