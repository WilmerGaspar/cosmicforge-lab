"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.chatCompletion = chatCompletion;
const groq_sdk_1 = require("groq-sdk");
const config_1 = require("../config");
const tools_1 = require("../tools");
const groq = new groq_sdk_1.Groq({
    apiKey: config_1.env.GROQ_API_KEY,
});
async function chatCompletion(messages) {
    try {
        const response = await groq.chat.completions.create({
            model: "llama-3.3-70b-versatile",
            messages: messages,
            tools: (0, tools_1.getToolsDefinitions)(),
            tool_choice: "auto",
            temperature: 0.2, // Baja temperatura para mantenerlo conciso y enfocado
            max_tokens: 1500,
        });
        return response.choices[0].message;
    }
    catch (error) {
        console.error("❌ Error en Groq LLM:", error);
        throw error;
    }
}
