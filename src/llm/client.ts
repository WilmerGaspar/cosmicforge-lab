import { Groq } from 'groq-sdk';
import { env } from '../config';
import { getToolsDefinitions } from '../tools';

const groq = new Groq({
  apiKey: env.GROQ_API_KEY,
});

export async function chatCompletion(messages: any[]) {
  try {
    const response = await groq.chat.completions.create({
      model: "llama-3.3-70b-versatile",
      messages: messages,
      tools: getToolsDefinitions(),
      tool_choice: "auto",
      temperature: 0.2, // Baja temperatura para mantenerlo conciso y enfocado
      max_tokens: 1500,
    });

    return response.choices[0].message;
  } catch (error) {
    console.error("❌ Error en Groq LLM:", error);
    throw error;
  }
}
