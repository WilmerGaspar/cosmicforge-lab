import { config } from 'dotenv';
import { z } from 'zod';
import path from 'path';

// Cargar variables de entorno desde .env
config();

// Definir el esquema de validación para las variables de entorno
const envSchema = z.object({
  TELEGRAM_BOT_TOKEN: z.string().min(1, "El token del bot de Telegram es obligatorio"),
  TELEGRAM_ALLOWED_USER_IDS: z.string().min(1, "Los IDs de usuario permitidos son obligatorios"),
  GROQ_API_KEY: z.string().min(1, "La API KEY de Groq es obligatoria"),
  OPENROUTER_API_KEY: z.string().optional(),
  OPENROUTER_MODEL: z.string().default("openrouter/free"),
  DB_PATH: z.string().default("./memory.db"), // Manteniendo por si hace falta, aunque no se usará
  GOOGLE_APPLICATION_CREDENTIALS: z.string().min(1, "La ruta de credenciales de Google es obligatoria para Firebase"),
});

// Validar y parsear las variables de entorno
const parsedEnv = envSchema.safeParse(process.env);

if (!parsedEnv.success) {
  console.error("❌ Error en la configuración de variables de entorno:");
  console.error(parsedEnv.error.format());
  process.exit(1);
}

// Exportar configuración ya parseada y procesada
export const env = {
  ...parsedEnv.data,
  // Convertir string de IDs separados por coma a un array de números
  allowedUserIds: parsedEnv.data.TELEGRAM_ALLOWED_USER_IDS
    .split(',')
    .map(id => parseInt(id.trim(), 10))
    .filter(id => !isNaN(id)),
  // Resolver ruta absoluta a la base de datos
  dbPath: path.resolve(process.cwd(), parsedEnv.data.DB_PATH)
};
