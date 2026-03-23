"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.env = void 0;
const dotenv_1 = require("dotenv");
const zod_1 = require("zod");
const path_1 = __importDefault(require("path"));
// Cargar variables de entorno desde .env
(0, dotenv_1.config)();
// Definir el esquema de validación para las variables de entorno
const envSchema = zod_1.z.object({
    TELEGRAM_BOT_TOKEN: zod_1.z.string().min(1, "El token del bot de Telegram es obligatorio"),
    TELEGRAM_ALLOWED_USER_IDS: zod_1.z.string().min(1, "Los IDs de usuario permitidos son obligatorios"),
    GROQ_API_KEY: zod_1.z.string().min(1, "La API KEY de Groq es obligatoria"),
    OPENROUTER_API_KEY: zod_1.z.string().optional(),
    OPENROUTER_MODEL: zod_1.z.string().default("openrouter/free"),
    DB_PATH: zod_1.z.string().default("./memory.db"), // Manteniendo por si hace falta, aunque no se usará
    GOOGLE_APPLICATION_CREDENTIALS: zod_1.z.string().min(1, "La ruta de credenciales de Google es obligatoria para Firebase"),
});
// Validar y parsear las variables de entorno
const parsedEnv = envSchema.safeParse(process.env);
if (!parsedEnv.success) {
    console.error("❌ Error en la configuración de variables de entorno:");
    console.error(parsedEnv.error.format());
    process.exit(1);
}
// Exportar configuración ya parseada y procesada
exports.env = {
    ...parsedEnv.data,
    // Convertir string de IDs separados por coma a un array de números
    allowedUserIds: parsedEnv.data.TELEGRAM_ALLOWED_USER_IDS
        .split(',')
        .map(id => parseInt(id.trim(), 10))
        .filter(id => !isNaN(id)),
    // Resolver ruta absoluta a la base de datos
    dbPath: path_1.default.resolve(process.cwd(), parsedEnv.data.DB_PATH)
};
