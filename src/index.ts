import { bot } from './bot/telegram';

async function bootstrap() {
  console.log("🚀 Iniciando OpenGravity...");
  
  // Limpiar/Cerrar la base de datos de SQLite o el bot de forma gracefully si es necesario
  process.once('SIGINT', () => {
    bot.stop();
    console.log("Bot detenido por SIGINT");
    process.exit(0);
  });
  process.once('SIGTERM', () => {
    bot.stop();
    console.log("Bot detenido por SIGTERM");
    process.exit(0);
  });

  bot.catch((err) => {
    console.error("❌ Error en el long polling del bot de Telegram:", err);
  });

  console.log("🤖 OpenGravity está escuchando en Telegram...");
  await bot.start();
}

bootstrap().catch(error => {
  console.error("❌ Error fatal al iniciar OpenGravity:", error);
  process.exit(1);
});
