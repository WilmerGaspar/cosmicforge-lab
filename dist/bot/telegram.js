"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.bot = void 0;
const grammy_1 = require("grammy");
const config_1 = require("../config");
const loop_1 = require("../agent/loop");
exports.bot = new grammy_1.Bot(config_1.env.TELEGRAM_BOT_TOKEN);
// Middleware de seguridad: Whitelist
exports.bot.use(async (ctx, next) => {
    if (ctx.from && ctx.from.id) {
        if (config_1.env.allowedUserIds.includes(ctx.from.id)) {
            return next(); // Usuario autorizado
        }
    }
    // Alternativa: ignorar silenciosamente en lugar de responder para evitar descubrir el bot
    console.log(`[Seguridad] Intento de acceso denegado del ID: ${ctx.from?.id}`);
});
exports.bot.command('start', async (ctx) => {
    await ctx.reply(`Hola 👋. Soy **OpenGravity**, tu asistente personal de IA seguro y local.\n\nPuedes empezar a hablarme de lo que necesites.`, { parse_mode: 'Markdown' });
});
exports.bot.on('message:text', async (ctx) => {
    const userId = ctx.from.id;
    const text = ctx.message.text;
    // Indicar acción de escritura en Telegram
    await ctx.replyWithChatAction('typing');
    try {
        const replyText = await (0, loop_1.processUserMessage)(userId, text);
        // El LLM podría generar una respuesta muy larga para un solo mensaje de Telegram (límite 4096 caracteres)
        if (replyText.length > 4000) {
            // Simplificación: dividir mensaje o enviar primeros 4000
            await ctx.reply(replyText.slice(0, 4000) + "\n\n...[Mensaje truncado]...", { parse_mode: 'Markdown' });
        }
        else {
            await ctx.reply(replyText, { parse_mode: 'Markdown' });
        }
    }
    catch (error) {
        console.error('[Bot] Error enviando respuesta:', error);
        await ctx.reply('❌ Ha ocurrido un error interno comunicándose con el motor LLM.', { parse_mode: 'Markdown' });
    }
});
