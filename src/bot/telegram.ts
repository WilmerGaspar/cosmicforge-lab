import { Bot } from 'grammy';
import { env } from '../config';
import { processUserMessage } from '../agent/loop';

export const bot = new Bot(env.TELEGRAM_BOT_TOKEN);

// Middleware de seguridad: Whitelist
bot.use(async (ctx, next) => {
  if (ctx.from && ctx.from.id) {
    if (env.allowedUserIds.includes(ctx.from.id)) {
      return next(); // Usuario autorizado
    }
  }
  // Alternativa: ignorar silenciosamente en lugar de responder para evitar descubrir el bot
  console.log(`[Seguridad] Intento de acceso denegado del ID: ${ctx.from?.id}`);
});

bot.command('start', async (ctx) => {
  await ctx.reply(`Hola 👋. Soy **OpenGravity**, tu asistente personal de IA seguro y local.\n\nPuedes empezar a hablarme de lo que necesites.`, { parse_mode: 'Markdown' });
});

bot.on('message:text', async (ctx) => {
  const userId = ctx.from.id;
  const text = ctx.message.text;

  // Indicar acción de escritura en Telegram
  await ctx.replyWithChatAction('typing');

  try {
    const replyText = await processUserMessage(userId, text);
    
    // El LLM podría generar una respuesta muy larga para un solo mensaje de Telegram (límite 4096 caracteres)
    if (replyText.length > 4000) {
      // Simplificación: dividir mensaje o enviar primeros 4000
      await ctx.reply(replyText.slice(0, 4000) + "\n\n...[Mensaje truncado]...", { parse_mode: 'Markdown' });
    } else {
      await ctx.reply(replyText, { parse_mode: 'Markdown' });
    }
    
  } catch (error) {
    console.error('[Bot] Error enviando respuesta:', error);
    await ctx.reply('❌ Ha ocurrido un error interno comunicándose con el motor LLM.', { parse_mode: 'Markdown' });
  }
});
