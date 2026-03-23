"use strict";
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    var desc = Object.getOwnPropertyDescriptor(m, k);
    if (!desc || ("get" in desc ? !m.__esModule : desc.writable || desc.configurable)) {
      desc = { enumerable: true, get: function() { return m[k]; } };
    }
    Object.defineProperty(o, k2, desc);
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
    Object.defineProperty(o, "default", { enumerable: true, value: v });
}) : function(o, v) {
    o["default"] = v;
});
var __importStar = (this && this.__importStar) || (function () {
    var ownKeys = function(o) {
        ownKeys = Object.getOwnPropertyNames || function (o) {
            var ar = [];
            for (var k in o) if (Object.prototype.hasOwnProperty.call(o, k)) ar[ar.length] = k;
            return ar;
        };
        return ownKeys(o);
    };
    return function (mod) {
        if (mod && mod.__esModule) return mod;
        var result = {};
        if (mod != null) for (var k = ownKeys(mod), i = 0; i < k.length; i++) if (k[i] !== "default") __createBinding(result, mod, k[i]);
        __setModuleDefault(result, mod);
        return result;
    };
})();
Object.defineProperty(exports, "__esModule", { value: true });
exports.addMessage = addMessage;
exports.getHistory = getHistory;
exports.clearHistory = clearHistory;
const admin = __importStar(require("firebase-admin"));
// Inicializar la app de Firebase Admin
// Requiere que la variable de entorno GOOGLE_APPLICATION_CREDENTIALS apunte al service-account.json
if (!admin.apps.length) {
    admin.initializeApp();
}
const db = admin.firestore();
/**
 * Agrega un mensaje al historial de un usuario en Firestore
 */
async function addMessage(userId, role, content = null, toolCalls = null, toolCallId = null) {
    try {
        const userDocRef = db.collection('users').doc(userId.toString());
        const messagesCollection = userDocRef.collection('messages');
        await messagesCollection.add({
            role,
            content,
            tool_calls: toolCalls,
            tool_call_id: toolCallId,
            created_at: admin.firestore.FieldValue.serverTimestamp()
        });
    }
    catch (error) {
        console.error("❌ Error guardando mensaje en Firestore:", error);
    }
}
/**
 * Obtiene el historial de mensajes de un usuario ordenados por fecha de creación
 */
async function getHistory(userId, limit = 30) {
    try {
        const userDocRef = db.collection('users').doc(userId.toString());
        const messagesSnapshot = await userDocRef
            .collection('messages')
            .orderBy('created_at', 'desc')
            .limit(limit)
            .get();
        // Invertir para que estén ordenados cronológicamente
        const docs = messagesSnapshot.docs.reverse();
        return docs.map(doc => {
            const data = doc.data();
            return {
                role: data.role,
                content: data.content || null,
                tool_calls: data.tool_calls || null,
                tool_call_id: data.tool_call_id || null,
            };
        });
    }
    catch (error) {
        console.error("❌ Error obteniendo historial de Firestore:", error);
        return [];
    }
}
/**
 * Borra el historial de un usuario
 */
async function clearHistory(userId) {
    try {
        const userDocRef = db.collection('users').doc(userId.toString());
        const messagesCollection = userDocRef.collection('messages');
        // Para borrar subcolecciones en Firebase, necesitamos buscar los docs y eliminarlos en lote
        const snapshot = await messagesCollection.get();
        if (snapshot.size === 0)
            return;
        const batch = db.batch();
        snapshot.docs.forEach((doc) => {
            batch.delete(doc.ref);
        });
        await batch.commit();
    }
    catch (error) {
        console.error("❌ Error limpiando historial de Firestore:", error);
    }
}
