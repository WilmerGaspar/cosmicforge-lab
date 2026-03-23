import { env } from '../config';
import * as admin from 'firebase-admin';

// Inicializar la app de Firebase Admin
// Requiere que la variable de entorno GOOGLE_APPLICATION_CREDENTIALS apunte al service-account.json
if (!admin.apps.length) {
  admin.initializeApp();
}

const db = admin.firestore();

export interface DBMessage {
  role: 'system' | 'user' | 'assistant' | 'tool';
  content: string | null;
  tool_calls: string | null;
  tool_call_id: string | null;
  created_at?: admin.firestore.Timestamp;
}

/**
 * Agrega un mensaje al historial de un usuario en Firestore
 */
export async function addMessage(
  userId: number, 
  role: 'system' | 'user' | 'assistant' | 'tool', 
  content: string | null = null, 
  toolCalls: string | null = null, 
  toolCallId: string | null = null
) {
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
  } catch (error) {
    console.error("❌ Error guardando mensaje en Firestore:", error);
  }
}

/**
 * Obtiene el historial de mensajes de un usuario ordenados por fecha de creación
 */
export async function getHistory(userId: number, limit: number = 30): Promise<DBMessage[]> {
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
        role: data.role as DBMessage['role'],
        content: data.content || null,
        tool_calls: data.tool_calls || null,
        tool_call_id: data.tool_call_id || null,
      };
    });
  } catch (error) {
    console.error("❌ Error obteniendo historial de Firestore:", error);
    return [];
  }
}

/**
 * Borra el historial de un usuario
 */
export async function clearHistory(userId: number) {
  try {
    const userDocRef = db.collection('users').doc(userId.toString());
    const messagesCollection = userDocRef.collection('messages');
    
    // Para borrar subcolecciones en Firebase, necesitamos buscar los docs y eliminarlos en lote
    const snapshot = await messagesCollection.get();
    
    if (snapshot.size === 0) return;

    const batch = db.batch();
    snapshot.docs.forEach((doc) => {
      batch.delete(doc.ref);
    });

    await batch.commit();
  } catch (error) {
    console.error("❌ Error limpiando historial de Firestore:", error);
  }
}
