import os
import logging
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

def entrenar_modelo(df: pd.DataFrame, valor_test: int):
    """
    Entrena un modelo para analizar la coherencia entre 
    la secuencia local y la referencia encontrada.
    """
    try:
        # 1. Limpieza de datos: Solo filas que tengan homología (no N/A)
        df_valid = df[df['id_homologo'] != 'N/A'].copy()

        if df_valid.empty:
            logging.warning("⚠️ No hay suficientes datos con homología para entrenar el modelo.")
            return 0, 0, []

        # 2. Nueva variable: Diferencia de longitud (Local vs Referencia si la tuviéramos)
        # Por ahora, mantengamos la lógica de predicción de GC pero aplicada al dataset validado
        X = df_valid[['longitud_adn']].values
        y = df_valid['contenido_gc'].values

        # 3. Entrenar el modelo "Estándar"
        modelo = LinearRegression()
        modelo.fit(X, y)
        r2 = modelo.score(X, y)

        # 4. Predicción para un valor dado
        prediccion = modelo.predict([[valor_test]])[0]

        # --- LOGICA DE BIO-VALIDACIÓN ---
        # Calculamos la desviación de cada gen respecto a lo "normal" en el dataset
        df_valid['prediccion_gc'] = modelo.predict(X)
        df_valid['desviacion'] = np.abs(df_valid['contenido_gc'] - df_valid['prediccion_gc'])
        
        # Si un gen se desvía más de 2 veces la desviación estándar, es una anomalía
        umbral = df_valid['desviacion'].std() * 2
        anomalias = df_valid[df_valid['desviacion'] > umbral]
        lista_anomalias = anomalias['gene_id'].tolist() if 'gene_id' in anomalias.columns else []

        if lista_anomalias:
            logging.warning(f"🚨 Se detectaron {len(lista_anomalias)} secuencias con %GC anómalo respecto al patrón.")

        logging.info(f"Modelo IA entrenado: R2={r2:.4f}. Predicción para {valor_test}bp: {prediccion:.2f}% GC")
        
        return r2, prediccion, lista_anomalias

    except Exception as e:
        logging.error(f"❌ Error en el entrenamiento del modelo: {e}")
        return 0, 0, []