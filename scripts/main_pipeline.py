import os
import sys
import logging
import argparse

# Configuración de rutas
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)

from src.procesado_sec import procesar_genes
from src.visualizer import generar_graficos
from src.model_ia import entrenar_modelo

def inicializar_logs():
    log_path = os.path.join(BASE_DIR, "logs", "pipeline.log")
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(
        filename=log_path, level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def mostrar_tabla_comparativa(df):
    """Muestra un resumen estético en la consola al finalizar."""
    print("\n" + "="*90)
    print(" 📊 RESUMEN COMPARATIVO: IDENTIDAD VS. MÉTRICAS PROPIAS")
    print("="*90)
    
    # Seleccionamos las columnas clave para la comparativa
    # Asegúrate de que estos nombres coincidan exactamente con el diccionario en procesado_sec.py
    vista_previa = df[['gene_id', 'id_homologo', 'longitud_adn', 'contenido_gc']].head(10)
    
    # Renombramos para una visualización limpia en consola
    vista_previa.columns = ['ID Local', 'Hit NCBI (Ref)', 'Long. (bp)', '% GC Local']
    
    print(vista_previa.to_string(index=False))
    print("="*90 + "\n")

def main():
    # --- CONFIGURACIÓN DE ARGUMENTOS ---
    parser = argparse.ArgumentParser(description="BioGen Predictive Pipeline")
    parser.add_argument(
        "--online", 
        action="store_true", 
        help="Activa la búsqueda de homología en bases de datos del NCBI (BLAST) y descarga de referencias"
    )
    args = parser.parse_args()

    inicializar_logs()
    print("\n" + "🧬" + " ="*20)
    print(" INICIANDO BIOGEN PREDICTIVE PIPELINE")
    print(" ="*20 + "\n")
    
    if args.online:
        print("🌐 MODO ONLINE: BLAST activado. Se identificarán homólogos y descargarán referencias.")
    else:
        print("🔌 MODO OFFLINE: Análisis local. Se omitirá la validación externa contra NCBI.")

    # Definición de rutas de trabajo
    RAW_DIR = os.path.join(BASE_DIR, "data", "raw")
    PROCESSED_FILE = os.path.join(BASE_DIR, "data", "processed", "features_genes.csv")
    GRAPHICS_DIR = os.path.join(BASE_DIR, "graficos")
    RESULTS_DIR = os.path.join(BASE_DIR, "resultados")

    try:
        # FASE 1: Procesamiento, Homología y Descarga de Referencias
        # Ahora procesar_genes se encarga de llamar a fetch_tool internamente si buscar_online es True
        df = procesar_genes(RAW_DIR, PROCESSED_FILE, buscar_online=args.online)

        if df is not None and not df.empty:
            # FASE 2: Visualización y Correlación
            print("[*] Generando análisis visual y matriz de correlación...")
            corr_matriz = generar_graficos(df, GRAPHICS_DIR)
            
            os.makedirs(RESULTS_DIR, exist_ok=True)
            corr_matriz.to_csv(os.path.join(RESULTS_DIR, "matriz_correlacion.csv"))

            # FASE 3: Inteligencia Artificial y Control de Calidad
            # El modelo ahora devuelve r2, la predicción y la lista de anomalías detectadas
            print("[*] Entrenando modelo de IA y ejecutando detector de anomalías...")
            r2, pred, lista_anomalias = entrenar_modelo(df, valor_test=5000)
            
            print(f"\n" + "-"*40)
            print(f"✅ RESULTADOS DEL PIPELINE")
            print(f"-"*40)
            print(f"📊 Precisión del modelo estadístico (R²): {r2:.4f}")
            print(f"🔮 Predicción teórica para 5000bp: {pred:.2f}% GC")

            # Notificación de anomalías detectadas por el modelo
            if lista_anomalias:
                print(f"\n🚨 ALERTA DE CALIDAD (IA):")
                print(f"Se han detectado desviaciones significativas en: {', '.join(lista_anomalias)}")
                print("Revisar posibles errores de secuenciación o mutaciones raras en estos IDs.")

            # FASE 4: Presentación de datos
            mostrar_tabla_comparativa(df)
            
            if args.online:
                REF_DIR = os.path.join(BASE_DIR, "data", "reference")
                refs_descargadas = len([
                    f for f in os.listdir(REF_DIR)
                    if f.endswith((".fasta", ".fa"))
                ]) if os.path.isdir(REF_DIR) else 0
                print(f"📥 Referencias descargadas en 'data/reference/': {refs_descargadas} archivo(s).")

            print(f"📁 Proceso finalizado. Resultados guardados en: {RESULTS_DIR}")
                
        else:
            print("⚠️ No se encontraron archivos válidos en 'data/raw' para analizar.")

    except Exception as e:
        logging.error(f"Error crítico en el flujo principal: {e}", exc_info=True)
        print(f"\n❌ ERROR CRÍTICO: {e}")
        print(f"Consulte 'logs/pipeline.log' para obtener más detalles técnicos.")

if __name__ == "__main__":
    main()
