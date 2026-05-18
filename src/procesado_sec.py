import os
import pandas as pd
import logging
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Importaciones de tus módulos personalizados
from src.homology_search import buscar_homologia_id
from src.fetch_tools import descargar_secuencia_homologa

# Configuración de Logs
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGS_DIR = os.path.join(BASE_DIR, "logs")
os.makedirs(LOGS_DIR, exist_ok=True)

for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig(
    filename=os.path.join(LOGS_DIR, 'pipeline.log'),
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def procesar_genes(ruta_entrada, ruta_salida, buscar_online=False):
    """
    Orquestador avanzado:
    1. Identifica homología.
    2. Descarga secuencia de referencia oficial (si online).
    3. Calcula métricas locales para comparación.
    """
    datos_procesados = []
    
    # Definimos dónde guardaremos las referencias descargadas
    RUTA_REFERENCIAS = os.path.join(BASE_DIR, "data", "reference")

    if not os.path.exists(ruta_entrada):
        logging.error(f"La ruta de entrada no existe: {ruta_entrada}")
        return None

    archivos = [f for f in os.listdir(ruta_entrada) if f.endswith((".fasta", ".fa"))]
    logging.info(f"🚀 Iniciando pipeline (Modo Online: {buscar_online}).")

    for archivo in archivos:
        path = os.path.join(ruta_entrada, archivo)
        try:
            for record in SeqIO.parse(path, "fasta"):
                adn = record.seq
                
                # --- PASO 1: Identificación y Descarga de Referencia ---
                id_homologo = "N/A"
                nombre_homologo = "Modo Offline"
                prot_ref_seq = "N/A"

                if buscar_online:
                    hit = buscar_homologia_id(str(adn))
                    id_homologo = hit['id']
                    nombre_homologo = hit['nombre']
                    prot_ref_seq = hit.get('prot_referencia', "No disponible")

                    # Si hay un hit válido, descargamos el FASTA de referencia usando fetch_tool
                    if id_homologo not in ["N/A", "Error"]:
                        descargar_secuencia_homologa(id_homologo, RUTA_REFERENCIAS)

                # --- PASO 2: Métricas de ADN locales ---
                gc = round(gc_fraction(adn) * 100, 2)
                
                # --- PASO 3: Traducción y métricas locales ---
                info_prot = realizar_traduccion(adn, record.id)
                
                # Registro final unificado para la tabla comparativa
                registro = {
                    "gene_id": record.id,
                    "fuente_homologia": nombre_homologo,
                    "id_homologo": id_homologo,
                    "longitud_adn": len(adn),
                    "contenido_gc": gc,
                    "peso_proteina_da": info_prot["peso_molecular_da"],
                    "secuencia_proteica_local": info_prot["secuencia_proteica"],
                    "secuencia_proteica_ref": prot_ref_seq
                }
                datos_procesados.append(registro)
                
            logging.info(f"✅ Procesado: {archivo}")

        except Exception as e:
            logging.error(f"❌ Error en {archivo}: {str(e)}")
            continue

    if datos_procesados:
        df = pd.DataFrame(datos_procesados)
        os.makedirs(os.path.dirname(ruta_salida), exist_ok=True)
        df.to_csv(ruta_salida, index=False)
        logging.info(f"💾 Dataset guardado en {ruta_salida}")
        return df
    
    return None

def realizar_traduccion(secuencia_adn, gene_id):
    """
    Traduce ADN a Proteína y calcula parámetros biofísicos.
    """
    try:
        logging.info(f"🧬 Traduciendo localmente: {gene_id}")
        proteina_seq = secuencia_adn.translate(to_stop=True)
        sec_str = str(proteina_seq)
        
        if len(sec_str) > 0:
            analisis = ProteinAnalysis(sec_str)
            peso_mol = round(analisis.molecular_weight(), 2)
            aromaticidad = round(analisis.aromaticity(), 3)
        else:
            peso_mol, aromaticidad = 0, 0
            logging.warning(f"⚠️ {gene_id} produjo secuencia vacía.")

        return {
            "longitud_proteina": len(sec_str),
            "peso_molecular_da": peso_mol,
            "aromaticidad": aromaticidad,
            "secuencia_proteica": sec_str
        }

    except Exception as e:
        logging.error(f"❗ Error en traducción de {gene_id}: {e}")
        return {
            "longitud_proteina": 0, 
            "peso_molecular_da": 0,
            "aromaticidad": 0, 
            "secuencia_proteica": None
        }