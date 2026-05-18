import os
import logging
from Bio import Entrez, SeqIO
from dotenv import load_dotenv

load_dotenv()

# Ruta para los Logs
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGS_DIR = os.path.join(BASE_DIR, "logs")

# Aseguramos que la carpeta de logs exista
os.makedirs(LOGS_DIR, exist_ok=True)

# Configuración del sistema de logs (evitando duplicados)
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig(
    filename=os.path.join(LOGS_DIR, 'pipeline.log'),
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def descargar_secuencia_homologa(gene_id, carpeta_destino):
    """
    Descarga una única secuencia de nucleótidos del NCBI basada en su ID.
    Útil para bajar la referencia encontrada en el proceso de homología.
    """
    # Aseguramos que la carpeta de destino exista (ej: data/reference)
    os.makedirs(carpeta_destino, exist_ok=True)
    filename = os.path.join(carpeta_destino, f"REF_{gene_id}.fasta")

    # Si el archivo ya existe, no lo descargamos de nuevo (ahorro de tiempo y red)
    if os.path.exists(filename):
        logging.info(f"El archivo de referencia {gene_id} ya existe localmente.")
        return filename

    # Solo necesitamos el email si vamos a hacer una llamada a NCBI
    Entrez.email = os.getenv("ENTREZ_EMAIL", "")
    if not Entrez.email:
        raise EnvironmentError("Define ENTREZ_EMAIL en el archivo .env antes de ejecutar el pipeline.")

    try:
        logging.info(f"Iniciando descarga de referencia NCBI: {gene_id}...")

        # Conexión con la BBDD nucleotide del NCBI
        with Entrez.efetch(db="nucleotide",
                           id=gene_id, 
                           rettype="fasta",
                           retmode="text") as handle:
            record = SeqIO.read(handle, "fasta")

        # Guardado en formato fasta
        SeqIO.write(record, filename, "fasta")

        file_size_kb = os.path.getsize(filename) / 1024
        
        if file_size_kb < 0.1:
            logging.warning(f"⚠️ {gene_id} es muy pequeño ({file_size_kb:.3f} KB)")
        else:
            logging.info(f"✅ Referencia {gene_id} guardada con éxito ({file_size_kb:.2f} KB)")
        
        return filename

    except Exception as e:
        logging.error(f"❌ Error al descargar referencia {gene_id}: {str(e)}")
        return None