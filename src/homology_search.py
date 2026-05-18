import os
import logging
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
from dotenv import load_dotenv

load_dotenv()
Entrez.email = os.getenv("ENTREZ_EMAIL", "")

def buscar_homologia_id(secuencia_adn):
    """
    Realiza BLAST, obtiene el ID y luego recupera la proteína de referencia.
    """
    try:
        logging.info("📡 Conectando con NCBI BLAST...")
        result_handle = NCBIWWW.qblast("blastn", "nt", secuencia_adn)
        blast_record = next(NCBIXML.parse(result_handle))
        
        if blast_record.alignments:
            mejor_hit = blast_record.alignments[0]
            accession_id = mejor_hit.accession
            titulo = mejor_hit.title
            
            logging.info(f"✅ Hit encontrado: {accession_id}. Recuperando proteína de referencia...")
            
            # --- LLAMADA INTERNA A LA SEGUNDA FUNCIÓN ---
            detalles = obtener_detalles_homologo(accession_id)
            
            return {
                "id": accession_id,
                "nombre": titulo,
                "prot_referencia": detalles["prot_ref"] if detalles else "No disponible",
                "producto_nombre": detalles["producto"] if detalles else "Desconocido"
            }
        
        logging.warning("⚠️ Sin homologías significativas.")
        return {"id": "N/A", "nombre": "Desconocido", "prot_referencia": "N/A", "producto_nombre": "N/A"}

    except Exception as e:
        logging.error(f"❌ Error en BLAST/Entrez: {e}")
        return {"id": "Error", "nombre": str(e), "prot_referencia": "Error", "producto_nombre": "Error"}

def obtener_detalles_homologo(accession_id):
    """
    Usa el ID para bajar el registro GenBank y extraer la proteína oficial.
    """
    try:
        # efetch busca el registro completo (GenBank es el formato más rico en información)
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        
        # El CDS (CoDing Sequence) contiene la traducción validada por curadores
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                return {
                    "prot_ref": feature.qualifiers["translation"][0],
                    "producto": feature.qualifiers.get("product", ["Desconocido"])[0]
                }
        return None
    except Exception as e:
        logging.warning(f"⚠️ No se pudieron obtener detalles para {accession_id}: {e}")
        return None