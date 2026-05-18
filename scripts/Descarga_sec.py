import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.fetch_tools import descargar_secuencia_homologa

# Ruta absoluta a la raíz del proyecto
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    mis_genes = ["NM_000546", "NM_007294", "NM_000041", "NM_001301717"]

    # Ruta absoluta: siempre apunta a BioGen-Predictive-Pipeline/data/raw/
    ruta_datos = os.path.join(BASE_DIR, "data", "raw")

    print("--- Iniciando descarga masiva ---")
    exitosos, fallidos = 0, []

    for gen in mis_genes:
        try:
            resultado = descargar_secuencia_homologa(gen, ruta_datos)
            if resultado:
                print(f"  ✅ {gen} descargado correctamente.")
                exitosos += 1
            else:
                print(f"  ⚠️  {gen} no se descargó (sin resultado).")
                fallidos.append(gen)
        except Exception as e:
            print(f"  ❌ {gen} falló: {e}. Continuando con el siguiente...")
            fallidos.append(gen)

    print(f"\n--- Proceso finalizado: {exitosos}/{len(mis_genes)} exitosos ---")
    if fallidos:
        print(f"    Genes con error: {', '.join(fallidos)}")