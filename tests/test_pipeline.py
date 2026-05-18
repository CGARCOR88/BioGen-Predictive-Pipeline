"""
Tests unitarios para BioGen-Predictive-Pipeline.
Cubre: procesado_sec, model_ia, visualizer y fetch_tools.
Ejecutar desde la raíz del proyecto:
    pytest tests/ -v
"""

import os
import sys
import tempfile
import shutil
import unittest
from unittest.mock import patch, MagicMock

import matplotlib
matplotlib.use('Agg')  # backend no interactivo: evita errores de Tk en tests

import pandas as pd
import numpy as np
from Bio.Seq import Seq

# Añadir raíz del proyecto al path para que los imports de src/ funcionen
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE_DIR)

# ─────────────────────────────────────────────────────────────
# DATOS DE PRUEBA COMPARTIDOS
# ─────────────────────────────────────────────────────────────
FASTA_CON_ORF = ">gene_test\nATGAAAGCAATTTCGACCGAGTTGTAA\n"
FASTA_SIN_ORF = ">gene_norf\nTAAAAAAAAAAAAAAAAAAAAAAAAAA\n"  # empieza con codón stop TAA
FASTA_GC_ALTO = ">gene_gc\nATGCGCGCGCGCGCGCGCGCGCGCTAA\n"

# DataFrame mínimo con homología para tests de model_ia
def make_df_con_homologia(n=6):
    return pd.DataFrame({
        "gene_id":      [f"G{i}" for i in range(n)],
        "id_homologo":  ["NM_000001"] * n,
        "longitud_adn": [1000, 2000, 3000, 4000, 5000, 6000][:n],
        "contenido_gc": [52.0, 50.0, 48.0, 46.0, 44.0, 42.0][:n],
    })

def make_df_sin_homologia(n=4):
    return pd.DataFrame({
        "gene_id":      [f"G{i}" for i in range(n)],
        "id_homologo":  ["N/A"] * n,
        "longitud_adn": [1000, 2000, 3000, 4000][:n],
        "contenido_gc": [50.0, 48.0, 46.0, 44.0][:n],
    })


# ═════════════════════════════════════════════════════════════
# 1. TESTS: procesado_sec — realizar_traduccion
# ═════════════════════════════════════════════════════════════
class TestRealizarTraduccion(unittest.TestCase):
    """Tests para la función realizar_traduccion de procesado_sec.py."""

    def setUp(self):
        # Importamos aquí para evitar que el import-level lance errors si falta .env
        from src.procesado_sec import realizar_traduccion
        self.realizar_traduccion = realizar_traduccion

    def test_traduccion_normal_devuelve_dict(self):
        """Una secuencia con ORF válido devuelve un dict con las 4 claves."""
        seq = Seq("ATGAAAGCAATTTCGACCGAGTTGTAA")
        resultado = self.realizar_traduccion(seq, "test_gene")
        self.assertIsInstance(resultado, dict)
        self.assertIn("longitud_proteina", resultado)
        self.assertIn("peso_molecular_da", resultado)
        self.assertIn("aromaticidad", resultado)
        self.assertIn("secuencia_proteica", resultado)

    def test_traduccion_orf_produce_proteina(self):
        """Un ORF válido produce secuencia proteica no vacía."""
        seq = Seq("ATGAAAGCAATTTCGACCGAGTTGTAA")
        resultado = self.realizar_traduccion(seq, "test_gene")
        self.assertGreater(resultado["longitud_proteina"], 0)
        self.assertIsNotNone(resultado["secuencia_proteica"])
        self.assertNotEqual(resultado["secuencia_proteica"], "")

    def test_traduccion_sin_orf_devuelve_longitud_cero(self):
        """Una secuencia que empieza con codón stop no produce proteína."""
        seq = Seq("TAATTTTTTTTTTTTTTTTTTTTTTTT")  # TAA = stop codon al inicio
        resultado = self.realizar_traduccion(seq, "gene_norf")
        self.assertEqual(resultado["longitud_proteina"], 0)
        self.assertEqual(resultado["peso_molecular_da"], 0)

    def test_traduccion_peso_molecular_positivo(self):
        """El peso molecular de una proteína válida es positivo."""
        seq = Seq("ATGAAAGCAATTTCGACCGAGTTGTAA")
        resultado = self.realizar_traduccion(seq, "test_gene")
        self.assertGreater(resultado["peso_molecular_da"], 0)

    def test_traduccion_aromaticidad_en_rango(self):
        """La aromaticidad debe estar entre 0 y 1."""
        seq = Seq("ATGAAAGCAATTTCGACCGAGTTGTAA")
        resultado = self.realizar_traduccion(seq, "test_gene")
        self.assertGreaterEqual(resultado["aromaticidad"], 0)
        self.assertLessEqual(resultado["aromaticidad"], 1)


# ═════════════════════════════════════════════════════════════
# 2. TESTS: procesado_sec — procesar_genes (modo offline)
# ═════════════════════════════════════════════════════════════
class TestProcesarGenes(unittest.TestCase):
    """Tests para procesar_genes en modo offline (sin NCBI)."""

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.output_csv = os.path.join(self.tmp_dir, "features.csv")

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def _escribir_fasta(self, nombre, contenido):
        ruta = os.path.join(self.tmp_dir, nombre)
        with open(ruta, "w") as f:
            f.write(contenido)
        return ruta

    def test_ruta_invalida_retorna_none(self):
        """Si la ruta de entrada no existe, retorna None."""
        from src.procesado_sec import procesar_genes
        resultado = procesar_genes("/ruta/inexistente", self.output_csv, buscar_online=False)
        self.assertIsNone(resultado)

    def test_procesa_fasta_con_orf_offline(self):
        """Con un FASTA válido en modo offline devuelve un DataFrame."""
        self._escribir_fasta("gen1.fasta", FASTA_CON_ORF)
        from src.procesado_sec import procesar_genes
        df = procesar_genes(self.tmp_dir, self.output_csv, buscar_online=False)
        self.assertIsNotNone(df)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertFalse(df.empty)

    def test_columnas_requeridas_presentes(self):
        """El DataFrame resultante contiene las columnas esenciales."""
        self._escribir_fasta("gen1.fasta", FASTA_CON_ORF)
        from src.procesado_sec import procesar_genes
        df = procesar_genes(self.tmp_dir, self.output_csv, buscar_online=False)
        for col in ["gene_id", "longitud_adn", "contenido_gc", "id_homologo"]:
            self.assertIn(col, df.columns, f"Falta columna: {col}")

    def test_guarda_csv(self):
        """El CSV de salida se crea en la ruta indicada."""
        self._escribir_fasta("gen1.fasta", FASTA_CON_ORF)
        from src.procesado_sec import procesar_genes
        procesar_genes(self.tmp_dir, self.output_csv, buscar_online=False)
        self.assertTrue(os.path.exists(self.output_csv))

    def test_id_homologo_es_na_en_offline(self):
        """En modo offline, id_homologo debe ser 'N/A'."""
        self._escribir_fasta("gen1.fasta", FASTA_CON_ORF)
        from src.procesado_sec import procesar_genes
        df = procesar_genes(self.tmp_dir, self.output_csv, buscar_online=False)
        self.assertTrue((df["id_homologo"] == "N/A").all())

    def test_carpeta_sin_fastas_retorna_none(self):
        """Una carpeta sin archivos .fasta retorna None."""
        from src.procesado_sec import procesar_genes
        resultado = procesar_genes(self.tmp_dir, self.output_csv, buscar_online=False)
        self.assertIsNone(resultado)

    def test_contenido_gc_entre_0_y_100(self):
        """El contenido GC calculado debe estar entre 0 y 100."""
        self._escribir_fasta("gen1.fasta", FASTA_CON_ORF)
        self._escribir_fasta("gen2.fasta", FASTA_GC_ALTO)
        from src.procesado_sec import procesar_genes
        df = procesar_genes(self.tmp_dir, self.output_csv, buscar_online=False)
        self.assertTrue((df["contenido_gc"] >= 0).all())
        self.assertTrue((df["contenido_gc"] <= 100).all())

    def test_longitud_adn_positiva(self):
        """La longitud del ADN debe ser mayor que cero."""
        self._escribir_fasta("gen1.fasta", FASTA_CON_ORF)
        from src.procesado_sec import procesar_genes
        df = procesar_genes(self.tmp_dir, self.output_csv, buscar_online=False)
        self.assertTrue((df["longitud_adn"] > 0).all())

    def test_procesar_multiples_fastas(self):
        """Procesar dos FASTAs produce dos registros."""
        self._escribir_fasta("gen1.fasta", FASTA_CON_ORF)
        self._escribir_fasta("gen2.fasta", FASTA_GC_ALTO)
        from src.procesado_sec import procesar_genes
        df = procesar_genes(self.tmp_dir, self.output_csv, buscar_online=False)
        self.assertEqual(len(df), 2)


# ═════════════════════════════════════════════════════════════
# 3. TESTS: model_ia — entrenar_modelo
# ═════════════════════════════════════════════════════════════
class TestEntrenarModelo(unittest.TestCase):
    """Tests para entrenar_modelo de model_ia.py."""

    def setUp(self):
        from src.model_ia import entrenar_modelo
        self.entrenar_modelo = entrenar_modelo

    def test_retorna_tres_valores(self):
        """entrenar_modelo debe retornar (r2, prediccion, lista_anomalias)."""
        df = make_df_con_homologia()
        resultado = self.entrenar_modelo(df, valor_test=5000)
        self.assertEqual(len(resultado), 3)

    def test_r2_entre_0_y_1(self):
        """El R² debe estar entre 0 y 1."""
        df = make_df_con_homologia()
        r2, _, _ = self.entrenar_modelo(df, valor_test=5000)
        self.assertGreaterEqual(r2, 0.0)
        self.assertLessEqual(r2, 1.0)

    def test_prediccion_es_float(self):
        """La predicción para valor_test debe ser un número."""
        df = make_df_con_homologia()
        _, prediccion, _ = self.entrenar_modelo(df, valor_test=5000)
        self.assertIsInstance(float(prediccion), float)

    def test_sin_homologia_retorna_ceros(self):
        """Sin filas con homología, retorna (0, 0, [])."""
        df = make_df_sin_homologia()
        r2, pred, anomalias = self.entrenar_modelo(df, valor_test=5000)
        self.assertEqual(r2, 0)
        self.assertEqual(pred, 0)
        self.assertEqual(anomalias, [])

    def test_df_vacio_retorna_ceros(self):
        """Un DataFrame vacío retorna (0, 0, [])."""
        df = pd.DataFrame(columns=["gene_id", "id_homologo", "longitud_adn", "contenido_gc"])
        r2, pred, anomalias = self.entrenar_modelo(df, valor_test=5000)
        self.assertEqual(r2, 0)
        self.assertEqual(pred, 0)
        self.assertEqual(anomalias, [])

    def test_lista_anomalias_es_lista(self):
        """lista_anomalias siempre debe ser una lista."""
        df = make_df_con_homologia()
        _, _, anomalias = self.entrenar_modelo(df, valor_test=5000)
        self.assertIsInstance(anomalias, list)

    def test_prediccion_lineal_coherente(self):
        """Con datos perfectamente lineales R² debe ser 1.0."""
        df = pd.DataFrame({
            "gene_id":      [f"G{i}" for i in range(5)],
            "id_homologo":  ["NM_000001"] * 5,
            "longitud_adn": [1000, 2000, 3000, 4000, 5000],
            "contenido_gc": [60.0, 55.0, 50.0, 45.0, 40.0],
        })
        r2, _, _ = self.entrenar_modelo(df, valor_test=3000)
        self.assertAlmostEqual(r2, 1.0, places=5)


# ═════════════════════════════════════════════════════════════
# 4. TESTS: visualizer — generar_graficos
# ═════════════════════════════════════════════════════════════
class TestGenerarGraficos(unittest.TestCase):
    """Tests para generar_graficos de visualizer.py."""

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        from src.visualizer import generar_graficos
        self.generar_graficos = generar_graficos

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def _make_df(self):
        return pd.DataFrame({
            "longitud_adn": [1000, 2000, 3000, 4000],
            "contenido_gc": [52.0, 50.0, 48.0, 46.0],
        })

    def test_crea_png(self):
        """Se genera el archivo heatmap_correlacion.png."""
        self.generar_graficos(self._make_df(), self.tmp_dir)
        self.assertTrue(os.path.exists(os.path.join(self.tmp_dir, "heatmap_correlacion.png")))

    def test_retorna_dataframe_correlacion(self):
        """Retorna un DataFrame de correlación."""
        corr = self.generar_graficos(self._make_df(), self.tmp_dir)
        self.assertIsInstance(corr, pd.DataFrame)

    def test_correlacion_es_cuadrada_2x2(self):
        """La matriz de correlación debe ser 2x2."""
        corr = self.generar_graficos(self._make_df(), self.tmp_dir)
        self.assertEqual(corr.shape, (2, 2))

    def test_diagonal_es_uno(self):
        """La diagonal de la correlación debe ser 1.0."""
        corr = self.generar_graficos(self._make_df(), self.tmp_dir)
        self.assertAlmostEqual(corr.loc["longitud_adn", "longitud_adn"], 1.0)
        self.assertAlmostEqual(corr.loc["contenido_gc", "contenido_gc"], 1.0)

    def test_crea_directorio_si_no_existe(self):
        """Crea el directorio de salida si no existe."""
        nuevo_dir = os.path.join(self.tmp_dir, "subdir_nuevo")
        self.generar_graficos(self._make_df(), nuevo_dir)
        self.assertTrue(os.path.isdir(nuevo_dir))


# ═════════════════════════════════════════════════════════════
# 5. TESTS: fetch_tools — descargar_secuencia_homologa
# ═════════════════════════════════════════════════════════════
class TestFetchTools(unittest.TestCase):
    """Tests para descargar_secuencia_homologa de fetch_tools.py."""

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        # Email de prueba para que el guard de fetch_tools no bloquee los tests
        os.environ["ENTREZ_EMAIL"] = "test@test.com"

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)
        os.environ.pop("ENTREZ_EMAIL", None)

    def test_no_descarga_si_archivo_ya_existe(self):
        """Si el archivo ya existe localmente, no llama a Entrez."""
        # Crear el archivo de antemano
        filename = os.path.join(self.tmp_dir, "REF_NM_000001.fasta")
        with open(filename, "w") as f:
            f.write(">NM_000001\nATGCCC\n")

        with patch("src.fetch_tools.Entrez.efetch") as mock_efetch:
            from src.fetch_tools import descargar_secuencia_homologa
            resultado = descargar_secuencia_homologa("NM_000001", self.tmp_dir)
            mock_efetch.assert_not_called()
            self.assertEqual(resultado, filename)

    def test_crea_carpeta_destino_si_no_existe(self):
        """Crea la carpeta destino si no existe antes de descargar."""
        nueva_carpeta = os.path.join(self.tmp_dir, "reference")

        mock_record = MagicMock()
        mock_record.id = "NM_000001"
        mock_record.seq = Seq("ATGCCCAAA")

        with patch("src.fetch_tools.Entrez.efetch") as mock_efetch, \
             patch("src.fetch_tools.SeqIO.read", return_value=mock_record), \
             patch("src.fetch_tools.SeqIO.write"):
            mock_efetch.return_value.__enter__ = MagicMock(return_value=MagicMock())
            mock_efetch.return_value.__exit__ = MagicMock(return_value=False)
            from src.fetch_tools import descargar_secuencia_homologa
            descargar_secuencia_homologa("NM_000001", nueva_carpeta)
            self.assertTrue(os.path.isdir(nueva_carpeta))

    def test_error_ncbi_retorna_none(self):
        """Si Entrez lanza una excepción, retorna None sin propagar el error."""
        with patch("src.fetch_tools.Entrez.efetch", side_effect=Exception("timeout NCBI")):
            from src.fetch_tools import descargar_secuencia_homologa
            resultado = descargar_secuencia_homologa("NM_INVALIDO", self.tmp_dir)
            self.assertIsNone(resultado)


if __name__ == "__main__":
    unittest.main(verbosity=2)
