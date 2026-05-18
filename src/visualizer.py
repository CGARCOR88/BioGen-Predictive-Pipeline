import os
import matplotlib
matplotlib.use('Agg')  # backend no interactivo (sin pantalla / CI)
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def generar_graficos(df: pd.DataFrame, output_dir: str):
    """Genera y guarda el heatmap de correlación."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Configurar el estilo visual
    sns.set_theme(style='whitegrid')
    plt.figure(figsize=(8, 6))

    # Calcular correlación
    corr = df[['longitud_adn', 'contenido_gc']].corr()

    # Dibujar heatmap
    sns.heatmap(corr, annot=True, cmap='coolwarm', fmt='.2f')
    plt.title('Análisis de Correlación: Longitud vs Contenido GC')

    # Guardar
    save_path = os.path.join(output_dir, 'heatmap_correlacion.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close() # Importante: cerrar para liberar memoria
    
    return corr