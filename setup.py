import os
import sys


def crear_proyecto(ruta_base="."):
    # 1. Definición de la estructura de carpetas
    carpetas = [
        "data/raw",
        "data/processed",
        "src",
        "scripts",
        "logs",
        "graficos",
        "resultados",
        "notebooks"
    ]

    for carpeta in carpetas:
        ruta_carpeta = os.path.join(ruta_base, carpeta)
        os.makedirs(ruta_carpeta, exist_ok=True)
        print(f"Creada: {ruta_carpeta}")

    # 2. Definicion del contenido del .gitignore
    gitignore_content = """# Entorno virtual
venv/
ENV/
.env

# Python compilado
__pycache__/
*.py[cod]
*$py.class

# Datos y Notebooks
.ipynb_checkpoints
data/raw/*
!data/raw/.gitkeep
"""

    # 3. Creacion del archivo .gitignore
    with open(os.path.join(ruta_base, '.gitignore'), 'w') as f:
        f.write(gitignore_content)
    print("Archivo .gitignore generado con éxito.")

    # 4. Creacion de un archivo .gitkeep en data/raw
    with open(os.path.join(ruta_base, 'data', 'raw', '.gitkeep'), 'w') as f:
        pass

# Usamos doble signo de igual para comparar
if __name__ == "__main__":
    ruta = sys.argv[1] if len(sys.argv) > 1 else "."
    crear_proyecto(ruta)