import os
import sys


def crear_proyecto(ruta_base="."):
    # 1. Definición de la estructura de carpetas (Actualizada y limpia)
    carpetas = [
        "data/raw",
        "data/processed",
        "src",           # Tu carpeta única de código fuente
        "logs",
        "graficos",
        "resultados"
    ]

    for carpeta in carpetas:
        ruta_carpeta = os.path.join(ruta_base, carpeta)
        os.makedirs(ruta_carpeta, exist_ok=True)
        print(f"Creada: {ruta_carpeta}")

    # 2. Definición del contenido del .gitignore (Incluye protección de entornos y VS Code)
    gitignore_content = """# Entorno virtual
venv/
ENV/
.env

# Editores de código e IDEs (¡Evita archivos basura!)
*.code-workspace
.vscode/

# Python compilado
__pycache__/
*.py[cod]
*$py.class

# Datos y Notebooks
.ipynb_checkpoints
data/raw/*
!data/raw/.gitkeep
"""

    # 3. Creación del archivo .gitignore
    with open(os.path.join(ruta_base, '.gitignore'), 'w') as f:
        f.write(gitignore_content)
    print("Archivo .gitignore generado con éxito.")

    # 4. Creación de un archivo .gitkeep en data/raw
    with open(os.path.join(ruta_base, 'data', 'raw', '.gitkeep'), 'w') as f:
        pass
    print("Archivo .gitkeep generado en data/raw.")

    # 5. Creación automática del __init__.py en src (¡El toque profesional!)
    with open(os.path.join(ruta_base, 'src', '__init__.py'), 'w') as f:
        pass
    print("Archivo src/__init__.py generado con éxito.")


if __name__ == "__main__":
    ruta = sys.argv[1] if len(sys.argv) > 1 else "."
    crear_proyecto(ruta)