# Calculadora de Transporte de Sedimentos (GEE + Campo)

Este proyecto es una aplicación web desarrollada con Flask que integra Google Earth Engine (GEE) para el cálculo de variables hidráulicas y el transporte de sedimentos, basándose en la metodología de Ramirez Quispe (2021).

## Tecnologías Utilizadas

- **Backend:** Python, Flask
- **Análisis Espacial:** Google Earth Engine API
- **Cálculo Numérico:** NumPy, Pandas
- **Frontend:** HTML5, Tailwind CSS
- **Servidor de Producción:** Gunicorn

## Instalación y Configuración Local

1.  **Clonar el repositorio:**
    ```bash
    git clone <url-del-repositorio>
    cd <nombre-del-directorio>
    ```

2.  **Instalar dependencias:**
    Se recomienda usar un entorno virtual.
    ```bash
    pip install -r requirements.txt
    ```

3.  **Autenticar Google Earth Engine:**
    Para ejecución local, necesitas autenticarte con tu cuenta de Google:
    ```bash
    earthengine authenticate
    ```

4.  **Ejecutar la aplicación:**
    ```bash
    python app.py
    ```
    La aplicación estará disponible en `http://127.0.0.1:5000`.

## Despliegue en Producción

Para desplegar esta aplicación de forma pública (por ejemplo, en Render, Railway o Heroku), sigue estos pasos:

### 1. Configuración de Credenciales de GEE

Dado que en producción no puedes realizar una autenticación interactiva, debes usar una **Cuenta de Servicio (Service Account)** de Google Cloud:

1.  Crea un proyecto en [Google Cloud Console](https://console.cloud.google.com/).
2.  Habilita la API de Google Earth Engine.
3.  Crea una **Service Account** y descarga la clave en formato **JSON**.
4.  Copia el contenido completo de ese archivo JSON.
5.  En tu plataforma de despliegue, crea una variable de entorno llamada `GEE_SERVICE_ACCOUNT_JSON` y pega el contenido del JSON como su valor.

### 2. Archivos de Configuración Incluidos

- **Dockerfile:** Permite crear una imagen de contenedor con todas las dependencias.
- **Procfile:** Utilizado por plataformas tipo PaaS para identificar el comando de inicio.
- **requirements.txt:** Lista de todas las librerías necesarias.

### 3. Comando de Inicio
El comando recomendado para producción es:
```bash
gunicorn --bind 0.0.0.0:$PORT app:app
```

## Pruebas

Para ejecutar las pruebas unitarias:
```bash
python3 -m pytest tests/test_models.py
export PYTHONPATH=$PYTHONPATH:. && python3 tests/test_gee.py
```

## Autor
Basado en "Transporte de Sedimentos con Python" por Ing. Robert Ramirez Quispe (2021).
