# Conexión directa con Google Earth Engine (cuenta de servicio)

La aplicación se conecta a Google Earth Engine (GEE) **desde el servidor**, sin
que el usuario tenga que copiar código a Colab. Para ello necesita las
credenciales de una **cuenta de servicio** de Google Cloud.

Cuando las credenciales están configuradas, los mapas del informe y de la
galería `/maps` se obtienen con **datos satelitales reales** (Sentinel-2, SRTM,
ESA WorldCover, JRC). Si no hay credenciales, la app sigue funcionando con
**datos sintéticos de demostración** (no se rompe nada).

## 1. Crear la cuenta de servicio (una sola vez)

1. Entra a [Google Cloud Console](https://console.cloud.google.com/) y crea o
   selecciona un proyecto.
2. Habilita la **Earth Engine API**:
   `APIs & Services → Library → "Earth Engine API" → Enable`.
3. Registra el proyecto en Earth Engine:
   [https://code.earthengine.google.com/register](https://code.earthengine.google.com/register)
4. Crea la cuenta de servicio:
   `IAM & Admin → Service Accounts → Create Service Account`.
   - Nombre, por ejemplo: `sedimentos-gee`.
   - No es necesario asignar roles especiales para usar GEE en modo lectura.
5. Genera la clave JSON:
   `(la cuenta creada) → Keys → Add Key → Create new key → JSON`.
   Se descargará un archivo `*.json`.
6. Autoriza la cuenta de servicio en Earth Engine (importante):
   [https://code.earthengine.google.com/](https://code.earthengine.google.com/)
   → en **Assets/Settings** o en la página de registro, agrega el
   `client_email` de la cuenta de servicio como usuario autorizado del proyecto
   de Earth Engine.

## 2. Configurar el secreto en el servidor

### Hugging Face Space
`Settings → Variables and secrets → New secret`:

| Nombre | Valor |
|--------|-------|
| `EE_SERVICE_ACCOUNT_JSON` | *(pega el contenido completo del archivo JSON)* |

> Pega el JSON **completo** (todo el contenido del archivo descargado, incluido
> `client_email`, `private_key`, `project_id`, etc.).

### Local / otro servidor
```bash
export EE_SERVICE_ACCOUNT_JSON="$(cat /ruta/a/clave.json)"
```

Alternativas soportadas:
- `EE_SERVICE_ACCOUNT_FILE` (ruta a un archivo .json) + `EE_SERVICE_ACCOUNT_EMAIL`
- Credenciales personales locales (`earthengine authenticate`) si no se define
  ninguna variable.

## 3. Verificar

Al iniciar, el log del servidor debe mostrar:

```
GEE inicializado con cuenta de servicio: sedimentos-gee@tu-proyecto.iam.gserviceaccount.com
```

En la galería `/maps` aparecerá el indicador **🛰 Datos satelitales reales —
Google Earth Engine**. Si en cambio ves **⚠ Datos de demostración**, revisa que
el secreto esté bien pegado y que la cuenta esté autorizada en Earth Engine.

### Diagnóstico en vivo: `/gee_status`

Abre en el navegador (sustituye por la URL de tu Space):

```
https://<tu-space>.hf.space/gee_status?probe=1
```

Devuelve un JSON con el estado exacto:

- `data_mode`: `"real"` o `"synthetic"`.
- `ready`: `true` si GEE se inicializó.
- `env_EE_SERVICE_ACCOUNT_JSON_present`: `true` si el secreto está definido.
- `init_error`: motivo si la inicialización falló (JSON mal pegado, cuenta no
  registrada, etc.).
- `probe` / `probe_error`: con `?probe=1` hace una consulta real al DEM SRTM;
  `"ok"` confirma acceso efectivo a datos. Si `probe` falla pero `ready` es
  `true`, normalmente la cuenta de servicio **no está autorizada en Earth
  Engine** (paso 6).

Casos típicos:

| Síntoma en `/gee_status` | Causa | Solución |
|--------------------------|-------|----------|
| `env_..._present: false` | Secreto no configurado | Añade `EE_SERVICE_ACCOUNT_JSON` (paso 2) |
| `init_error` con "JSON no es válido" | Pegaste mal el contenido | Pega el archivo `.json` completo |
| `ready: true`, `probe: failed` | Cuenta sin acceso a EE | Registra el proyecto y autoriza la cuenta (pasos 3 y 6) |
| `data_mode: real`, `probe: ok` | Todo correcto | ✅ |

## Seguridad

- **Nunca** subas el archivo JSON ni el contenido de `EE_SERVICE_ACCOUNT_JSON`
  al repositorio. Configúralo siempre como secreto/variable de entorno.
- El `.gitignore` del proyecto ignora archivos `*.json` de credenciales.
