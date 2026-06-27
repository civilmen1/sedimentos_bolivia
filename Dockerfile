FROM python:3.11-slim-bookworm

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcairo2 libpango-1.0-0 libpangocairo-1.0-0 \
    libgdk-pixbuf2.0-0 libffi-dev shared-mime-info \
    fonts-liberation fontconfig \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 7860

# Cachés escribibles para numba y matplotlib (en HF Spaces el HOME puede ser
# de solo lectura). Evita que la compilación JIT de numba falle o se corrompa.
ENV NUMBA_CACHE_DIR=/tmp/numba_cache \
    MPLCONFIGDIR=/tmp/mpl_cache \
    XDG_CACHE_HOME=/tmp/xdg_cache
RUN mkdir -p /tmp/numba_cache /tmp/mpl_cache /tmp/xdg_cache && chmod 777 /tmp/numba_cache /tmp/mpl_cache /tmp/xdg_cache

# Un solo worker: el geoproceso (pysheds/numba) es intensivo en memoria (~1.2 GB
# por cuenca) y la compilación JIT concurrente entre varios workers corrompe la
# caché de numba (KeyError en funciones generadas). Con 1 worker la compilación
# es secuencial y estable. Timeout amplio para el primer geoproceso (cacheado).
CMD ["gunicorn", "app:app", "--workers", "1", "--threads", "4", "--bind", "0.0.0.0:7860", "--timeout", "300"]
