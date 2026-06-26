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

# timeout amplio: la primera petición puede pagar la compilación JIT del
# geoproceso (mitigada por el warmup en segundo plano). Los mapas se cachean.
CMD ["gunicorn", "app:app", "--workers", "2", "--bind", "0.0.0.0:7860", "--timeout", "300"]
