FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 7860

CMD ["gunicorn", "app:app", "--workers", "2", "--bind", "0.0.0.0:7860", "--timeout", "120"]
