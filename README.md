---
title: Sedimentos Bolivia
emoji: 🏔️
colorFrom: blue
colorTo: green
sdk: docker
app_port: 7860
pinned: false
license: unlicense
short_description: Calculadora de transporte de sedimentos fluviales - Bolivia
---

# Calculadora de Transporte de Sedimentos — Bolivia

Herramienta técnica para el cálculo de transporte de sedimentos fluviales en Bolivia, basada en el libro **"Transporte de Sedimentos con Python – Tomo 1"** de Ing. Robert Ramirez Quispe (2021).

## Modelos implementados

- **Meyer-Peter & Müller (1948)** — Carga de fondo
- **Engelund-Hansen (1967)** — Carga total
- **Van Rijn (1984)** — Carga de fondo con parámetro T

## Parámetros calculados

- Velocidad de caída wₛ (Van Rijn 1993)
- Parámetro dimensional D*
- Criterio de Shields (θ₀ vs θ_c)
- Número de Froude y Rouse
- Modo de transporte (washload / mixto / fondo)
- Análisis de sensibilidad transporte vs tirante

## Uso

1. Ingrese las coordenadas del punto de muestreo (o haga clic en el mapa)
2. Ingrese datos granulométricos (d₅₀, d₉₀) y propiedades físicas
3. Ingrese datos hidráulicos de campo (tirante y, velocidad v)
4. Ingrese la pendiente del lecho S
5. Haga clic en **CALCULAR**

## Referencia

Ramirez Quispe, R. (2021). *Transporte de Sedimentos con Python – Tomo 1*. Bolivia.
