# Flujo 07: Re-planificación por Retrasos (Re-planning)

**Disparador:** SPI (Schedule Performance Index) cae por debajo de 0.85 (Alerta Roja).
**Agentes Principales:** `04_planning_schedule`, `01_pmo_orchestrator`

## Secuencia Paso a Paso
1. **Detección:** `10_monitoring_control` detecta el SPI bajo y lanza alerta.
2. **Análisis de Ruta Crítica:** `04_planning_schedule` identifica el cuello de botella o tarea retrasada que empuja la fecha fin.
3. **Escenarios:** El agente propone dos escenarios:
   - *Crashing:* Añadir recursos (consulta a `07_finance_agent` por costo y a `06_resources_agent` por disponibilidad).
   - *Fast-Tracking:* Paralelizar tareas (consulta dependencias en la WBS).
4. **Decisión:** Presenta las opciones al PM. El PM aprueba la mejor vía de acción y ajusta el forecast.