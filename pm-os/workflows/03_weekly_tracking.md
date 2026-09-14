# Flujo 03: Seguimiento Semanal (Weekly Tracking)

**Disparador:** Viernes 14:00 hrs (Programado)
**Agentes Principales:** `10_monitoring_control`, `08_stakeholders_comms`, `01_pmo_orchestrator`

## Secuencia Paso a Paso
1. **Recopilación:** `10_monitoring_control` consulta el avance (`schedule.json`), gastos (`budget.json`) y el RAID log.
2. **Cálculo de KPIs:** Calcula SPI, CPI, porcentaje de hitos completados y volatilidad de alcance.
3. **Determinación del Health Score:** Evalúa los semáforos de cada área y determina el color general del proyecto.
4. **Generación del Reporte:** Crea el archivo `status_report.md` en base al template.
5. **Redacción Ejecutiva:** `08_stakeholders_comms` extrae un resumen de 3 viñetas para el Dashboard Ejecutivo.
6. **Publicación:** `12_living_documentation` actualiza el repositorio. El Orquestador envía notificación (Slack/Email) al PM con el borrador del reporte para envío.