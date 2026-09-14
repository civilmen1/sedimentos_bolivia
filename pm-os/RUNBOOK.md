# Runbook Operativo (Cadencia de Ejecución)

Este documento describe cómo se utiliza Aura PM-OS en el día a día. Define qué hacen los agentes y qué hace el humano.

## 1. Rutina Diaria (Stand-up / Daily)
**Duración humana:** ~15 minutos.
**Nivel de automatización:** Alto (Event-driven).

*   **Sistema (08:00 AM):**
    *   `01_pmo_orchestrator` procesa todos los eventos nocturnos.
    *   `10_monitoring_control` revisa si hay `decisiones vencidas` o `issues estancados`.
    *   Se genera y envía el Dashboard Diario (Solo Alertas y Bloqueos) al PM.
*   **Project Manager:**
    *   Revisa alertas rojas/ámbar en el Dashboard.
    *   Responde a los `decision_request.md` generados por los agentes.
    *   Desbloquea al equipo técnico.

## 2. Rutina Semanal (Tracking & Status)
**Duración humana:** ~45 minutos (Sustituye 3-4 horas de recopilación manual).
**Nivel de automatización:** Medio (Requiere revisión de calidad humana).

*   **Sistema (Jueves PM):**
    *   `04_planning_schedule` y `07_finance_agent` solicitan el avance real (% Completado, Horas, Facturas).
    *   Si los datos están en herramientas integradas (ej. Jira/ERP), el sistema ingesta automáticamente.
*   **Sistema (Viernes AM):**
    *   Se ejecuta el Flujo 03 (Tracking Semanal).
    *   `10_monitoring_control` calcula el Health Score global.
    *   `05_risk_agent` evalúa el RAID log y sugiere mitigaciones vigentes.
    *   Se compila el `status_report.md` preliminar.
*   **Project Manager:**
    *   Revisa el `status_report.md`.
    *   Ajusta matices y valida las recomendaciones de los agentes (ej. acepta sugerencia de Fast-Tracking).
    *   Envía el reporte final consolidado.

## 3. Rutina Mensual (Steering / Forecast)
**Duración humana:** ~2 horas.
**Nivel de automatización:** Asistido.

*   **Sistema (Día 1 de cada mes):**
    *   `07_finance_agent` genera el cierre de Forecast (EAC, VAC).
    *   `06_resources_agent` proyecta la capacidad requerida para el siguiente mes vs. el Cronograma maestro.
    *   `09_quality_governance` audita la carpeta del proyecto y reporta artefactos huérfanos o desactualizados.
*   **Project Manager / PMO:**
    *   Revisión de tendencias (Snapshots históricos).
    *   Presentación de resultados al Sponsor/Comité de Dirección (basado en el output de `08_stakeholders_comms`).
    *   Ajustes metodológicos a las reglas de los agentes (si es necesario).