# Eventos y Triggers en Tiempo Real

Este documento define cómo Aura PM-OS reacciona a los cambios del entorno y del proyecto, orquestando la respuesta de los agentes.

## Estructura de Triggers

| Evento | Condición | Agente Activado | Artefactos Afectados | Salida Esperada | Severidad | SLA Respuesta |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Nuevo Proyecto** | Sistema recibe orden "Crear PRJ" | `02_initiation_charter` | Charter (Draft) | Borrador Charter | Baja | Inmediato |
| **Cambio Alcance** | Solicitud externa o nuevo req. | `03_scope_requirements` | WBS, Requisitos | Evaluación impacto alcance | Media | 4 horas |
| **Retraso Hito** | Hito supera fecha `baseline` | `04_planning_schedule` | Cronograma, Alerta | Propuesta replanificación, alerta de schedule | Alta | 2 horas |
| **Sobrecosto** | Gasto reportado > (Presupuesto Hito + Reserva) | `07_finance_agent` | Dashboard, Budget | Alerta financiera, forecast CPI ajustado | Alta | 2 horas |
| **Riesgo Nuevo** | Usuario o sistema loguea texto "Riesgo:..." | `05_risk_agent` | RAID Log | Score de severidad y sugerencia de mitigación | Media | 12 horas |
| **Issue Bloqueante** | Issue marcado `estado: Bloqueado` | `01_pmo_orchestrator` | Logs, Dashboards | Escalamiento inmediato al Sponsor (Decision Request) | Crítica | Inmediato |
| **Recurso Offline** | Disponibilidad reportada < 50% de plan | `06_resources_agent` | Plan Recursos | Alerta de cuello de botella, sugerencia reasignación | Media | 24 horas |
| **Decisión Vencida** | Fecha `resolucion` de Issue/CR expirada | `10_monitoring_control` | Dashboards, Status | Notificación automática al owner | Media | Diaria (batch) |
| **Pase de Fase** | Intento de cambiar `fase_actual` | `09_quality_governance` | Reporte Calidad | Validación de Go/No-Go | Alta | 1 hora |
| **Comms Externas** | Cierre de mes / Hito Mayor completado | `08_stakeholders_comms` | Plan Comms | Borrador de correo ejecutivo al comité | Baja | 24 horas |

## Manejo de Excepciones

Si un agente falla en generar la salida esperada dentro del SLA, el `01_pmo_orchestrator` toma el control y envía una notificación `Fallback` al Project Manager (humano) indicando:
> "SLA Expirado para [Evento]. Se requiere intervención humana para analizar [Condición]."