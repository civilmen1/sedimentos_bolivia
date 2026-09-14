# Catálogo de Agentes (Multi-Agent Architecture)

Este documento describe la red de 12 agentes especializados que operan Aura PM-OS, coordinados por el PMO Orchestrator.

## 1. Agente Orquestador PMO (`01_pmo_orchestrator`)
*   **Misión:** Coordinar a todos los agentes, consolidar el estado del proyecto, priorizar acciones y resolver dependencias.
*   **Responsabilidades:** Asignar tareas a otros agentes, compilar la vista integral, asegurar el flujo de trabajo continuo.
*   **Entradas:** Señales de eventos, outputs de otros agentes, directivas del PM.
*   **Salidas:** Órdenes de ejecución, consolidación de estado, alertas de sistema.
*   **Decisiones:** Recomendar acciones, pausar flujos no críticos. No puede aprobar cambios de línea base.
*   **Frecuencia:** Continua (Event-driven) / Diaria.
*   **Triggers:** Eventos de sistema, finalización de tareas de otros agentes.
*   **Dependencias:** Todos los agentes.
*   **Artefactos:** Dashboard Ejecutivo, Log de Orquestación.

## 2. Agente de Inicio y Charter (`02_initiation_charter`)
*   **Misión:** Transformar una idea o solicitud en un Project Charter estructurado.
*   **Responsabilidades:** Identificar objetivos, alcance inicial, restricciones y beneficios.
*   **Entradas:** Business Case, solicitudes de negocio, notas de kickoff.
*   **Salidas:** Project Charter en borrador.
*   **Decisiones:** Recomendar criterios de éxito. No aprueba el Charter final.
*   **Frecuencia:** Al inicio del proyecto o ante pivotes mayores.
*   **Triggers:** Creación de proyecto, nueva solicitud de sponsor.
*   **Dependencias:** Stakeholders, Orquestador.
*   **Artefactos:** Project Charter.

## 3. Agente de Alcance y Requisitos (`03_scope_requirements`)
*   **Misión:** Capturar y controlar requisitos, backlog y la WBS.
*   **Responsabilidades:** Versionar requisitos, detectar desalineación de alcance.
*   **Entradas:** Charter, reuniones de toma de requisitos, change requests.
*   **Salidas:** WBS, Registro de Requisitos actualizado.
*   **Decisiones:** Identificar scope creep. No aprueba cambios de alcance sin validación.
*   **Frecuencia:** Semanal o bajo demanda (cambio).
*   **Triggers:** Nuevo requisito, cambio solicitado.
*   **Dependencias:** Cronograma, Recursos.
*   **Artefactos:** Scope Statement, WBS, Registro de Requisitos.

## 4. Agente de Planificación y Cronograma (`04_planning_schedule`)
*   **Misión:** Gestionar hitos, ruta crítica y forecasts temporales.
*   **Responsabilidades:** Alertar sobre retrasos (slippages), proyectar escenarios.
*   **Entradas:** WBS, estimaciones, dependencias, actualizaciones de progreso.
*   **Salidas:** Cronograma maestro, Alertas de retraso.
*   **Decisiones:** Proponer replanificación. No modifica línea base sin aprobación.
*   **Frecuencia:** Continua (Tracking) / Semanal (Review).
*   **Triggers:** Tarea retrasada, actualización de avance.
*   **Dependencias:** Alcance, Recursos.
*   **Artefactos:** Cronograma Maestro, Registro de Hitos.

## 5. Agente de Riesgos (`05_risk_agent`)
*   **Misión:** Identificar, evaluar y monitorear riesgos e issues.
*   **Responsabilidades:** Mantener el RAID log, evaluar triggers, sugerir respuestas.
*   **Entradas:** Contexto del proyecto, lecciones aprendidas, alertas de cronograma/finanzas.
*   **Salidas:** Evaluación de riesgo, mitigaciones sugeridas.
*   **Decisiones:** Escalar riesgos a issues. No cierra riesgos críticos unilateralmente.
*   **Frecuencia:** Semanal / Event-driven (Nuevo riesgo).
*   **Triggers:** Desviación detectada, reporte externo.
*   **Dependencias:** Todos los agentes operativos.
*   **Artefactos:** Registro de Riesgos, Registro de Issues (RAID Log).

## 6. Agente de Recursos (`06_resources_agent`)
*   **Misión:** Gestionar capacidad, asignación y cuellos de botella de talento.
*   **Responsabilidades:** Identificar sobrecarga, recomendar staffing.
*   **Entradas:** Cronograma, perfiles, reportes de capacidad.
*   **Salidas:** Alertas de recursos, mapa de carga.
*   **Decisiones:** Sugerir reasignación. No contrata ni despide (obviamente).
*   **Frecuencia:** Semanal.
*   **Triggers:** Cambio en cronograma, indisponibilidad reportada.
*   **Dependencias:** Cronograma, Finanzas.
*   **Artefactos:** Plan de Recursos.

## 7. Agente Financiero (`07_finance_agent`)
*   **Misión:** Controlar presupuesto, burn rate y forecast.
*   **Responsabilidades:** Rastrear OPEX/CAPEX, emitir alertas por desviación.
*   **Entradas:** Presupuesto base, gastos reales, facturas, horas reportadas.
*   **Salidas:** Alertas de costo, forecast financiero.
*   **Decisiones:** Recomendar uso de reservas. No aprueba sobrecostos.
*   **Frecuencia:** Semanal / Mensual.
*   **Triggers:** Gasto ingresado, hito superado (pago).
*   **Dependencias:** Cronograma, Recursos.
*   **Artefactos:** Forecast Financiero.

## 8. Agente de Stakeholders y Comunicaciones (`08_stakeholders_comms`)
*   **Misión:** Gestionar mapa de stakeholders y estrategias de comunicación.
*   **Responsabilidades:** Redactar resúmenes ejecutivos, monitorear satisfacción.
*   **Entradas:** Análisis de actores, estado del proyecto.
*   **Salidas:** Plan de comunicación, borradores de correos/updates.
*   **Decisiones:** Sugerir intervención ante stakeholder insatisfecho.
*   **Frecuencia:** Semanal / Event-driven.
*   **Triggers:** Nuevo stakeholder, hito importante.
*   **Dependencias:** Orquestador, Riesgos.
*   **Artefactos:** Mapa de Stakeholders, Plan de Comunicaciones.

## 9. Agente de Calidad y Gobernanza (`09_quality_governance`)
*   **Misión:** Verificar cumplimiento metodológico (PMBOK) y Definition of Done.
*   **Responsabilidades:** Auditar consistencia, validar gates.
*   **Entradas:** Artefactos entregados, estándares de la PMO.
*   **Salidas:** Reporte de calidad, bloqueos por no cumplimiento.
*   **Decisiones:** Bloquear pase de fase si falta calidad.
*   **Frecuencia:** Por hito / Mensual.
*   **Triggers:** Intento de cierre de fase/entregable.
*   **Dependencias:** Orquestador.
*   **Artefactos:** Checklist de Calidad.

## 10. Agente de Seguimiento y Control (`10_monitoring_control`)
*   **Misión:** Consolidar KPIs, semáforos y varianzas (Health Score).
*   **Responsabilidades:** Generar el reporte de estado global y tendencias.
*   **Entradas:** Datos de Cronograma, Finanzas, Riesgos y Calidad.
*   **Salidas:** Dashboards operativos y ejecutivos, Health Score.
*   **Decisiones:** Cambiar el semáforo del proyecto (Verde/Ámbar/Rojo).
*   **Frecuencia:** Semanal.
*   **Triggers:** Cierre de semana, solicitud de reporte.
*   **Dependencias:** Todos los agentes operativos.
*   **Artefactos:** Status Report Semanal, Dashboard Operativo.

## 11. Agente de Cierre y Lecciones Aprendidas (`11_closure_lessons`)
*   **Misión:** Validar el cierre, capturar beneficios y registrar lecciones.
*   **Responsabilidades:** Preparar handover, documentar retrospectivas.
*   **Entradas:** Estado final, comentarios de stakeholders, RAID log histórico.
*   **Salidas:** Documento de cierre, repositorio de lecciones aprendidas.
*   **Decisiones:** Recomendar cierre formal. No libera el presupuesto sin firma.
*   **Frecuencia:** Al final del proyecto (o fase).
*   **Triggers:** Hito final completado.
*   **Dependencias:** Gobernanza, Finanzas.
*   **Artefactos:** Documento de Cierre, Lecciones Aprendidas.

## 12. Agente de Documentación Viva (`12_living_documentation`)
*   **Misión:** Mantener todos los artefactos sincronizados y versionados.
*   **Responsabilidades:** Actualizar metadatos, cruzar referencias.
*   **Entradas:** Cualquier cambio en cualquier artefacto.
*   **Salidas:** Archivos actualizados (Markdown/JSON).
*   **Decisiones:** Fusionar cambios menores. Requiere revisión para colisiones.
*   **Frecuencia:** Continua.
*   **Triggers:** Guardado/Actualización de un documento por otro agente.
*   **Dependencias:** Todos.
*   **Artefactos:** Todo el repositorio `projects/`.