# MANUAL DE USO: Aura PM-OS con Claude Code

Este manual explica cómo tú (el usuario) debes interactuar con **Claude Code** para que opere como el `01_pmo_orchestrator` y ejecute los flujos de Aura PM-OS en este repositorio.

## Concepto Core
En lugar de hacer clics en una interfaz gráfica, en Aura PM-OS **conversas con tu terminal usando Claude Code**. Tú le das el input inicial, y Claude Code (actuando como orquestador) lee las reglas en la carpeta `pm-os/`, asume las personalidades de los agentes, y genera/actualiza los archivos Markdown y JSON.

---

## 1. INICIAR UN NUEVO PROYECTO
Cuando tengas una nueva solicitud de negocio, ejecuta Claude Code en tu terminal y dale esta orden:

> **Usuario:** "Actúa como el Agente 02 (Initiation & Charter). Tengo este correo del Sponsor: 'Necesitamos migrar los servidores on-premise a AWS antes de fin de año con un presupuesto de $50k'. Inicia el Flujo 01 y genérame el borrador del Project Charter en una carpeta nueva en `/pm-os/projects/PRJ-AWS-001/`."

**Qué hará el sistema:**
Creará la carpeta, usará el template `project_charter.md`, llenará los datos y te lo mostrará.

---

## 2. PLANIFICACIÓN Y DESGLOSE
Una vez aprobado el Charter, debes ordenar la creación del desglose (WBS).

> **Usuario:** "Actúa como el Agente 03 y 04. Lee el Charter de PRJ-AWS-001. Necesito que generes el archivo `wbs.md` con las fases principales de la migración. Luego, genera un archivo `schedule.json` estimando los hitos."

**Qué hará el sistema:**
Creará el desglose estructurado y el JSON de tiempos, aplicando las reglas metodológicas de los prompts.

---

## 3. GESTIÓN DE RIESGOS (RAID LOG)
Para ingresar un riesgo nuevo al sistema, no tienes que abrir un Excel. Dile a Claude:

> **Usuario:** "Activa el Agente 05 (Risk). Hay un nuevo riesgo en PRJ-AWS-001: El equipo de seguridad podría demorar la revisión de arquitectura. Evalúa probabilidad e impacto, actualiza el RAID log (`raid_log.md`) y genera una salida en formato `risk_evaluation.md`."

**Qué hará el sistema:**
Calculará el score de severidad (P x I), actualizará el Log y te dará un reporte con recomendaciones de mitigación.

---

## 4. ACTUALIZACIÓN SEMANAL (STATUS REPORT)
Los viernes, para generar tu reporte de estado sin esfuerzo:

> **Usuario:** "Actúa como Agente 10 (Monitoring) y Agente 08 (Comms). Ejecuta el Flujo 03 para el PRJ-AWS-001. Asume que vamos en tiempo (SPI = 1.0) pero gastamos un poco más de lo esperado (CPI = 0.90). Genera el `status_report.md` semanal y redacta un correo ejecutivo usando la plantilla `executive_summary.md`."

**Qué hará el sistema:**
Consolidará los datos, calculará el Health Score (según las reglas de `KPI_HEALTH_SCORE.md`) y te dejará el reporte listo para enviar.

---

## 5. CONTROL DE CAMBIOS
Si el Sponsor pide algo nuevo que altera el alcance:

> **Usuario:** "Actúa como PMO Orchestrator. El Sponsor quiere añadir la migración de la base de datos secundaria al PRJ-AWS-001. Ejecuta el Flujo 04 de Control de Cambios. Analiza el impacto, dímelo y genera el documento `change_request.md`."

**Qué hará el sistema:**
Evaluará el impacto en tiempo y costo (actuando como los agentes 04 y 07), y te preparará el CR formal para firma.

---

## RESUMEN DE OPERACIÓN
1. **El repositorio es la base de datos:** Todo vive en `/pm-os/projects/`.
2. **Los Prompts son las leyes:** Claude Code lee `/pm-os/prompts/` para saber *cómo* pensar antes de hacer algo.
3. **Tú eres el Juez:** Claude te prepara recomendaciones (`recommendation.md`), solicitudes de decisión (`decision_request.md`) y alertas (`alert.md`). Tú apruebas y luego le dices: "Cambio aprobado, actualiza las líneas base."