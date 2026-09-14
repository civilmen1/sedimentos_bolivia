# Modelo de Gobernanza de Aura PM-OS

Este documento establece las reglas del juego para asegurar que la automatización provista por los agentes no genere caos, descontrol ni compromisos no autorizados.

## 1. Niveles de Autoridad del Sistema

El sistema opera bajo el principio de **autonomía controlada**. Se definen los siguientes niveles de autoridad para los agentes:
- **Total Autonomía (Informativa/Documental):** El sistema puede actualizar dashboards, logs (RAID, Decisiones), resúmenes ejecutivos, calcular varianzas y vincular artefactos sin intervención humana.
- **Autonomía Propositiva (Recomendaciones):** El sistema puede proponer planes de contingencia, sugerir replanificaciones de cronograma, redactar borradores de Change Requests o sugerir reasignación de recursos. *Requiere confirmación humana.*
- **Restricción Estricta (Aprobaciones Finales):** El sistema **NUNCA** puede:
  - Aprobar un aumento de presupuesto.
  - Extender la fecha de fin (Línea Base del Cronograma).
  - Cerrar un Riesgo/Issue clasificado como Crítico.
  - Eliminar alcance (Scope Baseline).

## 2. Tipos de Decisiones

Las decisiones procesadas por Aura PM-OS se clasifican en:
1.  **Informativas:** Cambios menores (ej. actualizar un estado de tarea a "En Progreso"). El agente las ejecuta y registra silenciosamente.
2.  **Recomendadas:** El sistema analiza datos (ej. tendencia de retraso) y recomienda una acción (ej. "Fast-tracking"). Requiere revisión del PM.
3.  **Sujetas a Confirmación:** El sistema redacta la actualización completa (ej. enviar minuta a stakeholders). Espera un "Aceptar" del humano.
4.  **Ejecutivas (Bloqueadas para la IA):** Decisiones de cambio de Baseline, Go/No-Go. Deben ser tomadas por el Sponsor/Steering Committee y registradas manualmente en el sistema por el PM.
5.  **De Cambio Controlado:** Ingresan por el Flujo de Change Control.

## 3. Mecanismo de Escalamiento

El sistema evalúa continuamente los indicadores y actúa según la severidad:
- **Informa:** Cambios menores (Log).
- **Alerta:** Desviación moderada (ej. SPI < 0.95, CPI < 0.95). Notifica al PM en el dashboard diario.
- **Eleva:** Si la alerta no se atiende en el SLA definido (ej. 48h) o si el KPI cae a zona roja, notifica al PMO/Sponsor.
- **Bloquea:** Si un Gate de calidad no se cumple, el Agente de Gobernanza impide que el proyecto pase a la siguiente fase metodológica.
- **Solicita Decisión:** Envía un formato estructurado (`decision_request.md`) requiriendo intervención humana explícita.

## 4. Flujo de Control de Cambios (Change Management)

Para cualquier alteración de las líneas base (Alcance, Tiempo, Costo):
1.  **Detección/Solicitud:** Inicia manual (Stakeholder) o automática (Agente detecta desalineación).
2.  **Evaluación de Impacto:** Los agentes (Alcance, Cronograma, Riesgo, Finanzas) analizan en paralelo cómo afecta el cambio.
3.  **Propuesta:** Se genera un `change_request.md` consolidado.
4.  **Revisión:** El PM revisa el impacto.
5.  **Aprobación/Rechazo:** El Sponsor (o Change Control Board) decide fuera del sistema. El PM registra la decisión.
6.  **Actualización de Línea Base:** Si es aprobado, los agentes actualizan la WBS, el Cronograma y el Presupuesto (Baseline V2).
7.  **Comunicación:** El Agente de Comms informa a los afectados.

## 5. Trazabilidad de Decisiones y Cambios

Toda alteración en un artefacto core debe incluir la siguiente metainformación (Log):
- `Timestamp:` (ISO 8601)
- `Origen:` (Usuario Humano o Nombre del Agente)
- `Evidencia:` (Dato o trigger que provocó el cambio)
- `Impacto Evaluado:` (Alto/Medio/Bajo en Costo/Tiempo)
- `Recomendación:` (Texto del sistema)
- `Decisión Humana:` (Aprobado/Rechazado/Modificado)
- `Responsable:` (ID de la persona que aprobó)
- `Estado:` (Implementado/Pendiente)