# Modelo de Datos y Memoria (Data Model)

Aura PM-OS persiste el estado del proyecto utilizando archivos estáticos (Markdown estructurado y JSON/YAML) para facilitar la lectura por humanos y por los LLMs subyacentes de los agentes.

## Esquema de Persistencia
- **Formato principal:** Markdown con Frontmatter (YAML) para metadatos clave, permitiendo que el cuerpo del archivo mantenga la narrativa y contexto.
- **Datos relacionales pesados:** Archivos JSON (ej. `schedule.json`, `budget.json`) dentro de la carpeta del proyecto.
- **Versionado:** Controlado vía Git. Cada *commit* representa un snapshot histórico del estado del proyecto.

## Entidades Principales (Esquemas Mínimos)

### 1. Proyecto (`project.md`)
*   `id_proyecto`: (String, ej. "PRJ-001")
*   `nombre`: (String)
*   `fase_actual`: (Enum: Iniciación, Planificación, Ejecución, Monitoreo, Cierre)
*   `estado_global`: (Enum: Verde, Ámbar, Rojo)
*   `health_score`: (Float 0-100)
*   `sponsor`: (ID Referencia)
*   `pm`: (ID Referencia)
*   `fecha_inicio`: (ISO 8601 Date)
*   `fecha_fin_estimada`: (ISO 8601 Date)
*   `presupuesto_total`: (Float)

### 2. Riesgo (`riesgo_[id].md`)
*   `id_riesgo`: (String)
*   `descripcion`: (String)
*   `categoria`: (String)
*   `probabilidad`: (Float 0.0 - 1.0)
*   `impacto`: (Float 1 - 5)
*   `severidad`: (Float, Calculado: P * I)
*   `estado`: (Enum: Abierto, Mitigado, Ocurrido, Cerrado)
*   `owner`: (ID Referencia)
*   `estrategia_respuesta`: (String)
*   `plan_contingencia`: (String)
*   `fecha_identificacion`: (Date)
*   `ultimo_review`: (Date)

### 3. Issue (`issue_[id].md`)
*   `id_issue`: (String)
*   `riesgo_origen`: (Referencia, opcional)
*   `descripcion`: (String)
*   `impacto_actual`: (String)
*   `urgencia`: (Alta, Media, Baja)
*   `estado`: (Abierto, En Resolución, Bloqueado, Resuelto)
*   `owner`: (ID Referencia)
*   `resolucion_esperada`: (Date)

### 4. Decisión (`decision_[id].md`)
*   `id_decision`: (String)
*   `contexto`: (String)
*   `opciones_evaluadas`: (List)
*   `decision_tomada`: (String)
*   `justificacion`: (String)
*   `impacto_estimado`: (String)
*   `aprobador`: (ID Referencia)
*   `fecha`: (Date)

### 5. Hito / Milestone (`hito_[id].md`)
*   `id_hito`: (String)
*   `nombre`: (String)
*   `fecha_baseline`: (Date)
*   `fecha_forecast`: (Date)
*   `estado`: (Futuro, En Riesgo, Retrasado, Completado)
*   `entregables_asociados`: (List of IDs)
*   `dependencias`: (List of IDs)

### 6. Recurso / Stakeholder (`actor_[id].md`)
*   `id_actor`: (String)
*   `nombre`: (String)
*   `rol_proyecto`: (String)
*   `tipo`: (Sponsor, Equipo, Proveedor, Cliente)
*   `interes`: (Alto, Bajo)
*   `poder`: (Alto, Bajo)
*   `disponibilidad_semanal`: (Horas)

## Relaciones y Links
Los agentes cruzan información referenciando los IDs en texto o metadata (ej. `[[PRJ-001]]` o `depende_de: HITO-003`).

## Snapshots y Bitácora
- `logs/bitacora.json`: Registro append-only de eventos menores ejecutados por los agentes.
- `reports/snapshots/`: Guardado semanal automatizado del estado completo (dashboard, kpis) para análisis de tendencias.