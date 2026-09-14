# Flujo 02: Planificación (Planning)

**Disparador:** Proyecto pasa a estado "Planificación".
**Agentes Principales:** `03_scope_requirements`, `04_planning_schedule`, `06_resources_agent`, `07_finance_agent`

## Secuencia Paso a Paso
1. **Definición de Alcance:** `03_scope_requirements` ingesta el Charter y, junto con el PM y el equipo (inputs manuales de talleres), genera la `wbs.md`.
2. **Cronograma Inicial:** `04_planning_schedule` toma la WBS, estima secuencias (apoyado por inputs del PM) y genera un borrador del cronograma (`schedule.json`).
3. **Asignación de Recursos:** `06_resources_agent` evalúa el cronograma contra la capacidad del pool de recursos. Alerta si hay roles sobrecargados.
4. **Presupuestación:** `07_finance_agent` toma las estimaciones de tiempo y perfiles (coste hora), más gastos fijos, y elabora el `budget.json` preliminar.
5. **Identificación de Riesgos Base:** `05_risk_agent` escanea el plan y sugiere riesgos estándar según la industria/metodología, creando los primeros registros en el RAID log.
6. **Aprobación de Líneas Base:** El PM solicita aprobación formal. Una vez aprobadas, las líneas base se congelan.