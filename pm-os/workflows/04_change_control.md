# Flujo 04: Control de Cambios (Change Control)

**Disparador:** Ingreso de solicitud de cambio por Stakeholder.
**Agentes Principales:** `03_scope_requirements`, `04_planning_schedule`, `07_finance_agent`

## Secuencia Paso a Paso
1. **Ingreso:** El PM registra la solicitud básica en el sistema.
2. **Evaluación 360:**
   - `03_scope_requirements` detalla cómo afecta la WBS.
   - `04_planning_schedule` simula el impacto en el camino crítico.
   - `07_finance_agent` calcula el coste adicional.
3. **Consolidación:** `01_pmo_orchestrator` agrupa los impactos en un borrador de `change_request.md`.
4. **Decisión:** El PM revisa y escala al CCB (Change Control Board).
5. **Aprobación:** Si se aprueba, `12_living_documentation` actualiza las Baselines y archiva el CR como aprobado.