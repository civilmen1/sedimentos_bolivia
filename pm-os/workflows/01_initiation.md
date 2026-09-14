# Flujo 01: Inicio del Proyecto (Initiation)

**Disparador:** Solicitud de nuevo proyecto por parte del Negocio.
**Agentes Principal:** `02_initiation_charter`, `01_pmo_orchestrator`

## Secuencia Paso a Paso
1. **Captura de Contexto:** El PM (humano) sube el Business Case inicial (o notas del Kickoff) al sistema.
2. **Generación de Borrador:** El agente `02_initiation_charter` procesa el texto, extrae objetivos, alcance preliminar y justificación. Redacta el borrador del Project Charter basado en `templates/project_charter.md`.
3. **Revisión de Calidad:** `09_quality_governance` revisa que el Charter tenga criterios SMART y métricas de éxito. Si falta, solicita corrección al agente 02.
4. **Validación Humana:** El sistema envía un `decision_request.md` al PM y Sponsor para revisar el borrador.
5. **Aprobación:** Sponsor aprueba externamente. PM ingresa "Aprobado" en el sistema.
6. **Creación de Baseline:** `12_living_documentation` guarda el Charter final en `/projects/[PRJ-ID]/charter.md` y cambia el estado del proyecto a "Planificación".