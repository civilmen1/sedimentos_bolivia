# Flujo 06: Gestión de Issues (Issue Management)

**Disparador:** Materialización de un riesgo o problema imprevisto bloqueante.
**Agentes Principales:** `05_risk_agent`, `01_pmo_orchestrator`

## Secuencia Paso a Paso
1. **Escalamiento:** Un riesgo se activa y se convierte en Issue, o se registra uno nuevo en el RAID log.
2. **Impacto:** `01_pmo_orchestrator` pausa los flujos dependientes afectados (ej. alerta al cronograma de un bloqueo inminente).
3. **Resolución:** El owner asignado trabaja en la solución fuera del sistema. El PM actualiza el estado.
4. **Cierre:** Al resolverse, el sistema documenta el tiempo de resolución (Issue Age) para métricas futuras.