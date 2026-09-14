# Flujo 05: Gestión de Riesgos (Risk Management)

**Disparador:** Nuevo riesgo detectado (por agente o manual).
**Agentes Principales:** `05_risk_agent`

## Secuencia Paso a Paso
1. **Identificación:** Se crea un borrador de `riesgo_[id].md`.
2. **Evaluación:** `05_risk_agent` asigna Probabilidad (P) e Impacto (I) basándose en contexto histórico.
3. **Estrategia:** Sugiere estrategia (Evitar, Mitigar, Transferir, Aceptar) y borrador de plan de contingencia.
4. **Revisión PM:** PM aprueba/modifica la evaluación.
5. **Monitoreo:** Cada semana, durante el Tracking, el agente re-evalúa si los triggers del riesgo se han activado, recomendando escalar a Issue si es inminente.