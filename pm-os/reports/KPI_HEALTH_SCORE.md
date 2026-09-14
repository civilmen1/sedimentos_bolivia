# Sistema de KPIs y Project Health Score

Este documento define cómo Aura PM-OS mide objetivamente la salud del proyecto (evitando el "síndrome del semáforo sandía": verde por fuera, rojo por dentro).

## 1. KPIs Fundamentales (Earned Value Management)

### Schedule Performance Index (SPI)
- **Fórmula:** `EV / PV` (Earned Value / Planned Value)
- **Fuente:** `04_planning_schedule`
- **Periodicidad:** Semanal
- **Umbrales:**
  - 🟢 Verde: `SPI >= 0.95`
  - 🟡 Ámbar: `0.85 <= SPI < 0.95`
  - 🔴 Rojo: `SPI < 0.85`
- **Acción Sugerida (Si Rojo):** Disparar flujo 07 (Re-planning). Evaluar Fast-tracking.

### Cost Performance Index (CPI)
- **Fórmula:** `EV / AC` (Earned Value / Actual Cost)
- **Fuente:** `07_finance_agent`
- **Periodicidad:** Semanal/Mensual
- **Umbrales:**
  - 🟢 Verde: `CPI >= 0.95`
  - 🟡 Ámbar: `0.90 <= CPI < 0.95`
  - 🔴 Rojo: `CPI < 0.90`
- **Acción Sugerida (Si Rojo):** Disparar flujo 08 (Budget Deviation). Calcular EAC (Estimate at Completion).

## 2. KPIs de Riesgos e Issues

### Índice de Riesgo Crítico (CRI)
- **Fórmula:** Suma de severidad de Riesgos Activos con Severidad > 15
- **Fuente:** `05_risk_agent` (RAID Log)
- **Umbrales:**
  - 🟢 Verde: `0 - 15`
  - 🟡 Ámbar: `16 - 30`
  - 🔴 Rojo: `> 30`
- **Acción Sugerida:** Escalar riesgos inminentes al Sponsor.

### Issue Age (Edad Promedio de Issues Críticos)
- **Fórmula:** Promedio de días desde la creación de issues `estado: Abierto` y `Urgencia: Alta`
- **Umbrales:**
  - 🟢 Verde: `< 3 días`
  - 🟡 Ámbar: `3 - 7 días`
  - 🔴 Rojo: `> 7 días`
- **Acción Sugerida:** El Orquestador envía recordatorio automático a los Owners de los issues bloqueados.

## 3. Health Score Consolidado

Aura PM-OS consolida las métricas anteriores en un Score global (0-100) para evitar que una buena métrica oculte un desastre.

### Fórmula de Ponderación Base:
- `SPI Score (0-40)`: SPI * 40 (Máx 40)
- `CPI Score (0-40)`: CPI * 40 (Máx 40)
- `Risk Deduction`: Restar 5 puntos por cada Riesgo Crítico Abierto (Máx -20)
- `Issue Deduction`: Restar 5 puntos por cada Issue Crítico con Edad > 7 días (Máx -20)

`Health Score = SPI_Score + CPI_Score - Risk_Deduction - Issue_Deduction`

### Umbrales del Health Score:
- 🟢 **VERDE (80 - 100):** Proyecto bajo control. Decisiones informativas.
- 🟡 **ÁMBAR (60 - 79):** Desviaciones detectadas. Requiere acciones preventivas por parte del PM. El PMO Orchestrator aumenta la frecuencia de monitoreo.
- 🔴 **ROJO (< 60):** Proyecto en crisis. Requiere intervención inmediata del Sponsor o Comité de Dirección (Flujo 09).

### Regla "Veto" (Kill Switch)
Si el SPI es `< 0.70` **O** el CPI es `< 0.70`, el Health Score se fuerza a **ROJO (<60)** independientemente del resto de las métricas. (Un proyecto que cuesta el doble o tarda el doble no puede estar "Ámbar").