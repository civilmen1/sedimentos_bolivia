# Flujo 08: Desviación Presupuestaria (Budget Deviation)

**Disparador:** CPI (Cost Performance Index) < 0.90 o Factura excede previsión.
**Agentes Principales:** `07_finance_agent`, `01_pmo_orchestrator`

## Secuencia Paso a Paso
1. **Detección:** Se ingresa un gasto real que rompe el umbral.
2. **Alerta Ejecutiva:** `07_finance_agent` emite una Alerta Financiera indicando el Overrun proyectado al final del proyecto (EAC).
3. **Mitigación:** Orquestador sugiere congelar contratación o usar Reserva de Contingencia.
4. **Aprobación Sponsor:** Se requiere autorización (vía CR) para liberar fondos de reserva.