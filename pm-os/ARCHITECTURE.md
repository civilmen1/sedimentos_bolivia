# Estructura del Repositorio y Arquitectura (ARCHITECTURE)

Aura PM-OS organiza su conocimiento e inteligencia en una estructura de carpetas predecible, diseñada para ser leída/escrita tanto por humanos como por LLMs.

## Árbol de Directorios (`/pm-os`)

```text
/pm-os
├── /agents/           # Lógica y scripts específicos de los agentes (si aplica Python/JS extra).
├── /dashboards/       # Vistas agregadas (UI/Markdown) consolidadas para humanos.
├── /governance/       # Reglas de negocio. (Source of truth para comportamiento del sistema).
│   ├── GOVERNANCE.md
│   └── EVENTS_AND_TRIGGERS.md
├── /knowledge/        # Modelos de datos y base de conocimiento del proyecto.
│   └── DATA_MODEL.md
├── /logs/             # Bitácoras, registros de orquestación (append-only por el sistema).
├── /projects/         # CARPETA PRINCIPAL DE DATOS VIVOS.
│   └── /PRJ-001/      # Cada proyecto tiene su propia subcarpeta.
│       ├── charter.md
│       ├── wbs.md
│       ├── schedule.json
│       ├── budget.json
│       ├── /risks/    # Documentos de riesgos individuales.
│       ├── /issues/
│       └── /decisions/
├── /prompts/          # (Source of truth de la IA). Prompts internos estandarizados por agente.
├── /reports/          # Reportes estáticos generados periódicamente (Snapshot/Read-only).
│   └── KPI_HEALTH_SCORE.md
├── /scripts/          # Automatizaciones auxiliares (CI/CD, webhooks).
├── /templates/        # Plantillas base (Source of truth metodológico).
│   ├── /outputs/      # Plantillas de salida de los agentes (Alertas, resúmenes).
│   └── [artefactos]
├── /workflows/        # Definición paso a paso de los procesos operativos.
├── README.md          # Visión del sistema.
├── AGENT_CATALOG.md   # Definición de la red.
└── RUNBOOK.md         # Operaciones diarias/semanales.
```

## Mecanismos de Actualización

1.  **Archivos Globales (`/prompts`, `/governance`, `/workflows`, `/templates`):**
    *   *Mantenidos por:* PMO / Ingeniero de Sistema.
    *   *Frecuencia:* Baja (Cambios metodológicos).
    *   *Permisos AI:* Solo lectura. Los agentes leen sus instrucciones de aquí.

2.  **Archivos por Proyecto (`/projects/...`):**
    *   *Mantenidos por:* Agentes Especializados (Escritura activa) y Humanos (Aprobaciones/Ediciones).
    *   *Frecuencia:* Alta (Diaria/Semanal).
    *   *Source of Truth:* Representa el estado real del proyecto.

3.  **Archivos de Salida (`/reports`, `/dashboards`, `/logs`):**
    *   *Mantenidos por:* Agente de Monitoreo (`10_monitoring_control`) y Orquestador.
    *   *Frecuencia:* Programada o Event-driven.
    *   *Permisos Humanos:* Solo lectura (No se editan a mano, se regeneran a partir de los datos en `/projects/`).