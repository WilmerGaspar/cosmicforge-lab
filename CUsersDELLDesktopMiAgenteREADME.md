# GNoME Material Validation Dashboard

Intelligent web dashboard for validating AI-discovered materials, especially those from GNoME (Graph Networks for Materials Exploration).

## Features

- **Chemical Screening**: Validate material formulas using SMACT and OpenBabel.
- **Thermodynamic Stability**: Check stability via Materials Project API.
- **Mechanical Properties**: Predict properties using ML models (MEGNet, ALIGNN).
- **Economic Analysis**: Estimate costs based on elemental compositions.
- **Synthesizability Scoring**: Semáforo system (green/yellow/red) for potential.
- **Batch Processing**: Validate multiple materials from uploaded files.
- **AI Agent Chat**: Get recommendations and insights from integrated LLM.
- **Visualizations**: Interactive charts with Plotly.

## Installation

1. Clone the repository.
2. Install dependencies: Collecting streamlit==1.55.0 (from -r requirements.txt (line 1))
  Downloading streamlit-1.55.0-py3-none-any.whl.metadata (9.8 kB)
3. Run the app: 

## Usage

1. Enter API keys in the sidebar (Materials Project, OpenAI, Groq).
2. Use the tabs for individual validation, batch processing, or agent chat.
3. For individual validation, enter a formula like 'LiFePO4' and click Validate.
4. View results, semáforo score, and visualizations.

## Requirements

See requirements.txt for dependencies.

## Contributing

Contributions welcome. Follow the plan in docs/superpowers/plans/.
