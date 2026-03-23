# GNoME Material Validation Dashboard Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Transform cosmicforge-lab into a Python/Streamlit dashboard for validating GNoME materials with AI agent integration.

**Architecture:** Web app with validation pipeline, ML predictions, API integrations, and agent chat. Uses modular Python scripts for components.

**Tech Stack:** Python, Streamlit, pymatgen, ML libraries, SQLite, APIs.

---

## Chunk 1: Project Setup and Migration

### Task 1: Create Python Environment

**Files:**
- Create: `requirements.txt`
- Create: `.gitignore` (Python)

- [ ] **Step 1: Write requirements.txt**

```
streamlit
pymatgen
matminer
ase
scikit-learn
mendeleev
openbabel
cgcnn
mp-api
plotly
pandas
numpy
requests
openai
groq-sdk
geneticalgorithm
tensorflow
torch
sqlalchemy
streamlit-authenticator
```

- [ ] **Step 2: Update .gitignore for Python**

```
__pycache__/
*.pyc
*.pyo
*.pyd
.Python
env/
venv/
.venv/
pip-log.txt
pip-delete-this-directory.txt
.tox/
.coverage
.coverage.*
.cache
nosetests.xml
coverage.xml
*.cover
*.log
.git
.mypy_cache
.pytest_cache
.hypothesis
```

- [ ] **Step 3: Commit setup**

```bash
git add requirements.txt .gitignore
git commit -m "feat: setup Python environment for dashboard"
```

### Task 2: Migrate to Streamlit App

**Files:**
- Create: `app.py` (main dashboard)
- Modify: `README.md`

- [ ] **Step 1: Copy provided app.py code**

Paste the entire GNoME validation code into `app.py`.

- [ ] **Step 2: Update README.md**

Replace content with:

```
# GNoME Material Validation Dashboard

Intelligent web dashboard for validating AI-discovered materials.

## Installation

pip install -r requirements.txt
streamlit run app.py
```

- [ ] **Step 3: Commit migration**

```bash
git add app.py README.md
git commit -m "feat: migrate to Streamlit dashboard with basic validation"
```

## Chunk 2: API Integrations

### Task 3: Add Materials Project API

**Files:**
- Modify: `app.py` (add MP integration)

- [ ] **Step 1: Add MP API key input**

In sidebar:

```python
mp_api_key = st.sidebar.text_input("Materials Project API Key", type="password")
```

- [ ] **Step 2: Update thermodynamic_stability function**

Use mp-api instead of simulation.

- [ ] **Step 3: Test MP connection**

Run app and validate a material with API key.

- [ ] **Step 4: Commit API integration**

```bash
git add app.py
git commit -m "feat: integrate Materials Project API for real stability data"
```

### Task 4: Add Other APIs (AFLOW, etc.)

**Files:**
- Modify: `app.py`

- [ ] **Step 1: Add API clients**

Import and initialize aflow-api, oqmd-api, etc.

- [ ] **Step 2: Extend validation pipeline**

Add API calls in thermodynamic_stability.

- [ ] **Step 3: Handle fallbacks**

If API fails, use simulation.

- [ ] **Step 4: Commit extended APIs**

```bash
git add app.py
git commit -m "feat: add AFLOW, OQMD, NOMAD API integrations"
```

## Chunk 3: ML Models and Agent

### Task 5: Integrate ML Models

**Files:**
- Create: `models.py`
- Modify: `app.py`

- [ ] **Step 1: Create models.py**

```python
from megnet.models import MEGNetModel
from alignn.models.alignn import ALIGNN

# Load pre-trained models
megnet = MEGNetModel.from_file('path/to/megnet.h5')
alignn = ALIGNN()
```

- [ ] **Step 2: Add to mechanical_properties**

Use models for predictions.

- [ ] **Step 3: Test predictions**

Validate output.

- [ ] **Step 4: Commit ML integration**

```bash
git add models.py app.py
git commit -m "feat: integrate MEGNet and ALIGNN for property predictions"
```

### Task 6: Add AI Agent

**Files:**
- Create: `agent.py`
- Modify: `app.py` (add chat tab)

- [ ] **Step 1: Create agent.py**

```python
import openai

def get_agent_response(prompt):
    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[{"role": "user", "content": prompt}]
    )
    return response.choices[0].message.content
```

- [ ] **Step 2: Add chat tab in app.py**

Use agent for recommendations.

- [ ] **Step 3: Integrate with validations**

Agent can query results.

- [ ] **Step 4: Commit agent**

```bash
git add agent.py app.py
git commit -m "feat: add AI agent for intelligent recommendations"
```

## Chunk 4: Database and Advanced Features

### Task 7: Add Database

**Files:**
- Create: `database.py`
- Modify: `app.py`

- [ ] **Step 1: Create database.py**

```python
import sqlite3

def save_validation(result):
    conn = sqlite3.connect('validations.db')
    # Save logic
```

- [ ] **Step 2: Integrate in app.py**

Save results after validation.

- [ ] **Step 3: Add query features**

View history.

- [ ] **Step 4: Commit DB**

```bash
git add database.py app.py
git commit -m "feat: add SQLite database for persistence"
```

### Task 8: Add Optimization and Feedback

**Files:**
- Create: `optimization.py`
- Modify: `app.py`

- [ ] **Step 1: Create optimization.py**

Genetic algorithm for composition optimization.

- [ ] **Step 2: Add feedback loop**

Learn from validations.

- [ ] **Step 3: Integrate in dashboard**

Tabs for optimization.

- [ ] **Step 4: Commit advanced features**

```bash
git add optimization.py app.py
git commit -m "feat: add genetic optimization and feedback loop"
```

## Chunk 5: Deployment and Security

### Task 9: Add Security

**Files:**
- Modify: `app.py` (add auth)

- [ ] **Step 1: Add streamlit-authenticator**

Require login.

- [ ] **Step 2: Secure API keys**

Environment variables.

- [ ] **Step 3: Test security**

Access control.

- [ ] **Step 4: Commit security**

```bash
git add app.py
git commit -m "feat: add authentication and secure API handling"
```

### Task 10: Deployment Setup

**Files:**
- Create: `Dockerfile`
- Create: `docker-compose.yml`

- [ ] **Step 1: Create Dockerfile**

```dockerfile
FROM python:3.9
COPY . /app
WORKDIR /app
RUN pip install -r requirements.txt
CMD streamlit run app.py
```

- [ ] **Step 2: Create docker-compose**

For easy deployment.

- [ ] **Step 3: Test deployment**

Run locally with Docker.

- [ ] **Step 4: Commit deployment**

```bash
git add Dockerfile docker-compose.yml
git commit -m "feat: add Docker deployment for remote access"
```