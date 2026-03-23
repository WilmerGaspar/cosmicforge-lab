---
description: GPT Researcher Methodology (Deep Research, Multi-Agent Collaboration, Objective Reporting)
---

# GPT Researcher Methodology

This workflow enforces the GPT Researcher methodology for all complex research, fact-checking, and data-gathering tasks.

## 1. Multi-Agent Architecture Simulation
When given a complex research task, adopt a multi-agent mindset:
- **Planner Agent Role:** First, analyze the research prompt. Generate a comprehensive set of sub-questions that must be answered to form a complete, objective view of the topic.
- **Execution Agent Role:** For each sub-question, perform targeted searches (web or local files). Do not rely on a single source.
- **Editor Agent Role:** Aggregate the findings from the execution phase into a cohesive, well-structured final report.

## 2. Objective and Unbiased Fact-Gathering
- Do not rely on assumptions or a single point of failure (e.g., just one website or document).
- Gather information from multiple sources to cross-reference facts and minimize bias or hallucinations.
- When conflicting information is found, present all viewpoints objectively.

## 3. Deep Research (Tree-based Exploration)
- Apply a recursive "Tree Exploration" pattern.
- Start with broad topics, and as you discover new subtopics, drill down into them (depth) while still covering the main branches (breadth).
- Maintain intelligent context management so that deep dives don't lose track of the original goal.

## 4. Comprehensive & Cited Reporting
- The final output must be a detailed, extensive research report.
- Include citations and references to the sources used (URLs for web searches, file paths for local documents).
- Provide summaries of the sources before jumping to conclusions.

## 5. Local Document Integration
- Treat the user's local workspace as a primary source of truth when applicable.
- Apply the same multi-agent querying and cross-referencing to local files (codebase, PDFs, Markdown docs, etc.) to generate internal reports.
