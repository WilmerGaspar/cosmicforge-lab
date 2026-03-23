---
description: Superpowers Methodology (TDD, Planning, Systematic Debugging)
---

# Superpowers Methodology

This workflow enforces the Superpowers methodology for all development tasks in this project.

## 1. Brainstorming & Planning (Before Coding)
- **Do not jump right into writing code**. 
- Ask the user what they are really trying to achieve. 
- Refine the design interactively.
- Once the design is approved, write an `implementation_plan.md`. The plan must be detailed enough for a junior engineer to follow.
- Split work into tiny tasks (2 to 5 minutes each).

## 2. Test-Driven Development (TDD)
- **Red-Green-Refactor**: Always write a failing test first.
- Observe the test fail (RED).
- Write the minimum code necessary to make it pass (GREEN).
- Commit / Refactor.
- Focus on YAGNI (You Aren't Gonna Need It) and DRY (Don't Repeat Yourself).

## 3. Sandboxed Execution
- Use isolated spaces (like git worktrees or separate branches) when exploring complex changes.
- Ensure a clean test baseline before starting.

## 4. Systematic Debugging
- Do not guess the error.
- Use a 4-phase root cause process.
- Verify before completing the fix.
- Add defense in depth.

## 5. Review and Completion
- Perform code review against the plan.
- Verify tests pass before completing a task or merging a branch.
