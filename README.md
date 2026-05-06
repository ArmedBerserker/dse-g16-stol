# dse-g16-stol

This repository contains the code, documentation, analyses, and supporting tools used by Group 16 for the 2026 Design Synthesis Exercise project on a STOL aircraft concept.

The goal of this repository is to keep all project work organized, reproducible, and easy for team members to understand and extend.

---

## Repository Overview

The repository may contain:

- Aircraft sizing and performance code
- Aerodynamic, propulsion, structural, and stability analysis tools
- Class definitions for aircraft components and configurations
- Input files and generated results
- Project documentation
- Trade studies and design notebooks
- Verification and validation material

Each subfolder should include its own short explanation if the contents are not immediately obvious.

---

## Getting Started

Clone the repository:

```bash
git clone <repository-url>
cd dse-g16-stol
```

Create and activate a virtual environment if needed:

```bash
python -m venv venv
source venv/bin/activate      # macOS/Linux
venv\Scripts\activate         # Windows
```

Install the required packages:

```bash
pip install -r requirements.txt
```

If `requirements.txt` does not yet exist, install the packages needed for the specific script or notebook you are running.

---

## Loading Aircraft Classes

To load a class from a file, use the project loader:

```python
from classes.aircraft_2 import loader

target = loader.load(file_path, target_class)
```

This approach should be preferred where possible, as it follows the existing structure of the repository.

---

## Coding Guidelines

Please follow the style already used in the repository. Consistent formatting and naming make it easier for the whole team to work on the same codebase.

### Class Names

Class names should:

- Start with a capital letter
- Use as few words as possible
- Use `PascalCase`

Examples:

```python
class Aircraft:
    pass

class Wing:
    pass

class PropulsionSystem:
    pass
```

Avoid unnecessarily long names such as:

```python
class CompleteAircraftConfigurationObject:
    pass
```

### Variable Names

Variable names should:

- Use lowercase letters
- Use underscores between words
- Be clear but not overly long

Examples:

```python
wing_area = 24.5
aspect_ratio = 8.2
fuel_mass = 120.0
```

Well-known aerodynamic variables may keep their standard compact notation.

Examples:

```python
cd0 = 0.028
clmax = 2.1
ld_ratio = 14.5
```

There is no need to write these as:

```python
c_d_0 = 0.028
c_l_max = 2.1
```

unless doing so improves clarity in a specific context.

### File Names

File names should be lowercase and use underscores where needed:

```text
aircraft_sizing.py
range_estimation.py
wing_loading.py
```

Avoid spaces, capital letters, and vague names.

Poor examples:

```text
New Code.py
FinalVersion2.py
Aircraft Stuff.py
```

---

## Branch Naming Guidelines

Use clear and consistent branch names so that everyone can understand the purpose of a branch without opening it.

Recommended branch format:

```text
<short-description>
```

Examples:

```text
range-estimation
landing-gear-sizing
propulsion-bug
update-readme
wing-loading-study
aircraft-class
```

### Branch Naming Rules

Branch names should:

- Be lowercase
- Use hyphens between words
- Be short but descriptive
- Avoid spaces
- Avoid personal names unless necessary
- Avoid vague names such as `new`, `update`, `final`, or `test2`

Good examples:

```text
hybrid-range-model
stol-performance
mass-fraction-equation
code-style-guide
```

Poor examples:

```text
MyBranch
new_stuff
final-version
testing
branch1
```

---

## Commit Guidelines

Commits should describe what changed and why. Keep each commit focused on one logical change.

Good commit messages:

```text
Add hybrid range estimation function
Fix fuel mass calculation in sizing script
Update README with branch naming rules
Refactor aircraft loading utilities
```

Poor commit messages:

```text
changes
stuff
fixed
final
work
```

Use present-tense verbs where possible, for example:

```text
Add
Fix
Update
Remove
Refactor
```

---

## Pull Request Guidelines

Before opening a pull request:

1. Make sure the code runs.
2. Check that your branch is up to date with the main branch.
3. Remove unnecessary temporary files.
4. Add comments where the logic is not obvious.
5. Include a short explanation of what was changed.

A good pull request description should include:

```text
## Summary
Briefly describe the purpose of the change.

## Changes
- Added ...
- Updated ...
- Fixed ...

## Notes
Mention assumptions, limitations, or things that still need work.
```

---

## Project Structure

A typical structure may look like:

```text
dse-g16-stol/
│
├── classes/             # Aircraft and component class definitions
├── data/                # Input data and reference values
├── docs/                # Documentation and reports
├── notebooks/           # Jupyter notebooks for analysis
├── scripts/             # Standalone scripts
├── tests/               # Verification and test scripts
├── results/             # Generated outputs and plots
├── requirements.txt     # Python package requirements
└── README.md
```

If a new folder is added, make sure its purpose is clear.

---

## Data and Results

Avoid committing large generated files unless they are required for the project record.

Prefer committing:

- Source code
- Input files
- Final plots or tables used in reports
- Small reference datasets
- Documentation

Avoid committing:

- Temporary files
- Large simulation outputs
- Cache files
- Personal copies of scripts
- Automatically generated files that can be recreated easily

---

## Collaboration Rules

To keep the repository manageable:

- Work on a separate branch for each task.
- Pull the latest changes before starting new work.
- Avoid directly committing to the main branch.
- Keep functions and files focused on one purpose.
- Document assumptions in the code or in a nearby markdown file.
- Do not delete or overwrite someone else’s work without checking first.
- Use pull requests for changes that affect shared code.

---

## Documentation Expectations

When adding a new model, script, or analysis, include enough information for another team member to understand:

- What the file does
- What inputs it needs
- What outputs it produces
- What assumptions are used
- How to run it

For example:

```python
"""
Calculates the required fuel mass for a given aircraft range.

Inputs:
    range_target: desired range [m]
    ld_ratio: lift-to-drag ratio [-]
    fuel_specific_energy: fuel specific energy [J/kg]

Outputs:
    fuel_mass: required fuel mass [kg]

Assumptions:
    Constant L/D during cruise.
    No reserve fuel included.
"""
```

---

## Units

Always make units clear. Prefer SI units unless there is a specific reason to use another system.

Recommended units:

| Quantity | Unit |
|---|---|
| Mass | kg |
| Weight | N |
| Distance | m or km |
| Time | s |
| Speed | m/s or kts| 
| Energy | J |
| Power | W or hp |
| Area | m² |
| Density | kg/m³ |

Include units in variable comments, function docstrings, or column names.

Examples:

```python
aircraft_mass = 1200.0  # kg
wing_area = 18.5        # m^2
cruise_speed = 62.0     # m/s
```
