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
