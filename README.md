# dse-g16-stol

This repository contains the code, documentation, analyses, and supporting tools used by Group 16 for the 2026 Design Synthesis Exercise project on a STOL aircraft concept.

The goal of this repository is to keep all project work organized, reproducible, and easy for team members to understand and extend.

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
---

## Loading Aircraft Classes

To load a class from a file, use the project loader:

```python
from classes.aircraft_2 import loader

target = loader.load(file_path, target_class)
```

This approach works on any class, as long as the .yaml files comply with the structure. Possible classes are
```python
Aircraft, Wing, Engine, Fuselage, Requirements, Mission
```
---

## General Example

The idea of the pre-existing code is that `.yaml` files are treated as a ground truth and are not changed by code. **Please do no overwrite a file once iterations start**. While building code, you may add any parameters you need to any file, but try to fit them properly into the existing category. Think: Is this a requirement or a parameter? for example. Any time a top-level heading is addded to a yaml, you must update the corresponding class in `aircraft_2.py` to contain that field. A yaml can be formatted in many ways. An example is shown below:

```yaml
field1 : thing1
field2 : thing2
field3 : 2
list1 : [0.8,    0.8,    0.8]
field4 :
    subfield1 : 0.1
    subfield2 : 0.2
```

The respective class then needs fields as:

_Note that subfields are interpreted as dictionaries under the corresponding field heading_


```python
field1 : str
field2 : str
field3 : int
list1 : list
field4 : dict
```


## Example: printing the L/D of a wing

### Method 1: Load an aircraft and then access its wing
```python
from classes.aircraft_2 import loader, Aircraft

aircraft = loader.load('aircraft.yaml', Aircraft)

print(aircraft.wing.ld)
```
The benefit of this method is that you can access other parameters of the aircraft if required, for example to get the propulsive efficiency you can run `ac.engine.eta_prop`. The downside however is that you need to have a defined `aircraft.yaml` file and corresponding input files. 

### Method 2: Load the wing only 
```python
from classes.aircraft_2 import loader, Wing

wing = loader.load('wing.yaml', Wing)

print(wing.ld)
```


The benefit of this method is that you do not need the whole aircraft defined, just the wing. This may be easier when you are working on a particular system. On the other hand, this does mean other aircraft values cannot be accessed (unless you define a seperate instance of that class as needed).

### Special case: 
In some classes, you have subfields (see that yaml example before). To access the values in that case, you must run `object.field1['subfield1']` to get the same values.
## Formatting Guidelines

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
<type>/<short-description>
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
feature/hybrid-range-model
analysis/stol-performance
fix/mass-fraction-equation
docs/code-style-guide
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

*Only commit generated images if they do no change (for eg. regression lines etc)*

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
| Speed | m/s |
| Energy | J |
| Power | W |
| Area | m² |
| Density | kg/m³ |

Include units in variable comments, function docstrings, or column names.

Examples:

```python
aircraft_mass = 1200.0  # kg
wing_area = 18.5        # m^2
cruise_speed = 62.0     # m/s
```

---

