# dse-g16-stol
This is the GitHub repository for G16 of the 2026 DSE session. It contains code and documentation utilised during the project.

## Guide for users
Please follow the trend established in the rest of the code. As a general guideline, classes should start with a capital letter, and have as few words as possible. Use all small letters for variables, with underscores between words, unless the variable is a well known aerodynamic property (for example no need to write cd0 as c_d_0). Finally name functions with camelCase (ie. a function is name with the first letter lowercase, with new words starting with uppercase)

### General Things
To implement a class from a file. The easiest way to do this is to do:
```
from classes.aircraft_2 import loader
target = loader.load(file_path, target_class)
```

Aircraft(Req)