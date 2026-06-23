import matplotlib.pyplot as plt
import numpy as np

# Categories
labels = [
    "Advanced R&D",
    "Wing Design",
    "Empennage Design",
    "Fuselage Design",
    "Landing Gear Design",
    "Auxiliary Systems Design",
    "Final Aircraft Design",
    "Software & Simulator",
    "Development Support",
    "Flight Test"
]

values = [3.09, 17.94, 5.78, 8.95, 2.61, 3.92, 4.36, 2.41, 10.16, 2.63]
sizes = [value / sum(values) for value in values]

# Create figure
fig, ax = plt.subplots(figsize=(10, 10))
ax.pie(sizes, labels=labels, autopct='%1.f%%', textprops={'fontsize': 16})

plt.tight_layout()
plt.savefig('RDTE.png')