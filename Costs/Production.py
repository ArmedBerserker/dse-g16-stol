import matplotlib.pyplot as plt
import numpy as np

# Categories
labels = [
    "Tooling Fabrication",
    "CNC Machines",
    "Manufacturing Jigs",
    "Wing Production",
    "Empennage Production",
    "Fuselage Production",
    "Landing Gear Production",
    "Auxiliary Systems Integration",
    "Final Aircraft Assembly"
]

values = [8.96, 0.6, 8.96, 20.62, 6.65, 10.29, 3.00, 4.51, 5.01]
sizes = [value / sum(values) for value in values]

# Create figure
fig, ax = plt.subplots(figsize=(12, 12))
ax.pie(sizes, labels=labels, autopct='%1.f%%', textprops={'fontsize': 17})

plt.tight_layout()
plt.savefig('Production.png')