import matplotlib.pyplot as plt
import numpy as np

# Categories
labels = [
    "Engineering Design",
    "Tooling",
    "Manufacturing",
    "Quality Control",
    "Development Support",
    "Flight Test",
    "Manufacturing Materials"
]

values = [49.07, 18.52, 39.87, 5.84, 10.16, 2.63, 4.37]
sizes = [value / sum(values) for value in values]

# Create figure
fig, ax = plt.subplots(figsize=(10, 10))
ax.pie(sizes, labels=labels, autopct='%1.f%%', textprops={'fontsize': 16})

plt.tight_layout()
plt.savefig('RDTEProd.png')