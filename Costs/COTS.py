import matplotlib.pyplot as plt
import numpy as np

# Categories
labels = [
    "Propulsion System Installation",
    "Avionics Installation",
    "De-icing System Installation",
    "Landing Gear Installation",
    "Interior Cabin Installation"
]

values = [35.62, 11.88, 8.91, 2.34, 6.57]
sizes = [value / sum(values) for value in values]

# Create figure
fig, ax = plt.subplots(figsize=(12, 12))
ax.pie(sizes, labels=labels, autopct='%1.f%%', textprops={'fontsize': 17})

plt.tight_layout()
plt.savefig('COTS.png')