import numpy as np 
import matplotlib.pyplot as plt

# Read data from the file
with open('output.txt', 'r') as file:
    lines = file.readlines()

x_values = []
y_values = []

for line in lines:
    # Skip empty lines
    if not line.strip():
        continue

    parts = line.split(',')  # Splitting based on comma

    # Check if there are at least three values in the line
    if len(parts) >= 3:
        x = float(parts[1])  # Second column
        y = float(parts[2])  # Third column
        x_values.append(x)
        y_values.append(y)

# Create the plot
plt.plot(x_values, y_values)
plt.xlabel('Second Column')
plt.ylabel('Third Column')
plt.title('Plot of Second Column vs. Third Column')
plt.grid(True)

# Show the plot (or save it as an image file)
plt.show()
