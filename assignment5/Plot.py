import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
import matplotlib.cm as cm
import imageio

delay = 100


# Open the input file
A = 1
try:
    with open('input.txt', 'r') as file:
        Nx = int(file.readline().strip())
        Ny = int(file.readline().strip())
        for _ in range(4):  # Skip the first 5 lines
            file.readline()
        A = float(file.readline())  # Read and store the value from the seventh line
        A = math.sqrt(A)
except FileNotFoundError:
    print("File not found.")
except Exception as e:
    print("Error:", e)

print(Nx,Ny,A)
# Step 1: Read the data from the file
data = np.genfromtxt('output.txt', delimiter=',', dtype=float)

# Step 2: Organize data into time steps
time_steps = {}
for row in data:
    t = row[0]
    x = int(float(row[1]))  # Convert to float and then to integer
    y = int(float(row[2]))  # Convert to float and then to integer
    value = row[3]
    if t not in time_steps:
        time_steps[t] = np.full((Nx, Ny), np.nan)  # Initialize with NaNs
    time_steps[t][x, y] = value

# Step 3: Fill missing values with NaN
# for t in time_steps:
#     time_steps[t] = np.nan_to_num(time_steps[t])

# Step 4: Create a single window to display the timelapse
fig, ax = plt.subplots(figsize=(5, 5)) 

cmap = cm.seismic  # Choose your desired colormap
cmap.set_bad(color='black')  # Set the color for NaN values


# Step 5: Plot the initial temperature map
im = ax.imshow(time_steps[0], cmap=cmap, vmin=-A, vmax=A, extent=[0, Ny, 0, Nx])
plt.colorbar(im, ax=ax)
ax.set_title('')
ax.set_xlabel('X')
ax.set_ylabel('Y')

# Step 6: Update function for animation
def update(frame):
    global current_timestep
    current_timestep = frame
    im.set_array(time_steps[frame])
    ax.set_title(f'Time Step {current_timestep}')
    return im,

# Step 7: Create the animation
animation = FuncAnimation(fig, update, frames=sorted(time_steps.keys()), interval=delay, blit=True, repeat=True)

gif_file = 'animation.gif'
with imageio.get_writer(gif_file, mode='I') as writer:
    for frame in sorted(time_steps.keys()):
        # Replace NaN values with a predefined fill value (e.g., 0)
        temp_frame = np.nan_to_num(time_steps[frame], nan=A+1)
        # Convert to uint8 for imageio
        writer.append_data((temp_frame * 255).clip(0, 255).astype(np.uint8))  # Scale and convert to uint8

plt.show()