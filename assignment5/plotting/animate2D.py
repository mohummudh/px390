import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap

# Replace this with the path to your file
file = 'output.txt'

columns = ["Time", "X", "Y", "Z"]
df = pd.read_csv(file, names=columns)

x, y, z, time = df.X.values, df.Y.values, df.Z.values, df.Time.values
grid = 18

fig, ax = plt.subplots()
squares = []
norm = plt.Normalize(z.min(), z.max())

# Custom colormap: white to black
# cmap = LinearSegmentedColormap.from_list('white_black', ['white', 'grey', 'black'], N=256)
cmap = plt.cm.viridis

# Function to create squares
def create_square(x, y, color):
    square = Rectangle((x-0.5, y-0.5), 1, 1, color=color)
    ax.add_patch(square)
    return square

# Initialize squares and timestamp
for i in range(grid):
    color = cmap(norm(z[i]))
    square = create_square(x[i], y[i], color)
    squares.append(square)

time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label='Z value')
ax.set_xlim(min(x) - 1, max(x) + 1)
ax.set_ylim(min(y) - 1, max(y) + 1)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_aspect('equal', adjustable='box')
ax.set_title('2D Grid Animation with White-Black Gradient and Timestamp')

# Update function for animation
def update(frame):
    for i, square in enumerate(squares):
        idx = frame*grid + i
        if idx < len(z):
            square.set_color(cmap(norm(z[idx])))
            square.set_xy((x[idx]-0.5, y[idx]-0.5))
    # Update timestamp
    time_text.set_text(f'Time: {time[frame*grid]:.2f}')

num_frames = len(z) // grid

ani = FuncAnimation(fig, update, frames=num_frames, interval=200)

animation_file = '2d_grid_visual.mp4'
ani.save(animation_file, writer='ffmpeg', fps=60)

plt.show()
