import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# Replace this with the path to your file
file = 'output.txt'

columns = ["Time", "X", "Y", "Z"]
df = pd.read_csv(file, names=columns)

t, x, y, z = df.Time.values, df.X.values, df.Y.values, df.Z.values

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

grid = 18
line1, = ax.plot3D(x[:grid], y[:grid], z[:grid]) 


ax.set_xlim(min(x) - 1, max(x) + 1)
ax.set_ylim(min(y) - 1, max(y) + 1)
ax.set_zlim(min(z) - 1, max(z) + 1)

def update(frame, line, v):
    line.set_data_3d(x[:grid], y[:grid], v[frame*grid : (frame + 1)*grid])
    return line

num_frames = len(z) // grid

ani = FuncAnimation(fig, update, frames=num_frames, fargs=(line1, z), interval=200)

animation_file = 'visual.mp4'
ani.save(animation_file, writer='ffmpeg', fps=60)

# Show the plot
plt.show()

