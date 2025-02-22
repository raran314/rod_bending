import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation


N=30
#node_positions = [np.array([np.sin(2 * np.pi * i / N), np.cos(2 * np.pi * i / N), 2]) for i in range(N)]
node_positions_0 = [np.array([2 * np.sin(4 * np.pi * i / N), i, 2]) for i in range(N)] #initial configuration of rod 


def edge_length(i,node_positions):
    ##length of the edge between nodes i+1 and i##
    if i >= len(node_positions) - 1:  # Prevent index out of bounds
        return 0  # Return zero if out of bounds
    else:
        return np.linalg.norm(np.array(node_positions[i+1])-np.array(node_positions[i]))

def delta_s(i,node_positions):
    return (1/2)*(edge_length(i,node_positions)+edge_length(i-1,node_positions))
    
def xi(i):
    if i <= 0 or i >= len(node_positions) - 1:
        return 0  
    delta_n = np.array(node_positions[i])-np.array(node_positions[i-1])
    delta_p = np.array(node_positions[i+1])-np.array(node_positions[i])
    return (np.dot(delta_n, delta_p))/((edge_length(i-1,node_positions)**2)*(edge_length(i,node_positions)**2))

def edge_angle(i):
    ##angle between vectors defining the edges i+1 and i, calculated with node i as the starting point##
    return np.arccos(xi(i)*edge_length(i-1,node_positions)*edge_length(i,node_positions))

def F_stretch(alpha_s,node_positions,node_positions_0,j):
    ##Bending force, F = -\partial E_stretch / \partial \vec{x}_{j}##
    force = -2*alpha_s*(((edge_length(j-1,node_positions)-edge_length(j-1,node_positions_0))/(edge_length(j-1,node_positions)))*(node_positions[j]-node_positions[j-1])-((edge_length(j,node_positions)-edge_length(j,node_positions_0))/(edge_length(j,node_positions)))*(node_positions[j+1]-node_positions[j]))
    return force

def F_bend(alpha_b,node_positions,node_positions_0,j):
    def A(d):
        return (2*edge_angle(d))/(-np.sin(edge_angle(d))*delta_s(d,node_positions))
    def B(a,b,c):
        return (xi(a)/(edge_length(b,node_positions)*(edge_length(c,node_positions)**3)))
    def C(d,K,Q):
        return A(d)*K-((edge_angle(d)**2)/(delta_s(d,node_positions))**2)*Q

    if j + 1 >= len(node_positions) or j + 2 >= len(node_positions):
        return np.zeros_like(node_positions[j])
    
    P1 = (1/(edge_length(j-1,node_positions)*edge_length(j,node_positions)))*(node_positions[j+1]-node_positions[j]+node_positions[j-1])-(B(j,j,j-1)*(node_positions[j]-node_positions[j-1])-B(j,j-1,j)*(node_positions[j+1]-node_positions[j]))
    P2 = (1/2)*((1/edge_length(j-1,node_positions))*(node_positions[j]-node_positions[j-1])-(1/edge_length(j,node_positions))*(node_positions[j+1]-node_positions[j]))
    P3 = (1/(edge_length(j-2,node_positions)*(edge_length(j-1,node_positions)))-B(j-1,j-1,j-2))*(node_positions[j-1]-node_positions[j-2])+B(j-1,j-2,j-1)*(node_positions[j]-node_positions[j-1])
    P4 = (1/2)*(1/(edge_length(j-1,node_positions)))*(node_positions[j]-node_positions[j-1])
    P5 = (1/(edge_length(j,node_positions)*(edge_length(j+1,node_positions)))-B(j+1,j,j+1))*(node_positions[j+1]-node_positions[j+2])+B(j+1,j+1,j)*(node_positions[j+1]-node_positions[j])
    P6 = -(1/2)*(1/(edge_length(j,node_positions)))*(node_positions[j+1]-node_positions[j])
          
    return -alpha_b*(C(j,P1,P2)+C(j-1,P3,P4)+C(j+1,P5,P6))

def iteration(node_positions,time_step, gamma, alpha_s, alpha_b)
    
#def iteration(node_positions, time_step, gamma, alpha_s, alpha_b, node_positions_0):
#    node_positions_new = []
#    for j in range(1, len(node_positions) - 1):
#        node_positions_new.append(node_positions[j] + (time_step/gamma) *
#                                 (F_bend(alpha_b, node_positions, node_positions_0, j) +
#                                  F_stretch(alpha_s, node_positions, node_positions_0, j)))
#    node_positions_new.insert(0,np.array([1,1,1]))
#    node_positions_new.append(np.array([2,2,2]))
#    return node_positions_new

#def apply_iterations(node_positions, n, time_step, gamma, alpha_s, alpha_b, node_positions_0):
#    history = [node_positions.copy()]
#    for _ in range(n):
#        node_positions = iteration(node_positions, time_step, gamma, alpha_s, alpha_b, node_positions_0)
#        history.append(node_positions.copy())
#    return history


# Example parameters
time_step = 0.001
gamma = 1.0
alpha_s = 100
alpha_b = 1
n_iterations = 5

# Apply iterations
history = apply_iterations(node_positions_0, n_iterations, time_step, gamma, alpha_s, alpha_b, node_positions_0)

# Create animation
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')
ax.set_title('3D Evolution of Points Over Time')
sc = ax.scatter([], [], [], c='b', marker='o')
lines = []
for _ in range(len(node_positions_0) - 1):
    lines.append(ax.plot([], [], [], c='gray')[0])

def update(frame):
    current_points = history[frame]
    x, y, z = np.array(current_points).T
    sc._offsets3d = (x, y, z)
    for i in range(len(current_points) - 1):
        lines[i].set_data([x[i], x[i+1]], [y[i], y[i+1]])
        lines[i].set_3d_properties([z[i], z[i+1]])
    return sc, *lines

ani = animation.FuncAnimation(fig, update, frames=len(history), interval=200, blit=False)
plt.show()
