import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

#data = np.loadtxt("states.txt",dtype=np.float32)
data = np.loadtxt("states.txt",dtype=np.float32,delimiter=',')

width = 2;

time = data[:,0];
position = data[:,1:4];
velocity = data[:,4:7];
euler_angle = data[:,7:10];
angular_rate = data[:,10:13];


fig = plt.figure()
ax = fig.add_subplot(221, projection='3d');
ax.plot(position[:,0], position[:,1], -position[:,2], linewidth = width)
ax.set_xlabel("x");
ax.set_ylabel("y");
ax.set_zlabel("z");
plt.xlim(-400, 400)
plt.ylim(-400, 1200)
ax.set_title("Position")
ax.legend(["Trajectory"]);
ax.set_box_aspect([1,2,1])

ax = fig.add_subplot(222);
ax.plot(time, velocity, linewidth = width)
ax.set_xlabel("Time (s)");
ax.set_ylabel("Velocity (m/s)");
ax.set_title("Velocity")
ax.legend(["Vx", "Vy", "Vz"])


ax = fig.add_subplot(223);
ax.plot(time, euler_angle, linewidth = width)
ax.set_xlabel("Time (s)");
ax.set_ylabel("Euler angle (deg)");
ax.set_title("Euler angle")
ax.legend(["roll", "pitch", "yaw"])


ax = fig.add_subplot(224);
ax.plot(time, angular_rate, linewidth = width)
ax.set_xlabel("Time (s)");
ax.set_ylabel("Angular Rate (deg/s)");
ax.set_title("Angular Rate")
ax.legend(["Wx", "Wy", "Wz"])

plt.show()

