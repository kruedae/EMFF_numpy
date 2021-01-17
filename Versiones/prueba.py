import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

ts=np.linspace(0,2,100)
x=np.cos(ts)
y=np.sin(ts)
print(x)

a=5+ \
        6
print(a)
fig=plt.figure()
ax=Axes3D(fig)

ax.plot(x,y,ts)

plt.show()
