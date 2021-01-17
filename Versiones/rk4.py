import numpy as np
import scipy.integrate as ode
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def twoBody(t, y):

    mu = 3.986004418 * 10**14

    r = np.sqrt(y[0]**2 + y[1]**2 + y[2]**2)

    ydot    = np.empty((6,1))
    ydot[:] = np.nan

    ydot[0] = y[3]             
    ydot[1] = y[4]             
    ydot[2] = y[5]             
    ydot[3] = (-mu/(r**3))*y[0]
    ydot[4] = (-mu/(r**3))*y[1]
    ydot[5] = (-mu/(r**3))*y[2]

    return ydot


# In m and m/s
# first three are the (x, y, z) position
# second three are the velocities in those same directions respectively
Y0 = np.array([-5614924.5443320004,
               -2014046.755686,
               2471050.0114869997,
               -673.03650300000004,
               582.41158099999996,
               1247.7034980000001])

#solution = ode.RK45(twoBody, t0 = 0.0, y0 = Y0, t_bound = 351.0)

solution=ode.solve_ivp(twoBody,[0., 35100.],Y0,method='LSODA', atol=1e-4, rtol=1e-6)
print(solution)
print(type(solution.y))
RR=solution.y

x1=RR[0]
y1=RR[1]
z1=RR[2]


#Imprimimos la funcion
ax=plt.axes

fig=plt.figure()
ax=Axes3D(fig)
plt.title('Posiciones de dos satelites EMFF')

ax.plot(x1,y1,z1)

plt.show()


