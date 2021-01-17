#IMPORTO LIBRERIAS
from astropy.constants import G, M_earth, R_earth
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Defino constantes astronomicas basicas
mu = G.value*M_earth.value
Re = R_earth.value

#Defino la funcion que convierte elementos orbitales en vectores cartesianos
def kep_2_cart(a,e,i,omega_AP,omega_LAN,T, EA):

    #1
    n = np.sqrt(mu/(a**3))
    M = n*(t - T)
    #2
#    MA = EA - e*np.sin(EA)
    MA=EA
    #3
    nu = 2*np.arctan(np.sqrt((1+e)/(1-e)) * np.tan(EA/2))
    #4
    r = a*(1 - e*np.cos(EA))
    #5
    h = np.sqrt(mu*a * (1 - e**2))
    #6
    Om = omega_LAN
    w =  omega_AP
    
    i=i*np.pi/180
    Om=Om*np.pi/180
    w=w*np.pi/180

    X = r*(np.cos(Om)*np.cos(w+nu) - np.sin(Om)*np.sin(w+nu)*np.cos(i))
    Y = r*(np.sin(Om)*np.cos(w+nu) + np.cos(Om)*np.sin(w+nu)*np.cos(i))
    Z = r*(np.sin(i)*np.sin(w+nu))

    #7
    p = a*(1-e**2)

    V_X = (X*h*e/(r*p))*np.sin(nu) - (h/r)*(np.cos(Om)*np.sin(w+nu) + \
    np.sin(Om)*np.cos(w+nu)*np.cos(i))
    V_Y = (Y*h*e/(r*p))*np.sin(nu) - (h/r)*(np.sin(Om)*np.sin(w+nu) - \
    np.cos(Om)*np.cos(w+nu)*np.cos(i))
    V_Z = (Z*h*e/(r*p))*np.sin(nu) - (h/r)*(np.cos(w+nu)*np.sin(i))

    return [X,Y,Z],[V_X,V_Y,V_Z]


#ELEMENTOS ORBITALES(ACA ES DONDE SE CAMBIAN LOS PARAMETROS SEGUN TLE CORRESPONDIENTE)

#a,e,i,omega_AP,omega_LAN,T,EA = 6795e3,0.0003099,51.6357,256.7529,198.7788,     0,103.3278#ESTACION ESPACIAL

#a,e,i,omega_AP,omega_LAN,T,EA=6795e3,0.0019298,0,111.0299,144.508,0,103.3278#FAC-SAT
a,e,i,omega_AP,omega_LAN,T,EA=6795e3,0.7,45,0,0,0,103.3278



t=np.linspace(0,2*np.pi,200)

X=[0]*len(t)
Y=[0]*len(t)
Z=[0]*len(t)
c=0
for j in t:
    ru_test2, v_test2 = kep_2_cart(a,e,i,omega_AP,omega_LAN,T, j)
    X[c] = ru_test2[0]*(1e-3)
    Y[c] = ru_test2[1]*(1e-3)
    Z[c] = ru_test2[2]*(1e-3)
    c+=1

#GRAFICA3-D
ax=plt.axes

fig=plt.figure()
ax=Axes3D(fig)
plt.title('Elipse')

ax.set_xlabel('X[km]')
ax.set_ylabel('Y[km]')
ax.set_zlabel('Z[km]')


ax.plot(X,Y,Z,color="orange")

u = np.linspace(0, 2 * np.pi, 200)
v = np.linspace(0, np.pi, 200)
x = 6371 * np.outer(np.cos(u), np.sin(v))
y = 6371 * np.outer(np.sin(u), np.sin(v))
z = 6371 * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b',alpha=0.5)




plt.show()


