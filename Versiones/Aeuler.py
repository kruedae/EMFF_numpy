#Importo las librerias
import numpy as np
import scipy.integrate as ode
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#Defino el vector con los momentos de inercia principales

I=np.array([1,1,1])

mu=1
A=1
#Defino el sistema de ecuaciones a solucionar funcion a solucionar tras haberlo separado en ecuaciones de primer orden

def ec2(t,RR,I):
    phi,the,psi, dphi,dthe,dpsi = RR[0], RR[1],RR[2],  RR[3], RR[4], RR[5], RR[6]
    
    #Defino el momento angular de la rueda y su derivada
    h=np.array([0,0,0])
    dh=np.array([0,0,0])
    
    A1 = np.array( [ [ (I[2]-I[1])/I[0],0,0], [0,(I[2]-I[0])/I[1],0], [0,0,0] ] )

    A2 = np.array( [ [0,(I[1]-I[2])/I[0],(I[1]-I[2])/I[0]], [(I[2]-I[0])/I[1], 0, (I[2]-I[0])/I[1]], [(I[0]-I[1])/I[2], (I[0]-I[1])/I[2], 0] ] )

    A3 = np.array( [ [0,1/I[0],-1/I[0]], [-1/I[1],0, 1/I[1] ], [1/I[2],-1/I[2],0] ] )

    A4 = np.array( [ [-1/I[0],0,0], [0,-1/I[1],0], [0,0,-1/I[2]] ] )

    dv = [dphi, dthe, dpsi]
    v = [phi, the, psi]


    #escribimos las fuerzas magneticas de cada satelite
    ddv = np.matmul(A1,v)+np.matmul(A2,dv)+np.matmul(A3,h)+np.matmul(A4,dh)

    


    return [dv[0],dv[1],dv[2],ddv[0],ddv[1],ddv[2]]




RR0=[0,np.pi,0,0,0,0]#vector de condiciones iniciales

#Solucionamos con las condiciones anteriores usando un Runge-Kutta de orden 3 (el metodo puede variar)
RR=ode.solve_ivp(ec2,[0., 100.],RR0,method='RK45', atol=1e-4, rtol=1e-6,max_step=2,args=[I])


#Guardamos las soluciones
phi=RR[0]
theta=RR[1]
psi=RR[2]


#Rotacion

x1=np.array([0.]*len(phi))
y1=np.array([1.]*len(phi))
z1=np.array([0.]*len(phi))

for i in range():
    MR=np.array([[np.cos(the[i])*np.cos(psi[i]), np.cos(the[i])*np.sin(psi[i]), np.sin(the)],
            [np.sin(the[i])*np.sin(phi[i])*np.cos(psi[i])-np.cos(phi[i])*np.sin(psi[i]),
                np.sin(the[i])*np.sin(phi[i])*np.sin(psi[i])+np.cos(phi[i])*np.cos(psi[i]), np.cos(the[i])*np.sin(phi[i])],
            [[np.sin(the[i])*np.cos(phi[i])*np.cos(psi[i])+np.sin(phi[i])*np.sin(psi[i]),
                np.sin(the[i])*np.cos(phi[i])*np.sin(psi[i])-np.sin(phi[i])*np.cos(psi[i]), np.cos(the[i])*np.cos(phi[i])]]
            ])
    v1=np.array([x1[i],y1[i],z1[i]])
    v1=np.matmul(MR,v1)
    x1[i],y1[i],z1[i]=v1[0],v1[1],v1[2]




#Imprimimos la funcion

#3-D
#Marco de la tierra
ax=plt.axes

fig=plt.figure()
ax=Axes3D(fig)
plt.title('Trayectoria relativa del segundo satelite')
ax.scatter(xr[0],yr[0],zr[0],marker='o', color="blue")
ax.scatter(xr[len(xr)-1],yr[len(xr)-1],zr[len(xr)-1],marker='o', color="blue")

ax.set_xlabel('X[m]')
ax.set_ylabel('Y[m]')
ax.set_zlabel('Z[m]')


ax.plot(x1,y1,z1)


ax.scatter(x1[0],y1[0],z1[0],marker='o', color="orange")
ax.scatter(x1[len(x1)-1],y1[len(x1)-1],z1[len(x1)-1],marker='o', color="orange")


plt.show()


