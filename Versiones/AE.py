
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
    phi,the,psi, dphi,dthe,dpsi = RR[0], RR[1],RR[2],  RR[3], RR[4], RR[5]
    
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




RR0=[0,0,0,2,2,2]#vector de condiciones iniciales

#Solucionamos con las condiciones anteriores usando un Runge-Kutta de orden 3 (el metodo puede variar)
RR1=ode.solve_ivp(ec2,[0., 2.],RR0,method='RK45', atol=1e-4, rtol=1e-6,max_step=0.05,args=[I])
RR=RR1.y
#Guardamos las soluciones
phi=RR[0]
the=RR[1]
psi=RR[2]
print(RR[0])
#Rotacion
def rota(phi,the,psi,x1,y1,y3):
    for i in range(1,len(x1)):
        MR=np.array([[np.cos(the[i])*np.cos(psi[i]), np.cos(the[i])*np.sin(psi[i]), np.sin(the[i])] ,
    
                [np.sin(the[i])*np.sin(phi[i])*np.cos(psi[i])-np.cos(phi[i])*np.sin(psi[i]) ,
                    np.sin(the[i])*np.sin(phi[i])*np.sin(psi[i])+np.cos(phi[i])*np.cos(psi[i]), np.cos(the[i])*np.sin(phi[i])] ,
    
                [np.sin(the[i])*np.cos(phi[i])*np.cos(psi[i])+np.sin(phi[i])*np.sin(psi[i]),
                    np.sin(the[i])*np.cos(phi[i])*np.sin(psi[i])-np.sin(phi[i])*np.cos(psi[i]), np.cos(the[i])*np.cos(phi[i])]
                ])
        v1=np.array([x1[i],y1[i],z1[i]])
        v1=np.matmul(MR,v1)
        x1[i],y1[i],z1[i]=v1[0],v1[1],v1[2]
    return [x1,y1,z1]

#Vectores a rotar
x1=np.array([1.]*len(phi))
y1=np.array([1.]*len(phi))
z1=np.array([1.]*len(phi))

x2=np.array([-1.]*len(phi))
y2=np.array([1.]*len(phi))
z2=np.array([1.]*len(phi))

x3=np.array([1.]*len(phi))
y3=np.array([-1.]*len(phi))
z3=np.array([1.]*len(phi))

x4=np.array([-1.]*len(phi))
y4=np.array([-1.]*len(phi))
z4=np.array([1.]*len(phi))

x5=np.array([1.]*len(phi))
y5=np.array([1.]*len(phi))
z5=np.array([2]*len(phi))

x6=np.array([-1.]*len(phi))
y6=np.array([1.]*len(phi))
z6=np.array([2.]*len(phi))

x7=np.array([1.]*len(phi))
y7=np.array([-1.]*len(phi))
z7=np.array([2.]*len(phi))

x8=np.array([-1.]*len(phi))
y8=np.array([-1.]*len(phi))
z8=np.array([2.]*len(phi))
print('R:',RR[0])
print(z5)
#Aplico la rotacion
print(rota(phi,the,psi,x5,y5,z5)[2])
x1,y1,z1=rota(phi,the,psi,x1,y1,z1)
x2,y2,z2=rota(phi,the,psi,x2,y2,z2)
x3,y3,z3=rota(phi,the,psi,x3,y3,z3)
x4,y4,z4=rota(phi,the,psi,x4,y4,z4)

x5,y5,z5=rota(phi,the,psi,x5,y5,z5)
print(z5)
x6,y6,z6=rota(phi,the,psi,x6,y6,z6)
x7,y7,z7=rota(phi,the,psi,x7,y7,z7)
x8,y8,z8=rota(phi,the,psi,x8,y8,z8)
#Imprimimos la funcion

#3-D
#Marco de la tierra
ax=plt.axes

fig=plt.figure()
ax=Axes3D(fig)
plt.title('Trayectoria relativa del segundo satelite')
ax.scatter(x1[0],y1[0],z1[0],marker='o', color="blue")
#ax.scatter(x1[len(x1)-1],y1[len(x1)-1],z1[len(z1)-1],marker='o', color="orange")
ax.scatter(x2[0],y2[0],z2[0],marker='o', color="blue")
#ax.scatter(x2[len(x1)-1],y2[len(x1)-1],z2[len(z2)-1],marker='o', color="orange")
ax.scatter(x3[0],y3[0],z3[0],marker='o', color="blue")
#ax.scatter(x3[len(x1)-1],y3[len(x1)-1],z3[len(z2)-1],marker='o', color="orange")
ax.scatter(x4[0],y4[0],z4[0],marker='o', color="blue")
#ax.scatter(x4[len(x1)-1],y4[len(x1)-1],z4[len(z2)-1],marker='o', color="orange")

ax.scatter(x5[0],y5[0],z5[0],marker='o', color="blue")
#ax.scatter(x5[len(x1)-1],y5[len(x1)-1],z5[len(z1)-1],marker='o', color="orange")
ax.scatter(x6[0],y6[0],z6[0],marker='o', color="blue")
#ax.scatter(x6[len(x1)-1],y6[len(x1)-1],z6[len(z2)-1],marker='o', color="orange")
ax.scatter(x7[0],y7[0],z7[0],marker='o', color="blue")
#ax.scatter(x7[len(x1)-1],y7[len(x1)-1],z7[len(z2)-1],marker='o', color="orange")
ax.scatter(x8[0],y8[0],z8[0],marker='o', color="blue")
#ax.scatter(x8[len(x1)-1],y8[len(x1)-1],z8[len(z2)-1],marker='o', color="orange")
ax.set_xlabel('X[m]')
ax.set_ylabel('Y[m]')
ax.set_zlabel('Z[m]')

print(x5[0],y5[0],z5[0])
#ax.plot(x1,y1,z1)
#ax.plot(x2,y2,z2)
#ax.plot(x3,y3,z3)
#ax.plot(x4,y4,z4)

#ax.plot(x5,y5,z5)
#ax.plot(x6,y6,z6)
#ax.plot(x7,y7,z7)
#ax.plot(x8,y8,z8)


plt.show()
