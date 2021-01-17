#Importo las librerias
import numpy as np
import scipy.integrate as ode
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#Defino los momentos magneticos de ambos satelites
m1=np.array( [1*10**4., 0., 0.] )
m2=np.array( [-1*10**4., 0., -1*10**4.] )

#Defino el sistema de ecuaciones a solucionar funcion a solucionar tras haberlo separado en ecuaciones de primer orden

def ec2(t,RR):
    x1,vx1,y1,vy1,z1,vz1,  x2,vx2,y2,vy2,z2,vz2 = RR[0], RR[1],RR[2], RR[3], RR[4], RR[5], RR[6], RR[7], RR[8],RR[9],RR[10],RR[11]

    #Definimos el vector posicion del primer y el segundo satelite
    
    m1=np.array( [1*10**4., 0., 0.] )
    m2=np.array( [-1*10**4., 0., 0.] )

    R1=np.array([x1,y1,z1])
    R2=np.array([x2,y2,z2])


    #Defino el vector posicion relativa del segundo respecto al primero
    r=R2-R1
    nr=np.sqrt(np.dot(r,r)) #norma del vector r
    Mew=4*np.pi*10**-7. #constante de permeabilidad

    #escribimos las funciones para la aceleracion de cada satelite
    fm1= 3.*Mew/( 4.*np.pi )*( -np.dot(m1,m2 )*r/( nr**5. )-np.dot( m1, r )*m2/(nr**5.)-np.dot(m2,r)*m1/(nr**5.)+5.*np.dot(m1,r)*np.dot(m2,r)*r/(nr**7.))
    fm2= 3.*Mew/(4.*np.pi)*(-np.dot(m2,m1)*(-r)/(nr**5.)-np.dot(m2,-r)*m1/(nr**5.)-np.dot(m1,-r)*m2/(nr**5.)+5.*np.dot(m2,-r)*np.dot(m1,-r)*(-r)/(nr**7.))

    #Terminos gravitacionales
    Me=3.986004*10**14.

    g1=Me*R1/( np.sqrt(np.dot(R1,R1))**3. )
    g2=Me*R2/( np.sqrt(np.dot(R2,R2))**3. )

    dv1=fm1-g1
    dv2=fm2-g2

    return [vx1,dv1[0],vy1,dv1[1],vz1,dv1[2], vx2,dv2[0],vy2,dv2[1],vz2,dv2[2]]

RR0=[35786000,0.,  0.,1075.,  0.,0.,    35786010,0., 0.,1075.,  0.,0.]#vector de condiciones iniciales

RR1=ode.solve_ivp(ec2,[0., 200000.],RR0,method='RK45', atol=1e-6, rtol=1e-10)
RR=RR1.y

#RR=ode.solve_ivp(ec2, t0=0.0, y0= RR0, t_bound=200.0)#solucionamos usando el paquete la funcion odeint

#print(RR)

#print(RR[0])
x1=RR[0]
y1=RR[2]
z1=RR[4]

x2=RR[6]
y2=RR[8]
z2=RR[10]

#Imprimimos la funcion
ax=plt.axes

fig=plt.figure()
ax=Axes3D(fig)
plt.title('Posiciones de dos satelites EMFF')

ax.plot(x1,y1,z1)

ax.plot(x2,y2,z2)

plt.show()


