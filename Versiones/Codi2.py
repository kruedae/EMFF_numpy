#Importo las librerias
import numpy as np
from scipy import linalg as lg
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#Defino los momentos magneticos de ambos satelites
m1=np.array( [1*10**4., 0] )
m2=np.array( [-1*10**4., 0] )

#Defino el sistema de ecuaciones a solucionar funcion a solucionar tras haberlo separado en ecuaciones de primer orden
def ec2(RR,t,m1,m2):
    x1,vx1,y1,vy1,x2,vx2,y2,vy2 = RR[0], RR[1],RR[2], RR[3], RR[4], RR[5], RR[6], RR[7]

    #Definimos el vector posicion del primer y el segundo satelite

    R1=np.array([x1,y1])
    R2=np.array([x2,y2])
    #Defino el vector posicion relativa del segundo respecto al primero
    r=R2-R1
    nr=np.sqrt(np.dot(r,r)) #norma del vector r
    Mew=4*np.pi*10**-7. #constante de permeabilidad

    #escribimos las funciones para la aceleracion de cada satelite
    dv1= 3.*Mew/( 4.*np.pi )*( -np.dot(m1,m2 )*r/( nr**5. )-np.dot( m1, r )*m2/(nr**5.)-np.dot(m2,r)*m1/(nr**5.)+5.*np.dot(m1,r)*np.dot(m2,r)*r/(nr**7.))
    dv2= 3.*Mew/(4.*np.pi)*(-np.dot(m2,m1)*(-r)/(nr**5.)-np.dot(m2,-r)*m1/(nr**5.)-np.dot(m1,-r)*m2/(nr**5.)+5.*np.dot(m2,-r)*np.dot(m1,-r)*(-r)/(nr**7.))
    return [vx1,dv1[0],vy1,dv1[1],vx2,dv2[0],vy2,dv2[1]]
RR0=[10.,0.,0.,0.,0.,0.,0.,0.]
ts=np.linspace(0,20,1000)
RR=odeint(ec2,RR0,ts,args=(m1,m2))
x1=RR[:,0]#Para cada ts, Zs toma un vector, aca selecciono la primer columna de cada uno de esos vecs
x2=RR[:,4]
print(x1)
print(x2)
plt.plot(ts,x1)
plt.plot(ts,x2)
plt.xlabel('Tiempo[s]')
plt.ylabel('Posicion[m]')
plt.title('Posiciones de dos satelites EMFF')
plt.show()

