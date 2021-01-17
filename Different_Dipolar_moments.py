import numpy as np
from scipy import linalg as lg
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#Defino el momento dipolar
m1=np.array([1*10**4.,0])
m2=np.array([-1*10**4.,0])

#Escribo la funcion a solucionar
def ec2(RR,t,m1,m2):
    x1,vx1,y1,vy1,x2,vx2,y2,vy2 = RR[0], RR[1],RR[2], RR[3], RR[4], RR[5], RR[6], RR[7]
    
    #Defino los vectores posicion
    R1=np.array([x1,y1])
    R2=np.array([x2,y2])
    r=R2-R1#Posicion relativa
    nr=np.sqrt(np.dot(r,r))

    Mew=4*np.pi*10**-7.#constante de permeatividad

    #Escribo las ecuaciones de movimiento
    dv1=3.*Mew/(4.*np.pi)*(-np.dot(m1,m2)*r/(nr**5.)-np.dot(m1,r)*m2/(nr**5.)-np.dot(m2,r)*m1/(nr**5.)+5.*np.dot(m1,r)*np.dot(m2,r)*r/(nr**7.))
    dv2=3.*Mew/(4.*np.pi)*(-np.dot(m2,m1)*(-r)/(nr**5.)-np.dot(m2,-r)*m1/(nr**5.)-np.dot(m1,-r)*m2/(nr**5.)+5.*np.dot(m2,-r)*np.dot(m1,-r)*(-r)/(nr**7.))
    return [vx1,dv1[0],vy1,dv1[1],vx2,dv2[0],vy2,dv2[1]]
#Escribo las condiciones iniciales de los dos satelites para cada uno de los cuatro casos
RR0=[10.,0.,0.,0.,0.,0.,0.,0.]
RR01=[5.,0.,0.,0.,0.,0.,0.,0.]
RR02=[20.,0.,0.,0.,0.,0.,0.,0.]
RR02=[20.,0.,0.,0.,0.,0.,0.,0.]
RR03=[1.,0.,0.,0.,0.,0.,0.,0.]

#Soluciono las ecuaciones usando la funcion odeint de la biblioteca Scipy
ts=np.linspace(0,30,1000)
RR=odeint(ec2,RR0,ts,args=(m1,m2))
RR1=odeint(ec2,RR01,ts,args=(m1,m2))
RR2=odeint(ec2,RR02,ts,args=(m1,m2))
RR3=odeint(ec2,RR03,ts,args=(m1,m2))

#Organizo las soluciones en vectores independientes
x1=RR[:,0]
x2=RR[:,4]
x11=RR1[:,0]
x21=RR1[:,4]
x12=RR2[:,0]
x22=RR2[:,4]
x13=RR3[:,0]
x23=RR3[:,4]

#Grafico las soluciones
plt.subplot(221)
plt.plot(ts,x1,'-b',label='10 m')
plt.plot(ts,x2)
plt.ylabel('Posicion[m]')
plt.legend()


plt.subplot(222)
plt.plot(ts,x11,'-b',label='5 m')
plt.plot(ts,x21)

plt.legend()


plt.subplot(223)
plt.plot(ts,x12,'-b',label='20 m')
plt.plot(ts,x22)
plt.legend()

plt.subplot(224)
plt.plot(ts,x13,'-b',label='1 m')
plt.plot(ts,x23)

plt.legend()

plt.xlabel('Tiempo[s]')
plt.subplots_adjust(hspace=0)
plt.suptitle('Posicion vs tiempo para diferentes distancias iniciales')
plt.show()

