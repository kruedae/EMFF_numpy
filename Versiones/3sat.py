#Importo las librerias
import numpy as np
import scipy.integrate as ode
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#Defino los momentos magneticos de ambos satelites
m1=np.array( [1*10**4.,1*10**4., 0.] )
m2=np.array( [-1*10**4., 0., 0.] )
m3=np.array( [-1*10**4.,0., 0.] )

#Defino el sistema de ecuaciones a solucionar funcion a solucionar tras haberlo separado en ecuaciones de primer orden

def ec2(t,RR):
    x1,vx1,y1,vy1,z1,vz1,  x2,vx2,y2,vy2,z2,vz2, x3,vx3,y3,vy3,z3,vz3  = RR[0], RR[1],RR[2], RR[3], RR[4], RR[5], RR[6], RR[7], RR[8],RR[9],RR[10],RR[11],RR[12], RR[13], RR[14],RR[15],RR[16],RR[17]

    #Definimos el vector posicion del primer y el segundo satelite
    m1=np.array( [-1*10**4., 0., 0.] )
    m2=np.array( [1*10**4., 0., 0.] )
    m3=np.array( [1*10**4., 0.,-1*10**4.] )


    R1=np.array([x1,y1,z1])
    R2=np.array([x2,y2,z2])
    R3=np.array([x3,y3,z3])


    #Defino el vector posicion relativa del segundo respecto al primero
    r12=R2-R1
    r13=R3-R1
    r23=R3-R2
    nr12=np.sqrt(np.dot(r12,r12)) #norma del vector r
    nr13=np.sqrt(np.dot(r13,r13)) #norma del vector r
    nr23=np.sqrt(np.dot(r23,r23)) #norma del vector r
    Mew=4*np.pi*10**-7. #constante de permeabilidad

    #escribimos las funciones para la aceleracion de cada satelite
    fm1= 3.*Mew/(4.*np.pi) * ( (-np.dot(m1,m2)*r12/(nr12**5.) - np.dot(m1,m3)*r13/(nr13**5.)  )  \
            - ( np.dot( m1, r12 )*m2/(nr12**5.)+np.dot( m1, r13 )*m3/(nr13**5.) )  \
            - ( np.dot(m2,r12)*m1/(nr12**5.)+ np.dot(m3,r13)*m1/(nr13**5.) ) \
            + 5.* ( np.dot(m1,r12)*np.dot(m2,r12)*r12/(nr12**7.)+np.dot(m1,r13)*np.dot(m3,r13)*r13/(nr13**7.))) 

    fm2= 3.*Mew/(4.*np.pi) * ( (-np.dot(m2,m1)*-r12/(nr12**5.)-np.dot(m2,m3)*r23/(nr23**5.)  )  \
            - ( np.dot( m2, -r12 )*m1/(nr12**5.)+np.dot( m2, r23 )*m3/(nr23**5.) )  \
            - ( np.dot(m1,-r12)*m2/(nr12**5.)+ np.dot(m3,r23)*m2/(nr23**5.) ) \
            + 5.* ( np.dot(m2,-r12)*np.dot(m1,-r12)*-r12/(nr12**7.)+np.dot(m2,r23)*np.dot(m3,r23)*r23/(nr23**7.))) 

    fm3= 3.*Mew/(4.*np.pi) * ( (-np.dot(m3,m2)*-r23/(nr23**5.)-np.dot(m3,m1)*-r13/(nr13**5.)  )  \
            - ( np.dot( m3, -r23 )*m2/(nr23**5.)+np.dot( m3, -r13 )*m1/(nr13**5.) )  \
            - ( np.dot(m2,-r23)*m3/(nr23**5.)+ np.dot(m1,-r13)*m3/(nr13**5.) ) \
            + 5.* ( np.dot(m3,-r23)*np.dot(m2,-r23)*-r23/(nr23**7.)+np.dot(m3,-r13)*np.dot(m1,-r13)*-r13/(nr13**7.))) 
 
    #Terminos gravitacionales
    Me=3.986004*10**14.

#    g1=Me*R1/( np.sqrt(np.dot(R1,R1))**3. )
#    g2=Me*R2/( np.sqrt(np.dot(R2,R2))**3. )
#    g3=Me*R3/( np.sqrt(np.dot(R3,R3))**3. )
   
    dv1=fm1#-g1
    dv2=fm2#-g2
    dv3=fm3#-g3

    return [vx1,dv1[0],vy1,dv1[1],vz1,dv1[2], vx2,dv2[0],vy2,dv2[1],vz2,dv2[2], vx3,dv3[0],vy3,dv3[1],vz3,dv3[2]]

#RR0=[35786000,0.,  0.,1075.,  0.,0.,    35786010,0., 0.,1075.,  0.,0.,    35786005,0., 0.,1075.,  0.,0.]

#vector de condiciones iniciales[Posx1,Velx1,Posy1,Vely1,...,Posx2,Velx2,...]
RR0=[0,0.,  0.,0.,  0.,0.,     -10,0., 0.,0,  0.,0.,    6,0., 0.,0.,  0.,0.]#vector de condiciones iniciales


RR1=ode.solve_ivp(ec2,[0., 30.],RR0,method='LSODA', atol=1e-6, rtol=1e-10)
RR=RR1.y

x1=RR[0]
y1=RR[2]
z1=RR[4]

x2=RR[6]
y2=RR[8]
z2=RR[10]

x3=RR[12]
y3=RR[14]
z3=RR[16]

#Imprimimos la funcion
ax=plt.axes

fig=plt.figure()
ax=Axes3D(fig)
plt.title('Posiciones de tres satelites EMFF')
ax.set_xlabel('x[m]')
ax.set_ylabel('y[m]')
ax.set_zlabel('z[m]')

ax.scatter(x1[len(x1)-1],y1[len(x1)-1],z1[len(x1)-1])
ax.scatter(x2[len(x1)-1],y2[len(x1)-1],z2[len(x1)-1])
ax.scatter(x3[len(x1)-1],y3[len(x1)-1],z3[len(x1)-1])

ax.plot(x1,y1,z1)

ax.plot(x2,y2,z2)

ax.plot(x3,y3,z3)

plt.show()


