#Importo las librerias
import numpy as np
import scipy.integrate as ode
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#Defino los momentos magneticos de ambos satelites

#Defino el sistema de ecuaciones a solucionar funcion a solucionar tras haberlo separado en ecuaciones de primer orden

def ec2(t,RR):
    x1,vx1,y1,vy1,z1,vz1,  x2,vx2,y2,vy2,z2,vz2,  xr,vxr,yr,vyr,zr,vzr   = RR[0], RR[1],RR[2], RR[3], RR[4], RR[5], RR[6], RR[7], RR[8],RR[9],RR[10],RR[11],   RR[12], RR[13], RR[14],RR[15],RR[16],RR[17]


    #Definimos los momentos magneticos

    m1=np.array( [1*10**4., 0., 0.] )
    m2=np.array( [-1*10**4., 0., 0.] )

    #Definimos el vector posicion del primer y el segundo satelite
    R1=np.array([x1,y1,z1])
    R2=np.array([x2,y2,z2])
    #Masas
    Ms1=1
    Ms2=1

    #Defino el vector posicion relativa del segundo respecto al primero
    r = R2-R1
    nr = np.sqrt(np.dot(r,r)) #norma del vector r
    Mew = 4*np.pi*10**-7. #constante de permeabilidad

    #escribimos las funciones para la aceleracion de cada satelite
    fm1 = 3.*Mew/( 4.*np.pi )*( -np.dot(m1,m2 )*r/( nr**5. )-np.dot( m1, r )*m2/(nr**5.)-np.dot(m2,r)*m1/(nr**5.)+5.*np.dot(m1,r)*np.dot(m2,r)*r/(nr**7.)) /Ms1
    fm2 = 3.*Mew/(4.*np.pi)*(-np.dot(m2,m1)*(-r)/(nr**5.)-np.dot(m2,-r)*m1/(nr**5.)-np.dot(m1,-r)*m2/(nr**5.)+5.*np.dot(m2,-r)*np.dot(m1,-r)*(-r)/(nr**7.)) /Ms2

    #Terminos gravitacionales
    Me=3.986004*10**14.

    g1 = Me*R1/( np.sqrt(np.dot(R1,R1))**3. )
    g2 = Me*R2/( np.sqrt(np.dot(R2,R2))**3. )

    #Describo las aceleraciones
    dv1 = fm1-g1
    dv2 = fm2-g2
    
    #RELATIVO
    Rcm = (R1*Ms1+R2*Ms2)/(Ms1+Ms2)
    rcm= Rcm/( np.sqrt(np.dot(Rcm,Rcm)) )


    Vcm = np.array( [Ms1*vx1+Ms2*vx2, Ms1*vy1+Ms2*vy2, Ms1*vz1+Ms2*vz2] )/(Ms1+Ms2)
    vcm= -Vcm/( np.sqrt( np.dot(Vcm,Vcm) ) )


    MeciRo = np.array( [-vcm, rcm, np.cross(-vcm,rcm)]  )

    #Velocidad angular:
    wro = np.sqrt( Me/np.sqrt(np.dot(Rcm,Rcm)**3.) )

    #Sistema relativo
    fmr1 = np.matmul(MeciRo,fm1)
#    print(fmr1)
    fmr2 = np.matmul(MeciRo,fm2)
    nRk = np.sqrt(np.dot(R1,R1))
    nRi = np.sqrt(np.dot(R2,R2))


    dvxr = (2*wro*vyr + xr*wro**2. - (Me*xr/(nRi**3.)) - fmr1[0] +fmr2[0])
    dvyr = (-2*wro*vxr + yr*wro**2. - Me*( nRk+yr )/nRi**3. + Me/nRk**2. - fmr1[1] +fmr2[1])
    dvzr = -Me*zr/(nRi**3.) - fmr1[2] + fmr2[2]



    return [vx1,dv1[0],vy1,dv1[1],vz1,dv1[2], vx2,dv2[0],vy2,dv2[1],vz2,dv2[2], vxr,dvxr,vyr,dvyr,vzr,dvzr]





RR0=[35786000,0.,  0.,1075.,  0.,0.,    35786010,0., 0.,1075., 0.,0.,  10,0,0,0,0,0]#vector de condiciones iniciales

RR1=ode.solve_ivp(ec2,[0., 25000.],RR0,method='RK45', atol=1e-4, rtol=1e-6,max_step=1)
RR=RR1.y

x1=RR[0]
y1=RR[2]
z1=RR[4]

x2=RR[6]
y2=RR[8]
z2=RR[10]

xr=RR[12]
yr=RR[14]
zr=RR[16]

#Imprimimos la funcion
ax=plt.axes

fig=plt.figure()
ax=Axes3D(fig)
plt.title('Posiciones de dos satelites EMFF')

#ax.plot(x1,y1,z1)

#ax.plot(x2,y2,z2)

ax.plot(xr,yr,zr)

plt.show()


#movimiento

ax=plt.axes

fig=plt.figure()
ax=Axes3D(fig)
plt.title('Posiciones de dos satelites EMFF')

ax.plot(x1,y1,z1)

ax.plot(x2,y2,z2)

#ax.plot(xr,yr,zr)

plt.show()

plt.plot(x1,y1)
plt.plot(x2,y2)
plt.show()

plt.plot(xr,yr)
plt.show()

