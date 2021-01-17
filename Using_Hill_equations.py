import numpy as np
import scipy.integrate as ode
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#Magnetic moments for both satellites

m1 = np.array([1*10**2,-1*10**2,0])
m2 = np.array([1*10**2,1*10**2.,-1*10**2])
m3 = np.array([1*10**3.,0,0])
m4 = np.array([1*10**3.,0,0])



#System of equations is defined according to hills equations and the forces described in Ashun's text

def eq_system(t,RR,m1,m2):
    x1,vx1,y1,vy1,z1,vz1,    x2,vx2,y2,vy2,z2,vz2,      xr,vxr,yr,vyr,zr,vzr =    RR[0],RR[1],RR[2],RR[3],RR[4],RR[5],     RR[6],RR[7],RR[8],RR[9],RR[10],RR[11],    RR[12],RR[13], RR[14],RR[15],RR[16],RR[17]



    #Definimos el vector posicion del primer y el segundo satelite
    #Position vector
    R1 = np.array([x1,y1,z1])
    R2 = np.array([x2,y2,z2])


    #Masas
    #Mass
    Ms1 = 20
    Ms2 = 20

    #Defino el vector posicion relativa del segundo respecto al primero
    r = R2-R1
    nr = np.sqrt(np.dot(r,r)) #norma del vector r
    Mew = 4*np.pi*10**-7. #constante de permeabilidad magnetica del vacio

    #escribimos las fuerzas magneticas de cada satelite
    fm1 = 3.*Mew/(4.*np.pi)*(-np.dot(m1, m2)*r/(nr**5.)
            -np.dot(m1, r)*m2/(nr**5.)
            -np.dot(m2, r)*m1/(nr**5.)
            +5.*np.dot(m1, r)*np.dot(m2, r)*r/(nr**7.))/Ms1
    fm2 = 3.*Mew/(4.*np.pi)*(-np.dot(m2,m1)*(-r)/(nr**5.)
            -np.dot(m2,-r)*m1/(nr**5.)
            -np.dot(m1,-r)*m2/(nr**5.)
            +5.*np.dot(m2,-r)*np.dot(m1,-r)*(-r)/(nr**7.)) /Ms2

    #Terminos gravitacionales
    Me = 3.986004*10**14.#Constante de gravitacion para la tierra

    g1 = Me*R1/(np.sqrt(np.dot(R1,R1))**3.)
    g2 = Me*R2/(np.sqrt(np.dot(R2,R2))**3.)

    #Describo las ecuaciones para las aceleraciones
    dv1 = fm1-g1
    dv2 = fm2-g2
    
    #RELATIVO

    Rcm = R1
    rcm = R1/np.sqrt(np.dot(R1,R1))

    Vcm = np.array([vx1,vy1,vz1])
    vcm = Vcm/( np.sqrt( np.dot(Vcm,Vcm) ) )
    
    
    MeciRo = np.array( [rcm, vcm, np.cross(vcm,rcm)]  )
#    MeciRo = np.array( [-vcm, rcm, np.cross(-vcm,rcm)]  )

    #Velocidad angular:
    wro = np.sqrt( Me/np.sqrt(np.dot(Rcm,Rcm)**3.) )

    #Sistema relativo
    fmr1 = np.matmul(MeciRo,fm1)
    fmr2 = np.matmul(MeciRo,fm2)
    nRk = np.sqrt(np.dot(R1,R1))
    nRi = np.sqrt(np.dot(R2,R2))


#    dvxr = (2*wro*vyr + xr*wro**2. - (Me*xr/(nRi**3.)) - fmr1[0] +fmr2[0])
#    dvyr = (-2*wro*vxr + yr*wro**2. - Me*( nRk+yr )/nRi**3. + Me/nRk**2. - fmr1[1] +fmr2[1])
#    dvzr = -Me*zr/(nRi**3.) - fmr1[2] + fmr2[2]
    dvxr = (3.*wro**2.*xr+2.*wro*vyr - fmr1[0] +fmr2[0])
    dvyr = (-2.*wro*vxr - fmr1[1] +fmr2[1])
    dvzr = -wro**2.*zr - fmr1[2] + fmr2[2]



    return [vx1,dv1[0],vy1,dv1[1],vz1,dv1[2], vx2,dv2[0],vy2,dv2[1],vz2,dv2[2], vxr,dvxr,vyr,dvyr,vzr,dvzr]




RR0 = [6917000.,0.,  0.,7591.19,  0.,0.,    6917005., 0., 0.,7591.19, 0.,0.,    5, 0,0,0,0,0]#vector de condiciones iniciales

#Solucionamos con las condiciones anteriores usando un Runge-Kutta de orden 3 (el metodo puede variar)
RR1 = ode.solve_ivp(eq_system,[0., 5800.],RR0,method='RK45', atol=1e-4, rtol=1e-6,max_step=2,args=[m1,m2])
RR = RR1.y

RR22 = ode.solve_ivp(eq_system,[0., 5800.],RR0,method='RK45', atol=1e-4, rtol=1e-6,max_step=2,args=[m3,m4])
RR2 = RR22.y

#Guardamos las soluciones
xr2 = RR2[12]
yr2 = RR2[14]
zr2 = RR2[16]

x1 = RR[0]
y1 = RR[2]
z1 = RR[4]

x2 = RR[6]
y2 = RR[8]
z2 = RR[10]

xr = RR[12]
yr = RR[14]
zr = RR[16]

#Imprimimos la funcion

#2-D
plt.show()
plt.title('posicion de ambos satelites vistos desde la tierra')
plt.axis('equal')
plt.plot(x1,y1)
plt.plot(x2,y2)
plt.plot(x1[len(xr)-1],y1[len(xr)-1], marker='o', markersize=3, color="blue")
plt.plot(x2[len(xr)-1],y2[len(xr)-1], marker='o', markersize=3, color="orange")
plt.show()

#2-D
plt.title('Posicion relativa')
plt.xlabel('x[m]')
plt.ylabel('y[m]')
plt.axis('equal')
plt.plot(xr,yr)
plt.plot(xr[0],yr[0], marker='o', markersize=2, color="blue")
plt.plot(xr[0],yr[0], marker='o', markersize=2, color="orange")

plt.show()




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


ax.plot(xr,yr,zr)

plt.show()

#Marco relativo

ax = plt.axes
fig = plt.figure()
ax = Axes3D(fig)

plt.title('Trayectoria para dos satelites')
ax.plot(x2,y2,z2)
ax.plot(x1,y1,z1)

ax.scatter(x1[0],y1[0],z1[0],marker='o', color="orange")
ax.scatter(x1[len(x1)-1],y1[len(x1)-1],z1[len(x1)-1],marker='o', color="orange")
ax.scatter(x2[len(x1)-1],y2[len(x1)-1],z2[len(xr)-1],marker='o', color="blue")
ax.scatter(x2[0],y2[0],z2[0],marker='o', color="blue")

ax.set_xlabel('X[m]')
ax.set_ylabel('Y[m]')
ax.set_zlabel('Z[m]')

plt.show()


