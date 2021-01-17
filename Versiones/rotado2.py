import numpy as np
from scipy import linalg as lg
from scipy.integrate import odeint
import matplotlib.pyplot as plt
m1=np.array([-10000.,0.,0.])
m2=np.array([1000.,0.,0.])
def ec2(RR,t,m1,m2):
    x1,vx1,y1,vy1,x2,vx2,y2,vy2 = RR[0], RR[1],RR[2], RR[3], RR[4], RR[5], RR[6], RR[7]
    R1=np.array([x1,y1,0])
    R2=np.array([x2,y2,0])
    r=R2-R1
    print('r',r)
    nr=np.sqrt(np.dot(r,r))
    Fmr=3.*(10**-7)/(nr**4.)*np.array([2.*m1[0]*m2[0]-m1[1]*m2[1]-m1[2]*m2[2],-m1[0]*m2[1]-m1[1]*m2[0],-m1[0]*m2[2]-m1[2]*m2[0]])
    print('Fmr',Fmr)
    Ang2=np.pi/2*0.
    Ang1=np.arcsin(r[2]/nr)
    M1rot=np.array([[np.cos(Ang1),0,-np.sin(Ang1)],[0,1,0],[np.sin(Ang1),0,np.cos(Ang1)]])
    M2rot=np.array([[np.cos(Ang2),np.sin(Ang2),0],[-np.sin(Ang2),np.cos(Ang2),0],[0,0,1]])
    Mrot=np.dot(M1rot,M2rot)
    print('Mrot',Mrot)
    FmI=np.dot(Mrot.transpose(),Fmr)
    print('FmI',FmI)
    dv1=FmI-(6.*10.**13.)*R1/np.sqrt(np.dot(R1,R1))**3.
    dv2=-FmI-(6.*10.**13.)*R2/np.sqrt(np.dot(R2,R2))**3.
    return [vx1,dv1[0],vy1,dv1[1],vx2,dv2[0],vy2,dv2[1]]
RR0=[35786000.,0.,0.,1075.,35786050.,0.,0.,1075.]
ts=np.linspace(0,200000,1000000)
RR=odeint(ec2,RR0,ts,args=(m1,m2))
x1=RR[:,0]#Para cada ts, Zs toma un vector, aca selecciono la primer columna de cada uno de esos vecs
y1=RR[:,2]
x2=RR[:,4]
print(x1)
print(x2)
plt.plot(x1,y1)
plt.show()

