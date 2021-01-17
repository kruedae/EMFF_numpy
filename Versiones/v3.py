import numpy as np
from scipy import linalg as lg
from scipy.integrate import odeint
import matplotlib.pyplot as plt
m1=np.array([1.*10.**2.,1.*10.**2.])
m2=np.array([1.*10.**4.,1.*10.**2.])
def ec2(RR,t,m1,m2):
    x1,vx1,y1,vy1,x2,vx2,y2,vy2 = RR[0], RR[1],RR[2], RR[3], RR[4], RR[5], RR[6], RR[7]
    R1=np.array([x1,y1])
    R2=np.array([x2,y2])
    r=R2-R1
    nr=np.sqrt(np.dot(r,r))
    dv1=-(6.*10.**13.)*R1/np.sqrt(np.dot(R1,R1))**3. -np.dot(m1,m2)*r/(nr**5.)-np.dot(m1,r)*m2/(nr**5.)-np.dot(m1,r)*m1/(nr**5.)+5*np.dot(m1,r)*np.dot(m2,r)*r/(nr**7.)-3.*(6.*10**13)*(6.3*10**6.)*0.001082/(2.*np.sqrt(np.dot(R1,R1)**5.))*R1

    dv2=-(6.*10.**13.)*R2/np.sqrt(np.dot(R2,R2))**3. -np.dot(m2,m1)*(-r)/(nr**5.)-np.dot(m2,-r)*m1/(nr**5.)-np.dot(m2,-r)*m2/(nr**5.)+5*np.dot(m2,-r)*np.dot(m1,-r)*(-r)/(nr**7.)-3.*(6.*10**13)*(6.3*10**6.)*0.001082/(2.*np.sqrt(np.dot(R2,R2)**5.))*R2
    return [vx1,dv1[0],vy1,dv1[1],vx2,dv2[0],vy2,dv2[1]]
RR0=[35786000.,0.,0.,1075.,35786050.,0.,0.,1075.]
ts=np.linspace(0,150000,3000000)
RR=odeint(ec2,RR0,ts,args=(m1,m2))
x1=RR[:,0]#Para cada ts, Zs toma un vector, aca selecciono la primer columna de cada uno de esos vecs
y1=RR[:,2]
print(x1)
x2=RR[:,4]#Para cada ts, Zs toma un vector, aca selecciono la primer columna de cada uno de esos vecs
y2=RR[:,6]
plt.plot(x1,y1)
plt.plot(x2,y2)
plt.show()


