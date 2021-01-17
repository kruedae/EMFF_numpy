#Programa para calcular la solucion a una ecuacion de orden superior
#Se transforma la ecuacion en un sistema de ecuaciones
import math
import pylab as p
#y''-2y'+2y=e^(2t)sent
def f(n,x,w1,w2):
    if n==0:
        return w2
    if n!=0:
        return math.exp(2*x)*math.sin(x)+2*w2-2*w1
#Funcion solucion
def sol(x):
    return math.exp(x)*((0.8-0.4*math.exp(x))*math.cos(x)+(0.2*math.exp(x)-0.8)*math.sin(x))

a=0.0
b=1
m=2
N=11
k1=[0.0]*m
k2=[0.0]*m
k3=[0.0]*m
k4=[0.0]*m
Y=[0.0]*m
#Vector Y, en la primer fila contiene las aproximaciones a la solucion y en la segunda las aproximaciones a su derivada
for k in range(0,m):
    Y[k]=[0.0]*N
#Valores iniciales
Y[0][0]=0.4
Y[1][0]=-0.6
XX=[0.0]*N
h=float((b-a)/N)
#N para el numero de puntos deseados
for i in range(1,N):
    XX[i]=a+h*float(i)
    # m para el orden de la ecuacion. En cada paso resuelvo tanto para w1 como para w2
    for j in range(0,m):
        k1[j]=h*f(j,XX[i-1],Y[0][i-1],Y[1][i-1])
    for j in range(0,m):
        k2[j]=h*f(j,XX[i-1]+h/2,Y[0][i-1]+k1[0]/2,Y[1][i-1]+k1[1]/2)
    for j in range(0,m):
        k3[j]=h*f(j,XX[i-1]+h/2,Y[0][i]+k2[0]/2,Y[1][i-1]+k2[1]/2)
    for j in range(0,m):
        k4[j]=h*f(j,XX[i-1]+h,Y[0][i-1]+k3[0],Y[1][i-1]+k3[1])
    for j in range(0,m):
        Y[j][i]=Y[j][i-1]+(1.0/6.0)*(k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j])
#Grafico la solucion
NS=101
XS=[0.0]*NS
YS=[0.0]*NS
hs=float((b-a)/NS)
for k in range (0,NS):
    XS[k]=a+hs*float(k)
    YS[k]=sol(XS[k])
p.plot(XX,Y[0],'p')
p.plot(XS,YS)
p.show()

