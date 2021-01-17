#CALCULO DEL NUMERO DE SATELITES A PARTIR DE UNA MODIFICACION DEL METODO WALKER MODIFICADO
import numpy as np

def walker(a,eps):
    Re=6371#radio de la tierra
    arch = open('walker.txt','wb')

    #CALCULOS 

    #Calculamos phi ()
    psi=np.arcsin(Re*np.cos(eps)/(Re+a))

    #calculamos theta()
    theta=np.pi/2-(psi+eps)
    print('theta',theta)
    arch.write('theta: '+str(theta)+'\n')

    #Calculamos r
    r=Re*theta
    print('r',r)
    arch.write('radio de cobertura: '+str(r)+'\n')

    #Calculo el area por satelite
    A=3*np.sqrt(3)*r**2./2

    print("Cobertura por satelite",A)
    arch.write('Cobertura por satelite: '+str(A)+'\n')

    #Inclinacion
    i=1.152-theta/2
    print('inclinacion max: ',i*180/np.pi)
    arch.write('inclinacion: '+str(i)+'\n')

    #Area de los polos
    rpolo=Re*(np.pi/2-1.152)
    Apolo=np.pi*rpolo**2.
    print('Area polo: ',Apolo)

    #Area total
    Atierra=510.1*10**6.
    Atotal=Atierra-Apolo
    print('Area total de cobertura: ',Atotal)
    arch.write('Area total de cobertura: '+str(Atotal)+'\n')

    #Numero minimo de satelites
    Nsatelites=int(Atotal/A)
    print('Numero de satelites',Nsatelites)
    arch.write('Numero de satelites: '+str(Nsatelites)+'\n')

    #Planos orbitales
    P=int(np.pi*Re/(r*np.cos(np.pi/6)))
    print('Planos orbitales: ',P)
    arch.write('Planos orbitales: '+str(P)+'\n')

    #Satelites por plano
    SP=Nsatelites/P
    aux=0
    Nfinalsat=Nsatelites
    while aux==0:
        if Nfinalsat%P==0:
            aux=1
        else:
            Nfinalsat=Nfinalsat+1
    print('# final de satelites: ',Nfinalsat)
    arch.write('Numero final de satelites: '+str(Nfinalsat)+'\n')
    print('satelite por plano: ',Nfinalsat/P)
    arch.write('Satelites por plano: '+str(Nfinalsat/P)+'\n')

    #Parametro de espacio relativo entre fases
    f=1
    print('Espacio relativo entre fases: ',f)

    #Diferencia entre fase

    df=360.*f/(Nfinalsat)
    print('Diferencia entre fases: ',df)
    arch.write('Diferencia entre fases: '+str(df)+'\n')
    
    arch.close()

#PARAMETROS DE ENTRADA
#Defino la altura
A=500
#Defino el angulo de elevacion minimo
ep=10*np.pi/180

walker(A,ep)
