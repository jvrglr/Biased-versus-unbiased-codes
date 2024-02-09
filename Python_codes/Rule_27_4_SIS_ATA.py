import numpy as np
import time
from numpy.random import binomial as iran_bin
from numpy.random import uniform as uniform
import matplotlib.pyplot as plt


def Gillespie(IC,tmax=20):
    """Compute I(t) from I(t=0)=IC using Gillespie algorithm"""
    """Returns density t_f,I/N(t_f). t_f=min(tmax,min_t{I(t)<=0})"""
    """beta,mu and N are implicit parameters that have to be defined outside this function"""
    I=IC;t=0.0
        
    while((t<tmax)and(I>0)):
        Wp=beta*I*(N-I)/N
        Wm=mu*I
        Wt=Wp+Wm
        t+=-np.log(uniform())/Wt
        p=Wp/Wt
        u=uniform()
        if (u<p):
            I=I+1
        else:
            I=I-1
    return t,I/N

def binomial(IC,Dt=1,tmax=20):
    """Compute I(t) from I(t=0)=IC using syncrhonous update with binomial algorithm"""
    """Returns density t_f,I/N(t_f). t_f=min(tmax,min_t{I(t)<=0})"""
    """beta,mu and N are implicit parameters that have to be defined outside this function"""
    I=IC;t=0.0
    mudt=mu*Dt
    q=1.0-np.exp(-mudt)
    bdt=beta/N*Dt

    while((t<tmax)and(I>0)):
        t+=Dt
        p=1.0-np.exp(-I*bdt)
        DS=iran_bin(N-I,p)
        DI=iran_bin(I  ,q)
        I+=DS-DI #Only can happen one thing

    return t,I/N


Ns=[1000]
beta=4.0;mu=1.0;gamma=1.0;tmax=20;IC=10
for N in Ns:
    print("N=",N)

    e0=0.25 #Scaling of errors binomial (err=e0*t)
    Ms=[2,10,30,50,100,300,1000,3000,5000,10000]
    
        
    tBs=np.zeros(len(Ms))
    tGs=np.zeros(len(Ms))
    eGs=np.zeros(len(Ms)) #Real precision measured
    eBs=np.zeros(len(Ms))
    MBs=np.zeros(len(Ms))
    dts=np.zeros(len(Ms))

    for jj,M in enumerate(Ms):
        #-----------------------Measure Gillespie
        I=0;I2=0
        start = time.time()
        for ll in range(int(M)):
            dumt,dumI= Gillespie(IC=IC,tmax=tmax)
            I+=dumI;I2+=dumI*dumI
        end = time.time()
        tG=end-start
        #-----------------------Prepare binomial
        I=I/M;I2=I2/M
        eG=np.sqrt(abs(I2-I**2)/M) #I had an error in this line!
        Dt=eG/e0/3
        Mdt=round(M*2.25)
        #-----------------------Measure binomial
        I=0;I2=0
        start = time.time()
        for ll in range(Mdt):
            dumt,dumI= binomial(IC=IC,Dt=Dt,tmax=tmax)
            I+=dumI;I2+=dumI*dumI #In order to have a fair comparison in CPU times
        end = time.time()
        tB=end-start
        I=I/Mdt;I2=I2/Mdt
        eB=np.sqrt(abs(I2-I**2)/Mdt)
        #----- print/save
        print("M=",M,"tG=",tG,"tB=",tB,"alpha=",tB/tG,"dt=",Dt,"MB=",Mdt)
        print("-----------------------------------------------------------------")
        tGs[jj]=tG
        tBs[jj]=tB
        MBs[jj]=Mdt
        dts[jj]=Dt
        eGs[jj]=eG
        eBs[jj]=eB
    tBs=np.array(tBs);tGs=np.array(tGs)
    alpha=tBs/tGs
    plt.plot(eGs,alpha)

    np.savetxt("data/ATA_TCPU_theoretical_1_realiz_N_"+str(N)+".dat",np.c_[tBs,tGs,eGs,eBs,MBs,dts]) 
