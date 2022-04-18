

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 12:54:55 2021

@author: u103424502
"""

import math
import numpy as np
from matplotlib.pylab import plt
import numba
@numba.jit


def E_applied(x):
    pow=math.pow
    exp=math.exp
    cos=math.cos
    PI=math.pi
    log=math.log
    atan=math.atan
    
    lamda=790
    width=18
    GDD=0
    beta=(pow(width,2))/(8*log(2.00))
    ganma=1+pow(GDD,2)/(4*pow(beta,2))
    a=GDD/(8*pow(beta,2)*ganma)
    ep=0.5*atan(GDD/(2*beta))
    return 1/(2*pow(ganma,0.25))*exp(-(pow(x,2))/(4*beta*ganma))*cos(2*PI*(2.9979/lamda)*100*x+a*pow(x,2))
def plasmon(x):
    cos=math.cos
    exp=math.exp
    PI=math.pi
    
    lamda=790
    lifetime=5
    return 0 if x<0 else exp(-x/lifetime)*cos(2*PI*(2.9979/lamda)*100*x)
def E_plasmon(x):
    a=-100
    b=100
    m=1000
    dt=(b-a)*1.0/m
    Sa=0.5*E_applied(x-a)*plasmon(a)
    Sb=0.5*E_applied(x-b)*plasmon(b)
    Ssigma=0
    for i in range (1,m):
        t=a+i*dt
        Ssigma=Ssigma+E_applied(x-t)*plasmon(t)
    return (Sa+Sb+Ssigma)*dt
def TPPL_acr(delay):
    pow=math.pow
    
    data=[[],[],[]]
    a=-100
    b=100
    m=1000
    dt=(b-a)*1.0/m
    for i in range (-1*delay,delay):
        Sa=0.5*pow((E_plasmon(a)+E_plasmon(a+0.1*i)),4)
        Sb=0.5*pow((E_plasmon(b)+E_plasmon(b+0.1*i)),4)
        Ssigma=0
        for j in range (1,m):
            t=a+j*dt
            Ssigma=Ssigma+pow((E_plasmon(t)+E_plasmon(t+0.1*i)),4)
        data[0].append(i/10)
        data[1].append((Sa+Sb+Ssigma)*dt)
        data[2].append(E_plasmon(i))
    return data
def E_applied_acr(delay):
    pow=math.pow
    
    data=[[],[]]
    a=-1000
    b=1000
    m=10000
    dt=(b-a)*1.0/m
    for i in range (-1*delay,delay):
        Sa=0.5*pow((E_applied(a)+E_applied(a+0.1*i)),4)
        Sb=0.5*pow((E_applied(b)+E_applied(b+0.1*i)),4)
        Ssigma=0
        for j in range (1,m):
            t=a+j*dt
            Ssigma=Ssigma+pow((E_applied(t)+E_applied(t+0.1*i)),4)
        data[0].append(i/10)
        data[1].append((Sa+Sb+Ssigma)*dt)
    return data
    
file='C:/Users/1000297123/Desktop/acr_E790_18fs.txt'
file_dp='C:/Users/1000297123/Desktop/acr_E790_18fs_p790_t5.txt'
def plt_TPPL_acr_write(delay,file_adr):
    data=TPPL_acr(delay)
    t=np.array(data[0])
    Y=np.array(data[1])
    Conv=np.array(data[2])
    nor=max(Y)-min(Y)
    Y=8*(Y/nor)
    with open(file_adr,"w") as f:
        for i in range(0,2*delay):
            str1=str(t[i])+' '+str(Y[i])+' '+str(Conv[i])+'\n'
            f.write(str1) 
    
    plt.plot(t,Y)
    plt.plot(t,Conv)
    plt.xlim(-delay/10,delay/10)
    
def plt_shg_acr_write(delay,file_adr):
    data=E_applied_acr(delay)
    t=np.array(data[0])
    Y=np.array(data[1])
    nor=max(Y)-min(Y)
    Y=8*(Y/nor)
    with open(file_adr,"w") as f:
        for i in range(0,2*delay):
            str1=str(t[i])+' '+str(Y[i])+'\n'
            f.write(str1) 
    
    plt.plot(t,Y)
    plt.xlim(-100,100)       
    
plt_TPPL_acr_write(1000,file_dp) 
#plt_shg_acr_write(1000,file)

plt.show()
