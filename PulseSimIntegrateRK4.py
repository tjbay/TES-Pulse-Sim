import random as rn
import math as math
import numpy as np
import matplotlib.pyplot as plt


def TES_resistance(T,I):
    return (Rn/2)*(1+math.tanh((T-Tc0+(I/A)**(2.0/3.0))/(dTc/(2))))

def TESdiffeq(T,I):        
    Rtes =  TES_resistance(T,I)    
    Tdot = (1/C)*(Rtes*(I**2) - kappa*(T**n-Tsub**n))
    Idot = -(1/L)*(I*(Rsh+Rp+Rtes)-Ib*Rsh)

    return [Tdot, Idot]   

def rungekutta4(T,I,dt):
    #calculates a single step of a runge-kutta 4th order solver

    #  k1 = dt * f(x)
    #  k2 = dt * f(x + k1/2)
    #  k3 = dt * f(x + k2/2)
    #  k4 = dt * f(x + k3)
    #  New X = x + (k0 + 2k1 + 2k2 + k3)/6
    
    k1_T, k1_I = TESdiffeq(T,I)
    k2_T, k2_I = TESdiffeq(T+dt*k1_T/2,I+dt*k1_I/2)
    k3_T, k3_I = TESdiffeq(T+dt*k2_T/2,I+dt*k2_I/2)
    k4_T, k4_I = TESdiffeq(T+dt*k3_T,I+dt*k3_I)

    newT = T + (k1_T*dt+2*k2_T*dt+2*k3_T*dt+k4_T*dt)/6
    newI = I + (k1_I*dt+2*k2_I*dt+2*k3_I*dt+k4_I*dt)/6
    
    return [newT, newI]

def addnoise():
    return rn.gauss(0,1)*1*10**(-9)



# constants
Tc0 = 110e-3
Tsub = 50e-3
dTc = 3e-3
k = 1.38e-23
h_bar = 1.054e-34
Rsh = 0.005
Rp = 0.002
Rn = 5.0
L = .35e-6
Sig = .4e-9
Vol = 24*24*0.035
gamma = 1.36e-16
C = gamma*Vol*Tc0*2.43

# values and initial conditions
count = 0
Ib = 20e-6
E0 = 3.1*1.602e-19
epsilon= 0.5
n=5.2
kappa = Vol*Sig
P0 = kappa*(Tc0**n-Tsub**n)
A = 3.52*math.sqrt(k*C/h_bar/Rn/Tc0)/1.5
Rbias = 0.018

I = Ib*Rsh/(Rbias+Rsh+Rp)
T = Tc0

dt = 1e-8
tstart = 0
tfinal = 500e-6
N = int((tfinal-tstart)/dt)


Tsave = np.zeros(N)
Rsave = np.zeros(N)
Isave = np.zeros(N)
tsave = np.zeros(N)


for i in range(N):
    tsave[i] = tstart + i*dt
    
    [T,I] = rungekutta4(T,I,dt)

    I += addnoise()

    if I < 0: I = 0

    if i == int(N/2):
        T = T + E0*epsilon/C

    Tsave[i] = T
    Isave[i] = I
    Rsave[i] = TES_resistance(T,I)

plt.figure(figsize=(12,8))
    
plt.subplot(2,2,1)
plt.plot(tsave*10**6,Tsave*10**3)
plt.xlim(200,300)
plt.ylim(105,109)
plt.ylabel('Temperature (mK)')

plt.subplot(2,2,2)
plt.plot(tsave*10**6,Isave*10**6)
plt.xlim(200,300)
plt.ylim(.4,1)
plt.ylabel('Current (microamps)')

plt.subplot(2,2,3)
plt.plot(tsave*10**6,Rsave)
plt.xlim(200,300)
plt.ylim(0,.3)
plt.xlabel('time (microseconds)')
plt.ylabel('Resistance (ohm)')

Psave = Isave*Isave*Rsave
plt.subplot(2,2,4)
plt.plot(tsave*10**6,Psave*10**15)
plt.xlim(200,300)
plt.ylim(0,175)
plt.xlabel('time (microseconds)')
plt.ylabel('Bias Power (femtowatts)')

plt.tight_layout()

plt.show()
