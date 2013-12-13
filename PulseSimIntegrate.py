import random as rn
import math as math
import numpy as np
import matplotlib.pyplot as plt
#import pickle

#  Calculates the TES resistance as a function of Temperature and Current
def TES_resistance(T,I):
    return (Rn/2)*(1+math.tanh((T-Tc0+(I/A)**(2.0/3.0))/(dTc/(2))))

#  The TES differential equations
def TESdiffeq(T,I):        
    Rtes =  TES_resistance(T,I)    
    Tdot = (1/C)*(Rtes*(I**2) - kappa*(T**n-Tsub**n))
    Idot = -(1/L)*(I*(Rsh+Rp+Rtes)-Ib*Rsh)
    return [Tdot, Idot]   
  
# calculates a single step of a runge-kutta 4th order solver

def rungekutta4(T,I,dt):    

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

# Simulates Johnson noise
def addnoise():
    return rn.gauss(0,1)*1*10**(-9)

#  Plots Graphs
def plotTESgraphs(tsave, Rsave, Isave, Tsave, Nphotons):

    plt.figure(figsize=(9,6))

    plt.hold(True)
    
    for i in range(Nphotons):
        plt.grid(True)
        plt.plot(tsave[i]*10**6,Isave[i]*10**6)
        plt.xlim(135,190)
        plt.ylim(0.4,0.9)
        plt.ylabel('Current ($\mu A$)', fontsize=15)
        plt.xlabel('time ($\mu s$)', fontsize=15)

    plt.tight_layout()
    plt.show()

# Constants
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

# Values and initial conditions
count = 0
Ib = 20e-6
E0 = 3.1*1.602e-19
epsilon= 0.5
n=5.2
kappa = Vol*Sig
P0 = kappa*(Tc0**n-Tsub**n)
A = 3.52*math.sqrt(k*C/h_bar/Rn/Tc0)/1.5
Rbias = 0.018
dt = 1e-8
tstart = 0
tfinal = 200e-6
steps = int((tfinal-tstart)/dt)
Nphotons = 5

# Initialize arrays to hold simulation data

Tsave = np.zeros((Nphotons, steps))
Rsave = np.zeros((Nphotons, steps))
Isave = np.zeros((Nphotons, steps))
tsave = np.zeros((Nphotons, steps))

# Set first element of arrays to initial conditions

Tsave[:,0] = Tc0
Isave[:,0] = Ib*Rsh/(Rbias+Rsh+Rp)
Rsave[:,0] = TES_resistance(Tc0, Ib*Rsh/(Rbias+Rsh+Rp))

#  Loop to calculate time steps

for i in range(steps-1):

    time = tstart+i*dt

    if i%2000 == 0: print (time*10**6)
    
    for j in range(Nphotons):

        T = Tsave[j,i]
        I = Isave[j,i]
        
        [newT,newI] = rungekutta4(T,I,dt)

    # Add some random noise.  Magnitude set to "look" right, not match 
    # any physical sources, such as Johnson noise.
        newI += addnoise()
        if newI < 0: newI = 0
        if i == int(steps/1.33):
            newT += E0*epsilon/C
        
        Tsave[j, i+1] = newT
        Isave[j, i+1] = newI
        Rsave[j, i+1] = TES_resistance(newT,newI)
        tsave[j, i+1] = time

plotTESgraphs(tsave, Rsave, Isave, Tsave, Nphotons)
    
#pickle.dump([tsave, Rsave, Isave, Tsave, Nphotons], open('300photonsSim.pkl', 'w'))

