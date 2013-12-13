import pickle
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams


def plotTESgraphs(tsave, Rsave, Isave, Tsave, Nphotons):

    rc('text',usetex=True)
    plt.figure(figsize=(12,8))
    plt.rcParams.update({'font.size':15})

    plt.hold(True)
    
    #axes fontsize, same font for math and text
    
    for i in range(Nphotons):
        plt.subplot(2,2,1)
        plt.grid(True)
        plt.plot(tsave[i]*10**6,Tsave[i]*10**3)
        plt.ylabel('Temperature ($mK$)', fontsize=15)
        plt.xlabel('time ($\mu s$)', fontsize=15)
        plt.xlim(135,190)
        plt.ylim(106.15,107.15)

        plt.subplot(2,2,2)
        plt.grid(True)
        plt.plot(tsave[i]*10**6,Isave[i]*10**6)
        plt.xlim(135,190)
        plt.ylim(0.4,0.9)
        plt.ylabel('Current ($\mu A$)', fontsize=15)
        plt.xlabel('time ($\mu s$)', fontsize=15)

        plt.subplot(2,2,3)
        plt.grid(True)
        plt.plot(tsave[i]*10**6,Rsave[i])
        plt.xlim(135,190)
        plt.ylim(.05,.25)
        plt.xlabel('time ($\mu s$)', fontsize=15)
        plt.ylabel('Resistance ($\Omega$)', fontsize=15)

        Psave = Isave[i]*Isave[i]*Rsave[i]
        plt.subplot(2,2,4)
        plt.grid(True)
        plt.plot(tsave[i]*10**6,Psave*10**15)
        plt.xlim(135,190)
        plt.ylim(0,150)
        plt.xlabel('time ($\mu s$)', fontsize=15)
        plt.ylabel('Bias Power ($fW$)', fontsize=15)

    plt.tight_layout()
    plt.suptitle('TES Sim Data', fontsize=20)
    plt.subplots_adjust(top=0.92)
    plt.show()


[tsave, Rsave, Isave, Tsave, Nphotons] = pickle.load(open('100photonsSim.pkl', 'r'))

plotTESgraphs(tsave, Rsave, Isave, Tsave, Nphotons)

