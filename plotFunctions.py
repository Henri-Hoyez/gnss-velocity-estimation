from parserRinex import *
import numpy as np

def plotSignalStrength(fileName,satName,t1,t2):
    df=obsToDataframeFinal(fileName+".obs")


    satIsolated=df.query('name=="'+satName+'"')
    print(satIsolated)
    satTimeIsolated=satIsolated.query('(second_of_week >='+str(t1)+' )and (second_of_week <='+str(t2)+')')
    print(len(satTimeIsolated))
    time=np.linspace(t1,t2,len(satTimeIsolated))
    print(len(time))
    plt.plot(time,satTimeIsolated['signal_strength'])
    plt.axis(x=[0,5])
    plt.grid()
    plt.show()

def plotNumberOfSatellites(fileName,t1,t2):
    numberOfSat=[]
    df=obsToDataframeFinal(fileName+".obs")
    timeDF=df.query('(second_of_week >='+str(t1)+' )and (second_of_week <='+str(t2)+')')
    ctr=1
    i=1
    while(i<len(timeDF.index)-1):
        while(timeDF.index[i][1]==timeDF.index[i-1][1] and i<len(timeDF.index)-1):
            ctr+=1
            i+=1
        numberOfSat.append(ctr)
        ctr=1
        i+=1
    time=np.linspace(t1,t2,len(numberOfSat))
    plt.plot(time,numberOfSat)
    plt.grid()
    plt.show()


    
plotSignalStrength("autoroute_apres_midi","G2",131150,133400)
plotNumberOfSatellites("autoroute_apres_midi",131150,133200)