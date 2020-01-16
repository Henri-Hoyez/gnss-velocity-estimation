from parserRinex import *
import numpy as np
from utils import nav_parser
import math

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

def getNumberOfSatellites(fileName,t1,t2):
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
    return numberOfSat



def getSatOnTime(fileName,t1):
    df=obsToDataframeFinal(fileName+".obs")
    timeDF=df.query('second_of_week =='+str(t1))
    ctr=1
    i=1
    while(timeDF.index[i][1]==timeDF.index[i-1][1] and i<len(timeDF.index)-1):
        ctr+=1
        i+=1
    numberOfSat=ctr
    return numberOfSat

def plotNumberOfSatellites(fileName,t1,t2):
    numberOfSat=getNumberOfSatellites(fileName,t1,t2)
    time=np.linspace(t1,t2,len(numberOfSat))
    plt.plot(time,numberOfSat)
    plt.grid()




def plotDOP(fileName,t1,t2):
    numberOfSat=getNumberOfSatellites(fileName,t1,t2)
    df=obsToDataframeFinal(fileName+".obs")
    df2=posToDataframe(fileName+".pos")
    dopTab=[]
    time=[]
    timeDF=df.query('(second_of_week >='+str(t1)+' )and (second_of_week <='+str(t2)+')')
    posDF=df2.query('(second_of_week >='+str(t1)+' )and (second_of_week <='+str(t2)+')')
    print(posDF,timeDF)
    i=1
    while(i<len(posDF)):
        print(i)
        print(posDF.index[i][1])
        time.append(int(posDF.index[i][1]))
        dopTab.append(calculateDOP(fileName,int(posDF.index[i][1])))
        i+=1
    
    
    plt.plot(time,dopTab)
    plt.grid()
    




def calculateDOP(fileName,t_obs):
    numberOfSat=getSatOnTime(fileName,t_obs)
    print(numberOfSat)
    df=obsToDataframeFinal(fileName+".obs")
    df2=posToDataframe(fileName+".pos")
    timeDF=df.query('second_of_week =='+str(t_obs))
    posDF=df2.query('second_of_week =='+str(t_obs))
    #print(timeDF)
    posX=posDF['X']
    posY=posDF['Y']
    posZ=posDF['Z']
    for j in range(len(timeDF.index)):
         if(not('G' in timeDF.index[j][2])):
             numberOfSat-=1

    matrixA=np.zeros((numberOfSat,4))
    print(matrixA.shape)
    for j in range(len(matrixA)):
        pseudoRange=timeDF.loc[(df.index[j],t_obs,timeDF.index[j]),'pseudo_range']
        if('G' in timeDF.index[j][2]):
            satPos=getSatPos(fileName,timeDF.index[j][2],t_obs)
        matrixA[j][0]=(satPos[0] - posX)/pseudoRange
        matrixA[j][1]=(satPos[1] - posY)/pseudoRange
        matrixA[j][2]=(satPos[2] - posZ)/pseudoRange
        matrixA[j][3]=-1
    
    transposeeA=matrixA.transpose()
    tmpA=np.dot(transposeeA,matrixA)
    if(not(np.linalg.det(tmpA)==0)):
        matrixQ=np.linalg.inv(tmpA)
        PDOP=math.sqrt(matrixQ[0][0]**2+matrixQ[1][1]**2+matrixQ[2][2]**2+matrixQ[3][3]**2)
        TDOP=math.sqrt(matrixQ[-1][-1]**2)
        GDOP=math.sqrt(PDOP**2+TDOP**2)
        return GDOP
    else:
        return 0


    PDOP=math.sqrt(matrixQ[0][0]**2+matrixQ[1][1]**2+matrixQ[2][2]**2+matrixQ[3][3]**2)
    TDOP=math.sqrt(matrixQ[-1][-1]**2)
    GDOP=math.sqrt(PDOP**2+TDOP**2)
    print("PDOP : ",PDOP)
    print("TDOP : ",TDOP)
    print("GDOP : ",GDOP)
def getSatPos(fileName,satName,t_obs):
    satelites=nav_parser.parse_nav_file('data/'+fileName+".nav")
    for i in satelites:
        nSat=i.name.replace(" ","")
        
        if(nSat==satName):
            return i.get_position(t_obs=t_obs)
plt.subplot(3,1,1)
plt.title("Signal Strength")
for i in range(30):
    plotSignalStrength("autoroute_apres_midi","G"+str(i),131700,132000)
plt.subplot(3,1,2)
plt.title("Number of satellites")
plotNumberOfSatellites("autoroute_apres_midi",131700,132000)

# satelites = nav_parser.parse_nav_file('data/autoroute_apres_midi.nav')
# print(satelites[0].name)
# print(satelites[0].get_position(t_obs=131120))
# getSatOnTime("autoroute_apres_midi",131120)
#calculateDOP("autoroute_apres_midi",131121)
plt.subplot(3,1,3)
plt.title("DOP")
plotDOP("autoroute_apres_midi",131700,132000)
plt.show()