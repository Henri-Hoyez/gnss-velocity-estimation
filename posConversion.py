import os 
import numpy as np
import math
import matplotlib.pyplot as plt

global windowSize
windowSize=10

def nmeaSpeed(positionFile):
    data=[]
    file=open(positionFile)
    ligne=file.readline()
    while(ligne!=""):
        if("$GPRMC" in ligne):
            dataLine=ligne.split(",")
            data.append(dataLine)
        ligne=file.readline()

    vitesse=[]
    for i in range(len(data)):
        vitesse.append(float(data[i][7])*1.852)
    return vitesse
    
def tableauData(positionFile):
    data=[]
    file=open(positionFile)
    compteur=0
    ligne=file.readline()

    while("%" in ligne):
        ligne=file.readline()
        compteur+=1

    ligne=file.readline()
    while(ligne !=""):
        dataLine=ligne.split("   ",4)
        data.append(dataLine[1:4])
        ligne=file.readline()
    print(len(data))
    return data


def velocityEstimation(positionFile):
    #Fonction qui extrait les données brutes du fichier .pos
    data=tableauData(positionFile)

    #Temps de 0 jusqu'au nombre de secondes d'enregistrement
    time=[]
    for i in range(len(data)):
        time.append(i)

    velocityWindow=[]
    velocityNorm=[]
    positionTime=[]
    velocityXmean=[]
    velocityYmean=[]
    velocityZmean=[]
    posX=[]
    posY=[]
    posZ=[]
    #Ajout du temps dans le tableau des positions
    for i in range(len(data)):
        positionTime.append([data[i],time[i]])
        
    #Calcul des vitesses sur une fenêtre de x secondes
    for i in range(len(positionTime)-windowSize):
        velocityX=[]
        velocityY=[]
        velocityZ=[]
        for j in range(i,i+windowSize):
            velocityX.append((float(positionTime[j+1][0][0])-float(positionTime[j][0][0]))/(float(positionTime[j+1][1])-float(positionTime[j][1])))
            velocityY.append((float(positionTime[j+1][0][1])-float(positionTime[j][0][1]))/(float(positionTime[j+1][1])-float(positionTime[j][1])))
            velocityZ.append((float(positionTime[j+1][0][2])-float(positionTime[j][0][2]))/(float(positionTime[j+1][1])-float(positionTime[j][1])))

        posX.append(float(positionTime[i][0][0]))
        posY.append(float(positionTime[i][0][1]))
        posZ.append(float(positionTime[i][0][2]))

        velocityXmean.append(sum(velocityX)/len(velocityX))
        velocityYmean.append(sum(velocityY)/len(velocityY))
        velocityZmean.append(sum(velocityZ)/len(velocityZ))

        print(velocityXmean[-1])
    for i in range(len(velocityXmean)):
        velocityNorm.append(math.sqrt(velocityXmean[i]**2+velocityYmean[i]**2+velocityZmean[i]**2)*3.6)
    return velocityNorm


myVel=velocityEstimation("autoroute_apres_midi.pos")
nmeaVel=nmeaSpeed("autoroute_apres_midi.nmea")
# plt.subplot(3,1,1)
# plt.plot(time2,velocityXmean)
# plt.subplot(3,1,2)
# plt.plot(time2,velocityYmean)
# plt.subplot(3,1,3)
# plt.plot(time2,velocityZmean)

# plt.subplot(2,1,1)
# plt.plot(posX,posY)

# plt.subplot(2,1,2)
time=np.linspace(0,2000,len(myVel))
time2=np.linspace(0,2000,len(nmeaVel))
plt.title("Vitesse mesurée sur l'autoroute")
plt.plot(time,myVel,label="Différence de position")
plt.plot(time2,nmeaVel,label="Ublox")
plt.xlabel("Time (s)")
plt.ylabel("Speed (km/h)")
plt.legend(prop={'size':12})
plt.grid()
plt.show()
