import os 
import numpy as np
import math
import matplotlib.pyplot as plt

def tableauData(positionFile):
    data=[]
    file=open(positionFile)
    compteur=0
    ligne=file.readline()
    while(ligne!=""):
        if("$GPRMC" in ligne):
            dataLine=ligne.split(",")
            data.append(dataLine)
        ligne=file.readline()

    vitesse=[]
    for i in range(len(data)):
        vitesse.append(float(data[i][7])*1.852)
    time=np.linspace(0,len(vitesse),len(vitesse))
    print(vitesse[0:10])
    plt.plot(time,vitesse)
    plt.grid()
    plt.show()
tableauData("autoroute_apres_midi.nmea")
