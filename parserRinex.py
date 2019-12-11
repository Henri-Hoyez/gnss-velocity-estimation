import pandas as pd
import os
from utils.gps_time import *
import time as ti
from datetime import date, timedelta, time, datetime


#formatFichier : en paramètre le fichier .obs direct 
#                   la fonction va le mettre au bon format et le sauvegarder sous un autre nom (ici tamponObs.txt)
#           UTILISATION UNIQUE

def formatFichier(file):
    f=open("data/"+file)
    f2=open("tamponObs.txt","w")
    
    ligne=f.readline()
    ligneTab=[]
    tabFinal=[]
    n=len(ligne)-1
    j=0
    while not(ligne==""):
        ligne2=""    
        ligneTab=[]
        for i in range(len(ligne)):
            ligneTab.append(ligne[i])

        for i in reversed(range(len(ligneTab)-1)):
            if(ligneTab[i]==" " and ligneTab[i+1]==" "):
                ligneTab.pop(i)
            if(ligneTab[i]==" " and (ligneTab[i-1]=="G" or ligneTab[i-1]=="R")):
                ligneTab.pop(i)
        for k in range(len(ligneTab)):
            ligne2+=ligneTab[k]
        tabFinal.append(ligne2)
        if(j%100==0):
            print(tabFinal[j])
            print(j)
        j+=1
        ligne=f.readline() 
        f2.write(ligne2)
    f.close()
    f2.close()
    return f2.name

# obsToDataframe : paramètre le fichier tamponObs au bon format 
# retourne le dataframe avec comme index : 
def obsToDataframe(file):
    df=pd.read_csv("data/"+file,delim_whitespace=True,header=None,skiprows=lambda x:x<=29,names=['name','pseudo_range','carrier_phase','doppler','signal_strength','4','5'])
    df['second_of_week']=df['pseudo_range']
    i=0
    df['n_week']=get_gps_week()
    #index des lignes de temps
    indexRef=df[df['name']==">"].index
    val=get_second_of_week(datetime(int(df.loc[indexRef[0],'pseudo_range']),int(df.loc[indexRef[0],'carrier_phase']),int(df.loc[indexRef[0],'doppler']),int(df.loc[indexRef[0],'signal_strength']),int(df.loc[indexRef[0],'4']),int(df.loc[indexRef[0],'5'])))-1
    

    df.loc[0:indexRef[0]-1,'second_of_week']=val
    for i in range(len(indexRef)-1):
        val+=1
        df.loc[indexRef[i]:indexRef[i+1]-1,'second_of_week']=val
        i+=1
    df.loc[indexRef[-1]:,'second_of_week']=val+1

    del df['4'],df['5'] 

    df=df.drop(df[df.name==">"].index)
    df=df.set_index(['n_week','second_of_week','name'])
    return df


# Fonction finale à utiliser : 
#       Paramètre : nom du fichier .obs (exemple avec "autoroute_apres_midi.obs")
#
def obsToDataframeFinal(file):
    goodFile=formatFichier(file)
    df=obsToDataframe(goodFile)
    return df
# formatFichier("autoroute_apres_midi.obs")
# df=obsToDataframe("tamponObs.txt")
# print(df)
df=obsToDataframeFinal("autoroute_apres_midi.obs")
print(df)