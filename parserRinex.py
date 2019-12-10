import pandas as pd
import os







def formatFichier(file):
    f=open(file)
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


def obsToDataframe(file):
    df=pd.read_csv(file,delim_whitespace=True,header=None,skiprows=lambda x:x<=29,names=['Sat','0','1','2','3','4','5'])
    print(df)
    i=0
    for index, row in df.iterrows():
        val=df.loc[index,'Sat']
    df['Time']=df['0']
    print(df)
# obsToDataframe("autoroute_apres_midi.obs")
# formatFichier("autoroute_apres_midi.obs")

obsToDataframe("tamponObs.txt")