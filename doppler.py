import numpy as np
# ftp://cddis.gsfc.nasa.gov/pub/gps/products/wwww/igswwwwd.sp3.Z | template ephemeris

class Doppler:
    # Compute basic parameters at request time
    def __init__(self):
        pass
    
    def __get_E_n(self,f_ti,Di,ri,ru,vi):
        c = 299792458
        norm_ri_ru = np.linalg.norm(ri-ru)
        K = -f_ti/(c*norm_ri_ru)
        R = vi[0]*(ri[0]-ru[0])+vi[1]*(-ri[1]-ru[1])
        return Di/K - R

    def get_usr_velocity(self):
        ru = []
        vi1 = []
        vi2 = []
        f_ti1 = 0
        f_ti2 = 0
        Di1 = 0
        Di2 = 0
        ri1 = []
        ri2 = [] 
        E1 = self.__get_E_n(f_ti1,Di1,ri1[:2],ru[:2],vi1[:2])
        E2 = self.__get_E_n(f_ti2,Di2,ri2[:2],ru[:2],vi2[:2])

        vuy = (E2-(E1*(-ri2[0]-ru[0]))/(-ri1[0]+ru[0]))/((ri2[0]-ru[0])/(ru[0]-ri1[0])+ru[1]-ri2[1])
        vux = (E1 - vuy*(-ri1[1]+ru[1]))/(-ri1[0]+ru[0])
        return (vux,vuy)



if __name__ == "__main__":
    # velocity = Doppler(208784,1574762406,1.28612756881,0.489020369661*10**-8,0.260362529662*10**-2,0.7931942133,0.622682263388*10**-5,0.1460313797*10**-5,0.515365547943*10**4,0.25971875*10**3,0.2853125*10**2,-0.204890966415*10**-7,0.577419996262*10**-7,0.108886242411*10,-0.815426822943*10**-8,0.963852456438,0.431089385164*10**-9)
    doppler = Doppler()





