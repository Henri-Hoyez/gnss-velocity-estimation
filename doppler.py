import numpy as np
from utils.nav_parser import parse_nav_file

# ftp://cddis.gsfc.nasa.gov/pub/gps/products/wwww/igswwwwd.sp3.Z | template ephemeris

class Doppler:
    # Compute basic parameters at request time
    def __init__(self):
        sats = parse_nav_file("data/Very_Bad_Trip/Belgique/autoroute_plus_tunnel.nav")

        self.ri1 = sats[10].get_position()
        self.ri2 = sats[9].get_position()
        print(self.ri1/1000)
        print(self.ri2/1000)
        self.vi1 = sats[10].get_velocity()
        self.vi2 = sats[9].get_velocity()
        

    
    def __get_K_n(self,f_ti,Di,ri,ru,vi):
        c = 299792458
        norm_ri_ru = np.linalg.norm(ri-ru)
        return (c*Di*norm_ri_ru)/f_ti + vi[0]*(ri[0]-ru[0]) + vi[1]*(ri[1]-ru[1])

    def get_usr_velocity(self):
        #G6 et G23 // 10 & 9
        ru = [4043743.6490,261011.8175]
        f_ti1 = 1575.42*10**6 #105690384.812
        f_ti2 = 1575.42*10**6 #101943589.013
        Di1 = 1319.955
        Di2 = -513.404
        K1 = self.__get_K_n(f_ti1,Di1,self.ri1[:2],ru[:2],self.vi1[:2])
        K2 = self.__get_K_n(f_ti2,Di2,self.ri2[:2],ru[:2],self.vi2[:2])
        a1 = self.ri1[0]-ru[0]
        a2 = self.ri2[0]-ru[0]
        b1 = self.ri1[1]-ru[1]
        b2 = self.ri2[1]-ru[1]

        K = np.transpose([K1,K2])
        matrix_to_inv = np.array([[a1,b1],[a2,b2]])
        m_inv = np.linalg.inv(matrix_to_inv)
        v = np.dot(m_inv,K)
        print(v)
        print(np.linalg.norm(v)*3.6)
        # vuy = (E2-(E1*(-self.ri2[0]-ru[0]))/(-self.ri1[0]+ru[0]))/((self.ri2[0]-ru[0])/(ru[0]-self.ri1[0])+ru[1]-self.ri2[1])
        # vux = (E1 - vuy*(-self.ri1[1]+ru[1]))/(-self.ri1[0]+ru[0])
        # print(np.linalg.norm([vux,vuy])*3.6)
        return v



if __name__ == "__main__":
    sats = parse_nav_file("data/Very_Bad_Trip/Belgique/autoroute_plus_tunnel.nav")
    # print(sats[1].name)
    # print(sats[1].get_position()/1000)
    # print(np.linalg.norm(sats[1].get_velocity())*3.6)
    # velocity = Doppler(208784,1574762406,1.28612756881,0.489020369661*10**-8,0.260362529662*10**-2,0.7931942133,0.622682273388*10**-5,0.1460313797*10**-5,0.515365547943*10**4,0.25971875*10**3,0.2853125*10**2,-0.204890966415*10**-7,0.577419996262*10**-7,0.108886242411*10,-0.815426822943*10**-8,0.963852456438,0.431089385164*10**-9)
    doppler = Doppler()
    doppler.get_usr_velocity()






