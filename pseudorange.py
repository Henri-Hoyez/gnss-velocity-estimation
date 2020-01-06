import numpy as np
from utils.nav_parser import parse_nav_file

# ftp://cddis.gsfc.nasa.gov/pub/gps/products/wwww/igswwwwd.sp3.Z | template ephemeris

class Doppler:
    # Compute basic parameters at request time
    def __init__(self):
        sats = parse_nav_file("data/Very_Bad_Trip/Belgique/autoroute_plus_tunnel.nav")
        print([sat.name for sat in sats])

        self.ri1 = sats[10].get_position()
        self.ri2 = sats[9].get_position()
        self.ri3 = sats[0].get_position()
        print(self.ri1)
        # print(self.ri2)
        self.vi1 = sats[10].get_velocity()
        self.vi2 = sats[9].get_velocity()
        self.vi3 = sats[0].get_velocity()

        self.cd1 = sats[10].clock_drift
        self.cd2 = sats[9].clock_drift
        self.cd3 = sats[0].clock_drift
        

    
    def __get_K_n(self,rho,clock_drift,ri,ru,vi):
        c = 299792458
        a = self.get_line_of_sight(ri,ru)
        print(a)
        return -rho - c*clock_drift + np.vdot(vi,a)

    def get_line_of_sight(self,ri,ru):
        return (ri-ru)/np.linalg.norm(ri-ru)

    def get_usr_velocity(self):
        #G06 G23 G03 // 10 & 9 & 0
        # origin obs : 4049390.9248   215421.6043  4906680.8162
        c = 299792458
        f_ti = 1575.42*10**6
        ru = np.array([4043730.5731   ,261041.2499  , 4909152.9334])*0
        rho1 = -c*1540.719/f_ti
        rho2 = -c*-207.044/f_ti
        rho3 = -c*386.173/f_ti

        k1 = self.__get_K_n(rho1,self.cd1,self.ri1,ru,self.vi1)
        k2 = self.__get_K_n(rho2,self.cd2,self.ri2,ru,self.vi2)
        k3 = self.__get_K_n(rho3,self.cd3,self.ri3,ru,self.vi3)

        a1 = self.get_line_of_sight(self.ri1,ru)
        a2 = self.get_line_of_sight(self.ri2,ru)
        a3 = self.get_line_of_sight(self.ri3,ru)

        K = np.array([k1,k2,k3])
        X = np.array([a1,a2,a3])
        v = np.linalg.lstsq(X,K,rcond=-1)
        # print(np.linalg.inv(X))
        print(K)

        print("////",np.linalg.norm(v[0]))
        
        return v



if __name__ == "__main__":
    sats = parse_nav_file("data/Very_Bad_Trip/Belgique/autoroute_plus_tunnel.nav")
    # print(sats[1].name)
    # print(sats[1].get_position()/1000)
    # print(np.linalg.norm(sats[1].get_velocity())*3.6)
    # velocity = Doppler(208784,1574762406,1.28612756881,0.489020369661*10**-8,0.260362529662*10**-2,0.7931942133,0.622682273388*10**-5,0.1460313797*10**-5,0.515365547943*10**4,0.25971875*10**3,0.2853125*10**2,-0.204890966415*10**-7,0.577419996262*10**-7,0.108886242411*10,-0.815426822943*10**-8,0.963852456438,0.431089385164*10**-9)
    doppler = Doppler()
    doppler.get_usr_velocity()
