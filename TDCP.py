import numpy as np
from utils.nav_parser import parse_nav_file

# ftp://cddis.gsfc.nasa.gov/pub/gps/products/wwww/igswwwwd.sp3.Z | template ephemeris

class TDCP:
    # Compute basic parameters at request time
    def __init__(self):
        sats = parse_nav_file("data/Very_Bad_Trip/Belgique/autoroute_plus_tunnel.nav")
        print([sat.name for sat in sats])
        print("blabla",-sats[0].get_position(sats[0].toe+1)/1000)
        self.x_sk = np.array([[sats[10].get_position(),sats[9].get_position(),-sats[0].get_position()],
        [sats[10].get_position(sats[10].toe+1),sats[9].get_position(sats[9].toe+1),-sats[0].get_position(sats[0].toe+1)]])
        self.x_rk = np.array([[4043743.6490 ,   261011.8175  , 4909156.8423],
        [4043730.5731  ,  261041.2499  , 4909152.9334]])
        print("groundtruth", np.linalg.norm(self.x_rk[1]-self.x_rk[0])*3.6)
        
        
    def get_phase_difference(self,phase1,phase2):
        return phase2 - phase1
        
    
    def __get_K_n(self,wavelength,phase1,phase2,delta_S,delta_G,delta_t):
        print("**** deltaS : ",delta_S ," deltaG : ", delta_G ," phase diff : ", wavelength * self.get_phase_difference(phase1,phase2)," result : ",delta_S - delta_G - wavelength * self.get_phase_difference(phase1,phase2))
        return (delta_S - delta_G - wavelength * self.get_phase_difference(phase1,phase2)) / (delta_t)

    def get_line_of_sight(self,ri,ru):
        return (ri-ru)/np.linalg.norm(ri-ru)


    def get_usr_velocity(self):
        #G06 G23 G03 // 10 & 9 & 0
        # origin obs : 4049390.9248   215421.6043  4906680.8162
        # c = 299792458
        # f_ti = 1575.42*10**6
        # ru = np.array([4043730.5731   ,261041.2499  , 4909152.9334])*0
        # rho1 = -c*1540.719/f_ti
        # rho2 = -c*-207.044/f_ti
        # rho3 = -c*386.173/f_ti
        phase1 = np.array([106543859.919,101729930.484,103940539.651])
        phase2 = np.array([106542319.194,101730136.853,103942948.132])

        c = 299792458
        f = 1575.42*10**6
        wavelength = c/f

        delta_S1 = np.inner(self.get_line_of_sight(self.x_sk[1,0],self.x_rk[1]),self.x_sk[1,0]) - np.inner(self.get_line_of_sight(self.x_sk[0,0],self.x_rk[0]),self.x_sk[0,0]) 
        delta_S2 = np.inner(self.get_line_of_sight(self.x_sk[1,1],self.x_rk[1]),self.x_sk[1,1]) - np.inner(self.get_line_of_sight(self.x_sk[0,1],self.x_rk[0]),self.x_sk[0,1])
        delta_S3 = np.inner(self.get_line_of_sight(self.x_sk[1,2],self.x_rk[1]),self.x_sk[1,2]) - np.inner(self.get_line_of_sight(self.x_sk[0,2],self.x_rk[0]),self.x_sk[0,2])
        delta_G1 = np.inner(self.get_line_of_sight(self.x_sk[1,0],self.x_rk[1]),self.x_rk[0]) -   np.inner(self.get_line_of_sight(self.x_sk[0,0],self.x_rk[0]),self.x_rk[0])
        delta_G2 = np.inner(self.get_line_of_sight(self.x_sk[1,1],self.x_rk[1]),self.x_rk[0]) -   np.inner(self.get_line_of_sight(self.x_sk[0,1],self.x_rk[0]),self.x_rk[0])
        delta_G3 = np.inner(self.get_line_of_sight(self.x_sk[1,2],self.x_rk[1]),self.x_rk[0]) -   np.inner(self.get_line_of_sight(self.x_sk[0,2],self.x_rk[0]),self.x_rk[0])
        
         
        k1 = self.__get_K_n(wavelength,phase1[0],phase2[0],delta_S1,delta_G1,1)
        k2 = self.__get_K_n(wavelength,phase1[1],phase2[1],delta_S2,delta_G2,1)
        k3 = self.__get_K_n(wavelength,phase1[2],phase2[2],delta_S3,delta_G3,1)

        #k1 = 25
        #k2 = 25
        #k3 = 25

        a1 = self.get_line_of_sight(self.x_sk[1,0],self.x_rk[1])
        a2 = self.get_line_of_sight(self.x_sk[1,1],self.x_rk[1])
        a3 = self.get_line_of_sight(self.x_sk[1,2],self.x_rk[1])

        K = np.array([k1,k2,k3])
        X = np.array([a1,a2,a3])
        v = np.linalg.lstsq(X,K,rcond=-1)
        # print(np.linalg.inv(X))
        print(X)
        # l = np.pi/6
        # phi = np.pi/3
        # to_ENU = np.array([[-np.sin(l),np.cos(l),0]
        # ,[-np.cos(l)*np.sin(phi),-np.sin(l)*np.sin(phi),np.cos(phi)]
        # ,[np.cos(l)*np.cos(phi),np.sin(l)*np.cos(phi),np.sin(phi)]])
        # v = np.dot(to_ENU,v)
        print("////",np.linalg.norm(v[0])*3.6)

        return v



if __name__ == "__main__":
    sats = parse_nav_file("data/Very_Bad_Trip/Belgique/autoroute_plus_tunnel.nav")
    # print(sats[1].name)
    # print(sats[1].get_position()/1000)
    # print(np.linalg.norm(sats[1].get_velocity())*3.6)
    # velocity = Doppler(208784,1574762406,1.28612756881,0.489020369661*10**-8,0.260362529662*10**-2,0.7931942133,0.622682273388*10**-5,0.1460313797*10**-5,0.515365547943*10**4,0.25971875*10**3,0.2853125*10**2,-0.204890966415*10**-7,0.577419996262*10**-7,0.108886242411*10,-0.815426822943*10**-8,0.963852456438,0.431089385164*10**-9)
    tdcp = TDCP()
    tdcp.get_usr_velocity()
