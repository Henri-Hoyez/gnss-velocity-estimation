import numpy as np
from utils.nav_parser import parse_nav_file

# ftp://cddis.gsfc.nasa.gov/pub/gps/products/wwww/igswwwwd.sp3.Z | template ephemeris

class Doppler:
    # Compute basic parameters at request time
    def __init__(self):
        sats = parse_nav_file("test_data/autoroute_plus_tunnel.nav")

        # set the good t_data ( 10 dec à 10:00:00 )

        sats[9].t_data = 208800.0
        sats[0].t_data = 208800.0

        print(sats[10].name,   ' - ', sats[9].name,   ' - ',sats[0].name )
        print(sats[10].t_data, ' - ', sats[9].t_data, ' - ',sats[0].t_data )

        self.ri1 = sats[10].get_position()
        self.ri2 = sats[9].get_position()
        self.ri3 = sats[0].get_position()
        print(self.ri1)
        #print(self.ri2)

        self.vi1 = sats[10].get_velocity()
        self.vi2 = sats[9].get_velocity()
        self.vi3 = sats[0].get_velocity()
        #print(self.vi2)
        







    def test_angles(self):
        
        ri12 = np.arccos(np.dot(self.ri1, self.ri2) / (np.linalg.norm(self.ri1) * np.linalg.norm(self.ri2)) )
        ri13 = np.arccos(np.dot(self.ri1, self.ri3) / (np.linalg.norm(self.ri1) * np.linalg.norm(self.ri3)) )
        ri23 = np.arccos(np.dot(self.ri2, self.ri3) / (np.linalg.norm(self.ri2) * np.linalg.norm(self.ri3)) )  

        print(ri12, ri13, ri23)


    def corelation_coeficiant(self, v1:np.ndarray, V2:np.ndarray):
        pass

    def get_k(self, D_i, f_ti, v_i, a_i):
        c = 3.00*10**8
        res = c*D_i/f_ti + np.dot(v_i, a_i)
        return res


    def azimuth_to_ECEF(self, sat_position:np.ndarray, user_position:np.ndarray):
        # sat_position  => [ x_sat , y_sat , z_sat]
        # user_position => [ x_usr , y_usr , z_usr]

        vect_norm = (sat_position - user_position) / np.linalg.norm(sat_position - user_position)
        # [vect_normx, vect_normy, vect_normz]


        theta = np.arctan( (vect_norm[0]) / ( vect_norm[2] ) )
        phi = np.arcsin( (vect_norm[2]) / ( np.sqrt(vect_norm[0]**2 + vect_norm[2]**2) ) )
        x, y, z = np.cos(phi) * np.cos(theta), np.cos(phi) * np.sin(theta), np.sin(phi) 

        return np.array([x, y, z])


    def get_usr_velocity(self):
        #G6 et G23 // 10 & 9

        ru = np.array([4043547.78553915, 254207.686387644, 4909623.02474359])

        f_ti1 = 1575.42*10**6 #105690384.812
        f_ti2 = 1575.42*10**6 #101943589.013
        f_ti3 = 1575.42*10**6

        Di1 = 1319.955
        Di2 = -513.404
        Di3 = -2687.413

        # convert azimits vectors to ECEF-translated 'base' beceause all velocity 
        # are axplained in this azimut base.

        a1 = self.azimuth_to_ECEF(self.ri1, ru)
        a2 = self.azimuth_to_ECEF(self.ri2, ru)
        a3 = self.azimuth_to_ECEF(self.ri3, ru)

        k1 = self.get_k(Di1, f_ti1, self.vi1, a1)
        k2 = self.get_k(Di2, f_ti2, self.vi2, a2)
        k3 = self.get_k(Di3, f_ti3, self.vi3, a3)

        K = np.array([k1, k2, k3]).T

        M = np.array([a1, a2, a3])

        v = np.dot(np.linalg.inv(M), K)
        print()
        print( np.linalg.inv(M) )
        print( K )

        print(np.linalg.norm(v) * 3.6 , 'km/h paupiette')


        



        



if __name__ == "__main__":
    sats = parse_nav_file("test_data/autoroute_plus_tunnel.nav")
    # print(sats[1].name)
    # print(sats[1].get_position()/1000)
    # print(np.linalg.norm(sats[1].get_velocity())*3.6)
    # velocity = Doppler(208784,1574762406,1.28612756881,0.489020369661*10**-8,0.260362529662*10**-2,0.7931942133,0.622682273388*10**-5,0.1460313797*10**-5,0.515365547943*10**4,0.25971875*10**3,0.2853125*10**2,-0.204890966415*10**-7,0.577419996262*10**-7,0.108886242411*10,-0.815426822943*10**-8,0.963852456438,0.431089385164*10**-9)
    doppler = Doppler()
    doppler.get_usr_velocity()






