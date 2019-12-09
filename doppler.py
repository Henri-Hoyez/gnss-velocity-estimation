import numpy as np
# ftp://cddis.gsfc.nasa.gov/pub/gps/products/wwww/igswwwwd.sp3.Z | template ephemeris

class Doppler:
    # Compute basic parameters at request time
    def __init__(self):
        pass
        

    def get_sat_pose(self,toe,t_data,mu0, delta_n, e, omega0, Cws, Cwc, a, Crc, Crs, Cic, Cis, omega_ascension0, omega_dot_ascension, i0, i_dot):
        a = a**2

        GM = 3.986004418*10**14
        omega_e = 7.2921151467*10**-5

        t = t_data-toe
        mu  = mu0 + (np.sqrt(GM/(a**3)) + delta_n) * t
        # E = mu + e * np.sin(E)

        # E = np.arange(0,2*np.pi,0.01)
        # sol = np.where(mu + e*np.sin(E) - E <= 0)[0][0]
        # E= E[sol]
        E = mu
        for i in range(7): 
            E = mu + e*np.sin(E)
        
        v = np.arctan((np.sqrt(1-e**2)*np.sin(E))/(np.cos(E)-e))

        # Correct for orbital perturbations
        omega = omega0 + Cwc * np.cos(2*(omega0+v)) + Cws * np.sin(2*(omega0+v))
        r = a*(1-e*np.cos(E))+ Crc * np.cos(2*(omega0+v)) + Crs*np.sin(2*(omega0+v))
        i = i0 + i_dot * t + Cic * np.cos(2*(omega0+v)) + Cis*np.sin(2*(omega0+v))
        
        # Compute  the  right  ascension
        omega_ascension = omega_ascension0 + (omega_dot_ascension-omega_e)*toe

        # Convert satellite position from orbital frame to ECEF frame
        r_vect = np.array([r*np.cos(v),r*np.sin(v),0])

        # Rotation matrix
        R = np.asarray([[np.cos(omega_ascension)*np.cos(omega)-np.sin(omega_ascension)*np.sin(omega)*np.cos(i), -np.cos(omega_ascension)*np.sin(omega)-np.sin(omega_ascension)*np.cos(omega)*np.cos(i), np.sin(omega_ascension)*np.sin(i)]
        ,[np.sin(omega_ascension)*np.cos(omega)+np.cos(omega_ascension)*np.sin(omega)*np.cos(i), -np.sin(omega_ascension)*np.sin(omega)+np.cos(omega_ascension)*np.cos(omega)*np.cos(i), -np.cos(omega_ascension)*np.sin(i)]
        ,[np.sin(omega)*np.sin(i), np.cos(omega)*np.sin(i), np.cos(i)]])
        

        print(np.dot(R,r_vect)/1000)
        return np.dot(R,r_vect)

    def get_sat_velocity(self):
        pos1 = self.get_sat_pose(208784, 208789, 0.128612756881*10,0.489020369661*10**-8,0.260362529662*10**-2,0.7931942133,0.622682273388*10**-5,0.1460313797*10**-5,0.515365547943*10**4,0.25971875*10**3,0.2853125*10**2,-0.204890966415*10**-7,0.577419996262*10**-7,0.108886242411*10,-0.815426822943*10**-8,0.963852456438,0.431089385175*10**-9)
        pos2 = self.get_sat_pose(208784, 209684, 0.128612756881*10,0.489020369661*10**-8,0.260362529662*10**-2,0.7931942133,0.622682273388*10**-5,0.1460313797*10**-5,0.515365547943*10**4,0.25971875*10**3,0.2853125*10**2,-0.204890966415*10**-7,0.577419996262*10**-7,0.108886242411*10,-0.815426822943*10**-8,0.963852456438,0.431089385175*10**-9)
        # pos2 = self.get_sat_pose(0.2160*10**6,216005,0.23388256954*10, 0.487698886045*10**-8, 0.260366254952*10**-2, 0.2081*10**4, 0.624172389507*10**-5, 0.157766044140*10**-5, 0.515365465546*10**4, 1, 0.31781250*10**2, -0.409781932831*10**-7, -0.2421438694*10**-7, 0.108880357969*10, -0.808605110220*10**-8, 0.963855990848, 0.468948104999*10**-9)
        print((np.linalg.norm(pos2-pos1)/(209684-208789)))
        # return np.linalg.norm(pos2-pos1)/(t_data2-t_data1)

if __name__ == "__main__":
    # velocity = Doppler(208784,1574762406,1.28612756881,0.489020369661*10**-8,0.260362529662*10**-2,0.7931942133,0.622682263388*10**-5,0.1460313797*10**-5,0.515365547943*10**4,0.25971875*10**3,0.2853125*10**2,-0.204890966415*10**-7,0.577419996262*10**-7,0.108886242411*10,-0.815426822943*10**-8,0.963852456438,0.431089385164*10**-9)
    doppler = Doppler()
    doppler.get_sat_velocity()





