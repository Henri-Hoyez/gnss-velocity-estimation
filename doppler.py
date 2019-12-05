import numpy as np
# ftp://cddis.gsfc.nasa.gov/pub/gps/products/wwww/igswwwwd.sp3.Z | template ephemeris

class Doppler:
    # Compute basic parameters at request time
    def __init__(self):
        pass
        

    def get_sat_pose(self,toe,t_data,mu0, delta_n, e, omega0, Cws, Cwc, a, Crc, Crs, Cic, Cis, omega_ascension0, omega_dot_ascension, i0, i_dot):
        a = a**2

        GM = 3.986004418*np.power(10,14)
        omega_e = 7.2921151467*10**-5

        t = t_data-toe
        mu  = mu0 + (np.sqrt(GM/np.power(a,3)) + delta_n) * t
        
        # E = mu + e * np.sin(E)
        E = np.arange(0,2*np.pi,0.01)
        sol = np.where(mu + e*np.sin(E) - E <= 0)[0][0]
        E= E[sol]
        E = np.roots([e/6,0,(1-e),-mu])
        E = np.real(np.select(np.imag(E)==0,E))
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
        

        print(np.dot(R,r_vect))
        return np.dot(R,r_vect)

    def get_sat_velocity(self,toe,t_data1,t_data2,mu0, delta_n, e, omega0, Cws, Cwc, a, Crc, Crs, Cic, Cis, omega_ascension0, omega_dot_ascension, i0, i_dot):
        pos1 = self.get_sat_pose(toe,t_data1,mu0, delta_n, e, omega0, Cws, Cwc, a, Crc, Crs, Cic, Cis, omega_ascension0, omega_dot_ascension, i0, i_dot)
        pos2 = self.get_sat_pose(toe,t_data2,mu0, delta_n, e, omega0, Cws, Cwc, a, Crc, Crs, Cic, Cis, omega_ascension0, omega_dot_ascension, i0, i_dot)
        print((np.linalg.norm(pos2-pos1)/(t_data2-t_data1))*3.6)
        # print((np.linalg.norm(pos2-pos1)))
        return np.linalg.norm(pos2-pos1)/(t_data2-t_data1)

if __name__ == "__main__":
    # velocity = Doppler(208784,1574762406,1.28612756881,0.489020369661*10**-8,0.260362529662*10**-2,0.7931942133,0.622682263388*10**-5,0.1460313797*10**-5,0.515365547943*10**4,0.25971875*10**3,0.2853125*10**2,-0.204890966415*10**-7,0.577419996262*10**-7,0.108886242411*10,-0.815426822943*10**-8,0.963852456438,0.431089385164*10**-9)
    doppler = Doppler()
    doppler.get_sat_velocity(208784,208789,215989,0.128612756881*10,0.489020369661*10**-8,0.260362529662*10**-2,0.7931942133,0.622682273388*10**-5,0.1460313797*10**-5,0.515365547943*10**4,0.25971875*10**3,0.2853125*10**2,-0.204890966415*10**-7,0.577419996262*10**-7,0.108886242411*10,-0.815426822943*10**-8,0.963852456438,0.431089385175*10**-9)





