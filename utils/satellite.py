import numpy as np
import datetime


class Satellite:

    def __init__(self,toe,t_data,mu0, delta_n, e, omega0, cws, cwc, sqrt_a, crc, crs, cic, cis, omega_ascension0, omega_dot_ascension, i0, i_dot):
        self.toe = toe
        self.t_data = t_data
        self.mu0 = mu0
        self.delta_n = delta_n
        self.e = e
        self.omega0 = omega0
        self.cws = cws
        self.cwc = cwc
        self.sqrt_a = sqrt_a
        self.crc = crc
        self.crs = crs
        self.cic = cic
        self.cis = cis
        self.omega_ascension0 = omega_ascension0
        self.omega_dot_ascension = omega_dot_ascension
        self.i0 = i0
        self.i_dot = i_dot



        self.name = None
        self.ephemeris_date = None
        self.type = None


    def set_type(self, type:str):
        self.type = type

    def set_name(self, name:str):
        self.name = name

    def set_ephemeris_date(self, ephemeris_date:datetime):
        self.ephemeris_date = ephemeris_date


    def get_position(self, t_obs = None):
        t_data = t_obs if t_obs != None else self.t_data

        a = self.sqrt_a**2

        GM = 3.986004418*10**14
        omega_e = 7.2921151467*10**-5

        t = t_data-self.toe
        mu  = self.mu0 + (np.sqrt(GM/(a**3)) + self.delta_n) * t
        # E = mu + e * np.sin(E)

        E = mu
        for i in range(7): 
            E = mu + self.e*np.sin(E)
        
        v = np.arctan((np.sqrt(1-self.e**2)*np.sin(E))/(np.cos(E)-self.e))

        # Correct for orbital perturbations
        omega = self.omega0 + self.cwc * np.cos(2*(self.omega0+v)) + self.cws * np.sin(2*(self.omega0+v))
        r = a*(1-self.e*np.cos(E))+ self.crc * np.cos(2*(self.omega0+v)) + self.crs*np.sin(2*(self.omega0+v))
        i = self.i0 + self.i_dot * t + self.cic * np.cos(2*(self.omega0+v)) + self.cis*np.sin(2*(self.omega0+v))
        
        # Compute  the  right  ascension
        omega_ascension = self.omega_ascension0 + (self.omega_dot_ascension-omega_e)*self.toe

        # Convert satellite position from orbital frame to ECEF frame
        r_vect = np.array([r*np.cos(v),r*np.sin(v),0])

        # Rotation matrix
        R = np.asarray([[np.cos(omega_ascension)*np.cos(omega)-np.sin(omega_ascension)*np.sin(omega)*np.cos(i), -np.cos(omega_ascension)*np.sin(omega)-np.sin(omega_ascension)*np.cos(omega)*np.cos(i), np.sin(omega_ascension)*np.sin(i)]
        ,[np.sin(omega_ascension)*np.cos(omega)+np.cos(omega_ascension)*np.sin(omega)*np.cos(i), -np.sin(omega_ascension)*np.sin(omega)+np.cos(omega_ascension)*np.cos(omega)*np.cos(i), -np.cos(omega_ascension)*np.sin(i)]
        ,[np.sin(omega)*np.sin(i), np.cos(omega)*np.sin(i), np.cos(i)]])
        
        return (-1)*np.dot(R,r_vect)

    def get_velocity(self,delta_t=1):
        pos1 = self.get_position()
        pos2 = self.get_position(self.t_data+delta_t)
        return (pos2-pos1)/(delta_t)
        

    
if __name__ == "__main__":
    sat = Satellite(208784, 208789, 0.128612756881*10,0.489020369661*10**-8,0.260362529662*10**-2,0.7931942133,0.622682273388*10**-5,0.1460313797*10**-5,0.515365547943*10**4,0.25971875*10**3,0.2853125*10**2,-0.204890966415*10**-7,0.577419996262*10**-7,0.108886242411*10,-0.815426822943*10**-8,0.963852456438,0.431089385175*10**-9)
    print(sat.get_position()/1000)
    print(np.linalg.norm(sat.get_velocity())*3.6)
    # pos2 = self.get_sat_pose(0.2160*10**6,216005,0.23388256954*10, 0.487698886045*10**-8, 0.260366254952*10**-2, 0.2081*10**4, 0.624172389507*10**-5, 0.157766044140*10**-5, 0.515365465546*10**4, 1, 0.31781250*10**2, -0.409781932831*10**-7, -0.2421438694*10**-7, 0.108880357969*10, -0.808605110220*10**-8, 0.963855990848, 0.468948104999*10**-9)
