import numpy as np
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Satellite:

    def __init__(self,toe,t_data,mu0, delta_n, e, omega0, cws, cwc, sqrt_a, crc, crs, cic, cis, omega_ascension0, omega_dot_ascension, i0, i_dot, clock_drift):
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
        self.clock_drift = clock_drift  


        self.name = None
        self.ephemeris_date = None
        self.type = None
        self.observation_error = None


    def set_type(self, type:str):
        self.type = type

    def set_name(self, name:str):
        self.name = name

    def set_ephemeris_date(self, ephemeris_date:datetime):
        self.ephemeris_date = ephemeris_date

    def set_observation_error(self, observation_error):
        self.observation_error = observation_error


    def show_trajetcory(self):
        t_in  = self.toe
        t_out = self.toe + 3*60**2

        positions = []

        for t in range(int(t_in), int(t_out)):
            positions.append(self.get_pos(t)[:3])

        positions = np.array(positions)
        fig = plt.figure()
        ax = Axes3D(fig)

        ax.scatter(positions[:,0],positions[:,1],positions[:,2])

        plt.show()


    def get_pos(self, t_obs:int=None):
        t_data = t_obs if t_obs != None else self.t_data

        bOMEGAE84 = 7.2921151467*10**-5
        bGM84 = 3.986005e14

        # self.sqrt_a     =  .515365547943e+04
        # self.toe        =  .208784000000e+06
        # self.mu0        =  .128612756881e+01
        # self.e          =  .260362529662E-02
        # self.delta_n    =  .489020369661e-08
        # self.omega0     =  .793194213300e+00
        # self.cws        =  .622682273388E-05
        # self.cwc        =  .146031379700E-05
        # self.crs        =  .285312500000E+02
        # self.crc        =  .285312500000E+03
        # self.cis        =  .577419996262E-07
        # self.cic        =  -.204890966415e-07
        # self.i_dot      =  .431089385175e-09
        # self.i0         =  .963852456438E+00
        # self.omega_ascension0  =  .108886242411e+01
        # self.omega_dot_ascension= -.815426822943E-08
        # t_data = 381600

        A = self.sqrt_a * self.sqrt_a
        
        n0 = np.sqrt(bGM84/(A*A*A))
        tk = t_data - self.toe
        n = n0 + self.delta_n
        mk = self.mu0 + n * tk

        mkdot = n
        ek = mk

        for _ in range(7):
            ek = mk + self.e*np.sin(ek)
        
        ekdot = mkdot/(1.0 - self.e*np.cos(ek))
        tak = np.arctan2( np.sqrt(1.0 - self.e * self.e)*np.sin(ek), np.cos(ek)-self.e)

        takdot = np.sin(ek)*ekdot*(1.0+self.e*np.cos(tak))/(np.sin(tak)*(1.0-self.e*np.cos(ek)))

        phik = tak + self.omega0

        corr_u = self.cws*np.sin(2.0*phik) + self.cwc*np.cos(2.0*phik)
        corr_r = self.crs*np.sin(2.0*phik) + self.crc*np.cos(2.0*phik)
        corr_i = self.cis*np.sin(2.0*phik) + self.cic*np.cos(2.0*phik)

        uk = phik + corr_u
        rk = A*(1.0 - self.e*np.cos(ek)) + corr_r
        ik = self.i0 + self.i_dot * tk + corr_i
        ukdot = takdot +2.0*(self.cws*np.cos(2.0*uk)-self.cwc*np.sin(2.0*uk))*takdot
        rkdot = A*self.e*np.sin(ek)*n/(1.0-self.e*np.cos(ek)) + 2.0*(self.crs*np.cos(2.0*uk)-self.crc*np.sin(2.0*uk))*takdot
        ikdot = self.i_dot + (self.cis*np.cos(2.0*uk)-self.cic*np.sin(2.0*uk))*2.0*takdot


        xpk = rk*np.cos(uk)
        ypk = rk*np.sin(uk)

        xpkdot = rkdot * np.cos(uk) - ypk * ukdot
        ypkdot = rkdot * np.sin(uk) + xpk * ukdot


        omegak = self.omega_ascension0 + (self.omega_dot_ascension-bOMEGAE84)*tk - bOMEGAE84*self.toe


        omegakdot = (self.omega_dot_ascension - bOMEGAE84)


        xk = xpk*np.cos(omegak) - ypk * np.sin(omegak)*np.cos(ik)
        yk = xpk*np.sin(omegak) + ypk * np.cos(omegak)*np.cos(ik)
        zk = ypk*np.sin(ik)


        xkdot = ( xpkdot-ypk*np.cos(ik)*omegakdot ) * np.cos(omegak) - ( xpk*omegakdot+ypkdot*np.cos(ik)-ypk*np.sin(ik)*ikdot )*np.sin(omegak)
        ykdot = ( xpkdot-ypk*np.cos(ik)*omegakdot ) * np.sin(omegak) + ( xpk*omegakdot+ypkdot*np.cos(ik)-ypk*np.sin(ik)*ikdot )*np.cos(omegak)
        zkdot = ypkdot * np.sin(ik) + ypk * np.cos(ik) * ikdot

        # self.convert_c_language(t_data)

        # print( "BCpos: t, xk, yk, zk:", t_data, xk, yk, zk )
        # print("BCvel: t, Vxk, Vyk, Vzk:", t_data, xkdot, ykdot, zkdot)
        # print("v= ", np.linalg.norm([ xkdot, ykdot, zkdot])*3.6 )
        # print()
        # exit()

        
        return np.array([xk, yk, zk, xkdot , ykdot , zkdot ])

    def get_velocity(self,delta_t=1):
        pos1 = self.get_pos(self.t_data - delta_t)[:3]

        pos2 = self.get_pos(self.t_data)[:3]

        return (pos2-pos1)/(delta_t)


    def convert_c_language(self, t_data):
        print("long double roota                    =  "+ str(self.sqrt_a) + ";")
        print("long double toe                      =  "+ str(self.toe) + ";")
        print("long double m0                       =  "+ str(self.mu0) + ";")
        print("long double e                        =  "+ str(self.e) + ";")
        print("long double delta_n                  =  "+ str(self.delta_n) + ";")
        print("long double smallomega               =  "+ str(self.omega0) + ";")
        print("long double cus                      =  "+ str(self.cws) + ";")
        print("long double cuc                      =  "+ str(self.cwc) + ";")
        print("long double crs                      =  "+ str(self.crs) + ";")
        print("long double crc                      =  "+ str(self.crc) + ";")
        print("long double cis                      =  "+ str(self.cis) + ";")
        print("long double cic                      =  "+ str(self.cic) + ";")
        print("long double idot                     =  "+ str(self.i_dot) + ";")
        print("long double i0                       =  "+ str(self.i0) + ";")
        print("long double bigomega0                =  "+ str(self.omega_ascension0) + ";")
        print("long double earthrate                =  bOMEGAE84;")
        print("long double bigomegadot              =  "+ str(self.omega_dot_ascension) +";")
        print("long double t                        =  "+ str(self.t_data) +";")
        print()


    def velocity_evolution(self):
        t_out = self.toe * 2*60*60
        
        velocities = [ np.linalg.norm(self.get_velocity(t)) for t in range(int(self.toe + 1), int(t_out)) ]

        plt.plot(list(range(self.toe + 1, t_out)), velocities, label=self.name)

        plt.grid()


    def point_satelite_angles(self, origin:np.ndarray, t:int):
        # soh cah toa
        sat_position = self.get_pos(t)[:3]
        orig_sat_vector = origin - sat_position
        r = np.linalg.norm(orig_sat_vector)

        phi = np.arcsin(orig_sat_vector[1]/r)
        theta = np.arctan(orig_sat_vector[2]/orig_sat_vector[0])

        return r, phi, theta


    def show_satelite_epheremide(self):
        print('*** --- ***')
        print("name", self.name)
        print("sqrt_a", self.sqrt_a)
        print("toe", self.toe)
        print("mu0", self.mu0)
        print("e", self.e)
        print("delta_n", self.delta_n)
        print("omega0 ", self.omega0)
        print("cws ", self.cws)
        print("cwc ",self.cwc)
        print("crs", self.crs)
        print('crc', self.crc )
        print("cis", self.cis)
        print("cic", self.cic)
        print("i_dot", self.i_dot)
        print("i0,", self.i0 )
        print("omega_ascension0", self.omega_ascension0   )
        print("omega_dot_ascension", self.omega_dot_ascension )
        print()
    
if __name__ == "__main__":
    pass
    # pos2 = self.get_sat_pose(0.2160*10**6,216005,0.23388256954*10, 0.487698886045*10**-8, 0.260366254952*10**-2, 0.2081*10**4, 0.624172389507*10**-5, 0.157766044140*10**-5, 0.515365465546*10**4, 1, 0.31781250*10**2, -0.409781932831*10**-7, -0.2421438694*10**-7, 0.108880357969*10, -0.808605110220*10**-8, 0.963855990848, 0.468948104999*10**-9)
