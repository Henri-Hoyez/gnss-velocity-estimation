import numpy as np
from utils.nav_parser import parse_nav_file
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from utils.pos_file_converter import parse_positions
from utils.parserRinex import obsToDataframeFinal
from utils.satelite_manager import get_satelites


import pandas as pd

# ftp://cddis.gsfc.nasa.gov/pub/gps/products/wwww/igswwwwd.sp3.Z | template ephemeris

class Doppler:
    # Compute basic parameters at request time
    def __init__(self):
        self.sats = parse_nav_file("test_data/autoroute_plus_tunnel.nav")

        # set the good t_data ( 26 nov à 10:00:00 )

        sats[9].t_data = 208800
        sats[0].t_data = 208800
        sats[10].t_data = 208800

        print(sats[10].name,   ' - ', sats[9].name,   ' - ',sats[0].name )
        print(sats[10].t_data, ' - ', sats[9].t_data, ' - ',sats[0].t_data )

        # stats_indx = [9, 0, 10]

        # for s in stats_indx:
        #     sats[s].show_info()


        self.ri1 = sats[10].get_position()
        self.ri2 = sats[9].get_position()
        self.ri3 = sats[0].get_position()
        #print(self.ri2)

        self.vi1 = sats[10].get_velocity()
        self.vi2 = sats[9].get_velocity()
        self.vi3 = sats[0].get_velocity()


    def test_angles(self):
        
        ri12 = np.arccos(np.dot(self.ri1, self.ri2) / (np.linalg.norm(self.ri1) * np.linalg.norm(self.ri2)) )
        ri13 = np.arccos(np.dot(self.ri1, self.ri3) / (np.linalg.norm(self.ri1) * np.linalg.norm(self.ri3)) )
        ri23 = np.arccos(np.dot(self.ri2, self.ri3) / (np.linalg.norm(self.ri2) * np.linalg.norm(self.ri3)) )  

        # print(ri12, ri13, ri23)


    def corelation_coeficiant(self, v1:np.ndarray, V2:np.ndarray):
        pass

    def get_k(self, D_i, f_ti, v_i, a_i):
        c = 3.00*10**8
        res =  c * D_i / f_ti + np.vdot(v_i, a_i)
        # print( c * D_i / f_ti )
        return res


    def azimuth_to_ECEF(self, vector:np.ndarray):
        # vector  => [ x_sat , y_sat , z_sat]
        # user_position => [ x_usr , y_usr , z_usr]

        r = np.linalg.norm(vector)

        theta = np.arctan( (vector[0]) / ( vector[2] ) )
        phi = np.arctan( vector[1] / np.sqrt( vector[0]**2 + vector[2]**2 ) )

        x, y, z = r*np.cos(phi) * np.cos(theta), r*np.sin(theta), r*np.cos(phi) * np.sin(theta) 

        return np.array([x, y, z])

    def project_on_ai(self, vi, ai):
        return np.dot(vi, ai)/np.linalg.norm(ai)

    def ecef_to_enu(self, origin, v):
        E = None
        N = None
        U = None





    def get_usr_velocity(self,t ,ru:np.ndarray, sats:list, D_i:list, f_ti:list ):

        a1 = (self.ri1 - ru) / (np.linalg.norm(self.ri1 - ru))
        a2 = (self.ri2 - ru) / (np.linalg.norm(self.ri2 - ru))
        a3 = (self.ri3 - ru) / (np.linalg.norm(self.ri3 - ru))

        k1 = self.get_k(D_i[0], f_ti[0], self.vi1, a1)
        k2 = self.get_k(D_i[1], f_ti[1], self.vi2, a2)
        k3 = self.get_k(D_i[2], f_ti[2], self.vi3, a3)

        K = np.array([k1, k2, k3])

        M = np.array([a1, a2, a3])

        v = np.dot(np.linalg.inv(M), K)

        vai = np.array([
            self.project_on_ai(self.vi1, a1),
            self.project_on_ai(self.vi2, a2),
            self.project_on_ai(self.vi3, a3)
        ])





        print('Before: ', v)
        print('after : ', vai)
        print('projection vai:', self.azimuth_to_ECEF(vai))
        # print('projection v  :', self.azimuth_to_ECEF(v))
        
        print('norme vai:', np.linalg.norm(self.azimuth_to_ECEF(vai)))
        print('norme v  :', np.linalg.norm(self.azimuth_to_ECEF(v)))
        
    
        return np.linalg.norm(self.azimuth_to_ECEF(vai))


    



    def draw_velocity_evolution(self):
        #G6 et G23 // G10 & G09

        # position_history = parse_positions('test_data/autoroute_plus_tunnel.pos')

        # position_history = position_history.set_index('gps_sec_of_week')

        position_history = pd.read_hdf('test_data/autoroute_plus_tunnel.hdf5').set_index('gps_sec_of_week')


        position_history =  position_history.reset_index().drop_duplicates(subset='gps_sec_of_week', keep='first').set_index('gps_sec_of_week')

        # rinex_dataframe = obsToDataframeFinal('test_data/autoroute_plus_tunnel.obs')

        # rinex_dataframe.to_hdf("test_data/rinex_data.hdf5", key='rinex')

        rinex_dataframe = pd.read_hdf('test_data/rinex_data.hdf5')


        # print(position_history)
 

        v_data = list()

        for time in list(position_history.index):
            try:
                reciver_informations = rinex_dataframe.loc[(2081, time)]
            except:
                print(time)
            sats = get_satelites(self.sats, time, reciver_informations)

            f_ti = [1.57542*10**9]*3 
            data = position_history.loc[time]

            # print(reciver_informations)
            ru = [np.array([data['x'], data['y'], data['z']])][0]
            d_i = [ reciver_informations.loc[s.name]['doppler'] for s in sats ]  
            
            if (len(d_i) < 3):
                continue

            v = self.get_usr_velocity(time, ru, sats, d_i, f_ti)
            v_data.append(np.linalg.norm(v)*3.6)

        plt.plot(list(position_history.index)[:len(v_data)], v_data)
        plt.grid()
        plt.show()
        
        
        
        # print('FINAL VELOCITY :', np.linalg.norm(v)*3.6, 'km/h')


    def test_velocity(self):
        di = [1319.955, -513.404, -2687.413]

        ru = np.array([4043547.78553915, 254207.686387644, 4909623.02474359])

        sats = [self.sats[i] for i  in [9, 0, 10]]
        
        doppler = Doppler()
        v = doppler.get_usr_velocity(208800, ru, sats,di, [1.57542*10**9]*3  )
       
       



if __name__ == "__main__":
    sats = parse_nav_file("test_data/autoroute_plus_tunnel.nav")

    doppler = Doppler() 

    doppler.test_velocity()

    # doppler = Doppler()
    # doppler.draw_velocity_evolution()
    # doppler.get_usr_velocity()







