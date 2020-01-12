import numpy as np
from utils.nav_parser import parse_nav_file
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

from utils.pos_file_converter import parse_positions
from utils.parserRinex import obsToDataframeFinal
from utils.satelite_manager import get_satelites
from sklearn.metrics import pairwise_distances
import seaborn as sns

from tqdm import tqdm



# ftp://cddis.gsfc.nasa.gov/pub/gps/products/wwww/igswwwwd.sp3.Z | template ephemeris

class Doppler:
    """ 
    This Class is aimed to compute the receiver velocity from the Doppler effect 
    """
    def __init__(self):
        # set the good t_data ( 26 nov à 10:00:00 )
        self.sats = parse_nav_file("test_data/autoroute_plus_tunnel.nav")


    def get_k(self, D_i, f_ti, v_i, a_i):
        c = 3.00*10**8
        res =  (c * (D_i / f_ti)) + np.vdot(v_i, a_i)% 2 * np.pi
        # print(  np.linalg.norm(v_i)*3.6 )
        return res


    def visualize(self, ru:np.ndarray, sats:list, t:int):
        """This function aimed to get a visualisation from the satelite data
        
        Arguments:
            ru {np.ndarray} -- The user position in the ECEF way
            sats {list} -- Satellites data which is used to compute satellites positions
            t {int} -- Time of computation, in gps time of week
        """

        indexes = [s.name for s in sats]
        indexes.append('usr')

        positions = []

        for s in sats:
            x, y, z, vx, vy, vz = s.get_pos(t)
            # vx, vy, vz = s.get_velocity()
            # print("//// V = ", np.linalg.norm([vx, vy, vz])*3.6)
            # exit()
            positions.append([ x, y, z, vx, vy, vz])

        positions.append([ru[0], ru[1], ru[2], 0, 0, 0])
        positions.append([0, 0, 0, 0, 0, 0])
        indexes.append("Orig")
        df = pd.DataFrame(np.array(positions), columns=['x','y','z', 'vx', 'vy', 'vz'], index=indexes)


        fig = plt.figure()
        ax = Axes3D(fig)

        for i in range( df.shape[0] ): #plot each point + it's index as text above
            x,y,z = df['x'][i], df['y'][i] , df['z'][i]
            vx, vy, vz = df['vx'][i]*3.6, df['vy'][i]*3.6 , df['vz'][i]*3.6

            ax.scatter(x, y, z, color='b') 
            ax.text( x, y, z,  '%s' % (str(indexes[i])), size=20, zorder=1,  
            color='k') 

            ax.quiver(x, y, z, vx, vy, vz  )

        plt.show()

    def get_available_satelites(self, t:int, sats:list, sats_obs:pd.DataFrame):
        """This function aimed to filter and return only the available satellites
        
        Arguments:
            t {int} -- time ine GPS week
            sats {list} -- Satelites data
            sats_obs {pd.DataFrame} -- The parsed RINEX data observation
        
        Returns:
            [type] -- a list of available satelites
        """
        obs_names = list(sats_obs.loc[2081, t].index)
        available_sats = list()


        for s in sats:

            if(t <= s.toe * 2*60*60 and t >= s.toe and s.name in obs_names):
                available_sats.append(s)


        return available_sats


    def get_usr_velocity(self,t:int ,ru:np.ndarray, sats:list, sats_obs:pd.DataFrame, f_ti:list ):
        """Compute the receiver velocity
        
        Arguments:
            t {int} -- Time in GPS time of week
            ru {np.ndarray} -- receiver position in ECEF
            sats {list} -- list of available satelites, only the first three sats are taken
            sats_obs {pd.DataFrame} -- parsed observation file
            f_ti {list} -- The emited frequencie of the satellite
        
        Returns:
            float -- velocity in meter peer second
        """

        # sats[0].show_trajetcory()

        obs = sats_obs.loc[2081, t]

        D_i = [obs.loc[s.name].doppler for s in sats[:3]]

        a1 = (sats[0].get_pos(t)[:3] - ru)/(np.linalg.norm(sats[0].get_pos(t)[:3] - ru))
        a2 = (sats[1].get_pos(t)[:3] - ru)/(np.linalg.norm(sats[1].get_pos(t)[:3] - ru))
        a3 = (sats[2].get_pos(t)[:3] - ru)/(np.linalg.norm(sats[2].get_pos(t)[:3] - ru))

        k1 = self.get_k(D_i[0], f_ti[0], sats[0].get_pos(t)[3:], a1)
        k2 = self.get_k(D_i[1], f_ti[1], sats[1].get_pos(t)[3:], a2)
        k3 = self.get_k(D_i[2], f_ti[2], sats[2].get_pos(t)[3:], a3)

        K = np.array([k1, k2, k3])

        M = np.array([a1, a2, a3])

        v = np.dot(np.linalg.inv(M), K)

        return np.linalg.norm(v)
    

    

    def draw_velocity_evolution(self, position_file_location:str, nav_file_location:str, obs_file_location:str ):
        """Compute, for all time in the position file location, the velocity with the doppler method
        
        Arguments:
            position_file_location {str} -- The position file location
            nav_file_location {str} -- The nav file location
            obs_file_location {str} -- The obs file location
        """
        
        sats = parse_nav_file(nav_file_location)

        # user_positions = parse_positions(position_file_location).set_index(["gps_sec_of_week"])

        user_positions = pd.read_hdf(position_file_location + '.hdf5', 'data')
        user_positions = user_positions.reset_index().drop_duplicates(subset ="gps_sec_of_week", inplace = False).set_index(["gps_sec_of_week"])

        # user_positions.to_hdf( position_file_location+'.hdf5', key='data')

        # sats_observation = obsToDataframeFinal(obs_file_location)

        # sats_observation.to_hdf(obs_file_location + '.hdf5', key='data')

        sats_observation = pd.read_hdf(obs_file_location + '.hdf5', 'data')

        time = list(sats_observation.index.get_level_values(1)[107000:107000+ 10000])

        vs = list()
        true_velocity = list()
        time_vel= list()

        for t in tqdm(time) : 
            tmp_sats = self.get_available_satelites(t, sats, sats_observation)
            
            try:
                ru = user_positions.loc[int(t)]
            except Exception :
                continue
            

            if len(tmp_sats) < 3:
                continue
            
            if len(tmp_sats) > 3: 
                tmp_sats = self.best_satelites(t, tmp_sats, ru)
            else:
                continue

           
            vs.append(self.get_usr_velocity( t, ru, tmp_sats, sats_observation, [1.57542*10**9] * 3))


            # try:
            #     true_velocity.append(np.linalg.norm([ru -user_positions.loc[t-1]]) * 3.6 )
            # except Exception:
            #     true_velocity.append(true_velocity[-1])
            

            time_vel.append(t)

        plt.plot(time_vel,vs, 'o', label="doppler")
        # plt.plot(time_vel, true_velocity, 'o', label="derivative")
        plt.legend()
        plt.grid()
        plt.show()



        

    def best_satelites(self, t:int, sats:list, user_position:np.ndarray):

        spheric_coordonate = np.array([s.point_satelite_angles(user_position, t) for s in sats])

        names = [s.name for s in sats]

        satelite_dataframe = pd.DataFrame(spheric_coordonate, columns=['r', 'phi', 'theta'], index=names).sort_values(['theta', 'phi'])

        best_names = list(satelite_dataframe.index)

        best_names = [best_names[0], 
                     best_names[int(len(best_names)/2 -1)], 
                     best_names[-1]]

        best_satelites = list(filter(lambda s: s.name in best_names, sats))
 
        return best_satelites
        
        

    def test_velocity(self):
        
        # di = [-2408.521, -2594.214, -2687.413]

        ru = np.array([4043547.78553915, 254207.686387644, 4909623.02474359])

        sats_obs = pd.read_hdf("test_data/satelites_observation.hdf5", 'data')

        sats_obs = sats_obs.reset_index()

        sats_obs.second_of_week = sats_obs.second_of_week - 26

        sats_obs =sats_obs.set_index(["n_week", "second_of_week", "name"])

        sats = self.get_available_satelites(208800, self.sats, sats_obs)

        print([s.name for s in sats[:3]])



        v = self.get_usr_velocity(208800, ru, sats, sats_obs, [1.57542*10**9] * 3)

        print('v1 :', v)

        sats = self.best_satelites(208800, sats, ru)

        # self.visualize(ru, sats, 208800)

        
        v = self.get_usr_velocity(208800, ru, sats , sats_obs, [1.57542*10**9] * 3)

        print('v2 :',v)

        # self.draw_velocity_evolution('test_data/autoroute_plus_tunnel.pos', 
        # 'test_data/autoroute_plus_tunnel.nav', 
        # 'test_data/autoroute_plus_tunnel.obs' )

        # print(v)



       
       



if __name__ == "__main__":
    sats = parse_nav_file("test_data/autoroute_plus_tunnel.nav")

    doppler = Doppler() 

    doppler.test_velocity()

    # doppler = Doppler()
    # doppler.draw_velocity_evolution()
    # doppler.get_usr_velocity()







