import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from utils.nav_parser import parse_nav_file
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

from utils.pos_file_converter import parse_positions
from utils.parserRinex import obsToDataframeFinal
from utils.satelite_manager import get_satelites
from sklearn.metrics import pairwise_distances
import seaborn as sns



# ftp://cddis.gsfc.nasa.gov/pub/gps/products/wwww/igswwwwd.sp3.Z | template ephemeris

class Doppler:
    # Compute basic parameters at request time
    def __init__(self):
        self.sats = parse_nav_file("data/Very_Bad_Trip/Belgique/autoroute_plus_tunnel.nav")
        print([sat.name for sat in sats])

        self.ri1 = sats[10].get_pos()[:3]
        self.ri2 = sats[9].get_pos()[:3]
        self.ri3 = sats[0].get_pos()[:3]

    def get_k(self, D_i, f_ti, v_i, a_i):
        c = 3.00*10**8
        res =  (c * (D_i / f_ti)) + np.vdot(v_i, a_i)
        # print(  np.linalg.norm(v_i)*3.6 )
        return res


    def visualize(self, ru:np.ndarray, sats:list, t:int):

        indexes = [s.name for s in sats]
        indexes.append('usr')

        positions = []
        for s in sats:
            x, y, z, vx, vy, vz = s.get_pos(t)
            positions.append([ x, y, z, vx, vy, vz])
            positions.append([ru[0], ru[1], ru[2], 0, 0, 0])

        positions.append([ru[0], ru[1], ru[2], 0, 0, 0])
        positions.append([0, 0, 0, 0, 0, 0])
        indexes.append("Orig")
        df = pd.DataFrame(np.array(positions), columns=['x','y','z', 'vx', 'vy', 'vz'], index=indexes)


        fig = plt.figure()
        ax = Axes3D(fig)

        for i in range( df.shape[0] ): #plot each point + it's index as text above
        # print(df)
        # print(np.linalg.norm(sats[0].get_velocity(t)))
            x,y,z = df['x'][i], df['y'][i] , df['z'][i]
            vx, vy, vz = df['vx'][i]*3.6, df['vy'][i]*3.6 , df['vz'][i]*3.6

            ax.scatter(x, y, z, color='b') 
            ax.text( x, y, z,  '%s' % (str(indexes[i])), size=20, zorder=1,  
            color='k') 

            ax.quiver(x, y, z, vx, vy, vz  )

        plt.show()

    def get_available_satelites(self, t:int, sats:list, sats_obs:pd.DataFrame):
        obs_names = list(sats_obs.loc[2081, t].index)
        available_sats = list()

    def get_usr_velocity(self,toe,ru,sats,Di,f_ti):
        #G6 et G23 // 10 & 9
        ru = np.array([4043743.6490  ,  261011.8175 ,  4909156.8423])
        f_ti = 1575.42*10**6

        for s in sats:

            if(t <= s.toe * 2*60*60 and t >= s.toe and s.name in obs_names):
                available_sats.append(s)

        return available_sats


    def get_usr_velocity(self,t:int ,ru:np.ndarray, sats:list, sats_obs:pd.DataFrame, f_ti:list ):

        # sats[0].show_trajetcory()

        obs = sats_obs.loc[2081, t]

        D_i = [obs.loc[s.name].doppler for s in sats]

        a1 = (sats[0].get_pos(t) [:3]- ru) / (np.linalg.norm(sats[0].get_pos(t)[:3]- ru))
        a2 = (sats[1].get_pos(t) [:3]- ru) / (np.linalg.norm(sats[1].get_pos(t)[:3]- ru))
        a3 = (sats[2].get_pos(t) [:3]- ru) / (np.linalg.norm(sats[2].get_pos(t)[:3]- ru))

        k1 = self.get_k(D_i[0], f_ti[0], sats[0].get_pos(t)[3:], a1)
        k2 = self.get_k(D_i[1], f_ti[1], sats[1].get_pos(t)[3:], a2)
        k3 = self.get_k(D_i[2], f_ti[2], sats[2].get_pos(t)[3:], a3)

        K = np.array([k1, k2, k3])

        M = np.array([a1, a2, a3])

        v = np.dot(np.linalg.inv(M), K)

        return np.linalg.norm(v)
    

    

    def draw_velocity_evolution(self, position_file_location:str, nav_file_location:str, obs_file_location:str ):
        sats = parse_nav_file(nav_file_location)

        # user_positions = parse_positions(position_file_location).set_index(["gps_sec_of_week"])

        user_positions = pd.read_hdf(position_file_location + '.hdf5', 'data')
        user_positions = user_positions.reset_index().drop_duplicates(subset ="gps_sec_of_week", inplace = False).set_index(["gps_sec_of_week"])

        # user_positions.to_hdf( position_file_location+'.hdf5', key='data')

        # sats_observation = obsToDataframeFinal(obs_file_location)

        # sats_observation.to_hdf(obs_file_location + '.hdf5', key='data')

        sats_observation = pd.read_hdf(obs_file_location + '.hdf5', 'data')

        time = list(sats_observation.index.get_level_values(1))

        vs = list()
        true_velocity = list()
        time_vel= list()

        for t in time : 
            tmp_sats = self.get_available_satelites(t, sats, sats_observation)
            if len(tmp_sats) < 3:
                continue
            
            try:
                ru = user_positions.loc[int(t)]
            except Exception :
                continue

            vs.append(self.get_usr_velocity( t, ru, tmp_sats, sats_observation ,[1.57542*10**9] * 3))
            try:
                true_velocity.append(np.linalg.norm([user_positions.loc[t]-user_positions.loc[t-1]])*3.6 )
            except Exception:
                true_velocity.append(true_velocity[-1])

            time_vel.append(t)

        plt.plot(time_vel,vs, label="doppler")
        plt.plot(time_vel, true_velocity, label="derivative")
        plt.legend()
        plt.grid()
        plt.show()



        

    def best_satelites(self, t:int, sats:list, user_position:np.ndarray):

        spheric_coordonate = np.array([s.point_satelite_angles(user_position, t) for s in sats])

        names = [s.name for s in sats]

        satelite_dataframe = pd.DataFrame(spheric_coordonate, columns=['r', 'phi', 'theta'], index=names).sort_values(['phi', 'theta'])

        

        print(satelite_dataframe)

        best_satelites = [sats[0], sats[ int(len(sats)/2) ], sats[-1]]

        
 
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

        self.visualize(ru, sats, 208800)

        
        v = self.get_usr_velocity(208800, ru, sats , sats_obs, [1.57542*10**9] * 3)

        print('v2 :',v)
        # self.draw_velocity_evolution('test_data/autoroute_plus_tunnel.pos', 
        # 'test_data/autoroute_plus_tunnel.nav', 
        # 'test_data/autoroute_plus_tunnel.obs' )

        # print(v)



       



if __name__ == "__main__":
    sats = parse_nav_file("data/Very_Bad_Trip/Belgique/autoroute_plus_tunnel.nav")

    doppler = Doppler() 

    doppler.test_velocity()

    # doppler = Doppler()
    # doppler.draw_velocity_evolution()
    # doppler.get_usr_velocity()



