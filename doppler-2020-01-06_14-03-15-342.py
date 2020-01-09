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



# ftp://cddis.gsfc.nasa.gov/pub/gps/products/wwww/igswwwwd.sp3.Z | template ephemeris

class Doppler:
    # Compute basic parameters at request time
    def __init__(self):
        self.sats = parse_nav_file("data/Very_Bad_Trip/Belgique/autoroute_plus_tunnel.nav")
        print([sat.name for sat in sats])

        self.ri1 = sats[10].get_pos()[:3]
        self.ri2 = sats[9].get_pos()[:3]
        self.ri3 = sats[0].get_pos()[:3]
        
        self.vi1 = sats[10].get_pos()[3:]
        self.vi2 = sats[9].get_pos()[3:]
        self.vi3 = sats[0].get_pos()[3:]
        print(self.vi3)
        exit()
        
    def __get_K_n(self,f_ti,Di,ri,ru,vi):
        c = 299792458
        ai = self.get_line_of_sight(ri,ru)
        return (c*Di)/f_ti - np.inner(vi,ai)
    def visualize(self, ru:np.ndarray, sats:list, t:int):

        indexes = [s.name for s in sats]
        indexes.append('usr')

        positions = []

        x, y, z, vx, vy, vz = s.get_pos(t)
        # vx, vy, vz = s.get_velocity()
        # print("//// V = ", np.linalg.norm([vx, vy, vz])*3.6)
        # exit()
        positions.append([ x, y, z, vx, vy, vz])

        for s in sats:
            positions.append([ru[0], ru[1], ru[2], 0, 0, 0])
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

    def get_line_of_sight(self,ri,ru):
        return (ri-ru)/np.linalg.norm(ri-ru)

    def get_usr_velocity(self):
        #G6 et G23 // 10 & 9
        ru = np.array([4043743.6490  ,  261011.8175 ,  4909156.8423])
        f_ti = 1575.42*10**6
        
        Di1 = 1319.955
        Di2 = -513.404
        Di3 = -2687.413


        k1 = self.__get_K_n(f_ti,Di1,self.ri1,ru,self.vi1)
        k2 = self.__get_K_n(f_ti,Di2,self.ri2,ru,self.vi2)
        k3 = self.__get_K_n(f_ti,Di3,self.ri3,ru,self.vi3)

        a1 = self.get_line_of_sight(self.ri1,ru)
        a2 = self.get_line_of_sight(self.ri2,ru)
        a3 = self.get_line_of_sight(self.ri3,ru)

        x = 0
        y = 0
        z = 0
        #print(x.shape,"||",a1.shape)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.quiver(x,y,z, a1, a2, a3,length=0.1, normalize=True)
        # ax.set_xlim([-1, 1])
        # ax.set_ylim([-1, 1])
        # ax.set_zlim([-1, 1])
        # plt.show()


        K = np.array([k1,k2,k3])
        X = np.array([a1,a2,a3])
        
        v = np.linalg.solve(X,K)

        print(v, "m/s")
        
        # print('FINAL VELOCITY :', np.linalg.norm(v)*3.6, 'km/h')


    def test_velocity(self):
        di = [1319.955, -513.404, -2687.413]

        ru = np.array([4043547.78553915, 254207.686387644, 4909623.02474359])

        sats = [self.sats[i] for i  in [9, 0, 10]]

        
        doppler = Doppler()
        v = doppler.get_usr_velocity(208800, ru, sats,di, [1.57542*10**9]*3  )

        # sats[0].show_trajetcory()

        # doppler.draw_velocity_evolution()
       
       



if __name__ == "__main__":
    sats = parse_nav_file("data/Very_Bad_Trip/Belgique/autoroute_plus_tunnel.nav")

    doppler = Doppler() 

    doppler.test_velocity()

    # doppler = Doppler()
    # doppler.draw_velocity_evolution()
    # doppler.get_usr_velocity()



