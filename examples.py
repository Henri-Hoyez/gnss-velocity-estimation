import pandas as pd
import numpy as np
from doppler import Doppler
from utils.nav_parser import parse_nav_file
from utils.pos_file_converter import parse_positions
from utils.parserRinex import obsToDataframeFinal


"""
This file aim to give somme examples on how to use our scripts.
"""

def doppler_import_data(obs_file_location, nav_file_location, position_file_location, load_hdf5=False):
    sats = parse_nav_file(nav_file_location)

    if load_hdf5:
        user_positions = pd.read_hdf(position_file_location + '.hdf5', 'data')
        sats_observation = pd.read_hdf(obs_file_location + '.hdf5', 'data')
    else:
        user_positions = parse_positions(position_file_location).set_index(["gps_sec_of_week"])
        user_positions = user_positions.reset_index().drop_duplicates(subset ="gps_sec_of_week", inplace = False).set_index(["gps_sec_of_week"])
        sats_observation = obsToDataframeFinal(obs_file_location)   

        user_positions.to_hdf( position_file_location+'.hdf5', key='data')
        sats_observation.to_hdf(obs_file_location + '.hdf5', key='data')
        

    sats_observation = sats_observation.reset_index().rename(columns={"second_of_week":'gps_sec_of_week'}).set_index(['gps_sec_of_week']) 
    user_positions = user_positions.merge(sats_observation, left_index=True, right_index=True).reset_index().set_index(['n_week','gps_sec_of_week', 'name']) 
    
    return sats, user_positions


def doppler_velocity_parameters_wrapper(time, gps_week, obs_file_location, 
                                        nav_file_location, position_file_location, load_hdf5=False):

    sats, rinex_dataframe = doppler_import_data(obs_file_location, nav_file_location, position_file_location, load_hdf5=True)
    
    # ru is the user position at a given time
    ru = [rinex_dataframe.loc[gps_week, time].iloc[0].x, 
          rinex_dataframe.loc[gps_week, time].iloc[0].y, 
          rinex_dataframe.loc[gps_week, time].iloc[0].z ]

    f_ti = [1.57542*10**9]*3

    return sats, rinex_dataframe, ru, f_ti




def doppler_velocity_V1():
    # import the data
    time = 208800
    gps_week = 2081

    sats, rinex_dataframe, ru, f_ti = doppler_velocity_parameters_wrapper(time, gps_week, 
                                                                          'test_data/autoroute_plus_tunnel.obs',
                                                                          'test_data/autoroute_plus_tunnel.nav',
                                                                          'test_data/autoroute_plus_tunnel.pos',)
    
    dp = Doppler()

    available_sats = dp.get_available_satelites(time, sats, rinex_dataframe)
    best_satelites = dp.best_satelites(time, available_sats, ru)

    velocity = dp.get_usr_velocity(time, ru, best_satelites, rinex_dataframe, f_ti)

    print("** doppler_velocity_V1 result: ", velocity, 'km/h')



def doppler_velocity_V2():
    time = 208800
    gps_week = 2081

    sats, rinex_dataframe, ru, _ = doppler_velocity_parameters_wrapper(time, gps_week, 
                                                                       'test_data/autoroute_plus_tunnel.obs',
                                                                        'test_data/autoroute_plus_tunnel.nav',
                                                                        'test_data/autoroute_plus_tunnel.pos',)

    dp= Doppler()

    available_sats = dp.get_available_satelites(time, sats, rinex_dataframe)
    best_satelites = dp.best_satelites(time, available_sats, ru)

    velocity_vector = dp.speed_for_the_win(time, ru, best_satelites, rinex_dataframe)[:3]

    velocity = np.linalg.norm(velocity_vector)

    print("** doppler_velocity_V2 result: ", velocity, 'km/h')
    


def main():
    doppler_velocity_V1()

    doppler_velocity_V2()



if __name__ == "__main__":
    main()