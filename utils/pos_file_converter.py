import numpy as np
import pandas as pd
from utils.gps_time import get_second_of_week
from datetime import datetime
import pyproj

from os import _exists

def spheric_to_ECEF(lat, lon, height):
    # return all values in meters
    
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    x, y, z = pyproj.transform(lla, ecef, lon, lat, height, radians=False)

    return x, y, z


def parse_positions(filename:str):

    # print( _exists(filename.split('.')[0]+'.hdf5') )
    # print( filename.split('.')[0]+'.hdf5') 
    # exit()

    if(_exists(filename.split('.')[0]+'.hdf5')):
        print('chargemnet')
        return pd.read_hdf(filename.split('.')[0]+'.hdf5')

    position_file = open(filename)

    lines = list(position_file)
    position_data = list()  # [ gps_time_of_week, x, y, z]

    for line in lines:
        line = line.rstrip()
        splited_line = line.split(" ")
        if splited_line[0] == '%':
            continue
        
        # take the date.
        year, month, day = splited_line[0].split('/')
        hour, minutes, seconds = splited_line[1].split(':')

        gps_sec_of_week = get_second_of_week(datetime(int(year), int(month), int(day), 
                                                      int(hour), int(minutes), int(seconds.split('.')[0]) ) )


        # print(year, month, day, hour, minutes, seconds)




        # print( get_second_of_week(datetime(2019, 11, 26, 10, 00, 0 ) ) )
        # print()

        # exit()
        splited_line = list(dict.fromkeys(splited_line))
        splited_line.remove('')
        
        lat = float(splited_line[2])
        lon = float(splited_line[3])
        height = float(splited_line[4])


        x, y, z = spheric_to_ECEF(lat, lon, height)

        position_data.append([gps_sec_of_week, x, y, z])

        df = pd.DataFrame(position_data, columns=['gps_sec_of_week', 'x', 'y', 'z'])

        # print(df.head())

    df.to_hdf(filename.split('.')[0]+'.hdf5', key='position')

    return df 
        
           


if __name__ == "__main__":
    df = parse_positions("test_data/autoroute_plus_tunnel.pos")