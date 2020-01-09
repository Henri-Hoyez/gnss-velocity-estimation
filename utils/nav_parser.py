import numpy as np
import pandas as pd
import datetime
from utils.satellite import Satellite
from utils.gps_time import get_second_of_week

GPS_SHIFT = 8
GLONASS_SHIFT = 4

def parse_date(message:str):
    cutted_message = message[4:24]
    year = cutted_message[:4]
    month = cutted_message[5:7]
    day = cutted_message[8:10]
    hour = cutted_message[11:13]
    minute = cutted_message[14:16]
    second = cutted_message[17:19]
    date = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))

    return date


def make_matrix_from_nav_message(message:list):

    message[0] = message[0][25:].split()

    for i in range(1, len(message)):
        message[i] = message[i].replace('\n', '').split()

    return message


def parse_gps_message(message:list):
    # print("GPS :")
    # print(message)
    # Quick and very dirty approach
    # parse the first line.
    gps_name = message[0][:3]
    nav_date = parse_date(message[0])

    nav_data = make_matrix_from_nav_message(message)
    toe = nav_data[3][0]
    t_data = toe                                     # I make t_data = t_toe beacause I don't know t_data
    mu0 = nav_data[1][3]                             
    delta_n = nav_data[1][2]
    e = nav_data[2][1]
    omega0 = nav_data[4][2]
    cws = nav_data[2][2]
    cwc = nav_data[2][0]
    sqrt_a = nav_data[2][3]
    crc = nav_data[4][1]
    crs = nav_data[1][1]
    cic = nav_data[3][1]
    cis = nav_data[3][3]
    omega_ascension0 = nav_data[3][2]
    omega_dot_ascension = nav_data[4][3]
    i0 = nav_data[4][0]
    i_dot = nav_data[5][0]

    clock_drift = nav_data[0][2]


    my_sat = Satellite(float(toe), float(t_data), float(mu0), float(delta_n), float(e), 
    float(omega0), float(cws), float(cwc), float(sqrt_a), float(crc), float(crs), float(cic), 
    float(cis), float(omega_ascension0), float(omega_dot_ascension), float(i0), float(i_dot), float(clock_drift))

    my_sat.set_name(gps_name.replace(" ","0"))
    my_sat.set_ephemeris_date(nav_date)
    my_sat.set_type("GPS")

    return my_sat
    

def parse_glonass_message(message:list):
    pass


def QZSS_message(message:list):
    # print("QZSS: ")
    # print(message)
    # print()
    pass

def parse_nav_file(file:str):
    """Aimed to clear the navigation file
    
    Arguments:
        file {str} -- Path to the rinex file
    Returns:
        cleaned file {str} -- The cleaned rinex file.
    """

    satellites = list()

    rinFile = open(file, "r")

    curr_line = rinFile.readline()

    # Skip the header TODO : make another method ? find the "END OF HEADER" in the list and shrink it ? 
    while curr_line[60:].rstrip() != "END OF HEADER":
        # print(curr_line)
        curr_line = rinFile.readline()

    rinex_input = list( rinFile )

    i = 0                                                  #Set Initial to Zero for the other while part

    while(i < len(rinex_input)):
        if (rinex_input[i][0] == "G"):
            parsed_sat = parse_gps_message(rinex_input[i:i+GPS_SHIFT])
            i += GPS_SHIFT
            satellites.append(parsed_sat)
            continue
        
        if rinex_input[i][0] == "R":
            parse_glonass_message(rinex_input[i:i+GLONASS_SHIFT])
            i += GLONASS_SHIFT
            # satellites.append(parsed_sat)
            continue

        if rinex_input[i][0] == "J":
            QZSS_message(rinex_input[i:i+GPS_SHIFT])
            i += GPS_SHIFT
            # satellites.append(parsed_sat)
            continue
        i += 1
    
    return satellites


if __name__ == "__main__":
    satelites = parse_nav_file('../test_data/autoroute_plus_tunnel.nav')
    print(satelites)
