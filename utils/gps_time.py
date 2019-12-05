import numpy as np
import matplotlib.pyplot as plt
from datetime import date, timedelta, time, datetime


def get_gps_week(time_of_week:date=None):
    if(time_of_week is None):
        time_of_week = date.today()

    ref = date(1980, 1, 6)

    return int((time_of_week - ref).days/7)
    



def get_second_of_week(time_of_week:date=None):
    if(time_of_week is None):
        time_of_week = datetime.now() 

    gps_day_second = (time_of_week.isoweekday()+1) % 7 * 24 * 60**2  # Get second of all day before

    gps_day_second += time_of_week.hour * 60**2 + time_of_week.minute * 60 + time_of_week.second

    return gps_day_second


    



if __name__ == "__main__":
    
    print(get_gps_week())
    print(get_second_of_week())



    
    










