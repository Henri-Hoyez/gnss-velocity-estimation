import numpy as np


class Doppler:

    def __init__():
        # Compute basic parameters at request time
        toe = None
        t_data = None
        mu0 = None
        delta_n = None
        e = None
        omega0 = None
        Cws = None
        Cwc = None
        a= None
        Crc = None
        Crs = None
        Cic = None
        Cis = None
        omega_ascension0 = None
        omega_dot_ascension = None
        i0 = None
        i_dot = None


        GM = 3.986004418*np.power(10,14)
        omega_e = 7.2921151467*np.power(10,-5)

        t = t_data-toe
        mu  = mu0 + (np.sqrt(GM/np.power(a,3)) + delta_n) * t
        # E = mu + e * np.sin(E) developpement limité
        v = np.arctan((np.sqrt(1-e**2)*np.sin(E))/(np.cos(E)-e))

        # Correct for orbital perturbations
        omega = omega0 + Cwc * np.cos(2*(omega0+v)) + Cws * np.sin(2*(omega0+v))
        r = a*(1-e*np.cos(E))+ Crc * np.cos(2*(omega0+v)) + Crs*np.sin(2*(omega0+v))
        i = i0 + i_dot * t + Cic * np.cos(2*(omega0+v)) + Cis*np.sin(2*(omega0+v))
        
        # Compute  the  right  ascension
        omega_ascension = omega_ascension0 + (omega_dot_ascension-omega_e)*toe

        # Convert satellite position from orbital frame to ECEF frame
        r_vect = np.array([r*np.cos(v),r*np.sin(v),0])

        # Rotation matrix
        R = [[np.cos(omega_ascension)*np.cos(omega)-np.sin(omega_ascension)*np.sin(omega)*np.cos(i), -np.cos(omega_ascension)*np.sin(omega)-np.sin(omega_ascension)*np.cos(omega)*np.cos(i), np.sin(omega_ascension)*np.sin(i)]
        ,[np.sin(omega_ascension)*np.cos(omega)+np.cos(omega_ascension)*np.sin(omega)*np.cos(i), -np.sin(omega_ascension)*np.sin(omega)+np.cos(omega_ascension)*np.cos(omega)*np.cos(i), -np.cos(omega_ascension)*np.sin(i)]
        ,[np.sin(omega)*np.sin(i), np.cos(omega)*np.sin(i), np.cos(i)]]





