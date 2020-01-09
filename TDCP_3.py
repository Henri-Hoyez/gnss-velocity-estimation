import numpy as np
from utils.nav_parser import parse_nav_file

# ftp://cddis.gsfc.nasa.gov/pub/gps/products/wwww/igswwwwd.sp3.Z | template ephemeris

class TDCP:
    # Compute basic parameters at request time
    def __init__(self):
        sats = parse_nav_file("data/Very_Bad_Trip/Belgique/autoroute_plus_tunnel.nav")
        print([sat.name for sat in sats])
        print("blabla",-sats[13].get_pos()[:3]/1000)
        self.x_sk = np.array([[sats[10].get_pos()[:3],sats[9].get_pos()[:3],-sats[0].get_pos()[:3],sats[13].get_pos()[:3]],
        [sats[10].get_pos(sats[10].toe+1)[:3],sats[9].get_pos(sats[9].toe+1)[:3],sats[0].get_pos(sats[0].toe+1)[:3],sats[13].get_pos(sats[13].toe+1)[:3]]])
        self.x_rk = np.array([[4043743.6490 ,   261011.8175  , 4909156.8423],
        [4043730.5731  ,  261041.2499  , 4909152.9334]])

        self.phase1 = np.array([106543859.919,101729930.484,103940539.651,110045178.0752])
        self.phase2 = np.array([106542319.194,101730136.853,103942948.132,110047781.4172])
        print("groundtruth", np.linalg.norm(self.x_rk[1]-self.x_rk[0])*3.6)
        
        
    def get_phase_difference(self,phase1,phase2):
        return phase2 - phase1
        
    
    def __get_K_n(self,wavelength,phase1,phase2,delta_S,delta_G,delta_t):
        print("**** deltaS : ",delta_S ," deltaG : ", delta_G ," phase diff : ", wavelength * self.get_phase_difference(phase1,phase2)," result : ",delta_S - delta_G - wavelength * self.get_phase_difference(phase1,phase2))
        return (delta_S - delta_G - wavelength * self.get_phase_difference(phase1,phase2)) / (delta_t)

    def get_line_of_sight(self,ri,ru):
        return (ri-ru)/np.linalg.norm(ri-ru)


    def get_usr_velocity(self):
        #G06 G23 G03 G19// 10 & 9 & 0 & 19
        
        phase1 = np.array([106543859.919,101729930.484,103940539.651,110045178.0752])
        phase2 = np.array([106542319.194,101730136.853,103942948.132,110047781.4172])

        c = 299792458
        f = 1575.42*10**6
        wavelength = c/f
        y = wavelength * self.get_phase_difference(phase1,phase2)

        a1 = np.append(self.get_line_of_sight(self.x_sk[1,0],self.x_rk[1]),1)
        a2 = np.append(self.get_line_of_sight(self.x_sk[1,1],self.x_rk[1]),1)
        a3 = np.append(self.get_line_of_sight(self.x_sk[1,2],self.x_rk[1]),1)
        a4 = np.append(self.get_line_of_sight(self.x_sk[1,3],self.x_rk[1]),1)

        X = np.array([a1,a2,a3,a4])
        W = np.diag([10+150**((-49.000)/10),10+150**((-50.000)/10),10+150**((-52.000)/10),10+150**((-45.000)/10)])
        W = np.linalg.inv(W)
        result = np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(X),W),X)),np.transpose(X)),W),y)
        print(W)
        print(np.linalg.norm(result[:3]))
        
        #v = np.linalg.lstsq(np.transpose(X),K,rcond=-1)
        ## print(np.linalg.inv(X))
        ## l = np.pi/6
        ## phi = np.pi/3
        ## to_ENU = np.array([[-np.sin(l),np.cos(l),0]
        ## ,[-np.cos(l)*np.sin(phi),-np.sin(l)*np.sin(phi),np.cos(phi)]
        ## ,[np.cos(l)*np.cos(phi),np.sin(l)*np.cos(phi),np.sin(phi)]])
        ## v = np.dot(to_ENU,v)
        #print("////",np.linalg.norm(v[0])*3.6)
#
        #return v


    def TDCP_matlab(self):
        ru_t = self.x_rk
        ri_jt = self.x_sk
        dphi = self.phase2-self.phase1
        SVcom = np.array([1023*10**6,1023*10**6,1023*10**6,1023*10**6])
        dt_data = 1
        el = np.array([np.pi/2,np.pi/4,np.pi/3,np.pi/6])
        
        num1 = np.array(ru_t[0]-ri_jt[0])
        den1 = np.diag(np.sqrt(np.dot(num1,np.transpose(num1))))
        e_jt1 = -np.array([np.divide(num1[:,0],den1),np.divide(num1[:,1],den1),np.divide(num1[:,2],den1)])
        print("aaaa",num1)
        print("bbbb",den1)
        print("cccc",e_jt1)

        num2 = ru_t[1]-ri_jt[1]
        den2 = np.diag(np.sqrt(np.dot(num2,np.transpose(num2))))
        e_jt2 = -np.array([np.divide(num2[:,0],den2),np.divide(num2[:,1],den2),np.divide(num2[:,2],den2)])

    
        SVDoppler = np.diag(e_jt2.T.dot(ri_jt[1].T)) - np.diag(e_jt1.T.dot(ri_jt[0].T))
        Dgeometry = np.dot(np.transpose(e_jt2),ru_t[0])-np.dot(np.transpose(e_jt1),ru_t[0])

        d_phi_adj = (3*10**8)/(1575.42*10**6)*dphi/100#-SVDoppler+Dgeometry
        print("******")
        print(dphi)
        print(SVDoppler)
        print(Dgeometry)
        print("******")
        dstd_el = np.sqrt(2)*np.divide(1,np.sin(el))
        W = np.diag(np.divide(1,dstd_el**2))
        H = np.vstack((e_jt2,np.ones((SVcom.shape[0],1)).T)).T
        # Db_dcdtu = np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(H),W),H)),np.transpose(H)),W),d_phi_adj)
        
        print(np.linalg.inv(H.T.dot(W).dot(H)).dot(H.T).dot(W))
        Db_dcdtu = np.linalg.inv(H.T.dot(W).dot(H)).dot(H.T).dot(W).dot(d_phi_adj)
        vu = (np.divide(Db_dcdtu[:3],dt_data),Db_dcdtu[3])
        print(vu)



if __name__ == "__main__":
    sats = parse_nav_file("data/Very_Bad_Trip/Belgique/autoroute_plus_tunnel.nav")
    # print(sats[1].name)
    # print(sats[1].get_position()/1000)
    # print(np.linalg.norm(sats[1].get_velocity())*3.6)
    # velocity = Doppler(208784,1574762406,1.28612756881,0.489020369661*10**-8,0.260362529662*10**-2,0.7931942133,0.622682273388*10**-5,0.1460313797*10**-5,0.515365547943*10**4,0.25971875*10**3,0.2853125*10**2,-0.204890966415*10**-7,0.577419996262*10**-7,0.108886242411*10,-0.815426822943*10**-8,0.963852456438,0.431089385164*10**-9)
    tdcp = TDCP()
    # tdcp.get_usr_velocity()
    tdcp.TDCP_matlab()
