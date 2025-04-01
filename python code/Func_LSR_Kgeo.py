# -*- coding: utf-8 -*-

import numpy as np
import math

def LiSparseReciproal_Kgeo(sza,vza,saa,vaa):
    sza= np.radians(sza)
    vza= np.radians(vza)
    n_saa= len(saa)

    saa= np.radians(saa)
    vaa= np.radians(vaa)
    raa= vaa-saa
    b2r= 1.0 
    theta_s= np.arctan(b2r*np.tan(sza))
    theta_v= np.arctan(b2r*np.tan(vza))
    cos_eta= np.cos(theta_s)*np.cos(theta_v)+np.sin(theta_s)*np.sin(theta_v)*np.cos(raa)
    D = ((np.tan(theta_s))**2+(np.tan(theta_v))**2-2.0*np.tan(theta_s)*np.tan(theta_v)*np.cos(raa))**0.5
    h2b= 2.0
    var1= h2b*(D**2+(np.tan(theta_s)*np.tan(theta_v)*np.sin(raa))**2)**0.5
    var2= 1/np.cos(theta_s)+1/np.cos(theta_v)
    cos_t= var1/var2
    for i_t in range(n_saa):
        if cos_t[i_t] > 1.0:
            cos_t[i_t]= 1.0
        if cos_t[i_t] < -1.0:
            cos_t[i_t]= -1.0
    t= np.arccos(cos_t)
    O= (t-np.sin(t)*np.cos(t))*(1/np.cos(theta_s)+1/np.cos(theta_v))
    O = O/math.pi
    sec_theta_s= 1/np.cos(theta_s)
    sec_theta_v= 1/np.cos(theta_v)
    K_geo= O-sec_theta_s-sec_theta_v+0.5*(1+cos_eta)*sec_theta_s*sec_theta_v
    return K_geo
