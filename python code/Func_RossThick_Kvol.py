# -*- coding: utf-8 -*-

import numpy as np
import math

def RossThick_Kvol(sza,vza,saa,vaa):
    sza= np.radians(sza)
    vza= np.radians(vza)

    saa= np.radians(saa)
    vaa= np.radians(vaa)
    raa= vaa-saa
    cos_kai= np.cos(sza)*np.cos(vza)+np.sin(sza)*np.sin(vza)*np.cos(raa)
    kai= np.arccos(cos_kai)
    numerator= (0.5*math.pi-kai)*cos_kai+np.sin(kai)
    denominator= np.cos(sza)+np.cos(vza)
    K_vol= numerator/denominator-0.25*math.pi
    return K_vol

