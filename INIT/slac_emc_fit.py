import numpy as np

def slac_emc_fit(x,A):
    alpha = -0.070+2.189*x-24.667*x**2+145.291*x**3-497.237*x**4+1013.129*x**5-1208.393*x**6+775.767*x**7-205.872*x**8

    C=np.exp(0.017+0.018*np.log(x)+0.005*(np.log(x))**2)

    emc=C*A**alpha
    
    return emc
