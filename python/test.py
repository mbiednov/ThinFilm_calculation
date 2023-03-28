import numpy as np

from pymultilayer import pymultilayer as ml

x = ml.Pymultilayer()

N = np.int(4)
n = np.array([1,2.22,1.45, 1.5],dtype=np.float64)
k = np.array([0,0,0,0],dtype=np.float64)
h = np.array([0,100,80,0],dtype=np.float64)
angle = np.float64(0)
wavelength = np.float64(632)

print ("Angle\tresult")


resR_p_all = np.array([])
resR_s_all = np.array([])
resT_p_all = np.array([])
resT_s_all = np.array([])
for i0 in range(90):
    resR_p = x.calReflectivity(N,n,k,h,np.float64(i0),wavelength,1) # p-polarise

    resR_s = x.calReflectivity(N,n,k,h,np.float64(i0),wavelength,0) # s-polarise

    resT_p = x.calTransmission(N,n,k,h,np.float64(i0),wavelength,1) # p-polarise
    resT_s = x.calTransmission(N,n,k,h,np.float64(i0),wavelength,0) # p-polarise
    print ("i0: {}, R_s: {},R_p: {}, T_p: {},T_s: {}".format(i0, resR_p, resR_s, resT_p, resT_s))

    resR_p_all = np.append(resR_p_all, resR_p)
    resR_s_all = np.append(resR_s_all, resR_s)
    resT_p_all = np.append(resT_p_all, resT_p)
    resT_s_all = np.append(resT_s_all, resT_s)

import matplotlib.pyplot as plt

plt.figure(dpi=150)
plt.plot(np.linspace(0,90,90),resR_p_all, label='R_p')
plt.plot(np.linspace(0,90,90),resR_s_all, label='R_s')
plt.plot(np.linspace(0,90,90),resT_p_all, label='T_p')
plt.plot(np.linspace(0,90,90),resT_s_all, label='T_s')
plt.legend()
plt.tight_layout(pad=0.5)
plt.savefig("results")
plt.close()


