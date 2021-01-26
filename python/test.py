import numpy as np

import pymultilayer as ml

x = ml.Pymultilayer()

N = np.int(4)
n = np.array([1,2.22,1.45, 1.5],dtype=np.float64)
k = np.array([0,0,0,0],dtype=np.float64)
h = np.array([0,100,80,0],dtype=np.float64)
angle = np.float64(0)
wavelength = np.float64(632)

print ("Angle\tp-result")

for i0 in range(90):
    res = x.calReflectivity(N,n,k,h,np.float64(i0),wavelength,1) # p-polarise
    print (i0, "\t", res)