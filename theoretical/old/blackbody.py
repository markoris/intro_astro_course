import matplotlib.pyplot as plt
import numpy as np

## constants
h = 6.626*10**-34 # [ J s]
c = 3.0*10**8 # [m s**-2]
k = 1.38*10**-23 # [J K]
sig = 5.67*10**-8 # [W m**-2 K**-4]
b = 2.897*10**-3 # [m K]

## Stellar types - temps                                                                                                                                       
spec_class = {'O':41000, 'B':31000, 'A':9500, 'F':7240, 'G':5290, 'K':5300, 'M':3850}

## blackbody/planck function
def planck(wav, temp):
    a = (2 * h * c**2) / wav**5
    b = h * c / (wav * k * temp)
    rad = a *  1 / (np.exp(b) - 1.0)
    return rad

## plotting function
def bb_plot(wav, temp):
    if len(temp) == 1:
        plt.plot(wav, planck(wav, temp))
    else:
        for i in range(len(temp)):
            plt.plot(wav, planck(wav, temp[i]))
    plt.xlabel('Wavelength')
    plt.ylabel('Intensity')
    plt.show()

## Wien's Law
def wien(temp):
    return b/temp

## Stefan-Boltazmann
def energy(temp):
    return sig * temp**4
