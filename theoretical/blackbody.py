import matplotlib.pyplot as plt
import numpy as np

## constants
h = 6.626*10**-34
c = 3.0*10**8
k = 1.38*10**-23

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
