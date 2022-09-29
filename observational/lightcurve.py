import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

## This assumes all data files will be same format from database
def load_data(file):
    data = np.loadtxt(file, skiprows=3, usecols=[1,2])
    time = [t[0] for t in data]
    flux = [f[1] for f in data]
    return time, flux

## Obtain folded time series
def fold(time, period):
    time_fold = []
    for i in range(len(time)):
        t_prime = (time[i] - time[0]) %  period
        time_fold.append(t_prime)
    return time_fold

## plot lightcurve
def plot_lc(time, flux):
    plt.plot(time,flux)
    plt.ylim([0.98,1.01])
    plt.xlabel('time')
    plt.ylabel('normalized flux')
    plt.show()

## Plot phase curve with folded time series    
def plot_phase(time, flux):
    plt.plot(time,flux,'k.')
    plt.ylim([0.98,1.01])
    plt.xlabel('phase')
    plt.ylabel('normalized flux')
    plt.show()










