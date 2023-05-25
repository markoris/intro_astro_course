#! /usr/bin/env python

# KJ Wagner 2023

import numpy as np
import argparse
import os
import sys

import lal
import lalsimulation as lalsim
from RIFT.misc.dag_utils import mkdir
import RIFT.lalsimutils as lalsimutils
from RIFT.misc.dag_utils import which
lalapps_path2cache = which('lal_path2cache')
from ligo.lw import lsctables, table, utils

from gwpy.timeseries import TimeSeries
from gwpy.plot import Plot

from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz
import h5py
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

# LIGO-specific readligo.py 
import readligo as rl

class blackHole:
    
    """
    Creates black hole object
    mass in solar units, distance in Mpc
    
    """
    
    def __init__(self, mass, distance):
        self.mass = mass
        self.distance = distance

class binaryBH:
    
    """
    Uses black hole objects to create BBH coalescence GW signals
    """

    # Change tref to be time of GW150914 for later ease
    tref = 1126259462.422-150
    
    def __init__(self, bh1, bh2, wave=None):
        self.bh1 = bh1
        self.bh2 = bh2
        self.wave = None
                
    def setParams(self, m1=10.*lal.MSUN_SI, m2=10.*lal.MSUN_SI, dist=1.e6*lal.PC_SI,
            s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., eccentricity=0.,
            tref=tref, fmin=18., approx='TaylorF2', 
            theta=0., phi=0., deltaF=1/16, deltaT=1./4096.            
            ):

        if os.path.isfile('params.xml.gz'):
            os.remove('params.xml.gz')
        
        # m1/m2/dist can be modified by user via blackHole object
        self.m1 = self.bh1.mass
        self.m2 = self.bh2.mass
        self.distance = self.bh1.distance
        
        # Use lalsimutils from RIFT to create waveform params xml
        thisWave = lalsimutils.ChooseWaveformParams()
        d_min = self.distance
        d_max = self.distance + 10
        thisWave.randomize(dMax=d_max,dMin=d_min,aligned_spin_Q=False,
                           volumetric_spin_prior_Q=False,sMax=1,
                           default_inclination=0.6407159166115236,
                           default_phase=2.5830560505242754,
                           default_polarization=2.078943813888743)
        thisWave.m1 = self.m1*lal.MSUN_SI
        thisWave.m2 = self.m2*lal.MSUN_SI
        
        # In case we want spin/ecc later
        thisWave.s1x = s1x
        thisWave.s1y = s1y
        thisWave.s1z = s1z
        thisWave.s2x = s2x
        thisWave.s2y = s2y
        thisWave.s2z = s2z
        thisWave.eccentricity = eccentricity
        
        # params relevant for waveform creation
        thisWave.tref = tref
        thisWave.fmin = fmin
        thisWave.approx = approx   # Using IMRPhenomXHM or TaylorF2 for now
        thisWave.theta = theta     # declination
        thisWave.phi = phi         # right ascension 
        thisWave.deltaF = deltaF
        thisWave.deltaT = deltaT
        thisWave.approx = lalsim.GetApproximantFromString(thisWave.approx)
        
        self.wave = thisWave
        thisWave = [thisWave]
        # Saves signal parameters in lal usable format
        return lalsimutils.ChooseWaveformParams_array_to_xml(thisWave,"params")


    def makeFrames(self):
        """
        As in pp_RIFT and util_LALWriteFrame
        Create waveform table and write to gwf frame file
        Zero pad to occupy length of axes
        """
        
        t_start = int(binaryBH.tref) - 150
        t_stop = int(binaryBH.tref) + 150
        ifos = ['H1','L1','V1']
        working_dir_full = os.getcwd()
        
        for ifo in ifos:

            if os.path.isfile(ifo+'-fake_strain.gwf'):
                os.remove(working_dir_full+'/'+ifo+'-fake_strain.gwf')
            
            # Create framework to write a signal with params from P
            filename = working_dir_full+"/params.xml.gz"
            instrument = ifo
            xmldoc = utils.load_filename(filename, verbose = True, contenthandler =lalsimutils.cthdler)
            sim_inspiral_table = lsctables.SimInspiralTable.get_table(xmldoc)
            self.wave.copy_sim_inspiral(sim_inspiral_table[0])
            self.wave.taper = lalsimutils.lsu_TAPER_START
            self.wave.detector = ifo
    
            # duration based on Newtonian inspiral [fmin, inf]
            T_est = lalsimutils.estimateWaveformDuration(self.wave)
            T_est = self.wave.deltaT*lalsimutils.nextPow2(T_est/self.wave.deltaT)
            if T_est < 16:
                T_est = 16

            self.wave.deltaF = 1./T_est

            # generate TD waveform based on params
            hoft = lalsimutils.hoft(self.wave)

            # Problem where tstart earlier than frame start
            # Solution as in util_LALWriteFrame
            if hoft.epoch > t_start:
                padStart = int((float(hoft.epoch)-t_start)/hoft.deltaT)
                # Remake hoft with padding
                pad_hoft = lal.CreateREAL8TimeSeries("Template h(t)", t_start , 0, hoft.deltaT,
                                               lalsimutils.lsu_DimensionlessUnit,
                                               hoft.data.length+padStart)
                pad_hoft.data.data = np.zeros(pad_hoft.data.length)
                pad_hoft.data.data[padStart:padStart+hoft.data.length] = hoft.data.data
                hoft = pad_hoft

            channel = instrument+":FAKE-STRAIN"

            tstart = int(hoft.epoch)
            duration = int(round(hoft.data.length*hoft.deltaT))
            fname = instrument.replace("1","")+"-fake_strain.gwf"
                        
            lalsimutils.hoft_to_frame_data(fname,channel,hoft)

    def getWaveform(self):
        """
        Save access to waveform later 
        Enables different plotting libs to be used
        """
        
        self.makeFrames()
        
        # Time series data
        # Save all in case want to explain H1/L1/V detectors
        # Particularly related to later plots with both detectors...
        working_dir_full = os.getcwd()
        dataL1 = TimeSeries.read(working_dir_full+'/L-fake_strain.gwf','L1:FAKE-STRAIN')
        dataH1 = TimeSeries.read(working_dir_full+'/H-fake_strain.gwf','H1:FAKE-STRAIN')
        dataV1 = TimeSeries.read(working_dir_full+'/V-fake_strain.gwf','V1:FAKE-STRAIN')

        # For speed, going to keep only H strain for later comparison
        dataH1.write('H-fake_strain.txt')

        return dataH1, dataL1, dataV1

    
def plotWaveform(binary, data):
    """
    Make the plot of signal strain vs time
    """

    # Create plot using GWPy
    plot = Plot(data, figsize=(12,8))
    ax = plot.gca()
    epoch = 1126259462.422-150
    ax.set_epoch(epoch)
    
    ax.set_xlim(epoch-0.15,epoch+0.025)
    ax.set_ylim(-1e-21, 1e-21)
    ax.set_title("Mass1 = {}, Mass2 = {}".format(binary.bh1.mass, binary.bh2.mass))
    ax.set_ylabel('GW Strain')

    
def loadLIGOdata():
    """
    Load given LIGO frames for GW150914
    Using readligo and gwosc tutorial
    https://gwosc.org/s/events/GW150914/GW150914_tutorial.html
    """
    
    fn_H1 = 'H-H1_LOSC_4_V1-1126259446-32.hdf5'
    strain_H1, time_H1, chan_dict_H1 = rl.loaddata(fn_H1, 'H1')

    fn_L1 = 'L-L1_LOSC_4_V1-1126259446-32.hdf5'
    strain_L1, time_L1, chan_dict_L1 = rl.loaddata(fn_L1, 'L1')

    return  strain_H1, time_H1, chan_dict_H1, strain_L1, time_L1, chan_dict_L1

        

def plotLIGOdata(ligo_data):
    """
    Plot strain data near time of event
    gwosc tutorial
    https://gwosc.org/s/events/GW150914/GW150914_tutorial.html
    """

    strain_H1, time_H1, chan_dict_H1, strain_L1, time_L1, chan_dict_L1 = ligo_data
    # both H1 and L1 will have the same time vector, so:
    time = time_H1
    # the time sample interval (uniformly sampled)
    dt = time[1] - time[0]

    # plot +- 5 seconds around the event:
    tevent = 1126259462.422         # Mon Sep 14 09:50:45 GMT 2015 
    deltat = 5.                     # seconds around the event
    # index into the strain time series for this time interval:
    indxt = np.where((time_H1 >= tevent-deltat) & (time_H1 < tevent+deltat))

    plt.figure(figsize=(12,8))
    plt.plot(time_H1[indxt]-tevent,strain_H1[indxt],'r',label='H1 strain')
    plt.plot(time_L1[indxt]-tevent,strain_L1[indxt],'g',label='L1 strain')
    plt.xlabel('time (s) since '+str(tevent))
    plt.ylabel('strain')
    plt.legend(loc='lower right')
    plt.title('Advanced LIGO strain data near GW150914');
    

def whiten(strain, interp_psd, dt):
    """
    Function to whiten data from gwosc tutorial
    whitening: transform to freq domain, divide by asd, then transform back, 
    https://gwosc.org/s/events/GW150914/GW150914_tutorial.html
    """
    
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)

   
    hf = np.fft.rfft(strain)
    white_hf = hf / (np.sqrt(interp_psd(freqs) /dt/2.))
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht


def return_whiten(ligo_data):
    """
    Return whitened strain for data and template
    Should probably generalize for template
    https://gwosc.org/s/events/GW150914/GW150914_tutorial.html
    """

    fs = 4096
    NFFT = 1*fs
    fmin = 10
    fmax = 2000
    
    strain_H1, time_H1, chan_dict_H1, strain_L1, time_L1, chan_dict_L1 = ligo_data
    wave_time, fake_H1 = np.genfromtxt('H-fake_strain.txt').transpose()
     
    time = time_H1
    dt = time[1] - time[0]
    
    Pxx_H1, freqs = mlab.psd(strain_H1, Fs = fs, NFFT = NFFT)
    Pxx_L1, freqs = mlab.psd(strain_L1, Fs = fs, NFFT = NFFT)
    
    # Interpolations of the ASDs computed above for whitening
    psd_H1 = interp1d(freqs, Pxx_H1)
    psd_L1 = interp1d(freqs, Pxx_L1)
    
    strain_H1_whiten = whiten(strain_H1,psd_H1,dt)
    strain_L1_whiten = whiten(strain_L1,psd_L1,dt)
    strain_wave_whiten = whiten(fake_H1,psd_H1,dt)

    # We need to suppress the high frequencies with some bandpassing:
    bb, ab = butter(4, [20.*2./fs, 300.*2./fs], btype='band')
    strain_H1_whitenbp = filtfilt(bb, ab, strain_H1_whiten)
    strain_L1_whitenbp = filtfilt(bb, ab, strain_L1_whiten)
    strain_wave_whitenbp =  filtfilt(bb, ab, strain_wave_whiten)

    return time, strain_H1_whitenbp, strain_L1_whitenbp, strain_wave_whitenbp


def plot_whiten(ligo_data, fit_template=False, binary=None):
    """
    Plot whitened and bandpassed data near GW150914
    Note new units for whitened data (sigma from mean vs strain units)
    https://gwosc.org/s/events/GW150914/GW150914_tutorial.html
    """

    time, strain_H1_whitenbp, strain_L1_whitenbp, strain_wave_whitenbp  = return_whiten(ligo_data)
    fs = 4096
    tevent = 1126259462.422
    strain_L1_shift = -np.roll(strain_L1_whitenbp,int(0.007*fs))
    
    plt.figure(figsize=(12,8))
    plt.plot(time-tevent,strain_H1_whitenbp,'r',label='Real H1 strain',alpha=0.75)
    plt.plot(time-tevent,strain_L1_shift,'g',label='Real L1 strain',alpha=0.75)

    # Option to plot whitened/bp data with and without waveform template
    if fit_template:
        if not binary:
            wave_time, wave_H1 = np.genfromtxt('H-fake_strain.txt').transpose()
            plt.plot(wave_time-tevent+150+0.018,strain_wave_whitenbp,'k',label='IMRPhenomXHM')
        else:
            wave_time, wave_H1 = np.genfromtxt('H-fake_strain.txt').transpose()
            plt.plot(wave_time-tevent+150+0.018,strain_wave_whitenbp,'k',label='IMRPhenomXHM')
            print("Mass1 = {}, Mass2 = {}, Distance = {}".format(binary.bh1.mass, binary.bh2.mass, binary.distance))
    
    plt.xlim([-0.1,0.05])
    plt.ylim([-4,4])
    plt.xlabel('time (s) since '+str(tevent))
    plt.ylabel('whitented strain')
    plt.legend(loc='lower left')
    plt.title('Advanced LIGO whitened strain data near GW150914')






