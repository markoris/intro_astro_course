#! /usr/bin/env python

# KJ Wagner 2023

## TODO : add in zero padding so signal occupies stretch of plot

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

    tref = 1000000000
    
    def __init__(self, bh1, bh2, wave=None):
        self.bh1 = bh1
        self.bh2 = bh2
        self.wave = None
                
    def setParams(self, m1=10.*lal.MSUN_SI, m2=10.*lal.MSUN_SI, dist=1.e6*lal.PC_SI,
            s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., eccentricity=0.,
            tref=1000000000, fmin=20., approx='TaylorF2', 
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
        d_max = self.distance + 100
        thisWave.randomize(dMax=d_max,dMin=d_min,aligned_spin_Q=False,
                           volumetric_spin_prior_Q=False,sMax=1)
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


    def makeWaveform(self):
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

            channel = instrument+":FAKE-STRAIN"

            ## TODO: add back in zero-padding to make signal occupy full xaxis

            tstart = int(hoft.epoch)
            duration = int(round(hoft.data.length*hoft.deltaT))
            fname = instrument.replace("1","")+"-fake_strain.gwf"
                        
            lalsimutils.hoft_to_frame_data(fname,channel,hoft)

    def plotWaveform(self):

        """
        Make the plot of signal strain vs time
        """
        self.makeWaveform()

        # Time series data
        working_dir_full = os.getcwd()
        dataL1 = TimeSeries.read(working_dir_full+'/L-fake_strain.gwf','L1:FAKE-STRAIN')
        dataH1 = TimeSeries.read(working_dir_full+'/H-fake_strain.gwf','H1:FAKE-STRAIN')
        dataV1 = TimeSeries.read(working_dir_full+'/V-fake_strain.gwf','V1:FAKE-STRAIN')

        # Calculates auto-spectral density
        # turns time series into frequency series
        dataL1_spec = dataL1.asd()
        dataH1_spec = dataH1.asd()
        dataV1_spec = dataV1.asd()

        # Create plot using GWPy
        plot = Plot(dataL1) # ,dataH1, dataV1)
        ax = plot.gca()

        ax.set_epoch(999999850)
        ax.set_xlim(999999850+145,999999850+151)
        ax.set_ylim(-1e-22, 1e-22)
        ax.set_title("Mass1 = {}, Mass2 = {}".format(self.bh1.mass, self.bh2.mass))
        ax.set_ylabel('GW Strain')
        #ax.set_ylabel(r'GW strain ASD [strain$/\sqrt{\mathrm{Hz}}$]')

        #plot.show()





