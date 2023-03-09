#! /usr/bin/env python

# KJ Wagner 2023

# TO DO
# - require that bh1 and bh2 have the same distance (user input specification in jpynb)

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
    
    def __init__(self, bh1, bh2, tref=None, wave=None):
        self.bh1 = bh1
        self.bh2 = bh2
        self.tref = None
        self.wave = None
        
    def setParams(self, m1=10.*lal.MSUN_SI, m2=10.*lal.MSUN_SI, dist=1.e6*lal.PC_SI,
            s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., eccentricity=0.,
            tref=1000000000, fmin=20., approx='TaylorF2', 
            theta=0., phi=0., deltaF=1/16, deltaT=1./4096.            
            ):
        
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
        self.tref = tref
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
        print("make a waveform")
        t_start = int(self.tref) - 150
        t_stop = int(self.tref) + 150
        ifos = ['H1','L1','V1']
        working_dir_full = os.getcwd()
        for ifo in ifos:
            # Create framework to write a signal with params from P
            filename = working_dir_full+"/params.xml.gz"
            instrument = ifo
            xmldoc = utils.load_filename(filename, verbose = True, contenthandler =lalsimutils.cthdler)
            sim_inspiral_table = lsctables.SimInspiralTable.get_table(xmldoc)
            self.wave.copy_sim_inspiral(sim_inspiral_table[0])
            self.wave.taper = lalsimutils.lsu_TAPER_START
            #self.wave.approx = lalsim.GetApproximantFromString(approx)
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

            tstart = int(hoft.epoch)
            duration = int(round(hoft.data.length*hoft.deltaT))
            fname = instrument.replace("1","")+"-fake_strain.gwf"

            lalsimutils.hoft_to_frame_data(fname,channel,hoft)

    def plotWaveform(self):
        print("plot a waveform")





