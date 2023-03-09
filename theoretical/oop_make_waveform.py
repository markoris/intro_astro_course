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
    
    def __init__(self, bh1, bh2):
        self.bh1 = bh1
        self.bh2 = bh2
        
    def setParams(self, m1=10.*lal.MSUN_SI, m2=10.*lal.MSUN_SI, dist=1.e6*lal.PC_SI,
            s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., eccentricity=0.,
            tref=0, fmin=40., approx='TaylorF2', 
            theta=0., phi=0., deltaF=None, deltaT=1./4096.            
            ):
        # m1/m2/dist can be modified by user via blackHole object
        self.m1 = self.bh1.mass
        self.m2 = self.bh2.mass
        self.distance = self.bh1.distance # require that bh1/bh2 have the same distance !
        # In case we want spin/ecc later
        #self.s1x = s1x
        #self.s1y = s1y
        #self.s1z = s1z
        #self.s2x = s2x
        #self.s2y = s2y
        #self.s2z = s2z
        #self.eccentricity = eccentricity
        # params relevant for waveform creation
        #self.tref = tref
        #self.fmin = fmin
        #self.approx = approx   # Using IMRPhenomXHM or TaylorF2 for now
        #self.theta = theta     # declination
        #self.phi = phi         # right ascension 
        #self.deltaF = deltaF
        #self.deltaT = deltaT

        thisWave = lalsimutils.ChooseWaveformParams()
        d_min = self.distance
        d_max = self.distance + 100
        thisWave.randomize(dMax=d_max,dMin=d_min,aligned_spin_Q=False,volumetric_spin_prior_Q=False,sMax=1)
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

        # get mchirp/eta from input masses - not used in orig ??
        #mc_val = thisWave.extract_param('mc')/lal.MSUN_SI
        #eta_val = thisWave.extract_param('eta')

        thisWave = [thisWave]
        # Saves signal parameters in lal usable format
        return lalsimutils.ChooseWaveformParams_array_to_xml(thisWave,"params")


    def makeWaveform(self):
        print("make a waveform")

    def plotWaveform(self):
        print("plot a waveform")





