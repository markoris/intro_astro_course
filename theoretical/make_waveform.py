#! /usr/bin/env python

# Uses pieces of util_LALWriteFrame.py from RIFT repository.

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

from gwpy.timeseries import TimeSeries
from gwpy.plot import Plot

# User inputs
parser = argparse.ArgumentParser()
parser.add_argument("--mass1", default=20,type=int,help="Mass of BH 1, input by user")
parser.add_argument("--mass2", default=10,type=int,help="Mass of BH 1, input by user")
parser.add_argument("--distance", default=1000,type=int,help="Distance of system away in Mpc, input by user")
opts =  parser.parse_args()

# Priors
n_events=1

mass1 = opts.mass1
mass2 = opts.mass2

# **** might need to be tweaked !!**** need a range such that waveforms appear diff
m_chirp = (mass1*mass2)**(3./5.)*(mass1+mass2)**(-1./5.)
mc_min=int(m_chirp)
mc_max= mc_min+5
mc_range = [mc_min,mc_max]
m_min=1

d_min = opts.distance
d_max = d_min + 100

eta_min=0.2
eta_max=0.249999
eta_range = [eta_min,eta_max]
ecc_min = 0.05
ecc_max = 0.3
chi_max = 1

no_spin=True
precessing_spin=False
aligned_spin=False
volumetric_spin=True
use_eccentric = False
use_matter=False

fix_sky_location=True
fix_event_time=True

fiducial_event_time=1000000000
fiducial_ra =0
fiducial_dec=0

#approx='TaylorF2'
approx='IMRPhenomXHM'
fmin_template=20

ifos = ['H1','L1', 'V1']
channels = {'H1': 'FAKE-STRAIN', 'L1':'FAKE-STRAIN','V1':'FAKE-STRAIN'}
flow= {'H1':20, 'L1':20, 'V1':20}
fmax = 2000
seglen_data = 16
seglen_analysis = 16

# for frame generation
srate_data=4096  
# for analysis
srate_analysis =4096 
use_LI_conditioning=True

### Create the signal parameter dictionary ###

P_list =[]; indx=0
while len(P_list) < n_events:
    P = lalsimutils.ChooseWaveformParams()
    # Randomize (sky location, etc)
    P.randomize(dMax=d_max,dMin=d_min,aligned_spin_Q=aligned_spin,volumetric_spin_prior_Q=volumetric_spin,sMax=chi_max)
    P.tref = fiducial_event_time
    P.fmin=fmin_template
    P.deltaF = 1./seglen_data
    P.deltaT = 1./srate_data
    # sky location
    if fix_sky_location:
        P.theta = fiducial_dec
        P.phi = fiducial_ra
    # some redundancy
    if no_spin:
        P.s1x=P.s1y=P.s1z=0
        P.s2x=P.s2y=P.s2z=0
    elif aligned_spin:
        P.s1x = P.s1y=0
        P.s2x = P.s2y=0
    if use_matter:
        P.lambda1 = np.random.uniform(0,lambda_max)
        P.lambda2 = np.random.uniform(0,lambda_max)
    if use_eccentric:
        P.eccentricity = np.random.uniform(ecc_min, ecc_max)
        P.fecc = 20.0
    P.approx=lalsim.GetApproximantFromString(approx)
    # Uniform in m1 and m2: 
    m1 = np.random.uniform(mc_range[0],mc_range[1]*2)
    m2 = np.random.uniform(m_min,mc_range[1]*1.5)
    m1,m2 = [np.maximum(m1,m2), np.minimum(m1,m2)]
    P.m1 = m1*lal.MSUN_SI
    P.m2 = m2*lal.MSUN_SI
    # ...but downselect in mchirp, eta
    mc_val = P.extract_param('mc')/lal.MSUN_SI
    eta_val = P.extract_param('eta')
    # Check rand vals are within specified ranges after converting inputs
    if mc_val < mc_range[0] or mc_val > mc_range[1]:
        continue
    if eta_val < eta_range[0] or eta_val > eta_range[1]:
        continue

    P_list.append(P)
    indx+=1
# Saves signal parameters in lal usable format
lalsimutils.ChooseWaveformParams_array_to_xml(P_list,"mdc")

### Make the signal frames to plot ###

t_start = int(fiducial_event_time)-150
t_stop = int(fiducial_event_time)+150

working_dir_full = os.getcwd()

mkdir('signal_frames')
for indx in np.arange(n_events):
    os.chdir(working_dir_full)
    print(working_dir_full)
    target_subdir = 'signal_frames/event_{}'.format(indx)
    # Test if directory already exists
    if os.path.exists(target_subdir):
        print(" Signal frames exist for event {}, skipping ".format(indx))
        continue
    else:
        print("Making new directory...")
        mkdir(target_subdir)
    print(" Writing ", indx)
    os.chdir(working_dir_full+"/"+target_subdir)

    for ifo in ifos:
        from ligo.lw import lsctables, table, utils # check all are needed

        # Create framework to write a signal with params from P
        filename = working_dir_full+"/mdc.xml.gz"
        event = indx
        instrument = ifo
        xmldoc = utils.load_filename(filename, verbose = True, contenthandler =lalsimutils.cthdler)
        sim_inspiral_table = lsctables.SimInspiralTable.get_table(xmldoc)
        P.copy_sim_inspiral(sim_inspiral_table[int(event)])
        P.taper = lalsimutils.lsu_TAPER_START
        P.approx = lalsim.GetApproximantFromString(approx)
        P.detector = ifo
        P.print_params()
        
        T_est = lalsimutils.estimateWaveformDuration(P)
        T_est = P.deltaT*lalsimutils.nextPow2(T_est/P.deltaT)
        if T_est < seglen_data:
            T_est = seglen_data
        P.deltaF = 1./T_est
        print(" Duration ", T_est)

        hoft = lalsimutils.hoft(P) 
        
        if seglen_data/hoft.deltaT > hoft.data.length:
            TDlenGoal = int(seglen_data/hoft.deltaT)
            hoft = lal.ResizeREAL8TimeSeries(hoft, 0, TDlenGoal)

        # zero pad some more on either side, to make sure the segment covers start to stop
        if t_start and hoft.epoch > t_start:
            nToAddBefore = int((float(hoft.epoch)-t_start)/hoft.deltaT)
            print(nToAddBefore, hoft.data.length)
            ht = lal.CreateREAL8TimeSeries("Template h(t)", t_start , 0, hoft.deltaT, lalsimutils.lsu_DimensionlessUnit, hoft.data.length+nToAddBefore)
            ht.data.data = np.zeros(ht.data.length)  # clear
            ht.data.data[nToAddBefore:nToAddBefore+hoft.data.length] = hoft.data.data
            hoft = ht

        if t_stop and hoft.epoch+hoft.data.length*hoft.deltaT < t_stop:
            nToAddAtEnd = int( (-(hoft.epoch+hoft.data.length*hoft.deltaT)+t_stop)/hoft.deltaT)
            print("Padding end ", nToAddAtEnd, hoft.data.length)
            hoft = lal.ResizeREAL8TimeSeries(hoft,0, int(hoft.data.length+nToAddAtEnd))

        channel = instrument+":FAKE-STRAIN"

        tstart = int(hoft.epoch)
        duration = int(round(hoft.data.length*hoft.deltaT))
        fname = instrument.replace("1","")+"-fake_strain-"+str(tstart)+"-"+str(duration)+".gwf"

        print("Writing signal with ", hoft.data.length*hoft.deltaT, " to file ", fname)
        lalsimutils.hoft_to_frame_data(fname,channel,hoft)


### Make the plot of signal strain vs time ###

# Time series data
dataL1 = TimeSeries.read(working_dir_full+'/signal_frames/event_0/L-fake_strain-999999850-300.gwf','L1:FAKE-STRAIN')
dataH1 = TimeSeries.read(working_dir_full+'/signal_frames/event_0/H-fake_strain-999999850-300.gwf','H1:FAKE-STRAIN')
dataV1 = TimeSeries.read(working_dir_full+'/signal_frames/event_0/V-fake_strain-999999850-300.gwf','V1:FAKE-STRAIN')

# Calculates auto-spectral density
# turns time series into frequency series
dataL1_spec = dataL1.asd()
dataH1_spec = dataH1.asd()
dataV1_spec = dataV1.asd()

# Create plot using GWPy
plot = Plot(dataL1) #,dataH1)
ax = plot.gca()

ax.set_epoch(999999850)
ax.set_xlim(999999850+145,999999850+151)
ax.set_ylim(-1e-22, 1e-22)
ax.set_title("Mass1 = {}, Mass2 = {}".format(mass1, mass2))
ax.set_ylabel('GW Strain')
#ax.set_ylabel(r'GW strain ASD [strain$/\sqrt{\mathrm{Hz}}$]')

plot.show()
