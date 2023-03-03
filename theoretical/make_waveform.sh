#!/bin/bash

if [ -d "signal_frames" ]
then
    rm -r signal_frames
    rm mdc.xml.gz
fi

echo "You will be asked for some numerical input variables. If you do not wish to supply an argument, simply hit enter when prompted."

echo "What is the mass of the first black hole?"
read mass1
mass1=$mass1

echo "What is the mass of the other black hole?"
read mass2
mass2=$mass2

echo "What is the distance to the system in Mpc"
read dist
distance=$dist

if [ -z "$mass1" ] && [ -z "$mass2" ] && [ -z "$distance" ]
then
    echo "Using defaults."
    python make_waveform.py
elif [ ! -z "$mass1" ] && [ ! -z "$mass2" ] && [ -z "$distance" ]
then
    echo "Using inputs m1=$mass1 and m2=$mass2"
    python make_waveform.py --mass1 $mass1 --mass2 $mass2
elif [ ! -z "$mass1" ] && [ -z "$mass2"]
then
    echo "If you are going to supply mass values, please specify both."
    exit 1
else
    python make_waveform.py --mass1 $mass1 --mass2 $mass2 --distance $distance
fi
