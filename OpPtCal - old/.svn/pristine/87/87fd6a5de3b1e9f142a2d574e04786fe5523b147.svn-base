# -*- coding: utf-8 -*-
"""
############################################################
#             trend_oppt.py
#
#  Mission:  MMS
#
#  Subsytem: Fast Plasma Instrument (FPI)
#
#  Author:   Daniel J. Gershman (daniel.j.gershman@nasa.gov)
#            Conrad Schiff (conrad.schiff-1@nasa.gov)
#
#  Purpose:  Trending script for Operating Point Calibration. 
#
#  Language: Python 2.7
#
#  Packages: matplotlib.pyplot
#            matplotlib.cm
#            matplotlib
#            numpy
#            sys
#            glob
#            os
#            spacepy.pycdf
#            opptcal_func
#            oppt_plot
#
#  History:  1) Initial version (D. Gershman) - 4/9/15
#            2) Update for observatory name in plot title (C. Schiff) - 2/16/16
#
#  $LastChangedBy: cschiff $
#  $LastChangedDate: 2016-02-16 22:36:45 +0000 (Tue, 16 Feb 2016) $
#  $LastChangedRevision: 7035 $
#
#  Copyright 2015 NASA Goddard Space Flight Center
"""
__contact__ = 'Dan Gershman (daniel.j.gershman@nasa.gov)'
import scipy.sparse as sp
from scipy.sparse import linalg
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as matplotlib
import numpy as np
import cPickle
import scipy.special as special
import gc as gc
import sys as sys
import scipy.optimize as optimize
import scipy.special as special
import scipy.io as io
import csv as csv
import glob as glob
import scipy.interpolate as interp
import pdb
import os
from spacepy import pycdf
import scipy as sci
from itertools import groupby
from xlrd import open_workbook
import opptcal_func as of
import oppt_plot as op

def serial_num_to_sc(serial_number):
    sc, bay = 'null', 'null'
    if serial_number == 'DES209':  sc,bay = 'mms1', 'DES0'
    if serial_number == 'DES211':  sc,bay = 'mms1', 'DES1'
    if serial_number == 'DES212':  sc,bay = 'mms1', 'DES2'
    if serial_number == 'DES207':  sc,bay = 'mms1', 'DES3'
    if serial_number == 'DES203':  sc,bay = 'mms2', 'DES0'
    if serial_number == 'DES204':  sc,bay = 'mms2', 'DES1'        
    if serial_number == 'DES201':  sc,bay = 'mms2', 'DES2'
    if serial_number == 'DES202':  sc,bay = 'mms2', 'DES3'
    if serial_number == 'DES206':  sc,bay = 'mms3', 'DES0'
    if serial_number == 'DES216':  sc,bay = 'mms3', 'DES1'
    if serial_number == 'DES213':  sc,bay = 'mms3', 'DES2'
    if serial_number == 'DES214':  sc,bay = 'mms3', 'DES3'
    if serial_number == 'DES215':  sc,bay = 'mms4', 'DES0'
    if serial_number == 'DES210':  sc,bay = 'mms4', 'DES1'
    if serial_number == 'DES208':  sc,bay = 'mms4', 'DES2'
    if serial_number == 'DES205':  sc,bay = 'mms4', 'DES3'        
#
    if serial_number == 'DIS002':  sc,bay = 'mms1', 'DIS0'
    if serial_number == 'DIS005':  sc,bay = 'mms1', 'DIS1'
    if serial_number == 'DIS006':  sc,bay = 'mms1', 'DIS2'
    if serial_number == 'DIS013':  sc,bay = 'mms1', 'DIS3'
    if serial_number == 'DIS014':  sc,bay = 'mms2', 'DIS0'
    if serial_number == 'DIS001':  sc,bay = 'mms2', 'DIS1'        
    if serial_number == 'DIS009':  sc,bay = 'mms2', 'DIS2'
    if serial_number == 'DIS012':  sc,bay = 'mms2', 'DIS3'
    if serial_number == 'DIS008':  sc,bay = 'mms3', 'DIS0'
    if serial_number == 'DIS015':  sc,bay = 'mms3', 'DIS1'
    if serial_number == 'DIS016':  sc,bay = 'mms3', 'DIS2'
    if serial_number == 'DIS007':  sc,bay = 'mms3', 'DIS3'
    if serial_number == 'DIS010':  sc,bay = 'mms4', 'DIS0'
    if serial_number == 'DIS011':  sc,bay = 'mms4', 'DIS1'
    if serial_number == 'DIS003':  sc,bay = 'mms4', 'DIS2'
    if serial_number == 'DIS004':  sc,bay = 'mms4', 'DIS3'
        
    return sc, bay        

macros = [[125,171,211,240,267,302],
          [125,171,211,240,269,303],
          [125,165,170,211,223,240,268,303],
          [125,170,210,240,268,303]]


base_dir = '/fpishare/Gershman/OpPtCal'
plt.close("all")
fig,axes = plt.subplots(nrows=4,ncols=1,figsize=(16,12), dpi= 100, sharex=True)

tr = 4.0e5
rt = 0.5

for obs in [1,2,3,4]:
 
    obs_macros = np.asarray(macros[obs-1])
    for unit in np.hstack( (np.arange(1,17),np.arange(201,217))):
        for head in range(2):
            prefix = '%sH%i' % (str(unit).zfill(3),head)   
    
            if unit > 200:
                sensor = 'des'                
            else:
                sensor = 'dis'

            data_path = '/data/ftp/mms%i/fpi/cal_sandbox/OpPtCal' % obs
            
            
            
            if os.path.exists(data_path) == True:
                dirs = os.listdir(data_path)
    
                files = []
                for cdir in dirs:
                    cfiles = glob.glob('%s/%s/analyzed/*%s_*_oppt.csv' % (data_path,cdir,prefix))
                    files = files + cfiles
                    
        
                if len(files) > 0:
        
                    
                    oppt_data = np.zeros((len(files),2,5,16),np.float)
                    orbits = np.zeros(len(files),np.float)
                    
                    for fii in range(len(files)):
                        
                        fname = files[fii]
                        
                        orbit = int(os.path.basename(fname).split('_')[1])
                        
                        datareader = open(fname,'r')
                        specFit = datareader.readlines()    # Load most recent entry in file  
                        
                        if len(specFit) < 80:
                            continue
                            
                        
                        skiporbit = 0
                        for pii in range(16):      
                            if len(np.array(map(float,specFit[78+pii].split()[2:]))) == 5:
                                oppt_data[fii,0,:,pii] = np.array(map(float,specFit[78+pii].split()[2:]))  
                                oppt_data[fii,1,:,pii] = np.array(map(float,specFit[97+pii].split()[2:]))                              
                                skiporbit = 0
                            else:
                                skiporbit = 1
                                
                        if skiporbit == 1:
                            continue
                        orbits[fii] = orbit  
                        print fname,orbit
                                           
                    
                    if np.max(orbits) == 0:
                        continue
                        
                    axes[0].cla()
                    axes[1].cla()
                    axes[2].cla() 
                    axes[3].cla() 

                    sii = orbits.argsort()
                                        
                    for pii in range(16):
                    
                        Qs00 = oppt_data[sii,0,0,pii]
                        Gs00 = 1.1*(3.*1.1-0.8)/(3.*1.1+0.2)*Qs00
                        SLs00 = special.gammainc(1.1,np.tile(tr,Qs00.size)/Qs00)  
                        Rs00 = oppt_data[sii,0,2,pii]
                        aaxt00 = oppt_data[sii,0,3,pii]
                        
                        if sensor == 'dis':
                            aaxt00 = 2.5e-3*np.ones(aaxt00.size,np.float)

                        
                        XT00 = of.cdf(tr/aaxt00,Qs00,1.1)/of.cdf(tr*np.ones(Qs00.size,np.float),Qs00,1.1)
                        


                        Qs50 = oppt_data[sii,1,0,pii]
                        Gs50 = 1.1*(3.*1.1-0.8)/(3.*1.1+0.2)*Qs50
                        SLs50 = special.gammainc(1.1,np.tile(tr,Qs50.size)/Qs50)  
                        Rs50 = oppt_data[sii,1,2,pii]
                        aaxt50 = oppt_data[sii,1,3,pii]
                        
                        if sensor == 'dis':
                            aaxt50 = 2.5e-3*np.ones(aaxt00.size,np.float)
                        
                        XT50 = of.cdf(tr/aaxt50,Qs50,1.1)/of.cdf(tr*np.ones(Qs50.size,np.float),Qs50,1.1)
                        
                        if pii >= 0 and pii <= 15:
                            gii = np.where(Rs00 < rt)[0]
                            
                            for ii in gii:
                                
                                cR = (rt-Rs00[ii])*0.5				
                                if Gs00[ii] > 0:	
                                                      
                                    axes[0].semilogy(orbits[sii][ii],Gs00[ii],'bs',alpha=cR,markeredgewidth=0.0)
                                    axes[1].semilogy(orbits[sii][ii],SLs00[ii],'bs',alpha=cR,markeredgewidth=0.0)
                                    axes[2].semilogy(orbits[sii][ii],XT00[ii],'bs',alpha=cR,markeredgewidth=0.0)
                                    axes[3].plot(orbits[sii][ii],aaxt00[ii],'bs',alpha=cR,markeredgewidth=0.0)


                            gii = np.where(Rs50 < rt)[0]
                                        
                            if np.max(orbits[sii][gii]) > 0:
                                for ii in gii:
                                    cR = (rt-Rs50[ii])*0.5					
                                    if Gs50[ii] > 0:					
                                        axes[0].semilogy(orbits[sii][ii],Gs50[ii],'rs',alpha=cR,markeredgewidth=0.0)
                                        axes[1].semilogy(orbits[sii][ii],SLs50[ii],'rs',alpha=cR,markeredgewidth=0.0)
                                        axes[2].semilogy(orbits[sii][ii],XT50[ii],'rs',alpha=cR,markeredgewidth=0.0)
                                        axes[3].plot(orbits[sii][ii],aaxt50[ii],'rs',alpha=cR,markeredgewidth=0.0)

                                          
                    if sensor == 'des':
                        axes[0].axhline(2e6,color='k',linestyle='--')                                    
                        axes[1].axhline(0.05,color='k',linestyle='--')        
                        axes[1].axhline(0.1,color='k',linestyle='--')                                
                        axes[2].axhline(0.1,color='k',linestyle='--')     
                        axes[3].axhline(0.022,color='k',linestyle='--')                   
                                  
                        axes[0].set_ylim([5e5,1e7])
                        axes[1].set_ylim([0.01,0.5])
                        axes[2].set_ylim([0.01,0.5])
                        axes[3].set_ylim([0,0.1])
                    else:
                        axes[0].axhline(5e6,color='k',linestyle='--')                                   
                        axes[1].axhline(0.05,color='k',linestyle='--')        
                        axes[1].axhline(0.1,color='k',linestyle='--')                                             
                        axes[2].axhline(0.05,color='k',linestyle='--')     
                        axes[3].axhline(0.0025,color='k',linestyle='--')                   
                                  
                        axes[0].set_ylim([2.5e6,5e7])
                        axes[1].set_ylim([0.01,0.5])
                        axes[2].set_ylim([0.01,0.5])
                        axes[3].set_ylim([0,0.1])                    
    
                    axes[0].set_xlim([np.min(orbits[orbits > 0])-1,np.max(orbits[orbits > 0])+60])
                    axes[3].set_xlabel('Orbit Number')
                    axes[0].set_ylabel('Gain (e-)')
                    axes[1].set_ylabel('Signal Loss')
                    axes[2].set_ylabel('Countrate Crosstalk')
                    axes[3].set_ylabel('Anode-Anode Crosstalk')
                    
                    if unit > 200:
                        serial_num = 'DES'+prefix[0:3]
                    else:
                        serial_num = 'DIS'+prefix[0:3]

                    sc, bay = serial_num_to_sc(serial_num)
                    axes[0].set_title('%s - %s' % (serial_num+prefix[-2:], sc.upper()) )
                    
                    for orb in obs_macros:
                        for aii in range(3):
                            axes[aii].axvline(orb,color='k')
                    
                    ofname = '/fpishare/Gershman/OpPtCal/output_test/%s%s_trend' % (sensor.upper(),prefix)
                    fig.savefig(ofname,dpi=fig.dpi)  
