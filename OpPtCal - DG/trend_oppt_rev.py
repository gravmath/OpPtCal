# -*- coding: utf-8 -*-
"""
Created on Thu Apr 09 12:26:14 2015

@author: djgersh
"""

"""
Analyze Operating Point Calibration Data (Commissioning)

Author:  Daniel J. Gershman (daniel.j.gershman@nasa.gov)

Created on 04/6/2015

Copyright 2015 NASA Goddard Space Flight Center

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

base_dir = '/fpishare/Gershman/OpPtCal'
plt.close("all")
fig,axes = plt.subplots(nrows=3,ncols=1,figsize=(16,9),dpi=100,sharex=True)

for obs in [1,2,3,4]:
    for unit in np.hstack( (np.arange(1,17),np.arange(201,217))):
        for head in range(2):
            prefix = '%sH%i' % (str(unit).zfill(3),head)   
    
            data_path = 'y:/mms%i/fpi/cal/OpPtCal' % obs
            
            if os.path.exists(data_path) == True:
                dirs = os.listdir(data_path)
    
                files = []
                for cdir in dirs:
                    cfiles = glob.glob('%s/%s/analyzed/*%s_*_oppt.csv' % (data_path,cdir,prefix))
                    files = files + cfiles
                    
        
                if len(files) > 0:
        
                    
                    oppt_data = np.zeros((len(files),4,4),np.float)
                    
                    
                    for fii in range(len(files)):
                        
                        fname = files[fii]
                        
                        orbit = int(os.path.basename(fname).split('_')[1])
                        
                        #print fname,orbit   
                        datareader = open(fname,'r')
                        specFit = datareader.readlines()    # Load most recent entry in file  
                        
                        (mcphv,mcpb,tr,trv,tau) = np.array(map(float,specFit[12].split()))    
                        p_Qs00 = np.array(map(float,specFit[23].split()))  
                        p_Qs50 = np.array(map(float,specFit[24].split()))  
                        lambdas = np.array(map(float,specFit[31].split()))  # current specfits aren't right       
                        
                        pxs = np.arange(0,16)   
                    
                        Qs00 = np.polyval(p_Qs00,pxs)
                        Qs50 = np.polyval(p_Qs50,pxs)    
                        
                        G00 = Qs00*1.1*(3.*1.1-0.8)/(3.*1.1+0.2)
                        G50 = Qs50*1.1*(3.*1.1-0.8)/(3.*1.1+0.2)
                    
                        SL00 = special.gammainc(1.1,np.tile(tr,16)/Qs00)  
                        SL50 = special.gammainc(1.1,np.tile(tr,16)/Qs50) 
                        
                        XT00 = of.cdf(tr/lambdas[0],Qs00,1.1)/of.cdf(tr,Qs00,1.1)
                        XT50 = of.cdf(tr/lambdas[0],Qs50,1.1)/of.cdf(tr,Qs50,1.1)
                        
                        oppt_data[fii,0,:] = [orbit,np.nanmax(G00),np.nanmax(SL00),np.nanmax(XT00)]
                        oppt_data[fii,1,:] = [orbit,np.nanmin(G00),np.nanmin(SL00),np.nanmin(XT00)] 
                    
                        oppt_data[fii,2,:] = [orbit,np.nanmax(G50),np.nanmax(SL50),np.nanmax(XT50)]
                        oppt_data[fii,3,:] = [orbit,np.nanmin(G50),np.nanmin(SL50),np.nanmin(XT50)] 
                    
                    
                    axes[0].cla()
                    axes[1].cla()
                    axes[2].cla()
                    
                    if unit > 200:
                        serial_num = 'DES'+prefix[0:3]
                    else:
                        serial_num = 'DIS'+prefix[0:3]
                    print serial_num
                    sc, bay = serial_num_to_sc(serial_num)
                    axes[0].set_title('%s - %s' % (serial_num+prefix[-2:], sc.upper()) )

                    
                    sii = oppt_data[:,0,0].argsort()
                    axes[0].semilogy(oppt_data[sii,0,0],oppt_data[sii,0,1],'r.',label='Op')
                    axes[0].semilogy(oppt_data[sii,0,0],oppt_data[sii,1,1],'r.')
                    axes[0].semilogy(oppt_data[sii,0,0],oppt_data[sii,2,1],'b.',label='Op+50')
                    axes[0].semilogy(oppt_data[sii,0,0],oppt_data[sii,3,1],'b.')
                    axes[0].fill_between(oppt_data[sii,0,0],oppt_data[sii,1,1],oppt_data[sii,0,1],facecolor='red',alpha=0.25)
                    axes[0].fill_between(oppt_data[sii,0,0],oppt_data[sii,3,1],oppt_data[sii,2,1],facecolor='blue',alpha=0.25)
        
                    if unit > 200:
                        axes[0].set_ylim([5e5,1e7])
                        axes[0].plot([0,2*oppt_data[-1,0,0]],[2e6,2e6],'k--',linewidth=2)
                    else:
                        axes[0].set_ylim([2.5e6,5e7])
                        axes[0].plot([0,2*oppt_data[-1,0,0]],[5e6,5e6],'k--',linewidth=2)
                        
                    axes[0].set_ylabel('Gain (e-)')
                    axes[0].legend(loc='lower left',fontsize=12)
                    
 
                    xs = oppt_data[sii[-40:],0,0]
                    ys = oppt_data[sii[-40:],0,1]
                    ps = np.polyfit(xs,np.log10(ys),1)                    
                    orbs = np.arange(xs[0],xs[-1]+40.)
#                    axes[0].plot(orbs,10**np.polyval(ps,orbs),'b--')
                    
                    xs = oppt_data[sii[-40:],0,0]
                    ys = oppt_data[sii[-40:],1,1]
                    ps = np.polyfit(xs,np.log10(ys),1)                    
                    orbs = np.arange(xs[0],xs[-1]+40.)
 #                   axes[0].plot(orbs,10**np.polyval(ps,orbs),'b--')                    
                    
                    xs = oppt_data[sii[-40:],0,0]
                    ys = oppt_data[sii[-40:],2,1]
                    ps = np.polyfit(xs,np.log10(ys),1)                    
                    orbs = np.arange(xs[0],xs[-1]+40.)
  #                  axes[0].plot(orbs,10**np.polyval(ps,orbs),'r--')
                    
                    xs = oppt_data[sii[-40:],0,0]
                    ys = oppt_data[sii[-40:],3,1]
                    ps = np.polyfit(xs,np.log10(ys),1)                  
                    orbs = np.arange(xs[0],xs[-1]+40.)
   #                 axes[0].plot(orbs,10**np.polyval(ps,orbs),'r--')       
                   
                                                          
                    axes[1].semilogy(oppt_data[sii,0,0],oppt_data[sii,0,2],'r.')
                    axes[1].semilogy(oppt_data[sii,0,0],oppt_data[sii,1,2],'r.')
                    axes[1].semilogy(oppt_data[sii,0,0],oppt_data[sii,2,2],'b.')
                    axes[1].semilogy(oppt_data[sii,0,0],oppt_data[sii,3,2],'b.')
                    axes[1].fill_between(oppt_data[sii,0,0],oppt_data[sii,1,2],oppt_data[sii,0,2],facecolor='red',alpha=0.25)
                    axes[1].fill_between(oppt_data[sii,0,0],oppt_data[sii,3,2],oppt_data[sii,2,2],facecolor='blue',alpha=0.25)
                    
                    xs = oppt_data[sii[-40:],0,0]
                    ys = oppt_data[sii[-40:],0,2]
                    ps = np.polyfit(xs,np.log10(ys),1)         
                    orbs = np.arange(xs[0],xs[-1]+40.)
#                    axes[1].plot(orbs,10**np.polyval(ps,orbs),'b--')
                    
                    xs = oppt_data[sii[-40:],0,0]
                    ys = oppt_data[sii[-40:],1,2]
                    ps = np.polyfit(xs,np.log10(ys),1)                 
                    orbs = np.arange(xs[0],xs[-1]+40.)
 #                   axes[1].plot(orbs,10**np.polyval(ps,orbs),'b--')                    
                    
                    xs = oppt_data[sii[-40:],0,0]
                    ys = oppt_data[sii[-40:],2,2]
                    ps = np.polyfit(xs,10**np.log10(ys),1)                     
                    orbs = np.arange(xs[0],xs[-1]+40.)
  #                  axes[1].plot(orbs,np.polyval(ps,orbs),'r--')
                    
                    xs = oppt_data[sii[-20:],0,0]
                    ys = oppt_data[sii[-20:],3,2]
                    ps = np.polyfit(xs,np.log10(ys),1)              
                    orbs = np.arange(xs[0],xs[-1]+40.)
   #                 axes[1].plot(orbs,10**np.polyval(ps,orbs),'r--')                    
                                        
                    
                    
                    
                    axes[1].set_ylim([1e-2,0.5])
                    axes[1].plot([0,2*oppt_data[-1,0,0]],[0.1,0.1],'k--',linewidth=2)
                    axes[1].set_ylabel('Signal Loss')
                    
                    axes[2].semilogy(oppt_data[sii,0,0],oppt_data[sii,0,3],'r.')
                    axes[2].semilogy(oppt_data[sii,0,0],oppt_data[sii,1,3],'r.')
                    axes[2].semilogy(oppt_data[sii,0,0],oppt_data[sii,2,3],'b.')
                    axes[2].semilogy(oppt_data[sii,0,0],oppt_data[sii,3,3],'b.')
                    axes[2].fill_between(oppt_data[sii,0,0],oppt_data[sii,1,3],oppt_data[sii,0,3],facecolor='red',alpha=0.25)
                    axes[2].fill_between(oppt_data[sii,0,0],oppt_data[sii,3,3],oppt_data[sii,2,3],facecolor='blue',alpha=0.25)
                    
                    if unit > 200:
                        axes[2].plot([0,2*oppt_data[-1,0,0]],[0.1,0.1],'k--',linewidth=2)
                    else:
                        axes[2].plot([0,2*oppt_data[-1,0,0]],[0.05,0.05],'k--',linewidth=2)
                    axes[2].set_ylim([1e-2,0.5])
                    axes[2].set_ylabel('Crosstalk')
                    axes[2].set_xlim([np.min(oppt_data[:,0,0])-0.5,np.max(oppt_data[:,0,0])+0.5+40])
                    axes[2].set_xlabel('Orbit')



                    xs = oppt_data[sii[-20:],0,0]
                    ys = oppt_data[sii[-20:],0,3]
                    ps = np.polyfit(xs,ys,1)                    
                    orbs = np.arange(xs[0],xs[-1]+40.)
    #                axes[2].plot(orbs,np.polyval(ps,orbs),'b--')
                    
                    xs = oppt_data[sii[-20:],0,0]
                    ys = oppt_data[sii[-20:],1,3]
                    ps = np.polyfit(xs,ys,1)                    
                    orbs = np.arange(xs[0],xs[-1]+40.)
     #               axes[2].plot(orbs,np.polyval(ps,orbs),'b--')                    
                    
                    xs = oppt_data[sii[-20:],0,0]
                    ys = oppt_data[sii[-20:],2,3]
                    ps = np.polyfit(xs,ys,1)                    
                    orbs = np.arange(xs[0],xs[-1]+40.)
      #              axes[2].plot(orbs,np.polyval(ps,orbs),'r--')
                    
                    xs = oppt_data[sii[-20:],0,0]
                    ys = oppt_data[sii[-20:],3,3]
                    ps = np.polyfit(xs,ys,1)                    
                    orbs = np.arange(xs[0],xs[-1]+40.)
       #             axes[2].plot(orbs,np.polyval(ps,orbs),'r--')                    
                                        
                    


                    
                    if unit > 200:
                        sensor = 'DES'
                    else:
                        sensor = 'DIS'
                    ofname = 'c:/users/cschiff/work/fpi/OpPtCal_trending_plots/%s%s_trend' % (sensor,prefix)
                    print ofname
                    fig.savefig(ofname,dpi=fig.dpi)
                    sys.exit()                    
