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
import datetime as datetime

debug = 1

if debug == 1:
    data_path   = 'C:/Users/cschiff/Work/Development/Python/FPI/OpPtCal/OpPtCal/'
    oppt_dir    = '0188_20150916_0316-20150917_0311'
    specfit_dir = 'spec_fit'
else:
    (data_path,oppt_dir,specfit_dir) = sys.argv[1:4]

plt.close("all")
fig0,axes0 = plt.subplots(nrows=4,ncols=4,figsize=(8,8))
plt.figure(fig0.number)
plt.subplots_adjust(hspace=0.,wspace=0.) 

fig1,axes1 = plt.subplots(nrows=2,ncols=3,figsize=(22,16))        
plt.figure(fig1.number)
plt.tight_layout()
    
fig2,axes2 = plt.subplots(nrows=4,ncols=4,figsize=(8,8))
plt.figure(fig2.number)
plt.subplots_adjust(hspace=0.,wspace=0.)   
  

output_dir = '%s/%s/analyzed' % (data_path,oppt_dir)   
path00 = '%s/%s/cdf/OpPt' % (data_path,oppt_dir)
path01 = '%s/%s/cdf/OpPt+50' % (data_path,oppt_dir)


data_paths = [path00,path01]

data = of.create_counts_matrix(data_paths)

counts = data[0]
serial = data[1]
threshs = data[2]
mcps = data[3]
aaxts = data[4]       

if os.path.exists(output_dir) == False:
    os.mkdir(output_dir)
    

if np.sum(counts) > 0:        

    for sensor in ['des','dis']:
  
       oppt_param = np.zeros((len(counts),2,16,5,2),np.float)
        
       op_Gammas = np.zeros(16,np.float)    
       op_Qs = np.zeros(16,np.float)
       pxs = np.arange(0,16)
        
       for head in [0,1]: 
           for cii in range(len(counts)):
               for eii in range(1):           
                    if (sensor == 'des' and serial[cii] > 200) or (sensor == 'dis' and serial[cii] < 200):
                        
                        
                        prefix = 'FM%sH%i' % (str(serial[cii]).zfill(3),head)   
                        print '%s: fitting %s MCP=%iV' % (oppt_dir,prefix,mcps[cii][head])
                        coppt_param = op.plot_tswpsum(fig0,axes0,counts,eii,threshs[cii][head],mcps[cii][head],aaxts[cii],serial[cii],head,cii,output_dir,oppt_dir)
                        oppt_param[cii,head,:,:,:] = coppt_param[cii,head,:,:,:]*1.
                        
                        Qs = oppt_param[cii,head,:,0,0]
                        gammas = oppt_param[cii,head,:,1,0]
    
                        use_px = np.where(np.abs(Qs-np.mean(Qs)) < np.std(Qs)*1.96)[0]                                         
                        
                        p_Gammas = np.polyfit(pxs[use_px],gammas[use_px],4)
                        p_Qs = np.polyfit(pxs[use_px],Qs[use_px],4)
    
                        if cii % 2 == 1:  
                            
                            print '%s: generating new SpecFit for %s' % (oppt_dir,prefix)

                            Qs_op = np.mean(np.polyval(p_Qs,pxs))
                            Qs_opplus = np.mean(np.polyval(op_Qs,pxs))
                            
                            lowMCP = mcps[cii][head]
                            highMCP = mcps[cii-1][head]
                            if Qs_op > Qs_opplus:
                                p_temp = p_Qs*1.
                                p_Qs = op_Qs*1.
                                op_Qs = p_temp*1.
                                
                                tempMCP = lowMCP*1.
                                lowMCP = highMCP*.1
                                highMCP = tempMCP*1.

                            Qs_op = np.mean(np.polyval(p_Qs,pxs))
                            Qs_opplus = np.mean(np.polyval(op_Qs,pxs))
                            
                            print '%s: low Q (%i) @ V=%i, high Q (%i) @ V=%i' % (oppt_dir,Qs_op,lowMCP,Qs_opplus,highMCP)
                                                                    
                            dlogQdV = (np.log10(Qs_opplus)-np.log10(Qs_op))/50.
                            
                            if serial[cii] < 200:
                                dlogQdV = -dlogQdV

                            fname = '%s/specfits/%s/%s%sH%i_%s.csv' % (of.base_dir,specfit_dir,sensor.upper(),str(serial[cii]).zfill(3),head,specfit_dir)
                            ofname = '%s/%s%sH%i_%s_oppt.csv' % (output_dir,sensor.upper(),str(serial[cii]).zfill(3),head,oppt_dir)
                            
                            with open(fname, 'r') as input_file, open(ofname, 'w') as output_file:
                                l = 0
                                
                                for line in input_file:                        
                                    if l == 23:
                                        output_file.write('%i\t%i\t%i\t%i\t%i\n' % (p_Qs[0],p_Qs[1],p_Qs[2],p_Qs[3],p_Qs[4]))
                                    elif l == 24:
                                        output_file.write('%i\t%i\t%i\t%i\t%i\n' % (op_Qs[0],op_Qs[1],op_Qs[2],op_Qs[3],op_Qs[4]))                                                                                                                          
                                    elif l == 26:
                                        output_file.write('0.000e+00\t0.000e+00\t0.000e+00\t0.000e+00\t1.100e+00\n')
                                    elif l == 28:
                                        output_file.write('%.3e\n' % dlogQdV)
                                         
                                    else:
                                        output_file.write(line)
                                    l+=1
                                                             
                        op_Gammas = p_Gammas
                        op_Qs = p_Qs                                
                            
        
       for unit in np.unique(serial):
          if (sensor =='des' and unit >200) or (sensor == 'dis' and unit < 200):
               prefix = 'FM%s' % str(unit).zfill(3)
               print '%s: plotting parameters for %s' % (oppt_dir,prefix)
               
               op.plot_phdparam(fig1,axes1,oppt_param,unit,serial,oppt_dir,specfit_dir,output_dir)                       
               op.plot_tswpnorm(fig2,axes2,counts,threshs,aaxts,unit,serial,oppt_param,oppt_dir,specfit_dir,output_dir)
               
else:
    print '%s: no data available' % oppt_dir

                   
                     