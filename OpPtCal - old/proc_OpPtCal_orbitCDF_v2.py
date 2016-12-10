# -*- coding: utf-8 -*-
"""
############################################################
#             proc_OpPtCal_orbitCDF_v2.py
#
#  Mission:  MMS
#
#  Subsytem: Fast Plasma Instrument (FPI)
#
#  Author:   Daniel J. Gershman (daniel.j.gershman@nasa.gov)
#            Conrad Schiff (conrad.schiff-1@nasa.gov)
#
#  Purpose:  Main workhorse for Operating Point Calibration. 
#            This script fits an individual days worth of 
#            electron and ion counts to characterize and 
#            calibrate each individual spectrometer.
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
#  History:  1) Initial version (D. Gershman) - 4/6/15
#            2) Update for flexi-paths (C. Schiff) - 1/06/2016
#
#  $LastChangedBy: dgershma $
#  $LastChangedDate: 2016-06-22 11:10:36 -0400 (Wed, 22 Jun 2016) $
#  $LastChangedRevision: 7159 $
#
#  Copyright 2015 NASA Goddard Space Flight Center
"""

__contact__ = 'Dan Gershman (daniel.j.gershman@nasa.gov)'

import opptcal_func as of
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys as sys
import os as os
import oppt_plot as op
import opptcal_paths as paths
import pdb

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


output_dir = '%s/%s/analyzed'    % (data_path,oppt_dir)   
path00     = '%s/%s/cdf/OpPt'    % (data_path,oppt_dir)
path01     = '%s/%s/cdf/OpPt+50' % (data_path,oppt_dir)


data_paths = [path00,path01]
print 'Loading data'
data = of.create_counts_matrix(data_paths)
print 'Data loaded'
counts = data[0]
serial = data[1]
threshs = data[2]
mcps = data[3]
aaxts = data[4]       

if os.path.exists(output_dir) == False:
    os.mkdir(output_dir)
    

if np.sum(counts) > 0:        

    for sensor in ['des','dis']:
  
       #changed the size of oppt_param to allow index 6 to have total counts - CS
       oppt_param = np.zeros((len(counts),2,16,8,2),np.float)
        
       op_Gammas = np.zeros(16,np.float)    
       op_Qs = np.zeros(16,np.float)
       pxs = np.arange(0,16)
        
       for head in [0,1]: 
           for cii in range(len(counts)):
               for eii in range(1):           
                    if (sensor == 'des' and serial[cii] > 200) or (sensor == 'dis' and serial[cii] < 200):
                        
                        
                        prefix = 'FM%sH%i' % (str(serial[cii]).zfill(3),head)   
                        print '%s: fitting %s MCP=%iV' % (oppt_dir,prefix,mcps[cii][head])

                        #pdb.set_trace()                        
                        coppt_param = op.plot_tswpsum(fig0,axes0,counts,eii,threshs[cii][head],mcps[cii][head],aaxts[cii],serial[cii],head,cii,output_dir,oppt_dir)
                        oppt_param[cii,head,:,:,:] = coppt_param[cii,head,:,:,:]*1.
                        
                        #aligned these assignments and unpacked the lim0crs to nonlinear fit debug purposes - CS, 6/17/16
                        Qs       = oppt_param[cii,head,:,0,0]
                        gammas   = oppt_param[cii,head,:,1,0]
                        lim0crs  = oppt_param[cii,head,:,2,0]
                        R2s      = oppt_param[cii,head,:,5,0]
                        

                        use_px   = np.where(np.abs(Qs-np.mean(Qs)) < np.std(Qs)*1.96)[0]                                         
                        
                        p_Gammas = np.polyfit(pxs[use_px],gammas[use_px],4)
                        p_Qs     = np.polyfit(pxs[use_px],Qs[use_px],4)
    
                        if cii % 2 == 1:  
                            
                            print '%s: generating new SpecFit for %s' % (oppt_dir,prefix)

                            Qs_op = np.mean(np.polyval(p_Qs,pxs))
                            Qs_opplus = np.mean(np.polyval(op_Qs,pxs))
                            
                            lowMCP          = mcps[cii][head]
                            highMCP         = mcps[cii-1][head]
                            specfit_lowMCP  = lowMCP
                            specfit_highMCP = highMCP
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

                            fname = '%s/%s/%s%sH%i_%s.csv' % (paths.OpPt_specfit_path,specfit_dir,sensor.upper(),str(serial[cii]).zfill(3),head,specfit_dir)
                            ofname = '%s/%s%sH%i_%s_oppt.csv' % (output_dir,sensor.upper(),str(serial[cii]).zfill(3),head,oppt_dir)
                            
                            with open(fname, 'rb') as input_file, open(ofname, 'wb') as output_file:
                                l = 0
                                for line in input_file:         
                                    if l == 12:
                                        #unpack the line
                                        if abs(specfit_lowMCP) < abs(specfit_highMCP):
                                            oppt_volt = specfit_lowMCP
                                        else:
                                            oppt_volt = specfit_highMCP
                                        ul    = line.split('\t')
                                        if sensor == 'des':
                                            ul[0] = str(int(oppt_volt*10.0/12.0))
                                        elif sensor == 'dis':
                                            ul[0] = str(int(oppt_volt+110))
                                        ul[1] = str(int(oppt_volt))
                                        output_file.write('%s\t%s\t%s\t%s\t%s' % (ul[0],ul[1],ul[2],ul[3],ul[4]))
                                    elif l == 23:
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
                                    
                                #added snippet to append the results of the nonlinear fit to the bottom
                                #of the specfit file for debugging purposes - CS, 6/17-20/16
                                Qs_OpPt             = oppt_param[cii-1,head,:,0,0]
                                lim0crs_OpPt        = oppt_param[cii-1,head,:,2,0]
                                R2s_OpPt            = oppt_param[cii-1,head,:,5,0]
                                total_counts_OpPt   = oppt_param[cii-1,head,:,6,0]
                                aaxt_OpPt           = oppt_param[cii-1,head,:,7,0]
                                
                                Qs_OpPt50           = oppt_param[cii,head,:,0,0]
                                lim0crs_OpPt50      = oppt_param[cii,head,:,2,0]
                                R2s_OpPt50          = oppt_param[cii,head,:,5,0]
                                total_counts_OpPt50 = oppt_param[cii,head,:,6,0]
                                aaxt_OpPt50         = oppt_param[cii,head,:,7,0]
                                
                                output_file.write('Non-linear fit parameters for the operating point:\n')
                                output_file.write('               Qs          lim0crs         R2s       AAXT       total counts\n')
                                for pixel in range(0,16):
                                    output_file.write('Pixel %02d:  %e  %e  %e  %e      %d\n' % (pixel, Qs_OpPt[pixel], lim0crs_OpPt[pixel], R2s_OpPt[pixel],aaxt_OpPt[pixel], total_counts_OpPt[pixel]) )

                                output_file.write('\n')
                                if sensor == 'des':
                                    output_file.write('Non-linear fit parameters for the operating point + 50 Volts:\n')
                                if sensor == 'dis':
                                    output_file.write('Non-linear fit parameters for the operating point - 50 Volts:\n')
                                output_file.write('               Qs          lim0crs         R2s       AAXT       total counts\n')                                    
                                for pixel in range(0,16):
                                    output_file.write('Pixel %02d:  %e  %e  %e  %e      %d\n' % (pixel, Qs_OpPt50[pixel], lim0crs_OpPt50[pixel], R2s_OpPt50[pixel],aaxt_OpPt50[pixel], total_counts_OpPt50[pixel]) )
                                    
                                output_file.close()
                            
                                                             
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

                   
                     