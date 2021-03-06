# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 15:53:46 2015

@author: fpi_egse
"""
import pdb
import numpy as np
import scipy.special as special
import glob as glob
import sys as sys
import os as os
from spacepy import pycdf
from xlrd import open_workbook

#base_dir = 'D:\\fpi\\FPIdataProcTools\\trunk\\CalAnalysis\\OpPtCal'
#lib_path = 'D:\\fpi\\FPIdataProcTools\\trunk\\CalAnalysis\\OpPtCal'
base_dir = 'c:/Users/cschiff/Work/Development/Python/FPI/OpPtCal/OpPtCal'
lib_path = 'c:/Users/cschiff/Work/Development/Python/FPI/OpPtCal/OpPtCal'

def cdf(trsb,Qb,gammasb):    
    # Gamma Distribution CDF
    gammasb = np.tile(1.1,trsb.size)
    fb = 1. - special.gammainc(gammasb,trsb/Qb)
    
    return fb

def Ns(trsb,lim0crsb,Qsb,gammasb,aaxts):  
    gammasb = np.tile(1.1,trsb.size)
    #cN = lim0crsb*cdf(trs*1.,Qsb*1.,gammasb*1.)
    #cNxt1 = lim0crsb*cdf(trs/aaxts[0],Qsb*1.,gammasb*1.)  
    #cNxt2 = lim0crsb*cdf(trs/aaxts[1],Qsb*1.,gammasb*1.) 

    cN = (1. - special.gammainc(gammasb,trsb/Qsb))*lim0crsb
    cNxt1 = (1. - special.gammainc(gammasb,trsb/Qsb/aaxts[0]))*lim0crsb
    cNxt2 = (1. - special.gammainc(gammasb,trsb/Qsb/aaxts[1]))*lim0crsb

    
    return (cN,cNxt1,cNxt2)
    
def pdf(trs,Q,gammab):
    f = 1./(special.gamma(gammab)*(Q**(gammab))) * trs**(gammab-1.) * np.exp(-trs/Q)
    return f

def Nspdf(trs,lim0cr,Q,gammab,aaxts):
    N = lim0cr*pdf(trs*1.,Q*1.,gammab*1.)
    
    Nxt1 = lim0cr*pdf(trs/aaxts[0],Q*1.,gammab*1.)       
    Nxt2 = lim0cr*pdf(trs/aaxts[1],Q*1.,gammab*1.)     

    return (N,Nxt1,Nxt2)    
    
def N_fitres(inputs,trs,Nd00,px,aaxts,gammab):
    
    (lim0cr,Q) = inputs

    if Q <= 0 or lim0cr <= 0 or gammab <= 0:
        return 1e100

    #(bN,bNxt1,bNxt2) = Ns(trs*1.,np.tile(lim0cr*1.,trs.size),np.tile(Q*1.,trs.size),np.tile(gammab*1.,trs.size),aaxts*1.)
    gammasb  = np.tile(gammab*1.,trs.size)
    Qsb      = np.tile(Q*1.,trs.size)
    lim0crsb = np.tile(lim0cr*1.,trs.size)
    bN       = (1. - special.gammainc(gammasb,trs/Qsb))*lim0crsb
    
    gammasb  = np.tile(gammab*1.,trs.size)
    Qsb      = np.tile(Q*1.,trs.size)
    lim0crsb = np.tile(lim0cr*1.,trs.size)    
    
    bNxt1    = (1. - special.gammainc(gammasb,trs/Qsb/aaxts[0]))*lim0crsb
    
    gammasb  = np.tile(gammab*1.,trs.size)
    Qsb      = np.tile(Q*1.,trs.size)
    lim0crsb = np.tile(lim0cr*1.,trs.size)    
    
    bNxt2    = (1. - special.gammainc(gammasb,trs/Qsb/aaxts[1]))*lim0crsb    
    
    Nfit     = bN*1.
    

    
    if px > 0:
        Nfit += bNxt1
    
    if px > 1:
        Nfit += bNxt2
    
    if px < 14:
        Nfit += bNxt2
        
    if px < 15:
        Nfit += bNxt1
   
  
    fit_res = (np.log10(Nfit)-np.log10(Nd00))/np.log10(Nd00)
    fit_res = fit_res[np.isfinite(fit_res)]
    fit_res = fit_res[np.isnan(fit_res)==0]

    if len(fit_res) == 0:
        return 1e100
        
    return fit_res
    
def get_serial_number(filename):


    st2        = os.path.basename(filename).split('_')
    obs        = st2[0].split('mms')[-1]
    spec       = st2[1].split('s')[-1]

    species    = st2[1].split('s')[0]

    serial_des = np.array([[209,211,212,207],[203,204,201,202],[206,216,213,214],[215,210,208,205]])
    serial_dis = np.array([[2,5,6,13],[14,1,9,12],[8,15,16,7],[10,11,3,4]])

    if species == 'de':
        serial = serial_des[int(obs)-1][int(spec)]
    else:
        serial = serial_dis[int(obs)-1][int(spec)]


    return serial

def threshconver(ser,h,xvalues):
    serial = str(ser)
    head   = str(h)

    if int(serial) < 200:
        sensor = 'DIS'
    else:
        sensor = 'DES'
        
    xlsfile   = lib_path + '/FPIthresholdConversions_%s.txt' % sensor

    serial_st = '%s S/N ' % sensor +serial #+' '+'Head '+head
    head_st   = 'Head '+head
        
    datafile  = np.genfromtxt(xlsfile,skip_header=2)

    thresh_volt = datafile[:,0]
    if int(serial) < 200:
        data = datafile[:,1]
    else:        
        data = datafile[:,(int(serial)-200)*2+int(head)-1]        

    threshs  = np.interp(xvalues,thresh_volt,data)

    return threshs

def create_counts_matrix(dirs):


    cdf_filenames0 = glob.glob(dirs[0] +'/*sigthrsh*.cdf')  # at operating point
    cdf_filenames1 = glob.glob(dirs[1] +'/*sigthrsh*.cdf')  # at operating point+50


    # make sure to get latest version
    cdf_filenames = []
    for fii in range(len(cdf_filenames0)):
        sensor_id = int(os.path.basename(cdf_filenames0[fii]).split('_')[1][-1])
        v_str = (os.path.basename(cdf_filenames0[fii]).split('_')[-1][1:-4]).split('.')
        
        skip = 0
        for gii in range(len(cdf_filenames0)):
            csensor_id = int(os.path.basename(cdf_filenames0[gii]).split('_')[1][-1])
            cv_str = (os.path.basename(cdf_filenames0[gii]).split('_')[-1][1:-4]).split('.')
            
            if sensor_id == csensor_id:
                if int(v_str[0]) < int(cv_str[0]):
                    skip = 1
                   
                    
                if int(v_str[0]) == int(cv_str[0]) and int(v_str[1]) < int(cv_str[1]):
                    
                    skip = 1
                    
                if int(v_str[0]) == int(cv_str[0]) and int(v_str[1]) == int(cv_str[1]) and int(v_str[2]) < int(cv_str[2]):
                    
                    skip = 1
                    
        if skip == 0:
            cdf_filenames = cdf_filenames + [cdf_filenames0[fii]]



    for fii in range(len(cdf_filenames1)):
        sensor_id = int(os.path.basename(cdf_filenames1[fii]).split('_')[1][-1])
        v_str = (os.path.basename(cdf_filenames1[fii]).split('_')[-1][1:-4]).split('.')
        
        skip = 0
        for gii in range(len(cdf_filenames1)):
            csensor_id = int(os.path.basename(cdf_filenames1[gii]).split('_')[1][-1])
            cv_str = (os.path.basename(cdf_filenames1[gii]).split('_')[-1][1:-4]).split('.')
            
            if sensor_id == csensor_id:
                if int(v_str[0]) < int(cv_str[0]):
                    skip = 1
                   
                    
                if int(v_str[0]) == int(cv_str[0]) and int(v_str[1]) < int(cv_str[1]):
                    
                    skip = 1
                    
                if int(v_str[0]) == int(cv_str[0]) and int(v_str[1]) == int(cv_str[1]) and int(v_str[2]) < int(cv_str[2]):
                    
                    skip = 1
                    
        if skip == 0:
            cdf_filenames = cdf_filenames + [cdf_filenames1[fii]]

            
    

    lii = []
    for fii in range(len(cdf_filenames)/2):
        lii.append(fii)
        lii.append(fii+len(cdf_filenames)/2)
    

        
    cdf_filenames = [cdf_filenames[lii[ii]] for ii in range(len(cdf_filenames))]
    

    serial=[]
    loop_index = 0
    counts = np.zeros(shape=(16,8,8,2,16,16))
    threshs = np.zeros(shape=(16,2,16))
    mcps = np.zeros(shape=(16,2))




    #pdb.set_trace()
    ####  ADD THRESH 01 AND MAKE ARRAYS LENTGH 16 [4DES,4DIS, 2MCPV = 16]
        
    for files in cdf_filenames:
        serial.append(get_serial_number(files))

        cdf = pycdf.CDF(files)

        cdf_keys = np.array(cdf.keys()[:])

        # extract the tag names for the data we need to construct the count matrix
        counts_head0_name = cdf_keys[np.char.find(cdf_keys,'counts_Head0') > -1]  #this does not assume a specific mms or dis/des number
        counts_head1_name = cdf_keys[np.char.find(cdf_keys,'counts_Head1') > -1]
        thresh_index_name = cdf_keys[np.char.find(cdf_keys,'thresh_index') > -1]
        it_index_name = cdf_keys[np.char.find(cdf_keys,'iteration_index') > -1]
        eg_index_name = cdf_keys[np.char.find(cdf_keys,'energy_index') > -1]
        pixel_index_name = cdf_keys[np.char.find(cdf_keys,'pixel_index') > -1]
        thresh_volt_name = cdf_keys[np.char.find(cdf_keys,'ThrshV') > -1]
        mcp_h0_name = cdf_keys[np.char.find(cdf_keys,'MCPvlt_Head0') > -1]
        mcp_h1_name = cdf_keys[np.char.find(cdf_keys,'MCPvlt_Head1') > -1]

        

        # use indices defined above to pull out data
        thresh_volt = np.array(cdf[thresh_volt_name[0]])
        counts_h0 = np.array(cdf[counts_head0_name[0]])
        counts_h1 = np.array(cdf[counts_head1_name[0]])
        thresh_index = np.array(cdf[thresh_index_name[0]])
        eg_index = np.array(cdf[eg_index_name[0]])
        pixel_index = np.array(cdf[pixel_index_name[0]])
        it_index = np.array(cdf[it_index_name[0]])
        mcp_h0 = np.array(cdf[mcp_h0_name[0]])
        mcp_h1 = np.array(cdf[mcp_h1_name[0]])
        

        cdf.close()

        # create an array based on the unique values in the parameter (probably use for looping)
        it_uniq = np.unique(it_index)
        pixel_uniq = np.unique(pixel_index)
        thresh_uniq = np.unique(thresh_index)
        eg_uniq = np.unique(eg_index)

        # need to remove double threshV at the first ESA voltage step, this step is not needed once the bug is fixed in the CDF files
        double_thresh_values = np.where(thresh_volt == np.roll(thresh_volt,-1)) #how many elements do we roll this??
        new_thresh_volt = np.delete(thresh_volt[:],double_thresh_values)
        new_mcp_h0 = np.delete(mcp_h0[:],double_thresh_values)
        new_mcp_h1 = np.delete(mcp_h1[:],double_thresh_values)

        mcp_h0_new_shape = new_mcp_h0.reshape(16,8,8)
        mcp_h1_new_shape = new_mcp_h1.reshape(16,8,8)
        thresh_volt_new_shape = new_thresh_volt.reshape(8,8,16)

        # now we need to restructure the data to match input for Dan's code

        for ll in range(0,2,1):

            threshs[loop_index,ll,:] = threshconver(serial[loop_index],ll,thresh_volt_new_shape[0,0,:])

            if ll == 0: mcps[loop_index,ll] = mcp_h0_new_shape[loop_index,0,0]
            if ll == 1: mcps[loop_index,ll] = mcp_h1_new_shape[loop_index,0,0]
           # pdb.set_trace()
            for pp in range(0,16,1):
               
                new_counts_h0 = np.roll(np.delete(counts_h0[:,pp],double_thresh_values),1) # here I remove the double thresh values and roll the array by one to the right, there is still a slight error (talk to george)
                new_counts_h1 = np.roll(np.delete(counts_h1[:,pp],double_thresh_values),1)

                #reshape the 1-D counts array
                counts_h0_new_shape = new_counts_h0.reshape(8,8,16) #[esa volts, iteration,thresh]
                counts_h1_new_shape = new_counts_h1.reshape(8,8,16) 
                # pdb.set_trace()
                for jj in range(len(it_uniq)):
                    for kk in range(len(thresh_uniq)):
                        for ii in range(len(eg_uniq)):

                            if ll == 0:
                                counts[loop_index,jj,ii,ll,pp,kk]=counts_h0_new_shape[jj,ii,kk] #[kk,jj,ii]
                            else:
                                counts[loop_index,jj,ii,ll,pp,kk]=counts_h1_new_shape[jj,ii,kk] #[kk,jj,ii]

        loop_index=loop_index+1
        print 'Loaded %s' % files

        
    aaxts = np.zeros((len(serial),2),np.float)
    for sii in range(len(serial)):
        if serial[sii] > 200:
            aaxts[sii,:] = [0.022,0.007]
        else:
            aaxts[sii,:] = [2.5e-3,1e-6]
    
    
    return counts, serial, threshs, mcps, aaxts 
