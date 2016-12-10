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


specfit_dir = 'specfit_000_0000000_0000-0000000_0000' 
output_dir = '%s/output_newrepository' % of.base_dir    


proc_script = open('%s/procOpPt.bat' % of.base_dir,'w')    
for obs in [1,2,3,4]:

    data_path = 'Z:\\mms%i\\fpi\\cal\\OpPtCal' % obs

    if os.path.exists(data_path):
    
        dirs = os.listdir(data_path)
        
        for oppt_dir in dirs:
            
            if '.' not in oppt_dir:
                #print oppt_dir
                orbitnum = int(oppt_dir.split('-')[0].split('_')[0])
                
                if orbitnum < 228:
                    continue
                oppt_files = glob.glob('%s/%s/analyzed/*_oppt.csv' % (data_path,oppt_dir))  

                if len(oppt_files) >= 16:
                    continue    
                else:
                    print '%s: adding to MMS%i processing list' % (oppt_dir,obs)
          
                proc_script.write('python D:\\fpi\\FPIdataProcTools\\trunk\\CalAnalysis\\OpPtCal\\proc_OpPtCal_orbitCDF_v2.py %s %s %s\n' % (data_path,oppt_dir,specfit_dir))
                
#proc_script.write('pause\n')
proc_script.close()
                           
                     