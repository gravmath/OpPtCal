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



def cdf(trs,Q,gamma):
    
    # Gamma Distribution CDF

    f = 1. - special.gammainc(gamma,trs/Q)

    return f

def Ns(trs,lim0cr,Q,gamma,aaxts):  

    N = lim0cr*cdf(trs,Q,gamma)
    
    Nxt1 = lim0cr*cdf(trs/aaxts[0],Q,gamma)   
    Nxt2 = lim0cr*cdf(trs/aaxts[1],Q,gamma) 

    return (N,Nxt1,Nxt2)
    
    
def pdf(trs,Q,gamma):
    f = 1./(special.gamma(gamma)*(Q**(gamma))) * trs**(gamma-1.) * np.exp(-trs/Q)
    return f

def Nspdf(trs,lim0cr,Q,gamma,aaxts):
    N = lim0cr*pdf(trs,Q,gamma)
    
    Nxt1 = lim0cr*pdf(trs/aaxts[0],Q,gamma)   
    Nxt2 = lim0cr*pdf(trs/aaxts[1],Q,gamma) 

    return (N,Nxt1,Nxt2)    
    

def N_fitres(inputs,trs,Nd00,px,aaxts,gamma):
    
    (lim0cr,Q) = inputs

    if Q < 0 or lim0cr < 0:
        return 1e100
        

    (N,Nxt1,Nxt2) = Ns(trs,np.tile(lim0cr,trs.size),np.tile(Q,trs.size),np.tile(gamma,trs.size),aaxts)

    Nfit = N
    
    if px > 0:
        Nfit += Nxt1
    
    if px > 1:
        Nfit += Nxt2
    
    if px < 14:
        Nfit += Nxt2
        
    if px < 15:
        Nfit += Nxt1
    

    fit_res = (np.log10(Nfit)-np.log10(Nd00))/np.log10(Nd00)
    fit_res = fit_res[np.isfinite(fit_res)]
        
    return fit_res


     

def plot_tswpsum(fig,axes,cs,eii,trs,ms,aaxt,unit,head,cii,tag):
    

    oppt_param = np.zeros((len(cs),2,16,5,2),np.float)

    prefix = 'FM%sH%i' % (str(unit).zfill(3),head)   
    
    
    print '%i, Plotting (RAW) %s @ MCPB=%i, %i' % (cii,prefix,ms,eii)
    
                
    #ofname = 'oppt_img/total_oppttswp_%s_%i_%i_%s' % (prefix,ms,eii,tag)
    ofname = '/Users/gbclark1/Desktop/oppt_img/total_oppttswp_%s_%i_%i_%s' % (prefix,ms,eii,tag)
    if eii > 0:
        return oppt_param
    else:
        counts_sub = np.sum(np.sum(cs[cii][:,:,head,:,:]*1.,axis=0),axis=0)
        


    counts_sub[counts_sub == 1] = 0.    


    Nmax = np.max(counts_sub)*2.
    Nmin = 0.5

    l10Max = 10**np.floor(np.log10(Nmax))                    

    trsf = np.logspace(3,8,100)
        
    Nmax = np.max(counts_sub)*2.
    Nmin = np.min(counts_sub[counts_sub > 0])*0.5
    
    l10Max = 10**np.floor(np.log10(Nmax))                    

    R2 = np.zeros((16,3),np.float)
    minr2s = np.ones(16,np.float)*100.
    l0s = np.zeros((16,3,2),np.float)
    cols = ['r','b','g']
    if Nmax > 0:                

        for px in range(16):
                
            pxx = np.floor((px)/4)
            pxy = px % 4

            Nd00 = counts_sub[:,px] #GBC changed this code perhaps because the counts cii,px indicies were transposed

            axes[pxx][pxy].cla()
            
            if(len(np.where(Nd00 > 0)[0]) > 3):
                p1 = Nd00[0:5]
                p2 = Nd00[5:10]
                p3 = Nd00[10:15]
                
                mms = [np.mean(p1[p1 > 0]),np.mean(p2[p2 > 0]),np.mean(p3[p3 > 0])]                             
                
                axes[pxx][pxy].loglog(trs,Nd00,'k.')
                for ii in range(3):
                    
                    gam = 1.1
                    
                    lim0cr = mms[ii]
                    cNd00 = Nd00*1.
                    cNd00[cNd00 > lim0cr*3.] = 0.

                    g_xt = np.where(Nd00[Nd00 > 0] > lim0cr*2.)[0]
                    g_sig = np.where(Nd00[Nd00 > 0] < lim0cr*0.5)[0]
                                  
                    if len(g_sig) > 0:
                      
                        g_sig = trs[g_sig[0]]
                        g_xt = g_sig*aaxt[0]
                        
                    elif len(g_sig) == 0 and len(g_xt) > 0:
                       
                        g_xt = trs[g_xt[-1]]
                        g_sig = g_xt/aaxt[0]
                        
                    else:
                        if ii > 0:    
                            g_xt = trs[0]
                            g_sig = g_xt/aaxt[0]
                        
                        else:
                         
                            g_sig = trs[-1]
                            
                        
                    Q = g_sig/(gam*(3.*gam-0.8)/(3.*gam+0.2))
                    
                    guesses = np.array([lim0cr,Q]) #,gam])
                    l0s[px,ii,0] = lim0cr
                    
                    if lim0cr > 0 and Q > 0 and gam > 0 and len(np.where(cNd00 > 0)[0]) > 3:
                    
                        kd,cov,infodict,mesg,ier = optimize.leastsq(N_fitres,guesses,args=(trs,cNd00,px,aaxt,gam),full_output=True)
                        s_sq = (infodict['fvec']**2).sum()/ (len(infodict['fvec'])-len(kd))     
                        cov = None
                        if(cov == None):
                            kderr = 0.*kd
                            pcov = np.zeros((kd.size,kd.size),np.float)
                        else:
                            pcov = cov * s_sq    
                            kderr = np.sqrt(pcov).diagonal() # 1sig confidence interval
                        
                        (lim0cr,Q) = kd#gamma) = kd
                        (lim0crerr,Qerr) = kderr#,gammaerr) = kderr
                        gamma = gam
                        gammaerr = 0.
                        
                        if lim0cr > 0:
                            (N,Nxt1,Nxt2) = Ns(trsf,lim0cr,Q,gamma,aaxt)
                            (Nt,Nxt1t,Nxt2t) = Ns(trs,lim0cr,Q,gamma,aaxt)
                                    
                            Nfit = N
                            Nfitt = Nt
 
                            
                            if px > 0:
                                Nfit += Nxt1
                                Nfitt += Nxt1t
                            
                            if px > 1:
                                Nfit += Nxt2
                                Nfitt += Nxt2t
                            
                            if px < 14:
                                Nfit += Nxt2
                                Nfitt += Nxt2t
                                
                            if px < 15:
                                Nfit += Nxt1
                                Nfitt += Nxt1t
                                                                                        
                            axes[pxx][pxy].loglog(trsf,Nfit,'-',color=cols[ii]) 
                            
                            Nfitc = Nfitt[cNd00 > 0]
                            Ndatac = cNd00[cNd00 > 0]

                            
                            
                            R2[px,ii] = 100.*np.sum(((Nfitc-Ndatac)/Ndatac)**2) / len(np.where(cNd00 > 0)[0])

                            Gain = Q*gamma*(3.*gamma-0.8)/(3.*gamma+0.2)
                            Gerr = np.abs(Qerr/Q)*Gain
    
                            if Qerr > Q:
                                Qerr = Q*0.
                            
                            if gammaerr > gamma:
                                gammaerr = gamma*0.
                                
                            if lim0crerr > lim0cr:
                                lim0crerr = lim0cr*0.
    
                            if Gerr > Gain:
                                Gerr = Gain*0.
     
                            l0s[px,ii,1] = lim0cr
                            R2s = R2[px,:]

                            lratio = l0s[px,ii,1]/l0s[px,ii,0]                            
                            
                            if len(R2s[R2s > 0]) > 0:
                                if ii == 0 or (ii == 1 and R2[px,ii]<minr2s[px] and np.abs(lratio-1.) < 0.3 ) or (ii == 2 and R2[px,ii]<minr2s[px] and np.abs(lratio-1.) < 0.3 and len(np.where(cNd00 > 0)[0]) > 8):
                                    minr2s[px] = R2[px,ii]
                                    oppt_param[cii,head,px,0,0] = Q
                                    oppt_param[cii,head,px,1,0] = gamma
                                    oppt_param[cii,head,px,2,0] = lim0cr
                                    oppt_param[cii,head,px,3,0] = Gain
        
                                    oppt_param[cii,head,px,0,1] = Qerr
                                    oppt_param[cii,head,px,1,1] = gammaerr
                                    oppt_param[cii,head,px,2,1] = lim0crerr
                                    oppt_param[cii,head,px,3,1] = Gerr                        
                                    oppt_param[cii,head,px,4,0] = ms
                                                    
                         
                if px == 0:
                    axes[pxx][pxy].set_title('%s @ MCPB=%iV %i' % (prefix,ms,eii))
                                                                           
                axes[pxx][pxy].set_xlim([3e4,3e7])
                axes[pxx][pxy].xaxis.set_ticks([1e5,1e6,1e7])        
       
                if pxx < 3:                  
                    axes[pxx][pxy].xaxis.set_ticklabels([])
                else:
                    axes[pxx][pxy].set_xlabel('e-')
            
                if pxy > 0:        
                    axes[pxx][pxy].yaxis.set_ticklabels([])
                else:
                    axes[pxx][pxy].set_ylabel('Counts')    
                    
        for px in range(16):
            pxx = np.floor((px)/4)
            pxy = px % 4

            axes[pxx][pxy].yaxis.set_ticks([1e-1*l10Max,1e-3*l10Max,1e-2*l10Max,1e-1*l10Max,l10Max])
            axes[pxx][pxy].set_ylim([Nmin,Nmax])
            axes[pxx][pxy].text(4e4,Nmax/2, '%.1f' % R2[px,0],color=cols[0])  
            axes[pxx][pxy].text(4e5,Nmax/2, '%.1f' % R2[px,1],color=cols[1])
            axes[pxx][pxy].text(4e6,Nmax/2, '%.1f' % R2[px,2],color=cols[2])
            axes[pxx][pxy].text(4e4,Nmin*1.1, '%.1f' % (l0s[px,0,1]/l0s[px,0,0]),color=cols[0])  
            axes[pxx][pxy].text(4e5,Nmin*1.1, '%.1f' % (l0s[px,1,1]/l0s[px,1,0]),color=cols[1])
            axes[pxx][pxy].text(4e6,Nmin*1.1, '%.1f' % (l0s[px,2,1]/l0s[px,2,0]),color=cols[2])
            axes[pxx][pxy].text(4e6,Nmin*3., '%.1f' % (minr2s[px]),color='k')
            
        fig.savefig(ofname,dpi=fig.dpi)

    return oppt_param 

def plot_phdparam(oppt_param,unit,serial):

    
    for head in range(2):
        
        prefix = 'FM%sH%i' % (str(unit).zfill(3),head)

        if unit > 200:
            fname = '/Users/gbclark1/Desktop/specfit/%s_comb.csv' % prefix
        else:
            fname = '/Users/gbclark1/Desktop/specfit/%s_specfit.csv' % prefix
            
        print 'Comparing %s to %s' % (prefix,fname)
                    
        datareader = open(fname,'r')
        specFit = datareader.readlines()    # Load most recent entry in file  
        
        (opAC,opDC) = np.array(map(float,specFit[4].split()))   
        
        defls = np.array(map(float,specFit[7].split()))
        gfacs = np.array(map(float,specFit[9].split()))  
        
        (mcphv,mcpb,tr,trv,tau) = np.array(map(float,specFit[12].split()))
        
        (sigma0,sigma1,kappa0,kappa1,ch_sig,ch_sp) = np.array(map(float,specFit[15].split()))
        
        (alpha,t_max,e_max) = np.array(map(float,specFit[18].split()))        
        p_ks = np.array(map(float,specFit[20].split()))      
    
        p_Qs = np.array(map(float,specFit[23].split()))  
        p_gammas = np.array(map(float,specFit[25].split()))  
        dQdV = np.array(map(float,specFit[27].split())) 
        
        lambdas = np.array(map(float,specFit[30].split())) 
            
        pxs = np.arange(0,16)   
        gammas = np.polyval(p_gammas,pxs)
        Qs = np.polyval(p_Qs,pxs)
        G = Qs*gammas*(3.*gammas-0.8)/(3.*gammas+0.2)
                                
        specfit_param = 0*oppt_param

        mcpY = np.array([])

        for mii in range(oppt_param.shape[0]):                                                    
            mcphv = oppt_param[mii,head,0,4,0] 
            
            specfit_param[mii,head,:,0,0] = Qs*10**(dQdV*(mcphv-mcpb))
            specfit_param[mii,head,:,1,0] = gammas
            
            if unit != serial[mii]:
                mcphv = 0
                
            mcpY = np.append(mcpY,mcphv)
        print mcpY              
        matplotlib.rcParams.update({'font.size': 14})
        plt.close("all")
        fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(22,16))        
        ofname = '/Users/gbclark1/Desktop/oppt_img/opptsummary_%s' % (prefix)
        plt.suptitle('Operating Point Calibration Summary for %s' % (prefix))
        cols = cm.rainbow(np.linspace(0, 1, 5)) 
 
        mcii = 0
        for mii in range(mcpY.size):
            if mcpY[mii] != 0:
                oppt_Q = oppt_param[mii,head,:,0,0]
                oppt_gam = oppt_param[mii,head,:,1,0]
                oppt_Gain = oppt_param[mii,head,:,3,0]
    
                oppt_Qerr = oppt_param[mii,head,:,0,1]
                oppt_gamerr = oppt_param[mii,head,:,1,1]
                oppt_Gainerr = oppt_param[mii,head,:,3,1]
    
                
                oppt_etaG = special.gammainc(oppt_gam,np.tile(tr,16)/oppt_Q)  
                oppt_xt = cdf(tr/lambdas[0],oppt_Q,oppt_gam)/cdf(tr,oppt_Q,oppt_gam)
                           
                if np.sum(oppt_Q) > 0:
                    axes[0][0].errorbar(pxs,oppt_Q,yerr=1.96*oppt_Qerr,fmt='o',color=cols[mcii])
                    axes[0][1].errorbar(pxs,oppt_gam,yerr=1.96*oppt_gamerr,fmt='o',color=cols[mcii])        
                    #axes[0][2].errorbar(pxs,oppt_sf,yerr=1.96*oppt_sferr,fmt='o',color=cols[mii])
                    axes[1][0].errorbar(pxs,oppt_Gain,yerr=1.96*oppt_Gainerr,fmt='o',color=cols[mcii])   
                    axes[1][1].semilogy(pxs,oppt_etaG,'o',color=cols[mcii])
                    if np.sum(oppt_xt) > 0:
                        axes[1][2].semilogy(pxs,oppt_xt,'o',color=cols[mcii])
                        
                specfit_Q = specfit_param[mii,head,:,0,0]
                specfit_gam = specfit_param[mii,head,:,1,0]
                specfit_Gain = specfit_Q*specfit_gam*(3.*specfit_gam-0.8)/(3.*specfit_gam+0.2)
                specfit_etaG = special.gammainc(specfit_gam,np.tile(tr,16)/specfit_Q)  
                specfit_xt = cdf(tr/lambdas[0],specfit_Q,specfit_gam)/cdf(tr,specfit_Q,specfit_gam)
    
    
                    
                axes[0][0].semilogy(pxs,specfit_Q,'-',color=cols[mcii],label='MCPB=%iV' % (mcpY[mii]))
                axes[0][1].plot(pxs,specfit_gam,'-',color=cols[mcii],label='MCPB=%iV' % (mcpY[mii]))        
                axes[1][0].semilogy(pxs,specfit_Gain,'-',color=cols[mcii],label='MCPB=%iV' % (mcpY[mii])) 
                axes[1][1].semilogy(pxs,specfit_etaG,'-',color=cols[mcii],label='MCPB=%iV' % (mcpY[mii]))
    
                if np.sum(specfit_xt) > 0:
                    axes[1][2].semilogy(pxs,specfit_xt,'-',color=cols[mcii],label='MCPB=%iV' % (mcpY[mii]))
                mcii+=1
            
            
                    
        axes[0][0].set_xlim([-0.5,15.5])        
        axes[0][0].set_ylim([1e5,1e8])
        axes[0][0].set_ylabel('Q (e-)')
        axes[0][0].set_xlabel('PX')
        axes[0][0].set_yscale('log')
    
        axes[0][1].set_xlim([-0.5,15.5])      
        axes[0][1].set_ylim([0.1,100])
        axes[0][1].set_ylabel('$\gamma$',fontsize=18)
        axes[0][1].set_xlabel('PX')
        axes[0][1].legend(loc='upper left')
        axes[0][1].set_yscale('log')


        fig.delaxes(axes[0][2])

    
        axes[1][0].set_xlim([-0.5,15.5])                      
        axes[1][0].set_ylim([1e4,1e9])
        axes[1][0].plot([-0.5,15.5],[2e6,2e6],'k--',linewidth=2)
        axes[1][0].set_ylabel('Gain (e-)')
        axes[1][0].set_xlabel('PX')
        axes[1][0].set_yscale('log')
    
        axes[1][1].set_xlim([-0.5,15.5])  
        axes[1][1].set_ylim([1e-4,1e0])
        axes[1][1].plot([-0.5,15.5],[0.05,0.05],'k--',linewidth=2)
        axes[1][1].set_ylabel('Signal Loss')
        axes[1][1].set_xlabel('PX')
        axes[1][1].set_yscale('log')
                  
        axes[1][2].set_xlim([-0.5,15.5])  
        axes[1][2].set_ylim([1e-4,1e0])
        axes[1][2].plot([-0.5,15.5],[0.05,0.05],'k--',linewidth=2)
        axes[1][2].set_ylabel('N(tr/AAXT,Q,$\gamma$)/N(tr,Q,$\gamma$)')
        axes[1][2].set_xlabel('PX')
        axes[1][2].set_yscale('log')
                                    
        plt.tight_layout()
       
        fig.savefig(ofname,dpi=fig.dpi)  





def plot_tswpnorm(cs,treshs,aaxts,unit,serial,tag,oppt_param):
    

    for head in range(2):
        plt.close("all")
        fig,axes = plt.subplots(nrows=4,ncols=4,figsize=(8,8))
        plt.subplots_adjust(hspace=0.,wspace=0.)   
        prefix = 'FM%sH%i' % (str(unit).zfill(3),head)   

        ofname = '/Users/gbclark1/Desktop/oppt_img/totalnorm_oppttswp_%s_%s' % (prefix,tag)


        if unit > 200:
            fname = '/Users/gbclark1/Desktop/specfit/%s_comb.csv' % prefix
        else:
            fname = '/Users/gbclark1/Desktop/specfit/%s_specfit.csv' % prefix
            
        print 'Comparing %s to %s' % (prefix,fname)
                    
        datareader = open(fname,'r')
        specFit = datareader.readlines()    # Load most recent entry in file  
        (mcphv,mcpb,tr,trv,tau) = np.array(map(float,specFit[12].split()))


        mcii = 0
        for cii in range(len(cs)):

            if serial[cii] == unit:
               

                  counts_sub = np.sum(np.sum(cs[cii][:,:,head,:,:]*1.,axis=0),axis=0)
                  cols = cm.rainbow(np.linspace(0, 1, 5)) 


                  trs = threshs[cii][head]
                  aaxt = aaxts[cii]
                  
                  #if unit < 200:
                  #    aaxt[1] = 8e-4
        
                  #counts_sub[:,-2:] = 0.
                  #counts_sub[:,0] = 0.
                  counts_sub[counts_sub == 1] = 0.    
    
                  Nmax = 20.
                  Nmin = 5e-3
                
                  l10Max = 10**np.floor(np.log10(Nmax))                    
                
                  trsf = np.logspace(3,8,100)
    
                  for px in range(16):
                
                      pxx = np.floor((px)/4)
                      pxy = px % 4
        
                      Nd00 = counts_sub[:,px]

                      Q = oppt_param[cii,head,px,0,0] 
                      gamma = oppt_param[cii,head,px,1,0] 
                      lim0cr = oppt_param[cii,head,px,2,0] 
                      mcphv = oppt_param[cii,head,px,4,0]             
                    
                      if(np.sum(Nd00) > 0 and lim0cr > 0):
     
                          
                          axes[pxx][pxy].loglog(trs,Nd00/np.tile(lim0cr,Nd00.size),'.',color=cols[mcii],label='MCPB=%iV' % mcphv)
     
                    
                          if lim0cr > 0:
                              (N,Nxt1,Nxt2) = Ns(trsf,lim0cr,Q,gamma,aaxt)
                           
                              Nfit = N
     
                              if px > 0:
                                  Nfit += Nxt1
                                 
                                
                              if px > 1:
                                  Nfit += Nxt2
                                  
                                
                              if px < 14:
                                  Nfit += Nxt2
                                
                                    
                              if px < 15:
                                  Nfit += Nxt1
                                
                                                                                            
                              axes[pxx][pxy].loglog(trsf,Nfit/np.tile(lim0cr,Nfit.size),'-',color=cols[mcii]) 
    
     
      
        
                      if px == 0:
                          axes[pxx][pxy].set_title('%s' % (prefix))                                                    
                      axes[pxx][pxy].set_xlim([3e4,3e7])
                      axes[pxx][pxy].xaxis.set_ticks([1e5,1e6,1e7])        
    
                      if px == 3:
                          axes[pxx][pxy].legend(loc='upper right',bbox_to_anchor=(1.,1.5),fontsize=10,ncol=3)
               
                      if pxx < 3:                  
                          axes[pxx][pxy].xaxis.set_ticklabels([])
                      else:
                          axes[pxx][pxy].set_xlabel('e-')
                
                      if pxy > 0:        
                          axes[pxx][pxy].yaxis.set_ticklabels([])
                      else:
                          axes[pxx][pxy].set_ylabel('CDF') 
                          
                      axes[pxx][pxy].yaxis.set_ticks([1e-1*l10Max,1e-3*l10Max,1e-2*l10Max,1e-1*l10Max,l10Max])
                      axes[pxx][pxy].set_ylim([Nmin,Nmax])
                      
                      if mcii == 0:
                          axes[pxx][pxy].text(4e4,Nmin*2., 'PX%i' % px)
                          axes[pxx][pxy].plot([tr,tr],[0.001,100.],'k--')
                    


                  
                  mcii+=1
        fig.savefig(ofname,dpi=fig.dpi) 

def plot_tswpnormPDF(cs,treshs,aaxts,unit,serial,tag,oppt_param):
    
 
    mcii = 0
    for cii in range(len(cs)):
        if serial[cii] == unit:    
   
            for head in range(2):
                plt.close("all")
                fig,axes = plt.subplots(nrows=4,ncols=4,figsize=(8,8))
                plt.subplots_adjust(hspace=0.,wspace=0.)   
                prefix = 'FM%sH%i' % (str(unit).zfill(3),head)   
                
                ofname = '/Users/gbclark1/Desktop/oppt_img/totalnorm_opptPDFtswp_%s_%s_%i' % (prefix,tag,mcii)


                if unit > 200:
                    fname = '/Users/gbclark1/Desktop/specfit/%s_comb.csv' % prefix
                else:
                    fname = '/Users/gbclark1/Desktop/specfit/%s_specfit.csv' % prefix
                    
                print 'Comparing %s to %s' % (prefix,fname)
                            
                datareader = open(fname,'r')
                specFit = datareader.readlines()    # Load most recent entry in file  
                (mcphv,mcpb,tr,trv,tau) = np.array(map(float,specFit[12].split()))
                  
                counts_sub = np.sum(np.sum(cs[cii][:,:,head,:,:]*1.,axis=0),axis=0)
                cols = cm.rainbow(np.linspace(0, 1, 5)) 
    
                trs = threshs[cii][head]
                aaxt = aaxts[cii]
              
                #if unit < 200:
                #    aaxt[1] = 8e-4
    
                #counts_sub[:,-2:] = 0.
                #counts_sub[:,0] = 0.
                counts_sub[counts_sub == 1] = 0.    

                Nmax = 5.
                Nmin = 5e-3
            
                l10Max = 10**np.floor(np.log10(Nmax))                    
            
                trsf = np.logspace(3,8,100)

                for px in range(16):
            
                    pxx = np.floor((px)/4)
                    pxy = px % 4
    
                    Nd00 = counts_sub[:,px]

                    Q = oppt_param[cii,head,px,0,0] 
                    gamma = oppt_param[cii,head,px,1,0] 
                    lim0cr = oppt_param[cii,head,px,2,0] 
                    mcphv = oppt_param[cii,head,px,4,0]             
           
                    if lim0cr > 0:
                        (N,Nxt1,Nxt2) = Nspdf(trsf,1.,Q,gamma,aaxt)

                       
                        Nxt = np.zeros(trsf.size,np.float)
 
                        if px > 0:
                            Nxt += Nxt1
                         
                        
                        if px > 1:
                            Nxt += Nxt2
                          
                        
                        if px < 14:
                            Nxt += Nxt2
                        
                            
                        if px < 15:
                            Nxt += Nxt1
                        
                        axes[pxx][pxy].loglog(trsf,N/np.max(N),'b-') 
                        axes[pxx][pxy].loglog(trsf,Nxt/np.max(N),'r-') 
    
                    if px == 0:
                        axes[pxx][pxy].set_title('%s @ MCPB=%iV' % (prefix,mcphv))                                                    
                    axes[pxx][pxy].set_xlim([3e4,3e7])
                    axes[pxx][pxy].xaxis.set_ticks([1e5,1e6,1e7])        

    
                  
                    if pxx < 3:                  
                        axes[pxx][pxy].xaxis.set_ticklabels([])
                    else:
                        axes[pxx][pxy].set_xlabel('e-')
            
                    if pxy > 0:        
                        axes[pxx][pxy].yaxis.set_ticklabels([])
                    else:
                        axes[pxx][pxy].set_ylabel('PDF') 
                      
                    axes[pxx][pxy].yaxis.set_ticks([1e-1*l10Max,1e-3*l10Max,1e-2*l10Max,1e-1*l10Max,l10Max])
                    axes[pxx][pxy].set_ylim([Nmin,Nmax])
                    
                    axes[pxx][pxy].text(5e6,2., 'PX%i' % px)
                    if px == 0:
                        axes[pxx][pxy].text(6e4,8e-3,'XT',color='red')
                        axes[pxx][pxy].text(5e6,8e-3,'SIG',color='blue')
                    axes[pxx][pxy].plot([tr,tr],[0.001,100.],'k--')
                
                        

                fig.savefig(ofname,dpi=fig.dpi)           
            mcii+=1
       
            
def get_serial_number(filename):

    st1 = filename.split('/')
    st2 = st1[-1].split('_')
    obs = st2[0].split('mms')[-1]
    spec = st2[1].split('s')[-1]

    species = st2[1].split('s')[0]

    serial_des = np.array([[209,211,212,207],[203,204,201,202],[206,216,213,214],[215,210,208,205]])
    serial_dis = np.array([[2,5,6,14],[14,1,9,12],[8,15,16,7],[10,11,3,4]])

    if species == 'de':
        serial = serial_des[int(obs)-1][int(spec)]
    else:
        serial = serial_dis[int(obs)-1][int(spec)]


    return serial

def threshconver(ser,h,xvalues):
    serial = str(ser)
    head = str(h)
    serial_st = 'DES S/N '+serial #+' '+'Head '+head
    head_st = 'Head '+head
    path = '/Users/gbclark1/svn/FPIdataProcTools/FPIthresholdConvert/data/'
    xlsfile = path+'FPIthresholdConversions.xls'

    book = open_workbook(xlsfile)
    name = book.sheet_names()
    sheet = book.sheet_by_name('Sheet1')

    # find the first column where the serial number matches
    for cell in range(sheet.ncols):
        dummy = sheet.col(cell)
        if dummy[0].value  == serial_st: break

    # determine the column for the correct head
    if sheet.col(cell)[1].value == head_st:
        cell = cell
    else:
        cell = cell+1


    #for electrons

    if int(serial) > 200:
        thresh_volt = []
        data = []
        for i in range(2,26,1):
            data.append(sheet.col(cell)[i].value)
            thresh_volt.append(sheet.col(0)[i].value)

    #for ions
    if int(serial) < 200:
        thresh_volt =[]
        data = []
        for i in range(29,47,1):
            data.append(sheet.col(1)[i].value)
            thresh_volt.append(sheet.col(0)[i].value)

    threshs = np.interp(xvalues,thresh_volt,data)

    return threshs


def create_counts_matrix(dir):


    cdf_filenames = glob.glob(dir+'/*sigthrsh*.cdf')
    serial=[]
    loop_index = 0
    counts = np.zeros(shape=(16,8,8,2,16,16))
    threshs = np.zeros(shape=(16,2,16))
    mcps = np.zeros(shape=(16,2))


    #pdb.set_trace()
    ####  ADD THRESH 01 AND MAKE ARRAYS LENTGH 16 [4DES,4DIS, 2MCPV = 16]

    for files in cdf_filenames:
        print files
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
                counts_h1_new_shape = new_counts_h1.reshape(8,8,16) # or 16,8,8 ....it's a mystery!
                # pdb.set_trace()
                for jj in range(len(it_uniq)):
                    for kk in range(len(thresh_uniq)):
                        for ii in range(len(eg_uniq)):

                            if ll == 0:
                                counts[loop_index,jj,ii,ll,kk,pp]=counts_h0_new_shape[jj,ii,kk] #[kk,jj,ii]
                            else:
                                counts[loop_index,jj,ii,ll,kk,pp]=counts_h1_new_shape[jj,ii,kk] #[kk,jj,ii]

        loop_index=loop_index+1

    aaxts = np.array([[0.022, 0.007],[0.022, 0.007],[0.022, 0.007],[0.022, 0.007],[0.022, 0.007],[0.022, 0.007],[0.022, 0.007],[0.022, 0.007],[2.49999994e-03,9.99999997e-07],[2.49999994e-03,9.99999997e-07],[2.49999994e-03,9.99999997e-07],[2.49999994e-03,9.99999997e-07],[2.49999994e-03,9.99999997e-07],[2.49999994e-03,9.99999997e-07],[2.49999994e-03,9.99999997e-07],[2.49999994e-03,9.99999997e-07]])
    return counts, serial, threshs, mcps, aaxts


Obs = 2
quarter=0
for FPImode in [0,1]:
    print '%i %i %i' % (Obs,quarter,FPImode)


    dir = '/Users/gbclark1/Desktop/cdf_files/'
    data = create_counts_matrix(dir)
    counts = data[0]
    serial = data[1]
    threshs = data[2]
    mcps = data[3]
    aaxts = data[4]

    # fname = 'IDLDataObs%i_FPI008_corrected' % (Obs)
    # oppt_data = io.readsav('%s.sav' % fname)
    # counts = oppt_data['data']['counts']
    # threshs = oppt_data['data']['thresh']
    # energy = oppt_data['data']['energy']
    # mcps = oppt_data['data']['mcpvlt']
    # aaxts = oppt_data['data']['anodxt']
    # serial = oppt_data['data']['serial']
    #
    
    plt.close("all")
    fig,axes = plt.subplots(nrows=4,ncols=4,figsize=(8,8))
    plt.subplots_adjust(hspace=0.,wspace=0.) 
    
    oppt_param = np.zeros((len(counts),2,16,5,2),np.float)
    
    for head in [0,1]: 
        for cii in range(len(counts)):
            for eii in range(1):           
                
                 if (FPImode == 0 and serial[cii] > 200) or (FPImode == 1 and serial[cii] < 200):

                     coppt_param = plot_tswpsum(fig,axes,counts,eii,threshs[cii][head],mcps[cii][head],aaxts[cii],serial[cii],head,cii,'all_%i' % (eii))
                     oppt_param[cii,head,:,:,:] = coppt_param[cii,head,:,:,:]
                 
                 
    for unit in np.unique(serial):
         
         if (FPImode == 0 and unit >200) or (FPImode == 1 and unit < 200):
             plot_phdparam(oppt_param,unit,serial)
             plot_tswpnorm(counts,threshs,aaxts,unit,serial,'',oppt_param)
         #plot_tswpnormPDF(counts,threshs,aaxts,unit,serial,'',oppt_param)
                 
                     