# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 15:54:50 2015

@author: fpi_egse
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as matplotlib
import numpy as np
import cPickle
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


def plot_tswpsum(fig,axes,cs,eii,trs,ms,aaxt,unit,head,cii,output_dir,oppt_dir):


    plt.figure(fig.number)
    
    oppt_param = np.zeros((len(cs),2,16,5,2),np.float)

    prefix = 'FM%sH%i' % (str(unit).zfill(3),head)   
    
    ofname = '%s/tswpfit_%s_%s_%i' % (output_dir,prefix,oppt_dir,ms)
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

            Nd00 = counts_sub[px,:] 

            axes[pxx][pxy].cla()
            
            if(len(np.where(Nd00 > 0)[0]) > 3):
                p1 = Nd00[0:5]
                p2 = Nd00[5:10]
                p3 = Nd00[10:15]
                
                mms = [np.mean(p1[p1 > 0]),np.mean(p2[p2 > 0]),np.mean(p3[p3 > 0])]                             
                
                axes[pxx][pxy].loglog(trs,Nd00,'k.')
                for ii in range(3):
                    print ii
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
                       
                        kd,cov,infodict,mesg,ier = optimize.leastsq(of.N_fitres,guesses,args=(trs[cNd00>0],cNd00[cNd00>0],px,aaxt,gam),full_output=True)

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
                            (N,Nxt1,Nxt2) = of.Ns(trsf,lim0cr,Q,gamma,aaxt)
                            (Nt,Nxt1t,Nxt2t) = of.Ns(trs,lim0cr,Q,gamma,aaxt)
                                    
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

def plot_phdparam(fig,axes,oppt_param,unit,serial,oppt_dir,specfit_dir,output_dir):

    plt.figure(fig.number)
    
    for head in range(2):
        
        prefix = 'FM%sH%i' % (str(unit).zfill(3),head)
        
        if unit > 200:
            sensor = 'des'
        else:
            sensor = 'dis'
            
        fname = '%s/specfits/%s/%s%sH%i_%s.csv' % (of.base_dir,specfit_dir,sensor.upper(),str(unit).zfill(3),head,specfit_dir)
        ofname = '%s/tswpparam_%s_%s' % (output_dir,prefix,oppt_dir)           
       
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
           
        p_gammas = np.array(map(float,specFit[26].split()))  
        dQdV = np.array(map(float,specFit[28].split())) 
        lambdas = np.array(map(float,specFit[31].split()))  # current specfits aren't right       
        
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

          
        matplotlib.rcParams.update({'font.size': 14})

        for xii in range(2):
            for yii in range(3):
                axes[xii][yii].cla()
                
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
                oppt_xt = of.cdf(tr/lambdas[0],oppt_Q,oppt_gam)/of.cdf(tr,oppt_Q,oppt_gam)
                           
                if np.sum(oppt_Q) > 0:
                    axes[0][0].errorbar(pxs,oppt_Q,yerr=1.96*oppt_Qerr,fmt='o',color=cols[mcii])
                    axes[0][1].errorbar(pxs,oppt_gam,yerr=1.96*oppt_gamerr,fmt='o',color=cols[mcii])  
                    axes[0][1].text(2,0.2,prefix,fontsize=20)
                    
                    #axes[0][2].errorbar(pxs,oppt_sf,yerr=1.96*oppt_sferr,fmt='o',color=cols[mii])
                    axes[1][0].errorbar(pxs,oppt_Gain,yerr=1.96*oppt_Gainerr,fmt='o',color=cols[mcii])   
                    axes[1][1].semilogy(pxs,oppt_etaG,'o',color=cols[mcii])
                    if np.sum(oppt_xt) > 0:
                        axes[1][2].semilogy(pxs,oppt_xt,'o',color=cols[mcii])
                        
                specfit_Q = specfit_param[mii,head,:,0,0]
                specfit_gam = specfit_param[mii,head,:,1,0]
                specfit_Gain = specfit_Q*specfit_gam*(3.*specfit_gam-0.8)/(3.*specfit_gam+0.2)
                specfit_etaG = special.gammainc(specfit_gam,np.tile(tr,16)/specfit_Q)  
                specfit_xt = of.cdf(tr/lambdas[0],specfit_Q,specfit_gam)/of.cdf(tr,specfit_Q,specfit_gam)
    
    
                    
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


        #fig.delaxes(axes[0][2])

    
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




def plot_tswpnorm(fig,axes,cs,threshs,aaxts,unit,serial,oppt_param,oppt_dir,specfit_dir,output_dir):
    
    plt.figure(fig.number)
    for head in range(2):
        
       
        prefix = 'FM%sH%i' % (str(unit).zfill(3),head)   

           
        ofname = '%s/tswpfit_norm_%s_%s' % (output_dir,prefix,oppt_dir)           
        if unit > 200:
            sensor = 'des'
        else:
            sensor = 'dis'
            
        fname = '%s/specfits/%s/%s%sH%i_%s.csv' % (of.base_dir,specfit_dir,sensor.upper(),str(unit).zfill(3),head,specfit_dir)
        datareader = open(fname,'r')
        specFit = datareader.readlines()    # Load most recent entry in file  
        (mcphv,mcpb,tr,trv,tau) = np.array(map(float,specFit[12].split()))

        for pxx in range(4):
            for pxy in range(4):
                axes[pxx][pxy].cla()            


        mcii = 0
        for cii in range(len(cs)):

            if serial[cii] == unit:
               

                  counts_sub = np.sum(np.sum(cs[cii][:,:,head,:,:]*1.,axis=0),axis=0)
                  cols = cm.rainbow(np.linspace(0, 1, 5)) 


                  trs = threshs[cii][head]
                  aaxt = aaxts[cii]
                  

                  counts_sub[counts_sub == 1] = 0.    
    
                  Nmax = 20.
                  Nmin = 5e-3
                
                  l10Max = 10**np.floor(np.log10(Nmax))                    
                
                  trsf = np.logspace(3,8,100)
    
                  for px in range(16):
                
                      pxx = np.floor((px)/4)
                      pxy = px % 4
        
                      Nd00 = counts_sub[px,:]

                      Q = oppt_param[cii,head,px,0,0] 
                      gamma = oppt_param[cii,head,px,1,0] 
                      lim0cr = oppt_param[cii,head,px,2,0] 
                      mcphv = oppt_param[cii,head,px,4,0]             
                    
                              
                    
                      if(np.sum(Nd00) > 0 and lim0cr > 0):
     
                          
                          axes[pxx][pxy].loglog(trs,Nd00/np.tile(lim0cr,Nd00.size),'.',color=cols[mcii],label='MCPB=%iV' % mcphv)
     
                    
                          if lim0cr > 0:
                              (N,Nxt1,Nxt2) = of.Ns(trsf,lim0cr,Q,gamma,aaxt)
                           
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
        #plt.close(fig)
        #del fig,axes
