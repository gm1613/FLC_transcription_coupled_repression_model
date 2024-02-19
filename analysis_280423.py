#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 09:07:55 2023

@author: menon
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

#%%
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

#%% Section 1: Plotting transcriptional activity over time 

type1b_out = np.loadtxt("20230117_autopathway_flc_type1b_tid1_only_analog_fca1_0124.txt",delimiter=' ')

t = type1b_out[:,0]*0.917 # Time t is multiplied by 0.917 to convert units from cell cycles (22h) to days (24h)

sense_prox = type1b_out[:,1]
sense_dist = type1b_out[:,2]
as_prox = type1b_out[:,3]
as_dist = type1b_out[:,4]
total_transc = sense_prox+sense_dist+as_prox+as_dist

fig, ax=plt.subplots(nrows=4,ncols=1)
ax[0].vlines(t,0,sense_prox,'grey',label='Sense Proximal')
ax[0].set_ylim([0,1.8])
ax[0].set_xlim([0,21])
ax[1].vlines(t,0,sense_dist,'grey',label='Sense Distal')
ax[1].set_ylim([0,1.8])
ax[1].set_xlim([0,21])
ax[2].vlines(t,0,as_prox,'grey',label='Antisense Proximal')
ax[2].set_ylim([0,0.5])
ax[2].set_xlim([0,21])
ax[2].set_yticks([0,0.5])
ax[3].vlines(t,0,as_dist,'grey',label='Antisense Distal')
ax[3].set_ylim([0,0.5])
ax[3].set_xlim([0,21])
ax[3].set_yticks([0,0.5])
fig.subplots_adjust(hspace=0.5)

fig, ax=plt.subplots(nrows=1,ncols=1)
ax.vlines(t,0,sense_dist,'grey')
ax.set_ylim([0,4])
ax.set_xlim([0,21])
fig.subplots_adjust(hspace=0.5)


# For 10, 14, 21 day timepoints, compute and print frequency averaged over previous  6hr
print("sense_dist:",np.mean(sense_dist[234:241]),np.mean(sense_dist[330:337]),np.mean(sense_dist[498:505]))



#%% Section 1b: Plotting differences in prox/dist ratio

type1b_out = np.loadtxt("20230117_autopathway_flc_type1b_tid1_fca1_1223.txt",delimiter=' ')


t = type1b_out[:,0]*0.917 # For the 6.75 days to 7 days window, take indices 162-168
                          # 0.25 day window is consistent with ~6hr lifetime of FLC mRNA
                          # Time t is multiplied by 0.917 to convert units from cell cycles (22h) to days (24h)


sense_prox_fca1 = np.mean(type1b_out[162:168,1])
sense_dist_fca1 = np.mean(type1b_out[162:168,2])
as_prox_fca1 = np.mean(type1b_out[162:168,3])
as_dist_fca1 = np.mean(type1b_out[162:168,4])

type1b_out = np.loadtxt("20230117_autopathway_flc_type1b_tid1_fca3_1223.txt",delimiter=' ')

sense_prox_fca3 = np.mean(type1b_out[162:168,1])
sense_dist_fca3 = np.mean(type1b_out[162:168,2])
as_prox_fca3 = np.mean(type1b_out[162:168,3])
as_dist_fca3 = np.mean(type1b_out[162:168,4])

type1b_out = np.loadtxt("20230117_autopathway_flc_type1b_tid1_clfnew_1223.txt",delimiter=' ')

sense_prox_clf = np.mean(type1b_out[162:168,1])
sense_dist_clf = np.mean(type1b_out[162:168,2])
as_prox_fca3 = np.mean(type1b_out[162:168,3])
as_dist_fca3 = np.mean(type1b_out[162:168,4])



type1b_out = np.loadtxt("20230117_autopathway_flc_type1b_tid1_Ler_1223.txt",delimiter=' ')

sense_prox_wt = np.mean(type1b_out[162:168,1])
sense_dist_wt = np.mean(type1b_out[162:168,2])
as_prox_wt = np.mean(type1b_out[162:168,3])
as_dist_wt = np.mean(type1b_out[162:168,4])

r0 = sense_prox_wt/sense_dist_wt
rfca3 = (sense_prox_fca3/sense_dist_fca3)/r0
rfca1 = (sense_prox_fca1/sense_dist_fca1)/r0
rclf = (sense_prox_clf/sense_dist_clf)/r0




#%% Section 2: Plotting changes over time (averaged over trajectories)
name_tag = "fca1_1223"
type1a_out = np.loadtxt("20230117_autopathway_flc_type1a_tid1_"+name_tag+".txt",delimiter=' ')

t = type1a_out[:,0]*0.917# Time t is multiplied by 0.917 to convert units from cell cycles (22h) to days (24h)

k27_whole = type1a_out[:,2]
k4 = type1a_out[:,3]
n_off = type1a_out[:,4]

fig, ax=plt.subplots()
ax.plot(t,k27_whole,'r',linewidth=2,label='H3K27me3')
ax.plot(t,k4,'b',linewidth=2,label='H3K4me1')
plt.title(name_tag)
plt.legend()

fig, ax=plt.subplots()
ax.plot(t,n_off/5000,'k--',linewidth=2,label='$N_{OFF}$')
plt.title(name_tag)
plt.legend()


#%% 
#Print average coverage of H3K4me1 at 7, 14, and 21 day timepoints
print("k4:",k4[382],k4[763],k4[1145])

#Print average coverage of H3K27me3 at 7, 14, and 21 day timepoints
print("k27_whole:",k27_whole[382],k27_whole[763],k27_whole[1145])

#%% Section 3: Plotting switchoff time histograms: WT vs fca3 vs fca1
type1c_out_wt = np.loadtxt("20230117_autopathway_flc_type1c_tid1_ler_1223.txt",delimiter=' ')
type1c_out_fca3 = np.loadtxt("20230117_autopathway_flc_type1c_tid1_fca3_1223.txt",delimiter=' ')
type1c_out_fca1 = np.loadtxt("20230117_autopathway_flc_type1c_tid1_fca1_1223.txt",delimiter=' ')

nswitched_wt = np.sum(type1c_out_wt>0)
nswitched_fca3 = np.sum(type1c_out_fca3>0)
nswitched_fca1 = np.sum(type1c_out_fca1>0)

wt_switchtimes = np.zeros(nswitched_wt)
fca3_switchtimes = np.zeros(nswitched_fca3)
fca1_switchtimes = np.zeros(nswitched_fca1)

n = type1c_out_wt.shape[0]

counter = 0;
for i in range(0,n):
    if (type1c_out_wt[i]>0.01):
        wt_switchtimes[counter]=type1c_out_wt[i]
        counter+=1

n = type1c_out_fca3.shape[0]

counter=0;
for i in range(0,n):
    if (type1c_out_fca3[i]>0.01):
        fca3_switchtimes[counter]=type1c_out_fca3[i]
        counter+=1
        
counter=0;
for i in range(0,n):
    if (type1c_out_fca1[i]>0.01):
        fca1_switchtimes[counter]=type1c_out_fca1[i]
        counter+=1
          

plt.figure()
sns.histplot(fca1_switchtimes*0.917,stat='density',bins=75,color='blue')
sns.histplot(fca3_switchtimes*0.917,stat='density',bins=50,color='green')
sns.histplot(wt_switchtimes*0.917,stat='density',bins=30,color='gray')

  


#%%  Section 4: Plotting steady state spatial profile of H3K27me3 at FLC

pf = np.loadtxt("20230117_autopathway_flc_type1d_tid1_Ler_short_0124.txt",delimiter=' ')

m=60 #Number of nucleosomes
ik27 = np.arange(0,m,1)
ik4 = np.arange(0,m,1)

pf0k27 = pf[ik27,1]
pf0k4 = pf[ik4,2]


N=30; #Number of histones

pfk27 = np.zeros((N,1))
pfk4 = np.zeros((N,1))


j=0
for i in range(0,m,2):
    pfk27[j]=(pf0k27[i]+pf0k27[i+1])/2
    pfk4[j]=(pf0k4[i]+pf0k4[i+1])/2
    j=j+1

window = 200*np.arange(0,30,1) # Assuming a total nucleosomal DNA length (including linker) of 200bp


avg_k27 = np.amax(pfk27,axis=0) 
avg_k4 = np.amax(pfk4,axis=0)



#%% Performing convolution with ChIP fragment size distribution (see Wu et al. 2016)

######################################################
# EXTRACTING ChIP FRAGMENT SIZE DISTRIBUTION FROM Wu et al. 2016
llim = 0          
ulim = 0.25      
lcd = 848
ucd = 196

pcd = np.array((671,318,574,406,530,623,637,686,690,730))

# pcd = np.array((671,318,350,406,500,545,580,600,630,680))

fdist = (ulim-llim)*(lcd-pcd)/(lcd-ucd) #PROBABILITY
x = np.arange(100,1001,100)             #FRAGMENT SIZE IN bp


########################################################

# SPECIFYING ChIP FRAGMENT SIZE DISTRIBUTION (FROM Wu et al. 2015)
#PROBABILITY
fdist= np.array((0.0678,0.2010,0.1050,0.1694,0.1219,0.0862,0.0809,0.0621,0.0605,0.0452))
x = np.arange(100,1001,100)#FRAGMENT SIZE IN bp


xc = np.arange(1,61,1)
coverage = np.zeros((100,2))
nuc_x = np.arange(1,100,2)
m_density = np.zeros((100,2))
p_density = np.zeros((100,2))


mwindow = 200*np.arange(0,50,1)
dwindow = 100*np.arange(0,100,1)
##########################################################
def density(x,coverage,fdist):
    ymin=1
    ymax=10
    f=0
    for y in range(ymin,ymax+1):
        h=0
        for z in range(x-y+1,x+1):
            for w in range(z,z+y-1+1):
                if (w>0 and w<=100):
                    h=h+coverage[w-1]
        f=f+fdist[y-1]*h;
    return f
###########################################################
for i in range(1,61):
    if np.mod(i,2)==1:
        coverage[i-1,1-1]=pfk27[int(np.floor(i/2))+1-1,1-1]
        coverage[i-1,2-1]=pfk4[int(np.floor(i/2))+1-1,1-1]
    else:
        coverage[i-1,1-1]=pfk27[int(i/2)-1,1-1]
        coverage[i-1,2-1]=pfk4[int(i/2)-1,1-1]      
        
for i in range(1,61):
    m_density[i-1,1-1] = density(xc[i-1],coverage[:,1-1],fdist)
    m_density[i-1,2-1] = density(xc[i-1],coverage[:,2-1],fdist)   
    



mk27 = np.amax(m_density[:,0])
fk27 = avg_k27/mk27



nm_density = fk27*m_density[:,0]

#np.savetxt('Ler_profile_0124.txt',nm_density,delimiter='\n')

#%% Plotting predicted ChIP profiles
nm1 = np.loadtxt('clfnew_profile_0124.txt')
nm2 = np.loadtxt('Ler_profile_0124.txt')
ratio = 1/5
fig, ax=plt.subplots()
plt.fill_between(dwindow/1000,0,nm2,facecolor='royalblue',label='WT')
plt.fill_between(dwindow/1000,0,nm1[0:100],facecolor='r',label='clf')
plt.ylim([0,1.6])
plt.xlim([-2,8])
ax.set_yticks([0,0.5,1,1.5])
ax.set_aspect(abs((8--2)/(1-0))*ratio)
plt.title('H3K27me3')
plt.legend()
plt.savefig('clf_H3K27me3_profile_new.pdf',format='pdf')

