import sys, glob

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

sys.path.append('../')
from growthrate import GrowthRate

file = 'bunchmonitor_qpv-0.2_nonlinearq26.h5'
file = 'bunchmonitor_qpv-1.25_nonlinearq26.h5'
file = 'bunchmonitor_qpv-0.65_nonlinearq26.h5'

qpv = file.split('qpv')[1].split('_')[0]
#1D-analysis: npersegs
fig, ax = plt.subplots(dpi=150)
npersegs = [32, 64, 128, 256, 512, 1024]
noverlap = 0.7
grvals = np.array([])
colors = ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e']
for i, nperseg in enumerate(npersegs):
    #analysis
    gr = GrowthRate(file=file, nperseg=nperseg, noverlap=nperseg*noverlap, grstart=1000)
    gr.MWFFT_analysis()
    #gr.plot()
    grvals = np.append(grvals, gr.slope*1e3)
    ax.plot(gr.points, gr.log_peak_amplitude_all, marker='.', lw=2, c=colors[i])
    ax.plot(gr.grx, gr.gry, lw=2, c=colors[i], label=f'size: {nperseg}')
    ax.text(0.7, 0.4-0.05*i, f'gr: {round(gr.slope*1e3,2)}'+' [$10^3/{n}_{turns}$]', transform=ax.transAxes, color=colors[i])

ax.legend()
ax.set_xlabel('MWFFT points')
ax.set_ylabel('Log(MWFFT Max. Amplitude)')
ax.set_title(f'Dependency with window size | $Q\'_v: {qpv}$ | overlap: 70%', fontweight='bold')
fig.savefig('MWFFTnperseg3.png')
fig.tight_layout()
plt.show()


#1D-analysis: noverlap
fig, ax = plt.subplots(dpi=150)
nperseg = 64
noverlaps = [0.9, 0.7, 0.5, 0.3]
grvals = np.array([])
colors = ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e']
for i, noverlap in enumerate(noverlaps):
    #analysis
    gr = GrowthRate(file=file, nperseg=nperseg, noverlap=nperseg*noverlap, grstart=1500)
    gr.MWFFT_analysis()
    #gr.plot()
    grvals = np.append(grvals, gr.slope*1e3)
    ax.plot(gr.points, gr.log_peak_amplitude_all, marker='.', lw=2, c=colors[i])
    ax.plot(gr.grx, gr.gry, lw=2, c=colors[i], label=f'overlap: {noverlap}')
    ax.text(0.7, 0.4-0.05*i, f'gr: {round(gr.slope*1e3,2)}'+' [$10^3/{n}_{turns}$]', transform=ax.transAxes, color=colors[i])

ax.legend()
ax.set_xlabel('MWFFT points')
ax.set_ylabel('Log(MWFFT Max. Amplitude)')
ax.set_title(f'Dependency with window overlap | $Q\'_v: {qpv}$ | window size: 64', fontweight='bold')
fig.savefig('MWFFTnoverlap3.png')
fig.tight_layout()
plt.show()

#2D-analysis
npersegs = np.array([32, 64, 128, 256, 512, 1024])
noverlaps = np.array([0.3, 0.5, 0.7, 0.9])
GR = np.zeros((len(npersegs), len(noverlaps)))
for nperseg in npersegs:
    for noverlap in noverlaps:
        #analysis
        gr = GrowthRate(file=file, nperseg=nperseg, noverlap=nperseg*noverlap, grstart=1500)
        gr.MWFFT_analysis()
        #gr.plot()
        gr.save(fname='MWFFTscan.csv', nperseg=nperseg, noverlap=noverlap)
        GR[npersegs == nperseg, noverlaps == noverlap] = gr.slope*1e3


# HEATMAP
# -------
fig, ax = plt.subplots(dpi=150)
GRdiff = GR - np.mean(GR)
vmin = np.min(GR[GR>0])*0.99
vmax = np.max(GR)*1.01
ax.imshow(GR.transpose(), vmin=vmin, vmax=vmax)
ax.set_xticks(np.arange(len(npersegs)), labels=npersegs)
ax.set_yticks(np.arange(len(noverlaps)), labels=noverlaps)
for i in range(len(npersegs)):
    for j in range(len(noverlaps)):
        ax.text(i, j, round(GR[i, j], 2), ha='center', va='center', color='k', alpha=0.5)

ax.set_xlabel('Points per window')
ax.set_ylabel('Window overlap')
ax.set_title(f'Heatmap of Growthrate for $Q\'_v={qpv}$')
ax.invert_yaxis()
fig.tight_layout()
fig.savefig('MWFFTheatmap.png')
fig.savefig('MWFFTheatmap.svg')
plt.show()

'''
#MWFFT colormaps
npersegs = np.array([32, 64, 128, 256, 512, 1024])
noverlaps = np.array([0.3, 0.5, 0.7, 0.9])
grvals = np.array([])
colors = ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e']
for nperseg in npersegs:
    for noverlap in noverlaps:
        #analysis
        gr = GrowthRate(file=file, nperseg=nperseg, noverlap=nperseg*noverlap, grstart=1000)
        gr.MWFFT_analysis()

        fig, ax = plt.subplots(dpi=100, figsize=(6,4))
        ax.imshow(np.log(np.abs(gr.MWFFT[2])), cmap='jet')
        ax.set_ylabel('Frequency')
        ax.set_xlabel('Window points')
        ax.set_title(f'log(|MWFFT|) with window size:${nperseg}$, overlap:${noverlap}$', fontweight='bold')
        fig.savefig(f'img/MWFFT_w{nperseg}_o{noverlap}.png')
        plt.show()
'''


