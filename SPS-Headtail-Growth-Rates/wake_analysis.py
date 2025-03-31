import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scienceplots


def impedance(time, wake):
    
    ft = np.fft.fft(wake, n=len(wake))
    freq = np.fft.fftfreq(len(wake), d=(time[1]-time[0]))
    ft = ft[freq>=0]
    freq = freq[freq>=0]

    return freq, np.real(ft) , np.imag(ft)


#Get wake data

plt.style.use(['science','no-latex', 'ieee'])

dold=np.transpose(np.loadtxt('old/SPS_Q26.wake', skiprows=1))
dnew_corr=np.transpose(np.loadtxt('old/new_model_corr.wake', skiprows=1))
dnew=np.transpose(np.loadtxt('old/SPS_wake_new.wake', skiprows=1))



tnew=dnew[0]
dipnew=dnew[2]
quadnew=dnew[4]

print(quadnew)

tcorr=dnew_corr[0]
dipcorr=dnew_corr[2]
quadcorr=dnew_corr[4]

told=dold[0]
dipold=dold[2]
quadold=dold[4]




#Print dipole wake function


fig, axes = plt.subplots(2,2, dpi=300, figsize=(3.25, 3.25))
fig.tight_layout()

#axes[0,0].set_xlabel(r"Time [ns]")
axes[0,0].set_ylabel(r'Wake function [V/pC/mm]')

axes[0,0].plot(tnew,dipnew, 'b', alpha=0.8, linestyle='-', ms=1.5, lw=0.8, label='new model: dip')
axes[0,1].plot(tnew,quadnew,'b', alpha=0.8,linestyle='-',ms=1.5, lw=0.8, label='new model: quad')


axes[0,0].plot(tcorr,dipcorr,'r', alpha=0.8, linestyle='-',ms=1.5, lw=0.8, label='corrected new model: dip')
axes[0,1].plot(tcorr,quadcorr,'r', alpha=0.8,linestyle='--',ms=1.5, lw=0.8, label='corrected new model: quad')



axes[0,0].plot(told,dipold,'g', alpha=0.8, linestyle='-',ms=1.5, lw=0.8, label='old model: dip')
axes[0,1].plot(told,quadold,'g', alpha=0.8, linestyle='--',ms=1.5, lw=0.8, label='old model: quad')




axes[0,0].set_xlim(left=-0.1, right=25)
axes[0,1].set_xlim(left=-0.1, right=25)


axes[0,0].title.set_text('Dipole wake function')
axes[0,1].title.set_text('Quadrupole wake function')

axes[0,0].legend()
axes[0,1].legend()
axes[1,0].legend()
axes[1,1].legend()




#ax.legend()


#printing dipole impedance

freq,real, imag =impedance(tnew, dipnew)
axes[1,0].plot(freq,real, '-b', label='new wall impedance' )

freq2,real2, imag2 =impedance(tcorr, dipcorr)
axes[1,0].plot(freq2,real2, '-r' ,label='corrected new wall impedance')

freq3,real3, imag3 =impedance(told, dipold)
axes[1,0].plot(freq3,real3, '-g', label='old wall impedance')

#freq4,real4, imag4 =impedance(tnew, dipnew)
#axes[1,0].plot(freq4,real4, '-m', label='impedance new')

#axx.plot(time,imag, ':')
axes[1,0].set_xlim(left=-0.1, right=5)
axes[1,0].set_xlabel('frequency [GHz]')
axes[1,0].set_ylabel('Transverse Impedance [Ohm/m]')


axes[1,0].title.set_text('Dipole Impedance')


#printing quadrupole impedance

freq,real, imag =impedance(tnew, quadnew)
axes[1,1].plot(freq,real, '-b',linestyle='-', label='new wall impedance')

freq2,real2, imag2 =impedance(tcorr, quadcorr)
axes[1,1].plot(freq2,real2, '-r',linestyle='--', label='corrected new wall impedance')

freq3,real3, imag3 =impedance(told, quadold)
axes[1,1].plot(freq3,real3, '-g', label='old wall impedance')

#freq4,real4, imag4 =impedance(tnew, quadnew)
#axes[1,1].plot(freq4,real4, '-m', label='impedance new')

#axx.plot(time,imag, ':')
axes[1,1].set_xlim(left=-0.1, right=5)
axes[1,1].set_xlabel('frequency [GHz]')

axes[1,1].title.set_text('Quadrupole Impedance')


plt.legend()
plt.show()


