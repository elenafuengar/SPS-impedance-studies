import numpy as np
import matplotlib.pyplot as plt
import scienceplots
from scipy.constants import c
from scipy.interpolate import interp1d
import scienceplots

plt.style.use(['science','no-latex', 'ieee'])
fig, ax = plt.subplots(1,1,dpi=300, figsize=(3.25, 3.25))
fig.tight_layout(rect=[0, 0.12, 1, 1])
def impedance(time, wake):
    
    ft = -1.j*np.fft.fft(wake, n=len(wake))
    freq = np.fft.fftfreq(len(wake), d=(time[1]-time[0])*1e-9)
    #ft=np.fft.fftshift(ft)
    #freq=np.fft.fftshift(freq)
    freq = freq
    return freq, np.real(ft) , np.imag(ft)

wake_new_corr=np.transpose(np.loadtxt('SPS_Q26.wake', skiprows=1))
tcorr=wake_new_corr[0]
dipcorr=wake_new_corr[2]
quadcorr=wake_new_corr[4]
freq,real, imag =impedance(tcorr, dipcorr)
freq,realquad, imagquad = impedance(tcorr,quadcorr)
imagtot=imag+imagquad
Z=interp1d(freq, imagtot)

Q=26
omega0=2*np.pi*43.4e3 #s^-1
omega_beta=omega0*Q #Tune

chromas=np.linspace(-0,1.8,50)

lenghts=[2.5,3,3.5,4,4.5, 5, 6]

for bl_val_ns in lenghts:
    y=[]
    for chroma in chromas:

        num=0
        den=0
        for p in np.linspace(-150000,150000):
            omega_xi=omega0*Q*chroma/(0.62e-3*2)
            omega_prime=omega0*p #+l*omega_s pero estamos en l=0
            
            sigma_z = bl_val_ns*1.e-9*c/4 # [m]

            omega=omega_prime+omega_beta-omega_xi
            bunchlenght=np.exp(-1*(omega*omega*sigma_z*sigma_z)/(c*c))
            num=num+Z(omega_prime+omega_beta)*bunchlenght
            
            den=den+bunchlenght

        eff_impedance=num/den
        y.append(eff_impedance)

    ax.plot(chromas,y, ms=1.5, c=plt.cm.jet([bl_val_ns/6]), lw=0.8, label=f'bunch length: {bl_val_ns} ns')
ax.legend()
#ax.invert_xaxis()
ax.set_xlabel(r"Chromaticity $Q_V'/ Q_V$")
ax.set_ylabel(r"Effective impedance")
plt.show()



