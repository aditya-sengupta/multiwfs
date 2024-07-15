'''
make a simple 1d toy model of an integrator on disturbence

building on toy_1d_ic_loop.py, here add subframes to allow a a frame delay that is an non-integer number of the frame rate... this takes longer to simulate
'''

# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from scipy.ndimage import shift
#utility functions to use throughout the simulation
import functions

nframes=10000
fpf=10 #number of subframes per frame
nsubframes=nframes*fpf

#GENERATE NOISELESS OPEN LOOP DISTURBANCE
pl=-11/3 #Kolmogorov?
white_noise=np.random.random(nsubframes)*2.*np.pi #phase
grid=np.mgrid[0:nsubframes]
xy=np.sqrt((grid-nsubframes/2.+0.5)**2.)
amplitude=(xy+1)**(pl/2.) #amplitude central value in xy grid is one, so max val is one in power law, everything else lower

amplitude[nsubframes//2]=0. #remove piston

amp=shift(amplitude,-nsubframes/2,mode='wrap')
complex_fourier_domain=functions.complex_amplitude(amp,white_noise)
real_time_domain=np.real(np.fft.fft(complex_fourier_domain))

#GENERATE A NOISE COMPONENT TO THE TURBULENCE
noise2signal=0.01 #std of noise relative to std of signal
noise=np.random.randn(nsubframes)
noise=noise*np.std(real_time_domain)/np.std(noise)*noise2signal
real_time_domain=real_time_domain+noise


#CREATE LOOP WITH AN NON-INTEGER FRAME DELAY OF THE WFS FRAME RATE
frame_delay_want=0.5 #requested frame delay
frame_delay = round((frame_delay_want-np.floor(frame_delay_want))*fpf)/fpf+np.floor(frame_delay_want) #the actual frame delay should be based on the integer number of subframes, so that the exact number of subframes can be calculated based in the exact frame delay
nsubframes_per_loop=int((frame_delay+2)*fpf) #how many subframes to save throughput the loop iteration
dmc_arr=np.zeros(nsubframes_per_loop) #in this subframe model the array of DM commands is the total number of subframes in between the first subframe at the start of the WFS exposure and the last subframe and the end of the DM commands.
ph_real_arr=np.zeros(fpf) #for open loop phase realizations each frame, you only care about what happens during the WFS exposure, so you only need fpf frames
gain, leak = 0.2,0.99
last_dmc=0
open_loop_time=np.zeros(nframes)
closed_loop_time=np.zeros(open_loop_time.shape)
dmc_time=np.zeros(open_loop_time.shape)
for i in range(nframes):
	for j in range(fpf):
		ph_real_arr[j]=real_time_domain[i*fpf+j]

	
	dmc_recon=np.mean(ph_real_arr-dmc_arr[:fpf],axis=0) #assume perfect measurement of the residual phase; just looking at latency effects of here

	#LQG STEPS:
	'''
	1. dmc_recon is y_t
	2. 1 [1x1] is dotted with K [1 x n_states]
	3. 2 [1 x n_states] is summed with dot[x_(t-1|t-1),(I-KC)A], producing an output with dimensions [1x n_states]
		- x_(t-1|t-1) is dimensions [1 x n_states]
		- (I-KC)A is dimensions [n_states x n_states]
	3a. once you have x_(t|t), 
	4. then sum output of 3 with u_(t-1) * (I-KC)D
		- (I-KC)D is [n_states x 1]
		- u_(t-1) is a scalar
	5. multiply the output of step 4 [n_states x 1] with G to give d_(t+1)
		- G is dimensions [1 x n_states]
		- d_(t+1) is a scalar
	
	'''
	
	dmc_apply=last_dmc*leak+dmc_recon*gain

	#save data
	open_loop_time[i]=np.mean(ph_real_arr,axis=0)
	closed_loop_time[i]=np.mean(ph_real_arr-dmc_arr[:fpf],axis=0)
	dmc_time[i]=-1*np.mean(dmc_arr[:fpf],axis=0)

	
	#reset for the next loop
	dmc_arr[-fpf:]=dmc_apply
	dmc_arr=np.roll(dmc_arr,-fpf) #shift the order of the dmc_arr so that eventually the applied DM commands generated in this frame are applied in a later frame based on the give frame delay 
	last_dmc=dmc_apply


#PLOT THE RESULTS

Twfs=1e-3
tau=frame_delay*Twfs
coldstart=100 #how many frames to discard due to closing the loop at first
freq,psdol=functions.genpsd(open_loop_time[coldstart:],Twfs)
freq,psdcl=functions.genpsd(closed_loop_time[coldstart:],Twfs)

tarr=np.linspace(0,nframes*Twfs,nframes)
size=15
font = {'family' : 'Times New Roman',
        'size'   : size}

mpl.rc('font', **font)
mpl.rcParams['image.interpolation'] = 'nearest'

colors=cm.viridis(np.linspace(0,0.5,2))

fig,axs=plt.subplots(nrows=1,ncols=2,figsize=(10,5))
ax1,ax2=axs.ravel()
ax1.plot(tarr,open_loop_time,color=colors[0],label='OL')
ax1.plot(tarr,closed_loop_time,color=colors[1],label='CL')
ax1.plot(tarr,dmc_time,color='k',ls='--',label='DMC')
ax1.set_ylabel('modal coefficient (arbitrary units)')
ax1.set_xlabel('time (s)')
ax1.legend(loc='best')

coldstart=10 #how many frames to discard due to closing the loop at first
ax2.plot(freq,psdol,color=colors[0],label='OL PSD')
ax2.plot(freq,psdcl,color=colors[1],label='CL PSD')
ax2.set_xlabel('temporal frequency (HZ)')
ax2.set_ylabel('PSD (arbitrary units)')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.legend(loc='upper left',framealpha=0.5)

axt2=ax2.twinx()
axt2.set_yscale('log')
axt2.set_xscale('log')
axt2.plot(freq,psdcl/psdol,color='k',label='(CL/OL) |RTF|$^2$')
axt2.plot(freq,np.abs(functions.Hrej(freq, Twfs, tau, gain, leak)**2),ls='--',color='k',label='model |RTF|$^2$')
axt2.set_ylabel('|rejection transfer function (RTF)|$^2$')
axt2.legend(loc='center left',framealpha=0.5)

fig.suptitle('$\\frac{\\sigma_\\mathrm{noise}}{\\sigma_\\mathrm{signal}}$='+str(noise2signal)+', T$_s$='+str(Twfs*1e3)+'ms, $\\tau$='+str(np.round(Twfs*frame_delay*1e3,1))+'ms, gain='+str(gain)+', leak='+str(leak),size=size)
plt.tight_layout()


#look at Hrej stability

colors=cm.viridis(np.linspace(0,0.75,4))

fig,ax=plt.subplots(figsize=(5,5))
ax.grid('on')
ax.axvline(0,color='k',ls='--')
ax.axhline(0,color='k',ls='--')
linfreq=np.linspace(np.min(freq),np.max(freq),1000)
Hol2plotp=functions.Hol(linfreq, Twfs, Twfs*frame_delay, gain, leak)
Hol2plotn=functions.Hol(-1*linfreq, Twfs, Twfs*frame_delay, gain, leak)
ax.plot(np.real(Hol2plotp),np.imag(Hol2plotp),color=colors[1])
ax.plot(np.real(Hol2plotn),np.imag(Hol2plotn),color=colors[1])
phasegrid=np.linspace(-np.pi,np.pi,500)
xunit,yunit=np.cos(phasegrid),np.sin(phasegrid)
ax.plot(xunit,yunit,color=colors[2],ls=':')
gain_margin_limit=2.5
ax.axvline(-1/gain_margin_limit,color=colors[0],ls='--')
#NEXT: FIND WHERE THE UNIT CIRCLE INTERSECTS THE BLUE CURVE
#np.where(np.logical_and)
ax.plot(-1,0,'-ro')
ax.set_xlim(-1.1,1.1)
ax.set_ylim(-1.1,1.1)
ax.plot(np.linspace(-2,0,10),np.linspace(-2,0,10),color=colors[3],ls='--')
ax.plot(np.linspace(-2,0,10),np.linspace(2,0,10),color=colors[3],ls='--')

phase_margin_point=np.abs(np.real(Hol2plotp)**2+np.imag(Hol2plotp)**2-1)
phase_margin_ind=np.where(phase_margin_point==np.min(phase_margin_point))[0][0]
ax.plot(np.real(Hol2plotp[phase_margin_ind]),np.imag(Hol2plotp[phase_margin_ind]),'-o',color=colors[3])
# %%
