'''
building on toy_1d_ic_loop_subframes.py, make a simple 1d toy model of an LQG Integral Controller LQG IC, using Lisa's matrix-based LQG formalism
'''
# %%
import os
import sys
import itertools
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.interpolation import rotate, shift
import scipy.ndimage.interpolation as inter
import scipy
import multiprocessing as mp
from glob import glob
#utility functions to use throughout the simulation
import functions
import matplotlib as mpl

p3i=lambda i: int(round(i)) #python2 to 3: change indicies that are floats to integers

nframes=1000

#CREATE LOOP WITH AN INTEGER FRAME DELAY
#free variables
frame_delay=1 #LQG is currently only setup for 1 frame
gain, leak = 0.43,0.999

alpha,noisevar=0.999,0.05 #LQG free parameters: AR model coeff of the temporal wavefront error's signal and noise terms


#LQG setup:
n_phase_past=1 #number of phase states
nd=n_phase_past + 1 #number of DM commands states. The "+1" is because you are predicting d_{t+1} into the future
n_ab=1 #number of AR coeffs describing phase
ns=n_phase_past+n_ab #number of states ignoring DM commands
nt=nd+ns #total number of states including DM commands

n_states=nt
Amat=np.zeros((ns,ns))
Amat[0,0]=alpha
Amat[1, 0] = 1.
nw = 1
Bmat = np.zeros((ns, nw))
Bmat[0,0] = 1
sigW=np.array([[1.]])

Cmat = np.zeros((1, ns))
Dmat = np.zeros((1, nd))    
f=frame_delay - 0 #this if for frame_delay <= 1
Cmat[0,(ns-2):(ns)] = [(1-f), f]
Dmat[0,(nd-2):(nd)] = [-(1-f), -f]

# Now convert to Looze form.... the prefix "l" denotes that this is a Looze form matrix
Cmat_pure = np.zeros((1, ns))
Cmat_pure[0,ns-1] = 1
Dmat_pure = np.zeros((1, nd))
Dmat_pure[0,nd-1] = -1
# like in Looze - actuators then aberration need two states
# https://numpy.org/doc/stable/reference/generated/numpy.tril.html
lowerdiag = np.tril(np.ones((nd, nd)), -1) - np.tril(np.ones((nd, nd)), -2)
lAmat = np.zeros((nt, nt))
lAmat[0:nd, 0:nd] = lowerdiag
lAmat[nd:nt, nd:nt] = Amat

lBmat = np.zeros((nt, nw))
lBmat[nd:nt,:] = Bmat

lCmat = np.zeros((1, nt))
lCmat[0,0:nd] = Dmat
lCmat[0,nd:nt] = Cmat

lDmat = np.zeros((nt, 1))    
lDmat[0,0] = 1

lsigW = sigW

lCmat_pure = np.append(Dmat_pure, Cmat_pure, axis=1)
lQmat = np.dot(np.transpose(lCmat_pure),lCmat_pure)

from scipy.linalg import eigvals, solve_discrete_are, solve
def mydare(A, B, Q, R):
        X = solve_discrete_are(A, B, Q, R)
        G = solve(B.T.dot(X).dot(B) + R, B.T.dot(X).dot(A))
        L = eigvals(A - B.dot(G))
        return X, L, G

def model2gainvecs(lAmat, lCmat, lDmat, lsigW, lsigv, lQmat, lRmat):
    [l1X, l1L, l1G] = mydare(lAmat, lDmat, lQmat, lRmat)
    [l2X, l2L, l2G] = mydare(np.conj(np.transpose(lAmat)), np.transpose(lCmat), lsigW, lsigv)

    lGmat = -np.linalg.inv(np.dot(np.dot(np.transpose(lDmat), l1X), lDmat) + lRmat)*(np.dot(np.dot(np.transpose(lDmat), l1X), lAmat))
    lKmat = np.dot(np.dot(l2X, np.transpose(lCmat)),np.linalg.inv(np.dot(np.dot(lCmat, l2X), np.transpose(lCmat)) + lsigv))
    return lGmat, lKmat

def helper_make_control(lqg_lAmat, lqg_lCmat, lqg_lDmat, lqg_lsigW, lqg_lQmat, noisevar):
        lRmat = np.zeros((1, 1))
        lEmat = np.zeros((nt, nt))
        for t in range(nt):
            lEmat[t,t] = 1

        [lGmat, lKmat] = model2gainvecs(lqg_lAmat, lqg_lCmat, lqg_lDmat, lqg_lsigW, noisevar, lqg_lQmat, lRmat)
        lqg_k = lKmat
        lqg_ikca = np.dot((lEmat - np.dot(lKmat, lqg_lCmat)), lqg_lAmat)
        lqg_ikcd = np.dot((lEmat - np.dot(lKmat, lqg_lCmat)), lqg_lDmat)
        lqg_g = lGmat
        return lqg_k,lqg_ikca,lqg_ikcd,lqg_g


mat_k, mat_ikca, mat_ikcd, mat_g = helper_make_control(lAmat, lCmat, lDmat,np.matmul(np.matmul(lBmat,lsigW), np.transpose(np.conj(lBmat))),lQmat, noisevar)

#IC version (non-LQG):
K=np.ones((n_states,1))*gain
ImKCtA=np.ones((n_states,n_states))*leak #(I-KC)A
ImKCtD=np.zeros((n_states,1)) #(I-KC)D
G=np.ones((1,n_states))

#LQG IC version
K=mat_k
ImKCtA= mat_ikca #(I-KC)A
ImKCtD= mat_ikcd#(I-KC)D
G=mat_g

#generate LQG ETFs:
def genlqg_tfs(freq):
	ldims = np.shape(mat_ikca)
	nstates = ldims[0]
	eye_n = np.zeros((nstates, nstates))
	for k1 in range(nstates):
	        eye_n[k1,k1] = 1
	cofz_lqg = np.zeros((len(freq)), complex)
	zinv = np.exp(-Twfs * 1j * 2 * np.pi * freq)
	for zz in range(len(freq)):
	        this_zinv = zinv[zz]           
	        mat3 = np.linalg.inv(eye_n - mat_ikca*this_zinv) 
	        leftside = 1. - np.dot(np.dot(mat_g, mat3), mat_ikcd)*this_zinv 
	        rightside = np.dot(np.dot(mat_g, mat3), mat_k)
	        cofz_lqg[zz] = rightside/leftside               

	Holsplane_nocont = lambda s: functions.Hwfs(s,Twfs)*functions.Hlag(s,Twfs*frame_delay)*functions.Hzoh(s,Twfs)
	Hol_nocont=Holsplane_nocont(functions.s2f(freq))
	Hol_lqg=Hol_nocont*cofz_lqg
	etf_lqg = np.abs(1 - Hol_lqg/(1. + Hol_lqg))**2
	return etf_lqg,Hol_lqg

def getbw(etf): #get bandwidth
        indsbw=np.where(np.diff(np.sign(etf-1)))[0] #find the zero crossings, including ones higher than the BW
        bw=np.min(freq[indsbw]) #the BW is the minimum frequency
        return bw


size=15
font = {'family' : 'Times New Roman',
        'size'   : size}

mpl.rc('font', **font)
mpl.rcParams['image.interpolation'] = 'nearest'

from matplotlib.lines import Line2D
from matplotlib import cm

fig=plt.figure(figsize=(5,5))
ax1t = plt.gca()

colors=cm.viridis(np.linspace(0,0.75,3))

Twfs=1e-3
tau=frame_delay*Twfs
coldstart=100 #how many frames to discard due to closing the loop at first
pl=-3#-11/3 #Kolmogorov?
white_noise=np.random.random(nframes)*2.*np.pi #phase
grid=np.mgrid[0:nframes]
xy=np.sqrt((grid-nframes/2.+0.5)**2.)
amplitude=(xy+1)**(pl/2.) #amplitude central value in xy grid is one, so max val is one in power law, everything else lower

amplitude[p3i(nframes/2)]=0. #remove piston

amp=shift(amplitude,-nframes/2,mode='wrap')
complex_fourier_domain=functions.complex_amplitude(amp,white_noise)
open_loop_time=np.real(np.fft.fft(complex_fourier_domain))
freq,psdol=functions.genpsd(open_loop_time[coldstart:],Twfs)

ax1t.set_yscale('log')
ax1t.set_xscale('log')
etf_ic=np.abs(functions.Hrej(freq, Twfs, tau, gain, leak))**2
bw_ic=getbw(etf_ic)
ax1t.plot(freq,etf_ic,ls='--',color=colors[1])
etf_lqg,Hol_lqg=genlqg_tfs(freq)
bw_lqg=getbw(etf_lqg)
ax1t.plot(freq,etf_lqg,ls='--',color=colors[2])
ax1t.set_ylabel('|rejection transfer function (RTF)|$^2$')
ax1t.legend(loc='lower right',framealpha=0.5)
ax1t.set_title('IC BW ='+str(np.round(bw_ic,1))+' Hz, LQG IC BW ='+str(np.round(bw_lqg,1))+'Hz',size=size)
plt.show()
# %%
ldims = np.shape(mat_ikca)
nstates = ldims[0]
eye_n = np.zeros((nstates, nstates))
for k1 in range(nstates):
        eye_n[k1,k1] = 1
cofz_lqg = np.zeros((len(freq)), complex)
zinv = np.exp(-Twfs * 1j * 2 * np.pi * freq)
zinv
# %%
for zz in range(len(freq)):
        this_zinv = zinv[zz]           
        mat3 = np.linalg.inv(eye_n - mat_ikca*this_zinv) 
        leftside = 1. - np.dot(np.dot(mat_g, mat3), mat_ikcd)*this_zinv 
        rightside = np.dot(np.dot(mat_g, mat3), mat_k)
        cofz_lqg[zz] = rightside/leftside     
# %%
