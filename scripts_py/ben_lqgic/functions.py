import itertools
import numpy as np
import scipy.ndimage.interpolation as inter
from scipy.ndimage.interpolation import rotate, shift
from scipy.ndimage import median_filter
from scipy.signal import welch, windows
from scipy.special import gamma,hyp2f1

p3i=lambda i: int(round(i)) #python2 to 3: change indicies that are floats to integers

def pupil_to_image(im):
	return np.fft.fft2(im,norm='ortho')

def image_to_pupil(im):
	return np.fft.ifft2(im,norm='ortho')
def mtf(im):
	return np.abs(np.fft.fftshift(image_to_pupil(im)))

def intensity(wavefront):
	return (np.abs(wavefront))**2.

def amplitude(wavefront):
	return np.abs(wavefront)

def phase(wavefront):
	return np.arctan2(np.imag(wavefront),np.real(wavefront))

def complex_amplitude(mag,phase):
	'''
	complex amplitude in terms of magnitude and phase
	'''
	return mag*np.cos(phase)+1j*mag*np.sin(phase)

def xy_plane(dim):
	'''
	define xy plane to use for future functions
	'''
	grid=np.mgrid[0:dim,0:dim]
	xy=np.sqrt((grid[0]-dim/2.+0.5)**2.+(grid[1]-dim/2.+0.5)**2.)
	return xy

def polar_grid(imagepix,pupilpix):
	'''
	make a polar image grid from a cartesian grid
	'''
	grid=np.mgrid[0:imagepix,0:imagepix]
	xy=np.sqrt((grid[0]-imagepix/2.+0.5)**2.+(grid[1]-imagepix/2.+0.5)**2.)
	xy[np.where(xy>pupilpix/2.)]=0.
	rad_norm=xy/np.max(xy)
	phi=np.arctan2(grid[1]-imagepix/2.+0.5,grid[0]-imagepix/2.+0.5)
	return rad_norm,phi

def lpf(im,s): 
	return median_filter(im,s)
def hpf(im,s):
	return im-lpf(im,s)

def mag_im(im_in,mag):
	'''
	magnify/demagnify image using cubic spline interpolation

	mag = magnification factor (numbers greater than 1 "zoom in," less than 1 "zoom out")
	'''
	imagepix=len(im_in) #assumes a square image
	grid=np.mgrid[0:imagepix,0:imagepix]
	maggrid=(grid-(imagepix/2))*1/mag+(imagepix/2)
	im_out=inter.map_coordinates(im_in,maggrid)
	return im_out

def make_kolmogorov_noise(wavefronterror,imagepix,wavelength):
	'''
	make kolmogorov noise with -5/3 power law:

	(1) take white noise as phase
	(2) the amplitude is a (-5/6) power law (i.e., weighting energy=amplitude**2. by -5/3 is the same as amplitude by -5/6 a la energy = amplitude**2)
	(3) combine (1) and (2) into a complex number and FFT image to pupil plane, take real part

	wavefronterror = rms WFE (nm)

	SOMEHOW THE -5/6 IN AMPLITUDE IS STILL WRONG, CHRISTIAN SAYS IT SHOULD BE -5/3, SO SOMETHING IS WRONG SOMEHWERE...
	'''

	white_noise=np.random.random((imagepix,imagepix))*2.*np.pi #phase
	#noisefft=np.fft.fftshift(np.fft.fft2(white_noise))
	xy=xy_plane(imagepix)
	#grid=np.mgrid[0:imagepix,0:imagepix]
	#xy=np.sqrt((grid[0])**2.+(grid[1])**2.)
	amplitude=(xy+1)**(-5./3.) #amplitude central value in xy grid is one, so max val is one in power law, everything else lower
	amp=shift(amplitude,(-imagepix/2.,-imagepix/2.),mode='wrap') #shift the peak amplitude from the center to all four corners of the image, which is where numpy weights the zero phase component in the fft
	image_wavefront=complex_amplitude(amp,white_noise) #combine random phase and power law amplitude into one complex number
	noise_wavefront=np.real(np.fft.fft2(image_wavefront))
	norm_factor=(wavefronterror/wavelength*2.*np.pi)/np.std(noise_wavefront) #normalization factor for phase error
	phase_out=noise_wavefront*norm_factor

	return phase_out

def make_noise(wavefronterror,imagepix,wavelength):
	'''
	make noise with -2 power law:

	(1) take white noise in image plane
	(2) IFFT to pupil plane
	(3) weight complex pupil plane (spatial frequency) by -1 power law (power~amplitude**2~x**(-2), so amplitude goes as (x**(-2))**(1/2)=x(-1)
	(4) FFT back to image plane

	wavefronterror = rms WFE (nm)
	'''

	white_noise=np.random.random((imagepix,imagepix))*2.*np.pi #phase
	#noisefft=np.fft.fftshift(np.fft.fft2(white_noise))
	xy=xy_plane(imagepix)
	#grid=np.mgrid[0:imagepix,0:imagepix]
	#xy=np.sqrt((grid[0])**2.+(grid[1])**2.)
	amplitude=(xy+1)**(-1.) #amplitude central value in xy grid is one, so max val is one in power law, everything else lower
	amplitude[p3i(imagepix/2),p3i(imagepix/2)]=0. #remove piston
	amp=shift(amplitude,(-imagepix/2.,-imagepix/2.),mode='wrap')
	image_wavefront=complex_amplitude(amp,white_noise)
	noise_wavefront=np.real(np.fft.fft2(image_wavefront))
	norm_factor=(wavefronterror/wavelength*2.*np.pi)/np.std(noise_wavefront) #normalization factor for phase error
	phase_out=noise_wavefront*norm_factor

	return phase_out

def zernike(n,m,rho,phi):
	'''
	make a zernike polynomial of specified n,m given input polar coordinate maps of rho (normalized to one; pupil coordinates only) and phi (radians)
	'''

	rad=gamma(n+1)*hyp2f1(-1./2.*(m+n),1./2.*(m-n),-n,rho**(-2))/gamma(1./2.*(2+n-m))/gamma(1./2.*(2+n+m))*rho**n
	if m>=0:
		cos=np.cos(m*phi)
		out=rad*cos
	else:
		sin=np.sin(-1*m*phi)
		out=rad*sin
	out[np.where(np.isnan(out)==True)]=0.
	return out

def remove_zernike(im_tar_phase,imagepix,pupilpix):
	'''
	remove zernikes from a phase screen up to n=5, which is the first 21 modes

	input:
	im_tar_phase = normal input phase screen generated from some power law
	imagepix = width of focal plane image in pixels
	pupilpix = width of pupil in pixels

	output: phase screen with first 21 zernikes removed from the pupil
	'''

	xy=xy_plane(imagepix)
	phase_screen=im_tar_phase-np.mean(im_tar_phase[np.where(xy<pupilpix/2.)]) #remove piston in the pupil just by changing to zero mean

	rho,phi=polar_grid(imagepix,pupilpix)

	zern_nm=[]
	for n in range(1,6): #remove zernikes up to n=5, which is the first 21 modes
		m=list(range(-1*n,n+2,2))
		for mm in m:
			zern_nm.append([n,mm])

	#reference array
	refarr=np.zeros((len(zern_nm),imagepix**2))
	for i in range(len(zern_nm)):
		z=zernike(zern_nm[i][0],zern_nm[i][1],rho,phi)
		refarr[i]=z.flatten()

	#covariance matrix:
	n=len(zern_nm)
	cov=np.zeros((n,n))
	for i in range(n):
		for j in range(i+1):
			if cov[j,i]==0.:
				cov[i,j]=np.sum(refarr[i,:]*refarr[j,:])
				cov[j,i]=cov[i,j]
			#print i*n+j,n**2-1
	covinv=np.linalg.pinv(cov,rcond=1e-7)

	#correlation image vector:
	tar=np.ndarray.flatten(phase_screen)
	cor=np.zeros((n,1))
	for i in range(n):
		cor[i]=np.sum(refarr[i]*tar)
		#print i, n-1

	coeffs=np.dot(covinv,cor)

	all_zern=np.dot(coeffs.T,refarr).reshape(imagepix,imagepix)
	out_phase=phase_screen-all_zern
	return out_phase

def remove_tt(im_tar_phase,imagepix,pupilpix):
	'''
	remove tip and tilt

	input:
	im_tar_phase = normal input phase screen generated from some power law
	imagepix = width of focal plane image in pixels
	pupilpix = width of pupil in pixels

	output: phase screen with tip and tilt removed from the pupil
	'''

	xy=xy_plane(imagepix)
	phase_screen=im_tar_phase-np.mean(im_tar_phase[np.where(xy<pupilpix/2.)]) #remove piston in the pupil just by changing to zero mean

	rho,phi=polar_grid(imagepix,pupilpix)

	zern_nm=[]
	for n in range(1,2): #remove tip,tilt zernikes
		m=list(range(-1*n,n+2,2))
		for mm in m:
			zern_nm.append([n,mm])

	#reference array
	refarr=np.zeros((len(zern_nm),imagepix**2))
	for i in range(len(zern_nm)):
		z=zernike(zern_nm[i][0],zern_nm[i][1],rho,phi)
		refarr[i]=z.flatten()

	#covariance matrix:
	n=len(zern_nm)
	cov=np.zeros((n,n))
	for i in range(n):
		for j in range(i+1):
			if cov[j,i]==0.:
				cov[i,j]=np.sum(refarr[i,:]*refarr[j,:])
				cov[j,i]=cov[i,j]
			#print i*n+j,n**2-1
	covinv=np.linalg.pinv(cov,rcond=1e-7)

	#correlation image vector:
	tar=np.ndarray.flatten(phase_screen)
	cor=np.zeros((n,1))
	for i in range(n):
		cor[i]=np.sum(refarr[i]*tar)
		#print i, n-1

	coeffs=np.dot(covinv,cor)

	all_zern=np.dot(coeffs.T,refarr).reshape(imagepix,imagepix)
	out_phase=phase_screen-all_zern
	return out_phase

def make_amp_err(percent,imagepix,pupilpix):
	'''
	make amplitude noise with -2 power law, but flat at low frequencies and zero at high frequencies:

	(1) take white noise in image plane
	(2) IFFT to pupil plane
	(3) weight complex pupil plane (spatial frequency) by -1 power law (power~amplitude**2~x**(-2), so amplitude goes as (x**(-2))**(1/2)=x(-1)
	(4) FFT back to image plane

	note that percent amplitude error is in intensity, so when I apply the percent variable, this is in intensity, but in the code it gets applied to amplitude, and amplitude = sqrt(intensity), hence the sqrt(percent) term

	'''

	white_noise=np.random.random((imagepix,imagepix))*2.*np.pi #phase
	#noisefft=np.fft.fftshift(np.fft.fft2(white_noise))
	xy=xy_plane(imagepix)
	#grid=np.mgrid[0:imagepix,0:imagepix]
	#xy=np.sqrt((grid[0])**2.+(grid[1])**2.)
	amplitude=(xy+1)**(-0.25) #amplitude central value in xy grid is one, so max val is one in power law, everything else lower
	
	#make all spatial frequencies below 1/Dpup/10 into a tophat function
	ind=np.where(np.abs(xy-(pupilpix/2.)/10.)==np.min(np.abs(xy-(pupilpix/2.)/10.)))
	amplitude=(xy+1)**(-1.)*xy[ind][0]
	amplitude[np.where(xy<(pupilpix/2.)/10.)]=1.
	amp_step_val=amplitude[np.where(np.abs(xy-(pupilpix/2.)/10.)==np.min(np.abs(xy-(pupilpix/2.)/10.)))][0]
	amplitude[np.where(xy<(pupilpix/2.)/10.)]=amp_step_val
	#amplitude[np.where(xy<(pupilpix/2.)/10.)]=0. #remove these spatial frequencies entirely instead of a tophat functions to see if this is limiting the SCC

	amplitude[p3i(imagepix/2),p3i(imagepix/2)]=0. #remove piston

	#remove alaising effects by cutting off power law just before the edge of the image
	amplitude[np.where(xy>imagepix/2.-1)]=0.

	amp=shift(amplitude,(-imagepix/2.,-imagepix/2.),mode='wrap')
	image_wavefront=complex_amplitude(amp,white_noise)
	noise_wavefront=np.real(np.fft.fft2(image_wavefront))
	norm_factor=np.sqrt(percent)/np.std(noise_wavefront[np.where(xy<pupilpix/2.)]) #normalize over the pupil, not the full image
	amp_out_ini=noise_wavefront*norm_factor

	amp_out=remove_tt(amp_out_ini,imagepix,pupilpix) #tip tilt removed phase screen
	#amp_out_gr1=remove_zernike(amp_out_ini,imagepix,pupilpix) #zernike removed phase screen; but, amp values are above one (when added to one below); change this so the max value is one
	#amp_out=amp_out_gr1-np.max(amp_out_gr1)


	return np.ones((imagepix,imagepix))+amp_out

def antialias(phin,imagepix,beam_ratio,Nact):
	'''
	anti-alias via a butterworth filter
	'''
	xy=xy_plane(imagepix)
	buttf = lambda rgrid,eps,r0,n: 1./np.sqrt(1+eps**2.*(xy/r0)**n) #butterworth filter
	phinput=phin-np.min(phin)
	phfilt=np.abs(pupil_to_image(np.fft.fftshift(image_to_pupil(phinput))*(buttf(xy,1,Nact/2*beam_ratio*0.99,100)))) #Nact actuators across the pupil
	phout=phfilt-np.mean(phfilt)
	return phout


def make_noise_pl(wavefronterror,imagepix,pupilpix,pl,Nact):
	'''
	make noise with a user input power law:

	(1) take white noise in image plane
	(2) IFFT to pupil plane
	(3) weight complex pupil plane (spatial frequency) by power law (power law is specified on the PSD, scaling is half that for amplitude)
	(4) FFT back to image plane and take the real part

	wavefronterror = rad rms WFE
	'''

	white_noise=np.random.random((imagepix,imagepix))*2.*np.pi #phase
	xy=xy_plane(imagepix)
	amplitude=(xy+1)**(pl/2.) #amplitude central value in xy grid is one, so max val is one in power law, everything else lower

	amplitude[p3i(imagepix/2),p3i(imagepix/2)]=0. #remove piston

	#remove alaising effects by cutting off power law just before the edge of the image
	amplitude[np.where(xy>imagepix/2.-1)]=0.

	amp=shift(amplitude,(-imagepix/2.,-imagepix/2.),mode='wrap')
	image_wavefront=complex_amplitude(amp,white_noise)
	noise_wavefront=np.real(np.fft.fft2(image_wavefront))

	beam_ratio=p3i(imagepix/pupilpix)

	norm_factor=wavefronterror/np.std(antialias(noise_wavefront,imagepix,beam_ratio,Nact)[np.where(xy<pupilpix/2.)]) #normalization factor for phase error over the pupil of modes within the DM control region
	phase_out_ini=noise_wavefront*norm_factor

	phase_out=phase_out_ini#remove_tt(phase_out_ini,imagepix,pupilpix) #tip tilt removed phase screen

	return phase_out

def make_wfe_atmosphere(r0,L0,imagepix,pupilpix,Dtel,rTTF=False):
	'''
	make open loop atmospheric WFE:

	(1) take white noise in image plane
	(2) IFFT to pupil plane
	(3) weight complex pupil plane (spatial frequency) by power law (power law is specified on the PSD, scaling is half that for amplitude)
	(4) FFT back to image plane and take the real part

	r0,L0 = Fried parameter (at observing wavelength) and outer scale, respectivley, in meters
	imagepix,pupilpix= image and pupil diameter, in pixels, respectively
	Dtel = telescope diameter, in meters
	if rTTF=True, the function also return the normalized phase screen WITH tip/tilt/focus and other terms < 0.5 c/p
	'''

	sfgrid=np.linspace(0.5,pupilpix/2,1000)/Dtel #spatial frequency grid over which to integrate to fine WFE; 0.5 c/p is the lowest spatial frequency the telescope can sample (sort of like focus or astigmatism). the highest spatial frequency (Nyquist limit) is pupilpix/2. c/p is converted into meters**(-1) by dividing by the telescope diameter
	wfe_norm=np.sqrt(0.49*r0**(-5/3)*np.trapz((sfgrid**2+(2.*np.pi/L0)**2)**(-11/6),x=sfgrid,dx=sfgrid[1]-sfgrid[0])) #WFE over the full pupil in rad rms (this is what the output WFE should be normalized to). paper reference: equation 6 in Johansson and Gavel, 1994, SPIE

	white_noise=np.random.random((imagepix,imagepix))*2.*np.pi #phase
	xy=xy_plane(imagepix)/Dtel #/Dtel converts xy radial grid into units of inverse meters

	amplitude=(xy**2+(2*np.pi/L0)**2)**((-11/6)/2) #divide by 2 in the outermost exponent to convert from PSD to amplitude power law
	#amplitude=xy**((-11/3)/2) #power law with no L0

	amplitude[np.where(xy<0.5)]=0. #remove piston and all spatial frequencies less than 0.5 c/p

	#remove alaising effects by cutting off power law just before the edge of the image
	amplitude[np.where(xy>imagepix/2.-1)]=0.

	amp=shift(amplitude,(-imagepix/2.,-imagepix/2.),mode='wrap')
	image_wavefront=complex_amplitude(amp,white_noise)
	noise_wavefront=np.real(np.fft.fft2(image_wavefront))

	beam_ratio=p3i(imagepix/pupilpix)
	noise_wavefront_rmtt=remove_tt(noise_wavefront,imagepix,pupilpix) #tip-tilt removed unnormalized phase screen
	noise_wavefront_rmttf=noise_wavefront_rmtt-antialias(noise_wavefront_rmtt,imagepix,beam_ratio,0.5) #remove spatial frequencies below 0.5 c/p
	norm_factor=wfe_norm/np.std(antialias(noise_wavefront_rmttf,imagepix,beam_ratio,pupilpix/2)[np.where(xy<pupilpix/2.)]) #normalization factor for phase error over the pupil of modes within the DM control region
	phase_out=noise_wavefront_rmttf*norm_factor
	if rTTF==False:
		return phase_out
	if rTTF==True:
		return phase_out, noise_wavefront*norm_factor

def translate_atm(phase_in,beam_ratio,Twfs,Dtel,vx=10,vy=0):
	'''
	specifying a non-default Twfs allows for a finer resolution in the bright star fast correction version below
	'''
	imagepix=phase_in.shape[0]
	grid=np.mgrid[0:imagepix,0:imagepix]
	fx,fy=(grid[1]-imagepix/2.)/beam_ratio/Dtel,(grid[0]-imagepix/2.)/beam_ratio/Dtel
	phase_alpha=-2.*np.pi*Twfs*(fx*vx+fy*vy)
	alpha=complex_amplitude(np.ones((imagepix,imagepix)),phase_alpha)
	alpha_ft=pupil_to_image(phase_in)*np.fft.fftshift(alpha)
	next_phase=np.real(image_to_pupil(alpha_ft))
	return next_phase

#modal gain optimization functions
def detrend(coeffs):
	x=np.linspace(0,1,len(coeffs))
	fit_coeffs=np.poly1d(np.polyfit(x,coeffs,1))
	coeffs_detrend=coeffs-fit_coeffs(x)
	return coeffs_detrend

#AO transfer functions
Hwfs = lambda s, Ts: (1. - np.exp(-Ts*s))/(Ts*s)
Hzoh=Hwfs
Hlag = lambda s,tau: np.exp(-tau*s)
Hint = lambda s, Ts, g: g/(1. - np.exp(-Ts*s)) #integrator
Hlint=  lambda s, Ts, g, l: g/(1. - l*np.exp(-Ts*s)) #leaky integrator
Hcont = Hlint
Holsplane = lambda s, Ts, tau, g, l:  Hwfs(s, Ts)*Hlag(s,tau)*Hcont(s, Ts, g, l)*Hzoh(s,Ts)
s2f = lambda f: 1.j*2.*np.pi*f
Hol = lambda f, Ts, tau, g, l:  Holsplane(s2f(f),Ts,tau,g,l) #open loop transfer function
Hrej = lambda f, Ts, tau, g, l: 1./(1. + Hol(f, Ts, tau, g, l)) #rejection transfer function
Hcl = lambda f, Ts, tau, g, l: Hol(f, Ts, tau, g, l)/(1. + Hol(f, Ts, tau, g, l))
Hn = lambda f, Ts, tau, g, l: Hcl(f, Ts, tau, g, l)/Hwfs(s2f(f), Ts)

def genpsd(tseries,dt,nseg=4):
	nperseg=2**int(np.log2(tseries.shape[0]/nseg)) #firstly ensures that nperseg is a power of 2, secondly ensures that there are at least nseg segments per total time series length for noise averaging
	window=windows.hann(nperseg)
	freq,psd=welch(tseries,fs=1./dt,window=window,noverlap=nperseg*0.25,nperseg=nperseg,detrend=False,scaling='density')
	freq,psd=freq[1:], psd[1:] #remove DC component (freq=0 Hz)
	return freq,psd