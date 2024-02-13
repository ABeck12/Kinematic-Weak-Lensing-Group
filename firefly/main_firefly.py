import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.modeling import models
from specutils import Spectrum1D, SpectralRegion
from specutils.analysis import correlation

import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit

from binning_1d import binning_1d
from load_templates import read_templates
dic_templates = read_templates()

# rest_value = 6000. * u.AA # what's this for?
obsfilename  = 'spec1d.m46.023.A2552.fits'
binfactor    = 16
templatename = 'VVDS Elliptical'

# **observed** spectrum (obs_spec)
hdulist = fits.open(obsfilename)
print('Opened...' + obsfilename)
hdr01data = hdulist[1].data
w = hdr01data["LAMBDA"][0]
f = hdr01data["FLUX"][0]
e = hdr01data["IVAR"][0]
arr_obs_wave = binning_1d(w, binfactor) * u.AA
arr_obs_flux = binning_1d(f, binfactor) * u.dimensionless_unscaled
arr_obs_erro = binning_1d(e, binfactor) * u.dimensionless_unscaled


# **template** spectrum (tem_spec)
arr_tem_wave = dic_templates[templatename][0] * u.AA
arr_tem_flux = dic_templates[templatename][1] * u.Unit("erg / (Angstrom s cm2)")
arr_tem_erro = 0.2e-19 * np.ones(len(arr_tem_flux))* u.Unit("erg / (Angstrom s cm2)")


# Apply the same small flux error for both spectra
uncertainty1 = StdDevUncertainty(arr_obs_erro)
uncertainty2 = StdDevUncertainty(arr_tem_erro)
arr_obs_uncty_obj = uncertainty1
arr_tem_uncty_obj = uncertainty2



# template rescale and shift up/down
def func_shift(x, k, b):
    y = k * x + b
    return y
arr_tem_wave_common = arr_tem_wave[(arr_tem_wave >= min(arr_obs_wave)) & 
                                   (arr_tem_wave <= max(arr_obs_wave))]
arr_tem_flux_common = arr_tem_flux[(arr_tem_wave >= min(arr_obs_wave)) & 
                                   (arr_tem_wave <= max(arr_obs_wave))]
interpo_func        = interpolate.interp1d(arr_obs_wave, arr_obs_flux, kind='linear') # f_obs(x_obs)
arr_obs_flux_resamp = interpo_func(arr_tem_wave_common) # resampled y_obs at x_tem; len=len(tem) 
h_ratio = np.max(arr_obs_flux_resamp) / np.max(arr_tem_flux_common).value
ver_off = np.min(arr_obs_flux_resamp) - np.min(arr_tem_flux_common).value
popt, pcov    = curve_fit(func_shift, 
                          arr_tem_flux_common, 
                          arr_obs_flux_resamp,
                          bounds=([0.8*h_ratio, .1*ver_off], 
                                  [5.0*h_ratio, 10*ver_off]))
best_k, best_b = popt[0]*u.dimensionless_unscaled, popt[1]*u.Unit("erg / (Angstrom s cm2)")
arr_tem_flux   = best_k * arr_tem_flux + best_b


# Assemble [wavelength, spec, error] into the spec objects
""" 
    spectral_axis: array of float, quantity with unit u.AA
    flux:          array of float, quantity with unit u.Jy or defined
    uncertainty:   array of float, quantity with unit u.Jy or defined
"""
obs_spec = Spectrum1D(spectral_axis = arr_obs_wave, 
                      flux          = arr_obs_flux, 
                      uncertainty   = arr_obs_uncty_obj)
tem_spec = Spectrum1D(spectral_axis = arr_tem_wave, 
                      flux          = arr_tem_flux, 
                      uncertainty   = arr_tem_uncty_obj)


# Template cross correlation
corr, lag = correlation.template_correlate(observed_spectrum = obs_spec, 
                                           template_spectrum = tem_spec,
                                           lag_units=u.dimensionless_unscaled)
                                     # or: lag_units=Unit('km / s')


z_peak = lag[np.argmax(corr)].value
z_guess = z_peak + 0.
arr_tem_wave *= (1+z_guess)




# All plottings
fig = plt.figure(figsize=(16,9),dpi=200)
plt.subplots_adjust(hspace=0.2, wspace=0.2)
gs = fig.add_gridspec(2, 2,
                      height_ratios=[3, 2],
                      width_ratios=[3, 1])
ax1 = fig.add_subplot(gs[0, :])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[1, 1])
ax1.set_title(obsfilename, fontsize=24, loc='left')

# Sub-plot for observed and template spectra
spec_axis1, spec_axis2 = arr_obs_wave.value, arr_tem_wave.value
flux1,      flux2      = arr_obs_flux.value, arr_tem_flux.value
uncertainty_arrayvalues= arr_obs_erro.value
ax1.grid(alpha=1, linestyle='dotted', zorder=-1)
ax1.plot(spec_axis1, 
        flux1, 
        label='Observed', 
        color='black', linewidth=1, zorder=11)
ax1.plot(spec_axis2, 
        flux2, 
        label='Template: '+templatename+' at z='+"{:.4f}".format(z_guess),
        color='red', linewidth=1,  zorder=11, linestyle=':')
ax1.fill_between(spec_axis1, 
                uncertainty_arrayvalues, 
                np.zeros(len(arr_obs_flux)), 
                color='black', linewidth=0, alpha=.15, label='Flux uncertainty', zorder=9)
ax1.set_ylabel('Flux')
ax1.minorticks_on()
ax1.tick_params(  which='minor', length=3,bottom=True,top=False,direction='in')
ax1.tick_params(  which='major', length=6,bottom=True,top=False,direction='in')
ax1.set_xlabel(r'Observed Wavelength ($\AA$)')
ax1.set_xlim(arr_obs_wave[0].value, arr_obs_wave[-1].value)
ax1.legend(loc='upper left', prop={'size': 15})

# Sub-plot for redshift probability distribution and guess
corr_normalized = corr/sum(corr) # normalization

ax2.grid(alpha=1, linestyle='dotted', zorder=-1)
ax3.grid(alpha=1, linestyle='dotted', zorder=-1)
ax2.plot(lag, 
        corr_normalized, 
        color='green', linewidth=2, label='Object 1', zorder=10)
ax3.plot(lag, 
        corr_normalized, 
        color='green', linewidth=3, zorder=10)
ax2.set_ylabel('Correlation')
ax3.set_ylabel('Correlation')
ax2.minorticks_on()
ax3.minorticks_on()
ax2.tick_params(which='minor', length=3, bottom=True,top=False,direction='in')
ax2.tick_params(which='major', length=9, bottom=True,top=False,direction='in')
ax3.tick_params(which='minor', length=3, bottom=True,top=True, direction='in')
ax3.tick_params(which='major', length=15,bottom=True,top=True, direction='inout',
                axis='x', labeltop='on', labelbottom='on', labelsize=12)
ax2.set_xlabel(r'Redshift ($z$)')
ax3.set_xlabel(r'Redshift ($z$)')
ax2.legend(loc='upper left',prop={'size': 12})

# Sub-plot for zoom-in redshift peak
try:
    z_peak = lag[np.argmax(corr_normalized)]
    cpdf   = np.array([0])
    for i in range(len(corr_normalized)): cpdf = np.append(cpdf, cpdf[i]+corr_normalized[i])
    cpdf   = cpdf[1:]
    i_left, i_mean, i_right = None, None, None
    for i in range(len(cpdf)):
        if cpdf[i] >= 0.1587 and  i_left==None:
            i_left  = i
        if cpdf[i] >= 0.5000 and  i_mean==None:
            i_mean  = i
        if cpdf[i] >= 0.8414 and i_right==None:
            i_right = i
    z_left, z_mean, z_right = lag[i_left], lag[i_mean], lag[i_right]
    ymin, ymax = np.min(corr_normalized), 1.25*np.max(corr_normalized)
    ax3.vlines(x=z_left,  ymin=ymin, ymax=ymax, 
               ls=':',  zorder=8, color='black', linewidth=1)
    ax3.vlines(x=z_peak,  ymin=0.8*ymax, ymax=ymax, 
               ls='--', zorder=8, color='green', linewidth=1)
    ax3.vlines(x=z_mean,  ymin=ymin, ymax=ymax, 
               ls='--', zorder=8, color='black', linewidth=1)
    ax3.vlines(x=z_right, ymin=ymin, ymax=ymax, 
               ls=':',  zorder=8, color='black', linewidth=1)
    ax3.text(z_left, 0,  "{:.4f}".format(z_left), size=8, 
             horizontalalignment='center', verticalalignment='top', 
             color='black', zorder=7)
    ax3.text(z_peak, ymax,  "{:.4f}".format(z_peak), size=10, 
             horizontalalignment='center', verticalalignment='bottom', 
             color='green', zorder=7)
    ax3.text(z_mean, 0,  "{:.4f}".format(z_mean), size=10, 
             horizontalalignment='center', verticalalignment='bottom', 
             color='black', zorder=7)
    ax3.text(z_right, 0, "{:.4f}".format(z_right), size=8, 
             horizontalalignment='center', verticalalignment='top', 
             color='black', zorder=7)
    ax3.set_xlim(lag[i_mean]-3*(lag[i_mean] -lag[i_left]),
                 lag[i_mean]+3*(lag[i_right]-lag[i_mean]))
    ax3.set_ylim(ymin, ymax)
    ax3.set_ylim(ymin, ymax)
except TypeError:
    print('A zoom-in for the redshift peak failed: indices of mean is invalid.')
    pass
plt.savefig("1dspec-redshift-guess.jpg", dpi=200, bbox_inches='tight')
