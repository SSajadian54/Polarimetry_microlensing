# Polarimetry_microlensing

There are several codes here, which simulate polarimetry in microlensing events.  

Some of them (depend_hot.py, depend_cool.py, and depend_medium.py) make the normalized limb-darkening and stokes profiles for different types of stars 
(diffrent values of logg and Teff) and in the standard filters (UBVRI) according to simulations of these profiles versus wavelength made by P. Harrington 
(http://www.astro.umd.edu/~jph/KURUCZ_CONTINUUM/).  

The main code "BPM150597.cpp" calculate the polarimetry microlensing curves using the precalculated limb-darkening and stokes profiles of stars.  

In the case that you use these codes for your academic publications, please cite this paper:  

https://ui.adsabs.harvard.edu/abs/2019MNRAS.487..908S/abstract


