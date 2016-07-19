#!/usr/bin/env python

#Python Script for performing HST photometry, schechter parameterization and size estimation
#v1.00 
#Author: Joseph T. W. McGuire, University of Leicester, Department of Physics and Astronomy, XROA Group 
#Contact: jtwm1@leicester.ac.uk
#Last Modified: Jun 2016

#---------------------------------------------------------------------------------------------------------------------#
#Imports
import os, sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy.io import fits
from scipy.interpolate import interp1d
from astropy.cosmology import Planck15

#---------------------------------------------------------------------------------------------------------------------#
#User Input
parser = argparse.ArgumentParser(description="Performs HST photometry, shechter paramitization and size estimation, \
									using input fits files. Please read code before first use as some functions may \
									be deactivated for authors ease of use. Assumes HST image counts are in 1e5 es^-1 \
									and use of GAIA input files.")
parser.add_argument("rootdir", help="Directory path of the program, ensure directory string ends with /")
parser.add_argument("input", help="File containing initial GRB host information")
parser.add_argument("-ZPT","--mABZPT",type=float,nargs=1,default=[26.4524],help="AB magnitude zero point, \
					left blank will default to infinite aperture WFC3 corrected value.")
parser.add_argument("-pta","--pixtoarc",type=float,nargs=1,default=[0.07],help="HST image pixel to arcsec scale, \
					left blank will default to 0.07arcsec per pixel.")
parser.add_argument("-fil","--photplam",type=float,nargs=1,default=[13922.907],help="HST Filter wavelength, \
					left blank will default to JH_140 filter. In units of A.")
args = parser.parse_args() 

rootdir    = args.rootdir
file_input = args.input

#User defined parameters
mABZPT      = args.mABZPT[0]
pix_to_arc  = args.pixtoarc[0] 
PhotPlam    = args.photplam[0]

print
print "Loading user defined parameters..."
print

#---------------------------------------------------------------------------------------------------------------------#
#Conversions and constants
PhotFnu    = 9.5291135e-08                 #Jy s e^-1. HST conversion rate from fits image e s^-1 to flux density in Jy 
Jy_to_erg  = 1e-23                         #Converts Jy to erg s^-1 Hz^-1 cm^-2
pc_to_cm   = 3.08567758e18                 #Converts parsecs to cm
rad_to_arc = 206264.806                    #Converts radians to arcsec
K_Jy       = 3631                          #AB magnitude Jy constant
K_fwhm     = 2. * np.sqrt(2. * np.log(2))  #Full-width half-maxiumum constant from FWHM = K * sigma

#---------------------------------------------------------------------------------------------------------------------#
#Load in txt file containing basic GRB parameters
Dictionary = {}
Dictionary.setdefault("Name", [])
Dictionary.setdefault("z", [])
Dictionary.setdefault("dL", [])
Dictionary.setdefault("dA", [])
Dictionary.setdefault("Av", [])
Dictionary.setdefault("PSF", [])
Dictionary.setdefault("EEF", [])

print "Loading GRB parameters..."
print

#Cosmology from Plank 2015 Results 
with open(rootdir + file_input, "r") as fi:
	next(fi)
	for line in fi:
		Dictionary["Name"].append(line.split()[0])
		Dictionary["z"].append(float(line.split()[1]))
		Dictionary["dL"].append(float(line.split()[2]))
		Dictionary["dA"].append(float(line.split()[3]))
		Dictionary["Av"].append(float(line.split()[4]))
		Dictionary["PSF"].append(float(line.split()[5]))
		Dictionary["EEF"].append(float(line.split()[6]))

#Output
FITs_Aperture    = []
FITs_Target      = []
FITs_Mean        = []
FITs_Std         = []
FITs_C_Rate      = []
FITs_app_mag     = []
FITs_abs_mag     = []
FITs_flux_den    = []
FITs_Lum         = []
FITs_B_UV        = []
FITs_m1600       = []
FITs_M1600       = []
FITs_f1600       = []
FITs_L1600       = []
FITs_m1500       = []
FITs_M1500       = []
FITs_f1500       = []
FITs_L1500       = []
FITs_alpha       = []
FITs_phi         = []
FITs_M_sch       = []
FITs_Lz3         = []
FITs_Lz6         = []
FITs_R_half_aper = []
FITs_R_half      = []
FITs_sigma_mag   = []
FITs_R_e         = []
FITs_P_cc        = []
			
#---------------------------------------------------------------------------------------------------------------------#
#Photometry Equations
def Counts(Target, Mean, Std):                           #Calculates total count rate in e s^-1
	CV   = (Target - Mean) / 1e5
 	pStd = CV + (Std / 1e5)    
	mStd = CV - (Std / 1e5)    
	return CV, pStd, mStd

def apmag(CV, pStd, mStd, GE):                           #Calculates apparent magnitude from count rate
	apmag = -2.5 * np.log10(CV)   + mABZPT - GE
	perr  = -2.5 * np.log10(pStd) + mABZPT - GE
	merr  = -2.5 * np.log10(mStd) + mABZPT - GE
	return apmag, perr, merr, apmag-perr, merr-apmag

def abmag(apmag, pStd, mStd, dL, z):                     #Calculates absolute magnitude from apparent and cosmology
	abmag = apmag - 5 * np.log10(dL/10) + 2.5 * np.log10(1+z)
	perr  = pStd  - 5 * np.log10(dL/10) + 2.5 * np.log10(1+z)
	merr  = mStd  - 5 * np.log10(dL/10) + 2.5 * np.log10(1+z)
	return abmag, perr, merr, abmag-perr, merr-abmag	

def Flux(apmag, pStd, mStd):                             #Calculates flux density from apparent magnitude
	Flux = 10**(apmag / -2.5) * K_Jy * Jy_to_erg
	perr = 10**(pStd  / -2.5) * K_Jy * Jy_to_erg
	merr = 10**(mStd  / -2.5) * K_Jy * Jy_to_erg
	return Flux, perr, merr, perr-Flux, Flux-merr

def Lum(Flux, pStd, mStd, dL, z):                        #Calculates luminosity from flux density and cosmology
	Lum  = (4 * np.pi * Flux * dL**2) / (1+z)
	perr = (4 * np.pi * pStd * dL**2) / (1+z)
	merr = (4 * np.pi * mStd * dL**2) / (1+z)
	return Lum, perr, merr, perr-Lum, Lum-merr

#---------------------------------------------------------------------------------------------------------------------#
#Radius equations
def radius(aperture, pinterp, minterp, dA):              #Calculates galaxy radius in kpc from aperture and cosmology
	radius = (dA * (aperture / rad_to_arc)) / 1000
	perr   = (dA * (pinterp  / rad_to_arc)) / 1000
	merr   = (dA * (minterp  / rad_to_arc)) / 1000
	return radius, perr, merr, perr-radius, radius-merr

def PSFcorr(aperture, pinterp, minterp, fwhm, wingPSF, rused):	 #Point Spread Function correction for radius estimation
	arcsec = np.sqrt(aperture**2 - (fwhm / K_fwhm)**2) - (((rused  - 0.24)/0.36) * wingPSF)
	perr   = np.sqrt(pinterp**2  - (fwhm / K_fwhm)**2) - (((rused  - 0.24)/0.36) * wingPSF)
	merr   = np.sqrt(minterp**2  - (fwhm / K_fwhm)**2) - (((rused  - 0.24)/0.36) * wingPSF)
	return arcsec, perr, merr, perr-arcsec, arcsec-merr

#---------------------------------------------------------------------------------------------------------------------#
#Conversion Equations
def abconv(abmag, perr, merr, beta, berr, z, newl):      #Converts absolute magnitude from one wavelength to another
	abconv = abmag - 2.5 * (beta + 2) * np.log10((1+z) * newl / PhotPlam)
	err    = np.sqrt(berr**2 + (np.average([perr, merr]))**2)
	return abconv, err, abconv-err, abconv+err

def apconv(abconv, pStd, mStd, dL, z):                   #Finds converted apparent magnitude from converted absolute
	apconv = abconv + 5 * np.log10(dL/10) - 2.5 * np.log10(1+z)
	perr   = pStd   + 5 * np.log10(dL/10) - 2.5 * np.log10(1+z)
	merr   = mStd   + 5 * np.log10(dL/10) - 2.5 * np.log10(1+z)
	aerr   = np.average([apconv-perr, merr-apconv])
	return apconv, aerr, apconv-aerr, apconv+aerr

#---------------------------------------------------------------------------------------------------------------------#
#Beta_UV linear parameterization from Duncan & Conselice 2015
def beta(abmag):                                         #Calculates beta from absolute magnitude
	beta = (-2.05) + (-0.13) * (abmag + 19.5)
	err  = 0.25#np.sqrt(0.04**2 + 0.04**2)
	return beta, err, beta+err, beta-err

#---------------------------------------------------------------------------------------------------------------------#
#Schecter parameterization from Bouwens et al. 2015
def alpha(z):                                            #Calculates alpha from redshift
	a   = (-1.87) + (-0.10) * (z - 6)
	err = np.sqrt(0.05**2 + 0.03**2)
	return a, err, a+err, a+err

def phi(z):                                              #Calculates phi* from redshift
	phi = (0.46) * 10**((-0.27) * (z - 6))
	err = np.sqrt(0.105**2 + 0.05**2)
	return phi, err, phi+err, phi-err

def M_sch(z):                                            #Calculates M* from redshift
	M   = (-20.95) + (0.01) * (z - 6)
	err = np.sqrt(0.10**2 + 0.06**2)
	return M, err, M-err, M+err

#---------------------------------------------------------------------------------------------------------------------#
#Other
def find_nearest(array, value):                          #Finds closest value given two sets of arrays to compare
    i = (np.abs(array - value)).argmin()
    return array[i], i

def ratio(mag1, mag2):                                   #Finds the ratio between two magnitudes
	ratio = 100**((mag1 - mag2) / 5.)
	return ratio

#---------------------------------------------------------------------------------------------------------------------#
#Photometry extraction

#Photometry FITs files contain aperture ID, x,y coords (pix), signal (1e5 e s^-1) and aperture size (pix)
#Photometry FITs files assumed to be made from gaiatofits.py or follow same convention
#ID 1 corresponds to galaxy aperture, ID 2-41 correspond to 40 background sky apertures

#All values are estimated with 1sigma standard deviation

for i in range(len(Dictionary["Name"])):
	print "Loading in ", Dictionary.get("Name")[i], " FITs contents..."
	print

	if os.path.isfile(Dictionary.get("Name")[i] + "_apr27_CenCorr_" + ".fits"):
		FIT_file     = fits.open(rootdir + Dictionary.get("Name")[i] + "_apr27_CenCorr_" + ".fits")
		FIT_contents = FIT_file[1].data

		#Sets up aperture list and arrays
		Apertures = [0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 
					 	0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 
					 	0.5, 0.525, 0.55, 0.575, 0.6]

		Target    = np.zeros((len(Apertures)))
		Mean      = np.zeros((len(Apertures)))
		Std       = np.zeros((len(Apertures)))

		C_Rate    = np.zeros((len(Apertures),3)) 

		#-----------------------------------------------------------------------------------------------------#	
		#Finds count rate, mean and standard diviation for each aperture
		for j in range(len(Apertures)):
			mask            = FIT_contents["APERTURE"] == np.around((Apertures[j] / pix_to_arc), 1)
			Contents_masked = FIT_contents[mask]

			Target[j] = Contents_masked["SIG"][0]  
			Mean[j]   = Contents_masked["SIG"][1:41].mean()
			Std[j]    = Contents_masked["SIG"][1:41].std()
	
		CV, pStd, mStd = Counts(Target, Mean, Std)

		C_Rate[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
				[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]] = CV
		C_Rate[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
				[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]] = pStd
		C_Rate[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
				[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]] = mStd

		#Finds the max count rate from apertures
		print "Finding max count rate..."
		print
		for k in range(len(Apertures)-1):
			if C_Rate[k][0] > C_Rate[k-1][0] and C_Rate[k][0] > C_Rate[k+1][0]:
				max_value = k

		#-----------------------------------------------------------------------------------------------------#	
		#Apparent and absolute magnitudes
		print "Calculating apparent and absolute magnitudes..."
		print
		#Using a 3 point averaged, EEF and AG corrected C_Rate
		C_Rate_avg      = np.average([C_Rate[max_value][0], C_Rate[max_value-1][0], C_Rate[max_value+1][0]])
		C_Rate_avg_perr = np.average([C_Rate[max_value][1], C_Rate[max_value-1][1], C_Rate[max_value+1][1]])
		C_Rate_avg_merr = np.average([C_Rate[max_value][2], C_Rate[max_value-1][2], C_Rate[max_value+1][2]])

		print C_Rate_avg

		C_Rate_EEF      = C_Rate_avg      / Dictionary.get("EEF")[i]
		C_Rate_EEF_perr = C_Rate_avg_perr / Dictionary.get("EEF")[i]
		C_Rate_EEF_merr = C_Rate_avg_merr / Dictionary.get("EEF")[i]

		print C_Rate_EEF

		app_mag = apmag(C_Rate_EEF, C_Rate_EEF_perr, C_Rate_EEF_merr, Dictionary.get("Av")[i])
		abs_mag = abmag(app_mag[0], app_mag[1], app_mag[2], Dictionary.get("dL")[i], Dictionary.get("z")[i])

		#Flux Density and Luminosity
		print "Calculating flux density and luminosity..."
		print
		flux_den = Flux(app_mag[0], app_mag[1], app_mag[2])
		Lumin    = Lum(flux_den[0], flux_den[1], flux_den[2], Dictionary.get("dL")[i] * pc_to_cm, 
						Dictionary.get("z")[i])

		#-----------------------------------------------------------------------------------------------------#	
		#Converts to standardized monochromatic wavelengths
		print "Converting to standardized monochromatic wavelengths..."
		print
		b_UV  = beta(abs_mag[0])

		M1600 = abconv(abs_mag[0], abs_mag[3], abs_mag[4], b_UV[0], b_UV[1], Dictionary.get("z")[i], 1600) 
		m1600 = apconv(M1600[0], M1600[2], M1600[3], Dictionary.get("dL")[i], Dictionary.get("z")[i])
		f1600 = Flux(m1600[0], m1600[2], m1600[3])
		L1600 = Lum(f1600[0], f1600[1], f1600[2], Dictionary.get("dL")[i] * pc_to_cm, Dictionary.get("z")[i])

		M1500 = abconv(abs_mag[0], abs_mag[3], abs_mag[4], b_UV[0], b_UV[1], Dictionary.get("z")[i], 1500) 
		m1500 = apconv(M1500[0], M1500[2], M1500[3], Dictionary.get("dL")[i], Dictionary.get("z")[i])
		f1500 = Flux(m1500[0], m1500[2], m1500[3])
		L1500 = Lum(f1500[0], f1500[1], f1500[2], Dictionary.get("dL")[i] * pc_to_cm, Dictionary.get("z")[i])

		#-----------------------------------------------------------------------------------------------------#	
		#Schechter parameters, parameterization from Bouwens et al 2015
		print "Calculating schecter parameters..."
		print
		a     = alpha(Dictionary.get("z")[i])
		p     = phi(Dictionary.get("z")[i])
		M     = M_sch(Dictionary.get("z")[i])

		#-----------------------------------------------------------------------------------------------------#	
		#L* Ratio
		print "Calculating L* ratios..."
		print		
		Lz3 = ratio(M_sch(3.0)[0], M1600[0])
		Lz6 = ratio(M[0], M1600[0])  

		#-----------------------------------------------------------------------------------------------------#	
		#Half-light radius estimation

		#Uses a cublic spline interpolation function to fill in values between apertures
		#Half-light radius is found when comparing the half-max count value to interpolated values 
		print "Calculating half-light radius and interpolation..."
		print
		interp_Aper       = np.linspace(Apertures[0], Apertures[-1], num=81, endpoint=True)

		interp_Counts     = interp1d(Apertures[0::2], C_Rate[0::2,0], kind='cubic')
		interp_C_perr     = interp1d(Apertures[0::2], C_Rate[0::2,1], kind='cubic')
		interp_C_merr     = interp1d(Apertures[0::2], C_Rate[0::2,2], kind='cubic')

		for l in range(len(interp_Aper)):
			if interp_Aper[l] == Apertures[max_value]:
				peak_aper = l

		interp_Count_mask = []
		interp_Cperr_mask = []
		interp_Cmerr_mask = []

		for m in range(len(interp_Aper)):
			if m < peak_aper:
				interp_Count_mask += [interp_Counts(interp_Aper)[m]]
				interp_Cperr_mask += [interp_C_perr(interp_Aper)[m]]
				interp_Cmerr_mask += [interp_C_merr(interp_Aper)[m]]

		half_max_C        = find_nearest(interp_Count_mask, C_Rate_EEF / 2.)
		half_max_C_perr   = find_nearest(interp_Cmerr_mask, C_Rate_EEF / 2.)
		half_max_C_merr   = find_nearest(interp_Cperr_mask, C_Rate_EEF / 2.)

		print half_max_C
		print half_max_C_perr
		print half_max_C_merr

		#Curtis-Lake wing PSF correction
		wingPSF = [0.025, 0.045, 0.047] 

		half_max_PSF      = PSFcorr(interp_Aper[half_max_C[1]], interp_Aper[half_max_C_perr[1]],
									interp_Aper[half_max_C_merr[1]], Dictionary.get("PSF")[i], 
									wingPSF[i], Apertures[max_value])

		half_max_R        = radius(half_max_PSF[0], half_max_PSF[1], half_max_PSF[2], Dictionary.get("dA")[i])

		#-----------------------------------------------------------------------------------------------------#
		#Astrometry offset
		offset = [0.06,0.13,0.21]
		offerr = [0.02,0.04,0.07]

		offkpc = radius(offset[i], offset[i] + offerr[i], offset[i] - offerr[i], Dictionary.get("dA")[i])

		#-----------------------------------------------------------------------------------------------------#
		#Probability of chance coincidence using H-count values from Metaclfe et al 2006
		print "Calculating probability of chance coincidence..."
		print
		H_mag    = [17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 
					22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0] 
		H_counts = [1349.0, 1199.0, 1798.0, 2398.0, 3672.0, 7868.0, 8167.0, 12663.0, 16934.0, 23528.0, 
					31396.0, 71800.0, 128000.0, 251000.0, 311000.0, 551000.0, 503000.0]

		print app_mag[0]

		for n in range(len(H_mag)): 
			if H_mag[n-1] < app_mag[0] < H_mag[n]:
				sigma_mag = (sum(H_counts[:n-1]) + ((app_mag[0] - H_mag[n-1]) * H_counts[n-1])) / (3600**2)

			elif H_mag[-1] < app_mag[0] < 29.0:
				sigma_mag = (sum(H_counts[:-1]) + ((app_mag[0] - H_mag[-1]) * H_counts[-1])) / (3600**2)

		print sigma_mag
	
		print half_max_PSF[0]

		R_e  = 2 * half_max_PSF[0]	
		P_cc = 1 - np.exp(-np.pi * R_e**2 * sigma_mag)

		print P_cc

		#-----------------------------------------------------------------------------------------------------#
		#Stores each loop into FITs input variable
		print "Storing calculated parameters..."
		print
		FITs_Aperture    += [Apertures[max_value]]
		FITs_Target      += [Target[max_value] /1e5]
		FITs_Mean        += [Mean[max_value] /1e5]
		FITs_Std         += [Std[max_value] /1e5]
		FITs_C_Rate      += [C_Rate_EEF, C_Rate_EEF_perr, C_Rate_EEF_merr] 
		FITs_app_mag     += [app_mag[0], app_mag[3], app_mag[4]]    
		FITs_abs_mag     += [abs_mag[0], abs_mag[3], abs_mag[4]]    
		FITs_flux_den    += [flux_den[0], flux_den[3], flux_den[4]] 
		FITs_Lum         += [Lumin[0], Lumin[3], Lumin[4]]                         
		FITs_B_UV        += [b_UV[0], b_UV[1]]             
		FITs_m1600       += [m1600[0], m1600[1]]          
		FITs_M1600       += [M1600[0], M1600[1]]          
		FITs_f1600       += [f1600[0], f1600[3], f1600[4]]          
		FITs_L1600       += [L1600[0], L1600[3], L1600[4]]                         
		FITs_m1500       += [m1500[0], m1500[1]]          
		FITs_M1500       += [M1500[0], M1500[1]]          
		FITs_f1500       += [f1500[0], f1500[3], f1500[4]]          
		FITs_L1500       += [L1500[0], L1500[3], L1500[4]]                 
		FITs_alpha       += [a[0], a[1]]                      
		FITs_phi         += [p[0], p[1]]                      
		FITs_M_sch       += [M[0], M[1]]    
		FITs_Lz3         += [Lz3]
		FITs_Lz6         += [Lz6]                  
		FITs_R_half_aper += [half_max_PSF[0], half_max_PSF[3], half_max_PSF[4]] 
		FITs_R_half      += [half_max_R[0], half_max_R[3], half_max_R[4]]       
		FITs_sigma_mag   += [sigma_mag]
		FITs_R_e         += [R_e]
		FITs_P_cc        += [P_cc]
			#-----------------------------------------------------------------------------------------------------#
		#Close FIT file
		FIT_file.close()

	else:
		print "Photometry FITs file was not found, skipping."
		print 

#---------------------------------------------------------------------------------------------------------------------#
#Creates and fills FITs Table for Photometry results
print "Creating FITs Table..."
print
tbhdu = fits.BinTableHDU.from_columns(
			[fits.Column(name='Name',                                 format='7A', array=Dictionary["Name"]),
			 fits.Column(name='z',                                    format='E',  array=Dictionary["z"]),
			 fits.Column(name='Aperture',     unit="arcsec",          format='E',  array=FITs_Aperture),
			 fits.Column(name='Target',       unit="e/s",             format='E',  array=FITs_Target),
			 fits.Column(name='Mean',         unit="e/s",             format='E',  array=FITs_Mean),
			 fits.Column(name='Std',          unit="e/s",             format='E',  array=FITs_Std),
			 fits.Column(name='C_Rate',       unit="e/s",             format='E',  array=FITs_C_Rate[0::3]),
			 fits.Column(name='C_pStd',       unit="e/s",             format='E',  array=FITs_C_Rate[1::3]),
			 fits.Column(name='C_mStd',       unit="e/s",             format='E',  array=FITs_C_Rate[2::3]),
			 fits.Column(name='J(AB)',        unit="mag",             format='E',  array=FITs_app_mag[0::3]),
			 fits.Column(name='M(AB)',        unit="mag",             format='E',  array=FITs_abs_mag[0::3]),
			 fits.Column(name='mag(AB)_perr', unit="mag",             format='E',  array=FITs_app_mag[1::3]),
			 fits.Column(name='mag(AB)_merr', unit="mag",             format='E',  array=FITs_abs_mag[2::3]),
			 fits.Column(name='Flux_Density', unit="ergs/s/Hz/cm2",   format='E',  array=FITs_flux_den[0::3]),
			 fits.Column(name='F_Den_err',    unit="ergs/s/Hz/cm2",   format='E',  array=FITs_flux_den[1::3]),
			 fits.Column(name='Lum',          unit="ergs/s/Hz",       format='E',  array=FITs_Lum[0::3]),
			 fits.Column(name='L_err',        unit="ergs/s/Hz",       format='E',  array=FITs_Lum[1::3]),
			 fits.Column(name='Beta',                                 format='E',  array=FITs_B_UV[0::2]),
			 fits.Column(name='Beta_err',                             format='E',  array=FITs_B_UV[1::2]),
			 fits.Column(name='J1600',        unit="mag",             format='E',  array=FITs_m1600[0::2]),
			 fits.Column(name='M1600',        unit="mag",             format='E',  array=FITs_M1600[0::2]),
			 fits.Column(name='mag(1600)_err',unit="mag",             format='E',  array=FITs_M1600[1::2]),
			 fits.Column(name='f1600',        unit="ergs/s/Hz/cm2",   format='E',  array=FITs_f1600[0::3]),
			 fits.Column(name='f1600_perr',   unit="ergs/s/Hz/cm2",   format='E',  array=FITs_f1600[1::3]),
			 fits.Column(name='f1600_merr',   unit="ergs/s/Hz/cm2",   format='E',  array=FITs_f1600[2::3]),
			 fits.Column(name='L1600',        unit="ergs/s/Hz",       format='E',  array=FITs_L1600[0::3]),
			 fits.Column(name='L1600_perr',   unit="ergs/s/Hz",       format='E',  array=FITs_L1600[1::3]),
			 fits.Column(name='L1600_merr',   unit="ergs/s/Hz",       format='E',  array=FITs_L1600[2::3]),
			 fits.Column(name='J1500',        unit="mag",             format='E',  array=FITs_m1500[0::2]),
			 fits.Column(name='M1500',        unit="mag",             format='E',  array=FITs_M1500[0::2]),
			 fits.Column(name='mag(1500)_err',unit="mag",             format='E',  array=FITs_m1500[1::2]),
			 fits.Column(name='f1500',        unit="ergs/s/Hz/cm2",   format='E',  array=FITs_f1500[0::3]),
			 fits.Column(name='f1500_perr',   unit="ergs/s/Hz/cm2",   format='E',  array=FITs_f1500[1::3]),
			 fits.Column(name='f1500_merr',   unit="ergs/s/Hz/cm2",   format='E',  array=FITs_f1500[2::3]),
			 fits.Column(name='L1500',        unit="ergs/s/Hz",       format='E',  array=FITs_L1500[0::3]),
			 fits.Column(name='L1500_perr',   unit="ergs/s/Hz",       format='E',  array=FITs_L1500[1::3]),
			 fits.Column(name='L1500_merr',   unit="ergs/s/Hz",       format='E',  array=FITs_L1500[2::3]),
			 fits.Column(name='alpha',                                format='E',  array=FITs_alpha[0::2]),
			 fits.Column(name='alpha_err',                            format='E',  array=FITs_alpha[1::2]),
			 fits.Column(name='phi',          unit="x10^-3 Mpc^-3",   format='E',  array=FITs_phi[0::2]),
			 fits.Column(name='phi_err',      unit="x10^-3 Mpc^-3",   format='E',  array=FITs_phi[1::2]),
			 fits.Column(name='M*UV',         unit="mag",             format='E',  array=FITs_M_sch[0::2]),
			 fits.Column(name='M*UV_err',     unit="mag",             format='E',  array=FITs_M_sch[1::2]),
			 fits.Column(name='L/L*(z=3)',                            format='E',  array=FITs_Lz3),
			 fits.Column(name='L/L*(z=6)',                            format='E',  array=FITs_Lz6),
			 fits.Column(name='R_half_aper',  unit="arcsec",          format='E',  array=FITs_R_half_aper[0::3]),
			 fits.Column(name='R_aper_perr',  unit="arcsec",          format='E',  array=FITs_R_half_aper[1::3]),
			 fits.Column(name='R_aper_merr',  unit="arcsec",          format='E',  array=FITs_R_half_aper[2::3]),
			 fits.Column(name='R_half_kpc',   unit="kpc",             format='E',  array=FITs_R_half[0::3]),
			 fits.Column(name='R_kpc_perr',   unit="kpc",             format='E',  array=FITs_R_half[1::3]),
			 fits.Column(name='R_kpc_merr',   unit="kpc",             format='E',  array=FITs_R_half[2::3]),
			 fits.Column(name='sigma_mag',                            format='E',  array=FITs_sigma_mag),
			 fits.Column(name='R_e',          unit="arcsec",          format='E',  array=FITs_R_e),
			 fits.Column(name='P_cc',                                 format='E',  array=FITs_P_cc)
			])
tbhdu.writeto(rootdir + "Host_photometry" + "_apr27_CenCorr_AGrm_" + ".fits", clobber=True)

#---------------------------------------------------------------------------------------------------------------------#
#Produces paper graphs

#Luminosity Function
print "Plotting UV LF Graph..."
print
#Values are taken from Table 5/6 in Bouwens et al 2015
M_SW   = [-22.52,-22.02,-21.52,-21.02,-20.52,-20.02,-19.52,-18.77,-17.77,-16.77]
phi_SW = [2e-6,15e-6,53e-6,176e-6,320e-6,698e-6,1246e-6,1900e-6,6680e-6,13640e-6]
phi_pm = [2e-6,6e-6,12e-6,25e-6,41e-6,83e-6,137e-6,320e-6,1380e-6,4200e-6]

phi_LF = 0.39e-3
Mlf_LF = -21.1
a_LF   = -1.9
M_LF  = np.linspace(-23., -16., num=161, endpoint=True)

#LF from Section 4.2 Bouwens et al 2015
LF = (phi_LF * (np.log(10)/2.5)) * 10**(-0.4 * (M_LF - Mlf_LF) * (a_LF + 1)) * np.exp(-10**(-0.4 * (M_LF - Mlf_LF)))

fig = plt.figure(1)
fig, ax = plt.subplots(1, figsize=(5.65685,4))

ax.errorbar(M_SW, phi_SW, yerr=phi_pm, ecolor='c', color='c', capsize=0, marker='.', ls="None")
ax.plot(M_LF, LF, 'c-')
ax.plot(-20.95, 2e-7, mec='c', mfc='c', marker=r'$\uparrow$', markersize=30)

ax.axvline(FITs_M1600[0], c=u"w", ls="--")
ax.axvline(FITs_M1600[2], c=u"w", ls="--")
ax.axvline(FITs_M1600[4], c=u"w", ls="--")

ax.axvline(FITs_M1600[2] + FITs_M1600[3], c=u"0.25", ls="-")

ax.axvspan(FITs_M1600[0] + FITs_M1600[1], FITs_M1600[0] - FITs_M1600[1], facecolor='k', alpha=0.1, lw=0)
ax.axvspan(FITs_M1600[2] + FITs_M1600[3], FITs_M1600[2] - FITs_M1600[3], facecolor='0.25', alpha=0.1, lw=0)
ax.axvspan(FITs_M1600[4] + FITs_M1600[5], FITs_M1600[4] - FITs_M1600[5], facecolor='0.50', alpha=0.1, lw=0)

ax.text(-17, 1e-6, "z~6", color='c', fontsize=12)
ax.text(FITs_M1600[0] - 0.23, 1e-5, "GRB 130606A", color='w', fontsize=12, rotation="vertical")
ax.text(FITs_M1600[2] - 0.23, 1e-5, "GRB 050904", color='w', fontsize=12, rotation="vertical")
ax.text(FITs_M1600[4] - 0.23, 1e-5, "GRB 140515A", color='w', fontsize=12, rotation="vertical")

plt.yscale("log")

ax.set_xlabel(r"$M_{1600}$", fontsize=20)
ax.set_ylabel(r"$\log(\rm{Number})\,\rm{mag}^{-1}\,\rm{Mpc}^{-3}$", fontsize=20)

plt.savefig(rootdir  + "Host_LF" + "_apr27_CenCorr_AGrm_" + ".eps", bbox_inches="tight")
plt.close(fig)

#GRB Host comparison
print "Plotting GRB Host Graph..."
print
GName  = []
zhosts = []
MUV    = []
Merr   = []
Paper  = []

#Values are taken from Table 2 in Tough VII 2015 and Table 2 in Griener 2015
with open(rootdir + "grb_hosts.txt", "r") as fs:
	next(fs)
	for line in fs:
		GName  += [line.split()[0]]
		zhosts += [float(line.split()[1])]
		MUV    += [float(line.split()[2])]
		Merr   += [float(line.split()[3])]
		Paper  += [line.split()[4]]

#Values are taken from Table 6 in Bouwens et al 2015
Mz_B   = [3.8,4.9,5.9,6.8,7.9]
M_B    = [-20.88,-21.17,-20.94,-20.77,-20.21]
Me_B   = [0.08,0.12,0.20,0.28,0.33]

MperrB = [-20.96,-21.29,-21.14,-21.05,-20.54]
MmerrB = [-20.80,-21.05,-20.74,-20.49,-19.88] 

#Values are taken from Table 1 in Oesch et al 2010
Mz_O   = [0.75,1.25,1.75]
M_O    = [-19.17,-20.08,-20.17]
Me_O   = [0.51,0.36,0.34]

Mz_O2  = [1.5,1.9,2.5]
MperrO2= [-20.33,-20.68,-20.86]
MmerrO2= [-19.31,-19.74,-20.52]

MperrO = [-19.68,-20.44,-20.51]
MmerrO = [-18.66,-19.72,-19.83]

#Values are taken from Table 3 in Reddy & Steidel et al 2015
Mz_R   = [2.3,3.05]
M_R    = [-20.70,-20.97]
Me_R   = [0.11,0.14]

MperrR = [-20.81,-21.11]
MmerrR = [-20.59,-20.83]

fig = plt.figure(1)
fig, ax = plt.subplots(1, figsize=(5.65685,4))

#Paper results
ax.errorbar(Dictionary.get("z")[0], FITs_M1600[0], yerr=[[FITs_M1600[1], FITs_M1600[1]]],
 	ecolor="r",marker="s",ls="None",capsize=0,mec="r",mfc="r",markersize=5)
ax.errorbar(Dictionary.get("z")[1], FITs_M1600[2], yerr=[[FITs_M1600[3], FITs_M1600[3]]],
 	ecolor="r",marker="s",ls="None",capsize=0,mec="r",mfc="r",markersize=5)
ax.errorbar(Dictionary.get("z")[2], FITs_M1600[4], yerr=[[FITs_M1600[5], FITs_M1600[5]]],
 	ecolor="r",marker="s",ls="None",capsize=0,mec="r",mfc="r",markersize=5)

#Other samples
for i in range(len(GName)):
	if Paper[i] == "T":
		if Merr[i] == 0.0:
			ax.errorbar(zhosts[i], MUV[i], yerr=Merr[i], marker="v", ls="None", capsize=0, mec="b", mfc="w", markersize=4)
		else:
			ax.errorbar(zhosts[i], MUV[i], yerr=Merr[i], ecolor='b', marker=".", ls="None", capsize=0, mec="b", mfc="b", markersize=6)
	elif Paper[i] == "G":
		if Merr[i] == 0.0:
			ax.errorbar(zhosts[i], MUV[i], yerr=Merr[i], marker="v", ls="None", capsize=0, mec="g", mfc="w", markersize=4)
		else:
			ax.errorbar(zhosts[i], MUV[i], yerr=Merr[i], ecolor="g", marker="D", ls="None", capsize=0, mec="g", mfc="g", markersize=3)
	elif Paper[i] == "Ta":
		if Merr[i] == 0.0:
			ax.errorbar(zhosts[i], MUV[i], yerr=Merr[i], marker="v", ls="None", capsize=0, mec="m", mfc="w", markersize=4)
		else:
			ax.errorbar(zhosts[i], MUV[i], yerr=Merr[i], ecolor='m', marker=".", ls="None", capsize=0, mec="m", mfc="m", markersize=6)

#M*
#z~1-2
ax.fill_between(Mz_O, MperrO, MmerrO, color='0.5', alpha=.25)
#z~2-3
ax.fill_between(Mz_R, MperrR, MmerrR, color='0.5', alpha=.25)
#z~3+
ax.fill_between(Mz_B, MperrB, MmerrB, color='0.5', alpha=.25)

ax.fill_between([0.0565,Mz_O[0]], [-17.061,MperrO[0]], [-17.001,MmerrO[0]], color='0.5', alpha=.25)
ax.fill_between([Mz_O[-1],Mz_R[0]], [MperrO[-1],MperrR[0]], [MmerrO[-1],MmerrR[0]], color='0.5', alpha=.25)
ax.fill_between([Mz_R[-1],Mz_B[0]], [MperrR[-1],MperrB[0]], [MmerrR[-1],MmerrB[0]], color='0.5', alpha=.25)

ax.text(4.5,-15.8, r"$M^{*}_{UV}$ - Various works", color='0.35', fontsize=11)
ax.text(4.5,-15.2, "TOUGH VII", color='b', fontsize=11)
ax.text(4.5,-14.6, "Greiner+ 2015*", color='g', fontsize=11)
ax.text(4.5,-14.0, "Tanvir+ 2012", color='m', fontsize=11)
ax.text(4.5,-13.4, "This work", color='r', fontsize=11)

ax.invert_yaxis()

ax.set_xlim(0, 7)

ax.set_xlabel(r"$z$", fontsize=20)
ax.set_ylabel(r"$M_{1600}$", fontsize=20)

plt.savefig(rootdir  + "Host_MZ" + "_apr27_CenCorr_AGrm_" + ".eps", bbox_inches="tight")
plt.close(fig)

#Half-light galaxy radius
print "Plotting Half-light radius Graph..."
print
ID  = []
r50 = []
Muv = []
z   = []
dis = []

#Data provided by Dr. Curtis-Lake from her submitted 2014 paper
with open(rootdir + "size_MUV_z6.txt", "r") as fs:
	next(fs)
	for line in fs:
		ID  += [int(line.split()[0])]
		r50 += [float(line.split()[1])]
		Muv += [float(line.split()[2])]
		z   += [float(line.split()[3])]
		dis += [int(line.split()[4])]

fig = plt.figure(1)
fig, ax = plt.subplots(1, figsize=(5.65685,4))

for i in range(len(ID)):
	if dis[i] == 0:
		ax.plot(Muv[i], r50[i], 'k.')
	else:
		ax.plot(Muv[i], r50[i], 'k.')

ax.errorbar(FITs_M1500[0], FITs_R_half[0],
	xerr=[[FITs_M1500[1], FITs_M1500[1]]],
	yerr=[[FITs_R_half[1], FITs_R_half[2]]],
	ecolor='r', mec='r', mfc='r', capsize=0, marker='s', markersize=5, ls="None")
ax.errorbar(FITs_M1500[2], FITs_R_half[3],
	xerr=[[FITs_M1500[3], FITs_M1500[3]]],
	yerr=[[FITs_R_half[4], FITs_R_half[5]]],
	ecolor='r', mec='r', mfc='r', capsize=0, marker='s', markersize=5, ls="None")
ax.errorbar(FITs_M1500[4], FITs_R_half[6],
	xerr=[[FITs_M1500[5], FITs_M1500[5]]],
	yerr=[[FITs_R_half[7], FITs_R_half[8]]],
	ecolor='r', mec='r', mfc='r', capsize=0, marker='s', markersize=5, ls="None")

ax.text(-22.5,0.0325,"Curtis-Lake+ 2015",color='k',fontsize=12)
ax.text(-22.5,0.02,"This work",color='r',fontsize=12)
ax.text(-19.0,0.02,"z~6",color='k',fontsize=12)

ax.set_yscale("log", nonposy="clip")
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

ax.set_xlabel(r"$M_{1500}$", fontsize=20)
ax.set_ylabel(r"$R_{\rm{half}}\,(\rm{kpc})$", fontsize=20)

plt.savefig(rootdir  + "Host_size" + "_apr27_CenCorr_AGrm_" + ".eps", bbox_inches="tight")
plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#End of program
print "Program executed correctly, ending process."
print
