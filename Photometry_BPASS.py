#!/usr/bin/env python

#Python Script for BPASS graphing and fitting
#v1.00 
#Author: Joseph T. W. McGuire, University of Leicester, Department of Physics and Astronomy, XROA Group 
#Contact: jtwm1@leicester.ac.uk
#Last Modified: June 2016

#---------------------------------------------------------------------------------------------------------------------#
#Imports
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import urllib2

#---------------------------------------------------------------------------------------------------------------------#
#File paths
rootdir         = "/home/jtwm1/HSTResearch/BPASS/Eliz/"

#---------------------------------------------------------------------------------------------------------------------#
#Conversion factors
lum_sol_to_erg_s = 3.846e33
SOL              = 299792458
pc_to_cm         = 3.08567758e18 
Jy_to_erg        = 1e-23
K_fwhm           = 2 * np.sqrt(2 * np.log(2))

#---------------------------------------------------------------------------------------------------------------------#
#Loads BPASS SED data"
print "Loading BPASS SED data..."
print
#100Myr continuous SF
Files = ["out_z001_bin_6.0_v1.cont", "out_z001_bin_6.1_v1.cont", "out_z001_bin_6.2_v1.cont", 
			"out_z001_bin_6.3_v1.cont", "out_z001_bin_6.4_v1.cont", "out_z001_bin_6.5_v1.cont",
			"out_z001_bin_6.6_v1.cont", "out_z001_bin_6.7_v1.cont", "out_z001_bin_6.8_v1.cont",
			"out_z001_bin_6.9_v1.cont", "out_z001_bin_7.0_v1.cont", "out_z001_bin_7.1_v1.cont",
			"out_z001_bin_7.2_v1.cont", "out_z001_bin_7.3_v1.cont", "out_z001_bin_7.4_v1.cont",
			"out_z001_bin_7.5_v1.cont", "out_z001_bin_7.6_v1.cont", "out_z001_bin_7.7_v1.cont",
			"out_z001_bin_7.8_v1.cont", "out_z001_bin_7.9_v1.cont", "out_z001_bin_8.0_v1.cont"]

#Loads in wavelength to one list, same for all files
with open(rootdir + Files[0], "r") as fi:
	next(fi)
	Lambda_A = [float(line.split()[0]) for line in fi] #In Angstrom

#Initializes spectral flux density list to be array[x] where x in number of files
LamSFD = [0] * len(Files) 

#Loads in spectral flux density of all files into array, which now has shape array[x][y] where y is number of values
for i in range(len(Files)):
	with open(rootdir + Files[i], "r") as fi:
		next(fi)
		LamSFD[i] = [float(line.split()[6]) for line in fi] #In A * ergs/s/A/cm2 at 10pc

#Expands y dimension to be same length for all x
for i in range(len(LamSFD)):
	LamSFD[i] += [0] * (len(Lambda_A) - len(LamSFD[i]))

#---------------------------------------------------------------------------------------------------------------------#
#Conversions
print "Conversions..."
print
#cm2 at 10pc
cm2 = 4 * np.pi * (10 * pc_to_cm)**2

#Hz at each wavelength
Hz = [((SOL * 1e10) / x) for x in Lambda_A]

#Converts LamSFD A * ergs/s/A/cm2 to SPD ergs/s/Hz
SPD = [0] * len(Files) 

for i in range(len(SPD)):
	SPD[i] = [((x * cm2) / y) for x, y in zip(LamSFD[i], Hz)]

#For the continuous star formation, each age bin must be weighted by
#(10**(age + 0.05) - 10**(age - 0.05)) / 1e6
Age = [6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 
			7.6, 7.7, 7.8, 7.9, 8.0]

Weighting    = [((10**(x + 0.05) - 10**(x - 0.05)) / 1e6) for x in Age]
Weighting[0] = 10**6.05 / 1e6 #Starting age has a different weighting

Weighted_SPD = [0] * len(Files)

for i in range(len(Weighted_SPD)):
	Weighted_SPD[i] = [(x * Weighting[i]) for x in SPD[i]]

#Sum Age bins up together for continuous
SPD_Con = [(a + b + c + d + e + f + g + h + i + j + k + l + m + n + o + p + q + r + s + t + u) 
			for a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u in
				zip(*Weighted_SPD)] 

#Instantaneous star formation at 10 Myr
SPD_Ins = [x for x in SPD[10]]

#Determines wavelength range of SED, cutting out values outside this range
SED_Range = [1000 <= x <= 30000 for x in Lambda_A]
SED_index = [x for x in range(len(SED_Range)) if SED_Range[x] == True]

Lam_range     = [Lambda_A[x] for x in SED_index]
SPD_Con_range = [SPD_Con[x]  for x in SED_index]
SPD_Ins_range = [SPD_Ins[x]  for x in SED_index]

#---------------------------------------------------------------------------------------------------------------------#
#Resolution scaling
print "Scaling to R~1000..."
print

#Finds the mean SPD within the range (continuum)
SPD_Con_mean = [np.mean([i, j, k, l, m]) for i, j, k, l, m in 
				zip(SPD_Con_range[0::5], SPD_Con_range[1::5], SPD_Con_range[2::5], 
					SPD_Con_range[3::5], SPD_Con_range[4::5])]
SPD_Ins_mean = [np.mean([i, j, k, l, m]) for i, j, k, l, m in 
				zip(SPD_Ins_range[0::5], SPD_Ins_range[1::5], SPD_Ins_range[2::5], 
					SPD_Ins_range[3::5], SPD_Ins_range[4::5])]

#Threshold for significant emission lines
Con_sig_K = 2.5
Ins_sig_K = 6.2

#Finds significant emission lines 
SPD_Con_sig = Con_sig_K * np.std(SPD_Con_mean) 
SPD_Ins_sig = Ins_sig_K * np.std(SPD_Ins_mean) 

Con_em   = [x >= SPD_Con_sig for x in SPD_Con_range]
Ins_em   = [x >= SPD_Ins_sig for x in SPD_Ins_range]
Con_em_i = [x for x in range(len(Con_em)) if Con_em[x] == True]
Ins_em_i = [x for x in range(len(Ins_em)) if Ins_em[x] == True]

#Finds emission lines within each range for wavelength and SPD
Lam_Con_em = [Lam_range[x] for x in Con_em_i]
Lam_Ins_em = [Lam_range[x] for x in Ins_em_i]
SPD_Con_em = [SPD_Con_range[x] for x in Con_em_i]
SPD_Ins_em = [SPD_Ins_range[x] for x in Ins_em_i]

#Significance detection testing
"""
#Plots emission lines
f, (ax1, ax2) = plt.subplots(2, figsize=(7,8))

ax1.set_xscale("log", nonposy="clip")
ax1.set_yscale("log", nonposy="clip")
ax2.set_xscale("log", nonposy="clip")
ax2.set_yscale("log", nonposy="clip")

ax1.plot(Lam_Con_em, SPD_Con_em, 'm.')
ax1.plot(Lam_range, SPD_Con_range, 'k-')
ax2.plot(Lam_Ins_em, SPD_Ins_em, 'm.')
ax2.plot(Lam_range, SPD_Ins_range, 'k-')

ax1.text(2000, 2e27, str(Con_sig_K) + "$\sigma$", color="k", fontsize=12)
ax2.text(2000, 2e25, str(Ins_sig_K) + "$\sigma$", color="k", fontsize=12)

ax1.set_xlim(1000, 30000)
ax1.set_ylim(1e27, 1e30)
ax2.set_xlim(1000, 30000)
ax2.set_ylim(1e25, 1e28)

ax2.set_xlabel(r"$\lambda\,(\rm{\AA})$", fontsize=20)
ax1.set_ylabel(r"$\rm{Power\,Density}\,(\rm{ergs\,s^{-1}\,Hz^{-1}})$", fontsize=20)
ax2.set_ylabel(r"$\rm{Power\,Density}\,(\rm{ergs\,s^{-1}\,Hz^{-1}})$", fontsize=20)

f.subplots_adjust(hspace=0.1)

plt.savefig(rootdir  + "emLines" + ".eps", bbox_inches="tight")
plt.close(f)
"""

#Finds emission line width resolution
Con_pos_wing_i = [x + 1 for x in Con_em_i]
Con_neg_wing_i = [x - 1 for x in Con_em_i]
Ins_pos_wing_i = [x + 1 for x in Ins_em_i]
Ins_neg_wing_i = [x - 1 for x in Ins_em_i]

Con_pos_wing_Lam = [Lam_range[x] for x in Con_pos_wing_i]
Con_neg_wing_Lam = [Lam_range[x] for x in Con_neg_wing_i]
Ins_pos_wing_Lam = [Lam_range[x] for x in Ins_pos_wing_i]
Ins_neg_wing_Lam = [Lam_range[x] for x in Ins_neg_wing_i]

Old_Con_del_Lam = [((y - x) / 2) for x, y in zip(Con_pos_wing_Lam, Con_neg_wing_Lam)]
Old_Ins_del_Lam = [((y - x) / 2) for x, y in zip(Ins_pos_wing_Lam, Ins_neg_wing_Lam)]

#Creates R~1000 emission line width
New_Con_del_Lam = [x / 1000 for x in Lam_Con_em]
New_Ins_del_Lam = [x / 1000 for x in Lam_Ins_em]

#Finds ratio between old and new line width resolutions
Con_del_ratio = [x / y for x, y in zip(Old_Con_del_Lam, New_Con_del_Lam)]
Ins_del_ratio = [x / y for x, y in zip(Old_Ins_del_Lam, New_Ins_del_Lam)]

#Gets new SPD emission line and range for new resolution
New_SPD_Con_em = [x * y for x, y in zip(SPD_Con_em, Con_del_ratio)]
New_SPD_Ins_em = [x * y for x, y in zip(SPD_Ins_em, Ins_del_ratio)]

New_SPD_Con_range = [x for x in SPD_Con_range]
New_SPD_Ins_range = [x for x in SPD_Ins_range]

#Scales new range to correct emission line positions
for i in range(len(New_SPD_Con_range)):
	for j in range(len(Con_em_i)):
		if i == Con_em_i0[j]:
			New_SPD_Con_range[i] = New_SPD_Con_em[j]
for i in range(len(New_SPD_Ins_range)):
	for j in range(len(Ins_em_i)):
		if i == Ins_em_i[j]:
			New_SPD_Ins_range[i] = New_SPD_Ins_em[j]

#Finds new line width wings 
New_Con_Lam_pos_wing = [x + y for x, y in zip(Lam_Con_em, New_Con_del_Lam)]
New_Con_Lam_neg_wing = [x - y for x, y in zip(Lam_Con_em, New_Con_del_Lam)]
New_Ins_Lam_pos_wing = [x + y for x, y in zip(Lam_Ins_em, New_Ins_del_Lam)]
New_Ins_Lam_neg_wing = [x - y for x, y in zip(Lam_Ins_em, New_Ins_del_Lam)]

#Finds new SPD wings 
New_Con_SPD_pos_wing = [New_SPD_Con_range[x - 1] for x in Con_em_i]
New_Con_SPD_neg_wing = [New_SPD_Con_range[x + 1] for x in Con_em_i]
New_Ins_SPD_pos_wing = [New_SPD_Ins_range[x - 1] for x in Ins_em_i]
New_Ins_SPD_neg_wing = [New_SPD_Ins_range[x + 1] for x in Ins_em_i]

#Adds in new scaled emission lines 
Lam_Con_range  = [x for x in Lam_range]
Lam_Ins_range  = [x for x in Lam_range]
Lam_Con_range += [x for x in New_Con_Lam_pos_wing]
Lam_Ins_range += [x for x in New_Ins_Lam_pos_wing]
Lam_Con_range += [x for x in New_Con_Lam_neg_wing]
Lam_Ins_range += [x for x in New_Ins_Lam_neg_wing]

#Sorts ranges into numerical order
Lam_Con_range.sort()
Lam_Ins_range.sort()
Lam_Con_range[:] = Lam_Con_range[::-1]
Lam_Ins_range[:] = Lam_Ins_range[::-1]

#Finds new emission line positions
New_Con_em_i = [x for x in range(len(Lam_Con_range)) 
					  for y in range(len(Lam_Con_em)) 
					  if Lam_Con_range[x] == Lam_Con_em[y]]
New_Ins_em_i = [x for x in range(len(Lam_Ins_range)) 
					  for y in range(len(Lam_Ins_em)) 
					  if Lam_Ins_range[x] == Lam_Ins_em[y]]

#Creates final SPD range with emisson lines scaled to R~1000 
Final_SPD_Con_range  = [x for x in New_SPD_Con_range[ : Con_em_i[0] ]]
for i in range(len(New_SPD_Con_em) - 1):
	Final_SPD_Con_range += [New_Con_SPD_pos_wing[i], New_SPD_Con_em[i], New_Con_SPD_neg_wing[i]]
	Final_SPD_Con_range += [x for x in New_SPD_Con_range[Con_em_i[i] + 1 : Con_em_i[i+1] ]]

Final_SPD_Con_range += [New_Con_SPD_pos_wing[-1], New_SPD_Con_em[-1], New_Con_SPD_neg_wing[-1]]
Final_SPD_Con_range += [x for x in New_SPD_Con_range[Con_em_i[-1] + 1 : ]]

Final_SPD_Ins_range  = [x for x in New_SPD_Ins_range[ : Ins_em_i[0] ]]
for i in range(len(New_SPD_Ins_em) - 1):
	Final_SPD_Ins_range += [New_Ins_SPD_pos_wing[i], New_SPD_Ins_em[i], New_Ins_SPD_neg_wing[i]]
	Final_SPD_Ins_range += [x for x in New_SPD_Ins_range[Ins_em_i[i] + 1 : Ins_em_i[i+1] ]]

Final_SPD_Ins_range += [New_Ins_SPD_pos_wing[-1], New_SPD_Ins_em[-1], New_Ins_SPD_neg_wing[-1]]
Final_SPD_Ins_range += [x for x in New_SPD_Ins_range[Ins_em_i[-1] + 1 : ]]

#Scaled significance detection testing
"""
#Plots scaled emission lines
f, (ax1, ax2) = plt.subplots(2, figsize=(7,8))

ax1.set_xscale("log", nonposy="clip")
ax1.set_yscale("log", nonposy="clip")
ax2.set_xscale("log", nonposy="clip")
ax2.set_yscale("log", nonposy="clip")

ax1.plot(Lam_Con_em, New_SPD_Con_em, 'm.')
ax1.plot(Lam_Con_range, New_SPD_Con_range2, 'k-')
ax2.plot(Lam_Ins_em, New_SPD_Ins_em, 'm.')
ax2.plot(Lam_Ins_range, New_SPD_Ins_range2, 'k-')

ax1.set_xlim(1000, 30000)
ax1.set_ylim(1e27, 1e30)
ax2.set_xlim(1000, 30000)
ax2.set_ylim(1e25, 1e28)

ax2.set_xlabel(r"$\lambda\,(\rm{\AA})$", fontsize=20)
ax1.set_ylabel(r"$\rm{Power\,Density}\,(\rm{ergs\,s^{-1}\,Hz^{-1}})$", fontsize=20)
ax2.set_ylabel(r"$\rm{Power\,Density}\,(\rm{ergs\,s^{-1}\,Hz^{-1}})$", fontsize=20)

f.subplots_adjust(hspace=0.1)

plt.savefig(rootdir  + "scaled_emLines" + ".eps", bbox_inches="tight")
plt.close(f)
"""

#---------------------------------------------------------------------------------------------------------------------#
"""
#Comparison to Elizabeth values
print "Comparing with Elizabeth values..."
print

#Comparison to the new_g_... files 
lamlines = [2014, 1908, 1900]

test   = [abs(x - lamlines[0]) < (3840 / (1+5.913) / 2) for x in Lambda_A]
testa  = [abs(x - lamlines[1]) < (3840 / (1+6.295) / 2) for x in Lambda_A]
testb  = [abs(x - lamlines[2]) < (3840 / (1+6.327) / 2) for x in Lambda_A]

test2  = [i for i in range(len(test))  if test[i] == True ]
test2a = [i for i in range(len(testa)) if testa[i] == True]
test2b = [i for i in range(len(testb)) if testb[i] == True]

test3  = [Lambda_A[x] for x in test2 ]
test3a = [Lambda_A[x] for x in test2a]
test3b = [Lambda_A[x] for x in test2b]

test4  = [SPD_Con[x] for x in test2 ]
test4a = [SPD_Con[x] for x in test2a]
test4b = [SPD_Con[x] for x in test2b]

test5  = [SPD_Ins[x] for x in test2 ]
test5a = [SPD_Ins[x] for x in test2a]
test5b = [SPD_Ins[x] for x in test2b]

testH  = [((SOL * 1e10) / x) for x in test3 ]
testHa = [((SOL * 1e10) / x) for x in test3a]
testHb = [((SOL * 1e10) / x) for x in test3b]

#ergs/s/Hz to ergs/s/A
test6  = [((x * y) / z) for x, y, z in zip(test4, testH, test3)]
test6a = [((x * y) / z) for x, y, z in zip(test4a, testHa, test3a)]
test6b = [((x * y) / z) for x, y, z in zip(test4b, testHb, test3b)]

test7  = [((x * y) / z) for x, y, z in zip(test5, testH, test3)]
test7a = [((x * y) / z) for x, y, z in zip(test5a, testHa, test3a)]
test7b = [((x * y) / z) for x, y, z in zip(test5b, testHb, test3b)]

print "Continuous"
print "2014A erg/s/A: ", np.log10(np.mean(test6))
print "1908A erg/s/A: ", np.log10(np.mean(test6a))
print "1900A erg/s/A: ", np.log10(np.mean(test6b))
print
print "Instant"
print "2014A erg/s/A: ", np.log10(np.mean(test7))
print "1908A erg/s/A: ", np.log10(np.mean(test7a))
print "1900A erg/s/A: ", np.log10(np.mean(test7b))
"""

#---------------------------------------------------------------------------------------------------------------------#
#Redshift scaling
print "Redshift scaling of old resolution..."
print
#Host cosmology
z  = [5.913, 6.295, 6.327]
dL = [58076.8897e6, 62455.5989e6, 62823.7951e6]
dL = [x * pc_to_cm for x in dL]

#cm2 at z
z_cm2 = [(4 * np.pi * x**2) for x in dL]

#Scales SPD ergs/s/Hz to SFD ergs/s/Hz/cm2 at host redshift
SFD_Con_13 = [((x * (1+z[0])) / z_cm2[0]) for x in SPD_Con]
SFD_Con_05 = [((x * (1+z[1])) / z_cm2[1]) for x in SPD_Con]
SFD_Con_14 = [((x * (1+z[2])) / z_cm2[2]) for x in SPD_Con]

SFD_Ins_13 = [((x * (1+z[0])) / z_cm2[0]) for x in SPD_Ins]
SFD_Ins_05 = [((x * (1+z[1])) / z_cm2[1]) for x in SPD_Ins]
SFD_Ins_14 = [((x * (1+z[2])) / z_cm2[2]) for x in SPD_Ins]

#Scales zSFD from ergs/s/Hz/cm2 to nJy
SFD_Con_13 = [((x / Jy_to_erg) * 1e9) for x in SFD_Con_13]
SFD_Con_05 = [((x / Jy_to_erg) * 1e9) for x in SFD_Con_05]
SFD_Con_14 = [((x / Jy_to_erg) * 1e9) for x in SFD_Con_14]

SFD_Ins_13 = [((x / Jy_to_erg) * 1e9) for x in SFD_Ins_13]
SFD_Ins_05 = [((x / Jy_to_erg) * 1e9) for x in SFD_Ins_05]
SFD_Ins_14 = [((x / Jy_to_erg) * 1e9) for x in SFD_Ins_14]

#Cut off zSFD prior Lyman-a break and scale wavelength to host redshifts
for i in range(len(Lambda_A)):
	if Lambda_A[i] <= 1215:
		SFD_Con_13[i] = 0.00001
		SFD_Con_05[i] = 0.00001
		SFD_Con_14[i] = 0.00001
		SFD_Ins_13[i] = 0.00001
		SFD_Ins_05[i] = 0.00001
		SFD_Ins_14[i] = 0.00001
  
Lambda_A_13 = [(x * (1+z[0])) for x in Lambda_A]
Lambda_A_05 = [(x * (1+z[1])) for x in Lambda_A]
Lambda_A_14 = [(x * (1+z[2])) for x in Lambda_A]

#---------------------------------------------------------------------------------------------------------------------#
#Redshift scaling with new resolution
print "Redshift scaling of new resolution..."
print
#Scales SPD ergs/s/Hz to SFD ergs/s/Hz/cm2 at host redshift
New_SFD_Con_13 = [((x * (1+z[0])) / z_cm2[0]) for x in Final_SPD_Con_range]
New_SFD_Con_05 = [((x * (1+z[1])) / z_cm2[1]) for x in Final_SPD_Con_range]
New_SFD_Con_14 = [((x * (1+z[2])) / z_cm2[2]) for x in Final_SPD_Con_range]

New_SFD_Ins_13 = [((x * (1+z[0])) / z_cm2[0]) for x in Final_SPD_Ins_range]
New_SFD_Ins_05 = [((x * (1+z[1])) / z_cm2[1]) for x in Final_SPD_Ins_range]
New_SFD_Ins_14 = [((x * (1+z[2])) / z_cm2[2]) for x in Final_SPD_Ins_range]

#Scales zSPD from ergs/s/Hz/cm2 to nJy
New_SFD_Con_13 = [((x / Jy_to_erg) * 1e9) for x in New_SFD_Con_13]
New_SFD_Con_05 = [((x / Jy_to_erg) * 1e9) for x in New_SFD_Con_05]
New_SFD_Con_14 = [((x / Jy_to_erg) * 1e9) for x in New_SFD_Con_14]

New_SFD_Ins_13 = [((x / Jy_to_erg) * 1e9) for x in New_SFD_Ins_13]
New_SFD_Ins_05 = [((x / Jy_to_erg) * 1e9) for x in New_SFD_Ins_05]
New_SFD_Ins_14 = [((x / Jy_to_erg) * 1e9) for x in New_SFD_Ins_14]

#Cut off zSFD prior Lyman-a break and scale wavelength to host redshifts
for i in range(len(Lam_Con_range)):
	if Lam_Con_range[i] <= 1215:
		New_SFD_Con_13[i] = 0.00001
		New_SFD_Con_05[i] = 0.00001
		New_SFD_Con_14[i] = 0.00001
for i in range(len(Lam_Ins_range)):
	if Lam_Ins_range[i] <= 1215:
		New_SFD_Ins_13[i] = 0.00001
		New_SFD_Ins_05[i] = 0.00001
		New_SFD_Ins_14[i] = 0.00001

Con_Lambda_A_13 = [(x * (1+z[0])) for x in Lam_Con_range]
Con_Lambda_A_05 = [(x * (1+z[1])) for x in Lam_Con_range]
Con_Lambda_A_14 = [(x * (1+z[2])) for x in Lam_Con_range]

Ins_Lambda_A_13 = [(x * (1+z[0])) for x in Lam_Ins_range]
Ins_Lambda_A_05 = [(x * (1+z[1])) for x in Lam_Ins_range]
Ins_Lambda_A_14 = [(x * (1+z[2])) for x in Lam_Ins_range]

#---------------------------------------------------------------------------------------------------------------------#
#Flux scaling
print "Flux scaling..."
print
#Host SFD in nJy and observed wavelength
Host_A       = 13922.807
Host_SFD     = [105.43, 34.33, 14.95]
Host_SFD_err = [14.58, 6.26, 4.5]

#Finds filter bandpass range 
BP_13 = [abs(x - Host_A) < (3840 / 2) for x in Lambda_A_13]
BP_05 = [abs(x - Host_A) < (3840 / 2) for x in Lambda_A_05]
BP_14 = [abs(x - Host_A) < (3840 / 2) for x in Lambda_A_14]

BP_13_index = [i for i in range(len(BP_13)) if BP_13[i] == True]
BP_05_index = [i for i in range(len(BP_05)) if BP_05[i] == True]
BP_14_index = [i for i in range(len(BP_14)) if BP_14[i] == True]

#Finds wavelengths and SFD within bandpass range
BP_A_13 = [Lambda_A_13[x] for x in BP_13_index]
BP_A_05 = [Lambda_A_05[x] for x in BP_05_index]
BP_A_14 = [Lambda_A_14[x] for x in BP_14_index]

BP_SFD_Con_13 = [SFD_Con_13[x] for x in BP_13_index]
BP_SFD_Con_05 = [SFD_Con_05[x] for x in BP_05_index]
BP_SFD_Con_14 = [SFD_Con_14[x] for x in BP_14_index]

BP_SFD_Ins_13 = [SFD_Ins_13[x] for x in BP_13_index]
BP_SFD_Ins_05 = [SFD_Ins_05[x] for x in BP_05_index]
BP_SFD_Ins_14 = [SFD_Ins_14[x] for x in BP_14_index]

#Finds mean of SFD within bandpass range
BP_SFD_Con_13_mean = np.mean(BP_SFD_Con_13)
BP_SFD_Con_05_mean = np.mean(BP_SFD_Con_05)
BP_SFD_Con_14_mean = np.mean(BP_SFD_Con_14)

BP_SFD_Ins_13_mean  = np.mean(BP_SFD_Ins_13)
BP_SFD_Ins_05_mean  = np.mean(BP_SFD_Ins_05)
BP_SFD_Ins_14_mean  = np.mean(BP_SFD_Ins_14)

#Finds scale factor between mean bandpass SFD and host SFD
Con_Scale_13 = Host_SFD[0] / BP_SFD_Con_13_mean
Con_Scale_05 = Host_SFD[1] / BP_SFD_Con_05_mean
Con_Scale_14 = Host_SFD[2] / BP_SFD_Con_14_mean

Ins_Scale_13 = Host_SFD[0] / BP_SFD_Ins_13_mean
Ins_Scale_05 = Host_SFD[1] / BP_SFD_Ins_05_mean
Ins_Scale_14 = Host_SFD[2] / BP_SFD_Ins_14_mean

#Scales SFD to host SFD
New_SFD_Con_13 = [x * Con_Scale_13 for x in New_SFD_Con_13]
New_SFD_Con_05 = [x * Con_Scale_05 for x in New_SFD_Con_05]
New_SFD_Con_14 = [x * Con_Scale_14 for x in New_SFD_Con_14]

New_SFD_Ins_13 = [x * Ins_Scale_13 for x in New_SFD_Ins_13]
New_SFD_Ins_05 = [x * Ins_Scale_05 for x in New_SFD_Ins_05]
New_SFD_Ins_14 = [x * Ins_Scale_14 for x in New_SFD_Ins_14]

print "Con CV: ", Con_Scale_13, Con_Scale_05, Con_Scale_14
print
print "Con Err: "
print ((Host_SFD[0] + Host_SFD_err[0]) / BP_SFD_Con_13_mean) - Con_Scale_13
print ((Host_SFD[1] + Host_SFD_err[1]) / BP_SFD_Con_05_mean) - Con_Scale_05
print ((Host_SFD[2] + Host_SFD_err[2]) / BP_SFD_Con_14_mean) - Con_Scale_14
print
print "Ins CV: ", Ins_Scale_13, Ins_Scale_05, Ins_Scale_14
print
print "Ins Err: "
print ((Host_SFD[0] + Host_SFD_err[0]) / BP_SFD_Ins_13_mean) - Ins_Scale_13
print ((Host_SFD[1] + Host_SFD_err[1]) / BP_SFD_Ins_05_mean) - Ins_Scale_05
print ((Host_SFD[2] + Host_SFD_err[2]) / BP_SFD_Ins_14_mean) - Ins_Scale_14
print

#---------------------------------------------------------------------------------------------------------------------#
#JWST Dictionaries
NIRSpec_27_Dict = {}
NIRSpec_27_Dict.setdefault("micron_I", [])
NIRSpec_27_Dict.setdefault("Jy_I", [])
NIRSpec_27_Dict.setdefault("flux_I", [])
NIRSpec_27_Dict.setdefault("micron_II", [])
NIRSpec_27_Dict.setdefault("Jy_II", [])
NIRSpec_27_Dict.setdefault("flux_II", [])
NIRSpec_27_Dict.setdefault("micron_III", [])
NIRSpec_27_Dict.setdefault("Jy_III", [])
NIRSpec_27_Dict.setdefault("flux_III", [])

NIRSpec_10_Dict = {}
NIRSpec_10_Dict.setdefault("micron_I", [])
NIRSpec_10_Dict.setdefault("Jy_I", [])
NIRSpec_10_Dict.setdefault("flux_I", [])
NIRSpec_10_Dict.setdefault("micron_II", [])
NIRSpec_10_Dict.setdefault("Jy_II", [])
NIRSpec_10_Dict.setdefault("flux_II", [])
NIRSpec_10_Dict.setdefault("micron_III", [])
NIRSpec_10_Dict.setdefault("Jy_III", [])
NIRSpec_10_Dict.setdefault("flux_III", [])

NIRSpec_Dict = {}
NIRSpec_Dict.setdefault("micron", [])
NIRSpec_Dict.setdefault("Jy", [])
NIRSpec_Dict.setdefault("flux", [])

#---------------------------------------------------------------------------------------------------------------------#
#JWST files
print "Loading JWST files..."
print
#NIRSpec
R2700_file = urllib2.urlopen("http://www.stsci.edu/jwst/instruments/nirspec/sensitivity/R2700_sens.txt") 
R1000_file = urllib2.urlopen("http://www.stsci.edu/jwst/instruments/nirspec/sensitivity/R1000_sens.txt") 
R100_file  = urllib2.urlopen("http://www.stsci.edu/jwst/instruments/nirspec/sensitivity/R100_sens.txt") 

R2700_temp_I   = []
R1000_temp_I   = []
R2700_temp_II  = []
R1000_temp_II  = []
R2700_temp_III = []
R1000_temp_III = []
R100_temp      = []

for line in R2700_file:
	if line.strip() == "Band I":
		break
for line in R2700_file:
	if line.strip() == "Band II":
		break
	R2700_temp_I.append(line.strip())
for line in R2700_file:
	if line.strip() == "Band III":
		break
	R2700_temp_II.append(line.strip())
for line in R2700_file:
	if line.strip() == "":
		break
	R2700_temp_III.append(line.strip())

for i in range(len(R2700_temp_I)):
	if not R2700_temp_I[i] == "":
		NIRSpec_27_Dict["micron_I"].append(float(R2700_temp_I[i].split()[0]))
		NIRSpec_27_Dict["Jy_I"].append(float(R2700_temp_I[i].split()[2]))
		NIRSpec_27_Dict["flux_I"].append(float(R2700_temp_I[i].split()[3]))
for i in range(len(R2700_temp_II)):
	if not R2700_temp_II[i] == "":
		NIRSpec_27_Dict["micron_II"].append(float(R2700_temp_II[i].split()[0]))
		NIRSpec_27_Dict["Jy_II"].append(float(R2700_temp_II[i].split()[2]))
		NIRSpec_27_Dict["flux_II"].append(float(R2700_temp_II[i].split()[3]))
for i in range(len(R2700_temp_III)):
	if not R2700_temp_III[i] == "":
		NIRSpec_27_Dict["micron_III"].append(float(R2700_temp_III[i].split()[0]))
		NIRSpec_27_Dict["Jy_III"].append(float(R2700_temp_III[i].split()[2]))
		NIRSpec_27_Dict["flux_III"].append(float(R2700_temp_III[i].split()[3]))

for line in R1000_file:
	if line.strip() == "Band I":
		break
for line in R1000_file:
	if line.strip() == "Band II":
		break
	R1000_temp_I.append(line.strip())
for line in R1000_file:
	if line.strip() == "Band III":
		break
	R1000_temp_II.append(line.strip())
for line in R1000_file:
	if line.strip() == "":
		break
	R1000_temp_III.append(line.strip())

for i in range(len(R1000_temp_I)):
	if not R1000_temp_I[i] == "":
		NIRSpec_10_Dict["micron_I"].append(float(R1000_temp_I[i].split()[0]))
		NIRSpec_10_Dict["Jy_I"].append(float(R1000_temp_I[i].split()[2]))
		NIRSpec_10_Dict["flux_I"].append(float(R1000_temp_I[i].split()[3]))
for i in range(len(R1000_temp_II)):
	if not R1000_temp_II[i] == "":
		NIRSpec_10_Dict["micron_II"].append(float(R1000_temp_II[i].split()[0]))
		NIRSpec_10_Dict["Jy_II"].append(float(R1000_temp_II[i].split()[2]))
		NIRSpec_10_Dict["flux_II"].append(float(R1000_temp_II[i].split()[3]))
for i in range(len(R1000_temp_III)):
	if not R1000_temp_III[i] == "":
		NIRSpec_10_Dict["micron_III"].append(float(R1000_temp_III[i].split()[0]))
		NIRSpec_10_Dict["Jy_III"].append(float(R1000_temp_III[i].split()[2]))
		NIRSpec_10_Dict["flux_III"].append(float(R1000_temp_III[i].split()[3]))

for line in R100_file:
	if line.strip() == "(micron)   R       (nJy)  (erg/s/cm2)":
		break
for line in R100_file:
	R100_temp.append(line.strip())

for i in range(len(R100_temp)):
	if not R100_temp[i] == "":
		NIRSpec_Dict["micron"].append(float(R100_temp[i].split()[0]))
		NIRSpec_Dict["Jy"].append(float(R100_temp[i].split()[2]))
		NIRSpec_Dict["flux"].append(float(R100_temp[i].split()[3]))

#MIRI
MIRI_um = [5.6, 7.7, 10.0]#, 11.3, 12.8, 15.0, 18.0, 21.0, 25.5]
MIRI_fl = [0.2e-6, 0.28e-6, 0.7e-6]#, 1.7e-6, 1.4e-6, 1.8e-6, 4.3e-6, 8.6e-6, 28e-6] 

#---------------------------------------------------------------------------------------------------------------------#
#JWST conversion
print "JWST conversions..."
print
#NIRSpec
NIRSpec_27_A_I   = []
NIRSpec_27_A_II  = []
NIRSpec_27_A_III = []
NIRSpec_10_A_I   = []
NIRSpec_10_A_II  = []
NIRSpec_10_A_III = []
NIRSpec_A    = []

for i in range(len(NIRSpec_27_Dict["micron_I"])):
	NIRSpec_27_A_I += [NIRSpec_27_Dict.get("micron_I")[i] * 10000]
for i in range(len(NIRSpec_27_Dict["micron_II"])):
	NIRSpec_27_A_II += [NIRSpec_27_Dict.get("micron_II")[i] * 10000]
for i in range(len(NIRSpec_27_Dict["micron_III"])):
	NIRSpec_27_A_III += [NIRSpec_27_Dict.get("micron_III")[i] * 10000]

for i in range(len(NIRSpec_10_Dict["micron_I"])):
	NIRSpec_10_A_I += [NIRSpec_10_Dict.get("micron_I")[i] * 10000]
for i in range(len(NIRSpec_10_Dict["micron_II"])):
	NIRSpec_10_A_II += [NIRSpec_10_Dict.get("micron_II")[i] * 10000]
for i in range(len(NIRSpec_10_Dict["micron_III"])):
	NIRSpec_10_A_III += [NIRSpec_10_Dict.get("micron_III")[i] * 10000]

for i in range(len(NIRSpec_Dict["micron"])):
	NIRSpec_A += [NIRSpec_Dict.get("micron")[i] * 10000]

NIRSpec_27_F_I   = []
NIRSpec_27_F_II  = []
NIRSpec_27_F_III = []
NIRSpec_10_F_I   = []
NIRSpec_10_F_II  = []
NIRSpec_10_F_III = []
NIRSpec_F    = []

for i in range(len(NIRSpec_27_Dict["Jy_I"])):
	NIRSpec_27_F_I += [NIRSpec_27_Dict.get("Jy_I")[i] / 3.3]
for i in range(len(NIRSpec_27_Dict["Jy_II"])):
	NIRSpec_27_F_II += [NIRSpec_27_Dict.get("Jy_II")[i] / 3.3]
for i in range(len(NIRSpec_27_Dict["Jy_III"])):
	NIRSpec_27_F_III += [NIRSpec_27_Dict.get("Jy_III")[i] / 3.3]

for i in range(len(NIRSpec_10_Dict["Jy_I"])):
	NIRSpec_10_F_I += [NIRSpec_10_Dict.get("Jy_I")[i] / 3.3]
for i in range(len(NIRSpec_10_Dict["Jy_II"])):
	NIRSpec_10_F_II += [NIRSpec_10_Dict.get("Jy_II")[i] / 3.3]
for i in range(len(NIRSpec_10_Dict["Jy_III"])):
	NIRSpec_10_F_III += [NIRSpec_10_Dict.get("Jy_III")[i] / 3.3]

for i in range(len(NIRSpec_Dict["Jy"])):
	NIRSpec_F += [NIRSpec_Dict.get("Jy")[i] / 3.3]

#MIRI
MIRI_A = [x * 10000 for x in MIRI_um]
MIRI_F = [(x * 1e9) / 2.8 for x in MIRI_fl]

#---------------------------------------------------------------------------------------------------------------------#
#Graphing
print "Plotting SED..."
print
f, (ax1, ax2) = plt.subplots(2, figsize=(5.65685,8))

ax1.set_xscale("log", nonposy="clip")
ax1.set_yscale("log", nonposy="clip")
ax2.set_xscale("log", nonposy="clip")
ax2.set_yscale("log", nonposy="clip")

ax1.plot(Con_Lambda_A_13, New_SFD_Con_13, color='k', linestyle='-', linewidth=0.5)
ax1.plot(Con_Lambda_A_05, New_SFD_Con_05, color='0.25', linestyle='-', linewidth=0.5)
ax1.plot(Con_Lambda_A_14, New_SFD_Con_14, color='0.50', linestyle='-', linewidth=0.5)

ax2.plot(Ins_Lambda_A_13, New_SFD_Ins_13, color='k', linestyle='-', linewidth=0.5)
ax2.plot(Ins_Lambda_A_05, New_SFD_Ins_05, color='0.25', linestyle='-', linewidth=0.5)
ax2.plot(Ins_Lambda_A_14, New_SFD_Ins_14, color='0.50', linestyle='-', linewidth=0.5)

ax1.errorbar(Host_A, Host_SFD[0], yerr=Host_SFD_err[0], ecolor='r', marker="s", 
				ls="None", capsize=0, mec="r", mfc="r", markersize=6)
ax1.errorbar(Host_A, Host_SFD[1], yerr=Host_SFD_err[1], ecolor='r', marker="o", 
				ls="None", capsize=0, mec="r", mfc="r", markersize=6)
ax1.errorbar(Host_A, Host_SFD[2], yerr=Host_SFD_err[2], ecolor='r', marker="D", 
				ls="None", capsize=0, mec="r", mfc="r", markersize=6)

ax2.errorbar(Host_A, Host_SFD[0], yerr=Host_SFD_err[0], ecolor='r', marker="s", 
				ls="None", capsize=0, mec="r", mfc="r", markersize=6)
ax2.errorbar(Host_A, Host_SFD[1], yerr=Host_SFD_err[1], ecolor='r', marker="o", 
				ls="None", capsize=0, mec="r", mfc="r", markersize=6)
ax2.errorbar(Host_A, Host_SFD[2], yerr=Host_SFD_err[2], ecolor='r', marker="D", 
				ls="None", capsize=0, mec="r", mfc="r", markersize=6)

ax1.plot(NIRSpec_10_A_I, NIRSpec_10_F_I, "c--", linewidth=2)
ax1.plot(NIRSpec_10_A_II, NIRSpec_10_F_II, "c--", linewidth=2)
ax1.plot(NIRSpec_10_A_III, NIRSpec_10_F_III, "c--", linewidth=2)

ax2.plot(NIRSpec_10_A_I, NIRSpec_10_F_I, "c--", linewidth=2)
ax2.plot(NIRSpec_10_A_II, NIRSpec_10_F_II, "c--", linewidth=2)
ax2.plot(NIRSpec_10_A_III, NIRSpec_10_F_III, "c--", linewidth=2)

ax1.text(6e4, 6e3, "z=5.913", color="k", fontsize=12)
ax1.text(6e4, 3.5e3, "z=6.295", color="0.25", fontsize=12)
ax1.text(6e4, 2e3, "z=6.327", color="0.50", fontsize=12)
ax1.text(1e4, 6e3, r"$100\,\rm{Myr}$ Continuous", color="k", fontsize=12)
ax1.text(1e4, 3.5e3, r"Model (A)", color="k", fontsize=12)

ax2.text(6e4, 6e3, "z=5.913", color="k", fontsize=12)
ax2.text(6e4, 3.5e3, "z=6.295", color="0.25", fontsize=12)
ax2.text(6e4, 2e3, "z=6.327", color="0.50", fontsize=12)
ax2.text(1e4, 6e3, r"$10\,\rm{Myr}$ Starburst", color="k", fontsize=12)
ax2.text(1e4, 3.5e3, r"Model (B)", color="k", fontsize=12)

ax1.set_xlim(1000 * (1+z[0]), 1e5)
ax1.set_ylim(2, 20000)
ax2.set_xlim(1000 * (1+z[0]), 1e5)
ax2.set_ylim(2, 20000)

ax1.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax2.yaxis.set_major_formatter(ticker.ScalarFormatter())

ax2.set_xlabel(r"$\lambda\,(\rm{\AA})$", fontsize=20)
ax1.set_ylabel(r"$\rm{Flux\,Density}\,(\rm{nJy})$", fontsize=20)
ax2.set_ylabel(r"$\rm{Flux\,Density}\,(\rm{nJy})$", fontsize=20)

f.subplots_adjust(hspace=0.1)

plt.savefig(rootdir  + "SED_jun08" + ".eps", bbox_inches="tight")
plt.close(f)

#---------------------------------------------------------------------------------------------------------------------#
#End of program
print "Program executed correctly, ending process."
print
