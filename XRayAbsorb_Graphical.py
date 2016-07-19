#!/usr/bin/env python

#Python Script for Graphical Presentation of data from XRayAbsorb.py
#v1.00 
#Author: Joseph T. W. McGuire, University of Leicester, Department of Physics and Astronomy, XROA Group 
#Contact: jtwm1@leicester.ac.uk
#Last Modified: Oct 2015

#---------------------------------------------------------------------------------------------------------------------#
#Imports
import os, sys, glob, shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import tarfile as tf
import re, argparse
from xspec import *
from tables import *
from astropy.io import fits
from astropy.table import Table, vstack 
from scipy.stats.kde import gaussian_kde

#---------------------------------------------------------------------------------------------------------------------#
#User Input
parser = argparse.ArgumentParser(description="Graphical Presentation of selected data, uses FITs format")
parser.add_argument("rootdir", help="Directory path of the program, ensure directory string ends with /")
parser.add_argument("input", help="Enter the 'output' string from FITs file you wish to load. \
		'XRAb_FITSTable_' + output + '.fits' is an example.")
parser.add_argument("Graphs", nargs=4, choices=["Y", "N"], help="Select which FITs file you wish to load \
		for illustration, choice is {Y, N} and the order is Ecut, XRAb, Z, MCMC")
parser.add_argument("-ls", "--lsmc", help="Compare XRAb data with L/SMC data.") 
args = parser.parse_args() 

rootdir = args.rootdir
output  = args.input 

#---------------------------------------------------------------------------------------------------------------------#
#User input checks

#args.rootdir
if re.split("[/]", rootdir)[-1] != "":
	print "Root Directory path did not end with /"
	sys.exit()

#args.Graphs
if args.Graphs[0] == "Y" and args.Graphs[1] == "Y":
	print "Ecut and XRAb cannot be envoked together"
	sys.exit()
elif args.Graphs[0] == "Y" and args.Graphs[2] == "Y":
	print "Ecut and Z cannot be envoked together"
	sys.exit()
elif args.Graphs[0] == "Y" and args.Graphs[3] == "Y":
	print "Ecut and MCMC cannot be envoked together"
	sys.exit()
if args.Graphs[1] == "Y" and args.Graphs[2] == "Y" and args.Graphs[3] == "Y":
	print "XRAb, Z and MCMC cannot be envoked together"
	sys.exit()
if args.Graphs[1] != "Y" and args.lsmc:
	print "lsmc cannot be envoked with any option except XRAb"
	sys.exit()

#---------------------------------------------------------------------------------------------------------------------#
#Loads selected FITs files
print
print
print "Loading FITs files"
print
print

if args.Graphs[0] == "Y":
	mainFITs     = fits.open(rootdir + "Ecut_FITSTable_" + args.input + ".fits")
	mainContents = mainFITs[1].data

if args.Graphs[1] == "Y" and args.Graphs[2] == "N":
	mainFITs     = fits.open(rootdir + "XRAb_FITSTable_" + args.input + ".fits")
	mainContents = mainFITs[1].data
elif args.Graphs[2] == "Y" and args.Graphs[1] == "N":
	mainFITs     = fits.open(rootdir + "Z_FITSTable_" + args.input + ".fits")
	mainContents = mainFITs[1].data
elif args.Graphs[1] == "Y" and args.Graphs[2] == "Y":
	mainFITs     = fits.open(rootdir + "XRAb_FITSTable_" + args.input + ".fits")
	ZFITs        = fits.open(rootdir + "Z_FITSTable_" + args.input + ".fits")
	mainContents = mainFITs[1].data
	ZContents    = ZFITs[1].data

if args.Graphs[1] == "Y" and args.Graphs[3] == "Y":
	mainFITs     = fits.open(rootdir + "MCMC_XRAb_FITSTable_" + args.input + ".fits")
	mainContents = mainFITs[1].data
elif args.Graphs[2] == "Y" and args.Graphs[3] == "Y":
	mainFITs     = fits.open(rootdir + "MCMC_Z_FITSTable_" + args.input + ".fits")
	mainContents = mainFITs[1].data

if args.lsmc:
	lsmcFITs     = fits.open(rootdir + "asu.fit")
	lsmcContents = lsmcFITs[1].data

#---------------------------------------------------------------------------------------------------------------------#
#Plot Functions
def NHz_CV_plot(AX, Colour, Table):
	AX.errorbar(Table["z"] + 1, Table["NH_Int"] * 1e22, yerr=[Table["NH_IntErr_low"] * 1e22, Table["NH_IntErr_hi"] * 1e22], \
			ecolor=Colour, marker=".", color=Colour, capsize=0, ls="None")
def NHz_UL_plot(AX, Colour, Table):
	AX.errorbar(Table["z"] + 1, Table["NH_Int"] * 1e22, yerr=[Table["NH_IntErr_low"] * 1e22, Table["NH_IntErr_hi"] * 1e22], \
			lolims=True, ecolor=Colour, color=Colour, capsize=2, ls="None")

def NHIG_CV_plot(AX, Colour, Table):
	AX.errorbar(Table["NH_Int"] * 1e22, Table["NH_Gal"] * 1e22, xerr=[Table["NH_IntErr_low"] * 1e22, Table["NH_IntErr_hi"] * 1e22], \
			ecolor=Colour, marker=".", color=Colour, ls="None", capsize=0)
def NHIG_UL_plot(AX, Colour, Table):	
	AX.errorbar(Table["NH_Int"] * 1e22, Table["NH_Gal"] * 1e22, xerr=[Table["NH_IntErr_low"] * 1e22, Table["NH_IntErr_hi"] * 1e22], \
			xlolims=True, ecolor=Colour, color=Colour, ls="None", capsize=2)

def NHIP_CV_plot(AX, Colour, Table):
	AX.errorbar(Table["NH_Int"] * 1e22, Table["Gamma"], xerr=[Table["NH_IntErr_low"] * 1e22, Table["NH_IntErr_hi"] * 1e22], \
			yerr=[Table["GammaErr_low"], Table["GammaErr_hi"]], ecolor=Colour, marker=".", color=Colour, ls="None", capsize=0)
def NHIP_UL_plot(AX, Colour, Table):
	AX.errorbar(Table["NH_Int"] * 1e22, Table["Gamma"], xerr=[Table["NH_IntErr_low"] * 1e22, Table["NH_IntErr_hi"] * 1e22], \
			yerr=[Table["GammaErr_low"], Table["GammaErr_hi"]], lolims=True, xlolims=True, ecolor=Colour, color=Colour, ls="None", capsize=2)

def NHIF_CV_plot(AX, Colour, Table):
	AX.errorbar(Table["NH_Int"] * 1e22, Table["cFx"], xerr=[Table["NH_IntErr_low"] * 1e22, Table["NH_IntErr_hi"] * 1e22], \
			yerr=[Table["cFxErr_low"], Table["cfxErr_hi"]], ecolor=Colour, marker=".", color=Colour, ls="None", capsize=0)
def NHIF_UL_plot(AX, Colour, Table):
	AX.errorbar(Table["NH_Int"] * 1e22, Table["cFx"], xerr=[Table["NH_IntErr_low"] * 1e22, Table["NH_IntErr_hi"] * 1e22], \
			yerr=[Table["cFxErr_low"], Table["cfxErr_hi"]], lolims=True, xlolims=True, ecolor=Colour, color=Colour, ls="None", capsize=2)

def NHIp_CV_plot(AX, Colour, Table):
	AX.errorbar(Table["NH_Int"] * 1e22, Table["meanpc"], xerr=[Table["NH_IntErr_low"] * 1e22, Table["NH_IntErr_hi"] * 1e22], \
			ecolor=Colour, marker=".", color=Colour, ls="None", capsize=0)
def NHIp_UL_plot(AX, Colour, Table):
	AX.errorbar(Table["NH_Int"] * 1e22, Table["meanpc"], xerr=[Table["NH_IntErr_low"] * 1e22, Table["NH_IntErr_hi"] * 1e22], \
			xlolims=True, ecolor=Colour, color=Colour, ls="None", capsize=2)

def gaussian(Data, mean, variance):
	fx = (1.0 / np.sqrt(2 * np.pi * variance)) * np.exp(- ((Data - mean)**2 / (2 * variance)))

	return fx

def NH_Bin_plot(AX, z, zerr, Colour, Table):
	AX.errorbar(z, np.median(Table["NH_Int"]) * 1e22, xerr=zerr, ecolor=Colour, marker=".", mfc=Colour, mec=Colour, ls="None", capsize=0)
#---------------------------------------------------------------------------------------------------------------------#
#Spectral Count Histogram
if args.Graphs[3] == "N":
	print "Creating Spectral Count Histogram..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, figsize=(7, 4))

	histbin = np.linspace(0, max(mainContents["SCounts"]), (int(max(mainContents["SCounts"])) / 20) + 1)
	n, bins, patches = ax.hist(mainContents["SCounts"], histbin, histtype="stepfilled", color="k", alpha=0.5)

	ax1a = plt.axes([.55, .5, .3, .3])

	histbin = np.linspace(0, 100, 21)
	n, bins, patches = ax1a.hist(mainContents["SCounts"], histbin, histtype="stepfilled", color="g", alpha=0.5)

	ax.set_xlabel(r"$cts$")
	ax.set_ylabel(r"Number")

	if args.Graphs[0] == "Y":
		plt.savefig(rootdir + "Ecut_N(SC)Hist_" + output + ".eps", bbox_inches="tight")
	elif args.Graphs[1] == "Y" and args.Graphs[2] == "N":
		plt.savefig(rootdir + "XRAb_N(SC)Hist_" + output + ".eps", bbox_inches="tight")
	elif args.Graphs[1] == "N" and args.Graphs[2] == "Y":
		plt.savefig(rootdir + "Z_N(SC)Hist_" + output + ".eps", bbox_inches="tight")

	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#Redshift Histogram
if args.Graphs[0] == "N" and args.Graphs[3] == "N":
	print "Creating Redshift Histogram..."
	print
	print	

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, figsize=(7, 4))

	histbin = np.linspace(0, 10.0, 41)
	n, bins, patches = ax.hist(mainContents["z"], histbin, histtype="stepfilled", color="k", alpha=0.5)

	ax.set_xlabel(r"$z$")
	ax.set_ylabel(r"Number")

	if args.Graphs[1] == "Y" and args.Graphs[2] == "N":
		plt.savefig(rootdir + "XRAb_N(z)Hist_" + output + ".eps", bbox_inches="tight")
	elif args.Graphs[1] == "N" and args.Graphs[2] == "Y":
		plt.savefig(rootdir + "Z_N(z)Hist_" + output + ".eps", bbox_inches="tight")

	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#Redshift Comparison Histogram
if args.Graphs[1] == "Y" and args.Graphs[2] == "N" and args.Graphs[3] == "N":
	print "Creating Redshift Comparison Histogram..."
	print
	print	

	S = []
	C = []
	with open(rootdir + "StarlingList.txt", "r") as fs:
		for line in fs:
			S += [float(line.split()[1])]
	with open(rootdir + "CampanaList.txt", "r") as fc:
		for line in fc:
			C += [float(line.split()[1])]

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, figsize=(7, 4))

	histbin = np.linspace(0, 10.0, 41)
	n, bins, patches = ax.hist(S, histbin, color="m")
	n, bins, patches = ax.hist(C, histbin, color="g", rwidth=0.5)
	n, bins, patches = ax.hist(mainContents["z"], histbin, color="r", histtype="step")

	ax.set_xlabel(r"$z$")
	ax.set_ylabel(r"Number")

	plt.savefig(rootdir + "XRAb_N(z)Hist_Comparison_" + output + ".eps", bbox_inches="tight")

	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#Subclass Histogram
if args.Graphs[1] == "Y" and args.Graphs[2] == "N" and args.Graphs[3] == "N":
	print "Creating Redshift Subclass Histogram..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, figsize=(7, 4))

	L   = mainContents["Type"] == "L"
	S   = mainContents["Type"] == "S"
	SE  = mainContents["Type"] == "SE"
	SEP = mainContents["Type"] == "SEP"
	LP  = mainContents["Type"] == "LP"
	LD  = mainContents["Type"] == "LD"
	UL  = mainContents["Type"] == "UL"
	X   = mainContents["Type"] == "X"

	temptb1 = mainContents[L]
	temptb2 = mainContents[S]
	temptb3 = mainContents[SE]
	temptb4 = mainContents[SEP]
	temptb5 = mainContents[LP]
	temptb6 = mainContents[LD]
	temptb7 = mainContents[UL]
	temptb8 = mainContents[X]

	histbin = np.linspace(0, 10.0, 41)
	n, bins, patches = ax.hist([temptb1["z"], temptb2["z"], temptb3["z"], temptb4["z"], temptb5["z"], temptb6["z"], temptb7["z"], temptb8["z"]], \
								histbin, stacked=True, color=["k", "b", "g", "c", "c", "m", "r", "y"])

	ax.set_xlabel(r"$z$")
	ax.set_ylabel(r"Number")

	plt.savefig(rootdir + "XRAb_N(z)Hist_Subclass_" + output + ".eps", bbox_inches="tight")

	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#NH G Histogram
if args.Graphs[1] == "Y" and args.Graphs[2] == "N" and args.Graphs[3] == "N":
	print "Creating NH and G Histograms..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, 2, figsize=(14, 4))

	min_NInt = np.log10(1e19 * np.power(mainContents["z"], 2.4))

	Real   = np.log10(mainContents["NH_Int"] * 1e22) >= min_NInt
	temptb = mainContents[Real]

	CVnh = temptb["NH_IntErr_stat"]      == "FFFFFFFFF"
	CVg  = temptb["GammaErr_stat"] == "FFFFFFFFF"
	
	temptbnh = temptb[CVnh]
	temptbg  = temptb[CVg]

	histbin = np.linspace(19, 24, 26)
	n, bins, patches = ax[0].hist(np.log10(temptbnh["NH_Int"] * 1e22), histbin, color="k", histtype="step")
	n, bins, patches = ax[0].hist(np.log10(temptb["NH_Int"] * 1e22), histbin, color="r", histtype="step", linestyle=("dashed"))

	histbin2 = np.linspace(1, 4, 51)
	n, bins, patches = ax[1].hist(temptbg["Gamma"], histbin2, color="k", histtype="step")
	n, bins, patches = ax[1].hist(temptb["Gamma"], histbin2, color="r", histtype="step", linestyle=("dashed"))

	ax[0].set_xlabel(r"$Log_{10}$ $N_{H_{X, Intrinsic}} (cm^{-2})$")
	ax[0].set_ylabel(r"Number")
	ax[1].set_xlabel(r"Spectral Index $\Gamma$")
	ax[1].set_ylabel(r"Number")

	plt.savefig(rootdir + "XRAb_N(NH_G)_Hist_" + output + ".eps", bbox_inches="tight")

	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#NH G redshift Histogram
if args.Graphs[1] == "Y" and args.Graphs[2] == "N" and args.Graphs[3] == "N":
	print "Creating NH and G Redshift Histograms..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, 2, figsize=(14, 4))

	min_NInt = np.log10(1e19 * np.power(mainContents["z"], 2.4))

	Real   = np.log10(mainContents["NH_Int"] * 1e22) >= min_NInt
	temptb = mainContents[Real]

	CVnh = temptb["NH_IntErr_stat"]      == "FFFFFFFFF"
	CVg  = temptb["GammaErr_stat"] == "FFFFFFFFF"
	
	temptbnh = temptb[CVnh]
	temptbg  = temptb[CVg]

	zl2  = temptb["z"] < mainContents["z"].mean()
	zge2 = temptb["z"] >= mainContents["z"].mean()

	temptbzl2  = temptb[zl2]
	temptbzge2 = temptb[zge2]

	histbin = np.linspace(19, 24, 26)
	n, bins, patches = ax[0].hist(np.log10(temptb["NH_Int"] * 1e22), histbin, color="k", histtype="step")
	n, bins, patches = ax[0].hist([np.log10(temptbzl2["NH_Int"] * 1e22), np.log10(temptbzge2["NH_Int"] * 1e22)], histbin, color=["c","g"], stacked=True, histtype="stepfilled")

	histbin2 = np.linspace(1, 4, 51)
	n, bins, patches = ax[1].hist(temptb["Gamma"], histbin2, color="k", histtype="step")
	n, bins, patches = ax[1].hist([temptbzl2["Gamma"], temptbzge2["Gamma"]], histbin2, color=["c","g"], stacked=True, histtype="stepfilled")

	ax[0].set_xlabel(r"$Log_{10}$ $N_{H_{X, Intrinsic}} (cm^{-2})$")
	ax[0].set_ylabel(r"Number")
	ax[1].set_xlabel(r"Spectral Index $\Gamma$")
	ax[1].set_ylabel(r"Number")

	plt.savefig(rootdir + "XRAb_N(NH_G)_Redshift_Hist_" + output + ".eps", bbox_inches="tight")

	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#NH G Subclass Histogram
if args.Graphs[1] == "Y" and args.Graphs[2] == "N" and args.Graphs[3] == "N":
	print "Creating NH and G subclass Histograms..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, 2, figsize=(14, 4))

	min_NInt = np.log10(1e19 * np.power(mainContents["z"], 2.4))

	Real   = np.log10(mainContents["NH_Int"] * 1e22) >= min_NInt
	temptb = mainContents[Real]

	L   = temptb["Type"] == "L"
	S   = temptb["Type"] == "S"
	SE  = temptb["Type"] == "SE"
	SEP = temptb["Type"] == "SEP"
	LP  = temptb["Type"] == "LP"
	LD  = temptb["Type"] == "LD"
	UL  = temptb["Type"] == "UL"
	X   = temptb["Type"] == "X"

	temptb1 = temptb[L]
	temptb2 = temptb[S]
	temptb3 = temptb[SE]
	temptb4 = temptb[SEP]
	temptb5 = temptb[LP]
	temptb6 = temptb[LD]
	temptb7 = temptb[UL]
	temptb8 = temptb[X]

	histbin = np.linspace(19, 24, 26)
	n, bins, patches = ax[0].hist([np.log10(temptb1["NH_Int"] * 1e22), np.log10(temptb2["NH_Int"] * 1e22), np.log10(temptb3["NH_Int"] * 1e22), \
									np.log10(temptb4["NH_Int"] * 1e22), np.log10(temptb5["NH_Int"] * 1e22), np.log10(temptb6["NH_Int"] * 1e22), \
									np.log10(temptb7["NH_Int"] * 1e22), np.log10(temptb8["NH_Int"] * 1e22)], \
								histbin, stacked=True, color=["k", "b", "g", "c", "c", "m", "r", "y"], histtype="stepfilled")

	histbin2 = np.linspace(1, 4, 51)
	n, bins, patches = ax[1].hist([temptb1["Gamma"], temptb2["Gamma"], temptb3["Gamma"], temptb4["Gamma"], temptb5["Gamma"], \
								temptb6["Gamma"], temptb7["Gamma"], temptb8["Gamma"]], \
								histbin2, stacked=True, color=["k", "b", "g", "c", "c", "m", "r", "y"], histtype="stepfilled")

	ax[0].set_xlabel(r"$Log_{10}$ $N_{H_{X, Intrinsic}} (cm^{-2})$")
	ax[0].set_ylabel(r"Number")
	ax[1].set_xlabel(r"Spectral Index $\Gamma$")
	ax[1].set_ylabel(r"Number")

	plt.savefig(rootdir + "XRAb_Subclass_Hist_" + output + ".eps", bbox_inches="tight")

	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#XRAb PDEs
if args.Graphs[1] == "Y" and args.Graphs[2] == "N" and args.Graphs[3] == "N":
	print "Creating XRAb NH and G PDEs..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, 2, figsize=(14, 4))

	min_NInt = np.log10(1e19 * np.power(mainContents["z"], 2.4))

	Real   = np.log10(mainContents["NH_Int"] * 1e22) >= min_NInt
	temptb = mainContents[Real]

	CVnh = temptb["NH_IntErr_stat"] == "FFFFFFFFF"
	CVg  = temptb["GammaErr_stat"]  == "FFFFFFFFF"
	temptba = temptb[CVnh]
	temptbb = temptb[CVg]

	fullnh_mean = np.log10(temptb["NH_Int"]* 1e22).mean()
	fullnh_var  = np.log10(temptb["NH_Int"]* 1e22).var()
	fullg_mean  = temptb["Gamma"].mean()
	fullg_var   = temptb["Gamma"].var()

	CVnh_mean   = np.log10(temptba["NH_Int"]* 1e22).mean()
	CVnh_var    = np.log10(temptba["NH_Int"]* 1e22).var()
	CVg_mean    = temptbb["Gamma"].mean()
	CVg_var     = temptbb["Gamma"].var()

	ax[0].plot(np.log10(temptb["NH_Int"] * 1e22), gaussian(np.log10(temptb["NH_Int"] * 1e22), fullnh_mean, fullnh_var), 'r.')
	ax[0].axvline(fullnh_mean, c=u"r", ls="-")
	ax[0].plot(np.log10(temptba["NH_Int"] * 1e22), gaussian(np.log10(temptba["NH_Int"] * 1e22), CVnh_mean, CVnh_var), 'k.')
	ax[0].axvline(CVnh_mean, c=u"k", ls="-")
	ax[1].plot(temptb["Gamma"], gaussian(temptb["Gamma"], fullg_mean, fullg_var), 'r.')
	ax[1].axvline(fullg_mean, c=u"r", ls="-")
	ax[1].plot(temptbb["Gamma"], gaussian(temptbb["Gamma"], CVg_mean, CVg_var), 'k.')
	ax[1].axvline(CVg_mean, c=u"k", ls="-")

	ax[0].set_xlabel(r"$Log_{10}$ $N_{H_{X, Intrinsic}} (cm^{-2})$")
	ax[0].set_ylabel(r"PDE")
	ax[1].set_xlabel(r"Spectral Index $\Gamma$")
	ax[1].set_ylabel(r"PDE")

	plt.savefig(rootdir + "XRAb_PDE_Hist_" + output + ".eps", bbox_inches="tight")

	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#Z NH Histogram
if args.Graphs[1] == "N" and args.Graphs[2] == "Y" and args.Graphs[3] == "N":
	print "Creating Z NH Histogram..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, 2, figsize=(14, 4))

	min_NInt = np.log10(1e19 * np.power(mainContents["z"], 2.4))

	Real   = np.log10(mainContents["NH_Int"] * 1e22) >= min_NInt
	temptb = mainContents[Real]

	CV = temptb["Author"] == "Ccv"
	temptba = temptb[CV]

	full_mean = np.log10(temptb["NH_Int"]* 1e22).mean()
	full_var  = np.log10(temptb["NH_Int"]* 1e22).var()
	CV_mean   = np.log10(temptba["NH_Int"]* 1e22).mean()
	CV_var    = np.log10(temptba["NH_Int"]* 1e22).var()

	histbin = np.linspace(19, 24, 26)
	n, bins, patches = ax[0].hist(np.log10(temptba["NH_Int"] * 1e22), histbin, color="k", histtype="step")
	n, bins, patches = ax[0].hist(np.log10(temptb["NH_Int"] * 1e22), histbin, color="r", histtype="step", linestyle=("dashed"))

	ax[1].plot(np.log10(temptb["NH_Int"] * 1e22), gaussian(np.log10(temptb["NH_Int"] * 1e22), full_mean, full_var), 'r.')
	ax[1].axvline(full_mean, c=u"r", ls="-")
	ax[1].plot(np.log10(temptba["NH_Int"] * 1e22), gaussian(np.log10(temptba["NH_Int"] * 1e22), CV_mean, CV_var), 'k.')
	ax[1].axvline(CV_mean, c=u"k", ls="-")

	ax[0].set_xlabel(r"$Log_{10}$ $N_{H_{X, Intrinsic}} (cm^{-2})$")
	ax[0].set_ylabel(r"Number")
	ax[1].set_xlabel(r"$Log_{10}$ $N_{H_{X, Intrinsic}} (cm^{-2})$")
	ax[1].set_ylabel(r"PDE")

	plt.savefig(rootdir + "Z_N(NH)Hist_" + output + ".eps", bbox_inches="tight")
	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#Z G Histogram
if args.Graphs[1] == "N" and args.Graphs[2] == "Y" and args.Graphs[3] == "N":
	print "Creating Z G Histogram..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, 2, figsize=(14, 4))

	CV     = mainContents["GammaErr_stat"] == "FFFFFFFFF"
	temptb = mainContents[CV]

	full_mean  = mainContents["Gamma"].mean()
	full_var   = mainContents["Gamma"].var()
	CV_mean    = temptb["Gamma"].mean()
	CV_var     = temptb["Gamma"].var()

	histbin = np.linspace(1, 4, 51)
	n, bins, patches = ax[0].hist(temptb["Gamma"], histbin, color="k", histtype="step")
	n, bins, patches = ax[0].hist(mainContents["Gamma"], histbin, color="r", histtype="step", linestyle=("dashed"))

	ax[1].plot(mainContents["Gamma"], gaussian(mainContents["Gamma"], full_mean, full_var), 'r.')
	ax[1].axvline(full_mean, c=u"r", ls="-")
	ax[1].plot(temptb["Gamma"], gaussian(temptb["Gamma"], CV_mean, CV_var), 'k.')
	ax[1].axvline(CV_mean, c=u"k", ls="-")

	ax[0].set_xlabel(r"Spectral Index $\Gamma$")
	ax[0].set_ylabel(r"Number")
	ax[1].set_xlabel(r"Spectral Index $\Gamma$")
	ax[1].set_ylabel(r"PDE")

	plt.savefig(rootdir + "Z_N(G)Hist_" + output + ".eps", bbox_inches="tight")
	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#L/SMC Histogram
if args.lsmc:
	print "Creating L/SMC Comparison Histogram..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(2, 2, figsize=(14, 8))

	LS_RV   = lsmcContents["N_HI_"] > 0
	temptb1 = lsmcContents[LS_RV]

	LS_CV    = temptb1["l_N_HI_"] != "<"
	LS_UL    = temptb1["l_N_HI_"] == "<"
	temptb1a = temptb1[LS_CV]
	temptb1b = temptb1[LS_UL]

	L_CV = temptb1a["MC"] == "L"
	S_CV = temptb1a["MC"] == "S"
	L_UL = temptb1b["MC"] == "L"
	S_UL = temptb1b["MC"] == "S"

	temptb2a = temptb1a[L_CV]
	temptb2b = temptb1b[L_UL]
	temptb3a = temptb1a[S_CV]
	temptb3b = temptb1b[S_UL]

	X_L     = mainContents["Type"] == "L"
	X_LP    = mainContents["Type"] == "LP"
	X_LD    = mainContents["Type"] == "LD"
	X_UL    = mainContents["Type"] == "UL"
	X_X     = mainContents["Type"] == "X"
	temptb4 = mainContents[X_L + X_LP + X_LD + X_UL + X_X]

	X_CV     = temptb4["NH_IntErr_stat"] == "FFFFFFFFF"
	X_UL     = temptb4["NH_IntErr_stat"] == "FFFTFFTFF"
	temptb4a = mainContents[X_CV]
	temptb4b = mainContents[X_UL]

	temptb4a["NH_Int"] = np.log10(temptb4a["NH_Int"] * 1e22)
	temptb4b["NH_Int"] = np.log10(temptb4b["NH_Int"] * 1e22)

	histbin_LCV = np.linspace(int(min(temptb2a["N_HI_"])) - 1, int(max(temptb4a["NH_Int"])) + 1, 21)
	histbin_LUL = np.linspace(int(min(temptb2b["N_HI_"])) - 1, int(max(temptb4b["NH_Int"])) + 1, 21)
	histbin_SCV = np.linspace(int(min(temptb3a["N_HI_"])) - 1, int(max(temptb4a["NH_Int"])) + 1, 21)
	histbin_SUL = np.linspace(int(min(temptb3b["N_HI_"])) - 1, int(max(temptb4b["NH_Int"])) + 1, 21)

	n, bins, patches = ax[0,0].hist(temptb4a["NH_Int"], histbin_LCV, color="k", histtype="stepfilled")
	n, bins, patches = ax[0,0].hist(temptb2a["N_HI_"], histbin_LCV, color="r", histtype="step")
	n, bins, patches = ax[0,1].hist(temptb4b["NH_Int"], histbin_LUL, color="k", histtype="stepfilled")
	n, bins, patches = ax[0,1].hist(temptb2b["N_HI_"], histbin_LUL, color="r", histtype="step")
	n, bins, patches = ax[1,0].hist(temptb4a["NH_Int"], histbin_SCV, color="k", histtype="stepfilled")
	n, bins, patches = ax[1,0].hist(temptb3a["N_HI_"], histbin_SCV, color="r", histtype="step")
	n, bins, patches = ax[1,1].hist(temptb4b["NH_Int"], histbin_SUL, color="k", histtype="stepfilled")
	n, bins, patches = ax[1,1].hist(temptb3b["N_HI_"], histbin_SUL, color="r", histtype="step")

	ax[0,0].set_xlabel(r"LMC vs $N_{H_{X, Int}}$ CV ($Log_{10}$ $N_{H}$ $cm^{-2}$)")
	ax[0,0].set_ylabel(r"Number")
	ax[0,1].set_xlabel(r"LMC vs $N_{H_{X, Int}}$ UL ($Log_{10}$ $N_{H}$ $cm^{-2}$)")
	ax[0,1].set_ylabel(r"Number")
	ax[1,0].set_xlabel(r"SMC vs $N_{H_{X, Int}}$ CV ($Log_{10}$ $N_{H}$ $cm^{-2}$)")
	ax[1,0].set_ylabel(r"Number")
	ax[1,1].set_xlabel(r"SMC vs $N_{H_{X, Int}}$ UL ($Log_{10}$ $N_{H}$ $cm^{-2}$)")
	ax[1,1].set_ylabel(r"Number")

	ax[0,0].set_ylim(0, 35)
	ax[0,1].set_ylim(0, 35)
	ax[1,0].set_ylim(0, 35)
	ax[1,1].set_ylim(0, 35)

	plt.savefig(rootdir + "XRAb_N(NH)Hist_LSMC_Comparison_" + output + ".eps", bbox_inches="tight")
	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#Flux Subclass Histogram
if args.Graphs[1] == "Y" and args.Graphs[2] == "N" and args.Graphs[3] == "N":
	print "Creating Flux subclass Histograms..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, 2, figsize=(14, 4))

	L   = mainContents["Type"] == "L"
	S   = mainContents["Type"] == "S"
	SE  = mainContents["Type"] == "SE"
	SEP = mainContents["Type"] == "SEP"
	LP  = mainContents["Type"] == "LP"
	LD  = mainContents["Type"] == "LD"
	UL  = mainContents["Type"] == "UL"
	X   = mainContents["Type"] == "X"

	temptb1 = mainContents[L]
	temptb2 = mainContents[S]
	temptb3 = mainContents[SE]
	temptb4 = mainContents[SEP]
	temptb5 = mainContents[LP]
	temptb6 = mainContents[LD]
	temptb7 = mainContents[UL]
	temptb8 = mainContents[X]

	mean = np.log10(mainContents["cFx"]).mean()
	var  = np.log10(mainContents["cFx"]).var()

	histbin = np.linspace(-14, -9, 26)
	n, bins, patches = ax[0].hist([np.log10(temptb1["cFx"]), np.log10(temptb2["cFx"]), np.log10(temptb3["cFx"]), np.log10(temptb4["cFx"]), \
								np.log10(temptb5["cFx"]), np.log10(temptb6["cFx"]), np.log10(temptb7["cFx"]), np.log10(temptb8["cFx"])], \
								histbin, stacked=True, color=["k", "b", "g", "c", "c", "m", "r", "y"], histtype="stepfilled")

	ax[1].plot(np.log10(mainContents["cFx"]), gaussian(np.log10(mainContents["cFx"]), mean, var), 'k.')
	ax[1].axvline(mean, c=u"k", ls="-")

	ax[0].set_xlabel(r"$Log_{10}$ $f_{X} (erg  cm^{-2}s^{-1})$")
	ax[0].set_ylabel(r"Number")
	ax[1].set_xlabel(r"$Log_{10}$ $f_{X} (erg  cm^{-2}s^{-1})$")
	ax[1].set_ylabel(r"PDE")

	plt.savefig(rootdir + "XRAb_Flux_Subclass_Hist_" + output + ".eps", bbox_inches="tight")

	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#NHz relation no frills
if args.Graphs[1] == "Y" and args.Graphs[2] == "N" and args.Graphs[3] == "N":
	print "Creating NHvsz No Frills plot..."
	print
	print

	nz=100
	zmin=0.01
	zmax=12
	lz1=np.linspace(np.log10(zmin+1),np.log10(zmax+1),num=nz)
	z=np.power(10,lz1)-1
	z1=z+1

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, figsize=(7, 4))

	ax.set_xscale("log", nonposx="clip")
	ax.set_yscale("log", nonposy="clip")
	
	CV = mainContents["NH_IntErr_stat"] == "FFFFFFFFF"
	UL = mainContents["NH_IntErr_stat"] == "FFFTFFTFF"
		
	temptb1 = mainContents[CV]
	temptb2 = mainContents[UL]

	minUL_NInt = 1e19 * np.power(temptb2["z"], 2.4)

	Gal = temptb2["NH_Int"] * 1e22 < minUL_NInt
	Int = temptb2["NH_Int"] * 1e22 >= minUL_NInt

	temptb2a = temptb2[Gal]
	temptb2b = temptb2[Int]

	NHz_CV_plot(ax, "k", temptb1)
	NHz_UL_plot(ax, "r", temptb2a)
	NHz_UL_plot(ax, "k", temptb2b)

	ax.plot(z1, 1e19 * np.power(z1,2.4), 'c--')

	ax.set_xlim(1.0, 12.0)
	ax.set_ylim(1e19, 1e24)

	ax.xaxis.set_major_formatter(FormatStrFormatter("%1i"))
	ax.xaxis.set_minor_formatter(FormatStrFormatter("%1i"))

	ax.set_xlabel(r"$z + 1$")
	ax.set_ylabel(r"$N_{H_{X, Intrinsic}} (cm^{-2})$")

	plt.savefig(rootdir + "XRAb_NHvsz_NoFrills_" + output + ".eps", bbox_inches="tight")
	plt.close(fig)

	with open(rootdir + "Graphics_ErrorLog.txt", "a") as fl:
		fl.write("XRAb GRBs with failed fitting: FFFFFFTTF\n")
		fl.write(temptb3["GRBname"] + "\n")
		fl.write("These were not plotted\n")
		fl.write("\n")

#---------------------------------------------------------------------------------------------------------------------#
#NHz relation subclass
if args.Graphs[1] == "Y" and args.Graphs[2] == "N" and args.Graphs[3] == "N":
	print "Creating NHvsz Subclass plot..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(2, 1, figsize=(7, 8))

	ax[0].set_xscale("log", nonposx="clip")
	ax[0].set_yscale("log", nonposy="clip")
	ax[1].set_xscale("log", nonposx="clip")
	ax[1].set_yscale("log", nonposy="clip")

	NF = mainContents["NH_IntErr_stat"] != "FFFFFFTTF"

	temptb1 = mainContents[NF]

	L   = temptb1["Type"] == "L"
	S   = temptb1["Type"] == "S"
	SE  = temptb1["Type"] == "SE"
	SEP = temptb1["Type"] == "SEP"
	LP  = temptb1["Type"] == "LP"
	LD  = temptb1["Type"] == "LD"
	UL  = temptb1["Type"] == "UL"
	X   = temptb1["Type"] == "X"

	temptb2 = temptb1[L]
	temptb3 = temptb1[S]
	temptb4 = temptb1[SE]
	temptb5 = temptb1[SEP]
	temptb6 = temptb1[LP]
	temptb7 = temptb1[LD]
	temptb8 = temptb1[UL]
	temptb9 = temptb1[X]
	
	L_CV   = temptb2["NH_IntErr_stat"] == "FFFFFFFFF"
	L_UL   = temptb2["NH_IntErr_stat"] == "FFFTFFTFF"
	S_CV   = temptb3["NH_IntErr_stat"] == "FFFFFFFFF"
	S_UL   = temptb3["NH_IntErr_stat"] == "FFFTFFTFF"
	SE_CV  = temptb4["NH_IntErr_stat"] == "FFFFFFFFF"
	SE_UL  = temptb4["NH_IntErr_stat"] == "FFFTFFTFF"
	SEP_CV = temptb5["NH_IntErr_stat"] == "FFFFFFFFF"
	SEP_UL = temptb5["NH_IntErr_stat"] == "FFFTFFTFF"
	LP_CV  = temptb6["NH_IntErr_stat"] == "FFFFFFFFF"
	LP_UL  = temptb6["NH_IntErr_stat"] == "FFFTFFTFF"
	LD_CV  = temptb7["NH_IntErr_stat"] == "FFFFFFFFF"
	LD_UL  = temptb7["NH_IntErr_stat"] == "FFFTFFTFF"
	UL_CV  = temptb8["NH_IntErr_stat"] == "FFFFFFFFF"
	UL_UL  = temptb8["NH_IntErr_stat"] == "FFFTFFTFF"
	X_CV   = temptb9["NH_IntErr_stat"] == "FFFFFFFFF"
	X_UL   = temptb9["NH_IntErr_stat"] == "FFFTFFTFF"
	
	temptb2a = temptb2[L_CV]
	temptb2b = temptb2[L_UL]
	temptb3a = temptb3[S_CV]
	temptb3b = temptb3[S_UL]
	temptb4a = temptb4[SE_CV]
	temptb4b = temptb4[SE_UL]
	temptb5a = temptb5[SEP_CV]
	temptb5b = temptb5[SEP_UL]
	temptb6a = temptb6[LP_CV]
	temptb6b = temptb6[LP_UL]
	temptb7a = temptb7[LD_CV]
	temptb7b = temptb7[LD_UL]
	temptb8a = temptb8[UL_CV]
	temptb8b = temptb8[UL_UL]
	temptb9a = temptb9[X_CV]
	temptb9b = temptb9[X_UL]

	NHz_CV_plot(ax[0], "k", temptb2a)
	NHz_UL_plot(ax[0], "k", temptb2b)
	NHz_CV_plot(ax[0], "b", temptb3a)
	NHz_UL_plot(ax[0], "b", temptb3b)
	NHz_CV_plot(ax[0], "g", temptb4a)
	NHz_UL_plot(ax[0], "g", temptb4b)
	NHz_CV_plot(ax[0], "c", temptb5a)
	NHz_UL_plot(ax[0], "c", temptb5b)
	NHz_CV_plot(ax[0], "c", temptb6a)
	NHz_UL_plot(ax[0], "c", temptb6b)
	NHz_CV_plot(ax[0], "m", temptb7a)
	NHz_UL_plot(ax[0], "m", temptb7b)
	NHz_CV_plot(ax[0], "r", temptb8a)
	NHz_UL_plot(ax[0], "r", temptb8b)
	NHz_CV_plot(ax[0], "y", temptb9a)
	NHz_UL_plot(ax[0], "y", temptb9b)

	NHz_CV_plot(ax[1], "b", temptb3a)
	NHz_UL_plot(ax[1], "b", temptb3b)
	NHz_CV_plot(ax[1], "g", temptb4a)
	NHz_UL_plot(ax[1], "g", temptb4b)
	NHz_CV_plot(ax[1], "c", temptb5a)
	NHz_UL_plot(ax[1], "c", temptb5b)
	NHz_CV_plot(ax[1], "c", temptb6a)
	NHz_UL_plot(ax[1], "c", temptb6b)
	NHz_CV_plot(ax[1], "m", temptb7a)
	NHz_UL_plot(ax[1], "m", temptb7b)
	NHz_CV_plot(ax[1], "r", temptb8a)
	NHz_UL_plot(ax[1], "r", temptb8b)
	NHz_CV_plot(ax[1], "y", temptb9a)
	NHz_UL_plot(ax[1], "y", temptb9b)

	ax[0].set_xlim(1.0, 12.0)
	ax[0].set_ylim(1e19, 1e24)
	ax[1].set_xlim(1.0, 12.0)
	ax[1].set_ylim(1e19, 1e24)

	ax[0].xaxis.set_major_formatter(FormatStrFormatter("%1i"))
	ax[0].xaxis.set_minor_formatter(FormatStrFormatter("%1i"))
	ax[1].xaxis.set_major_formatter(FormatStrFormatter("%1i"))
	ax[1].xaxis.set_minor_formatter(FormatStrFormatter("%1i"))

	ax[0].set_xlabel(r"$z + 1$")
	ax[0].set_ylabel(r"$N_{H_{X, Intrinsic}} (cm^{-2})$")
	ax[1].set_xlabel(r"$z + 1$")
	ax[1].set_ylabel(r"$N_{H_{X, Intrinsic}} (cm^{-2})$")

	plt.savefig(rootdir + "XRAb_NHvsz_SubClass_" + output + ".eps", bbox_inches="tight")
	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#NHz relation metallicity
if args.Graphs[1] == "Y" and args.Graphs[2] == "Y":
	print "Creating NHvsz Metallicity plot..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(2, 1, figsize=(7, 8))

	ax[0].set_xscale("log", nonposx="clip")
	ax[0].set_yscale("log", nonposy="clip")
	ax[1].set_xscale("log", nonposx="clip")
	ax[1].set_yscale("log", nonposy="clip")

	NF_X = mainContents["NH_IntErr_stat"] != "FFFFFFTTF"
	NF_Z = ZContents["NH_IntErr_stat"] != "FFFFFFTTF"
	F_Z  = ZContents["NH_IntErr_stat"] == "FFFFFFTTF"

	temptb1 = mainContents[NF_X]
	temptb2 = ZContents[NF_Z]
	temptb3 = ZContents[F_Z]

	X_CV = temptb1["NH_IntErr_stat"] == "FFFFFFFFF"
	X_UL = temptb1["NH_IntErr_stat"] == "FFFTFFTFF"
	Z_CV = temptb2["NH_IntErr_stat"] == "FFFFFFFFF"
	Z_UL = temptb2["NH_IntErr_stat"] == "FFFTFFTFF"

	temptb1a = temptb1[X_CV]
	temptb1b = temptb1[X_UL]
	temptb2a = temptb2[Z_CV]
	temptb2b = temptb2[Z_UL]

	CCVa = temptb2a["Author"] == "Ccv"
	CUPa = temptb2a["Author"] == "Cup"
	CCVb = temptb2b["Author"] == "Ccv"
	CUPb = temptb2b["Author"] == "Cup"

	temptb3a = temptb2a[CCVa] 
	temptb3b = temptb2b[CCVb]
	temptb4a = temptb2a[CUPa]
	temptb4b = temptb2b[CUPb]

	NHz_CV_plot(ax[0], "k", temptb1a)
	NHz_UL_plot(ax[0], "k", temptb1b)	
	NHz_CV_plot(ax[0], "r", temptb3a)
	NHz_UL_plot(ax[0], "r", temptb3b)		
	NHz_CV_plot(ax[0], "c", temptb4a)
	NHz_UL_plot(ax[0], "c", temptb4b)		

	NHz_CV_plot(ax[1], "r", temptb3a)
	NHz_UL_plot(ax[1], "r", temptb3b)		
	NHz_CV_plot(ax[1], "c", temptb4a)
	NHz_UL_plot(ax[1], "c", temptb4b)	

	ax[0].set_xlim(1.0, 12.0)
	ax[0].set_ylim(1e19, 1e24)
	ax[1].set_xlim(1.0, 12.0)
	ax[1].set_ylim(1e19, 1e24)

	ax[0].xaxis.set_major_formatter(FormatStrFormatter("%1i"))
	ax[0].xaxis.set_minor_formatter(FormatStrFormatter("%1i"))
	ax[1].xaxis.set_major_formatter(FormatStrFormatter("%1i"))
	ax[1].xaxis.set_minor_formatter(FormatStrFormatter("%1i"))

	ax[0].set_xlabel(r"$z + 1$")
	ax[0].set_ylabel(r"$N_{H_{X, Intrinsic}} (cm^{-2})$")
	ax[1].set_xlabel(r"$z + 1$")
	ax[1].set_ylabel(r"$N_{H_{X, Intrinsic}} (cm^{-2})$")

	plt.savefig(rootdir + "XRAb_NHvsz_Metallicity_" + output + ".eps", bbox_inches="tight")
	plt.close(fig)

	with open(rootdir + "Graphics_ErrorLog.txt", "a") as fl:
		fl.write("Z GRBs with failed fitting: FFFFFFTTF\n")
		fl.write(temptb3["GRBname"] + "\n")
		fl.write("These were not plotted\n")
		fl.write("\n")

#---------------------------------------------------------------------------------------------------------------------#
#NHz correlation tests
if args.Graphs[1] == "Y" or args.Graphs[2] == "Y":
	print "Creating NHvsz Correlation test plot..."
	print
	print

	fig = plt.figure(1)
	fig, ax = plt.subplots(2, 2, figsize=(14, 8))

	ax[0,0].set_xscale("log", nonposx="clip")
	ax[0,0].set_yscale("log", nonposy="clip")
	ax[0,1].set_xscale("log", nonposx="clip")
	ax[1,0].set_xscale("log", nonposx="clip")
	ax[1,0].set_yscale("log", nonposx="clip")
	ax[1,1].set_xscale("log", nonposx="clip")
	ax[1,1].set_yscale("log", nonposx="clip")

	NHI_CV = mainContents["NH_IntErr_stat"] == "FFFFFFFFF"
	NHI_UL = mainContents["NH_IntErr_stat"] == "FFFTFFTFF"

	temptb1a = mainContents[NHI_CV]
	temptb1b = mainContents[NHI_UL]

	NHIG_CV_plot(ax[0,0], "k", temptb1a)
	NHIG_UL_plot(ax[0,0], "r", temptb1b)

	NHIP_CV_plot(ax[0,1], "k", temptb1a)
	NHIP_UL_plot(ax[0,1], "r", temptb1b)

	NHIF_CV_plot(ax[1,0], "k", temptb1a)
	NHIF_UL_plot(ax[1,0], "r", temptb1b)

	NHIp_CV_plot(ax[1,1], "k", temptb1a)
	NHIp_UL_plot(ax[1,1], "r", temptb1b)

	ax[0,0].set_xlim(1e19, 1e24)
	ax[0,0].set_ylim(5e19, 1e22)
	ax[0,1].set_xlim(1e19, 1e24)
	ax[0,1].set_ylim(1.0, 3.5)
	ax[1,0].set_xlim(1e19, 1e24)
	ax[1,0].set_ylim(1e-13, 1e-10)
	ax[1,1].set_xlim(1e19, 1e24)
	ax[1,1].set_ylim(3e3, 3e5)

	ax[0,0].set_xlabel(r"$N_{H_{X, Intrinsic}} (cm^{-2})$")
	ax[0,1].set_xlabel(r"$N_{H_{X, Intrinsic}} (cm^{-2})$")
	ax[1,0].set_xlabel(r"$N_{H_{X, Intrinsic}} (cm^{-2})$")
	ax[1,1].set_xlabel(r"$N_{H_{X, Intrinsic}} (cm^{-2})$")

	ax[0,0].set_ylabel(r"$N_{H_{X, Gal}} (cm^{-2})$")
	ax[0,1].set_ylabel(r"Spectral Index, $\Gamma$")
	ax[1,0].set_ylabel(r"$f_{0.3-10.0 keV} (erg cm^{-2} s^{-1})$")
	ax[1,1].set_ylabel(r"Mean Photon Arrival Time")

	if args.Graphs[1] == "Y" and args.Graphs[2] == "N":
		plt.savefig(rootdir + "XRAb_NHvsz_Correlations_" + output + ".eps", bbox_inches="tight")
	if args.Graphs[1] == "N" and args.Graphs[2] == "Y":
		plt.savefig(rootdir + "Z_NHvsz_Correlations_" + output + ".eps", bbox_inches="tight")

	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#NHz binned distribution
if args.Graphs[1] == "Y" and args.Graphs[2] == "N" and args.Graphs[3] == "N":
	print "Creating Binned NHz plot..."
	print
	print

	min_NInt = 1e19 * np.power(mainContents["z"], 2.4)

	Real   = mainContents["NH_Int"] * 1e22 > min_NInt
	temptb = mainContents[Real]

	mask025 = temptb["z"] <= 0.25
	mask050 = temptb["z"] <= 0.5
	mask075 = temptb["z"] <= 0.75
	mask100 = temptb["z"] <= 1 
	mask125 = temptb["z"] <= 1.25
	mask150 = temptb["z"] <= 1.5
	mask175 = temptb["z"] <= 1.75
	mask200 = temptb["z"] <= 2 
	mask225 = temptb["z"] <= 2.25
	mask250 = temptb["z"] <= 2.50
	mask275 = temptb["z"] <= 2.75
	mask325 = temptb["z"] <= 3.25
	mask400 = temptb["z"] <= 4
	mask500 = temptb["z"] <= 5
	mask700 = temptb["z"] <= 7
	mask1000 = temptb["z"] <= 10

	temptb000025  = temptb[mask025]
	temptb025050  = temptb[mask050 - mask025]
	temptb050075  = temptb[mask075 - mask050]
	temptb075100  = temptb[mask100 - mask075]
	temptb100125  = temptb[mask125 - mask100]
	temptb125150  = temptb[mask150 - mask125]
	temptb150175  = temptb[mask175 - mask150]
	temptb175200  = temptb[mask200 - mask175]
	temptb200225  = temptb[mask225 - mask200]
	temptb225250  = temptb[mask250 - mask225]
	temptb250275  = temptb[mask275 - mask250]
	temptb275325  = temptb[mask325 - mask275]
	temptb325400  = temptb[mask400 - mask325]
	temptb400500  = temptb[mask500 - mask400]
	temptb500700  = temptb[mask700 - mask500]
	temptb7001000 = temptb[mask1000 - mask700]

	l = [len(temptb000025["z"]), len(temptb025050["z"]), len(temptb050075["z"]), len(temptb075100["z"]), len(temptb100125["z"]), \
		len(temptb125150["z"]), len(temptb150175["z"]), len(temptb175200["z"]), len(temptb200225["z"]), len(temptb225250["z"]), \
		len(temptb250275["z"]), len(temptb275325["z"]), len(temptb325400["z"]), len(temptb400500["z"]), len(temptb500700["z"]), \
		len(temptb7001000["z"])]

	x = [1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 4.000, 4.625, 5.500, 7.000, 9.500]

	y = [np.median(temptb000025["NH_Int"]) * 1e22, np.median(temptb025050["NH_Int"]) * 1e22, np.median(temptb050075["NH_Int"]) * 1e22, \
		np.median(temptb075100["NH_Int"]) * 1e22, np.median(temptb100125["NH_Int"]) * 1e22, np.median(temptb125150["NH_Int"]) * 1e22, \
		np.median(temptb150175["NH_Int"]) * 1e22, np.median(temptb175200["NH_Int"]) * 1e22, np.median(temptb200225["NH_Int"]) * 1e22, \
		np.median(temptb225250["NH_Int"]) * 1e22, np.median(temptb250275["NH_Int"]) * 1e22, np.median(temptb275325["NH_Int"]) * 1e22, \
		np.median(temptb325400["NH_Int"]) * 1e22, np.median(temptb400500["NH_Int"]) * 1e22, np.median(temptb500700["NH_Int"]) * 1e22, \
		np.median(temptb7001000["NH_Int"]) * 1e22]

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, figsize=(7, 4))

	ax.set_yscale("log", nonposy="clip")
	ax.set_xscale("log", nonposy="clip")
	ax.set_xlabel(r"$z + 1$")
	ax.set_ylabel(r"median $N_{H_{X, Intrinsic}}$ $(cm^{-2})$")

	NH_Bin_plot(ax, 1.125, 0.125, "k", temptb000025)
	NH_Bin_plot(ax, 1.375, 0.125, "k", temptb025050)
	NH_Bin_plot(ax, 1.625, 0.125, "k", temptb050075)
	NH_Bin_plot(ax, 1.875, 0.125, "k", temptb075100)
	NH_Bin_plot(ax, 2.125, 0.125, "k", temptb100125)
	NH_Bin_plot(ax, 2.375, 0.125, "k", temptb125150)
	NH_Bin_plot(ax, 2.625, 0.125, "k", temptb150175)
	NH_Bin_plot(ax, 2.875, 0.125, "k", temptb175200)
	NH_Bin_plot(ax, 3.125, 0.125, "k", temptb200225)
	NH_Bin_plot(ax, 3.375, 0.125, "k", temptb225250)
	NH_Bin_plot(ax, 3.625, 0.125, "k", temptb250275)
	NH_Bin_plot(ax, 4.000, 0.250, "k", temptb275325)
	NH_Bin_plot(ax, 4.625, 0.375, "k", temptb325400)
	NH_Bin_plot(ax, 5.500, 0.500, "k", temptb400500)
	NH_Bin_plot(ax, 7.000, 1.000, "k", temptb500700)
	NH_Bin_plot(ax, 9.500, 1.500, "k", temptb7001000)

	for i, txt in enumerate(l):
   		ax.annotate(txt, xy=(x[i],y[i]), xytext=(-3,5), textcoords="offset points", color="r", size="x-small")

	ax.set_xlim(1.0, 12.0)
	ax.xaxis.set_major_formatter(FormatStrFormatter("%1i"))
	ax.xaxis.set_minor_formatter(FormatStrFormatter("%1i"))

	plt.savefig(rootdir + "N(z)Bin_" + output + ".eps", bbox_inches="tight")
	plt.close(fig)

#---------------------------------------------------------------------------------------------------------------------#
#NHz relation with IGM models
if args.Graphs[1] == "Y" and args.Graphs[2] == "N" and args.Graphs[3] == "N":
	print "Creating NHvsz IGM Models plot..."
	print
	print

	subFITs      = fits.open(rootdir + "XRAb_FITSTable_IGMModel_1Z_0imind0inind_zlist_April2015_Run6_allpriors.fits") #+ output + ".fits")
	subContents  = subFITs[1].data

	fig = plt.figure(1)
	fig, ax = plt.subplots(1, figsize=(7, 4))

	ax.set_xscale("log", nonposx="clip")
	ax.set_yscale("log", nonposy="clip")
	
	CV = mainContents["NH_IntErr_stat"] == "FFFFFFFFF"
	UL = mainContents["NH_IntErr_stat"] == "FFFTFFTFF"
		
	temptb1 = mainContents[CV]
	temptb2 = mainContents[UL]

	minUL_NInt = 1e19 * np.power(temptb2["z"], 2.4)

	Gal = temptb2["NH_Int"] * 1e22 < minUL_NInt
	Int = temptb2["NH_Int"] * 1e22 >= minUL_NInt

	temptb2a = temptb2[Gal]
	temptb2b = temptb2[Int]

	NHz_CV_plot(ax, "k", temptb1)
	NHz_UL_plot(ax, "r", temptb2a)
	NHz_UL_plot(ax, "k", temptb2b)

	ax.errorbar(6.5, 5.2406 * 1e22,  
		yerr=[[12.644 * 1e22], [5.241 * 1e22]], 
		ecolor='g', mec='g', mfc='g', capsize=0, marker='D', markersize=3, ls="None", alpha=0.5)
	ax.errorbar(6.5, 1.46308 * 1e22,  
		yerr=[[0.687 * 1e22], [0.615 * 1e22]], 
		ecolor='c', mec='c', mfc='c', capsize=0, marker='D', markersize=3, ls="None")

#	ax.plot(subContents["z"] + 1, subContents["NH_Lyf"] * 1e22, 'b-')
	ax.plot(subContents["z"] + 1, subContents["NH_TCold"] * 1e22, 'k-')
	ax.plot(subContents["z"] + 1, subContents["NH_TWarm1"] * 1e22, 'm--')
	ax.plot(subContents["z"] + 1, subContents["NH_TWarm2"] * 1e22, 'm-.')
	ax.plot(subContents["z"] + 1, subContents["NH_TWarm3"] * 1e22, 'm:')

	ax.text(7.7, 2.5e20, "-   Cold IGM", color='k', fontsize=8) 
	ax.text(7.7, 1.6e20, "--  $10^{4}$ K IGM", color='m', fontsize=8) 
	ax.text(7.7, 1.0e20, "-.- $10^{6}$ K IGM", color='m', fontsize=8) 
	ax.text(7.7, 0.6e20, "..  $10^{7}$ K IGM", color='m', fontsize=8) 

	ax.set_xlim(1.0, 12.0)
	ax.set_ylim(1e19, 1e24)

	ax.xaxis.set_major_formatter(FormatStrFormatter("%1i"))
	ax.xaxis.set_minor_formatter(FormatStrFormatter("%1i"))

	ax.set_xlabel(r"$z + 1$")
	ax.set_ylabel(r"$N_{H_{X, Intrinsic}} (cm^{-2})$")

	plt.savefig(rootdir + "XRAb_NHvsz_IGM_1Z_0imind0inind_" + output + ".eps", bbox_inches="tight")
	plt.close(fig)

	subFITs.close()

#---------------------------------------------------------------------------------------------------------------------#
#NHz relation with IGM subtracted
if args.Graphs[1] == "Y" and args.Graphs[2] == "N" and args.Graphs[3] == "N":
	print "Creating NHvsz IGM subtracted plot..."
	print
	print

	subFITs      = fits.open(rootdir + "XRAb_FITSTable_ModelSub_1Z_0imind0inind_" + output + ".fits")
	subContents  = subFITs[1].data

	min_NInt = 1e19 * np.power(mainContents["z"], 2.4)
	Real = mainContents["NH_Int"] * 1e22 >= min_NInt

	temptb = mainContents[Real]

	CV = temptb["NH_IntErr_stat"] == "FFFFFFFFF"
	UL = temptb["NH_IntErr_stat"] == "FFFTFFTFF"

	tbcold_cv  = temptb[CV]
	tbcold_ul  = temptb[UL]
	tbwarm1_cv = temptb[CV]
	tbwarm1_ul = temptb[UL]
	tbwarm2_cv = temptb[CV]
	tbwarm2_ul = temptb[UL]
	tbwarm3_cv = temptb[CV]
	tbwarm3_ul = temptb[UL]

	for i in range(len(tbcold_cv["z"])):
		for j in range(len(subContents["z"])):
			if tbcold_cv["z"][i] == subContents["z"][j]: 
				tbcold_cv["NH_Int"][i]  = tbcold_cv["NH_Int"][i]  - subContents["NH_TCold"][j]
				tbwarm1_cv["NH_Int"][i] = tbwarm1_cv["NH_Int"][i] - subContents["NH_TWarm1"][j]	
				tbwarm2_cv["NH_Int"][i] = tbwarm2_cv["NH_Int"][i] - subContents["NH_TWarm2"][j]
				tbwarm3_cv["NH_Int"][i] = tbwarm3_cv["NH_Int"][i] - subContents["NH_TWarm3"][j]

		if tbcold_cv["NH_Int"][i] < 0:
			with open(rootdir + "XRAb_NHvsz_Models_ErrorLog.txt", "a") as fel:
				fel.write("TCold: " + tbcold_cv["GRBname"][i] + " now has negative NH: " + str(tbcold_cv["NH_Int"][i]) + "\n")
		if tbwarm1_cv["NH_Int"][i] < 0:
			with open(rootdir + "XRAb_NHvsz_Models_ErrorLog.txt", "a") as fel:
				fel.write("TWarm1: " + tbwarm1_cv["GRBname"][i] + " now has negative NH: " + str(tbwarm1_cv["NH_Int"][i]) + "\n")
		if tbwarm2_cv["NH_Int"][i] < 0:
			with open(rootdir + "XRAb_NHvsz_Models_ErrorLog.txt", "a") as fel:
				fel.write("TWarm2: " + tbwarm2_cv["GRBname"][i] + " now has negative NH: " + str(tbwarm2_cv["NH_Int"][i]) + "\n")
		if tbwarm3_cv["NH_Int"][i] < 0:
			with open(rootdir + "XRAb_NHvsz_Models_ErrorLog.txt", "a") as fel:
				fel.write("TWarm3: " + tbwarm3_cv["GRBname"][i] + " now has negative NH: " + str(tbwarm3_cv["NH_Int"][i]) + "\n")

	for i in range(len(tbcold_ul["z"])):
		for j in range(len(subContents["z"])):
			if tbcold_ul["z"][i] == subContents["z"][j]: 
				tbcold_ul["NH_Int"][i]  = tbcold_ul["NH_Int"][i]  - subContents["NH_TCold"][j]
				tbwarm1_ul["NH_Int"][i] = tbwarm1_ul["NH_Int"][i] - subContents["NH_TWarm1"][j]	
				tbwarm2_ul["NH_Int"][i] = tbwarm2_ul["NH_Int"][i] - subContents["NH_TWarm2"][j]
				tbwarm3_ul["NH_Int"][i] = tbwarm3_ul["NH_Int"][i] - subContents["NH_TWarm3"][j]

		if tbcold_ul["NH_Int"][i] < 0:
			with open(rootdir + "XRAb_NHvsz_Models_ErrorLog.txt", "a") as fel:
				fel.write("TCold: " + tbcold_ul["GRBname"][i] + " now has negative NH: " + str(tbcold_ul["NH_Int"][i]) + "\n")
		if tbwarm1_ul["NH_Int"][i] < 0:
			with open(rootdir + "XRAb_NHvsz_Models_ErrorLog.txt", "a") as fel:
				fel.write("TWarm1: " + tbwarm1_ul["GRBname"][i] + " now has negative NH: " + str(tbwarm1_ul["NH_Int"][i]) + "\n")
		if tbwarm2_ul["NH_Int"][i] < 0:
			with open(rootdir + "XRAb_NHvsz_Models_ErrorLog.txt", "a") as fel:
				fel.write("TWarm2: " + tbwarm2_ul["GRBname"][i] + " now has negative NH: " + str(tbwarm2_ul["NH_Int"][i]) + "\n")
		if tbwarm3_ul["NH_Int"][i] < 0:
			with open(rootdir + "XRAb_NHvsz_Models_ErrorLog.txt", "a") as fel:
				fel.write("TWarm3: " + tbwarm3_ul["GRBname"][i] + " now has negative NH: " + str(tbwarm3_ul["NH_Int"][i]) + "\n")

	fig = plt.figure(1)
	fig, ax = plt.subplots(4, 2, figsize=(13, 11))

	d = .015
	plt.subplots_adjust(hspace=0.1)

	#Top left Plot
	ax[0,0].set_xscale("log", nonposx="clip")
	ax[0,0].set_yscale("log", nonposy="clip")
	ax[1,0].set_xscale("log", nonposx="clip")
	ax[1,0].set_yscale("symlog", subsy=[2,3,4,5,6,7,8,9])

	for i in range(len(tbcold_cv["z"])):
		if tbcold_cv["NH_Int"][i] > 0:
			ax[0,0].errorbar(tbcold_cv["z"][i] + 1, tbcold_cv["NH_Int"][i] * 1e22, yerr=[[tbcold_cv["NH_IntErr_low"][i] *1e22], [tbcold_cv["NH_IntErr_hi"][i] *1e22]],
					capsize=0, marker='.', color='k', ecolor='k')
		else:
			ax[1,0].errorbar(tbcold_cv["z"][i] + 1, tbcold_cv["NH_Int"][i] * 1e22, yerr=[[tbcold_cv["NH_IntErr_low"][i] *1e22], [tbcold_cv["NH_IntErr_hi"][i] *1e22]],
					capsize=0, marker='.', color='k', ecolor='k')

	for i in range(len(tbcold_ul["z"])):
		if tbcold_ul["NH_Int"][i] > 0:
			ax[0,0].errorbar(tbcold_ul["z"][i] + 1, tbcold_ul["NH_Int"][i] * 1e22, yerr=[[tbcold_ul["NH_IntErr_low"][i] *1e22], [tbcold_ul["NH_IntErr_hi"][i] *1e22]],
					capsize=2, color='r', ecolor='r', lolims=True)
		else:
			ax[1,0].errorbar(tbcold_ul["z"][i] + 1, tbcold_ul["NH_Int"][i] * 1e22, yerr=[[tbcold_ul["NH_IntErr_low"][i] *1e22], [tbcold_ul["NH_IntErr_hi"][i] *1e22]],
					capsize=2, color='r', ecolor='r', lolims=True)

	ax[0,0].set_xlim(1.0, 12.0)
	ax[0,0].set_ylim(1e19, 1e24)
	ax[1,0].set_xlim(1.0, 12.0)
	ax[1,0].set_ylim(-1e24, -1e19)

	ax[0,0].spines['bottom'].set_visible(False)
	ax[0,0].xaxis.tick_top()
	ax[0,0].tick_params(labeltop='off')
	ax[1,0].spines['top'].set_visible(False)
	ax[1,0].xaxis.tick_bottom()

	ax[1,0].xaxis.set_major_formatter(FormatStrFormatter("%1i"))
	ax[1,0].xaxis.set_minor_formatter(FormatStrFormatter("%1i"))

	kwargs = dict(transform=ax[0,0].transAxes, color='k', clip_on=False)
	ax[0,0].plot((-d,d),(-d,+d), **kwargs)
	ax[0,0].plot((1-d,1+d),(-d,+d), **kwargs)

	kwargs.update(transform=ax[1,0].transAxes)
	ax[1,0].plot((1-d,1+d),(1-d,1+d), **kwargs)
	ax[1,0].plot((-d,d),(1-d,1+d), **kwargs)

	fig.text(0.075, 0.7, r"$N_{H_{X, Intrinsic}}$ $(cm^{-2})$ $T=Cold$ Model", ha='center', va='center', rotation='vertical')

	#Top Right Plot

	ax[0,1].set_xscale("log", nonposx="clip")
	ax[0,1].set_yscale("log", nonposy="clip")
	ax[1,1].set_xscale("log", nonposx="clip")
	ax[1,1].set_yscale("symlog", subsy=[2,3,4,5,6,7,8,9])

	for i in range(len(tbwarm1_cv["z"])):
		if tbwarm1_cv["NH_Int"][i] > 0:
			ax[0,1].errorbar(tbwarm1_cv["z"][i] + 1, tbwarm1_cv["NH_Int"][i] * 1e22, yerr=[[tbwarm1_cv["NH_IntErr_low"][i] *1e22], [tbwarm1_cv["NH_IntErr_hi"][i] *1e22]],
					capsize=0, marker='.', color='k', ecolor='k')
		else:
			ax[1,1].errorbar(tbwarm1_cv["z"][i] + 1, tbwarm1_cv["NH_Int"][i] * 1e22, yerr=[[tbwarm1_cv["NH_IntErr_low"][i] *1e22], [tbwarm1_cv["NH_IntErr_hi"][i] *1e22]],
					capsize=0, marker='.', color='k', ecolor='k')

	for i in range(len(tbwarm1_ul["z"])):
		if tbwarm1_ul["NH_Int"][i] > 0:
			ax[0,1].errorbar(tbwarm1_ul["z"][i] + 1, tbwarm1_ul["NH_Int"][i] * 1e22, yerr=[[tbwarm1_ul["NH_IntErr_low"][i] *1e22], [tbwarm1_ul["NH_IntErr_hi"][i] *1e22]],
					capsize=2, color='r', ecolor='r', lolims=True)
		else:
			ax[1,1].errorbar(tbwarm1_ul["z"][i] + 1, tbwarm1_ul["NH_Int"][i] * 1e22, yerr=[[tbwarm1_ul["NH_IntErr_low"][i] *1e22], [tbwarm1_ul["NH_IntErr_hi"][i] *1e22]],
					capsize=2, color='r', ecolor='r', lolims=True)

	ax[0,1].set_xlim(1.0, 12.0)
	ax[0,1].set_ylim(1e19, 1e24)
	ax[1,1].set_xlim(1.0, 12.0)
	ax[1,1].set_ylim(-1e24, -1e19)

	ax[0,1].spines['bottom'].set_visible(False)
	ax[0,1].xaxis.tick_top()
	ax[0,1].tick_params(labeltop='off')
	ax[1,1].spines['top'].set_visible(False)
	ax[1,1].xaxis.tick_bottom()

	ax[1,1].xaxis.set_major_formatter(FormatStrFormatter("%1i"))
	ax[1,1].xaxis.set_minor_formatter(FormatStrFormatter("%1i"))

	kwargs.update(transform=ax[0,1].transAxes)
	ax[0,1].plot((-d,d),(-d,+d), **kwargs)
	ax[0,1].plot((1-d,1+d),(-d,+d), **kwargs)

	kwargs.update(transform=ax[1,1].transAxes)
	ax[1,1].plot((1-d,1+d),(1-d,1+d), **kwargs)
	ax[1,1].plot((-d,d),(1-d,1+d), **kwargs)

	fig.text(0.5, 0.7, r"$N_{H_{X, Intrinsic}}$ $(cm^{-2})$ $T=10^{4} K$ Model", ha='center', va='center', rotation='vertical')

	#Bottom Left Plot

	ax[2,0].set_xscale("log", nonposx="clip")
	ax[2,0].set_yscale("log", nonposy="clip")
	ax[3,0].set_xscale("log", nonposx="clip")
	ax[3,0].set_yscale("symlog", subsy=[2,3,4,5,6,7,8,9])

	for i in range(len(tbwarm2_cv["z"])):
		if tbwarm2_cv["NH_Int"][i] > 0:
			ax[2,0].errorbar(tbwarm2_cv["z"][i] + 1, tbwarm2_cv["NH_Int"][i] * 1e22, yerr=[[tbwarm2_cv["NH_IntErr_low"][i] *1e22], [tbwarm2_cv["NH_IntErr_hi"][i] *1e22]],
					capsize=0, marker='.', color='k', ecolor='k')
		else:
			ax[3,0].errorbar(tbwarm2_cv["z"][i] + 1, tbwarm2_cv["NH_Int"][i] * 1e22, yerr=[[tbwarm2_cv["NH_IntErr_low"][i] *1e22], [tbwarm2_cv["NH_IntErr_hi"][i] *1e22]],
					capsize=0, marker='.', color='k', ecolor='k')

	for i in range(len(tbwarm2_ul["z"])):
		if tbwarm2_ul["NH_Int"][i] > 0:
			ax[2,0].errorbar(tbwarm2_ul["z"][i] + 1, tbwarm2_ul["NH_Int"][i] * 1e22, yerr=[[tbwarm2_ul["NH_IntErr_low"][i] *1e22], [tbwarm2_ul["NH_IntErr_hi"][i] *1e22]],
					capsize=2, color='r', ecolor='r', lolims=True)
		else:
			ax[3,0].errorbar(tbwarm2_ul["z"][i] + 1, tbwarm2_ul["NH_Int"][i] * 1e22, yerr=[[tbwarm2_ul["NH_IntErr_low"][i] *1e22], [tbwarm2_ul["NH_IntErr_hi"][i] *1e22]],
					capsize=2, color='r', ecolor='r', lolims=True)

	ax[2,0].set_xlim(1.0, 12.0)
	ax[2,0].set_ylim(1e19, 1e24)
	ax[3,0].set_xlim(1.0, 12.0)
	ax[3,0].set_ylim(-1e24, -1e19)

	ax[2,0].spines['bottom'].set_visible(False)
	ax[2,0].xaxis.tick_top()
	ax[2,0].tick_params(labeltop='off')
	ax[3,0].spines['top'].set_visible(False)
	ax[3,0].xaxis.tick_bottom()

	ax[3,0].xaxis.set_major_formatter(FormatStrFormatter("%1i"))
	ax[3,0].xaxis.set_minor_formatter(FormatStrFormatter("%1i"))

	kwargs.update(transform=ax[2,0].transAxes)
	ax[2,0].plot((-d,d),(-d,+d), **kwargs)
	ax[2,0].plot((1-d,1+d),(-d,+d), **kwargs)

	kwargs.update(transform=ax[3,0].transAxes)
	ax[3,0].plot((1-d,1+d),(1-d,1+d), **kwargs)
	ax[3,0].plot((-d,d),(1-d,1+d), **kwargs)

	ax[3,0].set_xlabel(r"$z + 1$")
	fig.text(0.075, 0.3, r"$N_{H_{X, Intrinsic}}$ $(cm^{-2})$ $T=10^{6} K$ Model", ha='center', va='center', rotation='vertical')

	#Bottom Right Plot 

	ax[2,1].set_xscale("log", nonposx="clip")
	ax[2,1].set_yscale("log", nonposy="clip")
	ax[3,1].set_xscale("log", nonposx="clip")
	ax[3,1].set_yscale("symlog", subsy=[2,3,4,5,6,7,8,9])

	for i in range(len(tbwarm3_cv["z"])):
		if tbwarm3_cv["NH_Int"][i] > 0:
			ax[2,1].errorbar(tbwarm3_cv["z"][i] + 1, tbwarm3_cv["NH_Int"][i] * 1e22, yerr=[[tbwarm3_cv["NH_IntErr_low"][i] *1e22], [tbwarm3_cv["NH_IntErr_hi"][i] *1e22]],
					capsize=0, marker='.', color='k', ecolor='k')
		else:
			ax[3,1].errorbar(tbwarm3_cv["z"][i] + 1, tbwarm3_cv["NH_Int"][i] * 1e22, yerr=[[tbwarm3_cv["NH_IntErr_low"][i] *1e22], [tbwarm3_cv["NH_IntErr_hi"][i] *1e22]],
					capsize=0, marker='.', color='k', ecolor='k')

	for i in range(len(tbwarm3_ul["z"])):
		if tbwarm3_ul["NH_Int"][i] > 0:
			ax[2,1].errorbar(tbwarm3_ul["z"][i] + 1, tbwarm3_ul["NH_Int"][i] * 1e22, yerr=[[tbwarm3_ul["NH_IntErr_low"][i] *1e22], [tbwarm3_ul["NH_IntErr_hi"][i] *1e22]],
					capsize=2, color='r', ecolor='r', lolims=True)
		else:
			ax[3,1].errorbar(tbwarm3_ul["z"][i] + 1, tbwarm3_ul["NH_Int"][i] * 1e22, yerr=[[tbwarm3_ul["NH_IntErr_low"][i] *1e22], [tbwarm3_ul["NH_IntErr_hi"][i] *1e22]],
					capsize=2, color='r', ecolor='r', lolims=True)

	ax[2,1].set_xlim(1.0, 12.0)
	ax[2,1].set_ylim(1e19, 1e24)
	ax[3,1].set_xlim(1.0, 12.0)
	ax[3,1].set_ylim(-1e24, -1e19)

	ax[2,1].spines['bottom'].set_visible(False)
	ax[2,1].xaxis.tick_top()
	ax[2,1].tick_params(labeltop='off')
	ax[3,1].spines['top'].set_visible(False)
	ax[3,1].xaxis.tick_bottom()

	ax[3,1].xaxis.set_major_formatter(FormatStrFormatter("%1i"))
	ax[3,1].xaxis.set_minor_formatter(FormatStrFormatter("%1i"))

	kwargs.update(transform=ax[2,1].transAxes)
	ax[2,1].plot((-d,d),(-d,+d), **kwargs)
	ax[2,1].plot((1-d,1+d),(-d,+d), **kwargs)

	kwargs.update(transform=ax[3,1].transAxes)
	ax[3,1].plot((1-d,1+d),(1-d,1+d), **kwargs)
	ax[3,1].plot((-d,d),(1-d,1+d), **kwargs)

	ax[3,1].set_xlabel(r"$z + 1$")
	fig.text(0.5, 0.3, r"$N_{H_{X, Intrinsic}}$ $(cm^{-2})$ $T=10^{7} K$ Model", ha='center', va='center', rotation='vertical')


	plt.savefig(rootdir + "NHvsz_IGMSubtraction_1Z_0imind0inind_" + output + ".pdf", bbox_inches="tight")
	plt.close(fig)

	subFITs.close()

#---------------------------------------------------------------------------------------------------------------------#
#Unloads FITs files
if args.Graphs[0] == "Y":
	mainFITs.close()
elif args.Graphs[1] == "Y" and args.Graphs[2] == "N":
	mainFITs.close()
elif args.Graphs[1] == "N" and args.Graphs[2] == "Y":
	mainFITs.close()
elif args.Graphs[1] == "Y" and args.Graphs[2] == "Y":
	mainFITs.close()
	ZFITs.close()
elif args.lsmc:
	lsmcFITs.close()

#---------------------------------------------------------------------------------------------------------------------#
print "Program executed correctly, ending process."
