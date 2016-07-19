#!/usr/bin/env python

#Python Script for X-Ray spectral fitting of GRBs using XSPEC
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
from astropy.io import fits

#---------------------------------------------------------------------------------------------------------------------#
#User Input
parser = argparse.ArgumentParser(description="Fits models to a GRB X-Ray Spectrum over a specified energy range \
		using XSPEC: ecut fits TBabs*powerlaw (for purpose of > 2keV spectral cutting, use in conjunction with \
		non-default energy limits), standard fits zTBabs*TBabs*powerlaw, Z fits TBvarabs*TBabs*powerlaw, mcmc \
		produces a distribution of each parameter (can't be used with cFlux) for standard or Z, cFlux fits \
		cflux* to standard or Z models (must have cflux option in conjunction with standard or Z). \
		Lightcurve option will create lightcurves of each GRB.")
parser.add_argument("rootdir", help="Directory path of the program, ensure directory string ends with /")
parser.add_argument("output", help="Specifying output name string for all produced objects")
parser.add_argument("-E", "--spectrumE", type=float, nargs=2, default=[0.3, 10.0], help="Spectrum energy limits in keV, \
		left blank defaults to 0.3 and 10.0")
parser.add_argument("-ci", "--confidence", type=float, nargs="?", default=2.706, help="Confidence range for fit error \
		calculation, left blank defaults to 2.706")
parser.add_argument("-par", "--parameters", type=float, nargs=2, default=[2.0, 0.0001], help="Initial Parameter values \
		for Gamma and norm respectively, left blank defaults to 2.0 and 0.0001")
parser.add_argument("-pri", "--prior", type=float, nargs=2, default=[1.98, 0.28], help="Prior mean and varience for Gamma, \
		left blank defaults to 0.98 and 0.28 based on Evans et al. 09 population paper")
parser.add_argument("-ch", "--chain", type=int, nargs=2, default=[1000, 10000], help="Burn and Length for MCMC chains, \
		left blank defaults to 1000 and 10000")
parser.add_argument("Function", nargs=6, choices=["Y", "N"], help="Select which function you wish to perfom, \
		choice is {Y, N} and the order is ecut, standard, Z, cFlux, mcmc, lightcurve")
parser.add_argument("-g", "--add", nargs=5, help="Fit a new GRB, adding to choden existing FITs file, format is \
		TargetID GRBname GalacticColumn Redshift Type")			
args = parser.parse_args() 

rootdir = args.rootdir
output = args.output

#---------------------------------------------------------------------------------------------------------------------#
#User Input checks

#args.rootdir
if re.split("[/]", rootdir)[-1] != "":
	print "Root Directory path did not end with /"
	sys.exit()

#args.Function
if args.Function[0] == "Y" and args.Function[1] == "Y":
	print "ecut and standard cannot be invoked together"
	sys.exit()
if args.Function[0] == "Y" and args.Function[2] == "Y":
	print "ecut and Z cannot be invoked together"
	sys.exit()
if args.Function[0] == "Y" and args.Function[3] == "Y":
	print "ecut and cFlux cannot be invoked together"
	sys.exit()
if args.Function[0] == "Y" and args.Function[4] == "Y":
	print "ecut and mcmc cannot be invoked together"
	sys.exit()
if args.Function[1] == "Y" and args.Function[2] == "Y":
	print "standard and Z cannot be invoked together"
	sys.exit()
if args.Function[1] == "Y" and args.Function[3] == "N" and args.Function[4] == "N":
	print "standard must be invoked with clfux or mcmc"
	sys.exit()
if args.Function[2] == "Y" and args.Function[3] == "N" and args.Function[4] == "N":
	print "Z must be invoked with clfux or mcmc"
	sys.exit()
if args.Function[3] == "Y" and args.Function[4] == "Y":
	print "cFlux and mcmc cannot be invoked together"
	sys.exit()
if args.Function[5] == "Y" and args.Function[0] == "Y":
	print "lightcurve and ecut cannot be invoked together"
	sys.exit()
if args.Function[5] == "Y" and args.Function[1] == "Y":
	print "lightcurve and standard cannot be invoked together"
	sys.exit()
if args.Function[5] == "Y" and args.Function[2] == "Y":
	print "lightcurve and Z cannot be invoked together"
	sys.exit()
if args.Function[5] == "Y" and args.Function[3] == "Y":
	print "lightcurve and cFlux cannot be invoked together"
	sys.exit()
if args.Function[5] == "Y" and args.Function[4] == "Y":
	print "lightcurve and mcmc cannot be invoked together"
	sys.exit()

#---------------------------------------------------------------------------------------------------------------------#
#Creates a global use dictionary
GRBDict = {}

GRBDict.setdefault("TargetID", [])
GRBDict.setdefault("GRBname", [])
GRBDict.setdefault("NHGal", [])
GRBDict.setdefault("z", [])

GRBDict.setdefault("Custom", [])
GRBDict.setdefault("Spectra", [])
GRBDict.setdefault("meanpc", [])

GRBDict.setdefault("start", [])
GRBDict.setdefault("stop", [])
GRBDict.setdefault("cstart", [])
GRBDict.setdefault("cstop", [])

GRBDict.setdefault("NHint", [])
GRBDict.setdefault("NHintErr_lo", [])
GRBDict.setdefault("NHintErr_hi", [])
GRBDict.setdefault("NHintErr_st", [])

GRBDict.setdefault("Gamma", [])
GRBDict.setdefault("GammaErr_lo", [])
GRBDict.setdefault("GammaErr_hi", [])
GRBDict.setdefault("GammaErr_st", [])

GRBDict.setdefault("norm", [])
GRBDict.setdefault("normErr_lo", [])
GRBDict.setdefault("normErr_hi", [])
GRBDict.setdefault("normErr_st", [])

GRBDict.setdefault("cstat", [])
GRBDict.setdefault("dof", [])
GRBDict.setdefault("cstatRed", [])

GRBDict.setdefault("CNRate", [])
GRBDict.setdefault("TRate", [])
GRBDict.setdefault("SCounts", [])
GRBDict.setdefault("SCBKGSub", [])
GRBDict.setdefault("Exposure", [])

GRBDict.setdefault("Type", [])
GRBDict.setdefault("Author", [])
GRBDict.setdefault("Z", [])
GRBDict.setdefault("M", [])

GRBDict.setdefault("cFx", [])
GRBDict.setdefault("cFxErr_lo", [])
GRBDict.setdefault("cFxErr_hi", [])

#---------------------------------------------------------------------------------------------------------------------#
#Populates dictionary with GRBs with z from Dr. Phil Evans GRB catalogue
print "Loading GRB list..."
print

if args.add:

	GRBDict["TargetID"].append(args.add[0])
	GRBDict["GRBname"].append(args.add[1])
	GRBDict["NHGal"].append(float(args.add[2]))
	GRBDict["z"].append(float(args.add[3]))
	GRBDict["Type"].append(args.add[4])

	GRBDict["Z"].append("0.5")
	GRBDict["M"].append("N/A")		

else:	

	if args.Function[0] == "Y":
		with open(rootdir + "StageFourList.txt", "r") as fi:
			next(fi)
			for line in fi:
				if len(line.split()) == 5:
					GRBDict["TargetID"].append(line.split()[0])
					GRBDict["GRBname"].append(line.split()[1])
					GRBDict["NHGal"].append(float(line.split()[2]))
					GRBDict["z"].append(float(line.split()[3]))
					GRBDict["Type"].append(line.split()[4])

					GRBDict["Author"].append("N/A")
					GRBDict["Z"].append("0.5")
					GRBDict["M"].append("N/A")				

	elif args.Function[1] == "Y":
		with open(rootdir + "test.list", "r") as fi:#"StageFiveList.txt", "r") as fi:
			next(fi)
			for line in fi:
				if len(line.split()) == 5:
					GRBDict["TargetID"].append(line.split()[0])
					GRBDict["GRBname"].append(line.split()[1])
					GRBDict["NHGal"].append(float(line.split()[2]))
					GRBDict["z"].append(float(line.split()[3]))
					GRBDict["Type"].append(line.split()[4])

					GRBDict["Z"].append("0.5")
					GRBDict["M"].append("N/A")			

	elif args.Function[2] == "Y":
		with open(rootdir + "CucchZupList.txt", "r") as fi:
			next(fi)
			for line in fi:
				if len(line.split()) == 6:
					GRBDict["TargetID"].append(line.split()[0])
					GRBDict["GRBname"].append(line.split()[1])
					GRBDict["NHGal"].append(float(line.split()[2]))
					GRBDict["z"].append(float(line.split()[3]))
					GRBDict["Z"].append(float(line.split()[4]))
					GRBDict["M"].append(line.split()[5])
					GRBDict["Author"].append("Cup")

					GRBDict["Type"].append("N/A")

		with open(rootdir + "CucchZcvList.txt", "r") as fi:
			next(fi)
			for line in fi:
				if len(line.split()) == 6:
					GRBDict["TargetID"].append(line.split()[0])
					GRBDict["GRBname"].append(line.split()[1])
					GRBDict["NHGal"].append(float(line.split()[2]))
					GRBDict["z"].append(float(line.split()[3]))
					GRBDict["Z"].append(float(line.split()[4]))
					GRBDict["M"].append(line.split()[5])
					GRBDict["Author"].append("Ccv")

					GRBDict["Type"].append("N/A")

	elif args.Function[5] == "Y":
		with open(rootdir + "UpdatedzList.txt", "r") as fi:
			next(fi)
			for line in fi:
				if len(line.split()) == 5:
					GRBDict["TargetID"].append(line.split()[0])
					GRBDict["GRBname"].append(line.split()[1])

#---------------------------------------------------------------------------------------------------------------------#
#Downloads XRT spectra from Swift data archive which fulfill criteria
#Ignores GRBs which are already within the directory
#Populates dictionary with spectra file locations
#Uses custom spectra if found within existing GRB folder
if not os.path.isdir(rootdir + "XRTdir"):
	os.mkdir(rootdir + "XRTdir")

	print "Creating directory for spectra..."
	print

os.chdir(rootdir + "XRTdir")

for i in range(len(GRBDict["TargetID"])):
	if os.path.isdir(GRBDict.get("GRBname")[i]):
		print "Loading " + GRBDict.get("GRBname")[i] + " spectra..."

		os.chdir(GRBDict.get("GRBname")[i])

		if args.Function[5] == "Y":
			if os.path.isfile("flux.qdp") == False:
				os.system("wget http://www.swift.ac.uk/xrt_curves/" + GRBDict.get("TargetID")[i] + "/flux.qdp")
			if os.path.isfile("curve2.qdp") == False:
				os.system("wget http://www.swift.ac.uk/xrt_curves/" + GRBDict.get("TargetID")[i] + "/curve2.qdp")
			if os.path.isfile("hardrat.qdp") == False:
				os.system("wget http://www.swift.ac.uk/xrt_curves/" + GRBDict.get("TargetID")[i] + "/hardrat.qdp")

		if os.path.isfile(GRBDict.get("GRBname")[i] + "custompc.pi") == True:
			GRBDict["Spectra"].append(rootdir + "XRTdir/" + GRBDict.get("GRBname")[i] + "/" + glob.glob("*custompc.pi")[0])

			with open(rootdir + "CustomSpectraTimes.txt", "r") as fcs:
				for line in fcs:
					if line.split()[0] == GRBDict.get("GRBname")[i]:
						GRBDict["meanpc"].append(float(re.split('[\W]',line)[8]))
						GRBDict["cstart"].append(float(re.split('[\W]',line)[4]))
						GRBDict["cstop"].append(float(re.split('[\W]',line)[5]))

			with open("donereg.txt", "r") as fi:
				if glob.glob("*pc.pi")[0] == "late_timepc.pi" or glob.glob("*pc.pi")[1] == "late_timepc.pi":
					next(fi)
					for line in fi:
						GRBDict["start"].append(float(line.split()[1]))
						GRBDict["stop"].append(float(line.split()[2]))
				elif glob.glob("*pc.pi")[0] == "interval0pc.pi" or glob.glob("*pc.pi")[1] == "interval0pc.pi":
					for line in fi:
						GRBDict["start"].append(float(line.split()[1]))
						GRBDict["stop"].append(float(line.split()[2]))

			GRBDict["Custom"].append("Y")

			os.chdir(rootdir + "XRTdir/")

		elif os.path.isfile(rootdir + "CustomSpectraFiles/" + GRBDict.get("GRBname")[i] + "custom.tar.gz") == True:
			os.system("cp " + rootdir + "CustomSpectraFiles/" + GRBDict.get("GRBname")[i] + "custom.tar.gz" + " " + \
					rootdir + "XRTdir/" + GRBDict.get("GRBname")[i] + "/" )

			tfile = tf.open(GRBDict.get("GRBname")[i] + "custom.tar.gz", "r:gz")
			tfile.extractall()
			os.remove(GRBDict.get("GRBname")[i] + "custom.tar.gz")

			GRBDict["Spectra"].append(rootdir + "XRTdir/" + GRBDict.get("GRBname")[i] + "/" + glob.glob("*custompc.pi")[0])

			with open(rootdir + "CustomSpectraTimes.txt", "r") as fcs:
				for line in fcs:
					if line.split()[0] == GRBDict.get("GRBname")[i]:
						GRBDict["meanpc"].append(float(re.split('[\W]',line)[8]))
						GRBDict["cstart"].append(float(re.split('[\W]',line)[4]))
						GRBDict["cstop"].append(float(re.split('[\W]',line)[5]))

			with open("donereg.txt", "r") as fi:
				if glob.glob("*pc.pi")[0] == "late_timepc.pi" or glob.glob("*pc.pi")[1] == "late_timepc.pi":
					next(fi)
					for line in fi:
						GRBDict["start"].append(float(line.split()[1]))
						GRBDict["stop"].append(float(line.split()[2]))
				elif glob.glob("*pc.pi")[0] == "interval0pc.pi" or glob.glob("*pc.pi")[1] == "interval0pc.pi":
					for line in fi:
						GRBDict["start"].append(float(line.split()[1]))
						GRBDict["stop"].append(float(line.split()[2]))

			GRBDict["Custom"].append("Y")

			os.chdir(rootdir + "XRTdir/")

		else:
			with open("donereg.txt", "r") as fi:
				if glob.glob("*pc.pi")[0] == "late_timepc.pi":
					next(fi)
					for line in fi:
						GRBDict["meanpc"].append(float(line.split()[4]))
						GRBDict["start"].append(float(line.split()[1]))
						GRBDict["stop"].append(float(line.split()[2]))
						GRBDict["cstart"].append("N/A")
						GRBDict["cstop"].append("N/A")
				elif glob.glob("*pc.pi")[0] == "interval0pc.pi":
					for line in fi:
						GRBDict["meanpc"].append(float(line.split()[4]))
						GRBDict["start"].append(float(line.split()[1]))
						GRBDict["stop"].append(float(line.split()[2]))
						GRBDict["cstart"].append("N/A")
						GRBDict["cstop"].append("N/A")
	
			GRBDict["Spectra"].append(rootdir + "XRTdir/" + GRBDict.get("GRBname")[i] + "/" + glob.glob("*pc.pi")[0])

			GRBDict["Custom"].append("N")

			os.chdir(rootdir + "XRTdir/")

	elif os.path.isdir(GRBDict.get("GRBname")[i]) == False:
		print
		print "Downloading XRT spectra..."
		print
		print
		print
		print 

		os.system("wget http://www.swift.ac.uk/xrt_spectra/" + GRBDict.get("TargetID")[i] + "/donereg.txt")

		with open("donereg.txt") as ft:
			rname  = []
			start  = []
			stop   = []
			meanpc = []

			for line in ft:
				rname  += [line.split()[0]]
				start  += [float(line.split()[1])]
				stop   += [float(line.split()[2])]
				meanpc += [float(line.split()[4])]

			B = 0
			for j in range(len(rname)):
				if (rname[j] == "late_time") & (meanpc[j] > 0.0):
					B = 1
					Name = "late_time"
					GRBDict["meanpc"].append(meanpc[j])
					GRBDict["start"].append(start[j])
					GRBDict["stop"].append(stop[j])
					GRBDict["cstart"].append("N/A")
					GRBDict["cstop"].append("N/A")
				elif (rname[j] == "interval0") & (start[j] > 3000.0):
					B = 1
					Name = "interval0"
					GRBDict["meanpc"].append(meanpc[j])
					GRBDict["start"].append(start[j])
					GRBDict["stop"].append(stop[j])
					GRBDict["cstart"].append("N/A")
					GRBDict["cstop"].append("N/A")

		if B == 1:
			os.mkdir(rootdir + "XRTdir/" + GRBDict.get("GRBname")[i])
			os.chdir(rootdir + "XRTdir/" + GRBDict.get("GRBname")[i])

			shutil.copyfile(rootdir + "XRTdir/donereg.txt", rootdir + "XRTdir/" + GRBDict.get("GRBname")[i] + "/donereg.txt")

			os.system("wget http://www.swift.ac.uk/xrt_spectra/" + GRBDict.get("TargetID")[i] + "/" + Name + ".tar.gz")

			if args.Function[5] == "Y":
				os.system("wget http://www.swift.ac.uk/xrt_curves/" + GRBDict.get("TargetID")[i] + "/flux.qdp")
				os.system("wget http://www.swift.ac.uk/xrt_curves/" + GRBDict.get("TargetID")[i] + "/curve2.qdp")
				os.system("wget http://www.swift.ac.uk/xrt_curves/" + GRBDict.get("TargetID")[i] + "/hardrat.qdp")

			tfile = tf.open(Name + ".tar.gz", "r:gz")
			tfile.extractall()

			for wt in glob.glob("*wt*"):
				os.remove(wt)

			GRBDict["Spectra"].append(rootdir + "XRTdir/" + GRBDict.get("GRBname")[i] + "/" + glob.glob("*pc.pi")[0])

			GRBDict["Custom"].append("N")

			os.chdir(rootdir + "XRTdir/")

		if B == 0:
			print
			print
			print
			print
			print "GRB ", GRBDict.get("GRBname")[i], " failed download criteria"
			print
			print
			print
			print 

			with open(rootdir + "XRAb_ErrorLog.txt", "a") as fel:
				fel.write(GRBDict.get("GRBname")[i] + " failed criteria and was not downloaded\n")

			GRBDict["Custom"].append("X")
			GRBDict["Spectra"].append("X")
			GRBDict["meanpc"].append("123456789")
			GRBDict["start"].append("X")
			GRBDict["stop"].append("X")
			GRBDict["cstart"].append("X")
			GRBDict["cstop"].append("X")

		os.remove("donereg.txt")		

#---------------------------------------------------------------------------------------------------------------------#
#Populates dictionary with dataset flag for Starling et al. (2013) and Campana et al. (2015) GRBs
#Converts Z value into solar if Z option is invoked
if args.Function[1] == "Y":
	for i in range(len(GRBDict["GRBname"])):
		GRBDict["Author"].append("M") 

	for i in range(len(GRBDict["GRBname"])): 
		with open(rootdir + "StarlingList.txt", "r") as fs:
			for line in fs:
				if line.split()[0] == GRBDict.get("GRBname")[i]:
					GRBDict["Author"][i] = "S"

		with open(rootdir + "CampanaList.txt", "r") as fc:
			for line in fc:
				if line.split()[0] == GRBDict.get("GRBname")[i]:
					if GRBDict["Author"][i] == "S":
						GRBDict["Author"][i] = "SC"
					else:
						GRBDict["Author"][i] = "C"

elif args.Function[2] == "Y":
	for i in range(len(GRBDict["GRBname"])):
		GRBDict["Z"][i] = 10**GRBDict["Z"][i]

#---------------------------------------------------------------------------------------------------------------------#
#Creates a lightcurve for each GRB, including: hardness ratio, automatic and custom start/stop times
#background cts and flux/background limits
if args.Function[5] == "Y":
	Time_pc    = []
	TErrpos_pc = []
	TErrneg_pc = []
	Flux_pc    = []
	FErrpos_pc = []
	FErrneg_pc = []
	Rate_pc    = []

	Time_up    = []
	TErrpos_up = []
	TErrneg_up = []
	Flux_up    = []
	FErrpos_up = []
	FErrneg_up = []
	Rate_up    = []

	Time_wt    = []
	TErrpos_wt = []
	TErrneg_wt = []
	Flux_wt    = []
	FErrpos_wt = []
	FErrneg_wt = []
	Rate_wt    = []

	Time_hr    = []
	TErrpos_hr = []
	TErrneg_hr = []
	Rate_hr    = []
	RErr_hr    = []

	BGRate_pc  = []
	BGRErr_pc  = []
	FracExp_pc = []

	datatemp  = []
	datatemp2 = []
	datatemp3 = []

	for i in range(len(GRBDict["TargetID"])):
		if GRBDict.get("start")[i] != "F":
			if os.path.isdir(rootdir + "XRTdir/" + GRBDict.get("GRBname")[i]): 
				os.chdir(rootdir + "XRTdir/" + GRBDict.get("GRBname")[i])
	
				if os.path.isfile("flux.qdp") == True and os.path.isfile("curve2.qdp") == True:

					print
					print "Creating " + GRBDict.get("GRBname")[i] + " lightcurve..."		
					print

					#Load Flux data
					with open("flux.qdp", "r") as ff:
						index    = 0
						wt_start = 0
						wt_end   = 0
	 					pc_start = 0
						pc_end   = 0
						up_start = 0
						up_end   = 0
						EOF      = 0
						for line in ff: 
							temp = line.strip()
							datatemp.append(temp)

							index += 1	

							if line.strip() == "! WT data":
								wt_start = index + 1
							if line.strip() == "! PC data":
								pc_start = index + 1
								wt_end   = index - 3
							if line.strip() == "! PC Upper limit":
								pc_end   = index - 3
								up_start = index 


						if ff.readline() == "":
							EOF = index
	
					if wt_start != 0: 
						for j in range(wt_end - wt_start + 1):
							Time_wt    += [float(datatemp[wt_start + j].split()[0])]
							TErrpos_wt += [float(datatemp[wt_start + j].split()[1])]
							TErrneg_wt += [0 - float(datatemp[wt_start + j].split()[2])]
							Flux_wt    += [float(datatemp[wt_start + j].split()[3])]
							FErrpos_wt += [float(datatemp[wt_start + j].split()[4])]
							FErrneg_wt += [0 - float(datatemp[wt_start + j].split()[5])]

					if pc_end == 0:
						pc_end = EOF - 1
						for k in range(pc_end - pc_start + 1):
							Time_pc    += [float(datatemp[pc_start + k].split()[0])]
							TErrpos_pc += [float(datatemp[pc_start + k].split()[1])]
							TErrneg_pc += [0 - float(datatemp[pc_start + k].split()[2])]
							Flux_pc    += [float(datatemp[pc_start + k].split()[3])]
							FErrpos_pc += [float(datatemp[pc_start + k].split()[4])]
							FErrneg_pc += [0 - float(datatemp[pc_start + k].split()[5])]
		
					else:
						for k in range(pc_end - pc_start + 1):
							Time_pc    += [float(datatemp[pc_start + k].split()[0])]
							TErrpos_pc += [float(datatemp[pc_start + k].split()[1])]
							TErrneg_pc += [0 - float(datatemp[pc_start + k].split()[2])]
							Flux_pc    += [float(datatemp[pc_start + k].split()[3])]
							FErrpos_pc += [float(datatemp[pc_start + k].split()[4])]
							FErrneg_pc += [0 - float(datatemp[pc_start + k].split()[5])]
	
					if up_start != 0: 
						up_end = EOF - 1
						for l in range(up_end - up_start + 1):
							Time_up    += [datatemp[up_start + l].split()[0]]
							TErrpos_up += [datatemp[up_start + l].split()[1]]
							TErrneg_up += [datatemp[up_start + l].split()[2]]
							Flux_up    += [datatemp[up_start + l].split()[3]]
							FErrpos_up += [datatemp[up_start + l].split()[4]]
							FErrneg_up += [datatemp[up_start + l].split()[5]]

					#Load ct Rate data for conversion
					with open("curve2.qdp", "r") as fct:
						index    = 0
						wt_start = 0
						wt_end   = 0
	 					pc_start = 0
						pc_end   = 0
						up_start = 0
						up_end   = 0
						EOF      = 0
						for line in fct: 
							temp = line.strip()
							datatemp2.append(temp)
							index += 1	
	
							if line.strip() == "! WT data":
								wt_start = index 
							if line.strip() == "! PC data":
								pc_start = index 
								wt_end   = index - 3
							if line.strip() == "! PC Upper limit":
								pc_end   = index - 3
								up_start = index 

						if fct.readline() == "":
							EOF = index - 1
		
					if wt_start != 0: 
						Rate_wt += [float(datatemp2[wt_start].split()[3])]
						Rate_wt += [float(datatemp2[wt_end].split()[3])]
	
					if pc_end == 0:
						pc_end   = EOF
						Rate_pc += [float(datatemp2[pc_start].split()[3])]
						Rate_pc += [float(datatemp2[pc_end].split()[3])]
						for k in range(pc_end - pc_start + 1):
							BGRate_pc  += [float(datatemp2[pc_start + k].split()[6])]
							BGRErr_pc  += [float(datatemp2[pc_start + k].split()[7])]
							FracExp_pc += [float(datatemp2[pc_start + k].split()[8])]
		
					else:
						Rate_pc += [float(datatemp2[pc_start].split()[3])]
						Rate_pc += [float(datatemp2[pc_end].split()[3])]
						for k in range(pc_end - pc_start + 1):
							BGRate_pc  += [float(datatemp2[pc_start + k].split()[6])]
							BGRErr_pc  += [float(datatemp2[pc_start + k].split()[7])]
							FracExp_pc += [float(datatemp2[pc_start + k].split()[8])]
	
					if up_start != 0: 
						up_end   = EOF 
						Rate_up += [float(datatemp2[up_start].split()[3])]
						Rate_up += [float(datatemp2[up_end].split()[3])]

					#Load HR data
					def file_len(fname):
						with open(fname) as f:
							for i, l in enumerate(f):
								pass
						return i + 1

					if os.path.isfile("hardrat.qdp") == True:
						if file_len("hardrat.qdp") != 9:
							with open("hardrat.qdp", "r") as fhr:
								index    = 0
	 							pc_start = 0
								pc_end   = 0
								EOF      = 0
								for line in fhr: 
									temp = line.strip()
									datatemp3.append(temp)

									index += 1	

									if line.strip() == "! PC -- hardness ratio":
										pc_start = index + 1

								if fhr.readline() == "":
									pc_end = index - 1 

							if pc_start != 0:
								for k in range(pc_end - pc_start + 1):
									Time_hr    += [float(datatemp3[pc_start + k].split()[0])]
									TErrpos_hr += [float(datatemp3[pc_start + k].split()[1])]
									TErrneg_hr += [0 - float(datatemp3[pc_start + k].split()[2])]
									Rate_hr    += [float(datatemp3[pc_start + k].split()[3])]
									RErr_hr    += [float(datatemp3[pc_start + k].split()[4])]
	
					#Plotting the Lightcurve
					if pc_start != 0 and file_len("hardrat.qdp") != 9:
						fig = plt.figure(figsize=(8, 5))
			
						ax1 = plt.subplot2grid((3,3), (0,0), colspan=3, rowspan=2)
						ax3 = plt.subplot2grid((3,3), (2,0), colspan=3, sharex=ax1)
		
						ax1.set_xscale("log", nonposx="clip")
						ax3.set_xscale("log", nonposx="clip")
						ax1.set_yscale("log", nonposy="clip")
		
						for j in range(len(Time_wt)):
							ax1.errorbar(Time_wt[j], Flux_wt[j], xerr=[[TErrneg_wt[j]], [TErrpos_wt[j]]], \
								yerr=[[FErrneg_wt[j]], [FErrpos_wt[j]]], capsize=0, c=u"b", marker="None", ls="None")

						for j in range(len(Time_pc)):
							ax1.errorbar(Time_pc[j], Flux_pc[j], xerr=[[TErrneg_pc[j]], [TErrpos_pc[j]]], \
								yerr=[[FErrneg_pc[j]], [FErrpos_pc[j]]], capsize=0, c=u"r", marker="None", ls="None")
					
						ax1.errorbar(Time_up, Flux_up, c=u"r", marker=u"v", ls="None")

						for l in range(len(Time_hr)):
							ax3.errorbar(Time_hr[l], Rate_hr[l], xerr=[[TErrneg_hr[l]], [TErrpos_hr[l]]], \
								yerr=[RErr_hr[l]], capsize=0, c=u"r", marker="None", ls="None")

						xmin, xmax = ax1.get_xlim()
						ymin, ymax = ax1.get_ylim()

						ConvRate = Flux_pc[0] / Rate_pc[0] 

						y2max = ymax / ConvRate
						y2min = ymin / ConvRate

						ax2 = ax1.twinx()
						ax2.set_yscale("log", nonposy="clip")
						ax2.set_ylim(y2min, y2max)
						ax2.set_ylabel(r"Count Rate $(0.3-10 keV) (s^{-1})$")
						ax2.set_xlim(xmin, xmax)

						for l in range(len(BGRate_pc)):
							ax2.errorbar(Time_pc[l], BGRate_pc[l], xerr=[[TErrneg_pc[l]], [TErrpos_pc[l]]], \
								yerr=[BGRErr_pc[l]], capsize=0, c=u"c", marker="None", ls="None")

						ax1.axhline(1e-13, c=u"k")
						ax2.axhline(1e-3, c=u"c", ls=":")
						ax1.axvline(GRBDict.get("start")[i], c=u"g", ls="--")
						ax1.axvline(GRBDict.get("stop")[i], c=u"g", ls="--")
						ax3.axvline(GRBDict.get("start")[i], c=u"g", ls="--")
						ax3.axvline(GRBDict.get("stop")[i], c=u"g", ls="--")
				
						if GRBDict.get("cstart")[i] != "N/A":
							ax1.axvline(GRBDict.get("cstart")[i], c=u"m", ls=":")
							ax1.axvline(GRBDict.get("cstop")[i], c=u"m", ls=":")
							ax3.axvline(GRBDict.get("cstart")[i], c=u"m", ls=":")
							ax3.axvline(GRBDict.get("cstop")[i], c=u"m", ls=":")
	
						plt.setp(ax1.get_xticklabels(), visible=False)
						fig.subplots_adjust(hspace=0.1)

						ax1.set_title(r"GRB " + GRBDict.get("GRBname")[i] + r" Lightcurve and Hardness Ratio")
						ax3.set_xlabel(r"Time since BAT trigger $(s)$")
						ax1.set_ylabel(r"Flux $(0.3-10 keV) (erg cm^{-2}s^{-1})$")
						ax3.set_ylabel(r"Hardness Ratio")

						plt.savefig(rootdir + "XRTdir/" + GRBDict.get("GRBname")[i] + "/" + \
						GRBDict.get("GRBname")[i] + "_Flux_Lightcurve.eps", bbox_inches="tight")
						plt.close(fig)	

					else:
						fig = plt.figure(1)
			
						fig, ax = plt.subplots(1, figsize=(8, 5))
		
						ax.set_xscale("log", nonposx="clip")
						ax.set_yscale("log", nonposy="clip")
						ax.set_ylim(1e-15,1e-8)
		
						for j in range(len(Time_wt)):
							ax.errorbar(Time_wt[j], Flux_wt[j], xerr=[[TErrneg_wt[j]], [TErrpos_wt[j]]], \
								yerr=[[FErrneg_wt[j]], [FErrpos_wt[j]]], capsize=0, c=u"b", marker="None", ls="None")

						for j in range(len(Time_pc)):
							ax.errorbar(Time_pc[j], Flux_pc[j], xerr=[[TErrneg_pc[j]], [TErrpos_pc[j]]], \
								yerr=[[FErrneg_pc[j]], [FErrpos_pc[j]]], capsize=0, c=u"r", marker="None", ls="None")
					
						ax.errorbar(Time_up, Flux_up, c=u"r", marker=u"v", ls="None")

						xmin, xmax = ax.get_xlim()
						ymin, ymax = ax.get_ylim()

						ConvRate = Flux_pc[0] / Rate_pc[0] 

						y2max = ymax / ConvRate
						y2min = ymin / ConvRate

						ax2 = ax.twinx()
						ax2.set_yscale("log", nonposy="clip")
						ax2.set_ylim(y2min, y2max)
						ax2.set_ylabel(r"Count Rate $(0.3-10 keV) (s^{-1})$")
						ax2.set_xlim(xmin, xmax)
	
						ax.axhline(1e-13, c=u"k")
						ax.axhline(1e-3, c=u"c")
						ax.axvline(GRBDict.get("start")[i], c=u"g", ls="--")
						ax.axvline(GRBDict.get("stop")[i], c=u"g", ls="--")
				
						if GRBDict.get("cstart")[i] != "N/A":
							ax.axvline(GRBDict.get("cstart")[i], c=u"m", ls=":")
							ax.axvline(GRBDict.get("cstop")[i], c=u"m", ls=":")

						ax.set_title(r"GRB " + GRBDict.get("GRBname")[i] + r" Lightcurve, No Hardness Ratio")
						ax.set_xlabel(r"Time since BAT trigger $(s)$")
						ax.set_ylabel(r"Flux $(0.3-10 keV) (erg cm^{-2}s^{-1})$")

						plt.savefig(rootdir + "Lightcurves/" + GRBDict.get("GRBname")[i] + "_Flux_Lightcurve.eps", bbox_inches="tight")
						plt.close(fig)						
			
				#Reset lists
				Time_pc    = []
				TErrpos_pc = []
				TErrneg_pc = []
				Flux_pc    = []
				FErrpos_pc = []
				FErrneg_pc = []
				Rate_pc    = []
	
				Time_up    = []
				TErrpos_up = []
				TErrneg_up = []
				Flux_up    = []
				FErrpos_up = []
				FErrneg_up = []
				Rate_up    = []

				Time_wt    = []
				TErrpos_wt = []
				TErrneg_wt = []
				Flux_wt    = []
				FErrpos_wt = []
				FErrneg_wt = []
				Rate_wt    = []

				Time_hr    = []
				TErrpos_hr = []
				TErrneg_hr = []
				Rate_hr    = []
				RErr_hr    = []

				BGRate_pc  = []
				BGRErr_pc  = []
				FracExp_pc = []
	
				datatemp  = []	
				datatemp2 = []
				datatemp3 = []

		else:
			print
			print GRBDict.get("GRBname")[i] + " failed and a lightcurve was not created"
			print
		
			with open(rootdir + "XRAb_ErrorLog.txt", "a") as fel:
				fel.write(GRBDict.get("GRBname")[i] + " failed and a lightcurve was not created\n")

	print "Program executed correctly, ending process."
	sys.exit()

#---------------------------------------------------------------------------------------------------------------------#
#PyXSpec code for fitting each XRT spectrum
#produces a spectral plot within each GRB subdirectory 
#populates dictionary with fit values
print "Entering XSPEC interface..."
print
print
print
print

for i in range(len(GRBDict["GRBname"])):
	if not GRBDict.get("Spectra")[i] == "F":
		infile = GRBDict.get("Spectra")[i]
		os.chdir(infile[:infile.rfind("/")])
	
		s = Spectrum(infile)

		AllData.ignore("bad")
		s.ignore("**-" + str(args.spectrumE[0]) + ", " + str(args.spectrumE[1]) + "-**")

		Xset.abund = "wilm"
		Xset.xsect = "vern"

		Fit.method = "leven 100"
		Fit.query = "yes"
		Fit.statMethod = "cstat"

		if args.Function[0] == "Y":
			m = Model("TBabs(powerlaw)", \
				setPars={1:GRBDict.get("NHGal")[i], 2:args.parameters[0], 3:args.parameters[1]})

		elif args.Function[1] == "Y":
			m = Model("TBabs(zTBabs(powerlaw))", \
				setPars={1:GRBDict.get("NHGal")[i], 3:GRBDict.get("z")[i], \
					4:args.parameters[0], 5:args.parameters[1]})

			nHI = m.zTBabs.nH

		elif args.Function[2] == "Y":
			m = Model("TBabs(TBvarabs(powerlaw))", setPars={1:GRBDict.get("NHGal")[i], \
					4:GRBDict.get("Z")[i], 5:GRBDict.get("Z")[i], \
					6:GRBDict.get("Z")[i], 7:GRBDict.get("Z")[i], 8:GRBDict.get("Z")[i], \
					9:GRBDict.get("Z")[i], 10:GRBDict.get("Z")[i], 11:GRBDict.get("Z")[i], \
					12:GRBDict.get("Z")[i], 13:GRBDict.get("Z")[i], 14:GRBDict.get("Z")[i], \
					15:GRBDict.get("Z")[i], 16:GRBDict.get("Z")[i], 17:GRBDict.get("Z")[i], \
					18:GRBDict.get("Z")[i], 19:GRBDict.get("Z")[i], \
					43:GRBDict.get("z")[i], 44:args.parameters[0], 45:args.parameters[1]})

			nHI = m.TBvarabs.nH

		nHG = m.TBabs.nH
		nHG.frozen = True

		G = m.powerlaw.PhoIndex
		n = m.powerlaw.norm

		Fit.bayes = "on"
		if args.Function[0] == "Y":
			G.prior = "gauss " + str(args.prior[0]) + " " + str(args.prior[1])
			n.prior = "jeff"
		if args.Function[1] == "Y" or args.Function[2] == "Y":
			nHI.prior = "jeff"
			G.prior = "gauss " + str(args.prior[0]) + " " + str(args.prior[1])
			n.prior = "jeff"

		#-----------------------------------------------------------------------------------------------------#	
		#Writes spectral information into Dictionary
		GRBDict["CNRate"].append(s.rate[0])
		GRBDict["TRate"].append(s.rate[2])
		GRBDict["SCounts"].append(int(s.rate[2] * s.exposure))
		GRBDict["SCBKGSub"].append(s.rate[0] * s.exposure)
		GRBDict["Exposure"].append(s.exposure)

		#-----------------------------------------------------------------------------------------------------#	
		#Performs the fit three times to (minorally) shake parameters
		Fit.perform()
		Fit.perform()
		Fit.perform()
	
		#-----------------------------------------------------------------------------------------------------#
		#Adds cflux into fitting if specified
		#cFx commands originally created by Dr. Andy Beardmore, modified by me
		if args.Function[3] == "Y":		
			best_fit_pars = []

			for component in m.componentNames:
				c = getattr(m, component)

				for parameter in c.parameterNames:
					p = getattr(c, parameter)

					comp_par = c.name + "." + p.name
					best_fit_pars.append((comp_par, p.values[0], p.values[1]))

			AllModels.calcFlux(str(args.spectrumE[0]) + " " + str(args.spectrumE[1]))
	   		Fx = s.flux[0]

			startpars = [args.spectrumE[0], args.spectrumE[1], np.log10(Fx)]

			for j in range(len(best_fit_pars)):
				if best_fit_pars[j][2] > 0.0:
					startpars.append(str(best_fit_pars[j][1]) + " " + str(best_fit_pars[j][1] * Fit.delta))
				else:
					startpars.append(str(best_fit_pars[j][1]) + " " + str(-best_fit_pars[j][1] * Fit.delta))					

			m = Model("cflux*" + m.expression, setPars = startpars) 

			if args.Function[1] == "Y":
				nHI = m.zTBabs.nH

			elif args.Function[2] == "Y":
				nHI = m.TBvarabs.nH

			cfl = m.cflux.lg10Flux
			G   = m.powerlaw.PhoIndex
			n   = m.powerlaw.norm
			n.frozen = True
			
			Fit.perform()
			Fit.perform()
			Fit.perform()

			#Converts confidence values into +/- values for matplotlib plotting
			GRBDict["cFx"].append(10.0**cfl.values[0])
			GRBDict["cFxErr_lo"].append(10.0**cfl.values[0] - 10.0**cfl.error[0])
			GRBDict["cFxErr_hi"].append(10.0**cfl.error[1] - 10.0**cfl.values[0])
	
		#-----------------------------------------------------------------------------------------------------#
		#If MCMC function is specified, fitting includes XSPEC mcmc production
		#Creates xcm file for MCMC stat extraction (pyxspec chain class lacks the functionality)
		#and then executes it
		if args.Function[4] == "Y":	
			os.chdir(rootdir + "XRTdir/" + GRBDict.get("GRBname")[i])
			if os.path.isfile(GRBDict.get("GRBname")[i] + "_" + output + "_chain1.fits"):
				AllChains += GRBDict.get("GRBname")[i] + "_" + output + "_chain1.fits"
				AllChains += GRBDict.get("GRBname")[i] + "_" + output + "_chain2.fits"
				AllChains += GRBDict.get("GRBname")[i] + "_" + output + "_chain3.fits"
				AllChains += GRBDict.get("GRBname")[i] + "_" + output + "_chain4.fits"
				AllChains += GRBDict.get("GRBname")[i] + "_" + output + "_chain5.fits"
				AllChains.show()
	
			else:
				AllChains.defBurn = args.chain[0]
				AllChains.defLength = args.chain[1]
				AllChains.show()
				c  = Chain(GRBDict.get("GRBname")[i] + "_" + output + "_chain1.fits")
				c2 = Chain(GRBDict.get("GRBname")[i] + "_" + output + "_chain2.fits")
				c3 = Chain(GRBDict.get("GRBname")[i] + "_" + output + "_chain3.fits")
				c4 = Chain(GRBDict.get("GRBname")[i] + "_" + output + "_chain4.fits")
				c5 = Chain(GRBDict.get("GRBname")[i] + "_" + output + "_chain5.fits")

			with open(rootdir + "mcmc.xcm", "w") as fw:
				fw.write("set xs_return_result 1\n")
				fw.write("set fileid [open " + rootdir + "MCMC_" + output + "_mcmcstat.txt a]\n")
				fw.write("chain load " + rootdir + "XRTdir/" + GRBDict.get("GRBname")[i] + \
						"/" + GRBDict.get("GRBname")[i] + "_" + output + "_chain1.fits\n")
				fw.write("chain load " + rootdir + "XRTdir/" + GRBDict.get("GRBname")[i] + \
						"/" + GRBDict.get("GRBname")[i] + "_" + output + "_chain2.fits\n")
				fw.write("chain load " + rootdir + "XRTdir/" + GRBDict.get("GRBname")[i] + \
						"/" + GRBDict.get("GRBname")[i] + "_" + output + "_chain3.fits\n")
				fw.write("chain load " + rootdir + "XRTdir/" + GRBDict.get("GRBname")[i] + \
						"/" + GRBDict.get("GRBname")[i] + "_" + output + "_chain4.fits\n")
				fw.write("chain load " + rootdir + "XRTdir/" + GRBDict.get("GRBname")[i] + \
						"/" + GRBDict.get("GRBname")[i] + "_" + output + "_chain5.fits\n")
				fw.write("chain stat 2\n")
				fw.write("tclout chain stat\n")
				fw.write("set nH [string trim $xspec_tclout]\n")
				if args.Function[1] == "Y":
					fw.write("chain stat 4\n")
					fw.write("tclout chain stat\n")
					fw.write("set G [string trim $xspec_tclout]\n")
					fw.write("chain stat 5\n")
					fw.write("tclout chain stat\n")
					fw.write("set n [string trim $xspec_tclout]\n")
				elif args.Function[2] == "Y":
					fw.write("chain stat 44\n")
					fw.write("tclout chain stat\n")
					fw.write("set G [string trim $xspec_tclout]\n")
					fw.write("chain stat 45\n")
					fw.write("tclout chain stat\n")
					fw.write("set n [string trim $xspec_tclout]\n")
				fw.write('puts $fileid "'+ GRBDict.get("GRBname")[i] + ' [lindex $nH] [lindex $G] [lindex $n]"\n')
				fw.write("close $fileid\n")
				fw.write("exit\n")

			os.system("xspec - "+ rootdir + "mcmc.xcm")

		#-----------------------------------------------------------------------------------------------------#
		#Fit confidence value estimation, if MCMC is envoked error will use MCMC to estimate
		if args.Function[0] == "Y":
		  	Fit.error(str(args.confidence) + " 2")
			Fit.error(str(args.confidence) + " 3")

		if args.Function[1] == "Y" and args.Function[3] == "Y":
		  	Fit.error(str(args.confidence) + " 7")
			Fit.error(str(args.confidence) + " 3")
			Fit.error(str(args.confidence) + " 5")

		elif args.Function[1] == "Y" and args.Function[4] == "Y":
		  	Fit.error(str(args.confidence) + " 4")
			Fit.error(str(args.confidence) + " 5")
			Fit.error(str(args.confidence) + " 2")

		if args.Function[2] == "Y" and args.Function[3] == "Y":
		  	Fit.error(str(args.confidence) + " 47")
			Fit.error(str(args.confidence) + " 3")
			Fit.error(str(args.confidence) + " 5")

		elif args.Function[2] == "Y" and args.Function[4] == "Y":
			Fit.error(str(args.confidence) + " 44")
			Fit.error(str(args.confidence) + " 45")
			Fit.error(str(args.confidence) + " 2")

		#-----------------------------------------------------------------------------------------------------#
		#Writes Fit values to Dictionary 
		GRBDict["Gamma"].append(G.values[0])
		GRBDict["norm"].append(n.values[0])
		GRBDict["cstat"].append(Fit.statistic)
		GRBDict["dof"].append(Fit.dof)
		GRBDict["cstatRed"].append((Fit.statistic / Fit.dof))
		GRBDict["NHint"].append(nHI.values[0])

		#Converts confidence values into +/- values for matplotlib plotting
		if args.Function[1] == "Y" or args.Function[2] == "Y":
			if nHI.error[0] > 0 and nHI.error[1] > 0:
				GRBDict["NHintErr_lo"].append(nHI.values[0] - nHI.error[0])
				GRBDict["NHintErr_hi"].append(nHI.error[1] - nHI.values[0])
				GRBDict["NHintErr_st"].append(nHI.error[2])

			elif nHI.error[0] <= 0 and nHI.error[1] > 0:
				GRBDict["NHintErr_lo"].append(nHI.error[0])
				GRBDict["NHintErr_hi"].append(nHI.error[1] - nHI.values[0])
				GRBDict["NHintErr_st"].append(nHI.error[2])

			elif nHI.error[0] > 0 and nHI.error[1] <= 0:
				GRBDict["NHintErr_lo"].append(nHI.values[0] - nHI.error[0])
				GRBDict["NHintErr_hi"].append(nHI.error[1])
				GRBDict["NHintErr_st"].append(nHI.error[2])

			else:
				GRBDict["NHintErr_lo"].append(nHI.error[0])
				GRBDict["NHintErr_hi"].append(nHI.error[1])
				GRBDict["NHintErr_st"].append(nHI.error[2])

		if G.error[0] > 0 and G.error[1] > 0:
			GRBDict["GammaErr_lo"].append(G.values[0] - G.error[0])
			GRBDict["GammaErr_hi"].append(G.error[1] - G.values[0])
			GRBDict["GammaErr_st"].append(G.error[2])

		elif G.error[0] <= 0 and G.error[1] > 0:
			GRBDict["GammaErr_lo"].append(G.error[0])
			GRBDict["GammaErr_hi"].append(G.error[1] - G.values[0])
			GRBDict["GammaErr_st"].append(G.error[2])

		elif G.error[0] > 0 and G.error[1] <= 0:
			GRBDict["GammaErr_lo"].append(G.values[0] - G.error[0])
			GRBDict["GammaErr_hi"].append(G.error[1])
			GRBDict["GammaErr_st"].append(G.error[2])

		else:
			GRBDict["GammaErr_lo"].append(G.error[0])
			GRBDict["GammaErr_hi"].append(G.error[1])
			GRBDict["GammaErr_st"].append(G.error[2])

		if args.Function[0] == "Y" or args.Function[4] == "Y":
			if n.error[0] > 0 and n.error[1] > 0:
				GRBDict["normErr_lo"].append(n.values[0] - n.error[0])
				GRBDict["normErr_hi"].append(n.error[1] - n.values[0])
				GRBDict["normErr_st"].append(n.error[2])

			elif n.error[0] <= 0 and n.error[1] > 0:
				GRBDict["normErr_lo"].append(n.error[0])
				GRBDict["normErr_hi"].append(n.error[1] - n.values[0])
				GRBDict["normErr_st"].append(n.error[2])

			elif n.error[0] > 0 and n.error[1] <= 0:
				GRBDict["normErr_lo"].append(n.values[0] - n.error[0])
				GRBDict["normErr_hi"].append(n.error[1])
				GRBDict["normErr_st"].append(n.error[2])

			else:
				GRBDict["normErr_lo"].append(n.error[0])
				GRBDict["normErr_hi"].append(n.error[1])
				GRBDict["normErr_st"].append(n.error[2])

		#-----------------------------------------------------------------------------------------------------#
		#Creates spectral plot for each GRB, placing them within each GRBs respective XRTdir folder
		Plot.xAxis = "keV"
		Plot.background = False
		
		if args.Function[0] == "Y":
			Plot.device = "Ecut_" + output + "/cps"
		elif args.Function[1] == "Y":
			Plot.device = "XRAb_" + output + "/cps"
		elif args.Function[2] == "Y":
			Plot.device = "Z_" + output + "/cps"
		elif args.Function[4] == "Y":
			Plot.device = "MCMC_" + output + "/cps"

		Plot.setRebin(5.0, 5)
		Plot.addCommand("label t X-Ray Spectrum of GRB" + GRBDict.get("GRBname")[i])
		Plot.addCommand("col 4 on 2")
		Plot("ldata resid ratio")
	
		#-----------------------------------------------------------------------------------------------------#
		#Clears PyXSPEC classes for next run		
		AllData.clear()
		AllModels.clear()
		AllChains.clear()

	#If GRB failed, dictionary is appended with failure values
	else:	
		GRBDict["NHint"].append("123456789")
		GRBDict["NHintErr_lo"].append("123456789")
		GRBDict["NHintErr_hi"].append("123456789")
		GRBDict["NHintErr_st"].append("X")
		GRBDict["Gamma"].append("123456789")
		GRBDict["GammaErr_lo"].append("123456789")
		GRBDict["GammaErr_hi"].append("123456789")
		GRBDict["GammaErr_st"].append("X")
		GRBDict["norm"].append("123456789")
		GRBDict["normErr_lo"].append("123456789")
		GRBDict["normErr_hi"].append("123456789")
		GRBDict["normErr_st"].append("X")
		GRBDict["cFx"].append("123456789")
		GRBDict["cFxErr_lo"].append("123456789")
		GRBDict["cFxErr_hi"].append("123456789")
		GRBDict["cstat"].append("123456789")
		GRBDict["dof"].append("123456789")
		GRBDict["cstatRed"].append("123456789")
		GRBDict["CNRate"].append("123456789")
		GRBDict["TRate"].append("123456789")
		GRBDict["SCounts"].append("123456789")
		GRBDict["SCBKGSub"].append("123456789")
		GRBDict["Exposure"].append("123456789")

#---------------------------------------------------------------------------------------------------------------------#
#Compiles all output data into FITS file, updating existing file if args.add is invoked

if args.add:

	tbhdu = fits.BinTableHDU.from_columns(
			[fits.Column(name='GRBname', format='7A', unit=None, array=GRBDict["GRBname"]),
			fits.Column(name='Author', format='3A', unit=None, array=GRBDict["Author"]),
			fits.Column(name='Type', format='2A', unit=None, array=GRBDict["Type"]),
			fits.Column(name='Custom', format='1A', unit=None, array=GRBDict["Custom"]),
			fits.Column(name='z', format='E', unit=None, array=GRBDict["z"]),
			fits.Column(name='NH_Gal', format='E', unit='10^22 cm^-2', array=GRBDict["NHGal"]),
			fits.Column(name='NH_Int', format='D', unit='10^22 cm^-2', array=GRBDict["NHint"]),
			fits.Column(name='NH_IntErr_low', format='D', unit='10^22 cm^-2', array=GRBDict["NHintErr_lo"]),
			fits.Column(name='NH_IntErr_hi', format='D', unit='10^22 cm^-2', array=GRBDict["NHintErr_hi"]),
			fits.Column(name='NH_IntErr_stat', format='9A', unit=None, array=GRBDict["NHintErr_st"]),
			fits.Column(name='Gamma', format='D', unit=None, array=GRBDict["Gamma"]),
			fits.Column(name='GammaErr_low', format='D', unit=None, array=GRBDict["GammaErr_lo"]),
			fits.Column(name='GammaErr_hi', format='D', unit=None, array=GRBDict["GammaErr_hi"]),
			fits.Column(name='GammaErr_stat', format='9A', unit=None, array=GRBDict["GammaErr_st"]),
			fits.Column(name='Norm', format='D', unit=None, array=GRBDict["norm"]),
			fits.Column(name='NormErr_low', format='D', unit=None, array=GRBDict["normErr_lo"]),
			fits.Column(name='NormErr_hi', format='D', unit=None, array=GRBDict["normErr_hi"]),
			fits.Column(name='NormErr_stat', format='9A', unit=None, array=GRBDict["normErr_st"]),
			fits.Column(name='cFx', format='D', unit='erg cm^-2 s^-1', array=GRBDict["cFx"]),
			fits.Column(name='cFxErr_low', format='D', unit='erg cm^-2 s^-1', array=GRBDict["cFxErr_lo"]),
			fits.Column(name='cFxErr_hi', format='D', unit='erg cm^-2 s^-1', array=GRBDict["cFxErr_hi"]),
			fits.Column(name='cstat', format='D', unit=None, array=GRBDict["cstat"]),
			fits.Column(name='dof', format='J', unit=None, array=GRBDict["dof"]),
			fits.Column(name='Red_cstat', format='D', unit=None, array=GRBDict["cstatRed"]),
			fits.Column(name='meanpc', format='D', unit='s', array=GRBDict["meanpc"]),
			fits.Column(name='CN_Rate', format='D', unit='cts s^-1', array=GRBDict["CNRate"]),
			fits.Column(name='T_Rate', format='D', unit='cts s^-1', array=GRBDict["TRate"]),
			fits.Column(name='SCounts', format='J', unit='cts', array=GRBDict["SCounts"]),
			fits.Column(name='SC_BKGSub', format='D', unit='cts', array=GRBDict["SCBKGSub"]),
			fits.Column(name='Exposure', format='D', unit='s', array=GRBDict["Exposure"])])

	tbhdu.writeto(rootdir + "XRAb_FITSTable_" + output + "_temp.fits", clobber=True)

	t1 = fits.open(rootdir + "XRAb_FITSTable_" + output + ".fits")
	t2 = fits.open(rootdir + "XRAb_FITSTable_" + output + "_temp.fits")

	nrows1 = t1[1].data.shape[0]
	nrows2 = t2[1].data.shape[0]
	nrows  = nrows1 + nrows2
	tbhdu  = fits.new_table(t1[1].columns, nrows=nrows)
	for name in t1[1].columns.names:
		tbhdu.data.field(name)[nrows1:] = t2[1].data.field(name)

	tbhdu.writeto(rootdir + "XRAb_FITSTable_" + output + ".fits", clobber=True)

	os.remove(rootdir + "XRAb_FITSTable_" + output + "_temp.fits")

else:

	if args.Function[0] == "Y":
		tbhdu = fits.BinTableHDU.from_columns(
				[fits.Column(name='GRBname', format='7A', unit=None, array=GRBDict["GRBname"]),
				fits.Column(name='Type', format='2A', unit=None, array=GRBDict["Type"]),
				fits.Column(name='Custom', format='1A', unit=None, array=GRBDict["Custom"]),
				fits.Column(name='z', format='E', unit=None, array=GRBDict["z"]),
				fits.Column(name='NH_Gal', format='E', unit='10^22 cm^-2', array=GRBDict["NHGal"]),
				fits.Column(name='Gamma', format='D', unit=None, array=GRBDict["Gamma"]),
				fits.Column(name='GammaErr_low', format='D', unit=None, array=GRBDict["GammaErr_lo"]),
				fits.Column(name='GammaErr_hi', format='D', unit=None, array=GRBDict["GammaErr_hi"]),
				fits.Column(name='GammaErr_stat', format='9A', unit=None, array=GRBDict["GammaErr_st"]),
				fits.Column(name='Norm', format='D', unit=None, array=GRBDict["norm"]),
				fits.Column(name='NormErr_low', format='D', unit=None, array=GRBDict["normErr_lo"]),
				fits.Column(name='NormErr_hi', format='D', unit=None, array=GRBDict["normErr_hi"]),
				fits.Column(name='NormErr_stat', format='9A', unit=None, array=GRBDict["normErr_st"]),
				fits.Column(name='cstat', format='D', unit=None, array=GRBDict["cstat"]),
				fits.Column(name='dof', format='J', unit=None, array=GRBDict["dof"]),
				fits.Column(name='Red_cstat', format='D', unit=None, array=GRBDict["cstatRed"]),
				fits.Column(name='meanpc', format='D', unit='s', array=GRBDict["meanpc"]),
				fits.Column(name='CN_Rate', format='D', unit='cts s^-1', array=GRBDict["CNRate"]),
				fits.Column(name='T_Rate', format='D', unit='cts s^-1', array=GRBDict["TRate"]),
				fits.Column(name='SCounts', format='J', unit='cts', array=GRBDict["SCounts"]),
				fits.Column(name='SC_BKGSub', format='D', unit='cts', array=GRBDict["SCBKGSub"]),
				fits.Column(name='Exposure', format='D', unit='s', array=GRBDict["Exposure"])])

		tbhdu.writeto(rootdir + "Ecut_FITSTable_" + output + ".fits", clobber=True)

	elif args.Function[1] == "Y" and args.Function[4] == "N":
		tbhdu = fits.BinTableHDU.from_columns(
				[fits.Column(name='GRBname', format='7A', unit=None, array=GRBDict["GRBname"]),
				fits.Column(name='Author', format='3A', unit=None, array=GRBDict["Author"]),
				fits.Column(name='Type', format='2A', unit=None, array=GRBDict["Type"]),
				fits.Column(name='Custom', format='1A', unit=None, array=GRBDict["Custom"]),
				fits.Column(name='z', format='E', unit=None, array=GRBDict["z"]),
				fits.Column(name='NH_Gal', format='E', unit='10^22 cm^-2', array=GRBDict["NHGal"]),
				fits.Column(name='NH_Int', format='D', unit='10^22 cm^-2', array=GRBDict["NHint"]),
				fits.Column(name='NH_IntErr_low', format='D', unit='10^22 cm^-2', array=GRBDict["NHintErr_lo"]),
				fits.Column(name='NH_IntErr_hi', format='D', unit='10^22 cm^-2', array=GRBDict["NHintErr_hi"]),
				fits.Column(name='NH_IntErr_stat', format='9A', unit=None, array=GRBDict["NHintErr_st"]),
				fits.Column(name='Gamma', format='D', unit=None, array=GRBDict["Gamma"]),
				fits.Column(name='GammaErr_low', format='D', unit=None, array=GRBDict["GammaErr_lo"]),
				fits.Column(name='GammaErr_hi', format='D', unit=None, array=GRBDict["GammaErr_hi"]),
				fits.Column(name='GammaErr_stat', format='9A', unit=None, array=GRBDict["GammaErr_st"]),
				fits.Column(name='Norm', format='D', unit=None, array=GRBDict["norm"]),
				fits.Column(name='NormErr_low', format='D', unit=None, array=GRBDict["normErr_lo"]),
				fits.Column(name='NormErr_hi', format='D', unit=None, array=GRBDict["normErr_hi"]),
				fits.Column(name='NormErr_stat', format='9A', unit=None, array=GRBDict["normErr_st"]),
				fits.Column(name='cFx', format='D', unit='erg cm^-2 s^-1', array=GRBDict["cFx"]),
				fits.Column(name='cFxErr_low', format='D', unit='erg cm^-2 s^-1', array=GRBDict["cFxErr_lo"]),
				fits.Column(name='cFxErr_hi', format='D', unit='erg cm^-2 s^-1', array=GRBDict["cFxErr_hi"]),
				fits.Column(name='cstat', format='D', unit=None, array=GRBDict["cstat"]),
				fits.Column(name='dof', format='J', unit=None, array=GRBDict["dof"]),
				fits.Column(name='Red_cstat', format='D', unit=None, array=GRBDict["cstatRed"]),
				fits.Column(name='meanpc', format='D', unit='s', array=GRBDict["meanpc"]),
				fits.Column(name='CN_Rate', format='D', unit='cts s^-1', array=GRBDict["CNRate"]),
				fits.Column(name='T_Rate', format='D', unit='cts s^-1', array=GRBDict["TRate"]),
				fits.Column(name='SCounts', format='J', unit='cts', array=GRBDict["SCounts"]),
				fits.Column(name='SC_BKGSub', format='D', unit='cts', array=GRBDict["SCBKGSub"]),
				fits.Column(name='Exposure', format='D', unit='s', array=GRBDict["Exposure"])])

		tbhdu.writeto(rootdir + "XRAb_FITSTable_" + output + ".fits", clobber=True)

	elif args.Function[2] == "Y" and args.Function[4] == "N":
		tbhdu = fits.BinTableHDU.from_columns(
				[fits.Column(name='GRBname', format='7A', unit=None, array=GRBDict["GRBname"]),
				fits.Column(name='Author', format='3A', unit=None, array=GRBDict["Author"]),
				fits.Column(name='Custom', format='1A', unit=None, array=GRBDict["Custom"]),
				fits.Column(name='z', format='E', unit=None, array=GRBDict["z"]),
				fits.Column(name='Z', format='E', unit='log10 Solar Metallicty', array=GRBDict["Z"]),
				fits.Column(name='M', format='2A', unit=None, array=GRBDict["M"]),
				fits.Column(name='NH_Gal', format='E', unit='10^22 cm^-2', array=GRBDict["NHGal"]),
				fits.Column(name='NH_Int', format='D', unit='10^22 cm^-2', array=GRBDict["NHint"]),
				fits.Column(name='NH_IntErr_low', format='D', unit='10^22 cm^-2', array=GRBDict["NHintErr_lo"]),
				fits.Column(name='NH_IntErr_hi', format='D', unit='10^22 cm^-2', array=GRBDict["NHintErr_hi"]),
				fits.Column(name='NH_IntErr_stat', format='9A', unit=None, array=GRBDict["NHintErr_st"]),
				fits.Column(name='Gamma', format='D', unit=None, array=GRBDict["Gamma"]),
				fits.Column(name='GammaErr_low', format='D', unit=None, array=GRBDict["GammaErr_lo"]),
				fits.Column(name='GammaErr_hi', format='D', unit=None, array=GRBDict["GammaErr_hi"]),
				fits.Column(name='GammaErr_stat', format='9A', unit=None, array=GRBDict["GammaErr_st"]),
				fits.Column(name='Norm', format='D', unit=None, array=GRBDict["norm"]),
				fits.Column(name='NormErr_low', format='D', unit=None, array=GRBDict["normErr_lo"]),
				fits.Column(name='NormErr_hi', format='D', unit=None, array=GRBDict["normErr_hi"]),
				fits.Column(name='NormErr_stat', format='9A', unit=None, array=GRBDict["normErr_st"]),
				fits.Column(name='cFx', format='D', unit='erg cm^-2 s^-1', array=GRBDict["cFx"]),
				fits.Column(name='cFxErr_low', format='D', unit='erg cm^-2 s^-1', array=GRBDict["cFxErr_lo"]),
				fits.Column(name='cFxErr_hi', format='D', unit='erg cm^-2 s^-1', array=GRBDict["cFxErr_hi"]),
				fits.Column(name='cstat', format='D', unit=None, array=GRBDict["cstat"]),
				fits.Column(name='dof', format='J', unit=None, array=GRBDict["dof"]),
				fits.Column(name='Red_cstat', format='D', unit=None, array=GRBDict["cstatRed"]),
				fits.Column(name='meanpc', format='D', unit='s', array=GRBDict["meanpc"]),
				fits.Column(name='CN_Rate', format='D', unit='cts s^-1', array=GRBDict["CNRate"]),
				fits.Column(name='T_Rate', format='D', unit='cts s^-1', array=GRBDict["TRate"]),
				fits.Column(name='SCounts', format='J', unit='cts', array=GRBDict["SCounts"]),
				fits.Column(name='SC_BKGSub', format='D', unit='cts', array=GRBDict["SCBKGSub"]),
				fits.Column(name='Exposure', format='D', unit='s', array=GRBDict["Exposure"])])

		tbhdu.writeto(rootdir + "Z_FITSTable_" + output + ".fits", clobber=True)

#---------------------------------------------------------------------------------------------------------------------#
#Loads chain stat file for compilation into FITS file
if args.Function[4] == "Y":
	GRB = []

	Ch1_par2_mean = []
	Ch2_par2_mean = []
	Ch3_par2_mean = []
	Ch4_par2_mean = []
	Ch5_par2_mean = []

	par2_allmean = []
	par2_varmean = []
	par2_allvar  = []
	par2_rubgel  = []

	Ch1_par4_mean = []
	Ch2_par4_mean = []
	Ch3_par4_mean = []
	Ch4_par4_mean = []
	Ch5_par4_mean = []

	par4_allmean = []
	par4_varmean = []
	par4_allvar  = []
	par4_rubgel  = []

	Ch1_par5_mean = []
	Ch2_par5_mean = []
	Ch3_par5_mean = []
	Ch4_par5_mean = []
	Ch5_par5_mean = []

	par5_allmean = []
	par5_varmean = []
	par5_allvar  = []
	par5_rubgel  = []

	Ch1_frac = []
	Ch2_frac = []
	Ch3_frac = []
	Ch4_frac = []
	Ch5_frac = []

	for i in range(len(GRBDict["GRBname"])):
		with open(rootdir + "MCMC_" + output + "_mcmcstat.txt", "r") as fr: 
			for line in fr:
				if line.split()[0] == GRBDict.get("GRBname")[i]:
					GRB += [line.split()[0]]

					Ch1_par2_mean += [line.split()[1]]
					Ch2_par2_mean += [line.split()[2]]
					Ch3_par2_mean += [line.split()[3]]
					Ch4_par2_mean += [line.split()[4]]
					Ch5_par2_mean += [line.split()[5]]

					par2_allmean += [line.split()[6]]
					par2_varmean += [line.split()[7]]
					par2_allvar  += [line.split()[8]]
					par2_rubgel  += [line.split()[9]]

					Ch1_par4_mean += [line.split()[15]]
					Ch2_par4_mean += [line.split()[16]]
					Ch3_par4_mean += [line.split()[17]]
					Ch4_par4_mean += [line.split()[18]]
					Ch5_par4_mean += [line.split()[19]]

					par4_allmean += [line.split()[20]]
					par4_varmean += [line.split()[21]]
					par4_allvar  += [line.split()[22]]
					par4_rubgel  += [line.split()[23]]

					Ch1_par5_mean += [line.split()[29]]
					Ch2_par5_mean += [line.split()[30]]
					Ch3_par5_mean += [line.split()[31]]
					Ch4_par5_mean += [line.split()[32]]
					Ch5_par5_mean += [line.split()[33]]

					par5_allmean += [line.split()[34]]
					par5_varmean += [line.split()[35]]
					par5_allvar  += [line.split()[36]]
					par5_rubgel  += [line.split()[37]]

					Ch1_frac += [line.split()[38]]
					Ch2_frac += [line.split()[39]]
					Ch3_frac += [line.split()[40]]
					Ch4_frac += [line.split()[41]]
					Ch5_frac += [line.split()[42]]
	
	tbhdu = fits.BinTableHDU.from_columns(
			[fits.Column(name='GRBname', format='7A', array=GRB),
			fits.Column(name='Chain_1_NH_Int_mean', format='E', array=Ch1_par2_mean),
			fits.Column(name='Chain_2_NH_Int_mean', format='E', array=Ch2_par2_mean),
			fits.Column(name='Chain_3_NH_Int_mean', format='E', array=Ch3_par2_mean),
			fits.Column(name='Chain_4_NH_Int_mean', format='E', array=Ch4_par2_mean),
			fits.Column(name='Chain_5_NH_Int_mean', format='E', array=Ch5_par2_mean),
			fits.Column(name='AllChain_NH_Int_mean', format='E', array=par2_allmean),
			fits.Column(name='AllChain_NH_Int_var', format='E', array=par2_allvar),
			fits.Column(name='AllChain_NH_Int_varofmean', format='E', array=par2_varmean),
			fits.Column(name='NH_Int_RubinGelmin', format='E', array=par2_rubgel),
			fits.Column(name='Chain_1_Gamma_mean', format='E', array=Ch1_par4_mean),
			fits.Column(name='Chain_2_Gamma_mean', format='E', array=Ch2_par4_mean),
			fits.Column(name='Chain_3_Gamma_mean', format='E', array=Ch3_par4_mean),
			fits.Column(name='Chain_4_Gamma_mean', format='E', array=Ch4_par4_mean),
			fits.Column(name='Chain_5_Gamma_mean', format='E', array=Ch5_par4_mean),
			fits.Column(name='AllChain_Gamma_mean', format='E', array=par4_allmean),
			fits.Column(name='AllChain_Gamma_var', format='E', array=par4_allvar),
			fits.Column(name='AllChain_Gamma_varofmean', format='E',  array=par4_varmean),
			fits.Column(name='Gamma_RubinGelmin', format='E', array=par4_rubgel),
	 		fits.Column(name='Chain_1_norm_mean', format='E', array=Ch1_par5_mean),
			fits.Column(name='Chain_2_norm_mean', format='E', array=Ch2_par5_mean),
			fits.Column(name='Chain_3_norm_mean', format='E', array=Ch3_par5_mean),
			fits.Column(name='Chain_4_norm_mean', format='E', array=Ch4_par5_mean),
			fits.Column(name='Chain_5_norm_mean', format='E', array=Ch5_par5_mean),
			fits.Column(name='AllChain_norm_mean', format='E', array=par5_allmean),
			fits.Column(name='AllChain_norm_var', format='E', array=par5_allvar),
			fits.Column(name='AllChain_norm_varofmean', format='E', array=par5_varmean),
			fits.Column(name='norm_RubinGelmin', format='E', array=par5_rubgel),
			fits.Column(name='Chain_1_frac_repval', format='E', array=Ch1_frac),
			fits.Column(name='Chain_2_frac_repval', format='E', array=Ch2_frac),
			fits.Column(name='Chain_3_frac_repval', format='E', array=Ch3_frac),
			fits.Column(name='Chain_4_frac_repval', format='E', array=Ch4_frac),
			fits.Column(name='Chain_5_frac_repval', format='E', array=Ch5_frac)])

	if args.Function[1] == "Y":
		tbhdu.writeto(rootdir + "MCMC_XRAb_FITSTable_" + output + ".fits", clobber=True)
	else:
		tbhdu.writeto(rootdir + "MCMC_Z_FITSTable_" + output + ".fits", clobber=True)

#---------------------------------------------------------------------------------------------------------------#
print "Program executed correctly, ending process."
