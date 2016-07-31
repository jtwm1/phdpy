#!/usr/bin/env python

#Python2.7 Script for GRB Ground Observatory Database
#v0.30 
#Authors: Joseph T. W. McGuire, University of Leicester, Department of Physics and Astronomy, XROA Group 
#		  Adam B. Higgins, University of Leicester, Department of Physics and Astronomy, XROA Group 
#Contact: jtwm1@leicester.ac.uk
#		  abh13@leicester.ac.uk

#---------------------------------------------------------------------------------------------------------------------#
#Imports
import os, sys, re, time
import numpy as np
import requests, argparse
from lxml import html
from bs4 import BeautifulSoup, SoupStrainer

#---------------------------------------------------------------------------------------------------------------------#
#User Input
parser = argparse.ArgumentParser(description="GRB Ground Observatory Database")
parser.add_argument("Input", nargs="+", help="Single or list of GRBs")
args = parser.parse_args()

GRBlist = args.Input

#User Input Flags
for i in range(len(GRBlist)):
	if not re.match(r'^[0-9][0-9][01][1-9][0-3][0-9]$|^[0-9][0-9][01][1-9][0-3][0-9][A-Za-z]$', GRBlist[i]):
		print("GRB names must be six characters in the format YYMMDD or seven with the addition of a letter between A-Z, can be lowercase.")
		sys.exit()

#---------------------------------------------------------------------------------------------------------------------#
#Creates directory for storing logs
if not os.path.isdir("WHT"):
	os.mkdir("WHT")
	print("Creating WHT directory...")
	print("")

#---------------------------------------------------------------------------------------------------------------------#
#Strips GRB list for web scraping
datelist = [re.sub(r'[A-Za-z]', '', i) for i in GRBlist]

yr = [int(i[:2]) for i in datelist]
yr = [str(i + 1900 if i >= 70 else i + 2000) for i in yr]
mm = [i[2:4] for i in datelist]
dd = [i[4:] for i in datelist]

#---------------------------------------------------------------------------------------------------------------------#
#Downloads observing logs of WHT for GRBlist
print("")
print("Downloading observing logs for WHT...")
print("")

for i in range(len(GRBlist)):
	if os.path.isfile(os.path.join("WHT", 
		os.path.basename("WHT Observing Log {}-{}-{}.txt".format(yr[i], mm[i], dd[i])))) == False:
		url = "http://www.ing.iac.es/astronomy/observing/inglogs.php?tel=wht&year={}&month={}&day={}&avance=next".format(yr[i], mm[i], dd[i])
	
		#Gets html link
		page = requests.get(url)

		#Gives an exception for a broken url link
		try:
			page.raise_for_status()
		except Exception as exc:
			print("")
			print("There was a problem: {}".format(exc))
			print("")

		#Parses html contents and separates out relevant log information
		obslog = BeautifulSoup(page.content, "lxml", parse_only=SoupStrainer(["title", "pre"]))
		header = obslog.title.string
		data   = obslog.select("pre")

		#Accounts for html change in WHT logs and then writes logs to txt file
		if len(data[0]) == 1:
			dataFile = open(os.path.join("WHT", os.path.basename(header + ".txt")), "wb")
			dataFile.write(data[0].getText().encode('utf-8'))
			dataFile.close()
		else:
			dataFile = open(os.path.join("WHT", os.path.basename(header + ".txt")), "wb")
			dataFile.write(data[1].getText().encode('utf-8'))
			dataFile.close()			

		print("Downloaded: " + header)
	else:
		print("WHT Observing Log {}-{}-{}".format(yr[i], mm[i], dd[i]) + " already exists and was not downloaded")

#---------------------------------------------------------------------------------------------------------------------#
#Grabs relelvant information from observing logs
print("")
print("Grabbing WHT log information...")
print("")

lostTime = []

for i in range(len(GRBlist)):
	with open(os.path.join("WHT", "WHT Observing Log " + yr[i] + "-" + mm[i] + "-" + dd[i] + ".txt"), "r") as logFile:
		lostTime = [line.strip().split()[3] for line in logFile if "TIME LOST" in line.strip()]

	print("On date of " + GRBlist[i])
	print("Time lost due to: Weather " + lostTime[0] + " Technical " + lostTime[1] + " Other " + lostTime[2]) 
	print("")
#---------------------------------------------------------------------------------------------------------------------#
#End of program
print("")
print("Program executed correctly, ending process.")
print("")
