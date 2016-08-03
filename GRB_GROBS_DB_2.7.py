#!/usr/bin/env python

#Python2.7 Script for GRB Ground Observatory Database
#v0.4 
#Authors: Joseph T. W. McGuire, University of Leicester, Department of Physics and Astronomy, XROA Group 
#		  Adam B. Higgins, University of Leicester, Department of Physics and Astronomy, XROA Group 
#Contact: jtwm1@leicester.ac.uk
#		  abh13@leicester.ac.uk

#---------------------------------------------------------------------------------------------------------------------#
#Imports
import os, sys, re
import numpy as np
import requests, argparse
from lxml import html
from bs4 import BeautifulSoup, SoupStrainer
from datetime import date, timedelta

#---------------------------------------------------------------------------------------------------------------------#
#User Input
parser = argparse.ArgumentParser(description="GRB Ground Observatory Database")
parser.add_argument("Input", nargs="+", help="Input a single GRB or list of GRBs")
parser.add_argument("Directory", nargs=1, help="Designate a folder to store observing logs")
parser.add_argument("Observatories", nargs="+", choices=["WHT"], help="Select which observatories you wish to scrape, \
					choices are [WHT]")
args = parser.parse_args()

GRBlist       = args.Input
rootDirectory = args.Directory[0]
Obslist       = args.Observatories

#User Input Flags
for i in range(len(GRBlist)):
	if not re.match(r'^[0-9][0-9][01][1-9][0-3][0-9]$|^[0-9][0-9][01][1-9][0-3][0-9][A-Za-z]$', GRBlist[i]):
		print "GRB names must be six characters in the format YYMMDD or seven with the addition of a letter \
				between A-Z, can be lowercase."
		sys.exit()

#---------------------------------------------------------------------------------------------------------------------#
#Creates root and child directories for storing logs
if not os.path.isdir(rootDirectory):
	os.mkdir(rootDirectory)
	print "Creating root directory..."
	print

for i in range(len(Obslist)):
	if not os.path.isdir(os.path.join(rootDirectory, Obslist[i])):
		os.mkdir(os.path.join(rootDirectory, Obslist[i]))
		print "Creating " + Obslist[i] + " child directory..."
		print

#---------------------------------------------------------------------------------------------------------------------#
#Strips GRB list into YY, MM and DD for web scraping
yr = [int(i[:2]) for i in GRBlist]
yr = [str(i + 1900 if i >= 70 else i + 2000) for i in yr]
mm = [i[2:4] for i in GRBlist]
dd = [i[4:6] for i in GRBlist]

#---------------------------------------------------------------------------------------------------------------------#
#Scrapes GRB information from GRB Online Index by Dr. Dan Perley
print
print "Grabbing GRB information from GRB Online Index (GRBOX) by Dr. Dan Perley..."
print

#Turns GRBs into date objects and increments by 1 day for refined search  
GRBdate    = [date(year=int(i), month=int(j), day=int(k)) for i, j, k in zip(yr, mm, dd)]
searchdate = [i + timedelta(days=1) for i in GRBdate]

starttime  = [re.sub('[A-Za-z]', '', i) for i in GRBlist]
endtime    = [i.strftime("%y%m%d") for i in searchdate]

GRBOXoutput = []
#Begins web scraping
for i in range(len(GRBlist)):
	#Url for GRBOX
	GRBOXurl = "http://www.astro.caltech.edu/grbox/grbox.php?starttime={}&endtime={}".format(starttime[i], endtime[i])

	#Grabs html link
	page = requests.get(GRBOXurl)

	#Gives an exception for a broken url link
	try:
		page.raise_for_status()
	except Exception as exc:
		print
		print "There was a problem: {}".format(exc)
		print 

	#Parses html contents and pulls GRB information
	GRBOXsoup = BeautifulSoup(page.content, "lxml", parse_only=SoupStrainer("tr", {"class": "normalrow"}))

	if len(GRBOXsoup) == 1:
		print "GRB", GRBlist[i], "was not found."
		continue
	
	GRBOXoutput += ["".join(i.findAll(text=True)) for i in GRBOXsoup.findAll("td")]

#Converts the unicode characters to ascii
GRBOXoutput = [i.encode("ascii", "ignore") for i in GRBOXoutput]

#---------------------------------------------------------------------------------------------------------------------#
#Downloads observing logs for GRBlist
print
print "Downloading observing logs..."
print

for i in range(len(GRBlist)):
	for j in range(len(Obslist)):
		#Scrapes WHT observing logs when selected
		if Obslist[j] == "WHT":
			if not os.path.isfile(os.path.join(rootDirectory, "WHT", 
							os.path.basename("WHT Observing Log {}-{}-{}.txt".format(yr[i], mm[i], dd[i])))):
				url = "http://www.ing.iac.es/astronomy/observing/inglogs.php?tel=wht&year={}&month={}&day={}&avance=next".format(yr[i], mm[i], dd[i])

				#Gets html link
				page = requests.get(url)

				#Gives an exception for a broken url link
				try:
					page.raise_for_status()
				except Exception as exc:
					print
					print "There was a problem: {}".format(exc)
					print

				#Parses html contents and separates out relevant log information
				WHTsoup = BeautifulSoup(page.content, "lxml", parse_only=SoupStrainer(["title", "pre"]))
				header = WHTsoup.title.string
				data   = WHTsoup.select("pre")

				#Accounts for html change in WHT logs and then writes logs to txt file
				if len(data[0]) == 1:
					dataFile = open(os.path.join(rootDirectory, "WHT", os.path.basename(header + ".txt")), "wb")
					dataFile.write(data[0].getText().encode('utf-8'))
					dataFile.close()
				else:
					dataFile = open(os.path.join(rootDirectory, "WHT", os.path.basename(header + ".txt")), "wb")
					dataFile.write(data[1].getText().encode('utf-8'))
					dataFile.close()			

				print "Downloaded:", header
			else:
				print "WHT Observing Log {}-{}-{}".format(yr[i], mm[i], dd[i]), "already exists and was not downloaded"

#---------------------------------------------------------------------------------------------------------------------#
#Grabs relevant information 
print
print "Collating GRB information..."
print

for i in range(len(GRBlist)):
	#Grabs data from GRBOX
	print "GRB", GRBlist[i]
	print
	
	for j in GRBOXoutput:
		if j == GRBlist[i]:
			index = GRBOXoutput.index(j)
			print "T90:", GRBOXoutput[index+1], "Redshift:", GRBOXoutput[index+5]
			print "Comments:", GRBOXoutput[index+2]
			print "RA:", GRBOXoutput[index+3], "DEC:", GRBOXoutput[index+4]
			print "XAG:", GRBOXoutput[index+6], "OAG:", GRBOXoutput[index+7], "RAG:", GRBOXoutput[index+8]

	print
	#Grabs data from each selected observatory
	for k in Obslist:
		if k == "WHT":
			with open(os.path.join(rootDirectory, "WHT", 
						"WHT Observing Log " + yr[i] + "-" + mm[i] + "-" + dd[i] + ".txt"), "r") as logFile:
				lostTime = [line.strip().split()[3] for line in logFile if "TIME LOST" in line.strip()]

			print "On date of", GRBlist[i], "for the WHT:"
			print "Time lost due to: Weather", lostTime[0], "Technical", lostTime[1], "Other", lostTime[2]
			print

	print
	print
	

#---------------------------------------------------------------------------------------------------------------------#
#End of program
print
print "Program executed correctly, ending process."
print
