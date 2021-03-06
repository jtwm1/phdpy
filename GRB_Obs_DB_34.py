#Python3.4 Script for GRB Ground Observatory Database
#v0.10 
#Authors: Joseph T. W. McGuire, University of Leicester, Department of Physics and Astronomy, XROA Group 
#		  Adam B. Higgins, University of Leicester, Department of Physics and Astronomy, XROA Group 
#Contact: jtwm1@leicester.ac.uk
#		  abh13@leicester.ac.uk

#---------------------------------------------------------------------------------------------------------------------#

import requests
from bs4 import BeautifulSoup
import numpy as np
import os
import argparse
import sys
import re
import warnings
from astropy import units as u
from astropy.coordinates import SkyCoord

### Allow Program to be read from the command line + flag up any incorrect
### GRB name formats and get GRB dates ready for scraping
parser = argparse.ArgumentParser(description='GRB Database of Observational Logs')
parser.add_argument('List',nargs='+',help='Input a single GRB or a list of GRBs')
parser.add_argument('Directory',nargs=1,help='Designate a folder to save the WHT Observing logs')
args = parser.parse_args()
grblist = args.List
folderpath = args.Directory

for i in range(len(grblist)):
    if not re.match(r'^[0-9][0-9][01][1-9][0-3][0-9]$|^[0-9][0-9][01][1-9][0-3][0-9][A-Za-z]$', grblist[i]):
        print("GRB names must be six characters in the format YYMMDD or seven with the addition of a letter between A-Z, can be lowercase.")
        sys.exit()

year = [int(i[:2]) for i in grblist]
year = [str(i + 1900 if i >= 70 else i + 2000) for i in year]
month = [i[2:4] for i in grblist]
day = [i[4:6] for i in grblist]    
    
### Gather GRB list from GRBOX
grburl = 'http://www.astro.caltech.edu/grbox/grbox.php?starttime=700101&endtime=161231'
r = requests.get(grburl).text
soup = BeautifulSoup(r,'lxml')
grbtable = soup.findAll('table')[1]
rows = grbtable.findAll('tr',attrs={'class':'normalrow'})
grb_info = {}
c = 0
for row in rows:
    cols = row.findAll('td')
    cols = [ele.text.strip() for ele in cols]
    if cols[1] == '':
        cols[1] = 'N/A'
    if cols[2] == '':
        cols[2] = 'No Extra Comments'
    try:
        coords = SkyCoord(cols[3],cols[4],unit=(u.hourangle, u.deg))
    except ValueError:
        cols[3] = 'N/A'
        cols[4] = 'N/A'
    if cols[5] == '':
        cols[5] = 'N/A'
    if cols[6] == '':
        cols[6] = 'No'
    else:
        cols[6] = 'Yes'
    if cols[7] == '':
        cols[7] = 'No'
    else:
        cols[7] = 'Yes'
    if cols[8] == '':
        cols[8] = 'No'
    else:
        cols[8] = 'Yes'
    grb_info[cols[0].upper()] = {'T90':cols[1],'RA':coords.ra.degree,
                                 'Dec':coords.dec.degree,'Redshift':cols[5],
                                 'X-ray AG':cols[6],'Optical AG':cols[7],
                                 'Radio AG':cols[8],'General Comments':cols[2]}                   
    c += 1     
    
### Create directory for WHT Log Files (if needed and change into it) 
if not os.path.exists(folderpath):
    os.makedirs(folderpath)
os.chdir(folderpath)

### Check URL works and scrape WHT log data and save to file
print('Downloading WHT log files...')
for i in range(0,len(grblist),1):
    filename = 'WHT Observing Log {}-{}-{}.txt'.format(year[i],month[i],day[i])    
    if not os.path.exists(os.path.join(folderpath,filename)
        url = 'http://www.ing.iac.es/astronomy/observing/inglogs.php?tel=wht&year={0}&month={1}&day={2}&avance=previous'.format(year[i],month[i],day[i])
        r = requests.get(url).text

        try:
            r.raise_for_status()
        except Exception as exc:
            print('\nThe following problem occurred: ',exc,'\n')

        soup = BeautifulSoup(r,'lxml')
        header = soup.title.string
        text = ''
        for ele in soup.findAll('pre'):
            text += ele.get_text()
    
        with open(filename,'w') as f:
            f.write(text)
        f.close()       
        print('Downloaded ',filename)
    else:
        print(filename,'already downloaded!')
