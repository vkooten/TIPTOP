# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Created on Tue Feb  7 14:06:58 2023

@author: Vankootenm
"""

from tiptop.tiptop import *
import numpy as np
import pandas as pd
import configparser
import shutil
import hcipy
def find_ngs(data,target_location):
    observing_zenith=target_location[0]
    observing_azimuth=target_location[1]
    #calculate the distance between two things in zenith and azimuth coordinates (polar coordinates)
    NGS=data[np.sqrt((observing_zenith/3600)**2+(data['zenith']/3600)**2+(observing_zenith/3600)*(data['zenith']/3600)*(abs(observing_azimuth-data['azimuth'])))<=1/60]
    ngs_zenith=NGS['zenith']-observing_zenith
    ngs_azimuth=NGS['azimuth']-observing_azimuth
    ngs_locations=[ngs_zenith,ngs_azimuth]
    return NGS, ngs_locations
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, np.rad2deg(phi))
config = configparser.RawConfigParser()
config.optionxform = str
rc("text", usetex=False)

path_config="C:\\Users\\Vankootenm\\Documents\\Python Scripts\\TIPTOP\\GNAO"





#%%%%%%%%%%%%%%%%%%%%#
if os.path.exists(path_config+"\\"+'prepped_starlist.csv'):
    data=pd.read_csv(path_config+"\\"+'prepped_starlist.csv')
else:
    #prep star catalogue
    catalogue= 'trilegal_l0_bm90_9deg2_UBVRIJHKLMN.dat'
    raw_data=pd.read_fwf(path_config+"\\"+catalogue)
    r_0=1.35e10
    max_mag=21
    data=raw_data[raw_data['R']<max_mag]
    
    #create new parameters we need
    zenith=np.random.uniform(0,3,len(data['R']))*3600 #just 3 degrees and we can then offset zenith to place anywhere on-sky
    azimuth=np.random.uniform(0,360,len(data['R'])) ##just 3 degrees and we can then offset azimuth to place anywhere on-sky
    r_flux=r_0*10**(-0.4*data['R'])
    data['R_flux']=r_flux
    data['zenith']=zenith
    data['azimuth']=azimuth
    
    #write file
    data.to_csv(path_config+"\\"+'prepped_starlist.csv')

name="GNAO_Zen25_1NGS"
n_observing=1
observing_zenith=np.random.uniform(0,3,n_observing)*3600 #arcseconds
observing_azimuth=np.random.uniform(0,360,n_observing) #degrees

#lets get our NGS stars
#1arcmin=0.0166667 degrees
NGS,ngs_locations=find_ngs(data,[observing_zenith, observing_azimuth])

#data[np.sqrt(((observing_zenith-data['zenith'])/3600)**2+(observing_azimuth-data['azimuth'])**2)<=1/60]


#setup for wide field mode which is 2'x2' or 0.0333333 x 0.0333333 degrees
#sample uniform grid across the field
n=60
evaluation_locations_zenith=np.linspace(0,85, n) #in arcseconds
evaluation_locations_azimuth=np.random.uniform(0,360, n) #degree


#lets setup the 4 LGS locations now for our system
#need to be in a circle with diameter=85" for wide field
delta=(85/2)
zenith=np.array([delta, delta, delta, delta]).flatten()
azimuth=np.array([45,135,225,315]).flatten()

num_NGS=1
count_ngs_fails=0
HO_do=True
if HO_do:
    #prep configuration files
    shutil.copy(path_config+"\\"+'GNAO_HO.ini',path_config+"\\"+name+'_HO.ini')
    config = configparser.RawConfigParser()
    config.optionxform = str

    config.read(path_config+"\\"+name+'_HO.ini')
    
    config.remove_option('sources_science', 'Zenith')
    config.remove_option('sources_science', 'Azimuth')
    config.set('sources_science', 'Zenith', np.array2string(evaluation_locations_zenith,separator=',', precision=5,max_line_width=1000))
    config.set('sources_science', 'Azimuth', np.array2string(evaluation_locations_azimuth,separator=',', precision=5,max_line_width=1000))

    config.remove_option('DM', 'OptimizationZenith')
    config.remove_option('DM', 'OptimizationAzimuth')
    config.remove_option('DM', 'OptimizationWeight')
    config.set('DM', 'OptimizationZenith', np.array2string(evaluation_locations_zenith, separator=',', precision=5,max_line_width=1000))
    config.set('DM', 'OptimizationAzimuth', np.array2string(evaluation_locations_azimuth,separator=',', precision=5,max_line_width=1000))
    config.set('DM', 'OptimizationWeight', str(n*[1]))

    config.remove_option('sources_HO', 'Zenith')
    config.remove_option('sources_HO', 'Azimuth')
    config.set('sources_HO', 'Zenith', np.array2string(zenith,separator=',', precision=5,max_line_width=1000))
    config.set('sources_HO', 'Azimuth', np.array2string(azimuth,separator=',', precision=5,max_line_width=1000))

    with open(path_config+"\\"+ name+'_HO.ini', 'w') as configfile:
        config.write(configfile)
    configfile.close()
    
    
    #def overallSimulation(path, parametersFile, outputDir, outputFile, doConvolve=False,doPlot=False, verbose=False, returnRes=False, addSrAndFwhm=False):
    file=name+'_HO'
    print(path_config+"\\"+file)
    overallSimulation(path_config, file, path_config, name+'_HO_Results', doPlot=True, doConvolve=True,addSrAndFwhm=True)

if len(NGS)<=num_NGS:
    print(NGS)
    print('Less than required NGS found...quitting')
    count_ngs_fails+=1
    pass
else:
    print('Sufficient NGS found...proceeding')

    #prep configuration files
    shutil.copy(path_config+"\\"+'GNAO_LO.ini',path_config+"\\"+name+'LO_'+str(i)+'.ini')
    config = configparser.RawConfigParser()
    config.optionxform = str

    config.read(path_config+"\\"+name+'HO.ini')
    
    config.remove_option('sources_science', 'Zenith')
    config.remove_option('sources_science', 'Azimuth')
    config.set('sources_science', 'Zenith', np.array2string(evaluation_locations_zenith,separator=',', precision=5,max_line_width=1000))
    config.set('sources_science', 'Azimuth', np.array2string(evaluation_locations_azimuth,separator=',', precision=5,max_line_width=1000))

    config.remove_option('DM', 'OptimizationZenith')
    config.remove_option('DM', 'OptimizationAzimuth')
    config.remove_option('DM', 'OptimizationWeight')
    config.set('DM', 'OptimizationZenith', np.array2string(evaluation_locations_zenith, separator=',', precision=5,max_line_width=1000))
    config.set('DM', 'OptimizationAzimuth', np.array2string(evaluation_locations_azimuth,separator=',', precision=5,max_line_width=1000))
    config.set('DM', 'OptimizationWeight', str(n*[1]))

    config.remove_option('sources_HO', 'Zenith')
    config.remove_option('sources_HO', 'Azimuth')
    config.set('sources_HO', 'Zenith', np.array2string(zenith,separator=',', precision=5,max_line_width=1000))
    config.set('sources_HO', 'Azimuth', np.array2string(azimuth,separator=',', precision=5,max_line_width=1000))

    config.remove_option('sources_LO', 'Zenith')
    config.remove_option('sources_LO', 'Azimuth')
    config.set('sources_HO', 'Zenith', np.array2string(ngs_locations[0],separator=',', precision=5,max_line_width=1000))
    config.set('sources_HO', 'Azimuth', np.array2string(ngs_locations[1],separator=',', precision=5,max_line_width=1000))

    with open(path_config+"\\"+name+'LO_'+str(i)+'.ini', 'w') as configfile:
        config.write(configfile)
    configfile.close()
    
    
    #def overallSimulation(path, parametersFile, outputDir, outputFile, doConvolve=False,doPlot=False, verbose=False, returnRes=False, addSrAndFwhm=False):
    file=name+'LO_'+str(i)
    overallSimulation(path_config, file, path_config, name+'_Results', doPlot=True, doConvolve=True,addSrAndFwhm=True)


