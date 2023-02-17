# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 07:25:08 2023

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

def plot_contour(data,pos,lgs_zenith,str_type,title):
    levels=np.array([np.percentile(data,25),np.percentile(data,50),np.percentile(data,75)])
    #level=np.flipud(levels)
    fig, (ax2) = plt.subplots(nrows=1)
    import matplotlib.cm as cm # matplotlib's color map library

    cpf = ax2.tricontourf(pos[0,:], pos[1,:], data,  cmap=cm.Reds)
    line_colors = ['black' for l in cpf.levels]
    linestyle = ['--' for l in cpf.levels]
    cp=ax2.tricontour(pos[0,:], pos[1,:], data, levels=levels,colors=line_colors,linestyle=linestyle)
    fmt = {}
    if str_type=='FWHM':
        strs = [75, 50, 25]
    else:
        strs = [25, 50, 75]

    for l, s in zip(cp.levels, strs):
        fmt[l] = s
    #plt.scatter(pos[0,:], pos[1,:],marker='*',s=50)

    #ax2.tricontour(X, Y, data, levels=3, linewidths=0., colors='k')
    ax2.clabel(cp,cp.levels,inline=2,fmt=fmt, fontsize=10, colors=line_colors)
    fig.colorbar(cpf, ax=ax2)
    ax2.plot(pos[0,:], pos[1,:], 'ko', ms=3)
    #ax2.set(xlim=(-40, 40), ylim=(-40, 40))
    #ax2.set_title('tricontour (%d points)' % npts)
    lgs_zenith = np.array([lgs_zenith, lgs_zenith, lgs_zenith, lgs_zenith])
    lgs_azimuth =np.array([45, 135, 225, 315])
    lgs_wide_x=lgs_zenith*np.cos(np.radians(lgs_azimuth))
    lgs_wide_y=lgs_zenith*np.sin(np.radians(lgs_azimuth))
    plt.scatter(lgs_wide_x,lgs_wide_y,marker='*',s=50)
    plt.subplots_adjust(hspace=0.5)
    plt.xlabel('')
    plt.text(-32,-38,"25-pctl: %.2f  50-pctl: %.2f 75-pctl: %.2f"% (cp.levels[0],cp.levels[1],cp.levels[2]),fontsize=10,color='white')
    plt.suptitle(str_type+' '+title)
    plt.savefig(path_config+"\\"+name+'_'+title+'_'+str_type+'.png',dpi=700)
    plt.show()
config = configparser.RawConfigParser()
config.optionxform = str
rc("text", usetex=False)

path_config="C:\\Users\\Vankootenm\\Documents\\Python Scripts\\TIPTOP\\GNAO"
######edit OOMAO HO simulation science fields first
#%%%%%%%%%%%%%%%%%%%%#
catalogue= 'OOMAO_sciencefield.csv'
sci_field=pd.read_csv(path_config+"\\"+catalogue,encoding= 'unicode_escape')
config = configparser.RawConfigParser()
config.optionxform = str
name='GNAO_HO_OOMAO'
config.read(path_config+"\\"+name+'.ini')

config.remove_option('sources_science', 'Zenith')
config.remove_option('sources_science', 'Azimuth')
config.set('sources_science', 'Zenith', np.array2string(np.array(sci_field['zen[arcsec]'].values,dtype=np.float64),separator=',', precision=5,max_line_width=1000))
config.set('sources_science', 'Azimuth', np.array2string(np.array(sci_field['azim[deg]'].values,dtype=np.float64),separator=',', precision=5,max_line_width=1000))
config.remove_option('DM', 'OptimizationZenith')
config.remove_option('DM', 'OptimizationAzimuth')
config.remove_option('DM', 'OptimizationWeight')
config.set('DM', 'OptimizationZenith', np.array2string(np.array(sci_field['zen[arcsec]'].values,dtype=np.float64), separator=',', precision=5,max_line_width=1000))
config.set('DM', 'OptimizationAzimuth', np.array2string(np.array(sci_field['azim[deg]'].values,dtype=np.float64),separator=',', precision=5,max_line_width=1000))
config.set('DM', 'OptimizationWeight', str(len(sci_field)*[1]))



with open(path_config+"\\"+name+'.ini', 'w') as configfile:
        config.write(configfile)
configfile.close()
file=name
sr, fwhm, pos = overallSimulation(path_config, file, path_config, name+'_Results', doPlot=True, doConvolve=True,addSrAndFwhm=True)
plot_contour(sr,pos,'SR','HO_widefield')
plot_contour(fwhm,pos,'FWHM','HO_widefield')

######## LO for wide field - NOT WORKING ####
catalogue= 'OOMAO_sciencefield.csv'
sci_field=pd.read_csv(path_config+"\\"+catalogue,encoding= 'unicode_escape')
sci_field=sci_field[0:8:2]
config = configparser.RawConfigParser()
config.optionxform = str
name='GNAO_LO_OOMAO'
config.read(path_config+"\\"+name+'.ini')

config.remove_option('sources_science', 'Zenith')
config.remove_option('sources_science', 'Azimuth')
config.set('sources_science', 'Zenith', np.array2string(np.array(sci_field['zen[arcsec]'].values,dtype=np.float64),separator=',', precision=5,max_line_width=1000))
config.set('sources_science', 'Azimuth', np.array2string(np.array(sci_field['azim[deg]'].values,dtype=np.float64),separator=',', precision=5,max_line_width=1000))
config.remove_option('DM', 'OptimizationZenith')
config.remove_option('DM', 'OptimizationAzimuth')
config.remove_option('DM', 'OptimizationWeight')
config.set('DM', 'OptimizationZenith', np.array2string(np.array(sci_field['zen[arcsec]'].values,dtype=np.float64), separator=',', precision=5,max_line_width=1000))
config.set('DM', 'OptimizationAzimuth', np.array2string(np.array(sci_field['azim[deg]'].values,dtype=np.float64),separator=',', precision=5,max_line_width=1000))
config.set('DM', 'OptimizationWeight', str(len(sci_field)*[1]))

#now lets do the 3 NGS case assuming 12 mag 
Hz=1000
r_ngs=(9.65e9*10**(-0.4*12))*0.8*(20*0.395)**2*(1./Hz)
config.remove_option('sensor_LO', 'NumberPhotons')
config.set('sensor_LO', 'NumberPhotons', np.array2string(np.array(1*[r_ngs]),separator=',', precision=5,max_line_width=1000))

with open(path_config+"\\"+name+'.ini', 'w') as configfile:
        config.write(configfile)
configfile.close()
file=name
sr, fwhm, pos = overallSimulation(path_config, file, path_config, name+'_Results', doPlot=True, doConvolve=True,addSrAndFwhm=True)

######## LO fornarrow field - NOT WORKING ####
catalogue= 'OOMAO_sciencefield.csv'
sci_field=pd.read_csv(path_config+"\\"+catalogue,encoding= 'unicode_escape')
sci_field=sci_field[0:8:2]
config = configparser.RawConfigParser()
config.optionxform = str
name='GNAO_LO_NF_OOMAO'
config.read(path_config+"\\"+name+'.ini')

config.remove_option('sources_science', 'Zenith')
config.remove_option('sources_science', 'Azimuth')
x=np.linspace(-10,10,5)
y=np.linspace(-10,10,5)
xv, yv = np.meshgrid(x, y)
z,t=cart2pol(xv.ravel(), yv.ravel())
config.set('sources_science', 'Zenith', np.array2string(z,separator=',', precision=5,max_line_width=1000))
config.set('sources_science', 'Azimuth', np.array2string(t,separator=',', precision=5,max_line_width=1000))
config.remove_option('DM', 'OptimizationZenith')
config.remove_option('DM', 'OptimizationAzimuth')
config.remove_option('DM', 'OptimizationWeight')
config.set('DM', 'OptimizationZenith', np.array2string(z, separator=',', precision=5,max_line_width=1000))
config.set('DM', 'OptimizationAzimuth', np.array2string(t,separator=',', precision=5,max_line_width=1000))
config.set('DM', 'OptimizationWeight', str(len(z)*[1]))

#now lets do the 3 NGS case assuming 12 mag 
Hz=1000
r_ngs=(9.65e9*10**(-0.4*12))*0.8*(20*0.395)**2*(1./Hz)
config.remove_option('sensor_LO', 'NumberPhotons')
config.set('sensor_LO', 'NumberPhotons', np.array2string(np.array(3*[r_ngs]),separator=',', precision=5,max_line_width=1000))

with open(path_config+"\\"+name+'.ini', 'w') as configfile:
        config.write(configfile)
configfile.close()
file=name
sr, fwhm, pos = overallSimulation(path_config, file, path_config, name+'_Results', doPlot=True, doConvolve=True,addSrAndFwhm=True)
plot_contour(sr,pos,10,'SR','LO_narrowfield')
plot_contour(fwhm,pos,10,'FWHM','LO_narrowfield')