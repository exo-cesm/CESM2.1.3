#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Thu Jul  9 09:30:19 2020

@author: Greg Cooke (based on https://github.com/nasapsg/globes exo_cam example)

See PSG handbook for more detailed information
https://psg.gsfc.nasa.gov/help.php#handbook

"""

#%% Initial imports

import xarray as xr
import sys
import numpy as np
from netCDF4 import Dataset as ncdf


'''
Instructions

Choose planetary and stellar parameters
Choose telescope and observation type
'''



#%% Initial inputs to define PSG setup
Object = 'Exoplanet'

'''
Star
'''
Star = 'TRAPPIST-1'#'Sun'

'''
Planetary parameter
'''
planet_distance = 10 #distance in pc
season = 180 # angle around orbit

'''
Upload locally or to website
Requires docker install and script changes
But only necessary if producing hundreds of PSG simulations

'''
upload = False 

'''
Choose telescope - may need updating

LUVOIR HDI, LUVOIR A-UV, LUVOIR A-VIS, LUVOIR A-NIR
LUVOIR B-UV, LUVOIR B-VIS, LUVOIR B-NIR
JWST NIRSpec, JWST MIRI
Keck HIRES, Keck NIRSpec
'''

Telescope_choice = 'Ideal' 
transit = True
albedo = False
albedo_test = False

'''
Choose observation output units
'''

obs_type = 'Jy'
obs_type = 'Wsrm2um'
if (albedo == True):
    obs_type = 'rif' #geometric albedo
    
'''
Observation exposure time
Only important if considering noise sources
'''    
exposure_time = 60 #in seconds
exposure_number = 1440

#GCM binning resolution
Binning = '10' 

'''
Specific parameters required for transit
'''
if (transit == True):
    planet_distance = 10
    Binning = '3'
    exposure_time = 60 #in seconds
    exposure_number = 780 
    if (Star == 'TRAPPIST-1'):
        exposure_time = 60 #in seconds
        exposure_number = 53 
    obs_type = 'rkm'
    season = 180

# Total exposure time    
total_T = exposure_time * exposure_number / (60**2)

'''
Specific parameters required for TRAPPIST-1 e transit
'''
if (Star == 'TRAPPIST-1'):
    total_T = round(total_T, 3)
    planet_distance = 12.47

'''
Selection of the number of stream pairs (NMAX) 
and Legendre terms (LMAX) 
Changing these changes simulation speed and accuracy
'''
LMAX = '41'
NMAX = '3'

'''
Cloud properties and whether they are included or not
First size is liquid cloud size, 
second size is ice cloud size.
'''
#include clouds or not (ice and liquid)
clouds = True
cloud_sizes = '5,100'
liquid_cloud = True
ice_cloud = True

'''
if True, Nitrogen only in atmosphere 
(removes other gas phase chemicals)
'''
N2_atmosphere = False

#%% Set star properties based on user input - other stars can be added

if (Star == 'Sun'):
    #Values for Earth
    stellar_type = 'G'
    stellar_temp = '5777'
    stellar_radius = '1'
    Metallicity = '0.0'
    period = '365'
    eccentricity = '0.01670'
    diameter = '12742'
    gravity = '9.8'
    
elif (Star == 'TRAPPIST-1'):
    #Values for TRAPPIST-1 e
    stellar_type = 'M'
    stellar_temp = '2566'
    stellar_radius = '0.1192'
    Metallicity = '0.04'
    period = '6.099'
    eccentricity = '0.00510'
    diameter = '11595'
    gravity = '9.15'

#%%  Wavelength ranges and spectral resolving power for instruments

HIRES_wav_start = '0.3'
HIRES_wav_end = '1.0'
HIRES_R = '50000'

Ideal_wav_start = '0.1'
Ideal_wav_end = '18'
Ideal_R = '250'

#%% Initialise some variables 

periapsis = 0
date = 'Dec02'
dm = '01-01'
sollon = '180.0'
sollat = '0.0'

#%% now define function to upload to PSG
def Convert_GCM(infile = 'b.e21.BWma1850.f19_g17.150pc_o2.001.cam.h1.0026-12-27-00000.nc',
                infile_swap = 'b.e21.BWma1850.f19_g17.150pc_o2.001.cam.h1.0026-12-27-00000.nc',
                fileout = 'Dec_27th_150pc_O2_psg.gcm', phase = 170, obs_type = 'ppm',
                sollon = '0', sollat = '0', date = '01-01', distance = '1', velocity = '0', season = '0',
                name = 'Baseline', Telescope = 'LUVOIR A-UV'):
    
    new_netcdf_name = 'Compressed_lon.nc'
    
    ncFile = xr.open_dataset(infile,decode_times=False) #open the file and decode time as false
    ncFile = ncFile.isel(time = np.arange(0,1))
    #compress longitude to fit upload limit
    avg_lon = 4
    lon_new = ncFile.isel(lon = np.arange(0,144,avg_lon))
    
    for i in range(0,int(144/avg_lon),1):
        lon_new.TS[:,:,i] = ncFile.TS.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.PS[:,:,i] = ncFile.PS.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.SSAVIS[:,:,i] = ncFile.SSAVIS.isel(lon =  np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.SSAVISdn[:,:,i] = ncFile.SSAVISdn.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.ALDIF[:,:,i] = ncFile.ALDIF.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.ALDIR[:,:,i] = ncFile.ALDIR.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.ASDIF[:,:,i] = ncFile.ASDIF.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.ASDIR[:,:,i] = ncFile.ASDIR.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.CH4[:,:,:,i] = ncFile.CH4.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.OH[:,:,:,i] = ncFile.OH.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.CLDICE[:,:,:,i] = ncFile.CLDICE.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.CLDLIQ[:,:,:,i] = ncFile.CLDLIQ.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.CO[:,:,:,i] = ncFile.CO.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.CO2[:,:,:,i] = ncFile.CO2.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.H2O[:,:,:,i] = ncFile.H2O.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.H2O2[:,:,:,i] = ncFile.H2O2.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.HNO3[:,:,:,i] = ncFile.HNO3.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.N2O[:,:,:,i] = ncFile.N2O.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.NO[:,:,:,i] = ncFile.NO.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.NO2[:,:,:,i] = ncFile.NO2.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.O[:,:,:,i] = ncFile.O.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.O1D[:,:,:,i] = ncFile.O1D.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.O2[:,:,:,i] = ncFile.O2.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
        lon_new.O3[:,:,:,i] = ncFile.O3.isel(lon = np.arange(i*avg_lon,i*avg_lon+avg_lon,1)).mean(dim ='lon')
 
    lon_new.to_netcdf(new_netcdf_name)
    
    itime = 0
    nfile = ncdf(new_netcdf_name)
    lat = (nfile.variables['lat'])[:]
    lon = (nfile.variables['lon'])[:];
    hyam = (nfile.variables['hyam'])[:];     hyam = np.flipud(hyam)
    hybm = (nfile.variables['hybm'])[:];     hybm = np.flipud(hybm)
    P0 = (nfile.variables['P0'])[:];         P0 = P0*1e-5
    PS = (nfile.variables['PS'])[:];         PS = PS[itime,:,:]*1e-5
    U = (nfile.variables['U'])[:];           U = U[itime,::-1,:,:]
    V = (nfile.variables['V'])[:];           V = V[itime,::-1,:,:]
    O3 = (nfile.variables['O3'])[:];           O3 = O3[itime,::-1,:,:]
    O2 = (nfile.variables['O2'])[:];           O2 = O2[itime,::-1,:,:]
    T = (nfile.variables['T'])[:];           T = T[itime,::-1,:,:]
    TS = (nfile.variables['TS'])[:];         TS = TS[itime,:,:]
    ASDIR = (nfile.variables['ASDIR'])[:];   ASDIR = ASDIR[itime,:,:]
    CLDICE = (nfile.variables['CLDICE'])[:]; CLDICE = CLDICE[itime,::-1,:,:]# + 1e-30
    CLDLIQ = (nfile.variables['CLDLIQ'])[:]; CLDLIQ = CLDLIQ[itime,::-1,:,:]# + 1e-30
    if (clouds == False): CLDLIQ = CLDLIQ * 0 + 1e-30; CLDICE = CLDICE * 0 + 1e-30
    CH4 = (nfile.variables['CH4'])[:];    CH4 = CH4[itime,::-1,:,:]
    CO2 = (nfile.variables['CO2'])[:];    CO2 = CO2[itime,::-1,:,:]
    N2O = (nfile.variables['N2O'])[:];    N2O = N2O[itime,::-1,:,:]
    H2O = (nfile.variables['H2O'])[:];         H2O = H2O[itime,::-1,:,:]
    OH = (nfile.variables['OH'])[:];         OH = OH[itime,::-1,:,:]
    HNO3 = (nfile.variables['HNO3'])[:];         HNO3 = HNO3[itime,::-1,:,:]
    O = (nfile.variables['O'])[:];         O = O[itime,::-1,:,:]
    nfile.close()
    
    # Fix variables
    ASDIR = np.where((ASDIR>=0) & (ASDIR<=1.0) & (np.isfinite(ASDIR)), ASDIR, 0.3)
    ASDIR = np.where((ASDIR>=0.981), 0, ASDIR)
    if (albedo_test == True):
        np.where(((ncFile.lat >= 70 or ncFile.lat <= -70) and  ASDIR>0.8), 0.06, ASDIR)
    CLDICE = np.where((CLDICE>0) & (np.isfinite(CLDICE)), CLDICE, 1e-30)
    CLDLIQ = np.where((CLDLIQ>0) & (np.isfinite(CLDLIQ)), CLDLIQ, 1e-30)
    CH4 = np.where((CH4>0) & (np.isfinite(CH4)), CH4, 1e-30)
    CO2 = np.where((CO2>0) & (np.isfinite(CO2)), CO2, 1e-30)
    N2O = np.where((N2O>0) & (np.isfinite(N2O)), N2O, 1e-30)
    H2O = np.where((H2O>0) & (np.isfinite(H2O)), H2O, 1e-30)
    

       
    sz = np.shape(T);
    sz = np.flip(sz);
    press3D = np.zeros((sz), dtype=np.float32); 
    for i in range(sz[0]):
    	for j in range(sz[1]):
    		press3D[i,j,:] = hyam * P0 + hybm * PS[j,i]; 
    N2 = 1.0 - (CH4 + CO2 + N2O + H2O + O2)
    if (N2_atmosphere == True): 
        N2 = N2 * 0 + 1;  
        O3 = O3 * 0 + 1e-30; 
        O2 = O2 * 0 + 1e-30;
        O = O * 0 + 1e-30; 
        if (clouds == True):
            print('Clouds included N2 atmosphere')
        else:
            CLDICE = CLDICE * 0 + 1e-30; 
            CLDLIQ = CLDLIQ * 0 + 1e-30; 
        CH4 = CH4 * 0 + 1e-30; 
        CO2 = CO2 * 0 + 1e-30; 
        N2O = N2O * 0 + 1e-30; 
        H2O = H2O * 0 + 1e-30; 
        OH = OH * 0 + 1e-30; 
        HNO3 = HNO3 * 0 + 1e-30;
    
    # Save object parameters
    newf = []
    newf.append('<OBJECT>'+Object)
    newf.append('<OBJECT-NAME>Earth')
    newf.append('<OBJECT-DATE>2020/'+date[:2]+'/'+date[-2:]+' 00:00')
    newf.append('<OBJECT-DIAMETER>'+diameter)
    newf.append('<OBJECT-GRAVITY>'+gravity)
    newf.append('<OBJECT-GRAVITY-UNIT>g')
    newf.append('<OBJECT-STAR-DISTANCE>'+distance)
    newf.append('<OBJECT-STAR-VELOCITY>'+velocity)
    newf.append('<OBJECT-SOLAR-LONGITUDE>'+sollon)
    newf.append('<OBJECT-SOLAR-LATITUDE>'+sollat)
    if (transit == True):
        newf.append('<OBJECT-SEASON>180')
    else:
        newf.append('<OBJECT-SEASON>'+str(season))
    newf.append('<OBJECT-STAR-TYPE>'+stellar_type)
    newf.append('<OBJECT-STAR-TEMPERATURE>'+stellar_temp)
    newf.append('<OBJECT-STAR-RADIUS>'+stellar_radius)
    newf.append('<OBJECT-OBS-LONGITUDE>'+sollon)
    newf.append('<OBJECT-OBS-LATITUDE>'+sollat)
    newf.append('<OBJECT-STAR-METALLICITY>'+Metallicity)
    newf.append('<OBJECT-PERIOD>'+period)
    newf.append('<OBJECT-PERIAPSIS>'+str(periapsis))
    newf.append('<OBJECT-OBS-VELOCITY>29.0')
    newf.append('<OBJECT-INCLINATION>90.00')
    newf.append('<OBJECT-ECCENTRICITY>'+eccentricity)
    # Save geometry parameters
    newf.append('<GEOMETRY>Observatory')
    newf.append('<GEOMETRY-OFFSET-NS>0.0')
    newf.append('<GEOMETRY-OFFSET-EW>0.0')
    newf.append('<GEOMETRY-OFFSET-UNIT>arcsec')
    newf.append('<GEOMETRY-OBS-ALTITUDE>'+str(planet_distance))
    newf.append('<GEOMETRY-ALTITUDE-UNIT>pc')
    newf.append('<GEOMETRY-USER-PARAM>0.0')
    newf.append('<GEOMETRY-STELLAR-TYPE>G')
    newf.append('<GEOMETRY-STELLAR-TEMPERATURE>5777')
    newf.append('<GEOMETRY-STELLAR-MAGNITUDE>0')
    newf.append('<GEOMETRY-SOLAR-ANGLE>90.000')
    newf.append('<GEOMETRY-OBS-ANGLE>90')
    newf.append('<GEOMETRY-PLANET-FRACTION>1.000e+00')
    if (transit == True): 
        newf.append('<GEOMETRY-STAR-DISTANCE>0')
    else:
        G_S_D = abs((1.49e11*np.sin(season*np.pi/180) / (planet_distance * 3.086e16))/4.8481e-6)
        newf.append('<GEOMETRY-STAR-DISTANCE>'+str(G_S_D))
    if (transit == True):
        newf.append('<GEOMETRY-STAR-FRACTION>'+str(8.38086e-05))
        newf.append('<GEOMETRY-PHASE>180')
    else:
        newf.append('<GEOMETRY-STAR-FRACTION>0')
        newf.append('<GEOMETRY-PHASE>'+str(season))
    newf.append('<GEOMETRY-REF>User')
    newf.append('<GEOMETRY-DISK-ANGLES>1')
    newf.append('<GEOMETRY-ROTATION>-0.00,0.00')
    newf.append('<GEOMETRY-AZIMUTH>0.000')
    newf.append('<GEOMETRY-BRDFSCALER>1.000')
    
    
    # Save atmosphere parameters
    newf.append('<ATMOSPHERE-DESCRIPTION>WACCM6 '+name+' simulation')
    newf.append('<ATMOSPHERE-STRUCTURE>Equilibrium')
    newf.append('<ATMOSPHERE-PRESSURE>1.0')
    newf.append('<ATMOSPHERE-PUNIT>bar')
    newf.append('<ATMOSPHERE-TEMPERATURE>300')
    newf.append('<ATMOSPHERE-WEIGHT>28.97')
    newf.append('<ATMOSPHERE-LAYERS>70')
    newf.append('<ATMOSPHERE-NMAX>'+NMAX)
    newf.append('<ATMOSPHERE-LMAX>'+LMAX)
    newf.append('<ATMOSPHERE-CONTINUUM>Rayleigh,Refraction,CIA_all,UV_all')
    newf.append('<ATMOSPHERE-NGAS>7')
    newf.append('<ATMOSPHERE-GAS>N2,H2O,O2,O3,CO2,N2O,CH4')
    newf.append('<ATMOSPHERE-TYPE>HIT[22],HIT[1],HIT[7],HIT[3],HIT[2],HIT[4],HIT[6]')
    newf.append('<ATMOSPHERE-ABUN>1,1,1,1,1,1,1')
    newf.append('<ATMOSPHERE-UNIT>scl,scl,scl,scl,scl,scl,scl')
    newf.append('<ATMOSPHERE-NAERO>2')
    newf.append('<ATMOSPHERE-AEROS>Water,WaterIce')
    newf.append('<ATMOSPHERE-ATYPE>AFCRL_WATER_HRI[0.20um-0.03m],AFCRL_ICE_HRI[0.20um-0.03m]')
    newf.append('<ATMOSPHERE-AABUN>1,1')
    newf.append('<ATMOSPHERE-AUNIT>scl,scl')
    newf.append('<ATMOSPHERE-ASIZE>'+cloud_sizes)
    newf.append('<ATMOSPHERE-ASUNI>um,um')
    # Save surface parameters
    newf.append('<SURFACE-TEMPERATURE>280')
    newf.append('<SURFACE-ALBEDO>0.3')
    newf.append('<SURFACE-EMISSIVITY>0.8')
    newf.append('<SURFACE-NSURF>0')
    
    if (Telescope == 'LUVOIR HDI'):
        newf.append('<GENERATOR-INSTRUMENT>LUVOIR_HDI: The HDI design provides a 2 x 3 arcminute field-of-view with two channels simultaneously (via beamsplitter), an ultraviolet/visible (UVIS) channel covering the range 0.2-1.0 um and a near-infrared (NIR) channel covering the range 0.8-2.5um. The UVIS channel employs a CMOS detector and the NIR a HgCdTe detector. Dozens of filters and several grisms provide spectroscopic capabilities up to RP=500 in both channels.')
        newf.append('<GENERATOR-RANGE1>0.2')
        newf.append('<GENERATOR-RANGE2>2.5')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>500')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>SINGLE')
        newf.append('<GENERATOR-DIAMTELE>15.0')
        newf.append('<GENERATOR-BEAM>1')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1')
        newf.append('<GENERATOR-TELESCOPE2>2.0')
        newf.append('<GENERATOR-TELESCOPE3>1.0')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>0.2@0.2,0.2@1,2.5@1.01,2.5@2.5')
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.5')
        newf.append('<GENERATOR-NOISEOTEMP>270')
        newf.append('<GENERATOR-NOISEOEFF>0.100@2.25e-1,0.112@2.75e-1,0.175@3.36e-1,0.211@4.75e-1,0.211@6.06e-1,0.211@7.75e-1,0.298@8.00e-1,0.342@1.26e+0,0.342@1.60e+0,0.335@2.22e+0,0.335@2.5')
        newf.append('<GENERATOR-NOISEOEMIS>0.1')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>8')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-TRANS>03-01')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.0')
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)
        
    if (Telescope == 'LUVOIR HDI 6m'):
        newf.append('<GENERATOR-INSTRUMENT>user')
        newf.append('<GENERATOR-RANGE1>0.2')
        newf.append('<GENERATOR-RANGE2>2.5')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>500')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>SINGLE')
        newf.append('<GENERATOR-DIAMTELE>6')
        newf.append('<GENERATOR-BEAM>1')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1')
        newf.append('<GENERATOR-TELESCOPE2>2.0')
        newf.append('<GENERATOR-TELESCOPE3>1.0')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>0.2@0.2,0.2@1,2.5@1.01,2.5@2.5')
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.5')
        newf.append('<GENERATOR-NOISEOTEMP>270')
        newf.append('<GENERATOR-NOISEOEFF>0.100@2.25e-1,0.112@2.75e-1,0.175@3.36e-1,0.211@4.75e-1,0.211@6.06e-1,0.211@7.75e-1,0.298@8.00e-1,0.342@1.26e+0,0.342@1.60e+0,0.335@2.22e+0,0.335@2.5')
        newf.append('<GENERATOR-NOISEOEMIS>0.1')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>8')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-TRANS>03-01')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.0')
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)

    
    if (Telescope == 'LUVOIR A-UV'):
        newf.append('<GENERATOR-INSTRUMENT>LUVOIR_A-UV: The Extreme Coronagraph for Living Planetary Systems (ECLIPS) delivers continuous spectral coverage from 200 nm to 2.5 um via three channels, UV (200 to 525 nm), VIS (515 nm to 1030 nm), and NIR (1 to 2 microns). The UV channel is effectively an imager and provides a maximum resolution of RP=7, while the VIS channel RP=140, and NIR=70. The core coronagraph throughput is practically twice for LUVOIR-B than A.')
        newf.append('<GENERATOR-RANGE1>0.2')
        newf.append('<GENERATOR-RANGE2>0.515')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>7')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>CORONA')
        newf.append('<GENERATOR-DIAMTELE>15.0')
        newf.append('<GENERATOR-BEAM>1')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1E-10')
        newf.append('<GENERATOR-TELESCOPE2>4.5')
        newf.append('<GENERATOR-TELESCOPE3>2.700e-11@0.000,2.700e-11@0.217,2.700e-11@1.410,2.110e-03@2.278,6.329e-03@2.767,1.688e-02@3.092,5.063e-02@3.418,1.097e-01@3.797,1.477e-01@4.014,1.751e-01@4.231,1.920e-01@4.557,2.004e-01@5.371,2.004e-01@6.130,2.131e-01@6.618,2.532e-01@7.215,2.679e-01@9.005,2.700e-01@12.640,2.700e-01@17.631,2.700e-01@20.561,2.700e-01@24.304,2.700e-01@27.288,2.700e-01@29.349')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>0@0.2,0@1,2.5@1.01,2.5@2.0')
        newf.append('<GENERATOR-NOISEOTEMP>270')
        newf.append('<GENERATOR-NOISEOEFF>0.0317@0.2000,0.0437@0.2261,0.0589@0.2580,0.0742@0.2986,0.0851@0.3377,0.0917@0.3667,0.0971@0.4029,0.1015@0.4493,0.1004@0.4971,0.1004@0.5140,0.1670@0.5150,0.1659@0.5377,0.1506@0.6304,0.1255@0.7087,0.0939@0.7986,0.0884@0.8435,0.1146@0.9058,0.1419@0.9594,0.1594@0.9942,0.1821@1.2200,0.1958@1.4100,0.2049@1.6200,0.2094@1.8700,0.2140@2.0000')
        newf.append('<GENERATOR-NOISEOEMIS>0.1')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>10')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-TRANS>03-01')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.0')
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)
        
    if (Telescope == 'LUVOIR B-UV'):
        newf.append('<GENERATOR-INSTRUMENT>LUVOIR_B-UV: The Extreme Coronagraph for Living Planetary Systems (ECLIPS) delivers continuous spectral coverage from 200 nm to 2.5 um via three channels, UV (200 to 525 nm), VIS (515 nm to 1030 nm), and NIR (1 to 2 microns). The UV channel is effectively an imager and provides a maximum resolution of RP=7, while the VIS channel RP=140, and NIR=70. The core coronagraph throughput is practically twice for LUVOIR-B than A.')
        newf.append('<GENERATOR-RANGE1>0.2')
        newf.append('<GENERATOR-RANGE2>0.515')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>7')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>CORONA')
        newf.append('<GENERATOR-DIAMTELE>6.0')
        newf.append('<GENERATOR-BEAM>1')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1E-10')
        newf.append('<GENERATOR-TELESCOPE2>4.5')
        newf.append('<GENERATOR-TELESCOPE3>4.578000e-11@0.000,4.578000e-11@0.216,4.578000e-11@0.649,2.110e-03@0.973,1.266e-02@1.459,8.228e-02@2.108,1.709e-01@2.973,2.658e-01@4.486,3.418e-01@6.757,3.945e-01@10.216,4.219e-01@14.108,4.409e-01@19.459,4.536e-01@23.622,4.578e-01@27.784,4.578e-01@29.459')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>0@0.2,0@1,2.5@1.01,2.5@2.0')
        newf.append('<GENERATOR-NOISEOTEMP>270')
        newf.append('<GENERATOR-NOISEOEFF>0.0317@0.2000,0.0437@0.2261,0.0589@0.2580,0.0742@0.2986,0.0851@0.3377,0.0917@0.3667,0.0971@0.4029,0.1015@0.4493,0.1004@0.4971,0.1004@0.5140,0.1670@0.5150,0.1659@0.5377,0.1506@0.6304,0.1255@0.7087,0.0939@0.7986,0.0884@0.8435,0.1146@0.9058,0.1419@0.9594,0.1594@0.9942,0.1821@1.2200,0.1958@1.4100,0.2049@1.6200,0.2094@1.8700,0.2140@2.0000')
        newf.append('<GENERATOR-NOISEOEMIS>0.1')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>10')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-TRANS>03-01')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.0')
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)
        
    if (Telescope == 'LUVOIR A-UV no coronagraph'):
        print('LUVOIR A-UV no coronagraph')
        newf.append('<GENERATOR-INSTRUMENT>LUVOIR_A-UV: The Extreme Coronagraph for Living Planetary Systems (ECLIPS) delivers continuous spectral coverage from 200 nm to 2.5 um via three channels, UV (200 to 525 nm), VIS (515 nm to 1030 nm), and NIR (1 to 2 microns). The UV channel is effectively an imager and provides a maximum resolution of RP=7, while the VIS channel RP=140, and NIR=70. The core coronagraph throughput is practically twice for LUVOIR-B than A.')
        newf.append('<GENERATOR-RANGE1>0.2')
        newf.append('<GENERATOR-RANGE2>0.515')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>7')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>SINGLE')
        newf.append('<GENERATOR-DIAMTELE>15.0')
        newf.append('<GENERATOR-BEAM>1')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE2>4.5')
        newf.append('<GENERATOR-TELESCOPE3>2.700e-11@0.000,2.700e-11@0.217,2.700e-11@1.410,2.110e-03@2.278,6.329e-03@2.767,1.688e-02@3.092,5.063e-02@3.418,1.097e-01@3.797,1.477e-01@4.014,1.751e-01@4.231,1.920e-01@4.557,2.004e-01@5.371,2.004e-01@6.130,2.131e-01@6.618,2.532e-01@7.215,2.679e-01@9.005,2.700e-01@12.640,2.700e-01@17.631,2.700e-01@20.561,2.700e-01@24.304,2.700e-01@27.288,2.700e-01@29.349')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>0@0.2,0@1,2.5@1.01,2.5@2.0')
        newf.append('<GENERATOR-NOISEOTEMP>270')
        newf.append('<GENERATOR-NOISEOEFF>0.0317@0.2000,0.0437@0.2261,0.0589@0.2580,0.0742@0.2986,0.0851@0.3377,0.0917@0.3667,0.0971@0.4029,0.1015@0.4493,0.1004@0.4971,0.1004@0.5140,0.1670@0.5150,0.1659@0.5377,0.1506@0.6304,0.1255@0.7087,0.0939@0.7986,0.0884@0.8435,0.1146@0.9058,0.1419@0.9594,0.1594@0.9942,0.1821@1.2200,0.1958@1.4100,0.2049@1.6200,0.2094@1.8700,0.2140@2.0000')
        newf.append('<GENERATOR-NOISEOEMIS>0.1')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>10')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-TRANS>03-01')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.0')
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)
        
    if (Telescope == 'LUVOIR A-VIS'):
        newf.append('<GENERATOR-INSTRUMENT>LUVOIR_A-VIS: The Extreme Coronagraph for Living Planetary Systems (ECLIPS) delivers continuous spectral coverage from 200 nm to 2.5 um via three channels, UV (200 to 525 nm), VIS (515 nm to 1030 nm), and NIR (1 to 2 microns). The UV channel is effectively an imager and provides a maximum resolution of RP=7, while the VIS channel RP=140, and NIR=70. The core coronagraph throughput is practically twice for LUVOIR-B than A.')
        newf.append('<GENERATOR-RANGE1>0.515')
        newf.append('<GENERATOR-RANGE2>1.0')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>140')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>CORONA')
        newf.append('<GENERATOR-DIAMTELE>15.0')
        newf.append('<GENERATOR-BEAM>1')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1E-10')
        newf.append('<GENERATOR-TELESCOPE2>4.5')
        newf.append('<GENERATOR-TELESCOPE3>2.700e-11@0.000,2.700e-11@0.217,2.700e-11@1.410,2.110e-03@2.278,6.329e-03@2.767,1.688e-02@3.092,5.063e-02@3.418,1.097e-01@3.797,1.477e-01@4.014,1.751e-01@4.231,1.920e-01@4.557,2.004e-01@5.371,2.004e-01@6.130,2.131e-01@6.618,2.532e-01@7.215,2.679e-01@9.005,2.700e-01@12.640,2.700e-01@17.631,2.700e-01@20.561,2.700e-01@24.304,2.700e-01@27.288,2.700e-01@29.349')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>0@0.2,0@1,2.5@1.01,2.5@2.0')
        newf.append('<GENERATOR-NOISEOTEMP>270')
        newf.append('<GENERATOR-NOISEOEFF>0.0317@0.2000,0.0437@0.2261,0.0589@0.2580,0.0742@0.2986,0.0851@0.3377,0.0917@0.3667,0.0971@0.4029,0.1015@0.4493,0.1004@0.4971,0.1004@0.5140,0.1670@0.5150,0.1659@0.5377,0.1506@0.6304,0.1255@0.7087,0.0939@0.7986,0.0884@0.8435,0.1146@0.9058,0.1419@0.9594,0.1594@0.9942,0.1821@1.2200,0.1958@1.4100,0.2049@1.6200,0.2094@1.8700,0.2140@2.0000')
        newf.append('<GENERATOR-NOISEOEMIS>0.1')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>10')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-TRANS>03-01')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.0')
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)
        
    if (Telescope == 'LUVOIR B-VIS'):
        newf.append('<GENERATOR-INSTRUMENT>LUVOIR_B-VIS: The Extreme Coronagraph for Living Planetary Systems (ECLIPS) delivers continuous spectral coverage from 200 nm to 2.5 um via three channels, UV (200 to 525 nm), VIS (515 nm to 1030 nm), and NIR (1 to 2 microns). The UV channel is effectively an imager and provides a maximum resolution of RP=7, while the VIS channel RP=140, and NIR=70. The core coronagraph throughput is practically twice for LUVOIR-B than A.')
        newf.append('<GENERATOR-RANGE1>0.515')
        newf.append('<GENERATOR-RANGE2>1.0')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>140')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>CORONA')
        newf.append('<GENERATOR-DIAMTELE>6.0')
        newf.append('<GENERATOR-BEAM>1')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1E-10')
        newf.append('<GENERATOR-TELESCOPE2>4.5')
        newf.append('<GENERATOR-TELESCOPE3>4.578000e-11@0.000,4.578000e-11@0.216,4.578000e-11@0.649,2.110e-03@0.973,1.266e-02@1.459,8.228e-02@2.108,1.709e-01@2.973,2.658e-01@4.486,3.418e-01@6.757,3.945e-01@10.216,4.219e-01@14.108,4.409e-01@19.459,4.536e-01@23.622,4.578e-01@27.784,4.578e-01@29.459')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>0@0.2,0@1,2.5@1.01,2.5@2.0')
        newf.append('<GENERATOR-NOISEOTEMP>270')
        newf.append('<GENERATOR-NOISEOEFF>0.0317@0.2000,0.0437@0.2261,0.0589@0.2580,0.0742@0.2986,0.0851@0.3377,0.0917@0.3667,0.0971@0.4029,0.1015@0.4493,0.1004@0.4971,0.1004@0.5140,0.1670@0.5150,0.1659@0.5377,0.1506@0.6304,0.1255@0.7087,0.0939@0.7986,0.0884@0.8435,0.1146@0.9058,0.1419@0.9594,0.1594@0.9942,0.1821@1.2200,0.1958@1.4100,0.2049@1.6200,0.2094@1.8700,0.2140@2.0000')
        newf.append('<GENERATOR-NOISEOEMIS>0.1')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>10')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-TRANS>03-01')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.0')
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)
        
    if (Telescope == 'LUVOIR A-NIR'):
        newf.append('<GENERATOR-INSTRUMENT>LUVOIR_A-VIS: The Extreme Coronagraph for Living Planetary Systems (ECLIPS) delivers continuous spectral coverage from 200 nm to 2.5 um via three channels, UV (200 to 525 nm), VIS (515 nm to 1030 nm), and NIR (1 to 2 microns). The UV channel is effectively an imager and provides a maximum resolution of RP=7, while the VIS channel RP=140, and NIR=70. The core coronagraph throughput is practically twice for LUVOIR-B than A.')
        newf.append('<GENERATOR-RANGE1>1.01')
        newf.append('<GENERATOR-RANGE2>2.0')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>70')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>CORONA')
        newf.append('<GENERATOR-DIAMTELE>15.0')
        newf.append('<GENERATOR-BEAM>1')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1E-10')
        newf.append('<GENERATOR-TELESCOPE2>4.5')
        newf.append('<GENERATOR-TELESCOPE3>2.700e-11@0.000,2.700e-11@0.217,2.700e-11@1.410,2.110e-03@2.278,6.329e-03@2.767,1.688e-02@3.092,5.063e-02@3.418,1.097e-01@3.797,1.477e-01@4.014,1.751e-01@4.231,1.920e-01@4.557,2.004e-01@5.371,2.004e-01@6.130,2.131e-01@6.618,2.532e-01@7.215,2.679e-01@9.005,2.700e-01@12.640,2.700e-01@17.631,2.700e-01@20.561,2.700e-01@24.304,2.700e-01@27.288,2.700e-01@29.349')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>0@0.2,0@1,2.5@1.01,2.5@2.0')
        newf.append('<GENERATOR-NOISEOTEMP>270')
        newf.append('<GENERATOR-NOISEOEFF>0.0317@0.2000,0.0437@0.2261,0.0589@0.2580,0.0742@0.2986,0.0851@0.3377,0.0917@0.3667,0.0971@0.4029,0.1015@0.4493,0.1004@0.4971,0.1004@0.5140,0.1670@0.5150,0.1659@0.5377,0.1506@0.6304,0.1255@0.7087,0.0939@0.7986,0.0884@0.8435,0.1146@0.9058,0.1419@0.9594,0.1594@0.9942,0.1821@1.2200,0.1958@1.4100,0.2049@1.6200,0.2094@1.8700,0.2140@2.0000')
        newf.append('<GENERATOR-NOISEOEMIS>0.1')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>10')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-TRANS>03-01')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.0')
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)
        
    if (Telescope == 'LUVOIR B-NIR'):
        newf.append('<GENERATOR-INSTRUMENT>LUVOIR_B-NIR: The Extreme Coronagraph for Living Planetary Systems (ECLIPS) delivers continuous spectral coverage from 200 nm to 2.5 um via three channels, UV (200 to 525 nm), VIS (515 nm to 1030 nm), and NIR (1 to 2 microns). The UV channel is effectively an imager and provides a maximum resolution of RP=7, while the VIS channel RP=140, and NIR=70. The core coronagraph throughput is practically twice for LUVOIR-B than A.')
        newf.append('<GENERATOR-RANGE1>1.01')
        newf.append('<GENERATOR-RANGE2>2.0')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>70')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>CORONA')
        newf.append('<GENERATOR-DIAMTELE>6.0')
        newf.append('<GENERATOR-BEAM>1')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1E-10')
        newf.append('<GENERATOR-TELESCOPE2>4.5')
        newf.append('<GENERATOR-TELESCOPE3>4.578000e-11@0.000,4.578000e-11@0.216,4.578000e-11@0.649,2.110e-03@0.973,1.266e-02@1.459,8.228e-02@2.108,1.709e-01@2.973,2.658e-01@4.486,3.418e-01@6.757,3.945e-01@10.216,4.219e-01@14.108,4.409e-01@19.459,4.536e-01@23.622,4.578e-01@27.784,4.578e-01@29.459')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>0@0.2,0@1,2.5@1.01,2.5@2.0')
        newf.append('<GENERATOR-NOISEOTEMP>270')
        newf.append('<GENERATOR-NOISEOEFF>0.0317@0.2000,0.0437@0.2261,0.0589@0.2580,0.0742@0.2986,0.0851@0.3377,0.0917@0.3667,0.0971@0.4029,0.1015@0.4493,0.1004@0.4971,0.1004@0.5140,0.1670@0.5150,0.1659@0.5377,0.1506@0.6304,0.1255@0.7087,0.0939@0.7986,0.0884@0.8435,0.1146@0.9058,0.1419@0.9594,0.1594@0.9942,0.1821@1.2200,0.1958@1.4100,0.2049@1.6200,0.2094@1.8700,0.2140@2.0000')
        newf.append('<GENERATOR-NOISEOEMIS>0.1')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>10')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-TRANS>03-01')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.0')
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)
    
    if (Telescope == 'HabEx SS-UV'):
        newf.append('<GENERATOR-INSTRUMENT>HabEx_SS-UV: The HabEx StarShade (SS) will provide extraordinary high-contrast capabilities from the UV (0.2 to 0.45 um), to the visible (0.45 to 1um), and to the infrared (0.975 to 1.8 um). By limiting the number of optical surfaces, this configuration provides high optical throughput (0.2 to 0.4) across this broad of wavelengths, while the quantum efficiency (QE) is expected to be 0.9 for the VU and visible detectors and 0.6 for the infrared detector. The UV channel provides a resolution (RP) of 7, visible channel a maximum of 140 and the infrared 40.')
        newf.append('<GENERATOR-RANGE1>0.2')
        newf.append('<GENERATOR-RANGE2>0.45')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>7')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>CORONA')
        newf.append('<GENERATOR-DIAMTELE>4.0')
        newf.append('<GENERATOR-BEAM>1.0')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1E-10')
        newf.append('<GENERATOR-TELESCOPE2>4.5')
        newf.append('<GENERATOR-TELESCOPE3>7e-11@-0.000e+00,7e-11@-7.483e-03,7e-11@-1.436e-02,3.544e-03@-2.040e-02,1.949e-02@-2.547e-02,3.367e-02@-2.788e-02,6.734e-02@-2.993e-02,1.241e-01@-3.210e-02,2.091e-01@-3.416e-02,2.818e-01@-3.572e-02,3.332e-01@-3.657e-02,3.987e-01@-3.802e-02,4.661e-01@-3.947e-02,5.352e-01@-4.031e-02,6.008e-01@-4.164e-02,6.344e-01@-4.236e-02,6.699e-01@-4.333e-02,6.911e-01@-4.441e-02,7.000e-01@-4.791e-02,7.000e-01@-5.853e-02,7.000e-01@-7.314e-02,7.000e-01@-8.340e-02,7.000e-01@-8.810e-02')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>0.008@0.2,0.008@0.975,0.32@0.976,0.32@1.8')
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@0.975,0.005@0.976,0.005@1.8')
        newf.append('<GENERATOR-NOISEOTEMP>270')
        newf.append('<GENERATOR-NOISEOEFF>0.2260@0.2000,0.2201@0.2112,0.2182@0.2224,0.2300@0.2587,0.2536@0.3035,0.2673@0.3399,0.2791@0.3874,0.2830@0.4294,0.2850@0.4490,0.1796@0.4500,0.1848@0.5161,0.1749@0.6503,0.1474@0.7734,0.1415@0.8434,0.1592@0.9021,0.1926@0.9608,0.1988@0.9750,0.1988@0.9760,0.2162@1.0587,0.2339@1.2462,0.2437@1.4503,0.2516@1.6657,0.2594@1.8000')
        newf.append('<GENERATOR-NOISEOEMIS>0.10')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>10')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>N')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-TRANS>02-01')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)  
        
    if (Telescope == 'HabEx SS-VIS'):
        newf.append('<GENERATOR-INSTRUMENT>HabEx_SS-VIS: The HabEx StarShade (SS) will provide extraordinary high-contrast capabilities from the UV (0.2 to 0.45 um), to the visible (0.45 to 1um), and to the infrared (0.975 to 1.8 um). By limiting the number of optical surfaces, this configuration provides high optical throughput (0.2 to 0.4) across this broad of wavelengths, while the quantum efficiency (QE) is expected to be 0.9 for the VU and visible detectors and 0.6 for the infrared detector. The UV channel provides a resolution (RP) of 7, visible channel a maximum of 140 and the infrared 40.')
        newf.append('<GENERATOR-RANGE1>0.45')
        newf.append('<GENERATOR-RANGE2>0.975')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>140')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>CORONA')
        newf.append('<GENERATOR-DIAMTELE>4.0')
        newf.append('<GENERATOR-BEAM>1.0')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1E-10')
        newf.append('<GENERATOR-TELESCOPE2>4.5')
        newf.append('<GENERATOR-TELESCOPE3>7e-11@-0.000e+00,7e-11@-1.113e-02,7e-11@-2.136e-02,3.544e-03@-3.033e-02,1.949e-02@-3.787e-02,3.367e-02@-4.146e-02,6.734e-02@-4.451e-02,1.241e-01@-4.774e-02,2.091e-01@-5.079e-02,2.818e-01@-5.313e-02,3.332e-01@-5.438e-02,3.987e-01@-5.654e-02,4.661e-01@-5.869e-02,5.352e-01@-5.995e-02,6.008e-01@-6.192e-02,6.344e-01@-6.300e-02,6.699e-01@-6.444e-02,6.911e-01@-6.605e-02,7.000e-01@-7.126e-02,7.000e-01@-8.705e-02,7.000e-01@-1.088e-01,7.000e-01@-1.240e-01,7.000e-01@-1.310e-01')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>0.008@0.2,0.008@0.975,0.32@0.976,0.32@1.8')
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@0.975,0.005@0.976,0.005@1.8')
        newf.append('<GENERATOR-NOISEOTEMP>270')
        newf.append('<GENERATOR-NOISEOEFF>0.2260@0.2000,0.2201@0.2112,0.2182@0.2224,0.2300@0.2587,0.2536@0.3035,0.2673@0.3399,0.2791@0.3874,0.2830@0.4294,0.2850@0.4490,0.1796@0.4500,0.1848@0.5161,0.1749@0.6503,0.1474@0.7734,0.1415@0.8434,0.1592@0.9021,0.1926@0.9608,0.1988@0.9750,0.1988@0.9760,0.2162@1.0587,0.2339@1.2462,0.2437@1.4503,0.2516@1.6657,0.2594@1.8000')
        newf.append('<GENERATOR-NOISEOEMIS>0.10')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>10')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>N')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-TRANS>02-01')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)  
        
    if (Telescope == 'HabEx SS-NIR'):
        newf.append('<GENERATOR-INSTRUMENT>HabEx_SS-NIR: The HabEx StarShade (SS) will provide extraordinary high-contrast capabilities from the UV (0.2 to 0.45 um), to the visible (0.45 to 1um), and to the infrared (0.975 to 1.8 um). By limiting the number of optical surfaces, this configuration provides high optical throughput (0.2 to 0.4) across this broad of wavelengths, while the quantum efficiency (QE) is expected to be 0.9 for the VU and visible detectors and 0.6 for the infrared detector. The UV channel provides a resolution (RP) of 7, visible channel a maximum of 140 and the infrared 40.')
        newf.append('<GENERATOR-RANGE1>0.975')
        newf.append('<GENERATOR-RANGE2>1.80')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>40')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>CORONA')
        newf.append('<GENERATOR-DIAMTELE>4.0')
        newf.append('<GENERATOR-BEAM>1.0')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1E-10')
        newf.append('<GENERATOR-TELESCOPE2>4.5')
        newf.append('<GENERATOR-TELESCOPE3>7e-11@-0.000e+00,7e-11@-1.995e-02,7e-11@-3.830e-02,3.544e-03@-5.439e-02,1.949e-02@-6.791e-02,3.367e-02@-7.434e-02,6.734e-02@-7.982e-02,1.241e-01@-8.561e-02,2.091e-01@-9.108e-02,2.818e-01@-9.526e-02,3.332e-01@-9.752e-02,3.987e-01@-1.014e-01,4.661e-01@-1.052e-01,5.352e-01@-1.075e-01,6.008e-01@-1.110e-01,6.344e-01@-1.130e-01,6.699e-01@-1.155e-01,6.911e-01@-1.184e-01,7.000e-01@-1.278e-01,7.000e-01@-1.561e-01,7.000e-01@-1.950e-01,7.000e-01@-2.224e-01,7.000e-01@-2.349e-01')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>0.008@0.2,0.008@0.975,0.32@0.976,0.32@1.8')
        newf.append('<GENERATOR-NOISE2>3e-5@0.2,3e-5@0.975,0.005@0.976,0.005@1.8')
        newf.append('<GENERATOR-NOISEOTEMP>270')
        newf.append('<GENERATOR-NOISEOEFF>0.2260@0.2000,0.2201@0.2112,0.2182@0.2224,0.2300@0.2587,0.2536@0.3035,0.2673@0.3399,0.2791@0.3874,0.2830@0.4294,0.2850@0.4490,0.1796@0.4500,0.1848@0.5161,0.1749@0.6503,0.1474@0.7734,0.1415@0.8434,0.1592@0.9021,0.1926@0.9608,0.1988@0.9750,0.1988@0.9760,0.2162@1.0587,0.2339@1.2462,0.2437@1.4503,0.2516@1.6657,0.2594@1.8000')
        newf.append('<GENERATOR-NOISEOEMIS>0.10')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>10')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>N')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-TRANS>02-01')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning) 
            
    if (Telescope == 'LIFE'):
        newf.append('<GENERATOR-INSTRUMENT>user')
        newf.append('<GENERATOR-RANGE1>3')
        newf.append('<GENERATOR-RANGE2>20')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>100')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>ARRAY')
        newf.append('<GENERATOR-DIAMTELE>3.5')
        newf.append('<GENERATOR-BEAM>0.006')
        newf.append('<GENERATOR-BEAM-UNIT>arcsec')
        newf.append('<GENERATOR-TELESCOPE1>4')
        newf.append('<GENERATOR-TELESCOPE2>2.0')
        newf.append('<GENERATOR-TELESCOPE3>1.0')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>20.0')
        newf.append('<GENERATOR-NOISE2>0.005')
        newf.append('<GENERATOR-NOISEOTEMP>40')
        newf.append('<GENERATOR-NOISEOEFF>0.4')
        newf.append('<GENERATOR-NOISEOEMIS>0.10')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>8')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-LOGRAD>Y')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-TRANS>02-01')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)
        
    if (Telescope == 'Keck HIRES'):
        newf.append('<GENERATOR-INSTRUMENT>Keck_HIRES: HIRES is the High Resolution Echelle Spectrometer at the Keck observatory on Manuakea (4200m), Hawaii. It is a grating cross-dispersed, echelle spectrograph capable of operating between 0.3 and 1.0 microns at high resolutions (up to RP 85,000). Several slits are available in width and length, since order separation varies between 6 and 43 arcseconds. Typical throughput is 0.15, yet it drops substantially (below 0.05) beyond 0.8 um.')
        newf.append('<GENERATOR-RANGE1>'+HIRES_wav_start)
        newf.append('<GENERATOR-RANGE2>'+HIRES_wav_end)
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>'+HIRES_R)
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>SINGLE')
        newf.append('<GENERATOR-DIAMTELE>10')
        newf.append('<GENERATOR-BEAM>0.5')
        newf.append('<GENERATOR-BEAM-UNIT>arcsec')
        newf.append('<GENERATOR-TELESCOPE1>1')
        newf.append('<GENERATOR-TELESCOPE2>2.0')
        newf.append('<GENERATOR-TELESCOPE3>1.0')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>3.0')
        newf.append('<GENERATOR-NOISE2>0.8')
        newf.append('<GENERATOR-NOISEOTEMP>273')
        newf.append('<GENERATOR-NOISEOEFF>0.15')
        newf.append('<GENERATOR-NOISEOEMIS>0.05')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>8')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>Y')
        newf.append('<GENERATOR-LOGRAD>Y')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-TRANS>02-01')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)
        
    if (Telescope == 'Keck NIRSpec'):
        newf.append('<GENERATOR-INSTRUMENT>Keck_NIRSpec: NIRSpec is the 2nd generation of the Near Infrared Spectrometer at the Keck observatory on Manuakea (4200m), Hawaii. It is a cryogenic cross-dispersed echelle spectrograph which features spectroscopy over the 0.95-5.5 micron range at resolutions up to 40,000 with a 2048x2048 H2RG Teledyne detector. The instrument samples this broad spectral range by employing an array of cross-dispersers and filter configurations.')
        newf.append('<GENERATOR-RANGE1>0.96')
        newf.append('<GENERATOR-RANGE2>5')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>25000')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>SINGLE')
        newf.append('<GENERATOR-DIAMTELE>10')
        newf.append('<GENERATOR-BEAM>0.5')
        newf.append('<GENERATOR-BEAM-UNIT>arcsec')
        newf.append('<GENERATOR-TELESCOPE1>1')
        newf.append('<GENERATOR-TELESCOPE2>2.0')
        newf.append('<GENERATOR-TELESCOPE3>1.0')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>10')
        newf.append('<GENERATOR-NOISE2>0.67')
        newf.append('<GENERATOR-NOISEOTEMP>273')
        newf.append('<GENERATOR-NOISEOEFF>0.18')
        newf.append('<GENERATOR-NOISEOEMIS>0.05')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>8')
        newf.append('<GENERATOR-TRANS-APPLY>Y')
        newf.append('<GENERATOR-TRANS-SHOW>Y')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-TRANS>02-01')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)  

    if (Telescope == 'E-ELT HIRES'):
        newf.append('<GENERATOR-INSTRUMENT>user')
        newf.append('<GENERATOR-RANGE1>'+HIRES_wav_start)
        newf.append('<GENERATOR-RANGE2>'+HIRES_wav_end)
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>'+HIRES_R)
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>SINGLE')
        newf.append('<GENERATOR-DIAMTELE>39.1')
        newf.append('<GENERATOR-BEAM>0.5')
        newf.append('<GENERATOR-BEAM-UNIT>arcsec')
        newf.append('<GENERATOR-TELESCOPE1>1')
        newf.append('<GENERATOR-TELESCOPE2>2.0')
        newf.append('<GENERATOR-TELESCOPE3>1.0')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>3.0')
        newf.append('<GENERATOR-NOISE2>0.8')
        newf.append('<GENERATOR-NOISEOTEMP>273')
        newf.append('<GENERATOR-NOISEOEFF>0.15')
        newf.append('<GENERATOR-NOISEOEMIS>0.05')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>8')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>Y')
        newf.append('<GENERATOR-LOGRAD>Y')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-TRANS>02-01')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)

    if (Telescope == 'E-ELT METIS'):
        newf.append('<GENERATOR-INSTRUMENT>user')
        newf.append('<GENERATOR-RANGE1>3')
        newf.append('<GENERATOR-RANGE2>13')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>50000')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>SINGLE')
        newf.append('<GENERATOR-DIAMTELE>39.1')
        newf.append('<GENERATOR-BEAM>0.5')
        newf.append('<GENERATOR-BEAM-UNIT>arcsec')
        newf.append('<GENERATOR-TELESCOPE1>1')
        newf.append('<GENERATOR-TELESCOPE2>2.0')
        newf.append('<GENERATOR-TELESCOPE3>1.0')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>10')
        newf.append('<GENERATOR-NOISE2>0.67')
        newf.append('<GENERATOR-NOISEOTEMP>273')
        newf.append('<GENERATOR-NOISEOEFF>0.18')
        newf.append('<GENERATOR-NOISEOEMIS>0.05')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>8')
        newf.append('<GENERATOR-TRANS-APPLY>Y')
        newf.append('<GENERATOR-TRANS-SHOW>Y')
        newf.append('<GENERATOR-LOGRAD>Y')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-TRANS>02-01')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)  
    
    if (Telescope == 'JWST NIRSpec'):
        newf.append('<GENERATOR-INSTRUMENT>JWST_NIRSpec-2700: JWST/NIRSpec is the infrared spectrometer onboard the James Webb Space Telescope. The high-resolution configuration (RP 2700) considers the maximum resolving power for each of the available instrument gratings/filters combinations (G140H+F070LP, G140H+F100LP, G235H+F170LP, G395H+F290LP). The three gratings cover the 1.0 to 5.3 um spectral region with an average resolving power of 2700.')
        newf.append('<GENERATOR-RANGE1>1')
        newf.append('<GENERATOR-RANGE2>5.3')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>2700')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>SINGLE')
        newf.append('<GENERATOR-DIAMTELE>5.64')
        newf.append('<GENERATOR-BEAM>1')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1')
        newf.append('<GENERATOR-TELESCOPE2>2.0')
        newf.append('<GENERATOR-TELESCOPE3>1.0')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>16.8')
        newf.append('<GENERATOR-NOISE2>0.005')
        newf.append('<GENERATOR-NOISEOTEMP>50')
        newf.append('<GENERATOR-NOISEOEFF>0.3')
        newf.append('<GENERATOR-NOISEOEMIS>0.10')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>8')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-LOGRAD>Y')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-TRANS>02-01')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)  
        
    if (Telescope == 'JWST MIRI'):
        newf.append('<GENERATOR-INSTRUMENT>JWST_MIRI-LRS: JWST/MIRI is the Mid-Infrared Instrument onboard the James Webb Space Telescope. The low resolution spectroscopy (LRS) configuration samples the 5 to 12 um spectral region with a resolving power of 100 (at 7.5um), and observations can be performed with a 0.5 arcsec slit or via slitless spectroscopy. Fringing and other systematic sources of noise have been identified in this instrument, yet these effects are not included in the noise simulator. Read noise is for FAST readout (32 e-/pixel/read) and for slitless measurements.')
        newf.append('<GENERATOR-RANGE1>5.00')
        newf.append('<GENERATOR-RANGE2>12.00')
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>100')
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>SINGLE')
        newf.append('<GENERATOR-DIAMTELE>5.64')
        newf.append('<GENERATOR-BEAM>1')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1')
        newf.append('<GENERATOR-TELESCOPE2>2.0')
        newf.append('<GENERATOR-TELESCOPE3>1.0')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>32')
        newf.append('<GENERATOR-NOISE2>0.2')
        newf.append('<GENERATOR-NOISEOTEMP>50')
        newf.append('<GENERATOR-NOISEOEFF>0.00@4.52,0.04@4.72,0.07@4.90,0.12@4.98,0.16@5.13,0.18@5.18,0.20@5.24,0.23@5.39,0.25@5.53,0.27@5.88,0.29@6.14,0.30@6.37,0.30@6.60,0.29@6.89,0.29@7.23,0.30@7.55,0.31@7.75,0.33@8.07,0.32@8.48,0.32@8.85,0.31@9.17,0.29@9.57,0.27@9.95,0.26@10.30,0.24@10.56,0.22@10.79,0.20@10.96,0.17@11.25,0.15@11.48,0.13@11.74,0.11@12.09,0.09@12.40,0.08@12.52')
        newf.append('<GENERATOR-NOISEOEMIS>0.10')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>0.76@5.01,1.03@5.65,1.44@6.57,1.87@7.54,2.17@8.24,2.49@8.95,2.76@9.55,3.08@10.30,3.44@11.03,3.72@11.69,4.03@12.35,4.30@13.02')
        newf.append('<GENERATOR-NOISEWELL>250000')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-TRANS>02-01')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)  
        
    if (Telescope == 'Ideal'):
        newf.append('<GENERATOR-INSTRUMENT>user')
        newf.append('<GENERATOR-RANGE1>'+Ideal_wav_start)
        newf.append('<GENERATOR-RANGE2>'+Ideal_wav_end)
        newf.append('<GENERATOR-RANGEUNIT>um')
        newf.append('<GENERATOR-RESOLUTION>'+Ideal_R)
        newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
        newf.append('<GENERATOR-TELESCOPE>SINGLE')
        newf.append('<GENERATOR-DIAMTELE>10')
        newf.append('<GENERATOR-BEAM>1')
        newf.append('<GENERATOR-BEAM-UNIT>diffrac')
        newf.append('<GENERATOR-TELESCOPE1>1')
        newf.append('<GENERATOR-TELESCOPE2>2.0')
        newf.append('<GENERATOR-TELESCOPE3>1.0')
        newf.append('<GENERATOR-NOISE>CCD')
        newf.append('<GENERATOR-NOISE1>20')
        newf.append('<GENERATOR-NOISE2>0.005')
        newf.append('<GENERATOR-NOISEOTEMP>288')
        newf.append('<GENERATOR-NOISEOEFF>0.4')
        newf.append('<GENERATOR-NOISEOEMIS>0.10')
        newf.append('<GENERATOR-NOISETIME>'+str(exposure_time))
        newf.append('<GENERATOR-NOISEFRAMES>'+str(exposure_number))
        newf.append('<GENERATOR-NOISEPIXELS>8')
        newf.append('<GENERATOR-NOISEWELL>250000')
        newf.append('<GENERATOR-TRANS-APPLY>N')
        newf.append('<GENERATOR-TRANS-SHOW>N')
        newf.append('<GENERATOR-LOGRAD>N')
        newf.append('<GENERATOR-GAS-MODEL>Y')
        newf.append('<GENERATOR-CONT-MODEL>Y')
        newf.append('<GENERATOR-CONT-STELLAR>Y')
        newf.append('<GENERATOR-RADUNITS>'+obs_type)
        newf.append('<GENERATOR-RESOLUTIONKERNEL>N')
        newf.append('<GENERATOR-TRANS>02-01')
        newf.append('<GENERATOR-GCM-BINNING>'+Binning)  
        
    # Save GCM parameters
    vars = '<ATMOSPHERE-GCM-PARAMETERS>'
    vars = vars + str("%d,%d,%d,%.1f,%.1f,%.2f,%.2f" %(sz[0],sz[1],sz[2],lon[0],lat[0],lon[1]-lon[0],lat[1]-lat[0]))
    vars = vars + ',Winds,Tsurf,Psurf,Albedo,Temperature,Pressure,N2,H2O,O2,O3,CO2,N2O,CH4,Water,WaterIce'
    newf.append(vars)
    
    with open(fileout,'w') as fw:
    	for i in newf: fw.write(i+'\n')
    with open(fileout,'ab') as fb:
        if sys.version_info>=(3,0,0): 
            bc=fb.write(bytes('<BINARY>',encoding = 'utf-8'))
        else: 
            bc=fb.write('<BINARY>')
        fb.write(np.asarray(U,order='C'))
        fb.write(np.asarray(V,order='C'))
        fb.write(np.asarray(TS,order='C'))
        fb.write(np.log10(np.asarray(PS,order='C')))
        fb.write(np.asarray(ASDIR,order='C'))   
        fb.write(np.asarray(T,order='C'))
        fb.write(np.log10(np.asarray(np.transpose(press3D),order='C')))
        fb.write(np.log10(np.asarray(N2,order='C')))
        fb.write(np.log10(np.asarray(H2O,order='C')))
        fb.write(np.log10(np.asarray(O2,order='C')))
        fb.write(np.log10(np.asarray(O3,order='C')))
        fb.write(np.log10(np.asarray(CO2,order='C')))
        fb.write(np.log10(np.asarray(N2O,order='C')))
        fb.write(np.log10(np.asarray(CH4,order='C'))) 
        fb.write(np.log10(np.asarray(CLDLIQ,order='C')))
        fb.write(np.log10(np.asarray(CLDICE,order='C')))
        if sys.version_info>=(3,0,0): bc=fb.write(bytes('</BINARY>',encoding = 'utf-8'))
        else: bc=fb.write('</BINARY>')
        
    fb.close()
    # print(U.shape, ASDIR.shape, press3D.shape)
        
#%% Now upload to PSG by calling function

infile = '/Users/gregcooke/b.e21.BWma1850.f19_g17.K_star_HZ.10pc_o2.0_obliq.001.cam.h1.0030-01-01-00000.nc'
fileout = '/Users/gregcooke/out.gcm'

if (transit == True):
    phase = 180 
else:
    phase = season

distance = 10

Name = 'test'
velocity = 0

Convert_GCM(infile = infile,
        fileout = fileout, phase = str(phase), obs_type = obs_type,
        sollon = sollon, sollat = sollat, date = date, 
        distance = str(distance), velocity = str(velocity),
        season = season, name = Name, Telescope = Telescope_choice)


#%% Calculations are finished
print('\nCalculations and/or uploads finished\n')
