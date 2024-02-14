#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 11:51:06 2022

@author: Greg Cooke
"""

#%% imports

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import cm
import cartopy.crs as ccrs
import matplotlib.colors as colors

#%% return subscript number or text
def sub(num):
    return r'$_{'+str(num)+'}$'

#%% return superscript number or text
def sup(num):
    return r'$^{'+str(num)+'}$'
#%% define colour bar function

Fontsize = 14
def cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]', orientation='vertical', tick = False, ticks = [2,3], shrink = 1):
    cbar = plt.colorbar(model, orientation = orientation, shrink = shrink)
    if (label == True):
        clabel = clabel
        cbar.set_label(clabel,size=15)
    cbar.ax.tick_params(labelsize=15)
    if (tick == True):
        cbar.set_ticks(ticks)

#%% import file for pressure and gaussian weights

File =  "/Users/gregcooke/Final_thesis/Data/b.e21.BWma1850.f19_g17.TRAPPIST1_e_MUSCLES.023.cam.h0.0270-0279.nc" #file name

Gw_file = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

gw = Gw_file.gw.values
Pressure = Gw_file.lev.values
#%%
def LWAV(DS, time = False):
    try:
        DS = DS.mean(dim='lon')
    except:
        DS = DS
    if (time == False):
        try:
            DS = DS.mean(dim='time')
        except:
            DS = DS
    DS = gw*DS #multiply by gaussian weights
    try:
        DS = np.sum(DS,axis=1)/gw.sum() #weighted mean
    except:
        DS = np.sum(DS,axis=0)/gw.sum() #weighted mean 
    return DS #return a numpy array that has been modified above

#%% define ozone column calculation function

def O3_col(ds, time = False, O2_mr = 0.21, lon = True, g = 9.145):
    
    g = g # gravity acceleration (cm/sec2)
    to_DU = 1./2.6867e20  # convert colunm density / cm^2 to Dobson units
    N2_mr = 1-O2_mr #nitrogen mixing ratio
    m_avg = ((28*N2_mr)+(O2_mr*32))*1.661e-27 # average weight of atmosphere
    
    #calcualte column integrals 
    o3 = ds['O3']
    
    # variables to calcualate hybrid pressure on model interfaces
    ps = ds['PS'].values # surface pressure
    p0 = ds['P0'].values
    hyai = ds['hyai'].values
    hybi = ds['hybi'].values
    
    #col = np.ndarray(ps.shape, np.float32)
    dp = np.ndarray(o3.shape, np.float32)
    
    nt = 1
    try:
        nt = ds['time'].size
    except:
        nt = 0
    ny = ds['lat'].size
    nz = ds['lev'].size
    nx = ds['lon'].size
    
    press = np.ndarray(nz, np.float32)
    if (nt > 0):
        for i in range(nt):
            for j in range(ny):
                for k in range(nx):
                    press = hyai * p0 + hybi * ps[i,j,k]
                    for iz in range(nz):
                        dp[i,iz,j,k] = (press[iz+1]-press[iz])
    else:
        for j in range(ny):
            for k in range(nx):
                press = hyai * p0 + hybi * ps[j,k]
                for iz in range(nz):
                    dp[iz,j,k] = (press[iz+1]-press[iz])
                    
    #calcualte column integrals 
    o3 = o3 * dp / m_avg
    
    if (nt > 0):
        for i in range(nt):
            for j in range(96):
                for k in range(144):
                    o3[i,:,j,k] = o3[i,:,j,k]/g    
    else:
        for j in range(96):
            for k in range(144):
                o3[:,j,k] = o3[:,j,k]/g
    
    if (nt > 0):
        o3col = np.sum(o3, axis = 1) * to_DU # convert to DU
    else:
        o3col = np.sum(o3, axis = 0) * to_DU # convert to DU
    if (lon == False):
        o3col  = o3col.mean(dim='lon')
    if (time == False):
        try:
            o3col  = o3col.mean(dim='time')
        except:
            o3col  = o3col
    print('O3 columns calculated')
    return o3col

#%% Define general column functions

def general_col(ds, Var = 'H2O', time = False, O2_mr = 0.21, lon = True, g = 9.145):
    
    g = g # gravity acceleration (cm/sec2)
    #to_DU = 1./2.6867e20  # convert colunm density / cm^2 to Dobson units
    N2_mr = 1-O2_mr #nitrogen mixing ratio
    m_avg = ((28*N2_mr)+(O2_mr*32))*1.661e-27 # average weight of atmosphere

    #calcualte column integrals 
    var = ds[Var]

    # variables to calcualate hybrid pressure on model interfaces
    ps = ds['PS'].values # surface pressure
    p0 = ds['P0'].values
    hyai = ds['hyai'].values
    hybi = ds['hybi'].values

    #col = np.ndarray(ps.shape, np.float32)
    dp = np.ndarray(var.shape, np.float32)

    nt = 1
    try:
        nt = ds['time'].size
    except:
        nt = 0
    ny = ds['lat'].size
    nz = ds['lev'].size
    nx = ds['lon'].size

    press = np.ndarray(nz, np.float32)
    if (nt > 0):
        for i in range(nt):
            for j in range(ny):
                for k in range(nx):
                    press = hyai * p0 + hybi * ps[i,j,k]
                    for iz in range(nz):
                        dp[i,iz,j,k] = (press[iz+1]-press[iz])
    else:
        for j in range(ny):
            for k in range(nx):
                press = hyai * p0 + hybi * ps[j,k]
                for iz in range(nz):
                    dp[iz,j,k] = (press[iz+1]-press[iz])
                    
    #calcualte column integrals 
    var = var * dp / (m_avg*g)

    varcol = np.sum(var, axis = 0) # convert to DU
    lon = True
    if (lon == False):
        varcol  = varcol.mean(dim='lon')
    if (time == False):
        try:
            varcol  = varcol.mean(dim='time')
        except:
            varcol  = varcol
    print('Columns calculated')
    return varcol

#%%General robinson plot for columns

def Surface_rob_plot_gen_col(cmap = 'coolwarm', mmw = 18, var = 'H2O',
                             clabel = 'H'+ sub(2) +'O column [kg m'+sup(-2)+']',
                             levels = 12, save_name = 'Surface_temp', 
                             log = False, vmax = 10, vmin = 0.1):

    Figure = plt.figure(figsize = (12,7))
    #clabel = 'Zonal winds [ms'+sup(-1)+']'
    clabel = clabel
    cmap = cmap
    gs = gridspec.GridSpec(2,2)
    gs.update(hspace = 0.12)
    gs.update(wspace = 0.05)
    
    ax = plt.subplot(gs[0,0], projection=ccrs.Robinson(central_longitude=180))
    ax.coastlines('110m')
    
    P19_col = (general_col(P19.mean(dim = 'time'), Var = var, time = False, O2_mr = 0.21, lon = True, g = 9.145))
    P19_col = P19_col*1.66e-27*mmw
    
    if(log == True):
        model = P19_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, norm = colors.LogNorm(vmax = vmax, vmin = vmin), transform=ccrs.PlateCarree())
    else:
        model = P19_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, transform=ccrs.PlateCarree(), vmax = vmax, vmin =  vmin)
    plt.tick_params(labelsize = Fontsize, length = 4, width = 2, right = True)
    plt.title('P19 PI', fontsize = 15)
    #cbar(model, label = True, clabel = clabel, orientation='vertical')
    
    W21_col = (general_col(W21.mean(dim = 'time'), Var = var, time = False, O2_mr = 0.21, lon = True, g = 9.145))
    W21_col = W21_col*1.66e-27*mmw
    
    ax = plt.subplot(gs[0,1], projection=ccrs.Robinson(central_longitude=180))
    ax.coastlines('110m')
    if(log == True):
        model = W21_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, norm = colors.LogNorm(vmax = vmax, vmin = vmin), transform=ccrs.PlateCarree())
    else:
        model = W21_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, transform=ccrs.PlateCarree(), vmax = vmax, vmin =  vmin)
    plt.tick_params(labelsize = Fontsize, length = 4, width = 2, right = True, labelleft = False)
    plt.title('W21 PI', fontsize = 15)
    
    P19_01pc_col = (general_col(P19_01pc.mean(dim = 'time'), Var = var, time = False, O2_mr = 0.00021, lon = True, g = 9.145))
    P19_01pc_col = P19_01pc_col*1.66e-27*mmw
    
    ax = plt.subplot(gs[1,0], projection=ccrs.Robinson(central_longitude=180))
    ax.coastlines('110m')
    
    if(log == True):
        model = P19_01pc_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, norm = colors.LogNorm(vmax = vmax, vmin = vmin), transform=ccrs.PlateCarree())
    else:
        model = P19_01pc_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, transform=ccrs.PlateCarree(), vmax = vmax, vmin =  vmin)
    plt.tick_params(labelsize = Fontsize, length = 4, width = 2, right = True)
    plt.title('P19 0.1% PAL', fontsize = 15)
    #cbar(model, label = True, clabel = clabel, orientation='vertical')
    
    W21_01pc_col = (general_col(W21_01pc.mean(dim = 'time'), Var = var, time = False, O2_mr = 0.00021, lon = True, g = 9.145))
    W21_01pc_col = W21_01pc_col*1.66e-27*mmw
    
    ax = plt.subplot(gs[1,1], projection=ccrs.Robinson(central_longitude=180))
    ax.coastlines('110m')
    if(log == True):
        model = W21_01pc_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, norm = colors.LogNorm(vmax = vmax, vmin = vmin), transform=ccrs.PlateCarree())
    else:
        model = W21_01pc_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, transform=ccrs.PlateCarree(), vmax = vmax, vmin =  vmin)
    plt.tick_params(labelsize = Fontsize, length = 4, width = 2, right = True, labelleft = False)
    plt.title('W21 0.1% PAL', fontsize = 15)
    #cbar(model, label = True, clabel = clabel, orientation='vertical')

    cax = Figure.add_axes([0.92, 0.135, 0.025, 0.74])#start left, vertical start, width, height
    Cbar = plt.colorbar(model, orientation='vertical', cax=cax)
    Cbar.set_label(clabel, size=15)
    Cbar.ax.tick_params(labelsize=15) 
    
    plt.savefig('/Users/gregcooke/python_output/'+save_name + '.png', dpi = 200, bbox_inches = 'tight')
    plt.savefig('/Users/gregcooke/python_output/'+save_name + '.eps', format = 'eps', dpi = 200, bbox_inches = 'tight')


#%% color choices

blues = cm.get_cmap('Blues')
W21_no_TL_color = 'k'
W21_color = '#2980C9'
W21_10pc_color = 'b'
W21_1pc_color = '#ADB0FC'
W21_01pc_color = 'slategrey'

reds = cm.get_cmap('Reds')
P19_no_TL_color = reds(1.0)
P19_color = '#E37D0F'
P19_10pc_color = '#994A05'
P19_1pc_color = '#E0D713'
P19_01pc_color = 'crimson'

greens = cm.get_cmap('Greens')
darkgreen = greens(0.8)
green = greens(0.6)
lightgreen = greens(0.4)

#%% set tick parameters

plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['text.usetex'] = False

#%%

def define_time(i = 0, starting_date = '0120-01-01-00000', year = '0120', month = '01', day = '01', time_of_day= '0000'):
    if (i == 0 or i == 48 or i == 96 or i == 144 or i == 192 or i == 240 or i == 288):
        time_of_day = '00000'    
        if(i == 48 or i == 96 or i == 144 or i == 192 or i == 240 or i == 288):
            day = str(int(day)+int(i/48))
    else:
        time_of_day = str(i*1800)   
        
    if (int(time_of_day) >= 86400):
        day = str(int(day)+(int(time_of_day)//86400))    
        time_of_day =  str(int(time_of_day)%86400)
        
    while (len(time_of_day) < 5):
        time_of_day = '0' + time_of_day
        
    if (len(day) < 2):
        day = '0'+day
        
    time = year+'-'+month+'-'+day+'-'+time_of_day
    
    return time

#%% define cumulative ozone column function
def cumulative_o3_col(ds, time = False, O2_mr = 0.21, lon = True, g = 9.80616):
    
    
    g = g# gravity acceleration (cm/sec2)
    to_DU = 1./2.6867e20  # convert colunm density / cm^2 to Dobson units
    N2_mr = 1-O2_mr #nitrogen mixing ratio
    m_avg = ((28*N2_mr)+(O2_mr*32))*1.661e-27 # average weight of atmosphere
    
    #calcualte column integrals 
    o3 = ds['O3']
    
    # variables to calcualate hybrid pressure on model interfaces
    ps = ds['PS'].values # surface pressure
    p0 = ds['P0'].values
    hyai = ds['hyai'].values
    hybi = ds['hybi'].values
    

    #col = np.ndarray(ps.shape, np.float32)
    dp = np.ndarray(o3.shape, np.float32)
    
    nt = 1
    try:
        nt = ds['time'].size
    except:
        nt = 0
    ny = ds['lat'].size
    nz = ds['lev'].size
    nx = ds['lon'].size
    
    press = np.ndarray(nz, np.float32)
    if (nt > 0):
        for i in range(nt):
            for j in range(ny):
                for k in range(nx):
                    press = hyai * p0 + hybi * ps[i,j,k]
                    for iz in range(nz):
                        dp[i,iz,j,k] = (press[iz+1]-press[iz])
    else:
        for j in range(ny):
            for k in range(nx):
                press = hyai * p0 + hybi * ps[j,k]
                for iz in range(nz):
                    dp[iz,j,k] = (press[iz+1]-press[iz])
                    
    #calcualte column integrals 
    o3 = o3 * dp / m_avg
    
    if (nt > 0):
        for i in range(nt):
            for j in range(96):
                for k in range(144):
                    o3[i,:,j,k] = o3[i,:,j,k]/g    
    else:
        for j in range(96):
            for k in range(144):
                o3[:,j,k] = (o3[:,j,k]/g) 
    o3col = o3 * to_DU
    #o3col = np.sum(o3, axis = 0)  # convert to DU
    if (lon == False):
        o3col  = o3col.mean(dim='lon')
    if (time == False):
        try:
            o3col  = o3col.mean(dim='time')
        except:
            o3col  = o3col
    print('O3 columns calculated')
    return o3col

#%% define number density function

def n_dens(ds, var = 'H2O'):
    k = 1.381e-23
    ndens = ds[var]*ds.lev*100/(ds.T*k)
    return ndens

#%%

Blues = cm.get_cmap('Blues')
darkblue = Blues(0.99)
blue = Blues(0.7)
lightblue = Blues(0.4)

Reds = cm.get_cmap('Reds')
darkred = Reds(0.99)
red = Reds(0.7)
lightred = Reds(0.4)

Reds = cm.get_cmap('Blues')
darkred = Reds(0.99)
red = Reds(0.7)
lightred = Reds(0.3)

Oranges = cm.get_cmap('Oranges')
darkorange = Oranges(0.99)
orange = Oranges(0.7)
lightorange = Oranges(0.4)

#%% Read in tidally locked cases
path = '/Users/gregcooke/TP1e_ApJ_paper_data/'

File =  path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.W21.cam.h0.0293-0302.nc" #file name
W21 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.W21.cam.h1.0293-0302.nc"  
W21_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

#checked and updated
File =  path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.W21.no_T_lock.cam.h0.0310-0319.nc" #file name
W21_no_TL = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.W21.no_T_lock.cam.h1.0310-0319.nc"  
W21_no_TL_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

#checked and updated
File =  path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.W21.10pc_o2.cam.h0.0300-0309.nc" #file name
W21_10pc = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.W21.10pc_o2.cam.h1.0300-0309.nc"  
W21_10pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

#checked and updated
File =  path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.W21.1pc_o2.cam.h0.0271-0280.nc" #file name
W21_1pc = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.W21.1pc_o2.cam.h1.0271-0280.nc"  
W21_1pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

#W21 0.1% PAL case files
File =  path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.W21.0.1pc_o2.cam.h0.0250-0259.nc" #file name
W21_01pc = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.W21.0.1pc_o2.cam.h1.0250-0259.nc"  
W21_01pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

#P19 case files
File =  path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.P19.cam.h0.0320-0329.nc" #file name
P19 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.P19.cam.h1.0320-0329.nc"  
P19_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

#P19 noTL case files
File =  path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.P19.no_T_lock.cam.h0.0280-0289.nc" #file name
P19_no_TL = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.P19.no_T_lock.cam.h1.0280-0289.nc"  
P19_no_TL_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

#P19 10% PAL case files
File =  path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.P19.10pc_o2.cam.h0.0350-0359.nc" #file name
P19_10pc = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.P19.10pc_o2.cam.h1.0350-0359.nc"  
P19_10pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

#P19 1% PAL case files
File =  path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.P19.1pc_o2.cam.h0.0310-0319.nc" #file name
P19_1pc = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.P19.1pc_o2.cam.h1.0310-0319.nc"  
P19_1pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

#P19 0.1% PAL case files
File =  path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.P19.0.1pc_o2.cam.h0.0300-0309.nc" #file name
P19_01pc = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.TRAPPIST1_e.P19.0.1pc_o2.cam.h1.0300-0309.nc"  
P19_01pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false


#%% O2 cases

path = '/Users/gregcooke/Final_Thesis/Data/'

File = path+"b.e21.BWma1850.f19_g17.baseline.cam.h0.0014-0017.nc"
PI = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# Pre-industrial baseliine run (PI; 100% PAL of O2)
File =  path+"b.e21.BWma1850.f19_g17.baseline.cam.h0.0014-0017.nc" #file name
Pre_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 50% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.50pc_o2.002.cam.h0.0031-0034.nc" #file name
Fifty_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 10% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.10pc_o2.001.cam.h0.0034-0037.nc" #file name
Ten_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 5% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.5pc_o2.002.cam.h0.0032-0035.nc" #file name
Five_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 1% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.1pc_o2.cam.h0.0028-0031.nc" #file name
One_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 0.5% PAL of O2
File = path+"b.e21.BWma1850.f19_g17.0.5pc_o2.001.cam.h0.0046-0049.nc" #file name
Zero5_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 0.1% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.0.1pc_o2.001.cam.h0.0080-0083.nc" #file name
Zero1_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

File =  path+"b.e21.BWma1850.f19_g17.150pc_o2.002.cam.h0.0034-0037.nc" #file name
One50_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

#%% Calculate cumulative ozone columns

W21_o3 = cumulative_o3_col(W21, time = False, O2_mr = 0.21, lon = True, g = 9.1454)
W21_no_TL_o3 = cumulative_o3_col(W21_no_TL, time = False, O2_mr = 0.21, lon = True, g = 9.1454)
W21_10pc_o3 = cumulative_o3_col(W21_10pc, time = False, O2_mr = 0.021, lon = True, g = 9.1454)
W21_1pc_o3 = cumulative_o3_col(W21_1pc, time = False, O2_mr = 0.0021, lon = True, g = 9.1454)
W21_01pc_o3 = cumulative_o3_col(W21_01pc, time = False, O2_mr = 0.00021, lon = True, g = 9.1454)

P19_o3 = cumulative_o3_col(P19, time = False, O2_mr = 0.21, lon = True, g = 9.1454)
P19_no_TL_o3 = cumulative_o3_col(P19_no_TL, time = False, O2_mr = 0.21, lon = True, g = 9.1454)
P19_10pc_o3 = cumulative_o3_col(P19_10pc, time = False, O2_mr = 0.021, lon = True, g = 9.1454)
P19_1pc_o3 = cumulative_o3_col(P19_1pc, time = False, O2_mr = 0.0021, lon = True, g = 9.1454)
P19_01pc_o3 = cumulative_o3_col(P19_01pc, time = False, O2_mr = 0.00021, lon = True, g = 9.1454)

PI_o3 = cumulative_o3_col(PI, time = False, O2_mr = 0.21, lon = True, g = 9.81)

#%% cumulative ozone column plot

# W21 data
W21_cumu_o3 = LWAV(W21_o3.cumsum(dim='lev'))
W21_no_TL_cumu_o3 = LWAV(W21_no_TL_o3.cumsum(dim='lev'))
W21_10pc_cumu_o3 = LWAV(W21_10pc_o3.cumsum(dim='lev'))
W21_1pc_cumu_o3 = LWAV(W21_1pc_o3.cumsum(dim='lev'))
W21_01pc_cumu_o3 = LWAV(W21_01pc_o3.cumsum(dim='lev'))

# P19 data
P19_cumu_o3 = LWAV(P19_o3.cumsum(dim='lev'))
P19_no_TL_cumu_o3 = LWAV(P19_no_TL_o3.cumsum(dim='lev'))
P19_10pc_cumu_o3 = LWAV(P19_10pc_o3.cumsum(dim='lev'))
P19_1pc_cumu_o3 = LWAV(P19_1pc_o3.cumsum(dim='lev'))
P19_01pc_cumu_o3 = LWAV(P19_01pc_o3.cumsum(dim='lev'))

#PI data
PI_cumu_o3 = np.cumsum(LWAV(PI_o3))
      
plt.figure(figsize = (7,5))
plt.plot(W21_cumu_o3, LWAV(W21.Z3)/1e3, color = W21_color, label = 'W21 PI')
plt.plot(W21_no_TL_cumu_o3, LWAV(W21_no_TL.Z3)/1e3, color = W21_no_TL_color, label = 'W21 PI noTL')
plt.plot(W21_10pc_cumu_o3, LWAV(W21_10pc.Z3)/1e3, color = W21_10pc_color, label = 'W21 10% PAL')
plt.plot(W21_1pc_cumu_o3, LWAV(W21_1pc.Z3)/1e3, color = W21_1pc_color, label = 'W21 1% PAL')
plt.plot(W21_01pc_cumu_o3, LWAV(W21_01pc.Z3)/1e3, color = W21_01pc_color, label = 'W21 0.1% PAL')

plt.plot(P19_cumu_o3, LWAV(P19.Z3)/1e3, color = P19_color, label = 'P19 PI')
plt.plot(P19_no_TL_cumu_o3, LWAV(P19_no_TL.Z3)/1e3, color = P19_no_TL_color, label = 'P19 PI noTL')
plt.plot(P19_10pc_cumu_o3, LWAV(P19_10pc.Z3)/1e3, color = P19_10pc_color, label = 'P19 10% PAL')
plt.plot(P19_1pc_cumu_o3, LWAV(P19_1pc.Z3)/1e3, color = P19_1pc_color, label = 'P19 1% PAL')
plt.plot(P19_01pc_cumu_o3, LWAV(P19_01pc.Z3)/1e3, color = P19_01pc_color, label = 'P19 0.1% PAL')

plt.plot(PI_cumu_o3, LWAV(P19_01pc.Z3)/1e3, color = 'k', label = 'PI')

plt.ylim(0,80)
plt.xlim(0.01)
plt.xscale('log')
plt.ylabel('Pressure [hPa]', fontsize = 15)
plt.xlabel('Cumulative ozone column [DU]', fontsize = 15)
plt.tick_params(labelsize = 15, right = True)

plt.legend(loc = 0)

#%% Calculate ozone columns

# W21 ozone columns
W21_O3_col = O3_col(W21.mean(dim = 'time'), time = True, lon = True, O2_mr = 0.21, g = 9.1454)
W21_10pc_O3_col = O3_col(W21_10pc.mean(dim = 'time'), time = True, lon = True, O2_mr = 0.21/10, g = 9.1454)
W21_1pc_O3_col = O3_col(W21_1pc.mean(dim = 'time'), time = True, lon = True, O2_mr = 0.21/100, g = 9.1454)
W21_01pc_O3_col = O3_col(W21_01pc.mean(dim = 'time'), time = True, lon = True, O2_mr = 0.21/1000, g = 9.1454)
W21_no_TL_O3_col = O3_col(W21_no_TL.mean(dim = 'time'), time = True, lon = True, O2_mr = 0.21, g = 9.1454)

# P19 ozone columns
P19_O3_col = O3_col(P19.mean(dim = 'time'), time = True, lon = True, O2_mr = 0.21, g = 9.1454)
P19_10pc_O3_col = O3_col(P19_10pc.mean(dim = 'time'), time = True, lon = True, O2_mr = 0.21/10, g = 9.1454)
P19_1pc_O3_col = O3_col(P19_1pc.mean(dim = 'time'), time = True, lon = True, O2_mr = 0.21/100, g = 9.1454)
P19_01pc_O3_col = O3_col(P19_01pc.mean(dim = 'time'), time = True, lon = True, O2_mr = 0.21/1000, g = 9.1454)
P19_no_TL_O3_col = O3_col(P19_no_TL.mean(dim = 'time'), time = True, lon = True, O2_mr = 0.21, g = 9.1454)

# Earth ozone columns
O3_col_One50 = O3_col(One50_pc_h0.mean(dim='time'), time = True, lon = True, O2_mr = 0.21*1.5, g = 9.81)
O3_col_Pre = O3_col(Pre_pc_h0.mean(dim='time'), time = True, lon = True, O2_mr = 0.21, g = 9.81)
O3_col_Fifty = O3_col(Fifty_pc_h0.mean(dim='time'), time = True, lon = True, O2_mr = 0.21*0.5, g = 9.81)
O3_col_Ten = O3_col(Ten_pc_h0.mean(dim='time'), time = True, lon = True, O2_mr = 0.21*0.1, g = 9.81)
O3_col_Five = O3_col(Five_pc_h0.mean(dim='time'), time = True, lon = True, O2_mr = 0.21*0.05, g = 9.81)
O3_col_One = O3_col(One_pc_h0.mean(dim='time'), time = True, lon = True, O2_mr = 0.21*0.01, g = 9.81)
O3_col_Zero5 = O3_col(Zero5_pc_h0.mean(dim='time'), time = True, lon = True, O2_mr = 0.21*0.005, g = 9.81)
O3_col_Zero1 = O3_col(Zero1_pc_h0.mean(dim='time'),time = True, lon = True, O2_mr = 0.21*0.001, g = 9.81)

#%% print mean, minimum and maximum ozone columns

print('\nP19 min, mean, max:')
print(P19_O3_col.min().values)
print(LWAV(P19_O3_col).values)
print(P19_O3_col.max().values)

print('\nW21 min, mean, max:')
print(W21_O3_col.min().values)
print(LWAV(W21_O3_col).values)
print(W21_O3_col.max().values)

print('\nW21 0.1% PAL min, mean, max:')
print(W21_01pc_O3_col.min().values)
print(LWAV(W21_01pc_O3_col).values)
print(W21_01pc_O3_col.max().values)

print('\nP19 no TL min, mean, max:')
print(W21_no_TL_O3_col.min().values)
print(LWAV(W21_no_TL_O3_col).values)
print(W21_no_TL_O3_col.max().values)

print('\nP19 0.1% PAL min, mean, max:')
print(P19_01pc_O3_col.min().values)
print(LWAV(P19_01pc_O3_col).values)
print(P19_01pc_O3_col.max().values)

print('\nP19 no TL min, mean, max:')
print(P19_no_TL_O3_col.min().values)
print(LWAV(P19_no_TL_O3_col).values)
print(P19_no_TL_O3_col.max().values)

#%% create effect for shading surrounding molecular features

def gradient_image(ax, direction=0.3, cmap_range=(0, 1), **kwargs):
    """
    Draw a gradient image based on a colormap.

    Parameters
    ----------
    ax : Axes
        The axes to draw on.
    direction : float
        The direction of the gradient. This is a number in
        range 0 (=vertical) to 1 (=horizontal).
    cmap_range : float, float
        The fraction (cmin, cmax) of the colormap that should be
        used for the gradient, where the complete colormap is (0, 1).
    **kwargs
        Other parameters are passed on to `.Axes.imshow()`.
        In particular, *cmap*, *extent*, and *transform* may be useful.
    """
    phi = direction * np.pi / 2
    v = np.array([np.cos(phi), np.sin(phi)])
    X = np.array([[v @ [1, 0], v @ [1, 1]],
                  [v @ [0, 0], v @ [0, 1]]])
    a, b = cmap_range
    X = a + (b - a) / X.max() * X
    im = ax.imshow(X, interpolation='bicubic', clim=(0, 1),
                   aspect='auto', **kwargs)
    return im


def gradient_bar(ax, x, y, width=0.5, bottom=0):
    for left, top in zip(x, y):
        right = left + width
        gradient_image(ax, extent=(left, right, bottom, top),
                       cmap=plt.cm.Blues_r, cmap_range=(0, 0.8))


fig, ax = plt.subplots()
ax.set(xlim=(0, 10), ylim=(0, 1))

#%% Direct imaging

markercolor = 'grey'
lw = 2

sp_path = '/Users/gregcooke/TP_1e_UV_paper_data/spectra_files/'

P19_EmSp = np.genfromtxt(sp_path+'TP-1_e_SSPO_Ideal_24.0hr_12.47pc_Wsrm2um-0320-01-01-00000.txt')
P19_10pc_EmSp = np.genfromtxt(sp_path+'TP-1_e_SSPO_10pc_o2_Ideal_24.0hr_12.47pc_Wsrm2um-0350-01-01-00000.txt')
P19_1pc_EmSp = np.genfromtxt(sp_path+ 'TP-1_e_SSPO_1pc_o2_Ideal_24.0hr_12.47pc_Wsrm2um-0310-01-01-00000.txt')
P19_01pc_EmSp = np.genfromtxt(sp_path+'TP-1_e_SSPO_0_1pc_o2_Ideal_24.0hr_12.47pc_Wsrm2um-0300-01-01-00000.txt')
P19_noTL_EmSp = np.genfromtxt(sp_path+'TP-1_e_SSPO_1day_Ideal_24.0hr_12.47pc_Wsrm2um-0280-01-01-00000.txt')

N2_P19_EmSp = np.genfromtxt(sp_path+'N2_TP-1_e_SSPO_Ideal_24.0hr_12.47pc_Wsrm2um-0320-01-01-00000.txt')
N2_P19_10pc_EmSp = np.genfromtxt(sp_path+'N2_TP-1_e_SSPO_10pc_o2_Ideal_24.0hr_12.47pc_Wsrm2um-0350-01-01-00000.txt')
N2_P19_1pc_EmSp = np.genfromtxt(sp_path+ 'N2_TP-1_e_SSPO_1pc_o2_Ideal_24.0hr_12.47pc_Wsrm2um-0310-01-01-00000.txt')
N2_P19_01pc_EmSp = np.genfromtxt(sp_path+'N2_TP-1_e_SSPO_0_1pc_o2_Ideal_24.0hr_12.47pc_Wsrm2um-0300-01-01-00000.txt')
N2_P19_noTL_EmSp = np.genfromtxt(sp_path+'N2_TP-1_e_SSPO_1day_Ideal_24.0hr_12.47pc_Wsrm2um-0280-01-01-00000.txt')

W21_EmSp = np.genfromtxt(sp_path+'TP-1_e_SSPO_MUSC_Ideal_24.0hr_12.47pc_Wsrm2um-0293-01-01-00000.txt')
W21_10pc_EmSp = np.genfromtxt(sp_path+'TP-1_e_MUSC_10pc_Ideal_24.0hr_12.47pc_Wsrm2um-0300-01-01-00000.txt')
W21_1pc_EmSp = np.genfromtxt(sp_path+'TP-1_e_MUSC_1pc_Ideal_24.0hr_12.47pc_Wsrm2um-0271-01-01-00000.txt')
W21_01pc_EmSp = np.genfromtxt(sp_path+'TP-1_e_MUSC_01pc_Ideal_24.0hr_12.47pc_Wsrm2um-0242-01-01-00000.txt')
W21_noTL_EmSp = np.genfromtxt(sp_path+'TP-1_e_MUSC_noTL_Ideal_24.0hr_12.47pc_Wsrm2um-0310-01-01-00000.txt')

N2_W21_EmSp = np.genfromtxt(sp_path+'N2_TP-1_e_SSPO_MUSC_Ideal_24.0hr_12.47pc_Wsrm2um-0293-01-01-00000.txt')
N2_W21_10pc_EmSp = np.genfromtxt(sp_path+'N2_TP-1_e_MUSC_10pc_Ideal_24.0hr_12.47pc_Wsrm2um-0300-01-01-00000.txt')
N2_W21_1pc_EmSp = np.genfromtxt(sp_path+'N2_TP-1_e_MUSC_1pc_Ideal_24.0hr_12.47pc_Wsrm2um-0271-01-01-00000.txt')
N2_W21_01pc_EmSp = np.genfromtxt(sp_path+'N2_TP-1_e_MUSC_01pc_Ideal_24.0hr_12.47pc_Wsrm2um-0242-01-01-00000.txt')
N2_W21_noTL_EmSp = np.genfromtxt(sp_path+'N2_TP-1_e_MUSC_noTL_Ideal_24.0hr_12.47pc_Wsrm2um-0310-01-01-00000.txt')

gs = gridspec.GridSpec(1,2)
plt.figure(figsize = (13,4.5))
gs.update(wspace = 0.2)
ax = plt.subplot(gs[0,0])
ax.spines['top'].set_visible(False)

plt.plot(P19_EmSp[:,0], P19_EmSp[:,4], lw = lw, color = P19_color, label = 'P19 PI')
plt.plot(P19_10pc_EmSp[:,0], P19_10pc_EmSp[:,4], lw = lw, color = P19_10pc_color, label = 'P19 10% PAL')
plt.plot(P19_1pc_EmSp[:,0], P19_1pc_EmSp[:,4], lw = lw, color = P19_1pc_color, label = 'P19 1% PAL')
plt.plot(P19_01pc_EmSp[:,0], P19_01pc_EmSp[:,4], lw = lw, color = P19_01pc_color, label = 'P19 0.1% PAL')
plt.plot(P19_noTL_EmSp[:,0], P19_noTL_EmSp[:,4], lw = lw, color = P19_no_TL_color, label = 'P19 noTL')

plt.plot(W21_EmSp[:,0], W21_EmSp[:,4], lw = lw, color = W21_color, label = 'W21 PI')
plt.plot(W21_10pc_EmSp[:,0], W21_10pc_EmSp[:,4], lw = lw, color = W21_10pc_color, label = 'W21 10% PAL')
plt.plot(W21_1pc_EmSp[:,0], W21_1pc_EmSp[:,4], lw = lw, color = W21_1pc_color, label = 'W21 1% PAL')
plt.plot(W21_01pc_EmSp[:,0], W21_01pc_EmSp[:,4], lw = lw, color = W21_01pc_color, label = 'W21 0.1% PAL')
plt.plot(W21_noTL_EmSp[:,0], W21_noTL_EmSp[:,4], lw = lw, color = W21_no_TL_color, label = 'W21 noTL')

plt.xlim(4.6,4.9)
plt.ylim(0, 0.3e-8)

x1 = 4.65; x2 = 0.3

gradient_image(ax, direction=0, extent=((4.67-4.6)/0.3, (4.86-4.6)/0.3, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.2)

gradient_image(ax, direction=0, extent=((4.67-4.6)/0.3, (5-4.6)/0.3, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.2)


plt.text(4.7,3.05e-9,'O'+r'$_3$', color = markercolor, fontsize = Fontsize, horizontalalignment = 'center')
plt.text(4.88,3.05e-9,'CO'+r'$_2$', color = markercolor, fontsize = Fontsize, horizontalalignment = 'center')

plt.tick_params(labelsize = 15, labelleft = True, width = 2, length = 4, left = True, right = True)
plt.xlabel('Wavelength ['+r'$\rm \mu$'+'m]', fontsize = 15)
plt.ylabel('Spectral radiance [Wsr'+sup(-1)+'m'+sup(-2)+r'$\rm \mu$'+'m'+sup(-1)+']', fontsize = 15)


JWST_color1 = '#FFC500'
JWST_color2 = '#FF4BEC'
plt.axvspan(4.6,4.9,0,0.1,color = JWST_color1 , alpha = 0.3)
plt.text(4.6+(0.1*(4.9-4.6)),0.5e-10,'LIFE', color = 'k', fontsize = 15, horizontalalignment = 'center')
plt.text(4.6+(0.5*(4.9-4.6)),0.5e-10,'JWST NIRSpec', color = 'k', fontsize = 15, horizontalalignment = 'center')
plt.text(4.6+(0.85*(4.9-4.6)),0.5e-10,'ELT METIS', color = 'k', fontsize = 15, horizontalalignment = 'center')


ax = plt.subplot(gs[0,1])
ax.spines['top'].set_visible(False)

a, = plt.plot(P19_EmSp[:,0], P19_EmSp[:,4], lw = lw, color = P19_color, label = 'P19 PI')
b, = plt.plot(P19_10pc_EmSp[:,0], P19_10pc_EmSp[:,4], lw = lw, color = P19_10pc_color, label = 'P19 10% PAL')
c, = plt.plot(P19_1pc_EmSp[:,0], P19_1pc_EmSp[:,4], lw = lw, color = P19_1pc_color, label = 'P19 1% PAL')
d, = plt.plot(P19_01pc_EmSp[:,0], P19_01pc_EmSp[:,4], lw = lw, color = P19_01pc_color, label = 'P19 0.1% PAL')
e, = plt.plot(P19_noTL_EmSp[:,0], P19_noTL_EmSp[:,4], lw = lw, color = P19_no_TL_color, label = 'P19 noTL')

f, = plt.plot(W21_EmSp[:,0], W21_EmSp[:,4], lw = lw, color = W21_color, label = 'W21 PI')
g, = plt.plot(W21_10pc_EmSp[:,0], W21_10pc_EmSp[:,4], lw = lw, color = W21_10pc_color, label = 'W21 10% PAL')
h, = plt.plot(W21_1pc_EmSp[:,0], W21_1pc_EmSp[:,4], lw = lw, color = W21_1pc_color, label = 'W21 1% PAL')
i, = plt.plot(W21_01pc_EmSp[:,0], W21_01pc_EmSp[:,4], lw = lw, color = W21_01pc_color, label = 'W21 0.1% PAL')
j, = plt.plot(W21_noTL_EmSp[:,0], W21_noTL_EmSp[:,4], lw = lw, color = W21_no_TL_color, label = 'W21 noTL')

plt.xlim(9,10.2)
plt.ylim(0, 2e-9)
x1 = 9; x2 = 1.2
plt.text(9.65,2.05e-9,'O'+r'$_3$', color = markercolor, fontsize = Fontsize, horizontalalignment = 'center')
plt.tick_params(labelsize = 15, labelleft = True, width = 2, length = 4, left = True, right = True)

legend1 = plt.legend(loc = (0.01,0.1),handles=[a,b,c,d,e], handlelength = 0.5, fontsize = 13, frameon = False, ncol = 1)
text_colors = [P19_color,P19_10pc_color,P19_1pc_color,P19_01pc_color, P19_no_TL_color]  # Define text colors for each label
for text, color_ in zip(legend1.get_texts(), text_colors):
    text.set_color(color_)

legend = plt.legend(loc = (0.01,0.75), handles=[f,g,h,i,j],  handlelength = 0.5, fontsize = 13, frameon = False, ncol = 2)
text_colors = [W21_color,W21_10pc_color,W21_1pc_color,W21_01pc_color, W21_no_TL_color]  # Define text colors for each label
for text, color_ in zip(legend.get_texts(), text_colors):
    text.set_color(color_)
    
plt.gca().add_artist(legend)      

plt.gca().add_artist(legend1)       
    

plt.xlabel('Wavelength ['+r'$\rm \mu$'+'m]', fontsize = 15)

plt.axvspan(9,10.2,0,0.1,color = JWST_color2, alpha = 0.3)
plt.text(9+(0.1*(10.2-9)),0.05e-9,'LIFE', color = 'k', fontsize = 15, horizontalalignment = 'center')
plt.text(9+(0.5*(10.2-9)),0.05e-9,'JWST MIRI', color = 'k', fontsize = 15, horizontalalignment = 'center')
plt.text(9+(0.85*(10.2-9)),0.05e-9,'ELT METIS', color = 'k', fontsize = 15, horizontalalignment = 'center')

gradient_image(ax, direction=0, extent=((9-9)/1.2, (10.2-9)/1., 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.1)

plt.savefig('/Users/gregcooke/python_output/Direct_imaging.png', bbox_inches = 'tight', dpi = 400)

#%% Figure out amount of radiation in each band
#%% Plot ozone columns

levels = 21

factor = 2.6867e20
factor = 1
cmap  ='inferno'
plt.figure(figsize = (12, 10* 1.3*1.2))
gs = gridspec.GridSpec(5 ,2)
gs.update(hspace = 0.12)
gs.update(wspace = 0.1)
ax = plt.subplot(gs[0,0], projection=ccrs.Robinson(central_longitude=180))
model = (P19_O3_col*factor).plot.contourf(cmap = cmap, add_labels = False, add_colorbar = False, vmin = 600*factor, vmax = 3000*factor, levels = levels, transform=ccrs.PlateCarree())
plt.tick_params(labelsize = 15, labelbottom = False, width = 2, length = 4, top = True)
cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]')
plt.title('P19 PI', y = 0.5, fontsize = Fontsize, color = 'w')

ax = plt.subplot(gs[1,0], projection=ccrs.Robinson(central_longitude=180))
model = (P19_10pc_O3_col*factor).plot.contourf(cmap = cmap, add_labels = False, add_colorbar = False, vmin = 200*factor, vmax = 1500*factor, levels = levels, transform=ccrs.PlateCarree())
plt.tick_params(labelsize = 15, labelbottom = False, width = 2, length = 4, top = True)
cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]')
plt.title('P19 10% PAL', y = 0.5, fontsize = Fontsize, color = 'w')

ax = plt.subplot(gs[2,0], projection=ccrs.Robinson(central_longitude=180))
model = (P19_1pc_O3_col*factor).plot.contourf(cmap = cmap, add_labels = False, add_colorbar = False, vmin = 200*factor, vmax = 700*factor, levels = levels, transform=ccrs.PlateCarree())
plt.tick_params(labelsize = 15, labelbottom = False, width = 2, length = 4, top = True)
cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]')
plt.title('P19 1% PAL', y = 0.5, fontsize = Fontsize, color = 'w')

ax = plt.subplot(gs[3,0], projection=ccrs.Robinson(central_longitude=180))
model = (P19_01pc_O3_col*factor).plot.contourf(cmap = cmap, add_labels = False, add_colorbar = False, vmin = 200*factor, vmax = 700*factor, levels = levels, transform=ccrs.PlateCarree())
plt.tick_params(labelsize = 15, labelbottom = False, width = 2, length = 4, top = True)
cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]')
plt.title('P19 0.1% PAL', y = 0.5, fontsize = Fontsize, color = 'w')

ax = plt.subplot(gs[4,0], projection=ccrs.Robinson(central_longitude=180))
model = (P19_no_TL_O3_col*factor).plot.contourf(cmap = cmap, add_labels = False, add_colorbar = False, vmin = 600*factor, vmax = 3000*factor, levels = levels, transform=ccrs.PlateCarree())
plt.tick_params(labelsize = 15, labelbottom = False, width = 2, length = 4, top = True)
cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]')
plt.title('P19 noTL', y = 0.5, fontsize = Fontsize, color = 'w')

ax = plt.subplot(gs[0,1], projection=ccrs.Robinson(central_longitude=180))
model = (W21_O3_col*factor).plot.contourf(cmap = cmap, add_labels = False, add_colorbar = False, vmin = 20*factor, vmax = 120*factor, levels = levels, transform=ccrs.PlateCarree())
plt.tick_params(labelsize = 15, labelbottom = False, width = 2, length = 4, top = True)
cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]')
plt.title('W21 PI', y = 0.5, fontsize = Fontsize, color = 'w')

ax = plt.subplot(gs[1,1], projection=ccrs.Robinson(central_longitude=180))
model = (W21_10pc_O3_col*factor).plot.contourf(cmap = cmap, add_labels = False, add_colorbar = False, vmin = 200*factor, vmax = 700*factor, levels = levels, transform=ccrs.PlateCarree())
plt.tick_params(labelsize = 15, labelbottom = False, width = 2, length = 4, top = True)
cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]')
plt.title('W21 10% PAL', y = 0.5, fontsize = Fontsize, color = 'w')

ax = plt.subplot(gs[2,1], projection=ccrs.Robinson(central_longitude=180))
model = (W21_1pc_O3_col*factor).plot.contourf(cmap = cmap, add_labels = False, add_colorbar = False, vmin = 600*factor, vmax = 2000*factor, levels = levels, transform=ccrs.PlateCarree())
plt.tick_params(labelsize = 15, labelbottom = False, width = 2, length = 4, top = True)
cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]')
plt.title('W21 1% PAL', y = 0.5, fontsize = Fontsize, color = 'w')

ax = plt.subplot(gs[3,1], projection=ccrs.Robinson(central_longitude=180))
model = (W21_01pc_O3_col*factor).plot.contourf(cmap = cmap, add_labels = False, add_colorbar = False, vmin = 800*factor, vmax = 2000*factor, levels = levels, transform=ccrs.PlateCarree())
plt.tick_params(labelsize = 15, labelbottom = False, width = 2, length = 4, top = True)
cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]')
plt.title('W21 0.1% PAL', y = 0.5, fontsize = Fontsize, color = 'w')

ax = plt.subplot(gs[4,1], projection=ccrs.Robinson(central_longitude=180))
model = (W21_no_TL_O3_col*factor).plot.contourf(cmap = cmap, add_labels = False, add_colorbar = False, levels =levels, transform=ccrs.PlateCarree())
plt.tick_params(labelsize = 15, width = 2, length = 4, top = True)
cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]')
plt.xlabel('Longitude ['+sup('\circ')+']', fontsize = Fontsize)
plt.ylabel('Latitude ['+sup('\circ')+']', fontsize = Fontsize)
plt.title('W21 noTL', y = 0.5, fontsize = Fontsize, color = 'w')

plt.savefig('/Users/gregcooke/python_output/O3_col.png', dpi = 400, bbox_inches = 'tight')

#%% O2 photolyis plot

k = 1.381e-23
ylim = 6e-6

lw = 2

def JO2(DS, color = 'k'):
    plt.plot(2*LWAV((DS.jo2_a+DS.jo2_b)*DS.O2*DS.lev*100/(k*DS.T)), P19.lev.values, color = color, lw = lw)

def JO3(DS, color = 'k'):
    plt.plot(2*LWAV((DS.jo3_a+DS.jo3_b)*DS.O3*DS.lev*100/(k*DS.T)), P19.lev.values, color = color, lw = lw)

def n_dens_plot(DS, var, color = 'k'):
    plt.plot(LWAV(DS[var]*DS.lev*100/(k*DS.T)), P19.lev.values, color = color, lw = lw)

def n_dens_plot_Ox(DS, var, color = 'k'):
    plt.plot(LWAV((DS['O3']+DS['O'])*DS.lev*100/(k*DS.T)), P19.lev.values, color = color, lw = lw)

plt.figure(figsize = (12,4))
gs = gridspec.GridSpec(1,3)

plt.subplot(gs[0,0])
JO2(P19, color = P19_color)
JO2(P19_10pc, color = P19_10pc_color)
JO2(P19_1pc, color = P19_1pc_color)
JO2(P19_01pc, color = P19_01pc_color)
JO2(P19_no_TL, color = P19_no_TL_color)

JO2(W21, color = W21_color)
JO2(W21_10pc, color = W21_10pc_color)
JO2(W21_1pc, color = W21_1pc_color)
JO2(W21_01pc, color = W21_01pc_color)
JO2(W21_no_TL, color = W21_no_TL_color)

plt.xlabel('Production of O', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, labelleft = True)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xlim(1e6, 3e12)
plt.xscale('log')

plt.subplot(gs[0,1])

n_dens_plot(P19, var = 'O', color = P19_color)
n_dens_plot(P19_10pc, var = 'O', color = P19_10pc_color)
n_dens_plot(P19_1pc, var = 'O', color = P19_1pc_color)
n_dens_plot(P19_01pc, var = 'O', color = P19_01pc_color)
n_dens_plot(P19_no_TL, var = 'O', color =P19_no_TL_color)

n_dens_plot(W21, var = 'O', color = W21_color)
n_dens_plot(W21_10pc, var = 'O', color = W21_10pc_color)
n_dens_plot(W21_1pc, var = 'O', color = W21_1pc_color)
n_dens_plot(W21_01pc, var = 'O', color = W21_01pc_color)
n_dens_plot(W21_no_TL, var = 'O', color = W21_no_TL_color)

plt.xlabel('O number density', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, labelleft = False)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xlim(1e6, 3e18)
plt.xscale('log')

plt.subplot(gs[0,2])

n_dens_plot(P19, var = 'O3', color = P19_color)
n_dens_plot(P19_10pc, var = 'O3', color = P19_10pc_color)
n_dens_plot(P19_1pc, var = 'O3', color = P19_1pc_color)
n_dens_plot(P19_01pc, var = 'O3', color = P19_01pc_color)
n_dens_plot(P19_no_TL, var = 'O3', color =P19_no_TL_color)

n_dens_plot(W21, var = 'O3', color = W21_color)
n_dens_plot(W21_10pc, var = 'O3', color = W21_10pc_color)
n_dens_plot(W21_1pc, var = 'O3', color = W21_1pc_color)
n_dens_plot(W21_01pc, var = 'O3', color = W21_01pc_color)
n_dens_plot(W21_no_TL, var = 'O3', color = W21_no_TL_color)

plt.xlabel('O'+sub(3)+' number density', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, labelleft = False)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xlim(1e16, 1e20)
plt.xscale('log')

fs = 12
x = 1.5e16
plt.text(x, 2e-5, 'P19',color = P19_color, fontsize = fs)
plt.text(x, 7e-5, 'P19 10% PAL',color = P19_10pc_color, fontsize = fs)
plt.text(x, 2e-4, 'P19 1% PAL',color = P19_1pc_color, fontsize = fs)
plt.text(x, 7e-4, 'P19 0.1% PAL',color = P19_01pc_color, fontsize = fs)
plt.text(x, 2e-3, 'P19 noTL',color = P19_no_TL_color, fontsize = fs)
x = 1e18
plt.text(x, 2e-5, 'W21',color = W21_color, fontsize = fs)
plt.text(x, 7e-5, 'W21 10% PAL',color = W21_10pc_color, fontsize = fs)
plt.text(x, 2e-4, 'W21 1% PAL',color = W21_1pc_color, fontsize = fs)
plt.text(x, 7e-4, 'W21 0.1% PAL',color = W21_01pc_color, fontsize = fs)
plt.text(x, 2e-3, 'W21 noTL',color = W21_no_TL_color, fontsize = fs)

plt.figure(figsize = (12,4))
gs = gridspec.GridSpec(1,3)

plt.subplot(gs[0,0])
JO2(P19, color = P19_color)
JO2(P19_10pc, color = P19_10pc_color)
JO2(P19_1pc, color = P19_1pc_color)
JO2(P19_01pc, color = P19_01pc_color)
JO2(P19_no_TL, color = P19_no_TL_color)

JO2(W21, color = W21_color)
JO2(W21_10pc, color = W21_10pc_color)
JO2(W21_1pc, color = W21_1pc_color)
JO2(W21_01pc, color = W21_01pc_color)
JO2(W21_no_TL, color = W21_no_TL_color)

plt.xlabel('Production of O', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, labelleft = True)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xlim(1e6, 3e12)
plt.xscale('log')

plt.subplot(gs[0,1])

JO3(P19, color = P19_color)
JO3(P19_10pc, color = P19_10pc_color)
JO3(P19_1pc, color = P19_1pc_color)
JO3(P19_01pc, color = P19_01pc_color)
JO3(P19_no_TL, color = P19_no_TL_color)

JO3(W21, color = W21_color)
JO3(W21_10pc, color = W21_10pc_color)
JO3(W21_1pc, color = W21_1pc_color)
JO3(W21_01pc, color = W21_01pc_color)
JO3(W21_no_TL, color = W21_no_TL_color)

plt.xlabel('Destruction of ozone', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, labelleft = False)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xlim(1e10, 1e15)
plt.xscale('log')

fs = 12

x = 2e10
plt.text(x, 2e-5, 'P19',color = P19_color, fontsize = fs)
plt.text(x, 7e-5, 'P19 10% PAL',color = P19_10pc_color, fontsize = fs)
plt.text(x, 2e-4, 'P19 1% PAL',color = P19_1pc_color, fontsize = fs)
plt.text(x, 7e-4, 'P19 0.1% PAL',color = P19_01pc_color, fontsize = fs)
plt.text(9e10, 2e-3, 'P19 noTL',color = P19_no_TL_color, fontsize = fs)
x = 5e12
plt.text(x, 2e-5, 'W21',color = W21_color, fontsize = fs)
plt.text(x, 7e-5, 'W21 10% PAL',color = W21_10pc_color, fontsize = fs)
plt.text(x, 2e-4, 'W21 1% PAL',color = W21_1pc_color, fontsize = fs)
plt.text(x, 7e-4, 'W21 0.1% PAL',color = W21_01pc_color, fontsize = fs)
plt.text(x, 2e-3, 'W21 noTL',color = W21_no_TL_color, fontsize = fs)

plt.subplot(gs[0,2])

n_dens_plot_Ox(P19, var = 'O3', color = P19_color)
n_dens_plot_Ox(P19_10pc, var = 'O3', color = P19_10pc_color)
n_dens_plot_Ox(P19_1pc, var = 'O3', color = P19_1pc_color)
n_dens_plot_Ox(P19_01pc, var = 'O3', color = P19_01pc_color)
n_dens_plot_Ox(P19_no_TL, var = 'O3', color =P19_no_TL_color)

n_dens_plot_Ox(W21, var = 'O3', color = W21_color)
n_dens_plot_Ox(W21_10pc, var = 'O3', color = W21_10pc_color)
n_dens_plot_Ox(W21_1pc, var = 'O3', color = W21_1pc_color)
n_dens_plot_Ox(W21_01pc, var = 'O3', color = W21_01pc_color)
n_dens_plot_Ox(W21_no_TL, var = 'O3', color = W21_no_TL_color)

plt.xlabel('O'+sub(3)+' number density', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, labelleft = False)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xlim(1e16, 1e20)
plt.xscale('log')
#%% Temperature and mixing ratio plot

def MR(var = 'H2O'):
    plt.plot(LWAV(P19[var]), P19.lev.values, color = P19_color, lw = lw, label = 'P19')
    plt.plot(LWAV(P19_10pc[var]), P19.lev.values, color = P19_10pc_color, lw = lw, label = 'P19 10% PAL')
    plt.plot(LWAV(P19_1pc[var]), P19.lev.values, color = P19_1pc_color, lw = lw, label = 'P19 1% PAL')
    plt.plot(LWAV(P19_01pc[var]), P19.lev.values, color = P19_01pc_color, lw = lw,label = 'P19\n0.1% PAL')
    plt.plot(LWAV(P19_no_TL[var]), P19_no_TL.lev.values, lw = lw,color = P19_no_TL_color, label = 'P19 noTL')
    plt.plot(LWAV(W21[var]), P19.lev.values, color = W21_color,lw = lw, label = 'W21')
    plt.plot(LWAV(W21_10pc[var]), P19.lev.values, color = W21_10pc_color,lw = lw, label = 'W21 10% PAL')
    plt.plot(LWAV(W21_1pc[var]), P19.lev.values, color = W21_1pc_color,lw = lw, label = 'W21 1% PAL')
    plt.plot(LWAV(W21_01pc[var]), P19.lev.values, color = W21_01pc_color,lw = lw, label = 'W21\n0.1% PAL')
    plt.plot(LWAV(W21_no_TL[var]), P19.lev.values, color = W21_no_TL_color,lw = lw, label = 'W21\nno TL')

lw = 2
ylim = 6e-6

plt.figure(figsize = (11*1.2,11*1.3*1.2))
gs = gridspec.GridSpec(4,2)
gs.update(hspace = 0.28)
gs.update(wspace = 0.12)
plt.subplot(gs[0,0])
plt.axhline(55, lw = 2, ls = ':', color = 'k')
plt.text(270, 0.06, 'Thermosphere', fontsize = 14, horizontalalignment = 'center')
plt.text(270, 10, 'Middle\natmosphere', fontsize = 14, horizontalalignment = 'center')
plt.text(270, 400, 'Troposphere', fontsize = 14, horizontalalignment = 'center')
plt.axhline(0.1, lw = 2, ls = ':', color = 'k')
plt.plot(LWAV(P19.T), P19.lev.values, color = P19_color, lw = lw, label = 'P19')
plt.plot(LWAV(P19_10pc.T), P19.lev.values, color = P19_10pc_color, lw = lw,label = 'P19 10% PAL')
plt.plot(LWAV(P19_1pc.T), P19.lev.values, color = P19_1pc_color, lw = lw,label = 'P19 1% PAL')
plt.plot(LWAV(P19_01pc.T), P19.lev.values, color = P19_01pc_color, lw = lw,label = 'P19 0.1% PAL')
plt.plot(LWAV(P19_no_TL.T), P19_no_TL.lev.values, lw = lw,color = P19_no_TL_color, label = 'P19 noTL')
plt.plot(LWAV(W21.T), P19.lev.values, color = W21_color,lw = lw, label = 'W21')
plt.plot(LWAV(W21_10pc.T), P19.lev.values, color = W21_10pc_color, lw = lw, label = 'W21 10% PAL')
plt.plot(LWAV(W21_1pc.T), P19.lev.values, color = W21_1pc_color, lw = lw, label = 'W21 1% PAL')
plt.plot(LWAV(W21_01pc.T), P19.lev.values, color = W21_01pc_color, lw = lw, label = 'W21 0.1% PAL')
plt.plot(LWAV(W21_no_TL.T), P19_no_TL.lev.values, lw = lw,color = P19_no_TL_color, label = 'P19 noTL')
plt.ylabel('Pressure [hPa]', fontsize = 15)
plt.xlabel('Temperature [K]', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, width = 2, length = 4)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xlim(170,300)

plt.text(175, 4e-5, 'a', color = 'k', fontsize = 15)

k = 1.381e-23

plt.subplot(gs[0,1])
MR(var = 'O')
plt.xlabel('O mixing ratio', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, labelleft = False, width = 2, length = 4)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xscale('log')
plt.text(1e-19, 4e-5, 'b', color = 'k', fontsize = 15)
x = 8e-1
plt.text(x, 0.8e0, 'W21 PI',color = W21_color, fontsize = 14, horizontalalignment = 'right')
plt.text(x, 4e0, 'W21 10% PAL',color = W21_10pc_color, fontsize = 14, horizontalalignment = 'right')
plt.text(x, 2e1, 'W21 1% PAL',color = W21_1pc_color, fontsize = 14, horizontalalignment = 'right')
plt.text(x, 1.2e2, 'W21 0.1% PAL',color = W21_01pc_color, fontsize = 14, horizontalalignment = 'right')
plt.text(x, 8e2, 'W21 noTL',color = W21_no_TL_color, fontsize = 14, horizontalalignment = 'right')
x = 1e-17
plt.text(x, 0.0009, 'P19 PI',color = P19_color, fontsize = 15)
plt.text(x, 0.005, 'P19 10% PAL',color = P19_10pc_color, fontsize = 15)
plt.text(x, 0.025, 'P19 1% PAL',color = P19_1pc_color, fontsize = 15)
plt.text(x, 0.15, 'P19 0.1% PAL',color = P19_01pc_color, fontsize = 15)
plt.text(x, 1, 'P19 noTL',color = P19_no_TL_color, fontsize = 15)


plt.subplot(gs[1,0])
MR(var = 'O2')
plt.xlabel('O'+sub(2)+' mixing ratio', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, width = 2, length = 4)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e-6)
plt.ylabel('Pressure [hPa]', fontsize = 15)
plt.text(1.5e-6, 500, 'c', color = 'k', fontsize = 15)

plt.subplot(gs[1,1])
MR(var = 'O3')
plt.xlabel('O'+sub(3)+' mixing ratio', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, labelleft = False, width = 2, length = 4)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e-9, 2e-4)
plt.text(1e-4, 500, 'd', color = 'k', fontsize = 15)

plt.subplot(gs[2,0])
MR(var = 'H2O')
plt.ylabel('Pressure [hPa]', fontsize = 15)
plt.xlabel('H'+sub(2)+'O mixing ratio', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, width = 2, length = 4)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xscale('log')
plt.text(1e-12, 500, 'e', color = 'k', fontsize = 15)

k = 1.381e-23

plt.subplot(gs[2,1])
MR(var = 'CH4')
plt.xlabel('CH'+sub(4)+' mixing ratio', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, labelleft = False, width = 2, length = 4)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xscale('log')
#plt.legend(loc = (0.01,0.01), fontsize = 12, handlelength = 0.5, frameon = False, ncol = 2)
plt.xlim(1e-22)
plt.text(1.4e-6, 4e-5, 'f', color = 'k', fontsize = 15)

plt.subplot(gs[3,0])
MR(var = 'N2O')
plt.xlabel('N'+sub(2)+'O mixing ratio', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, width = 2, length = 4)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xscale('log')
plt.ylabel('Pressure [hPa]', fontsize = 15)
plt.text(5e-14, 500, 'g', color = 'k', fontsize = 15)

plt.subplot(gs[3,1])
plt.plot(LWAV(P19_h2.CO2), P19.lev.values, color = P19_color,lw = lw, label = 'P19')
plt.plot(LWAV(P19_10pc_h2.CO2), P19_01pc.lev.values, color = P19_10pc_color,lw = lw, label = 'P19 10% PAL')
plt.plot(LWAV(P19_1pc_h2.CO2), P19_01pc.lev.values, color = P19_1pc_color,lw = lw, label = 'P19 1% PAL')
plt.plot(LWAV(P19_01pc_h2.CO2), P19_01pc.lev.values, color = P19_01pc_color,lw = lw, label = 'P19 0.1% PAL')
plt.plot(LWAV(P19_no_TL_h2.CO2), P19_01pc.lev.values, color = P19_no_TL_color, lw = lw,label = 'P19 noTL')
plt.plot(LWAV(W21_h2.CO2), P19_no_TL.lev.values, color = W21_color,lw = lw,label = 'W21')
plt.plot(LWAV(W21_10pc_h2.CO2), P19_no_TL.lev.values, color = W21_10pc_color,lw = lw,label = 'W21 10% PAL')
plt.plot(LWAV(W21_1pc_h2.CO2), P19_no_TL.lev.values, color = W21_1pc_color,lw = lw,label = 'W21 1% PAL')
plt.plot(LWAV(W21_01pc_h2.CO2), P19_no_TL.lev.values, color = W21_01pc_color,lw = lw,label = 'W21 0.1% PAL')
plt.plot(LWAV(W21_no_TL_h2.CO2), P19_no_TL.lev.values, color = W21_no_TL_color,lw = lw,label = 'W21 noTL')

plt.text(5e-4, 500, 'h', color = 'k', fontsize = 15)

plt.xlabel('CO'+sub(2)+' mixing ratio', fontsize = 15)
plt.tick_params(labelsize = 15, right = True, labelleft = False, width = 2, length = 4)
plt.ylim(1e3,ylim)
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e-8, 1e-3)

plt.savefig('/Users/gregcooke/python_output/Temp_mixing_ratio.eps', format = 'eps', bbox_inches = 'tight', dpi = 400)
plt.savefig('/Users/gregcooke/python_output/Temp_mixing_ratio.png', bbox_inches = 'tight', dpi = 400)

#%%  Cumulative ozone column function definition

def cumulative_col(ds, time = False, var = 'O2', O2_mr = 0.21, lon = True, g = 9.80616):
    
    g = g # gravity acceleration (cm/sec2)
    N2_mr = 1-O2_mr #nitrogen mixing ratio
    m_avg = ((28*N2_mr)+(O2_mr*32))*1.661e-27 # average weight of atmosphere
    
    #calcualte column integrals 
    var = ds[var]
    
    # variables to calcualate hybrid pressure on model interfaces
    ps = ds['PS'].values # surface pressure
    p0 = ds['P0'].values
    hyai = ds['hyai'].values
    hybi = ds['hybi'].values
    

    #col = np.ndarray(ps.shape, np.float32)
    dp = np.ndarray(var.shape, np.float32)
    
    nt = 1
    try:
        nt = ds['time'].size
    except:
        nt = 0
    ny = ds['lat'].size
    nz = ds['lev'].size
    nx = ds['lon'].size
    
    press = np.ndarray(nz, np.float32)
    if (nt > 0):
        for i in range(nt):
            for j in range(ny):
                for k in range(nx):
                    press = hyai * p0 + hybi * ps[i,j,k]
                    for iz in range(nz):
                        dp[i,iz,j,k] = (press[iz+1]-press[iz])
    else:
        for j in range(ny):
            for k in range(nx):
                press = hyai * p0 + hybi * ps[j,k]
                for iz in range(nz):
                    dp[iz,j,k] = (press[iz+1]-press[iz])
                    
    #calcualte column integrals 
    var_col = var * dp / (m_avg * g)
    
    #o3col = np.sum(o3, axis = 0)  # convert to DU
    if (lon == False):
        var_col  = var_col.mean(dim='lon')
    if (time == False):
        try:
            var_col  = var_col.mean(dim='time')
        except:
            var_col  = var_col
    print('Columns calculated')
    return var_col

#%% Plot surface temperature
levels = 12
yplace = 1

TS_min_values = np.array([P19.TS.values.min(), P19_no_TL.TS.values.min(), 
                 P19_10pc.TS.values.min(), P19_1pc.TS.values.min(), 
                 P19_01pc.TS.values.min(), W21.TS.values.min(), 
                 W21_no_TL.TS.values.min(), W21_10pc.TS.values.min(),
                 W21_1pc.TS.values.min(), W21_01pc.TS.values.min()]) 

print(str(TS_min_values.min()) + ' --> ' + str(TS_min_values.max()))

TS_max_values = np.array([P19.TS.values.max(), P19_no_TL.TS.values.max(), 
                 P19_10pc.TS.values.max(), P19_1pc.TS.values.max(), 
                 P19_01pc.TS.values.max(), W21.TS.values.max(), 
                 W21_no_TL.TS.values.max(), W21_10pc.TS.values.max(),
                 W21_1pc.TS.values.max(), W21_01pc.TS.values.max()]) 

print(str(TS_max_values.min()) + ' --> ' + str(TS_max_values.max()))

TS_mean_values = np.array([LWAV(P19.TS).values, 
                           LWAV(P19_no_TL.TS).values,
                           LWAV(P19_10pc.TS).values,
                           LWAV(P19_1pc.TS).values,
                           LWAV(P19_01pc.TS).values,
                           LWAV(W21.TS).values,
                           LWAV(W21_no_TL.TS).values,
                           LWAV(W21_10pc.TS).values,
                           LWAV(W21_1pc.TS).values,
                           LWAV(W21_01pc.TS).values]) 

print(str(TS_mean_values.min()) + ' --> ' + str(TS_mean_values.max()))

Fontsize = 15
Figure = plt.figure(figsize = (12, 10* 1.3*1.2))
gs = gridspec.GridSpec(5 ,2)
gs.update(hspace = 0.15)
gs.update(wspace = 0.0)
ax = plt.subplot(gs[0,0], projection=ccrs.Robinson(central_longitude=180))
model = P19.TS.mean(dim='time').plot.contourf(cmap = 'inferno', add_labels = False, add_colorbar = False, levels = levels, vmax = 280,vmin=170, transform=ccrs.PlateCarree())
plt.title('P19 PI', y = yplace, fontsize = Fontsize, color = 'k')

ax = plt.subplot(gs[1,0], projection=ccrs.Robinson(central_longitude=180))
model = P19_no_TL.TS.mean(dim='time').plot.contourf(cmap = 'inferno', add_labels = False, add_colorbar = False, levels = levels, vmax = 280,vmin=170, transform=ccrs.PlateCarree())
plt.title('P19 noTL', y = yplace, fontsize = Fontsize, color = 'k')

ax = plt.subplot(gs[2,0], projection=ccrs.Robinson(central_longitude=180))
model = P19_10pc.TS.mean(dim='time').plot.contourf(cmap = 'inferno', add_labels = False, add_colorbar = False, levels = levels, vmax = 280,vmin=170,transform=ccrs.PlateCarree())
plt.title('P19 10% PAL', y = yplace, fontsize = Fontsize, color = 'k')

ax = plt.subplot(gs[3,0], projection=ccrs.Robinson(central_longitude=180))
model = P19_1pc.TS.mean(dim='time').plot.contourf(cmap = 'inferno', add_labels = False, add_colorbar = False, levels = levels, vmax = 280,vmin=170,transform=ccrs.PlateCarree())
plt.title('P19 1% PAL', y = yplace, fontsize = Fontsize, color = 'k')

ax = plt.subplot(gs[4,0], projection=ccrs.Robinson(central_longitude=180))
model = P19_01pc.TS.mean(dim='time').plot.contourf(cmap = 'inferno', add_labels = False, add_colorbar = False, levels = levels, vmax = 280,vmin=170,transform=ccrs.PlateCarree())
plt.title('P19 0.1% PAL', y = yplace, fontsize = Fontsize, color = 'k')

ax = plt.subplot(gs[0,1], projection=ccrs.Robinson(central_longitude=180))
model = W21.TS.mean(dim='time').plot.contourf(cmap = 'inferno', add_labels = False, add_colorbar = False, levels = levels, vmax = 280,vmin=170,transform=ccrs.PlateCarree())
plt.title('W21 PI', y = yplace, fontsize = Fontsize, color = 'k')

ax = plt.subplot(gs[1,1], projection=ccrs.Robinson(central_longitude=180))
model = W21_no_TL.TS.mean(dim='time').plot.contourf(cmap = 'inferno', add_labels = False, add_colorbar = False, levels = levels, vmax = 280,vmin=170,transform=ccrs.PlateCarree())
plt.title('W21 noTL', y = yplace, fontsize = Fontsize, color = 'k')

ax = plt.subplot(gs[2,1], projection=ccrs.Robinson(central_longitude=180))
model = W21_10pc.TS.mean(dim='time').plot.contourf(cmap = 'inferno', add_labels = False, add_colorbar = False, levels = levels, vmax = 280,vmin=170,transform=ccrs.PlateCarree())
plt.title('W21 10% PAL', y = yplace, fontsize = Fontsize, color = 'k')

ax = plt.subplot(gs[3,1], projection=ccrs.Robinson(central_longitude=180))
model = W21_1pc.TS.mean(dim='time').plot.contourf(cmap = 'inferno', add_labels = False, add_colorbar = False, levels = levels, vmax = 280,vmin=170,transform=ccrs.PlateCarree())
plt.title('W21 1% PAL', y = yplace, fontsize = Fontsize, color = 'k')

ax = plt.subplot(gs[4,1], projection=ccrs.Robinson(central_longitude=180))
model = W21_01pc.TS.mean(dim='time').plot.contourf(cmap = 'inferno', add_labels = False, add_colorbar = False, levels = levels, vmax = 280,vmin=170,transform=ccrs.PlateCarree())
plt.title('W21 0.1% PAL', y = yplace, fontsize = Fontsize, color = 'k')

#cax = Figure.add_axes([0.125,0.08, 0.775, 0.02]) #start left, vertical start, width, height
cax = Figure.add_axes([0.925,0.125, 0.03, 0.755]) #start left, vertical start, width, height
Cbar = plt.colorbar(model, orientation='vertical', cax=cax)
Cbar.set_label('Temperature [K]', size=15)
Cbar.ax.tick_params(labelsize=15) 
plt.savefig('/Users/gregcooke/python_output/Surface_temperatures.png', bbox_inches = 'tight', dpi = 400)

#%% Transmission spectra final

dpi = 400


sp_path = '/Users/pygjc/TL_UV_paper/spectra_output/'

P19_TrSp = np.genfromtxt(sp_path+'TP-1_e_SSPO_Ideal_0.883hr_12.47pc_rkm-0309-01-01-00000.txt')
P19_10pc_TrSp = np.genfromtxt(sp_path+'TP-1_e_SSPO_10pc_o2_Ideal_0.883hr_12.47pc_rkm-0309-01-01-00000.txt')
P19_1pc_TrSp = np.genfromtxt(sp_path+'TP-1_e_SSPO_1pc_o2_Ideal_0.883hr_12.47pc_rkm-0265-01-01-00000.txt')
P19_01pc_TrSp = np.genfromtxt(sp_path+'TP-1_e_SSPO_0_1pc_o2_Ideal_0.883hr_12.47pc_rkm-0274-01-01-00000.txt')
P19_noTL_TrSp = np.genfromtxt(sp_path+'TP-1_e_SSPO_1day_Ideal_0.883hr_12.47pc_rkm-0228-01-01-00000.txt')

W21_TrSp = np.genfromtxt(sp_path+'TP-1_e_SSPO_MUSC_Ideal_0.883hr_12.47pc_rkm-0294-01-01-00000.txt')
W21_10pc_TrSp = np.genfromtxt(sp_path+'TP-1_e_MUSC_10pc_Ideal_0.883hr_12.47pc_rkm-0296-01-01-00000.txt')
W21_1pc_TrSp = np.genfromtxt(sp_path+'TP-1_e_MUSC_1pc_Ideal_0.883hr_12.47pc_rkm-0270-01-01-00000.txt')
W21_01pc_TrSp = np.genfromtxt(sp_path+'TP-1_e_MUSC_01pc_Ideal_0.883hr_12.47pc_rkm-0277-01-01-00000.txt')
W21_noTL_TrSp = np.genfromtxt(sp_path+'TP-1_e_MUSC_noTL_Ideal_0.883hr_12.47pc_rkm-0256-01-01-00000.txt')


alpha = 1
ls = ':'
texty = 72
texty1 = 74
texty2 = 71
fontsize = 12

cmap_test = plt.cm.binary

ylim = 70

fig = plt.figure(facecolor='w', figsize = (10,4.5))

#plt.style.use('dark_background')
gs = gridspec.GridSpec(1,4000)
markercolor = 'grey'
Fontsize = 15

lw = 2

hatch_color = 'grey'
color = 'k'

x1 = 0.2; x2 = 0.8
ax = plt.subplot(gs[0,:1000])
ax.spines['left'].set_color(color)
ax.spines['bottom'].set_color(color)
ax.spines['top'].set_color(color)
ax.spines['right'].set_color(color)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_linewidth(0)

plt.xlim(0.2,1)

x = np.linspace(0,100)
y = np.linspace(0,100)
# Create a vertical gradient between x=9 and x=10

gradient_image(ax, direction=0, extent=((0.2-0.2)/0.8, (0.34-0.2)/0.8, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)
gradient_image(ax, direction=0, extent=((0.41-0.2)/0.8, (1-0.2)/0.8, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)
gradient_image(ax, direction=0, extent=((0.68-0.2)/0.8, (0.7-0.2)/0.8, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)
gradient_image(ax, direction=0, extent=((0.751-0.2)/0.8, (0.771-0.2)/0.8, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)

plt.text(0.25,texty,'O'+r'$_3$', color = markercolor, fontsize = fontsize, horizontalalignment = 'center')
plt.text(0.55,texty,'O'+r'$_3$', color = markercolor, fontsize = fontsize, horizontalalignment = 'center')
plt.text(0.77,texty,'O'+r'$_2$', color = markercolor, fontsize = fontsize, horizontalalignment = 'center')
plt.text(0.69,texty,'O'+r'$_2$', color = markercolor, fontsize = fontsize, horizontalalignment = 'center')

plt.errorbar([0.45],[15],[(5e-6*((0.1192*696340e3)**2)/(2*0.91*6371e3))/1000], lw=lw, capsize = 5, color = 'k')
plt.text(0.52, 9, '5\nppm', fontsize = fontsize, color = 'k', horizontalalignment = 'center')

plt.plot(P19_TrSp[:,0],P19_TrSp[:,5], lw = lw, color = P19_color, label = 'P19 PI', alpha = alpha)
plt.plot(P19_10pc_TrSp[:,0],P19_10pc_TrSp[:,5], lw = lw, color = P19_10pc_color, label = 'P19 10% PAL', alpha = alpha)
plt.plot(P19_1pc_TrSp[:,0],P19_1pc_TrSp[:,5], lw = lw, color = P19_1pc_color, label = 'P19 1% PAL', alpha = alpha)
plt.plot(P19_01pc_TrSp[:,0],P19_01pc_TrSp[:,5], lw = lw, color = P19_01pc_color, label = 'P19 0.1% PAL', alpha = alpha)
plt.plot(W21_TrSp[:,0],W21_TrSp[:,5], lw = lw, color = W21_color, label = 'W21 PI', alpha = alpha)
plt.plot(W21_10pc_TrSp[:,0],W21_10pc_TrSp[:,5], lw = lw, color = W21_10pc_color, label = 'W21 10% PAL', alpha = alpha)
plt.plot(W21_1pc_TrSp[:,0],W21_1pc_TrSp[:,5],lw = lw,  color = W21_1pc_color, label = 'W21 1% PAL', alpha = alpha)
plt.plot(W21_01pc_TrSp[:,0],W21_01pc_TrSp[:,5], lw = lw, color = W21_01pc_color, label = 'W21 0.1% PAL', alpha = alpha)

plt.ylabel('Effective altitude [km]', fontsize = 15, color = color)
plt.ylim(0,ylim)
plt.text(0.3, 1.5, 'HWO', fontsize = fontsize, horizontalalignment = 'left', color = 'k')

plt.xticks(np.arange(0.2, 0.8, 0.2,))
plt.tick_params(labelsize = 15, top = False,  width = 2, length = 4, right = False, colors = color)


plt.axvspan(0.6,1,0,5/70,color = JWST_color1, alpha = 0.2, lw=0)
HabWorlds_color = '#008BFF'
plt.axvspan(0,1,0,5/70,color = HabWorlds_color, alpha = 0.2, lw=0)

ax = plt.subplot(gs[0,1000:])

ax.spines['left'].set_color(color)
ax.spines['bottom'].set_color(color)
ax.spines['right'].set_color(color)
ax.spines['top'].set_color(color)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['left'].set_linewidth(0)

plt.plot(P19_TrSp[:,0],P19_TrSp[:,5], lw = lw, color = P19_color, label = 'P19 PI', alpha = alpha)
plt.plot(P19_10pc_TrSp[:,0],P19_10pc_TrSp[:,5], lw = lw, color = P19_10pc_color, label = 'P19\n10% PAL', alpha = alpha)
plt.plot(P19_1pc_TrSp[:,0],P19_1pc_TrSp[:,5], lw = lw, color = P19_1pc_color, label = 'P19\n1% PAL', alpha = alpha)
plt.plot(P19_01pc_TrSp[:,0],P19_01pc_TrSp[:,5], lw = lw, color = P19_01pc_color, label = 'P19\n0.1% PAL', alpha = alpha)
plt.plot(W21_TrSp[:,0],W21_TrSp[:,5], lw = lw, color = W21_color, label = 'W21 PI', alpha = alpha)
plt.plot(W21_10pc_TrSp[:,0],W21_10pc_TrSp[:,5], lw = lw, color = W21_10pc_color, label = 'W21\n10% PAL', alpha = alpha)
plt.plot(W21_1pc_TrSp[:,0],W21_1pc_TrSp[:,5],lw = lw,  color = W21_1pc_color, label = 'W21\n1% PAL', alpha = alpha)
plt.plot(W21_01pc_TrSp[:,0],W21_01pc_TrSp[:,5], lw = lw, color = W21_01pc_color, label = 'W21\n0.1% PAL', alpha = alpha)

legend = plt.legend(loc = (1.15,0.05), handlelength = 0.5, ncol = 1, fontsize = 12, frameon = False)
text_colors = [P19_color,P19_10pc_color,P19_1pc_color,P19_01pc_color,
               W21_color,W21_10pc_color,W21_1pc_color,W21_01pc_color]  # Define text colors for each label

for text, color_ in zip(legend.get_texts(), text_colors):
    text.set_color(color_)
    
x1 = 1; x2 = 10

plt.ylim(0,ylim)
plt.xlim(1,11)

#
x1 = 1; x2 = 10

plt.xlim(1,11)
plt.xticks(np.arange(1,12))

plt.text(7.3, 1.5, 'JWST MIRI', fontsize = fontsize, color = 'k')

plt.text(2.8, 1.5, 'JWST NIRSpec', fontsize = fontsize, color = 'k')

# Plot the gradients using imshow
gradient_image(ax, direction=0, extent=((8.5-1)/10, (10.35-1)/10, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)
gradient_image(ax, direction=0, extent=((7.4-1)/10, (8-1)/10, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)
gradient_image(ax, direction=0, extent=((5.2-1)/10, (7.1-1)/10, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)
gradient_image(ax, direction=0, extent=((6-1)/10, (7-1)/10, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)
gradient_image(ax, direction=0, extent=((1.25-1)/10, (1.3-1)/10, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)
gradient_image(ax, direction=0, extent=((1.95-1)/10, (2.08-1)/10, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)
gradient_image(ax, direction=0, extent=((2.6-1)/10, (2.85-1)/10, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)
gradient_image(ax, direction=0, extent=((3.2-1)/10, (3.4-1)/10, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)
gradient_image(ax, direction=0, extent=((4.1-1)/10, (4.5-1)/10, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)
gradient_image(ax, direction=0, extent=((4.7-1)/10, (4.9-1)/10, 0.0, 1), transform=ax.transAxes,
               cmap=plt.cm.binary, cmap_range=(0, 1), alpha=0.3)



plt.text(5.5,texty,'H'+r'$_2$'+'O', color = 'grey', fontsize = fontsize, horizontalalignment = 'center')
plt.text(6.4,texty,'O'+r'$_2$'+'-X', color = 'grey', fontsize = fontsize, horizontalalignment = 'center')
plt.text(7.7,texty2,'CH'+sub(4), color = 'grey', fontsize = fontsize, horizontalalignment = 'center')
plt.text(7.7,texty1+0.5,'N'+sub(2)+'O', color = 'grey', fontsize = fontsize, horizontalalignment = 'center')
plt.text(9.3,texty,'O'+sub(3), color = 'grey', fontsize = fontsize, horizontalalignment = 'center')
plt.text(6.4,texty,'O'+r'$_2$'+'-X', color = 'grey', fontsize = fontsize, horizontalalignment = 'center')
plt.text(1.27,texty,'O'+r'$_2$', color = 'grey', fontsize = fontsize, horizontalalignment = 'center')
plt.text(2,texty,'CO'+r'$_2$', color = markercolor, fontsize = fontsize, horizontalalignment = 'center')
plt.text(2.7,texty,'CO'+r'$_2$', color = markercolor, fontsize = fontsize, horizontalalignment = 'center')
plt.text(3.3,texty1+0.5,'CH'+r'$_4$', color = markercolor, fontsize = fontsize, horizontalalignment = 'center')
plt.text(3.3,texty2,'O'+r'$_3$', color = markercolor, fontsize = fontsize, horizontalalignment = 'center')
plt.text(4.25,texty,'CO'+r'$_2$', color = markercolor, fontsize = fontsize, horizontalalignment = 'center')
plt.text(4.8,texty1+0.5,'CO'+r'$_2$', color = markercolor, fontsize = fontsize, horizontalalignment = 'center')
plt.text(4.8,texty2,'O'+r'$_3$', color = markercolor, fontsize = fontsize, horizontalalignment = 'center')

plt.xlabel('Wavelength ['+r'$\rm \mu$'+'m]', fontsize = 15, x =  0.33, color = color)
plt.tick_params(which = 'both',labelsize = 15,  top = False, labelleft = False, left = False, right = False, width = 2, length = 4, colors = color)
ax.set_yticks([])

ax2 = ax.twinx()

ax2.spines['left'].set_color(color)
ax2.spines['bottom'].set_color(color)
ax2.spines['right'].set_color(color)
ax2.spines['top'].set_color(color)
ax2.spines['left'].set_linewidth(0)
ax2.spines['right'].set_visible(True)
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)

#plt.errorbar([0.52],[4900],[5], lw=lw, capsize = 5, color = 'k', zorder = 1)

#ax2.plot(P19_TrSp[:,0], 1e6*((P19_TrSp[:,5]+(6371*0.91))/(695500*0.1192))**2, ls = '-',lw = 2, color = P19_color, alpha = alpha)
plt.ylim(1e6*((6371*0.91)/(695500*0.1192))**2,1e6*(((6371*0.91)+70)/(695500*0.1192))**2)

plt.ylabel('Transit depth [ppm]', fontsize = 15, color = color)
plt.tick_params(which = 'both', labelsize = 15,  top = False, labelleft = False, width = 2, length = 4, left = False, right = True, colors = color)

plt.axvspan(1,5,0,5/70,color = JWST_color1, alpha = 0.2, lw=0)
plt.axvspan(5,11,0,5/70,color = JWST_color2, alpha = 0.2, lw=0)
HabWorlds_color = '#008BFF'
plt.axvspan(1,2.5,0,5/70,color = HabWorlds_color, alpha = 0.2, lw=0)

plt.savefig('/Users/gregcooke/python_output/TP1_e_transmission_spectra.png', facecolor='white', dpi = dpi, bbox_inches = 'tight', transparent = True)
plt.savefig('/Users/gregcooke/python_output/TP1_e_transmission_spectra.svg', format='svg', dpi = dpi, bbox_inches = 'tight', transparent = True)

#%%O2 vs O3 column curve

color = 'k'
dpi =300


plt.figure(figsize = (12,4*1.2))

# Data
x = [1, 2, 3, 4, 5]
y = [2, 4, 6, 8, 10]
positive_error = [0.5, 0.8, 0.2, 1.0, 0.4]
negative_error = [0.3, 0.6, 0.4, 0.7, 0.2]

# Plot
lw = 3
capsize = 6
markersize = 7

plt.errorbar([1.5], [LWAV(O3_col_One50)], yerr=[[LWAV(O3_col_One50)-O3_col_One50.min()],[O3_col_One50.max() - LWAV(O3_col_One50)]], fmt='', ls = '', marker = 'o', markersize = markersize, capsize = capsize, lw = lw, color = color, label = 'Earth')
plt.errorbar([1], [LWAV(O3_col_Pre)], yerr=[[LWAV(O3_col_Pre)-O3_col_Pre.min()],[O3_col_Pre.max() - LWAV(O3_col_Pre)]], fmt='-',  marker = 'o', markersize = markersize, capsize = capsize, lw = lw, color = color)
plt.errorbar([0.5], [LWAV(O3_col_Fifty)], yerr=[[LWAV(O3_col_Fifty)-O3_col_Fifty.min()],[O3_col_Fifty.max() - LWAV(O3_col_Fifty)]], fmt='-',  marker = 'o', markersize = markersize, capsize = capsize, lw = lw, color = color)
plt.errorbar([0.1], [LWAV(O3_col_Ten)], yerr=[[LWAV(O3_col_Ten)-O3_col_Ten.min()],[O3_col_Ten.max() - LWAV(O3_col_Ten)]], fmt='-',  marker = 'o', markersize = markersize, capsize = capsize, lw = lw, color = color)
plt.errorbar([5/100], [LWAV(O3_col_Five)], yerr=[[LWAV(O3_col_Five)-O3_col_Five.min()],[O3_col_Five.max() - LWAV(O3_col_Five)]], fmt='-',  marker = 'o', markersize = markersize, capsize = capsize, lw = lw, color = color)
plt.errorbar([1/100], [LWAV(O3_col_One)], yerr=[[LWAV(O3_col_One)-O3_col_One.min()],[O3_col_One.max() - LWAV(O3_col_One)]], fmt='-',  marker = 'o', markersize = markersize, capsize = capsize, lw = lw, color = color)
plt.errorbar([0.5/100], [LWAV(O3_col_Zero5)], yerr=[[LWAV(O3_col_Zero5)-O3_col_Zero5.min()],[O3_col_Zero5.max() - LWAV(O3_col_Zero5)]], fmt='-',  marker = 'o', markersize = markersize, capsize = capsize, lw = lw, color = color)
plt.errorbar([0.1/100], [LWAV(O3_col_Zero1)], yerr=[[LWAV(O3_col_Zero1)-O3_col_Zero1.min()],[O3_col_Zero1.max() - LWAV(O3_col_Zero1)]], fmt='-',  marker = 'o', markersize = markersize, capsize = capsize, lw = lw, color = color)

plt.errorbar([1], [LWAV(P19_O3_col)], yerr=[[LWAV(P19_O3_col)-P19_O3_col.min()],[P19_O3_col.max() - LWAV(P19_O3_col)]], fmt='-', capsize = capsize, lw = lw,ls = '',  marker = 'o', markersize = markersize, color = P19_color, label = 'TRAPPIST-1 e P19 (high UV)')
plt.errorbar([1], [LWAV(W21_O3_col)], yerr=[[LWAV(W21_O3_col)-W21_O3_col.min()],[W21_O3_col.max() - LWAV(W21_O3_col)]], fmt='-', capsize = capsize,ls = '',  lw = lw, marker = 'o', markersize = markersize, color = W21_color, label = 'TRAPPIST-1 e W21 (low UV)')
plt.errorbar([10/100], [LWAV(W21_10pc_O3_col)], yerr=[[LWAV(W21_10pc_O3_col)-W21_10pc_O3_col.min()],[W21_10pc_O3_col.max() - LWAV(W21_10pc_O3_col)]], fmt='-', capsize = capsize, lw = lw,  marker = 'o', markersize = markersize,color = W21_color)
plt.errorbar([1/100], [LWAV(W21_1pc_O3_col)], yerr=[[LWAV(W21_1pc_O3_col)-W21_1pc_O3_col.min()],[W21_1pc_O3_col.max() - LWAV(W21_1pc_O3_col)]], fmt='-', capsize = capsize, lw = lw, marker = 'o', markersize = markersize,color = W21_color)
plt.errorbar([1/1000], [LWAV(W21_01pc_O3_col)], yerr=[[LWAV(W21_01pc_O3_col)-W21_01pc_O3_col.min()],[W21_01pc_O3_col.max() - LWAV(W21_01pc_O3_col)]], fmt='-', capsize = capsize, lw = lw,  marker = 'o', markersize = markersize,color = W21_color)
plt.errorbar([10/100], [LWAV(P19_10pc_O3_col)], yerr=[[LWAV(P19_10pc_O3_col)-P19_10pc_O3_col.min()],[P19_10pc_O3_col.max() - LWAV(P19_10pc_O3_col)]], fmt='-', capsize = capsize, lw = lw,  marker = 'o', markersize = markersize,color = P19_color)
plt.errorbar([1/100], [LWAV(P19_1pc_O3_col)], yerr=[[LWAV(P19_1pc_O3_col)-P19_1pc_O3_col.min()],[P19_1pc_O3_col.max() - LWAV(P19_1pc_O3_col)]], fmt='-', capsize = capsize, lw = lw,  marker = 'o', markersize = markersize,color = P19_color)
plt.errorbar([1/1000], [LWAV(P19_01pc_O3_col)], yerr=[[LWAV(P19_01pc_O3_col)-P19_01pc_O3_col.min()],[P19_01pc_O3_col.max() - LWAV(P19_01pc_O3_col)]], fmt='-', capsize = capsize, lw = lw,  marker = 'o', markersize = markersize,color = P19_color)

#plt.legend(loc = 0, fontsize = 15, frameon = False)
plt.text(3.2e-3, 200, 'TRAPPIST-1 e P19\n(high UV)', color = P19_color, fontsize = 15, horizontalalignment = 'center')
plt.text(3.2e-3, 900, 'TRAPPIST-1 e W21\n(low UV)', color = W21_color, fontsize = 15, horizontalalignment = 'center')
plt.text(3.2e-3, 40, 'Earth', color = 'k', fontsize = 15, horizontalalignment = 'center')

plt.ylabel('O'+sub(3)+' column [DU]', fontsize = 15, color = color)
plt.xlabel('O'+sub(2)+' concentration [PAL]', fontsize = 15, color = color)
plt.tick_params(labelsize = 15,  top = False, labelleft = True, width = 2, length = 4, left = True, right = False, colors = color)
plt.tick_params(which = 'minor', colors = color)

plt.yscale('log')
plt.ylim(10,10000)

plt.xscale('log')
plt.savefig('/Users/gregcooke/python_output/O2_O3_column_curve.png', dpi = dpi, bbox_inches = 'tight')
plt.savefig('/Users/gregcooke/python_output/O2_O3_column_curve.eps', format='eps', dpi = dpi, bbox_inches = 'tight')

#%% Figures for appendix

Surface_rob_plot_gen_col(cmap = 'Blues', mmw = 18, var = 'H2O',
                             clabel = 'H'+ sub(2) +'O column [kg m'+sup(-2)+']',
                             levels = 11, save_name = 'H2O_col', 
                             log = True, vmax = 10, vmin = 0.1)

Surface_rob_plot_gen_col(cmap = 'Blues_r', mmw = 18, var = 'CLDLIQ',
                             clabel = 'H'+ sub(2) +'O cloud liquid column [kg m'+sup(-2)+']',
                             levels = 11, save_name = 'CLDLIQ_col', 
                             log = True, vmax = 1, vmin = 0.01)

Surface_rob_plot_gen_col(cmap = 'Blues_r', mmw = 18, var = 'CLDICE',
                             clabel = 'H'+ sub(2) +'O cloud ice column [kg m'+sup(-2)+']',
                             levels = 11, save_name = 'CLDLIQ_col', 
                             log = True, vmax = 0.005, vmin = 0.0001)


#%%
Fontsize = 15
save_name = 'High_cloud_fraction'
clabel = 'High cloud fraction'
log = False
levels = 11
vmax = 1
vmin = 0

Figure = plt.figure(figsize = (12,7))
#clabel = 'Zonal winds [ms'+sup(-1)+']'
cmap = 'Blues_r'
gs = gridspec.GridSpec(2,2)
gs.update(hspace = 0.12)
gs.update(wspace = 0.05)
    
ax = plt.subplot(gs[0,0], projection=ccrs.Robinson(central_longitude=180))
ax.coastlines('110m')

P19_col = P19.CLDHGH

if(log == True):
    model = P19_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, norm = colors.LogNorm(vmax = vmax, vmin = vmin), transform=ccrs.PlateCarree())
else:
    model = P19_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, transform=ccrs.PlateCarree(), vmax = vmax, vmin =  vmin)
plt.tick_params(labelsize = Fontsize, length = 4, width = 2, right = True)
plt.title('P19 PI', fontsize = 15)

W21_col = W21.CLDHGH

ax = plt.subplot(gs[0,1], projection=ccrs.Robinson(central_longitude=180))
ax.coastlines('110m')
if(log == True):
    model = W21_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, norm = colors.LogNorm(vmax = vmax, vmin = vmin), transform=ccrs.PlateCarree())
else:
    model = W21_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, transform=ccrs.PlateCarree(), vmax = vmax, vmin =  vmin)
plt.tick_params(labelsize = Fontsize, length = 4, width = 2, right = True, labelleft = False)
plt.title('W21 PI', fontsize = 15)

P19_01pc_col = P19_01pc.CLDHGH

ax = plt.subplot(gs[1,0], projection=ccrs.Robinson(central_longitude=180))
ax.coastlines('110m')

if(log == True):
    model = P19_01pc_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, norm = colors.LogNorm(vmax = vmax, vmin = vmin), transform=ccrs.PlateCarree())
else:
    model = P19_01pc_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, transform=ccrs.PlateCarree(), vmax = vmax, vmin =  vmin)
plt.tick_params(labelsize = Fontsize, length = 4, width = 2, right = True)
plt.title('P19 0.1% PAL', fontsize = 15)

W21_01pc_col = W21_01pc.CLDHGH

ax = plt.subplot(gs[1,1], projection=ccrs.Robinson(central_longitude=180))
ax.coastlines('110m')
if(log == True):
    model = W21_01pc_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, norm = colors.LogNorm(vmax = vmax, vmin = vmin), transform=ccrs.PlateCarree())
else:
    model = W21_01pc_col.plot(cmap = cmap, levels = levels, add_labels = False, add_colorbar = False, transform=ccrs.PlateCarree(), vmax = vmax, vmin =  vmin)
plt.tick_params(labelsize = Fontsize, length = 4, width = 2, right = True, labelleft = False)
plt.title('W21 0.1% PAL', fontsize = 15)

cax = Figure.add_axes([0.92, 0.135, 0.025, 0.74])#start left, vertical start, width, height
Cbar = plt.colorbar(model, orientation='vertical', cax=cax)
Cbar.set_label(clabel, size=15)
Cbar.ax.tick_params(labelsize=15) 

plt.savefig('/Users/gregcooke/python_output/'+save_name + '.png', dpi = 200, bbox_inches = 'tight')
plt.savefig('/Users/gregcooke/python_output/'+save_name + '.eps', format = 'eps', dpi = 200, bbox_inches = 'tight')

#%% Ox production

cmap  ='inferno'
vmax = 1e12; vmin = 1e8

levels = 21

Figure = plt.figure(figsize = (12,9))
gs = gridspec.GridSpec(2,2)

plt.subplot(gs[0,0])
DS = P19; plt.title('P19 PI', fontsize = 15)
O2_prod_SSPO = n_dens(DS,'O2')*(DS.jo2_a+DS.jo2_b)
model = (2*O2_prod_SSPO).mean(dim = 'lon').plot(cmap = cmap, norm = colors.LogNorm(vmax = vmax, vmin = vmin), levels = levels, add_labels = False, add_colorbar = False)
plt.ylim(1000,1e-4); plt.yscale('log')
plt.tick_params(labelsize = 15, left = True, right = True, labelleft = True,labelbottom = False,  width = 2, length = 4)
plt.ylabel('Pressure [hPa]', fontsize = 15)
plt.xticks(np.arange(-90, 90+30, 30))
print(((2*O2_prod_SSPO).mean(dim = 'lon')).max())
 
plt.subplot(gs[1,0])
DS = P19_01pc; plt.title('P19 0.1% PAL', fontsize = 15)
O2_prod_SSPO = n_dens(DS,'O2')*(DS.jo2_a+DS.jo2_b)
model = (2*O2_prod_SSPO).mean(dim = 'lon').plot(cmap = cmap, norm = colors.LogNorm(vmax = vmax, vmin = vmin), levels = levels, add_labels = False, add_colorbar = False)
plt.ylim(1000,1e-4); plt.yscale('log')
plt.tick_params(labelsize = 15, left = True, right = True, labelleft = True, width = 2, length = 4)
plt.ylabel('Pressure [hPa]', fontsize = 15)
plt.xlabel('Latitude ['+sup('\circ')+']', fontsize = 15); plt.xticks(np.arange(-90, 90+30, 30))
print(((2*O2_prod_SSPO).mean(dim = 'lon')).max())

plt.subplot(gs[0,1])
DS = W21; plt.title('W21 PI', fontsize = 15)
O2_prod_SSPO = n_dens(DS,'O2')*(DS.jo2_a+DS.jo2_b)
model = (2*O2_prod_SSPO).mean(dim = 'lon').plot(cmap = cmap, norm = colors.LogNorm(vmax = vmax, vmin = vmin), levels = levels, add_labels = False, add_colorbar = False)
plt.ylim(1000,1e-4); plt.yscale('log')
plt.tick_params(labelsize = 15, left = True, right = True, labelleft = False,labelbottom = False,  width = 2, length = 4)
plt.xticks(np.arange(-90, 90+30, 30))
print(((2*O2_prod_SSPO).mean(dim = 'lon')).max())

plt.subplot(gs[1,1])
DS = W21_01pc; plt.title('W21 0.1% PAL', fontsize = 15)
O2_prod_SSPO = n_dens(DS,'O2')*(DS.jo2_a+DS.jo2_b)
model = (2*O2_prod_SSPO).mean(dim = 'lon').plot(cmap = cmap, norm = colors.LogNorm(vmax = vmax, vmin = vmin), levels = levels, add_labels = False, add_colorbar = False)
plt.ylim(1000,1e-4); plt.yscale('log')
plt.tick_params(labelsize = 15, left = True, right = True, labelleft = False, width = 2, length = 4)
plt.xlabel('Latitude ['+sup('\circ')+']', fontsize = 15); plt.xticks(np.arange(-90, 90+30, 30))
print(((2*O2_prod_SSPO).mean(dim = 'lon')).max())

cax = Figure.add_axes([0.95,0.125, 0.03, 0.755]) #start left, vertical start, width, height
Cbar = plt.colorbar(model, orientation='vertical', cax=cax, extend = 'neither')
Cbar.set_label('O production from O'+sub(2)+' photolysis [molecules m'+sup(-3)+' s'+sup(-1)+']', size=15)
Cbar.ax.tick_params(labelsize=15) 
plt.tick_params(labelsize = 15, left = True, right = True, labelleft = False, width = 2, length = 4)

plt.savefig('/Users/gregcooke/python_output/Ox_prod.png', dpi = 400, bbox_inches='tight')
  
#%% Ozone number density zonal mean plot

levels = 21

vmax = 1e20; vmin = 1e16
ylim = 1e-2

Figure = plt.figure(figsize = (12,9))
gs = gridspec.GridSpec(2,2)

plt.subplot(gs[0,0])
DS = P19; plt.title('P19 PI', fontsize = 15)
model = n_dens(DS,'O3').mean(dim = 'lon').plot(cmap = cmap, norm = colors.LogNorm(vmax = vmax, vmin = vmin), levels = levels, add_labels = False, add_colorbar = False)
print((n_dens(DS,'O3').mean(dim = 'lon')).max())
plt.ylim(1000,ylim); plt.yscale('log')
plt.tick_params(labelsize = 15, left = True, right = True, labelleft = True, labelbottom = False, width = 2, length = 4)
plt.ylabel('Pressure [hPa]', fontsize = 15)
plt.xticks(np.arange(-90, 90+30, 30))

plt.subplot(gs[1,0])
DS = P19_01pc; plt.title('P19 0.1% PAL', fontsize = 15)
model = model = n_dens(DS,'O3').mean(dim = 'lon').plot(cmap = cmap, norm = colors.LogNorm(vmax = vmax, vmin = vmin), levels = levels, add_labels = False, add_colorbar = False)
print((n_dens(DS,'O3').mean(dim = 'lon')).max())
plt.ylim(1000,ylim); plt.yscale('log')
plt.tick_params(labelsize = 15, left = True, right = True, labelleft =  True, width = 2, length = 4)
plt.ylabel('Pressure [hPa]', fontsize = 15)
plt.xlabel('Latitude ['+sup('\circ')+']', fontsize = 15); plt.xticks(np.arange(-90, 90+30, 30))

plt.subplot(gs[0,1])
DS = W21; plt.title('W21 PI', fontsize = 15)
model = n_dens(DS,'O3').mean(dim = 'lon').plot(cmap = cmap, norm = colors.LogNorm(vmax = vmax, vmin = vmin), levels = levels, add_labels = False, add_colorbar = False)
print((n_dens(DS,'O3').mean(dim = 'lon')).max())
plt.ylim(1000,ylim); plt.yscale('log')
plt.tick_params(labelsize = 15, left = True, right = True, labelleft = False, labelbottom = False, width = 2, length = 4)
plt.xticks(np.arange(-90, 90+30, 30))

plt.subplot(gs[1,1])
DS = W21_01pc; plt.title('W21 0.1% PAL', fontsize = 15)
model = n_dens(DS,'O3').mean(dim = ['lon','time']).plot(cmap = cmap, norm = colors.LogNorm(vmax = vmax, vmin = vmin), levels = levels, add_labels = False, add_colorbar = False)
print((n_dens(DS,'O3').mean(dim = 'lon')).max())
plt.ylim(1000,ylim); plt.yscale('log')
plt.tick_params(labelsize = 15, left = True, right = True, labelleft = False, width = 2, length = 4)
plt.xlabel('Latitude ['+sup('\circ')+']', fontsize = 15); plt.xticks(np.arange(-90, 90+30, 30))

cax = Figure.add_axes([0.95,0.125, 0.03, 0.755]) #start left, vertical start, width, height
Cbar = plt.colorbar(model, orientation='vertical', cax=cax)
Cbar.set_label('O'+sub(3)+' number density [molecules m'+sup(-3)+']', size=15)
Cbar.ax.tick_params(labelsize=15) 
plt.tick_params(labelsize = 15, left = True, right = True, labelleft = False, width = 2, length = 4)

plt.savefig('/Users/gregcooke/python_output/O3_dens.png', dpi = 400, bbox_inches='tight')
#%% Show TRAPPIST-1 stellar spectra
  
solar_dir = '/Users/gregcooke/stellar_files/'
solar_file = 'TRAPPIST1_flux_at_e_W21.nc'#read in file
ds = xr.open_dataset(solar_dir+solar_file)#attahc file to dataset
ssi = ds['ssi'].isel(time=0) #define dataset from file

W21_flux = ssi.values
W21_wavelength = ssi.wavelength.values

solar_dir = '/Users/gregcooke/stellar_files/'
solar_file = 'TRAPPIST1_flux_at_e_P19.nc'#read in file
ds = xr.open_dataset(solar_dir+solar_file)#attahc file to dataset
ssi = ds['ssi'].isel(time=0) #define dataset from file

P19_flux = ssi.values
P19_wavelength = ssi.wavelength.values

plt.figure(figsize = (11*1.1,4.5*1.1))
ax = plt.subplot()
markercolor = 'whitesmoke'
plt.plot(P19_wavelength, P19_flux/1000, lw = 2, color = P19_color, label = 'TRAPPIST-1 at e (P19)')
plt.plot(P19_wavelength, W21_flux/1000, lw = 2, color = W21_color, label = 'TRAPPIST-1 at e (W21)')
plt.ylabel('Top of atmosphere irradiance\nper unit wavelength [W m'+r'$^{-2}$'+' nm'+r'$^{-1}$'+']',fontsize=15,color='k')
plt.xscale('log')
plt.yscale('log')
plt.tick_params(labelsize = 15, length = 4, width = 2, right = True, colors = 'k')#, labelbottom = False) #change the size of the ticks and labels  
plt.tick_params(which = 'minor', labelsize = 15, length = 2, width = 1, right = True, colors = 'k')#, labelbottom = False) #change the size of the ticks and labels  
plt.xlim(100,5e3)#WACCM_wavelength.max())

plt.legend(frameon = False,  labelcolor='linecolor', loc = (0.6,0.2), fontsize = 15)

plt.axvspan(100, 400, alpha = 0.25, facecolor = 'slategrey', edgecolor = 'slategrey', lw = 0)
plt.ylim(1e-6,1)
plt.yticks([1e-6, 1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1])
plt.text(200,2e-5,'UV',color='#434343',fontsize=20, horizontalalignment = 'center')
plt.text(121,0.14,'L-'+r'$\alpha$',color='#434343',fontsize=15, horizontalalignment = 'center')

xtick_positions = [100, 200, 400, 600, 1000, 2000, 4000]
xtick_labels = [100, 200, 400, 600, 1000, 2000, 4000]
plt.xticks(xtick_positions, xtick_labels)

plt.xlabel('Wavelength [nm]', fontsize = 15,color='k')
plt.savefig('/Users/gregcooke/python_output/TP_1_spectra.png', dpi = 200, bbox_inches = 'tight')#, transparent = True)

#%% Check flux in each band

from scipy.integrate import trapz

W21_TSI = trapz(W21_flux, W21_wavelength)/1000
P19_TSI = trapz(P19_flux, P19_wavelength)/1000

# Schumann Runge continum
W21_HC = trapz(W21_flux[189:189+40], W21_wavelength[189:189+40])/1000
P19_HC = trapz(P19_flux[189:189+40], P19_wavelength[189:189+40])/1000

# Schumann Runge bands
W21_HC = trapz(W21_flux[189:189+40], W21_wavelength[189:189+40])/1000
P19_HC = trapz(P19_flux[189:189+40], P19_wavelength[189:189+40])/1000

# Herzberg continuum
W21_HC = trapz(W21_flux[189:189+40], W21_wavelength[189:189+40])/1000
P19_HC = trapz(P19_flux[189:189+40], P19_wavelength[189:189+40])/1000

# Hartley bands
W21_HC = trapz(W21_flux[189:189+40], W21_wavelength[189:189+40])/1000
P19_HC = trapz(P19_flux[189:189+40], P19_wavelength[189:189+40])/1000

# Huggins bands
W21_HC = trapz(W21_flux[189:189+40], W21_wavelength[189:189+40])/1000
P19_HC = trapz(P19_flux[189:189+40], P19_wavelength[189:189+40])/1000

# Chappuis bands
W21_HC = trapz(W21_flux[189:189+40], W21_wavelength[189:189+40])/1000
P19_HC = trapz(P19_flux[189:189+40], P19_wavelength[189:189+40])/1000


