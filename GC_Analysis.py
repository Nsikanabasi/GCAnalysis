#!/apps/python3/3.5.2/bin/python3
# -*- coding: utf-8 -*-
################################################################################
#########################  GEOS-Chem Data Analysis Codes #######################
################################################################################

#Import all the required libraries and applications
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt #for plotting a georeferenced data
from mpl_toolkits.basemap import Basemap #for plotting a georeferenced data
# importing the needed file for analysis
my_example_nc_file = '/home/574/nsu574/geosfp_4x5_tropchem/coards.20130701.nc'
fh = Dataset(my_example_nc_file,mode='r')
#Reading in the variables from the .nc file earlier imported.
lons = fh.variables['lon'][:]
lats = fh.variables['lat'][:]
time = fh.variables['time'][:]
tmax_units = fh.variables['time'].units
#closing the files after Reading
fh.close()
# Get some parameters for the Stereographic Projection
m = Basemap(width=5000000,height=3500000,
            resolution='l',projection='stere',\
            lat_ts=40,lat_0=lat_0,lon_0=lon_0)

# If the lon and lat variables are 1D,use meshgrid to create 2D arrays.
# This is not necessary if coordinates are already in 2D arrays.
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)
"""
Put todays scripts from the tutorial here and write out the details
"""
"""
################################################################################
########################## PLOTTING DATA OVER A MAP#############################
################################################################################
cs = m.pcolor(xi,yi,np.squeeze(tmax))

# Add Grid Lines
m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# Add Colorbar
cbar = m.colorbar(cs, location='bottom', pad="10%")
cbar.set_label(tmax_units)

# Add Title
plt.title('DJF Maximum Temperature') #Put the appropriate title of the plot here

plt.show() #shows the plot
"""
"""
from netCDF4 import MFDataset



from mpl_toolkits.basemap import Basemap,cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from datetime import datetime,timedelta
import math

import sys
from pyhdf.SD import SD, SDC
################################################################################
##### FILL IN SETTINGS  ########################################################
################################################################################
Year = '2013' # The year of interest e.g. 2010, 2011, 2012, ......
Tracer = 'o3'  # The species of interest e.g. co, co2, c2h6, hcho, etc
################################################################################
##### SETTING UP OTHER VARIABLES  ##############################################
################################################################################
g = 9.8               # gravitational force (m/s)
NAv = 6.0221415e+23		   # Avogadro number
MWair = 28.9644          # Molar mass of air (g/mole)
psurf = 1013.25		#Atmospheric pressure at the surface (mb)
H = (8.314*240)/(28.9644E-03*9.8)
pcol_const = (NAv* 10)/(MWair*g)   # scaling factor for turning vertical...
#...mixing ratio (vmr) into pcol (using hPa levels)
################################################################################
##### STATIONS OF INTEREST   ###################################################
################################################################################
#2x25		       Real					IDL/Python
#Station           lon      lat          I    	 J
#Wollongong      : 150.88   -34.41        132   28
#Darwin          : 130.89   -12.43        124   39
#Lauder          : 169.68   -45.04        140   22
#Cape Grim 		 : 144.7    -40.7         130   25
#Reunion         : 55.49    -20.90        94    35
#Cape Fergusson	 : 147.06   -19.28        131   35

#Pulls out the Tracer, the Molecular weight of the Tracer (g/mol)...
#... and the gc_field
if Tracer == 'co':
    MWtracer = 28.01 #Molecular weight of the Tracer
	gc_field = 'IJ_AVG_S__CO' # as represented in GEOS-Chem
if Tracer == 'c2h6':
	MWtracer = 30.07
	gc_field = 'IJ_AVG_S__C2H6'
if Tracer == 'h2co':
	MWtracer = 30.031
	gc_field = 'IJ_AVG_S__CH2O'
if Tracer == 'o3'
	MWtracer = 48
	gc_field = 'IJ_AVG_S__O3'
################################################################################
##### GEOS-Chem   ##############################################################
################################################################################
#[time,lev,lat,lon]
print("Processing GEOS-Chem data")
#select GC File(s)
dataset1 = MFDataset('/home/574/nsu574/geosfp_4x5_tropchem/coards.'+Year+'*.nc')

####################################
#Monthly work
####################################
#Extract from Monthly files
gc_times=dataset1.variables['time'][:]	# Extract GC time (TAU format - hours since 1985)
NAIR=dataset1.variables['BXHGHT_S__N_AIR_'][:] # Extract air density (molecules per m3)
BXHT=dataset1.variables['BXHGHT_S__BXHEIGHT'][:] # Extract grid box height (metres)
gc_area=dataset1.variables['DXYP__DXYP'][:]	# Extract grid column surface area
gc_vmr=dataset1.variables[gc_field][:]	# Extract volume molar ratio (ppb Carbon atoms)

#Create long / lat grid
gc_lon = np.arange(-182.5,182.5,5)	# Create grid edges for longitudes (res 2.5 deg)
gc_lat = np.arange(-91,95,4)	# Create grid edges for latitude (res 2 deg)

#convert gc time
d0=datetime(1985,1,1,0)			# Tau 0 = 1/1/1985
gc_datetimes = []	# Create list for new time array
for time in gc_times:
	hrs=timedelta(hours=int(time))	# DateTime package extracts 'hours' since 1985 date format
	gc_datetimes.append(d0+hrs)	# Fill list with formated dates

#Create box centres heights
BXHT = np.mean(BXHT, axis=(0,2,3))	# Remove the time, long and lat dimensions
gc_lvl = []						# Create list for level altitude
gc_lvl.append((BXHT[0])/2)	# Fill first level
for i in range(1,47):
	 gc_lvl.append((gc_lvl[i-1]+BXHT[i-1]/2)+(BXHT[i])/2)	# Fill in all other levels

gc_vmr_surf = gc_vmr[:,0,:,:]			# Extract surface vmr
gc_vmr_surf_avg = np.mean(gc_vmr_surf, axis=(1,2))	# Average to global surface vmr
gc_vmr_surf_map = np.mean(gc_vmr_surf, axis=(0))	# Average to surface vmr / year
gc_vmr_vert_avg = np.mean(gc_vmr, axis=(0,2,3))		# Average to global and yearly vertical profile

NAIR = np.mean(NAIR, axis=(0,2,3))			# Remove the time, long and lat dimensions
gc_molec_vert_avg = gc_vmr_vert_avg*1E-9*NAIR		# Convert to vmr to C/


#Select stations from Monthly files

#surf_stations=[]
#surf_stations.append(gc_datetimes)
#surf_stations.append(gc_vmr_surf[:,28,132])		# Extract surface vmr for wollongong
#surf_stations.append(gc_vmr_surf[:,27,132])		# Extract surface vmr for grid bellow wollongong
#surf_stations.append(gc_vmr_surf[:,39,124])		# Extract surface vmr for darwin
#surf_stations.append(gc_vmr_surf[:,22,140])		# Extract surface vmr for lauder
#surf_stations.append(gc_vmr_surf[:,25,130])		# Extract surface vmr for cape grim
#surf_stations.append(gc_vmr_surf[:,35,131])		# Extract surface vmr for capeferguson
#surf_stations=np.transpose(np.array(surf_stations))
#np.savetxt("Monthly_GC_Stations_Surface_"+Tracer+"_"+Inventory+"_"+Year+"_"+Tracer.upper()+".csv",surf_stations, fmt="%s", delimiter=",", header="Month,Wollongong,Below_Wollongong,Darwin,Lauder,Cape_Grim,Cape_Ferguson")


if False:
	#Find matching FTS/GC time stamps and extrac GC data
	station_matched_vmr = []
	station_matched_airden = []
	for i in np.arange(0,len(nd_datetimes),1):
		for j in np.arange(0,len(station_datetimes),1):
			if nd_datetimes[i] == station_datetimes[j]:
				station_matched_vmr.append(station_vmr[j,:])
				station_matched_airden.append(gc_airden[j,:])
				break
	station_matched_vmr = np.array(station_matched_vmr)
	station_matched_airden = np.array(station_matched_airden)

	station_matched_pcol = []
	for i in np.arange(0,len(nd_datetimes),1):
		station_matched_pcol.append(station_matched_vmr[i,:] / 1e9 * station_matched_airden[i,:] * BXHT[:] * 100)
	station_matched_pcol =np.array(station_matched_pcol)					## !!!!!! ##

	#Interpolate GC data to FTS levels
	rev_nd_lvl = nd_lvl[::-1]
	interp_station_vmr = []
	interp_station_airden = []
	for i in np.arange(0,len(nd_datetimes),1):
		interp_station_vmr.append(np.interp(rev_nd_lvl,gc_lvl,station_matched_vmr[i,:]))
		interp_station_airden.append(np.interp(rev_nd_lvl,gc_lvl,station_matched_airden[i,:]))
	interp_station_vmr=np.array(interp_station_vmr)[:,::-1]
	interp_station_airden=np.array(interp_station_airden)[:,::-1]

	nd_BXHGT=np.array(nd_BXHGT)

	monthly_station_vmr=gc_vmr[:,:,J,I]
	interp_monthly_station_vmr=[]
	for i in np.arange(0,12,1):
		interp_monthly_station_vmr.append(np.interp(rev_nd_lvl,gc_lvl,monthly_station_vmr[i,:]))
	interp_monthly_station_vmr =np.array(interp_monthly_station_vmr)[:,::-1]

	#Calculate interpolated GC partial columns
	station_pcol =[]
	for i in np.arange(0,len(nd_datetimes),1):
		station_pcol.append(interp_station_vmr[i,:] / 1e9 * interp_station_airden[i,:] * nd_BXHGT[:] * 100)
	station_pcol =np.array(station_pcol)						## !!!!!! ##

	monthly_station_pcol =[]
	for i in np.arange(0,12,1):
		monthly_station_pcol.append(2.12e13 * interp_monthly_station_vmr[i,:] * nd_edgePressdif[:])
	monthly_station_pcol =np.array(monthly_station_pcol)

	#Apply averaging kernels and aprioris on GC partial columns
	station_ap_pcol =[]
	for i in np.arange(0,len(nd_datetimes),1):
		station_ap_pcol.append(nd_apr[i,:] + nd_avk[i,:]*(station_pcol[i,:]-nd_apr[i,:]))
	station_ap_pcol =np.array(station_ap_pcol)					## !!!!!! ## and nd_pcol

	monthly_station_ap_pcol =[]
	for i in np.arange(0,12,1):
		monthly_station_ap_pcol.append(avg_apr[i,:] + avg_avk[i,:]*(monthly_station_pcol[i,:]-avg_apr[i,:]))
	monthly_station_ap_pcol =np.array(monthly_station_ap_pcol)

	#Sum partial columns into total columns
	station_tcol=np.sum(station_ap_pcol, axis=(1))
	monthly_station_tcol=np.sum(monthly_station_ap_pcol, axis=(1))

	toprint = []
	toprint.append(nd_datetimes)
	toprint.append(station_tcol)
	toprint.append(nd_xp)
	toprint=np.transpose(np.array(toprint))
	np.savetxt(Station+"_"+Year+"_"+Inventory+"_"+Tracer+"_GC_vs_FTS.csv",toprint, delimiter=",", fmt="%s", header="Date_Time,"+Station+"_"+Inventory+"_GC_"+Tracer+"_tcol,FTS_"+Tracer+'_'+Station+"_tcol")

	toprint2 = []
	toprint2.append(gc_datetimes)
	toprint2.append(monthly_station_tcol)
	toprint2.append(avg_nd_xp)
	toprint2=np.transpose(np.array(toprint2))
	np.savetxt("Monthly_"+Station+"_"+Year+"_"+Inventory+"_"+Tracer+"_GC_vs_FTS.csv",toprint2, delimiter=",", fmt="%s", header="Date_Time,"+Station+"_"+Inventory+"_GC_"+Tracer+"_tcol,FTS_"+Tracer+'_'+Station+"_tcol")



# Plotting
	#plot vertical profile
#   	f = plt.figure(1,figsize=(10,10))
#   	plt.plot(range(0,12),darwin_gc_vmr_surf,'-k', linewidth=3, label='FTS')
#    	plt.ylabel("CO level")
#   	plt.xlabel("Month ppb")
#	plt.ylim((35,300))

	#plot emissions over map
#	f = plt.figure(2,figsize=(12,8))
#	m=Basemap(llcrnrlat=-49,  urcrnrlat=-8.9, 			# basemap over area of interest
#         	llcrnrlon=111.25, urcrnrlon=178.75,
#        	resolution='l',projection='merc',
#	       	lat_0=0, lon_0=0)
#	gc_lon,gc_lat = np.meshgrid(gc_lon,gc_lat) 			# lat lon are 1D, basemap uses 2D mesh
#	gc_xi, gc_yi = m(gc_lon,gc_lat)
#	cmap=plt.get_cmap('OrRd')
#	norm=colors.LogNorm(vmin=0.01, vmax=50)
#	cs = m.pcolormesh(gc_xi,gc_yi,np.squeeze(gc_tot_bbemissions),cmap=cmap,norm=norm) # draw the data onto the map
#	m.drawcoastlines()  						# add coastlines
#	m.drawparallels([0], labels=[0,0,0,0]) 				# draw equator, label nothing
#	m.drawparallels(np.arange(-89.,89.,8.),labels=[1,0,0,0]) 	# draw parallels
#	m.drawmeridians(np.arange(-178.75,178.75,10.),labels=[0,0,0,1]) # draw meridians
#	cb=m.colorbar(cs,"right",size="5%", pad="2%")			#add colorbar
#	cb.ax.tick_params(labelsize=15)
#	cb.set_label('Kg of '+Tracer.upper()+' emmited', fontsize=20,labelpad=10)
#	plt.suptitle(Inventory +' '+ Tracer.upper()+' Emission map', fontsize=20) 	#add Title
#	plt.savefig(Inventory+' '+Tracer.upper()+' emission.png', dpi = 100)

	#plot emission difference over map
#	f = plt.figure(3,figsize=(12,8))
#	m=Basemap(llcrnrlat=-49,  urcrnrlat=-8.9, 			# basemap over area of interest
#         	llcrnrlon=111.25, urcrnrlon=178.75,
#       	resolution='l',projection='merc',
#          	lat_0=0, lon_0=0)
#	cmap=plt.get_cmap('bwr')
#	#norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,vmin=-2, vmax=2)
#	cs = m.pcolormesh(gc_xi,gc_yi,np.squeeze(gc_emmission_dif/1e+3),cmap=cmap,vmin=-10, vmax=10) #,norm=norm)
#	m.drawcoastlines()  						# add coastlines
#	m.drawparallels([0], labels=[0,0,0,0]) 				# draw equator, label nothing
#	m.drawparallels(np.arange(-89.,89.,8.),labels=[1,0,0,0]) 	# draw parallels
#	m.drawmeridians(np.arange(-178.75,178.75,10.),labels=[0,0,0,1]) # draw meridians
#	cb=m.colorbar(cs,"right",size="5%", pad="2%")			#add colorbar
#	cb.ax.tick_params(labelsize=15)
#	cb.set_label('Difference in Kg of '+Tracer.upper()+' emmited', fontsize=20,labelpad=10)
#	plt.suptitle(Inventory +' vs GFED '+ Tracer.upper()+' Emissions', fontsize=20) 	#add Title
#	#plt.savefig(Inventory+' vs GFED '+Tracer.upper()+' emission.png', dpi = 100)
"""
