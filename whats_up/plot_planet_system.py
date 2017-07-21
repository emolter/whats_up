#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from ephemeris.ephemeris import Ephemeris
from matplotlib import text

'''This standalone function uses the framework provided by the Ephemeris class
to plot the position of moons in planetary systems at a given date and time.'''

planet = 'Uranus'
obs_code = 568
tstart = '2017-06-23 05:30'
tend = '2017-06-23 06:00' #this is not actually used, but needs to get put into call to Ephemeris. Should just put 30 mins after tstart
stepsize = '30 minutes' #this not used either
outfile = 'images/system_diagram.png'
fov = 1.0/60.0 #field of view of image in decimal degrees 

# sample planet systems of interest. Arbitrary target lists are possible
if planet == 'Mars':
    targetlist = ['Mars','Phobos','Deimos']
elif planet == 'Jupiter':
    targetlist = ['Jupiter','Io','Europa','Ganymede','Callisto']
elif planet == 'Saturn':
    targetlist = ['Saturn','Titan']
elif planet == 'Uranus':
    targetlist = ['Uranus','Oberon','Titania','Ariel','Umbriel','Miranda']
elif planet == 'Neptune':
    targetlist = ['Neptune','Triton']
elif planet == 'Pluto':
    targetlist = ['Pluto','Charon']
else:
    targetlist = [planet]

colors = ['k' for val in targetlist]


def coord_to_deg(coord):
    [d,m,s] = [float(val) for val in coord.split()]
    return np.sign(d)*(np.abs(d) + m/60.0 + s/3600.0)
    
def deg_to_coord(deg):
    d = int(deg)
    m = int(np.abs(deg - d)*60)
    s = (np.abs(deg - d)*60 - m)*60.0
    return '{: 03d}'.format(d) + ' ' + '{:02d}'.format(m) + ' ' + '{:04.1f}'.format(s)
    

fig,ax = plt.subplots(1,1, figsize = (10,10))
ax.set_aspect('equal')
for i in range(len(targetlist)):
    target = targetlist[i]
    ephem = Ephemeris(target,obs_code,tstart,tend,stepsize)
    ra = ephem.ra[0]
    dec = ephem.dec[0]
    coords = (15.0*coord_to_deg(ra),coord_to_deg(dec))
    r = ephem.ang_diam[0]/(2.0*3600.)
    
    if target == planet:
        xlim = [coords[0] - fov/2.0,coords[0]+fov/2.0]
        ylim = [coords[1] - fov/2.0,coords[1]+fov/2.0]
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        disk = plt.Circle(coords,r,color=colors[i], hatch = 'xxx', fill = False)
        ax.add_artist(disk)
    else:
        plt.plot(coords[0],coords[1],label=target,color=colors[i],linestyle = '', marker = '.')
        
    ax.text(coords[0] - r - fov/200., coords[1],target,color = colors[i], va='top',ha='left', fontsize = 12)
    

ax.set_xticklabels([deg_to_coord(val/15.0) for val in ax.get_xticks()]) #convert labels to hms
ax.set_yticklabels([deg_to_coord(val) for val in ax.get_yticks()]) #convert labels to dms
ax.set_xlabel('Right Ascension (HH MM SS)', fontsize = 14)
ax.set_ylabel('Declination (DD MM SS)', fontsize = 14)
ax.set_title(planet + ' system UTC '+tstart,y=0.93,fontsize=16)
plt.gca().invert_xaxis()

ax.annotate('N',xy = (0.15,0.05), 
                xytext = (0.15,0.15),
                textcoords = 'axes fraction', 
                xycoords='axes fraction', 
                color = 'k', 
                fontsize = 16,
                va = 'center', ha = 'center',
                arrowprops = dict(facecolor='k', arrowstyle = '<|-', linewidth = 3, fill=True))
ax.annotate('E',xy = (0.15,0.05), 
                xytext = (0.05,0.05),
                textcoords = 'axes fraction', 
                xycoords='axes fraction', 
                color = 'k', 
                fontsize = 16,
                va = 'center', ha = 'center',
                arrowprops = dict(facecolor='k', arrowstyle = '<|-', linewidth = 3, fill=True))

ax.get_yaxis().set_tick_params(which='both', direction='in',length=8, right=True)
ax.get_xaxis().set_tick_params(which='both', direction='in',length=8, top = True)

plt.tight_layout()
plt.savefig(outfile)
plt.show()
        