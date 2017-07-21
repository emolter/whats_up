#!/usr/bin/env python
from urllib.request import urlopen
import numpy as np
from time import strftime, gmtime, time
from datetime import datetime,timedelta
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import sys

from .get_ephem import get_ephemerides, naif_lookup
from .alt_az_pa import alt_az_pa
from .read_scope_limits import read_scope_limits

##################################
######## helper functions ########
################################## 

def coord_to_deg(coord):
    [d,m,s] = [float(val) for val in coord.split()]
    return np.sign(d)*(np.abs(d) + m/60.0 + s/3600.0)
    
def deg_to_coord(deg):
    d = int(deg)
    m = int(np.abs(deg - d)*60)
    s = (np.abs(deg - d)*60 - m)*60.0
    return '{: 03d}'.format(d) + ' ' + '{:02d}'.format(m) + ' ' + '{:04.1f}'.format(s)

def spherical_cosines(phi1,theta1,phi2,theta2):
    '''Returns the angular distance between two points in lat/lon or alt/az
    (phi1, theta1) and (phi2, theta2).'''
    phi1, theta1, phi2, theta2 = np.radians(phi1), np.radians(theta1), np.radians(phi2), np.radians(theta2)
    dtheta = theta2 - theta1
    ans = np.abs(np.arccos(np.sin(phi1)*np.sin(phi2) + np.cos(phi1)*np.cos(phi2)*np.cos(dtheta)))
    return np.degrees(ans)   

def message_lookup(mess_code, elev, airmass, ang_sep, moon_dist):
    '''Helper function to observability(); decodes message codes
    o = elevation ok, m = elevation below min, x = elevation above max
    o = airmass ok, x = airmass above max
    o = separation ok, m = separation below min, z = this is parent body, or no data
    o = visibility ok, e = eclipse, c = occultation, t = transit
    example: 'oozc' -> elevation ok, airmass ok, separation not specified, body is occulted by parent
    '''
    #elevation
    message = []
    if mess_code[0] == 'o':
        message.append('Elevation is ok: '+str(elev)+' degrees.')
    elif mess_code[0] == 'm':
        message.append('Elevation is too low: '+str(elev)+' degrees.')
    elif mess_code[0] == 'x':
        message.append('Elevation is too near zenith: '+str(elev)+' degrees.')

    #airmass
    if mess_code[1] == 'o':
        message.append('Airmass is ok: '+str(airmass))
    elif mess_code[1] == 'x':
        message.append('Airmass is too high: '+str(airmass))

    #angular sep
    if mess_code[2] == 'o':
        message.append('Angular separation from host planet is ok: angular separation = '+str(ang_sep)+' arcsec.')
    elif mess_code[2] == 'm':
        message.append('Too close to host planet: angular separation = '+str(ang_sep)+' arcsec.')
    elif mess_code[2] == 'x':
        message.append('Minimum angular separation from parent body not specified (or target is itself a parent body)')
    
    #occultations and eclipses
    if mess_code[3] == 'o':
        message.append('No occultations or eclipses.')
    elif mess_code[3] == 'e':
        message.append('Target is in total or partial eclipse!')
    elif mess_code[3] == 'c':
        message.append('Target is occulted by parent!')
    elif mess_code[3] == 't':
        message.append('Target is transiting parent!')

    #distance from moon
    if mess_code[4] == 'o':
        message.append('Moon distance is ok.')
    elif mess_code[4] == 'm':
        message.append('Target is too close to the moon: angular distance = '+str(moon_dist)+' degrees.')
    elif mess_code[4] == 'x':
        message.append('No minimum moon distance specified.')
        
    return message
    
def read_ephem_line(arr):
    '''Helper to ephemeris.__init__. Converts ephemeris data to float, putting np.nan for "n.a."'''
    arr_out = []
    for s in arr:
        if s.strip() == 'n.a.':
            arr_out.append(np.nan)
        else:
            arr_out.append(float(s))
    return np.asarray(arr_out) 



### Class ###

class Ephemeris():
    '''Functions relevant to an ephemeris for a single target'''
    
    def __init__(self,target,obs_code,tstart,tend,stepsize):
        '''Immediately run get_ephemerides, then set a bunch of 
        class variables corresponding to different information found in
        the ephemeris.
        '''
        self.target = target
        self.obs_code = obs_code
        self.ephem, self.observatory_coords = get_ephemerides(naif_lookup(self.target),self.obs_code,tstart,tend,stepsize)
        self.times = self.ephem[:,0]
        self.sun = [s.strip() for s in self.ephem[:,1]]
        self.moon = [s.strip() for s in self.ephem[:,2]]
        self.ra = self.ephem[:,3] #dms
        self.dec = self.ephem[:,4] #dms
        self.dra = np.asarray([float(s) for s in self.ephem[:,5]]) #arcsec hr-1 
        self.ddec = np.asarray([float(s) for s in self.ephem[:,6]]) #arcsec hr-1
        self.azimuth = np.asarray([float(s) for s in self.ephem[:,7]]) #degrees, North = 0 = 360
        self.elevation = np.asarray([float(s) for s in self.ephem[:,8]]) #degrees, above horizon
        self.airmass = read_ephem_line(self.ephem[:,9]) 
        self.extinction = read_ephem_line(self.ephem[:,10]) #magnitudes; currently not used
        self.vmag = read_ephem_line(self.ephem[:,11]) 
        self.sbrt = read_ephem_line(self.ephem[:,12]) #currently not used
        self.ang_sep = read_ephem_line(self.ephem[:,13]) #arcsec
        self.visibility = [s.strip(' ') for s in self.ephem[:,14]]
        self.ang_diam = read_ephem_line(self.ephem[:,15]) #arcsec
        
        ##planet orientation information only used for plot_planet_system
        self.ob_lon = read_ephem_line(self.ephem[:,16]) #degrees, positive to west
        self.ob_lat = read_ephem_line(self.ephem[:,17]) #degrees
        self.np_ang = read_ephem_line(self.ephem[:,18]) #degrees
        self.np_dist = read_ephem_line(self.ephem[:,19]) #arcsec
        
    def observability(self,limits_file=None, limits_pad = 1.0, moon_loc = None, min_moon_dist = None, min_ang_sep = None, max_airmass = None):
        '''Evaluate at what times over the given time range
        the target can be observed. Compares location on sky at each time step
        to the telescope limits and airmass limits; evaluates whether object 
        is in an eclipse or occultation or simply too close to parent body
        scope_limits, if provided, given as [alt_lower_interp_func, alt_upper_interp_func] as output by read_scope_limits.py
        moon_loc given as [[moon azimuth list],[moon elev list]] as output by observability.__init__
        '''
        
        print('--------------------------------------')
        self.min_ang_sep = min_ang_sep
        self.min_moon_dist = min_moon_dist
        if moon_loc:
            moon_azi = np.asarray(moon_loc[0])
            moon_elev = np.asarray(moon_loc[1])
            self.moon_dist = spherical_cosines(self.elevation,self.azimuth,moon_elev,moon_azi)
        
        observable = []
        message_list = [] 
        
        print('Start of ephemeris: '+self.times[0])
        
        for i in range(len(self.times)):
            truths = []
            mess_code = ''
            
            elevation = self.elevation[i]
            azimuth = self.azimuth[i]
            ang_sep = self.ang_sep[i]
            visibility = self.visibility[i]
            moon_dist = self.moon_dist[i]
            airmass = self.airmass[i]
            
            #check azimuth to initialize elevation limits
            if limits_file:
                [interp_l, interp_u] = read_scope_limits(limits_file)
                min_elev = interp_l(azimuth) + limits_pad
                max_elev = interp_u(azimuth) - limits_pad
            else:
                min_elev = 0.0
                max_elev = 90.0
                
            #check elevation
            if min_elev < elevation < max_elev:
                mess_code += 'o'
                truths.append(True)  
            elif elevation <= min_elev:
                truths.append(False)
                mess_code += 'm'
            elif elevation >= max_elev:
                truths.append(False)
                mess_code += 'x'                

            #check airmass
            if not max_airmass:
                max_airmass = 999
            if airmass < max_airmass:
                mess_code += 'o'
                truths.append(True)
            else:
                mess_code += 'x'
                truths.append(False)
            
            #check distance from parent body
            if self.min_ang_sep != None:
                if ang_sep > self.min_ang_sep:
                    truths.append(True)
                    mess_code += 'o'
                else:
                    truths.append(False)
                    mess_code += 'm'
            else:
                mess_code += 'x'
            
            #check for occultations, eclipses by parent body
            if visibility != '-' and visibility != 'n.a.': #check whether it even has a parent body
                if visibility == '*':
                    truths.append(True)
                    mess_code += 'o'
                else:
                    truths.append(False)
                    if visibility == 'p' or visibility == 'u':
                        mess_code += 'e'
                    if visibility == 'O' or visibility == 'P' or visibility == 'U':
                        mess_code += 'c'
                    if visibility == 't':
                        mess_code += 't'
            else:
                mess_code += 'x'
                
            ##check location of moon is farther than minimum
            if min_moon_dist:
                if not moon_loc:
                    sys.exit('ERROR: A minimum moon distance was specified, but the location of the moon was not passed into observability()!')
                
                if moon_dist <= self.min_moon_dist:
                    truths.append(False)
                    mess_code += 'm'
                else:
                    mess_code += 'o'
            else:
                mess_code += 'x'
            
            #Assess overall observability
            message_list.append(mess_code)
            if np.all(truths):
                observable.append(True)
            else:
                observable.append(False)
        

        #figure out when things go from not observable to observable
        observable = np.asarray(observable)
        switches = np.where(observable[:-1] != observable[1:])[0]
        switches = [val + 1 for val in switches] #easier if it's the first index of the new state
        
        #print statements for start of time range
        if observable[0]:
            print('********** '+self.target+' can be observed at start of time range ('+self.times[0]+')! **********')
            message = message_lookup(message_list[0],self.elevation[0],self.airmass[0],self.ang_sep[0],self.moon_dist[0])
            for mess in message:
                print(mess)          
        else:
            print('* '+self.target+' CANNOT be observed at start of time range ('+self.times[0]+') *')
            message = message_lookup(message_list[0],self.elevation[0],self.airmass[0],self.ang_sep[0],self.moon_dist[0])
            for mess in message:
                print(mess)
        
        #print statements for changes between observable and not
        for val in switches:
            if observable[val]:
                print('********** '+self.target+' can be observed starting at: '+self.times[val]+' **********')
                message = message_lookup(message_list[val],self.elevation[val],self.airmass[val],self.ang_sep[val],self.moon_dist[val])
                for mess in message:
                    print(mess)
            else:
                print('* '+self.target+' can NO LONGER be observed starting at: '+self.times[val]+' *')
                message = message_lookup(message_list[val],self.elevation[val],self.airmass[val],self.ang_sep[val],self.moon_dist[val])
                for mess in message:
                    print(mess)
        
        print('End of ephemeris: '+self.times[-1])
        
        self.observable = observable
        self.message_list = message_list
        self.switches = switches #indices where it switches from observable to not or vice versa. useful when plotting
       
    def generate_starlist(self, observatory_params, standards_params = [False,None], starlist_save_path=''):
        '''Put ephemeris info into format readable by Keck, output as .txt.
        Note this will make list for entire time range, not just observable times.
        This is by design - more flexible this way - but be careful!'''
        
        date = self.times[0][:12].strip()
        outfname = 'starlist_'+self.target+'_'+date+'.txt'
        self.starlist_name = outfname
        
        #pad target name with whitespace, then truncate to be correct length
        handle = self.target + '                 '
        handle = handle[:9]
        handle = handle[0].upper() + handle[1:] #target name must be capitalized for AO acquisition software to recognize extended sources automatically
        
        firstline = 'record divider  00 00 00.00 +00 00 00.0 2000.0 #####################\n'
        
        with open(starlist_save_path+outfname,'w') as f:
            f.write(firstline)
            for i in range(len(self.times)):
                t = self.times[i]
                name = handle + ' ' + t[-5:] #identifier includes time
                ra = self.ra[i]
                dec = self.dec[i]
                dra = str(self.dra[i]/15.0)[:6] #positive implies moving east, divide by 15 because that's the format Keck takes
                ddec = str(self.ddec[i])[:5]
                vmag = str(self.vmag[i])[:4]
                
                rmag = vmag #for now. all planets are ~zero color, and the rmag is just used for the waveguide sensor
                              
                line = name + ' ' + ra + ' ' + dec + ' 2000.0 dra='+dra+' ddec='+ddec+' rotmode=pa rotdest=0'+' vmag='+vmag+' rmag='+rmag+'\n'
                f.write(line)
            if standards_params[0]:
                standard_lines = self.find_standard_stars(observatory_params,standards_params[1])
                for ln in standard_lines:
                    f.write(ln)
                
    def find_standard_stars(self,observatory_params, standard_star_file ,moon_loc = None, min_moon_dist = None):
        '''Query database to find standard stars that are near enough to target,
        append those to the starlist in Keck-readable format'''
        if observatory_params[1]:
            scope_limits_file = observatory_params[2]
            limits_pad = observatory_params[3]
        else:
            scope_limits_file = None
            limits_pad = 1.0
            
        loc_target = (15.0*coord_to_deg(self.ra[-1]),coord_to_deg(self.dec[-1]))
        dists = []
        not_up = []
        standardlist = []
        with open(standard_star_file,'r') as f:
            f.readline()
            for line in f:
                standardlist.append(line)
                #name = line[:15]
                ra = line[15:27]
                dec = line[28:39]
                elev, azi, pa = alt_az_pa(self.times[-1],ra,dec,self.observatory_coords[0],self.observatory_coords[1])
                
                #check azimuth to initialize elevation limits
                if scope_limits_file:
                    [interp_l, interp_u] = read_scope_limits(scope_limits_file)
                    min_elev = interp_l(azi) + limits_pad
                    max_elev = interp_u(azi) - limits_pad
                else:
                    min_elev = 0.0
                    max_elev = 90.0
            
                #check elevation
                if min_elev < elev < max_elev:
                    not_up.append(False)
                else:
                    not_up.append(True)
                
                #moon avoid
                if min_moon_dist:
                    if not moon_loc:
                        sys.exit('ERROR: A minimum moon distance was specified, but the location of the moon was not passed into find_standard_stars()!')
                
                    if moon_dist <= self.min_moon_dist:
                        not_up[-1] = False                    
                
                loc = (15.0*coord_to_deg(ra),coord_to_deg(dec))
                dist = spherical_cosines(loc[1],loc[0],loc_target[1],loc_target[0])
                dists.append(dist)
        
        dists = np.asarray(dists)
        not_up = np.asarray(not_up)
        dists[not_up] = 100000 #make them really far away so min doesn't find them
        
        if np.min(dists) < 100000:
            closest_three = []
            for j in range(3):
                closest_i = np.argmin(dists)
                closest_star = standardlist[closest_i]
                closest_three.append(closest_star)
                dists[closest_i] = 100000
                if not np.min(dists) < 100000:
                    print('Only found %d stars!' %float(j+1))
                    return closest_three
            return closest_three
        else:
            sys.exit('ERROR: Standard star find failed!  Could not find any observable stars.')