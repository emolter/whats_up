#!/usr/bin/env python

from urllib.request import urlopen
import numpy as np
from time import strftime, gmtime, time
from datetime import datetime,timedelta
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import sys

from ephemeris.ephemeris import Ephemeris
from ephemeris.read_scope_limits import read_scope_limits
from get_twilight import get_twilight


helpstr = """
    -----------------------------------------------------------------------
    whats_up -- Displays ephemeris information and creates Keck starlists for solar system targets
    
    Purpose:
            To streamline observation planning for non-sidereal targets.
    
    Limitations:
            Does not accept comets at this time.  Will hopefully update later.
            Does not accept time resolutions finer than 1 minute right now. Will hopefully update later.
            Does not accept J1950 times.
    
    Usage:
           	./whats_up.py input_file.txt
            
    Parameters:
            input file must be .txt or .dat
            See documentation and sample input file for parameter descriptions
     
    Output:
            Text to terminal.
            If you want text to a log file instead, use ./whats_up.py input_file.txt > logfile_name.txt
            Plots saved as .png
            (optional) Keck starlist as .txt
    
    Example:
            ./whats_up.py inputs.txt
    
    Help:
            ./whats_up.py -h
    
    Modification history:
           	2017-Jul-18     Ned Molter    Original version.
    -----------------------------------------------------------------------
    """

def main(argv):
    
    if argv == '-h' or argv == '-help':
        print(helpstr)
        sys.exit()
        
    else:
        infile = argv
    
    args = {}    
    with open(infile,'r') as f:
        for line in f:
            line = line.split('#',1)[0]
            if len(line) < 2:
                continue
            l = line.split()
            param = l[0].strip(', \n')
            arg = ' '.join(l[1:]).strip(', \n')
            args[param] = arg
    
    
    ###########################################
    ### Parse arguments and validate inputs ###
    ###### Error handling is not perfect ######    
    ###########################################
    
    try:
        targetlist = args['targetlist']
        targetlist = [s.strip(', \n') for s in targetlist.split()]
        if len(targetlist) == 0:
            raise
    except:
        sys.exit('ERROR: "targetlist" parameter not found or formatted incorrectly. See sample input file for help.')
    
    try:
        tstart = args['tstart']
        datetime.strptime(tstart,'%Y-%m-%d %H:%M')
    except:
        sys.exit('ERROR: "tstart" parameter not found or formatted incorrectly. See sample input file for help.')
    
    try:
        tend = args['tend']
        datetime.strptime(tend,'%Y-%m-%d %H:%M')
    except:
        sys.exit('ERROR: "tstart" parameter not found or formatted incorrectly. See sample input file for help.')
    
    dt = datetime.strptime(tend,'%Y-%m-%d %H:%M') - datetime.strptime(tstart,'%Y-%m-%d %H:%M')
    if dt.total_seconds() <= 0:
        sys.exit('ERROR: End time is before start time!')
    
    try:
        stepsize = args['stepsize']
    except:
        print('stepsize not specified. Using step size of 30 minutes.')
        stepsize = '30 minutes'
        
    try:
        min_ang_sep = args['min_ang_sep'].split()
        min_ang_sep = [float(s.strip(', ')) for s in min_ang_sep]
        if len(min_ang_sep) == 0:
            raise
        if len(min_ang_sep) == 1:
            min_ang_sep = [min_ang_sep[0] for val in targetlist]
    except:
        print('min_ang_sep not specified. Using min_ang_sep = 0.')
        min_ang_sep = [0.0 for val in targetlist]
        
    try:
        min_moon_dist = args['min_moon_dist'].split()
        min_moon_dist = [float(s.strip(', ')) for s in min_moon_dist]
        if len(min_moon_dist) == 0:
            raise
        if len(min_moon_dist) == 1:
            min_moon_dist = [min_moon_dist[0] for val in targetlist]
    except:
        print('moon_dist not specified or specified incorrectly. Using min_moon_dist = 5.0 degrees')
        min_moon_dist = [5.0 for val in targetlist]
        
    try:
        max_airmass = args['max_airmass'].split()
        max_airmass = [float(s.strip(', ')) for s in max_airmass]
        if len(max_airmass) == 0:
            raise
        if len(max_airmass) == 1:
            max_airmass = [max_airmass[0] for val in targetlist]
    except:
        print('max_airmass not specified or specified incorrectly. Using no max airmass')
        max_airmass = [None for val in targetlist]
        
    try:
        observatory_code = args['observatory_code']
        if len(observatory_code) == 0:
            raise
    except:
        sys.exit('ERROR: Observatory code could not be read. Should be valid JPL Horizons observatory ID. See sample input file for help.')
    
    try:
        ingest_limits = args['ingest_limits']
        if ingest_limits.lower() == 'true' or ingest_limits.lower() == 't' or ingest_limits == '1':
            ingest_limits = True
        elif ingest_limits.lower() == 'false' or ingest_limits.lower() == 'f' or ingest_limits == '0':
            ingest_limits = False
        else:
            raise
    except:
        print('ingest_limits not specified or could not be read. Setting ingest_limits = False.')
        ingest_limits = False
        
    if ingest_limits:
        try:
            limits_file = args['limits_file']
            if len(limits_file) == 0:
                raise
        except:
            print('ERROR: limits_file not found! Required if ingest_limits set to True.  See sample input file for help.')
        try:
            limits_pad = float(args['limits_pad'])
        except:
            print('limits_pad not specified or could not be read. Setting limits_pad = 1.0 degrees')
            limits_pad = 1.0
    else:
        limits_file = None
        limits_pad = None

    try:
        make_starlist = args['make_starlist']
        if make_starlist.lower() == 'true' or make_starlist.lower() == 't' or make_starlist == '1':
            make_starlist = True
        elif make_starlist.lower() == 'false' or make_starlist.lower() == 'f' or make_starlist == '0':
            make_starlist = False
        else:
            raise
    except:
        print('make_starlist not specified or could not be read. Setting make_starlist = False.')
        make_starlist = False
    if make_starlist:
        try:
            starlist_save_path = args['starlist_save_path']
            if len(starlist_save_path) == 0:
                raise 
        except:
            starlist_save_path = ''
        try:
            include_standards = args['include_standards']
            if include_standards.lower() == 'true' or include_standards.lower() == 't' or include_standards == '1':
                include_standards = True
            elif include_standards.lower() == 'false' or include_standards.lower() == 'f' or include_standards == '0':
                include_standards = False
            else:
                raise
        except:
            print('include_standards not specified or could not be read. Setting include_standards = False.')
            include_standards = False
        if include_standards:
            try:
                standards_file = args['standards_file']
                if len(standards_file) == 0:
                    raise
            except:
                print('ERROR: standards_file not found! Required if include_standards set to True.  See sample input file for help.')
        else:
            standards_file = None
    else:
        include_standards = None
        standards_file = None
        
    try:
        show_plots = args['show_plots']
        if show_plots.lower() == 'true' or show_plots.lower() == 't' or show_plots == '1':
            show_plots = True
        elif show_plots.lower() == 'false' or show_plots.lower() == 'f' or show_plots == '0':
            show_plots = False
        else:
            raise
    except:
        print('show_plots not specified or could not be read. Setting show_plots = False.')
        show_plots = False
        
    try:
        plot_save_path = args['plot_save_path']
        if len(plot_save_path) == 0:
            raise 
    except:
        plot_save_path = ''   

    print('-----------------------------') 
    
    #consolidate to make things easier to read   
    observation_params = [targetlist,tstart,tend,stepsize,min_ang_sep,max_airmass,min_moon_dist,make_starlist,starlist_save_path]
    observatory_params = [observatory_code,ingest_limits,limits_file,limits_pad]
    standards_params = [include_standards,standards_file]
    plot_params = [show_plots,plot_save_path]
     
    
    ########################################
    ### run the actual steps of the code ###
    ########################################
    
    obs = Observation(observation_params,observatory_params)
    obs.make_ephemerides()

    if make_starlist:
        for ephem in obs.ephemlist:
            ephem.generate_starlist(observatory_params, standards_params, starlist_save_path)
    
    obs.eval_observability(plot_params)
    print(' ')

 



class Observation:
    '''Handle ephemerides for multiple targets at once.
    Print pretty text to terminal about whether things are observable.
    Make pretty plots of whether things are observable.'''
    
    def __init__(self,observation_params,observatory_params):
        '''Night-specific, but not object-specific, parameters are set.
        observation_params = [targetlist,tstart,tend,stepsize,min_ang_sep,max_airmass,min_moon_dist,make_starlist,starlist_save_path]
        observatory_params = [observatory_code,ingest_limits,limits_file,limits_pad]
        see documentation for descriptions of each parameter'''
        self.targetlist = observation_params[0]
        self.tstart = observation_params[1]
        self.tend = observation_params[2]
        self.stepsize = observation_params[3]
        self.min_ang_sep = observation_params[4]
        self.max_airmass = observation_params[5]
        self.min_moon_dist = observation_params[6]
        self.make_starlist = observation_params[7]
        self.obs_code = observatory_params[0]
        self.ingest_limits = observatory_params[1]
        self.limits_file = observatory_params[2]
        self.limits_pad = observatory_params[3]
        
        moon = Ephemeris('Moon',self.obs_code,self.tstart,self.tend,self.stepsize)
        self.moon_loc = [moon.azimuth, moon.elevation]
        
    def make_ephemerides(self):        
        '''make ephemeris objects for each target for the specified time range'''
        ephemlist = [] 
        for obj in self.targetlist:
            planet_ephem = Ephemeris(obj,self.obs_code,self.tstart,self.tend,self.stepsize)
            ephemlist.append(planet_ephem)
        self.ephemlist = ephemlist
        
    def eval_observability(self, plot_params = [False,'']):
        '''Print whether each target observable over the specified time range.
        Plot two nice visualizations of this and save as .png, optionally show that
        plot to screen.
        plot_params = [show_plot,plot_save_path]
        show_plot: boolean, whether to show plot to screen
        plot_save_path: string, the directory into which go the plots. Defaults to 
            current working directory. Must end with a /'''
        show_plots = plot_params[0]
        plot_save_path = plot_params[1]
        
        for j in range(len(self.ephemlist)):
            ephem = self.ephemlist[j]
            min_moon_dist = self.min_moon_dist[j]
            min_ang_sep = self.min_ang_sep[j]
            max_airmass = self.max_airmass[j]
            
            print('--------------------------------------')
            print('----------------- '+ephem.target+' -----------------')
            ephem.observability(limits_file = self.limits_file, limits_pad = self.limits_pad, 
                    moon_loc=self.moon_loc, min_moon_dist = min_moon_dist, min_ang_sep=min_ang_sep, max_airmass = max_airmass)

        fontsize = 14
        alpha = 0.6
        params =   {'font.size':fontsize,
                    'legend.fontsize': fontsize,
                    'axes.labelsize': fontsize + 2,
                    'axes.titlesize': fontsize + 4,
                    'axes.titlepad': 10, 
                    'xtick.labelsize': fontsize,
                    'ytick.labelsize': fontsize,
                    'lines.linewidth': 2,
                    'lines.markersize': 12
                    }
        plt.rcParams.update(params)
        
        fig, ax = plt.subplots(1,1,subplot_kw={'projection':'polar'},figsize=(8,8))
        cmap = get_cmap('plasma')
        
        for j in range(len(self.ephemlist)):
            body = self.ephemlist[j]
            color = cmap(float(j)/len(self.ephemlist))
        
            #plot ranges where it's observable
            switches = body.switches
            if body.observable[0]:
                switches.insert(0,0)
            if body.observable[-1]:
                switches.append(len(body.observable)-1)
            
            dummy_bool = True
            for i in range(len(switches)):
                s = switches[i]
                t = body.times[s][-5:]
                azi = body.azimuth[s]
                elev = body.elevation[s]
                #label start and end points with time
                ax.plot(np.radians(azi),90 - elev,marker = 'o',markersize = 6,color = color)
                if i % 2 == 0.: #just spreading out text a bit to make easier to read
                    if j % 2 == 0.:
                        ax.text(np.radians(azi),90 - elev,t, color = color, va='top', ha = 'right')
                    else:
                        ax.text(np.radians(azi),90 - elev,t, color = color, va='bottom', ha = 'right')
                else:
                    if j % 2 == 0.:
                        ax.text(np.radians(azi),90 - elev,t, color = color, va='top', ha = 'left')
                    else:
                        ax.text(np.radians(azi),90 - elev,t, color = color, va='bottom', ha = 'left')
            
                if body.observable[s]:
                    try:
                        az = body.azimuth[s:switches[i+1]+1]
                        el = body.elevation[s:switches[i+1]+1]
                        if dummy_bool:
                            label = body.target
                        else:
                            label = ''
                        ax.plot(np.radians(az),90 - el, color=color, label=label)
                        dummy_bool = False
                    except:
                        pass 
        
        if self.ingest_limits:
            #plot the telescope limits
            [interp_l,interp_u] = read_scope_limits(self.limits_file)
            az = np.linspace(0,360,91)
            el_l = 90 - interp_l(az) #hack to make 90 (the zenith) in the center
            el_u = 90 - interp_u(az)
            az = np.radians(az)
            
            ax.plot(az,el_l,color = 'darkred', linestyle = '--', label = 'telescope limits', alpha = alpha)
            ax.plot(az,el_u,color = 'darkred', linestyle = '--', alpha = alpha)
            
        #hack to make elevation = 90 (the zenith) in the center    
        ax.set_yticks(range(-1, 90, 10))                   # Define the yticks
        ax.set_yticklabels(map(str, range(90, -1, -10)))   # Change the labels
        
        #make things pretty
        ax.set_theta_zero_location('N')
        ax.text(np.radians(10),100,'Ephemeris '+body.times[0][0:-5], fontsize = fontsize + 4, ha = 'right', bbox=dict(facecolor='none', edgecolor='k'))
        ax.text(np.radians(0-3), 90, 'N', fontsize = fontsize + 10, ha = 'left', va = 'bottom')
        ax.text(np.radians(90-3), 90, 'E', fontsize = fontsize + 10, ha = 'right', va = 'bottom')
        ax.text(np.radians(180-4), 90, 'S', fontsize = fontsize + 10, ha = 'right', va = 'top')
        ax.text(np.radians(270-3), 90, 'W', fontsize = fontsize + 10, ha = 'left', va = 'top')
        ax.legend()
        
        print(' ')
        plt.tight_layout()
        plt.savefig(plot_save_path+'az_alt.png',bbox = None)
        if show_plots:
            plt.show()
        plt.close()    
         
            
        ### The second, simpler plot - just a chart of what's observable when
        
        fig, ax = plt.subplots(1,1, figsize = (8,3))
        cmap = get_cmap('plasma')
        
        ax.set_ylim([-1,len(self.ephemlist)])
        ax.set_yticks(range(len(self.ephemlist)))
        ax.set_yticklabels(self.targetlist + [''])
        
        for j in range(len(self.ephemlist)):
            body = self.ephemlist[j]
            switches = body.switches
            switches.insert(0,0)
            switches.append(len(body.observable)-1)
            
            color = cmap(float(j)/len(self.ephemlist))
            t = body.times 
            hr = [float(val[-5:-3]) for val in t]
            mn = [float(val[-2:])/60.0 for val in t]
            t_float = np.add(hr,mn)

            for i in range(len(switches)):
                t_start = t_float[switches[i]]

                try:
                    t_end = t_float[switches[i+1]]
                except:
                    t_end = t_float[-1]
                    
                if body.observable[switches[i]]:
                    ax.fill_between([t_start,t_end],[j-0.35,j-0.35],[j+0.35,j+0.35],color=color)
                else:
                    ax.plot([t_start,t_end],[j,j],lw=2,color='k')
                
                if i != 0 and i != len(switches) - 1:
                    if body.observable[switches[i]] and switches[i] != switches[-1]:
                        ax.text(t_start,j+0.5, t[switches[i]][-5:], color = color, va='top', ha = 'right',fontsize=12)
                    else:
                        ax.text(t_start,j-0.1, t[switches[i]][-5:], color = color, va='top', ha = 'left',fontsize=12)
        '''        
        #get twilight times from Keck website and plot
        twilight_times = get_twilight(self.tstart[2:11])
        for val in ['sunset','dusk_18deg','dawn_18deg','sunrise']:
            t = twilight_times[val]
            label = val + '       '
            label = label[:7].strip(' ')
            hr = float(t[:2])
            mn = float(t[3:5])/60.0
            t_float = np.add(hr,mn)
            ax.axvline(t_float, linestyle = '--', color = 'steelblue')
            if val == 'sunset' or val == 'dawn_18deg':
                ax.text(t_float, 3.65, t[:-3], rotation = 0, color = 'steelblue', fontsize=10, ha = 'right')
                ax.text(t_float, -0.95, label, rotation = 0, color = 'steelblue', fontsize=10, ha = 'right') 
            else:
                ax.text(t_float, 3.65, t[:-3], rotation = 0, color = 'steelblue', fontsize=10, ha = 'left')
                ax.text(t_float, -0.95, label, rotation = 0, color = 'steelblue', fontsize=10, ha = 'left')             
        '''            
        ax.set_title('Observable Targets '+self.tstart[:-6])
        ax.set_xlabel('UT Time (Hours)')
        plt.tight_layout()
        plt.savefig(plot_save_path+'observability.png')
        if show_plots:
            plt.show()
        plt.close()
    
    
               
if __name__ == "__main__":
    try:
        sys.argv[1]
    except:        
        sys.exit('ERROR: Input file required. ./whats_up.py -h for help.')
    main(sys.argv[1])
     