import numpy as np
from datetime import datetime, timedelta

def coord_to_deg(coord):
    [d,m,s] = [float(val) for val in coord.split()]
    return np.sign(d)*(np.abs(d) + m/60.0 + s/3600.0)

def find_lst(time,long_east):
    '''time in UTC, longitude to the east, positive or negative.
    Does not include Equation of Time so error is a few minutes'''
    j2000 = datetime.strptime('2000-01-01 12:00:00','%Y-%m-%d %H:%M:%S')
    time = time.strip()
    try:
        ut_time = datetime.strptime(time,'%Y-%m-%d %H:%M')
    except:
        try:
            ut_time = datetime.strptime(time,'%Y-%b-%d %H:%M')
        except:
            sys.exit('Could not read input time. Ensure format correct')
    
    diff_days = (ut_time - j2000)/timedelta(days=1)

    #from http://aa.usno.navy.mil/faq/docs/GAST.php
    gmst_hrs = (18.697374558 + 24.06570982441908*diff_days)
    lst = (gmst_hrs + long_east/15)%24
    
    return lst
    
def alt_az_pa(time,ra,dec,lat,long_east):
    '''Calculates altitude, azimuth, parallactic angle.
    Input RA as string 'hh mm ss.ss', dec as string 'dd mm ss.s', rest as degrees
    Does not include altitude of telescope, so some error'''
    lst = find_lst(time,long_east)
    ra = coord_to_deg(ra)
    dec = coord_to_deg(dec)
    ha = 15*(lst - ra)

    ra = np.radians(ra)
    dec = np.radians(dec)
    ha = np.radians(ha)
    lat = np.radians(lat)
    
    #from http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm
    sinh = np.sin(lat)*np.sin(dec) + np.cos(lat)*np.cos(dec)*np.cos(ha)
    alt = np.arcsin(sinh)
    cos_az = (np.sin(dec) - np.sin(lat)*np.sin(alt))/(np.cos(alt)*np.cos(lat))
    if alt < 0 :
        az = np.arccos(-cos_az) + np.pi
    else:
        az = np.arccos(cos_az)
    
    sinb = -np.sin(az)*np.cos(lat)/np.cos(dec)
    cosb = -np.cos(ha)*np.cos(az) - np.sin(ha)*np.sin(az)*np.sin(lat)
    pa = np.arctan(cosb/sinb)
    
    alt = np.degrees(alt)
    az = np.degrees(az)
    pa = np.degrees(pa)
    
    return alt, az, pa

## test cases
#t_lst = lst('2017-06-16 02:00',360-155.4747)  
#print(alt_az_pa('12 50 13.72','-03 53 46.3',19.8260,t_lst)) #should be 32,108
#print(alt_az_pa('23 02 31.27','-07 07 02.0',19.8260,t_lst)) #should be -61, 292
