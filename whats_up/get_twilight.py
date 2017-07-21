#!/usr/bin/env python

from urllib.request import urlopen

def get_twilight(date):
    '''Queries Keck website to get twilight time data'''
    twilight_times = {}
    '''date in format YY-MM-DD'''
    base_url = 'https://www.keck.hawaii.edu/observing/schedule/ws/telsched.php?field=Twilight&verbosity=0&date='
    url = base_url + date
    dataline = urlopen( url ).readlines()[0]
    for entry in dataline.split(','):
        entry = entry.strip(''' ` }{ ' ''')
        key, val = entry.split(":'")
        twilight_times[key] = val
    return twilight_times