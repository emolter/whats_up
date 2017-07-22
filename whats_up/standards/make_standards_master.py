import numpy as np
from urllib.request import urlopen
from astroquery.simbad import Simbad
from astropy.table import Table

do_brights, do_faints = True, True


def generate_starlist(table,outfname):
    '''Put astropy table info into format readable by Keck, output as .txt.
    '''
    
    firstline = 'record divider  00 00 00.00 +00 00 00.0 2000.0 #####################\n'
    
    with open(outfname,'w') as f:
        f.write(firstline)
        for i in range(len(table['MAIN_ID'])):
            #pad target name with whitespace, then truncate to be correct length
            handle = table['MAIN_ID'][i].decode('utf-8').replace(' ','') + '               '
            name = handle[:15]

            ra = table['RA'][i]+'00'
            dec = table['DEC'][i]+'00'
            ra = ra[:11]
            dec = dec[:11]
            rmag = str(table['FLUX_R'][i])
            vmag = str(table['FLUX_V'][i])
            #rmag must exist, but need only be approximate (used to close AO loops)
            if rmag == '' or rmag == '--' or rmag == np.nan:
                rmag = vmag
            hmag = str(table['FLUX_H'][i])
            jmag = str(table['FLUX_J'][i])
            kmag = str(table['FLUX_K'][i])
            try:
                lmag = str(table['FLUX_L'][i])
                mmag = str(table['FLUX_M'][i])
            except:
                lmag = '--'
                mmag = '--'
            sptype = str(table['SP_TYPE'][i].decode("utf-8"))
            
            line = name + ' ' + ra + ' ' + dec + ' 2000.0 rotmode=pa rotdest=0 vmag='+vmag+' rmag='+rmag+' hmag='+hmag+' jmag='+jmag+' kmag='+kmag+' lmag='+lmag+' mmag='+mmag+' sptype='+sptype+'\n'
            f.write(line)
                


#gemini stuff from text file gemini_phot_standards.txt
def gemini_bright():
    with open('gemini_phot_standards.txt','r') as f:
        header = f.readline()
    
        dist, name, ra, dec, y_mag, j_mag, h_mag, k_mag, l_mag, m_mag, catalog = [],[],[],[],[],[],[],[],[],[],[]
        for line in f:
            l = line.split()
            dist.append(float(l[0].strip(', \n')))
            name.append(l[1].strip(', \n'))
            ra.append(l[2].strip(', \n'))
            dec.append(l[3].strip(', \n'))
            y_mag.append(float(l[4].strip(', \n')))
            j_mag.append(float(l[5].strip(', \n')))
            h_mag.append(float(l[6].strip(', \n')))
            k_mag.append(float(l[7].strip(', \n')))
            l_mag.append(float(l[8].strip(', \n')))
            m_mag.append(float(l[9].strip(', \n')))
            catalog.append(l[10].strip(', \n'))
            
    name = np.asarray(name)
    l_mag = np.asarray(l_mag)
    m_mag = np.asarray(m_mag)
    l_min = 8.0

    ok = np.where( (l_mag > 0.0) & (m_mag > 0.0) & (l_mag < l_min) )

    name_ok = name[ok]
    l_mag_ok = l_mag[ok]
    m_mag_ok = m_mag[ok]
     
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('sptype','flux(U)', 'flux(V)', 'flux(B)', 'flux(R)', 'flux(H)', 'flux(J)', 'flux(K)')
    result_table = custom_simbad.query_objects(name_ok, wildcard=True)
    
    result_table['FLUX_L'] = l_mag_ok
    result_table['FLUX_M'] = m_mag_ok
    
    #ensure there's an A star
    sptypes = np.asarray(result_table['SP_TYPE'])
    sptruths = [val.startswith(b'A') for val in sptypes]
    return result_table[sptruths]

def gemini_faint():
    
    with open('gemini_phot_standards.txt','r') as f:
        header = f.readline()
    
        dist, name, ra, dec, y_mag, j_mag, h_mag, k_mag, l_mag, m_mag, catalog = [],[],[],[],[],[],[],[],[],[],[]
        for line in f:
            l = line.split()
            dist.append(float(l[0].strip(', \n')))
            name.append(l[1].strip(', \n'))
            ra.append(l[2].strip(', \n'))
            dec.append(l[3].strip(', \n'))
            y_mag.append(float(l[4].strip(', \n')))
            j_mag.append(float(l[5].strip(', \n')))
            h_mag.append(float(l[6].strip(', \n')))
            k_mag.append(float(l[7].strip(', \n')))
            l_mag.append(float(l[8].strip(', \n')))
            m_mag.append(float(l[9].strip(', \n')))
            catalog.append(l[10].strip(', \n'))    
    
    name = [string.replace('-','') for string in name]
    name = np.asarray(name)
    j_mag = np.asarray(j_mag)
    h_mag = np.asarray(h_mag)
    k_mag = np.asarray(k_mag)
    
    ok = np.where( (k_mag > 0 ) & (j_mag > 0 ) & (h_mag > 0 ) & ([string.startswith('FS') for string in name]) )[0]
    name_ok = name[ok]
    k_ok = k_mag[ok]
    h_ok = h_mag[ok]
    j_ok = j_mag[ok]
    
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('sptype','flux(U)', 'flux(V)', 'flux(B)', 'flux(R)')
    result_table = custom_simbad.query_objects(name_ok, wildcard = True)

    ## MAKE SURE THAT THE SIMBAD QUERY DOES NOT RANDOMIZE ORDER!
    #ids = np.asarray(result_table['MAIN_ID'])
    #for j in range(len(ids)):
    #    print(ids[j],name_ok[j])
    # this test shows that it retains its order - good
    
    #add gemini fluxes back in
    result_table['FLUX_K'] = k_ok
    result_table['FLUX_H'] = h_ok
    result_table['FLUX_J'] = j_ok
    
    #ensure v mag exists
    vmags = np.asarray(result_table['FLUX_V'])
    goodones = np.logical_and( np.invert((np.isnan(vmags))), vmags < 14.0) #warning from evaluating < on the nans, but doesn't break anything

    return result_table[goodones]


if do_brights:
    brightfname = 'standards_bright.txt' #have L and M mags
    brights = gemini_bright() 
    generate_starlist(brights,brightfname)
    brights.write('fancy_brights.txt', format='ascii', overwrite=True)  
    
    
if do_faints:
    faintfname = 'standards_faint.txt'
    faints = gemini_faint()
    generate_starlist(faints,faintfname)
    faints.write('fancy_faints.txt', format='ascii', overwrite=True)
     
'''    
custom_simbad = Simbad()
custom_simbad.add_votable_fields('sptype','flux(U)', 'flux(V)', 'flux(B)', 'flux(R)', 'flux(I)', 'flux(J)', 'flux(K)')
result_table = custom_simbad.query_objects(name, wildcard=True)
'''