###script to get vlass observing and QL imaging status

import requests, numpy as np
from io import BytesIO
#from astropy.io import ascii
from astropy.table import Table

outputfile = '../vlass_status.fits'

def parse_row(row, names):
    'assumes row is a list of column values for a row'
    rdata = []
    for i in range(len(names)):
        if i > len(row)-1:
            rdata.append('')
        else:
            rdata.append(row[i])
    if len(row)>len(names):
        rdata[-1] = ' '.join(row[len(names)-1:])
    return rdata




def get_vlass_status(url = 'https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php',
                     startrow=3, round2=2,
                     colinfo={'Tile': {'dtype': str, 'unit': None},
                              'Dec_min': {'dtype': float, 'unit': 'deg'},
                              'Dec_max': {'dtype': float, 'unit': 'deg'},
                              'RA_min': {'dtype': float, 'unit': 'hour'},
                              'RA_max': {'dtype': float, 'unit': 'hour'},
                              'epoch': {'dtype': str, 'unit': None},
                              'obsdate': {'dtype': str, 'unit': None},
                              'status': {'dtype': str, 'unit': None}}):
    'obtain latest status of VLASS observing and QL imaging'
    names = list(colinfo.keys())
    response = requests.get(url)
    if response.status_code==200:
        bytes = BytesIO(response.content)
        dstring = bytes.read().decode()
        dlines = dstring.split('\n')
        dlines = dlines[startrow:]
        datarows = []
        for line in dlines:
            if len(line)>0:
                datarows.append(parse_row(line.split(), names=names))
        datarows = np.array(datarows)
        data = Table(datarows, names=names)
        for col in data.colnames:
            if colinfo[col]['dtype'] == float:
                data[col] = np.array(data[col]).astype(colinfo[col]['dtype'])
                data[col] = np.round(np.array(data[col]), round2)
            data[col].unit = colinfo[col]['unit']
        
    else:
        print(f'Warning: unable to obtain updated VLASS status from url: {url}')
        data = None
        
    return data


vdat = get_vlass_status()

vdat.write(outputfile, overwrite=True) 
