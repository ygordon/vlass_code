###tidy up cirada csv files and add metadata -- makes a more user friendly file

import numpy as np, argparse
from astropy.table import Table, MaskedColumn
from astropy import units as u

###########################################################
###########################################################
###setup/parameters

sunit = u.Unit('mJy')
sbunit = u.Unit('mJy / beam')
posunit = u.Unit('deg')
sepsizeunit = u.Unit('arcsec')
pxunit = u.Unit('pix')


coldict = {'Component_name': {'dtype': '<U31',
                              'unit': None},
           'Component_id': {'dtype': 'int64',
                            'unit': None},
           'Isl_id': {'dtype': 'int64',
                      'unit': None},
           'RA': {'dtype': 'float64',
                  'unit': posunit,
                  'round2': 6},
           'DEC': {'dtype': 'float64',
                   'unit': posunit,
                   'round2': 6},
           'E_RA': {'dtype': 'float64',
                    'unit': posunit,
                    'round2': 6},
           'E_DEC': {'dtype': 'float64',
                     'unit': posunit,
                     'round2': 6},
           'Total_flux': {'dtype': 'float64',
                          'unit': sunit,
                          'round2': 2},
           'E_Total_flux': {'dtype': 'float64',
                            'unit': sunit,
                            'round2': 2},
           'Peak_flux': {'dtype': 'float64',
                         'unit': sbunit,
                         'round2': 2},
           'E_Peak_flux': {'dtype': 'float64',
                           'unit': sbunit,
                           'round2': 2},
           'Maj': {'dtype': 'float64',
                   'unit': sepsizeunit,
                   'round2': 2},
           'E_Maj': {'dtype': 'float64',
                     'unit': sepsizeunit,
                     'round2': 2},
           'Min': {'dtype': 'float64',
                   'unit': sepsizeunit,
                   'round2': 2},
           'E_Min': {'dtype': 'float64',
                     'unit': sepsizeunit,
                     'round2': 2},
           'PA': {'dtype': 'float64',
                  'unit': posunit,
                  'round2': 1},
           'E_PA': {'dtype': 'float64',
                    'unit': posunit,
                    'round2': 1},
           'Isl_Total_flux': {'dtype': 'float64',
                              'unit': sunit,
                              'round2': 2},
           'E_Isl_Total_flux': {'dtype': 'float64',
                                'unit': sunit,
                                'round2': 2},
           'Isl_rms': {'dtype': 'float64',
                       'unit': sbunit,
                       'round2': 2},
           'Isl_mean': {'dtype': 'float64',
                        'unit': sbunit,
                        'round2': 2},
           'Resid_Isl_rms': {'dtype': 'float64',
                             'unit': sbunit,
                             'round2': 2},
           'Resid_Isl_mean': {'dtype': 'float64',
                              'unit': sbunit,
                              'round2': 2},
           'RA_max': {'dtype': 'float64',
                      'unit': posunit,
                      'round2': 6},
           'DEC_max': {'dtype': 'float64',
                       'unit': posunit,
                       'round2': 6},
           'E_RA_max': {'dtype': 'float64',
                        'unit': posunit,
                        'round2': 6},
           'E_DEC_max': {'dtype': 'float64',
                         'unit': posunit,
                         'round2': 6},
           'S_Code': {'dtype': '<U1',
                      'unit': None},
           'Xposn': {'dtype': 'float64',
                     'unit': pxunit,
                     'round2': 2},
           'E_Xposn': {'dtype': 'float64',
                       'unit': pxunit,
                       'round2': 2},
           'Yposn': {'dtype': 'float64',
                     'unit': pxunit,
                     'round2': 2},
           'E_Yposn': {'dtype': 'float64',
                       'unit': pxunit,
                       'round2': 2},
           'Xposn_max': {'dtype': 'float64',
                         'unit': pxunit,
                         'round2': 2},
           'E_Xposn_max': {'dtype': 'float64',
                           'unit': pxunit,
                           'round2': 2},
           'Yposn_max': {'dtype': 'float64',
                         'unit': pxunit,
                         'round2': 2},
           'E_Yposn_max': {'dtype': 'float64',
                           'unit': pxunit,
                           'round2': 2},
           'Maj_img_plane': {'dtype': 'float64',
                             'unit': sepsizeunit,
                             'round2': 2},
           'E_Maj_img_plane': {'dtype': 'float64',
                               'unit': sepsizeunit,
                               'round2': 2},
           'Min_img_plane': {'dtype': 'float64',
                             'unit': sepsizeunit,
                             'round2': 2},
           'E_Min_img_plane': {'dtype': 'float64',
                               'unit': sepsizeunit,
                               'round2': 2},
           'PA_img_plane': {'dtype': 'float64',
                            'unit': posunit,
                            'round2': 1},
           'E_PA_img_plane': {'dtype': 'float64',
                              'unit': posunit,
                              'round2': 1},
           'DC_Maj': {'dtype': 'float64',
                      'unit': sepsizeunit,
                      'round2': 2},
           'E_DC_Maj': {'dtype': 'float64',
                        'unit': sepsizeunit,
                        'round2': 2},
           'DC_Min': {'dtype': 'float64',
                      'unit': sepsizeunit,
                      'round2': 2},
           'E_DC_Min': {'dtype': 'float64',
                        'unit': sepsizeunit,
                        'round2': 2},
           'DC_PA': {'dtype': 'float64',
                     'unit': posunit,
                     'round2': 1},
           'E_DC_PA': {'dtype': 'float64',
                       'unit': posunit,
                       'round2': 1},
           'DC_Maj_img_plane': {'dtype': 'float64',
                                'unit': sepsizeunit,
                                'round2': 2},
           'E_DC_Maj_img_plane': {'dtype': 'float64',
                                  'unit': sepsizeunit,
                                  'round2': 2},
           'DC_Min_img_plane': {'dtype': 'float64',
                                'unit': sepsizeunit,
                                'round2': 2},
           'E_DC_Min_img_plane': {'dtype': 'float64',
                                  'unit': sepsizeunit,
                                  'round2': 2},
           'DC_PA_img_plane': {'dtype': 'float64',
                               'unit': posunit,
                               'round2': 1},
           'E_DC_PA_img_plane': {'dtype': 'float64',
                                 'unit': posunit,
                                 'round2': 1},
           'Tile': {'dtype': '<U6',
                    'unit': None},
           'Subtile': {'dtype': '<U14',
                       'unit': None},
           'QL_image_RA': {'dtype': 'float64',
                           'unit': posunit,
                           'round2': 6},
           'QL_image_DEC': {'dtype': 'float64',
                            'unit': posunit,
                            'round2': 6},
           'NVSS_distance': {'dtype': 'float64',
                             'unit': sepsizeunit,
                             'round2': 2},
           'FIRST_distance': {'dtype': 'float64',
                              'unit': sepsizeunit,
                              'round2': 2},
           'Peak_to_ring': {'dtype': 'float64',
                            'unit': None,
                            'round2': 2},
           'Duplicate_flag': {'dtype': 'int64',
                              'unit': None},
           'Quality_flag': {'dtype': 'int64',
                            'unit': None},
           'Source_name': {'dtype': '<U19',
                           'unit': None},
           'Source_type': {'dtype': '<U2',
                           'unit': None},
           'NN_dist': {'dtype': 'float64',
                       'unit': sepsizeunit,
                       'round2': 2},
           'BMAJ': {'dtype': 'float64',
                    'unit': sepsizeunit,
                    'round2': 2},
           'BMIN': {'dtype': 'float64',
                    'unit': sepsizeunit,
                    'round2': 2},
           'BPA': {'dtype': 'float64',
                   'unit': posunit,
                   'round2': 2},
           'Best_neuron_y': {'dtype': 'int64',
                             'unit': None},
           'Best_neuron_x': {'dtype': 'int64',
                             'unit': None},
           'Neuron_dist': {'dtype': 'float64',
                           'unit': None,
                           'round2': 2},
           'P_sidelobe': {'dtype': 'float64',
                          'unit': None,
                          'round2': 1}}


###########################################################
###########################################################
###functionality

def tidy_table(data, coldict, maskfill=-99, sortcol='RA'):
    'cleans up data table and adds units'
    for col in data.colnames:
        if col in list(coldict.keys()):
            colmasked = False
            cdtype = coldict[col]['dtype']
            if type(data[col]) == MaskedColumn:
                colmasked = True
                cmask = data[col].mask
                data[col].fill_value = maskfill
                data[col] = data[col].filled()
            cdata = np.array(data[col]).astype(cdtype)
            if 'round2' in list(coldict[col].keys()):
                cdata = np.round(cdata, coldict[col]['round2'])
            data[col] = cdata
            data[col].unit = coldict[col]['unit']
            if colmasked == True:
                data[col] = MaskedColumn(data[col])
                data[col].mask = cmask
    
    data.sort(sortcol)
    
    return


def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="tidy CIRADA VLASS catalog")
    parser.add_argument("input_file",
                        help="CIRADA csv file of VLASS catalog")
    parser.add_argument("--outname", action='store',
                        type=str, default=None,
                        help="output filename")
    
    args = parser.parse_args()
    
    ###split input images into list
    if args.outname is None:
        args.outname = '.'.join([args.input_file.rsplit('.', 1)[0], 'fits'])
    
    return args

###########################################################
###########################################################
###main

if __name__ == '__main__':
    args = parse_args()
    data = Table.read(args.input_file)
    tidy_table(data=data, coldict=coldict)
    data.write(args.outname, format='fits')


