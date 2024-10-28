###command line script to stack VLASS images based on Michael's jupyter notebook

###point to multiple images separated by comma to stack
##e.g. >python vlass_image_stacking.py vlass_ep1.fits,vlass_ep2.fits,vlass_ep3.fits

import numpy as np, argparse
from astropy.io import fits
from reproject.mosaicking import *
from reproject import reproject_exact
from astropy.wcs import WCS
from astropy.convolution import Gaussian2DKernel, convolve_fft
import radio_beam as rb
from spectral_cube import SpectralCube
from warnings import filterwarnings
filterwarnings('ignore')


############################################################################
############################################################################
###parameters

############################################################################
############################################################################
###functions
    
def make_vlass_coad(input_files, outname='vlass_coad.fits',
                    overwrite_old=False):
    'make coadded VLASS image from inputs'
    
    ###load, reduce to 2D and extract necessary info
    cubes = [SpectralCube.read(i) for i in input_files]
    imdat = [(i.hdu.data[0], i.wcs.dropaxis(2)) for i in cubes]
    wcs_out, shape_out = find_optimal_celestial_wcs(imdat)
    
    ###find common beam and convolve images with
    beams = rb.Beams(beams = [rb.Beam.from_fits_header(i.header) for i in cubes])
    common = rb.commonbeam.common_manybeams_mve(beams)
    convolved = [i.convolve_to(common, convolve=convolve_fft) for i in cubes]
    indata = [(i.hdu.data[0], i.wcs.dropaxis(2)) for i in convolved]
    
    ###coad convolved images
    outdata, footprint = reproject_and_coadd(indata, wcs_out,
                                             shape_out=shape_out,
                                             reproject_function=reproject_exact)
    
    ###make output fits hdu and write to file
    hdu_out = fits.PrimaryHDU(outdata, header=wcs_out.to_header())
    hdu_out.header = common.attach_to_header(hdu_out.header) ##adds beam info
    ###add additional info to header (obs dates [median/range], comment on input images)?
    hdu_out = fits.HDUList(hdu_out)
    hdu_out.writeto(outname, overwrite=overwrite_old)
    
    return


def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="coad VLASS images")
    parser.add_argument("input_files",
                        help="image files to coad")
    parser.add_argument("--indir", action='store',
                        type=str, default='.',
                        help="input directory")
    parser.add_argument("--outdir", action='store',
                        type=str, default='.',
                        help="output directory")
    parser.add_argument("--outname", action='store',
                        type=str, default='vlass_coad.fits',
                        help="output filename")
    
    args = parser.parse_args()
    
    ###split input images into list
    args.input_files = args.input_files.split(',')
    args.input_files = ['/'.join([args.indir.removesuffix('/'), i])
                        for i in args.input_files]
    args.outname = '/'.join([args.outdir.removesuffix('/'), args.outname])
    
    return args

############################################################################
############################################################################
###main

if __name__ == '__main__':
    args = parse_args()
    make_vlass_coad(input_files=args.input_files,
                    outname=args.outname)



