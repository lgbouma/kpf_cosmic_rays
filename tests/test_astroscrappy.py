import os
from glob import glob
from astropy.io import fits
from kpf_cosmic_rays.paths import DATADIR, RESULTSDIR
from kpf_cosmic_rays.plotting import plot_dark, plot_kpf_L0

include_bias = 0
do_darks = 1
do_SoCal = 0
do_LFC = 0

# darks from 20220518
fitsdir = os.path.join(DATADIR, '20220518')
fitsnames = [
    "KP.20220518.06621.03.fits",
    "KP.20220518.08534.24.fits",
    "KP.20220518.10447.44.fits",
    "KP.20220518.12360.74.fits",
]
darkpaths = [os.path.join(fitsdir, f) for f in fitsnames]

# 20220519 SoCal frames
fitsdir = os.path.join(DATADIR, '20220519')
fitsnames = [
  "KP.20220519.00090.23.fits",
  "KP.20220519.00218.29.fits",
  "KP.20220519.00346.34.fits",
  "KP.20220519.00474.59.fits",
  "KP.20220519.00602.83.fits",
  "KP.20220519.00731.11.fits",
]
socalpaths = [os.path.join(fitsdir, f) for f in fitsnames]
if do_SoCal:
    for socalpath in socalpaths:
        socaldir = os.path.join(RESULTSDIR, 'SoCal')
        plot_kpf_L0(socalpath, outdir=socaldir, typestr='SoCal')

# 20220519 LFC frames
fitsdir = os.path.join(DATADIR, '20220519')
fitsnames = [
  "KP.20220519.00849.87.fits",
  "KP.20220519.00967.90.fits",
  "KP.20220519.01085.93.fits",
]
lfcpaths = [os.path.join(fitsdir, f) for f in fitsnames]
if do_LFC:
    for lfcpath in lfcpaths:
        lfcdir = os.path.join(RESULTSDIR, 'LFC')
        plot_kpf_L0(lfcpath, outdir=lfcdir, typestr='LFC')


if include_bias:
    fitsdir = os.path.join(DATADIR, '20220518')
    biaspath = os.path.join(fitsdir, "KP.20220518.04508.75.fits")
    biasdir = os.path.join(RESULTSDIR, 'bias')
    plot_kpf_L0(biaspath, outdir=biasdir, typestr='bias')
else:
    biaspath = None

if do_darks:
    for darkpath in darkpaths:
        if include_bias:
            plot_dark(darkpath, get_cosmics=0, biaspath=biaspath)
            plot_dark(darkpath, get_cosmics=1, biaspath=biaspath)
        plot_dark(darkpath, get_cosmics=0)
        plot_dark(darkpath, get_cosmics=1)
