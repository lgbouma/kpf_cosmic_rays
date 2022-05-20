"""
plot_dark
plot_kpf_L0
"""
import os, pickle, re, subprocess, itertools
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter

from datetime import datetime
from numpy import array as nparr
from astropy.io import fits
from itertools import product
from matplotlib import patches

from kpf_cosmic_rays.helpers import read_fits
from kpf_cosmic_rays.paths import RESULTSDIR

from aesthetic.plot import savefig, format_ax, set_style


def plot_dark(darkpath, biaspath=None, get_cosmics=0):

    outdir = os.path.join(RESULTSDIR, 'darks')
    outname = os.path.basename(darkpath).replace('.fits','')
    s = ''
    if biaspath is not None:
        s+='_L0-Bias'
    else:
        s+='_L0'
    if get_cosmics:
        s+='_cosmicmask'
    outpath = os.path.join(outdir, f'{outname}{s}.png')

    amp_names = [
        'GREEN_AMP1', 'GREEN_AMP2', 'GREEN_AMP3', 'GREEN_AMP4',
        'RED_AMP1', 'RED_AMP2'
    ]

    plt.close('all')

    fig = plt.figure(figsize=(9,16))
    axd = fig.subplot_mosaic(
        """
        01
        23
        45
        45
        """
    )

    for ix, amp_name in enumerate(amp_names):

        ax = axd[str(ix)]

        img, hdr = read_fits(darkpath, ext=amp_name)
        gain = hdr['CCDGAIN']

        # Subtract out the bias image to help mitigate fixed-pattern noise.
        # Include a "zero-point" to avoid negative values confusing the
        # laplacian.
        if biaspath is not None:
            bias_img, bias_hdr = read_fits(biaspath, ext=amp_name)
            orig_1pct = int(np.ceil(np.nanpercentile(img, 1)))
            img -= bias_img
            img += orig_1pct

        if get_cosmics:
            # https://astroscrappy.readthedocs.io/en/latest/api/astroscrappy.detect_cosmics.html
            from astroscrappy import detect_cosmics
            readnoise = 4.0 # e-; Ashley Baker, priv comm 2022/05/16

            if biaspath is not None:
                # adding bias subtraction for darks introduces major changes to
                # the laplacian
                sigclip = 4
                sigfrac = 0.3
            else:
                sigclip = 4
                sigfrac = 0.1

            crmask, cleanarr = detect_cosmics(
                img, gain=gain, readnoise=readnoise,
                sigclip=sigclip, sigfrac=sigfrac
            )

        if not get_cosmics:
            vmin, vmax = (
                int(np.nanpercentile(img, 1)), int(np.nanpercentile(img, 99))
            )
            norm = colors.Normalize(vmin=vmin, vmax=vmax)
            cset = ax.imshow(img.astype(int), cmap='binary_r', norm=norm, origin='lower',
                             interpolation='none')
        elif get_cosmics:
            vmin, vmax = 0, 1
            norm = colors.Normalize(vmin=vmin, vmax=vmax)
            cset = ax.imshow(crmask.astype(int), cmap='binary', origin='lower',
                            interpolation='none')

            assert np.all(np.unique(crmask.astype(int)) == np.array([0,1]))

        ax.set_title(amp_name)

        # colorbars
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        cb = fig.colorbar(cset, ax=ax, cax=cax, extend='both')

        ticks = np.arange(int(vmin), int(vmax)+1, 1)
        cb.set_ticks(ticks)
        cb.set_ticklabels(np.array(ticks).astype(str))

        cb.ax.tick_params(direction='in')
        cb.ax.tick_params(labelsize='xx-small')

    fig.text(-0.01,0.5, 'CCD Column', va='center', rotation=90)
    fig.text(0.5,-0.01, 'CCD Row', ha='center')

    titlestr = outname
    if biaspath is None:
        titlestr += " (L0 dark), "
    else:
        titlestr += " (L0 dark - bias + ZP), "
    if not get_cosmics:
        titlestr += "1st to 99th pctile linear stretch"
    else:
        titlestr += f"cosmics boolean mask (σ={sigclip}, f={sigfrac})"

    fig.suptitle(titlestr, y=1.0)

    fig.tight_layout()

    savefig(fig, outpath, dpi=400, writepdf=False)


def plot_kpf_L0(fitspath, outdir=None, typestr='', overwrite=0):

    s = '' if outdir is None else "_"+os.path.basename(outdir)
    if outdir is None: outdir = os.path.join(RESULTSDIR, 'L0')
    if not os.path.exists(outdir): os.mkdir(outdir)
    outname = os.path.basename(fitspath).replace('.fits','')
    outpath = os.path.join(outdir, f'{outname}{s}.png')

    if not overwrite and os.path.exists(outpath):
        print(f"Found {outpath}. Continue.")
        return 1

    amp_names = [
        'GREEN_AMP1', 'GREEN_AMP2', 'GREEN_AMP3', 'GREEN_AMP4',
        'RED_AMP1', 'RED_AMP2'
    ]

    plt.close('all')

    fig = plt.figure(figsize=(9,16))
    axd = fig.subplot_mosaic(
        """
        01
        23
        45
        45
        """
    )

    for ix, amp_name in enumerate(amp_names):

        ax = axd[str(ix)]

        img, hdr = read_fits(fitspath, ext=amp_name)

        vmin, vmax = (
            int(np.nanpercentile(img, 1)), int(np.nanpercentile(img, 99))
        )
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cset = ax.imshow(img.astype(int), cmap='binary_r', norm=norm, origin='lower',
                         interpolation='none')

        ax.set_title(amp_name)

        # colorbars
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        cb = fig.colorbar(cset, ax=ax, cax=cax, extend='both')

        if vmax - vmin < 20:
            dtick = 1
        elif vmax - vmin < 100:
            dtick = 10
        elif vmax - vmin < 1000:
            dtick = 50
        elif vmax - vmin < 3000:
            dtick = 100
        elif vmax - vmin < 10000:
            dtick = 500
        ticks = np.arange(int(vmin), int(vmax)+dtick, dtick)
        cb.set_ticks(ticks)
        cb.set_ticklabels(np.array(ticks).astype(str))

        cb.ax.tick_params(direction='in')
        cb.ax.tick_params(labelsize='xx-small')

    fig.text(-0.01,0.5, 'CCD Column', va='center', rotation=90)
    fig.text(0.5,-0.01, 'CCD Row', ha='center')

    fig.suptitle(outname + f' (L0,{typestr}) , 1st to 99th pctile linear stretch', y=1.0)

    fig.tight_layout()

    savefig(fig, outpath, dpi=400, writepdf=False)
