"""
Given a science frame (e.g., SoCal), inject cosmic rays from the 30-minute
darks, and test astroscrappy recovery fraction.

NOTE: this assumes that the "cr_cache" has been created, by running
test_astroscrappy.py with do_darks = 1.
"""
import os, pickle
from glob import glob
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.io import fits
from kpf_cosmic_rays.paths import DATADIR, RESULTSDIR
from kpf_cosmic_rays.plotting import plot_dark, plot_kpf_L0
from kpf_cosmic_rays.helpers import read_fits

# https://astroscrappy.readthedocs.io/en/latest/api/astroscrappy.detect_cosmics.html
from astroscrappy import detect_cosmics

# https://github.com/lgbouma/aesthetic ; pip install aesthetic
from aesthetic.plot import savefig, format_ax, set_style
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

OUTDIR = os.path.join(RESULTSDIR, "cosmic_injection_recovery")
if not os.path.exists(OUTDIR): os.mkdir(OUTDIR)

verbose = True

def run_cosmic_injection_recovery(
    infitspath, plot_interim=1,
    cosmic_args={'sigclip':4, 'sigfrac':0.1, 'readnoise':4.0, 'bkg_err':5}
):
    """
    Given a L0 science frame, inject cosmics into each channel, and try to
    recover them.  Output result .CSV files to OUTDIR.

    Note: this procedure assumes the "background image" (cosmic free) is known!
    This assumption might benefit from some more serious poking.

    It also injects 30 minutes worth of cosmic rays into the target frame,
    regardless of the target frame exptime.

    Default arguments for CR detection:
        sigclip = 4 # passed to detect_cosmics
        sigfrac = 0.1  # passed to detect_cosmics
        bkg_err = 5 # noise injected into the image to make CR identification harder
        readnoise = 4.0 # e-; Ashley Baker, priv comm 2022/05/16
    """

    outname = (
        'cr-inj-recov_'+
        os.path.basename(infitspath).replace(".fits","")+
        f"_sigclip{cosmic_args['sigclip']}_sigfrac{cosmic_args['sigfrac']}"+
        f"_bkgerr{cosmic_args['bkg_err']}"
        '.csv'
    )
    outpath = os.path.join(OUTDIR, outname)
    if os.path.exists(outpath):
        print(f"Found {outpath}")
        return

    amp_names = [
        'GREEN_AMP1', 'GREEN_AMP2', 'GREEN_AMP3', 'GREEN_AMP4',
        'RED_AMP1', 'RED_AMP2'
    ]

    if verbose:
        print(infitspath)

    recov_fracs = {}

    if plot_interim:
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

        target_img, hdr = read_fits(infitspath, ext=amp_name)
        gain = hdr['CCDGAIN']
        hdul = fits.open(infitspath)
        target_exptime = hdul[0].header['EXPTIME']
        hdul.close()

        cachedir = os.path.join(DATADIR, 'cr_cache')
        cachepkls = glob(os.path.join(cachedir, "crcache*.pkl"))
        assert len(cachepkls) > 0

        np.random.seed(ix)
        cachepkl = np.random.choice(cachepkls)

        if verbose:
            print(f'{amp_name}: got {cachepkl}')

        with open(cachepkl, 'rb') as f:
            d = pickle.load(f)

        # metadata; currently not used
        cr_exptime = 1800 # The cosmic ray cache images had 1800 sec exptimes.
        exptime_mult_factor = cr_exptime / target_exptime

        # construct the cosmic ray image
        _img = d[f'{amp_name}_img']
        cr_mask = d[f'{amp_name}_crmask']
        cr_img = np.zeros_like(_img)
        cr_img[cr_mask] = _img[cr_mask]

        err_img = np.round(np.random.normal(
            loc=0, scale=cosmic_args['bkg_err'], size=_img.shape
        )).astype(int)

        # add the cosmic ray image to the target image
        img = target_img + cr_img

        # assume the background image (CR-free) has some error
        inbkg = target_img + err_img

        # detect cosmic rays in the resulting image 
        recov_mask, clean_img = detect_cosmics(
            img, inbkg=inbkg, gain=gain,
            readnoise=cosmic_args['readnoise'],
            sigclip=cosmic_args['sigclip'], sigfrac=cosmic_args['sigfrac']
        )

        # calculate the fraction of injected cosmic ray pixels that were
        # recovered as such.
        n_cr_pixels = np.sum(cr_mask)
        n_incorrect = np.sum(~(cr_mask == recov_mask))
        recov_frac = (n_cr_pixels - n_incorrect)/(n_cr_pixels)

        recov_fracs[amp_name+"_n_cr_pixels"] = n_cr_pixels
        recov_fracs[amp_name+"_n_incorrect"] = n_incorrect
        recov_fracs[amp_name+"_recov_frac"] = np.round(recov_frac, 5)

        if plot_interim:

            ax = axd[str(ix)]
            vmin, vmax = (
                int(np.nanpercentile(img, 1)), int(np.nanpercentile(img, 99))
            )
            norm = colors.Normalize(vmin=vmin, vmax=vmax)
            cset = ax.imshow(img.astype(int), cmap='binary_r', norm=norm,
                             origin='lower', interpolation='none')

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

    out_df = pd.DataFrame(recov_fracs, index=[0])
    out_df.to_csv(outpath, index=False)
    print(f"Wrote {outpath}")

    if plot_interim:
        fig.text(-0.01,0.5, 'CCD Column', va='center', rotation=90)
        fig.text(0.5,-0.01, 'CCD Row', ha='center')
        titlestr = outname.replace(".csv","") + f" exptime={target_exptime}"
        fig.suptitle(titlestr, y=1.0)
        fig.tight_layout()
        savefig(fig, outpath.replace(".csv",".png"), dpi=400, writepdf=False)


def main():
    # 20220519 SoCal frames
    fitsdir = os.path.join(DATADIR, '20220519')
    fitsnames = [
      #"KP.20220519.00090.23.fits",
      #"KP.20220519.00218.29.fits",
      #"KP.20220519.00346.34.fits",
      #"KP.20220519.00474.59.fits",
      "KP.20220519.00602.83.fits",
      "KP.20220519.00731.11.fits",
    ]
    socalpaths = [os.path.join(fitsdir, f) for f in fitsnames]

    for socalpath in socalpaths:

        for bkg_err in range(0,5):
            cosmic_args={'sigclip':4, 'sigfrac':0.1,
                         'readnoise':4.0, 'bkg_err':bkg_err}
            run_cosmic_injection_recovery(socalpath, cosmic_args=cosmic_args)

    csvpaths = glob(os.path.join(OUTDIR, "cr-inj-recov_KP*csv"))
    df = pd.concat([pd.read_csv(f) for f in csvpaths], ignore_index=True)
    df.index = [os.path.basename(c).replace("cr-inj-recov_","") for c in csvpaths]
    outcsv = os.path.join(OUTDIR, "cr-inj-recov-metadata.csv")
    df.T.to_csv(outcsv)
    print(f"Wrote {outcsv}")


if __name__ == "__main__":
    main()
