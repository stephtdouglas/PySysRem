import logging

import numpy as np
import astropy.io.fits as fits

class source:
    def __init__(self, mags, magerrors, epochs, flags, filename):
        self.mags = mags
        self.magerrors = magerrors
        self.epochs = epochs
        self.flags = flags
        self.filename = filename

        self.stats()

    def clean_up(self):
        # First, save the original arrays just in case
        self.orig_mags = np.copy(self.mags)
        self.orig_magerrors = np.copy(self.magerrors)
        self.orig_epochs = np.copy(self.epochs)
        self.orig_flags = np.copy(self.flags)

        # Remove anything with a flag > 0
        good = self.flags==0
        self.mags = self.mags[good]
        self.magerrors = self.magerrors[good]
        self.epochs = self.epochs[good]
        self.flags = self.flags[good]
        self.good_mask = good

        self.stats()

    def stats(self):
        self.median, self.std = self._calc_stats(self.mags)

    def _calc_stats(self,flux):
        return np.median(flux), np.std(flux)

    def remove_epoch(self, epoch):
        if epoch in self.epochs:
            epoch_n = np.where(self.epochs == epoch)[0][0]
            self.xschip = np.array(list(self.xschip[0:epoch_n]) + list(self.xschip[epoch_n+1:]),dtype=np.float32)
            self.yschip = np.array(list(self.yschip[0:epoch_n]) + list(self.yschip[epoch_n+1:]),dtype=np.float32)
            self.xs = np.array(list(self.xs[0:epoch_n]) + list(self.xs[epoch_n+1:]),dtype=np.float32)
            self.ys = np.array(list(self.ys[0:epoch_n]) + list(self.ys[epoch_n+1:]),dtype=np.float32)
            self.mags = np.array(list(self.mags[0:epoch_n]) + list(self.mags[epoch_n+1:]),dtype=np.float32)
            self.magerrors = np.array(list(self.magerrors[0:epoch_n]) + list(self.magerrors[epoch_n+1:]),dtype=np.float32)
            self.epochs = np.array(list(self.epochs[0:epoch_n]) + list(self.epochs[epoch_n+1:]),dtype=np.float32)
            self.flags = np.array(list(self.flags[0:epoch_n]) + list(self.flags[epoch_n+1:]),dtype=np.float32)
            self.images = np.array(list(self.images[0:epoch_n]) + list(self.images[epoch_n+1:]))
        else:
            print "*"*20, "Tried to remove epoch", epoch, "but not present in list", self.epochs

    @classmethod
    def from_ptf(cls,filename):
        lc = np.loadtxt(filename,usecols=np.arange(8))

        epochs = lc[:,0]
        mags = lc[:,1]
        magerrors = lc[:,2]
        xskys = lc[:,3]
        yskys = lc[:,4]
        xschip = lc[:,5]
        yschip = lc[:,6]
        flags = lc[:,7]
        xc = np.median(xskys)
        yc = np.median(yskys)
        # Right now it's not saving the image filenames, which should be fine

        ptf_lc = cls(mags, magerrors, epochs, flags, filename)
        return ptf_lc

    @classmethod
    def from_k2sff(cls,filename,ext=1):
        """ Create a light curve object from a K2SFF light curve.
        If no extension is specified, the "BESTAPER" extension (1) will be used.
        """

        if ext==0:
            logging.warning("Ext 0 is the header. Using BESTAPER (ext=1)")
            ext=1

        with fits.open(filename) as hdu:
            table = hdu[ext].data

            epochs = table["T"]
            mags = table["FCOR"]
            magerrors = np.ones_like(mags)
            xskys = np.zeros_like(mags)
            yskys = np.zeros_like(mags)
            xschip = np.zeros_like(mags)
            yschip = np.zeros_like(mags)
            flags = table["MOVING"]
            xc, yc = 0,0

        k2sff_lc = cls(mags, magerrors, epochs, flags, filename)
        return k2sff_lc

def fix_epochs(star_list):
    """ Correct a set of light curves so that all have the same time axis."""

    stars_dim = len(star_list)
    epoch_lengths = np.zeros(stars_dim,int)
    # Get all the epoch array lengths
    for i, star in enumerate(star_list):
        epoch_lengths[i] = len(star.epochs)

    if np.all(epoch_lengths==epoch_lengths[0]):
        print("No need for length correction.")
        return star_list

    # Find the longest
    longest = np.argmax(epoch_lengths)
    base_epochs = star_list[longest].epochs

    # I'm assuming the others are just missing points, and don't have extras.
    # It will get more complicated if we have to cross them all
    cadence = 0.02043
    half_cadence = cadence / 2.0
    for i, star in enumerate(star_list):
        len_diff = len(base_epochs) - len(star.epochs)
        if ((len_diff==0) and
            (np.all(np.abs(star.epochs-base_epochs)<half_cadence))):
            # If the arrays match, don't worry about it
            continue

        missing_epochs = np.zeros(len_diff)
        mcount = 0
        for i, ep in enumerate(base_epochs):
            ep_diff = np.abs(star.epochs - ep)
            if np.any(ep_diff<half_cadence):
                continue
            else:
                missing_epochs[mcount] = ep
                mcount += 1

        new_epochs = np.append(star.epochs,missing_epochs)
        extras = np.zeros_like(missing_epochs)
        new_flux = np.append(star.mags,extras)
        new_errors = np.append(star.magerrors,extras)
        new_flags = np.append(star.flags,np.ones(len(missing_epochs),int)*9999)

        sort_eps = np.argsort(new_epochs)
        star.epochs = new_epochs[sort_eps]
        star.mags = new_flux[sort_eps]
        star.magerrors = new_errors[sort_eps]
        star.flags = new_flags[sort_eps]

    return star_list
