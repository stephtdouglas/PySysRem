import numpy

class source:
    def __init__(self, xc, yc=0, xs=[], ys=[], xschip=[], yschip=[], mags=[], magerrors=[], epochs=[], flags=[], images=[], filename=""):
        if isinstance(xc, float) or isinstance(xc, int):
            self.xc = xc
            self.yc = yc
            self.xschip = xschip
            self.yschip = yschip
            self.xs = xs
            self.ys = ys
            self.mags = mags
            self.magerrors = magerrors
            self.epochs = epochs
            self.flags = flags
            self.images = images
            self.filename = filename
            
        if isinstance(xc, str): # a filename is specified
            epochs = []
            mags = []
            mag_errs = []
            xskys = []
            yskys = []
            xchips = []
            ychips = []
            flags = []
            images = []

            for l in open(xc,"r"):
                if xc.split('/')[-1].find("sysrem") == -1:
                  if len(l.strip()) > 0 and l[0] != '#':
                    s = l.split()
                    epochs.append(float(s[0]))
                    mags.append(float(s[1]))
                    mag_errs.append(float(s[2]))
                    xskys.append(float(s[3]))
                    yskys.append(float(s[4]))
                    xchips.append(float(s[5]))
                    ychips.append(float(s[6]))
                    flags.append(float(s[7]))
                    if len(s) == 9:
                        images.append(s[8])
                    else:
                        images.append("")
                else:
                  if len(l.strip()) > 0 and l[0] != '#':
                    s = l.split()
                    epochs.append(float(s[0]))
                    mags.append(float(s[2]))
                    mag_errs.append(float(s[3]))
                    xskys.append(float(s[4]))
                    yskys.append(float(s[5]))
                    xchips.append(float(s[6]))
                    ychips.append(float(s[7]))
                    flags.append(float(s[8]))
                    if len(s) == 10:
                        images.append(s[9])
                    else:
                        images.append("")
                    

            self.epochs = numpy.array(epochs)
            self.mags = numpy.array(mags)
            self.magerrors = numpy.array(mag_errs)
            self.xs = numpy.array(xskys)
            self.ys = numpy.array(yskys)
            self.xschip = numpy.array(xchips)
            self.yschip = numpy.array(ychips)
            self.flags = numpy.array(flags)
            self.xc = numpy.median(xskys)
            self.yc = numpy.median(yskys)
            self.images = numpy.array(images)
            self.filename = xc
        
    def good_epochs_only(self):
        f = self.flags
        return source(self.xc, self.yc, self.xs[f == 0], self.ys[f == 0], self.xschip[f == 0], self.yschip[f == 0], self.mags[f == 0], self.magerrors[f == 0], self.epochs[f == 0], self.flags[f == 0],self.images[f == 0], self.filename)

    def n_good_epochs(self):
        return len(self.flags[self.flags == 0])
    
    def remove_epoch(self, epoch):
        if epoch in self.epochs:
            epoch_n = numpy.where(self.epochs == epoch)[0][0]
            self.xschip = numpy.array(list(self.xschip[0:epoch_n]) + list(self.xschip[epoch_n+1:]),dtype=numpy.float32)
            self.yschip = numpy.array(list(self.yschip[0:epoch_n]) + list(self.yschip[epoch_n+1:]),dtype=numpy.float32)
            self.xs = numpy.array(list(self.xs[0:epoch_n]) + list(self.xs[epoch_n+1:]),dtype=numpy.float32)
            self.ys = numpy.array(list(self.ys[0:epoch_n]) + list(self.ys[epoch_n+1:]),dtype=numpy.float32)
            self.mags = numpy.array(list(self.mags[0:epoch_n]) + list(self.mags[epoch_n+1:]),dtype=numpy.float32)
            self.magerrors = numpy.array(list(self.magerrors[0:epoch_n]) + list(self.magerrors[epoch_n+1:]),dtype=numpy.float32)
            self.epochs = numpy.array(list(self.epochs[0:epoch_n]) + list(self.epochs[epoch_n+1:]),dtype=numpy.float32)
            self.flags = numpy.array(list(self.flags[0:epoch_n]) + list(self.flags[epoch_n+1:]),dtype=numpy.float32)
            self.images = numpy.array(list(self.images[0:epoch_n]) + list(self.images[epoch_n+1:]))
        else:
            print "*"*20, "Tried to remove epoch", epoch, "but not present in list", self.epochs
            
