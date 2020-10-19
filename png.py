from astropy.io import fits as pyfits
from PIL import Image, ImageDraw
import numpy as np
import os
import sys
import pyplz_rgbtools


# Marshall, P. et al. (2016) "Space Warps I"

f = open('0simulations/filter/sim0_lens.cat', 'r')
lines = f.readlines()
f.close()

ncol = 5
nlens = len(lines)
nrow = nlens/ncol
if nlens%ncol > 0:
    nrow += 1

rgbbands = ('i', 'r', 'g')
datadir = '0simulations/data/'

# marshall16 parameters
scales = [0.4, 0.6, 1.8]
alpha = 1.
Q = 1.

def make_rgb(rootname, rgbbands, datadir='/', outputdir='./'):

    rname = datadir + 'sim0_%s_%s_sci.fits'%(rootname, rgbbands[0])

    ny, nx = pyfits.open(rname)[0].data.shape

    im = Image.new('RGB', (nx, ny), 'black')

    data = []

    for i in range(3):

        band = rgbbands[i]
        dname = datadir + 'sim0_%s_%s_sci.fits'%(rootname, band)
        d = pyfits.open(dname)[0].data.copy()
        data.append(d)

    im.putdata(pyplz_rgbtools.marshall16_pil_format(data, scales=scales, alpha=alpha, Q=Q))

    return im

npix = 101

fullimg = Image.new('RGB', (ncol*npix, nrow*npix), 'black')

for i in range(nlens):

    name = lines[i+24].split()[0]

    colno = i%ncol
    rowno = i/ncol

    im = make_rgb(name, rgbbands, datadir=datadir)

    fullimg.paste(im, (npix*colno, npix*rowno))

fullimg.save('0simulations/sample_collage_%dx%d.png'%(nrow, ncol))
