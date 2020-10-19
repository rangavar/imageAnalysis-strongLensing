from astropy.io import fits as pyfits
from PIL import Image, ImageDraw
import numpy as np
import os
import sys
import pyplz_rgbtools
import h5py

f = open('0simulations/filter/sim0_lens.cat', 'r')
lines = f.readlines()
f.close()
lines = lines[1:]

data = h5py.File('0simulations/filter/sim0_parameters.hdf5', 'r')
i_mags = {}
count = 0
for line in lines:
    ID = line.split(' ')[0]
    i_mags[ID] = data['source_mag_i'][count]
    count = count + 1

temp_sorted_list = sorted(i_mags, key=i_mags.get)
delete_list = [0,15,21,26,27,28,29,30,31,32,41,42,43,45,47,48,49,53,54,55,59,61,64,66]
#delete_list = []
delete_id = []
for key in delete_list:
    delete_id.append(temp_sorted_list[key])
clean_sorted_list = []
for key in temp_sorted_list:
    if key not in delete_id:
        clean_sorted_list.append(key)

sorted_list = {}
for name in clean_sorted_list:
    sorted_list[name] = i_mags[name]
print(len(sorted_list))

ncol = 5
nlens = len(sorted_list)
nrow = nlens/ncol
nrow = int(nrow)
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

count = 0
for name in sorted_list.keys():
    print(name,sorted_list[name])

    colno = int(count%ncol)
    rowno = int(count/ncol)

    im = make_rgb(name, rgbbands, datadir=datadir)

    fullimg.paste(im, (npix*colno, npix*rowno))

    count = count + 1

fullimg.save('0simulations/sample_collage.png')
