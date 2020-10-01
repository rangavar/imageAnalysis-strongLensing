import matplotlib
matplotlib.use('Agg')
from astropy.io import fits
import h5py
import numpy
import os
import matplotlib.pyplot as plt

os.makedirs('sources')

NBach = 1
directoryNumbers = []
directoryNames = []
for i in range(NBach):
    directoryNumbers.append(i)
    j = str(10*i+1).zfill(3)
    k = str(10*i+10).zfill(3)
    directoryNames.append('%s..%s'%(j,k))

idAddress = {}
names = []
error = []
bad = []
for i in directoryNumbers:
    f = open('../../../artemis/41-1-%s/progress/names.txt'%directoryNames[directoryNumbers[i]])
    lines = f.readlines()
    f.close()

    for line in lines:
        idAddress['%s'%line.split('\n')[0]] = i
        names.append('%s'%line.split('\n')[0])

    f = open('../../../artemis/41-1-%s/plots/colors.txt'%directoryNames[directoryNumbers[i]],'r')
    lines = f.readlines()
    f.close()

    for line in lines[1:]:
        error.append(float(line.split(' ')[-1].split('\n')[0]))

    f = open('../../../artemis/41-1-%s/bad.txt'%directoryNames[directoryNumbers[i]],'r')
    lines = f.readlines()
    f.close()

    for line in lines:
        line = line.split(' - ')[0]
        bad.append(line)

#this doesn't seem to be the right place for the source centroid position (rather we want to use the pixel with maximum pixle value as the centroid
x = []
y = []
chain = h5py.File('../41makeSims1/0simulations/filter/sim0_parameters.hdf5', 'r')
for i in range(len(names)):
    x.append(chain['x_source'][i])
    y.append(chain['y_source'][i])

#f = open('../../suicide.badmasks.may6/suicide3/colors.txt','r')
#lines = f.readlines()
#f.close()
#chain = h5py.File('../4makeSims/3simulations/sm0.hdf5', 'r')
#for i in range(17):
#    names.append(chain['object_ids'][i])
#    x.append(chain['x_source'][i])
#    y.append(chain['y_source'][i])
#    error.append(float(lines[i+1].split(' ')[-1].split('\n')[0]))

ratio = numpy.empty(len(names))
lines = []
lines.append('name ratio error\n')
for i in range(len(names)):
    noSource = fits.open('../3quarry/data/%s_i_sci.fits'%names[i])
    #this definition is wrong and will be changed shortly
    lens = float(noSource[0].data[int(y[i]),int(x[i])])
    imageFile = fits.open('../41makeSims1/0simulations/data/sim0_%s_i_sci.fits'%names[i])
    image = imageFile[0].data
    noSource[0].data = image - noSource[0].data

    #finding the source maximum pixel
    xRange,yRange = numpy.shape(noSource[0].data)
    iMax = 0.
    for x_pixel in range(xRange):
        for y_pixel in range(yRange):
            if noSource[0].data[x_pixel][y_pixel] >= iMax:
                iMax = noSource[0].data[x_pixel][y_pixel]
                sourcex = x_pixel
                sourcey = y_pixel

    #lens again
    newnoSource = fits.open('../3quarry/data/%s_i_sci.fits'%names[i])
    lens = float(newnoSource[0].data[sourcex,sourcey])
    source = float(noSource[0].data[sourcex,sourcey])
    noSource.writeto('sources/%s_source.fits'%names[i])
    ratio[i] = float(lens/source)
    lines.append('%s %f %f\n'%(names[i],ratio[i],error[i]))

#for i in range(60,77):
#    noSource = fits.open('../3quarry/data/%s_i_sci.fits'%names[i])
#    lens = float(noSource[0].data[int(y[i]),int(x[i])])
#    imageFile = fits.open('../4makeSims/3simulations/data0/sim0_%s_i_sci.fits'%names[i])
#    image = imageFile[0].data
#    noSource[0].data = image - noSource[0].data
#    #source = float(numpy.max(noSource[0].data))
#    source = float(noSource[0].data[int(y[i]),int(x[i])])
#    noSource.writeto('sources/%s_source.fits'%names[i])
#    ratio[i] = float(lens/source)
#    lines.append('%s %f %f\n'%(names[i],ratio[i],error[i]))

f = open('ratio.txt','w')
f.writelines(lines)
f.close()

#bad = [12,15,26,33,43,44,46,58]

#c = 0
#d = 0
#for i in range(60):
#    if i in bad:
#        if c == 0:
#            plt.scatter(ratio[i],error[i],c='r',s=3 , label = 'third light model is needed')
#            c = c+1
#        else:
#            plt.scatter(ratio[i],error[i],c='r',s=3)

#    else:
#        if d == 0:
#            plt.scatter(ratio[i],error[i],c='b',s=3 , label = 'no need for third light model')
#            d = d + 1
#        else:
#            plt.scatter(ratio[i],error[i],c='b',s=3)

newRatio = ratio
newError = error
for i in range(len(names)):
    if names[i] in bad:
        for j in range(len(newRatio)):
            if newRatio[j] == ratio[i]:
                ratioKey = j
        for k in range(len(newError)):
            if newError[k] == error[i]:
                errorKey = k
        newRatio = numpy.delete(newRatio,ratioKey)
        newError = numpy.delete(newError,errorKey)

plt.scatter(newRatio,newError,c='b',s=3)
plt.xscale('log')
plt.yscale('log')
plt.axhline(0.2,linestyle = '--')
plt.axvline(2.0,linestyle = ':')

plt.xlabel('surface brightness ratio of lens to source,  $\lambda_{i}$',fontsize=13)
plt.ylabel('color measurement error, $\Delta C$',fontsize=13)
#plt.set_xscale('log')
#plt.set_yscale('log')
#plt.legend(loc='upper left')
plt.savefig('ratio.eps',format='eps', dpi=1000)

#newRatio = []
#newError = []
#for i in range(20):
#    newRatio.append(ratio[i])
#    newError.append(error[i])

#plt.close()
#plt.scatter(newRatio,newError,s=5)
#plt.xlabel('lens to source ratio at the source position')
#plt.ylabel('error')
#plt.title('0Sim')
#plt.savefig('0ratio_maskingModified.pdf',dpi=300)

#newRatio = []
#newError = []
#for i in range(20,40):
#    newRatio.append(ratio[i])
#    newError.append(error[i])

#plt.close()
#plt.scatter(newRatio,newError,s=5)
#plt.xlabel('lens to source ratio at the source position')
#plt.ylabel('error')
#plt.title('1Sim')
#plt.ylim(0,0.8)
#plt.savefig('1ratio_maskingModified.pdf',dpi=300)

#newRatio = []
#newError = []
#for i in range(40,60):
#    newRatio.append(ratio[i])
#    newError.append(error[i])

#plt.close()
#plt.scatter(newRatio,newError,s=5)
#plt.xlabel('lens to source ratio at the source position')
#plt.ylabel('error')
#plt.title('2Sim')
#plt.savefig('2ratio_maskingModified.pdf',dpi=300)

#newRatio = []
#newError = []
#for i in range(60,77):
#    newRatio.append(ratio[i])
#    newError.append(error[i])

#plt.close()
#plt.scatter(newRatio,newError,s=5)
#plt.xlabel('lens to source ratio at the source position')
#plt.ylabel('error')
#plt.title('3Sim')
#plt.savefig('3ratio_maskingModified.pdf',dpi=300)

#newRatio = []
#newError = []
#for i in range(60):
#    newRatio.append(ratio[i])
#    newError.append(error[i])

#plt.close()
#plt.scatter(newRatio,newError,s=5)
#plt.xlabel('ratio of lens to source surface brightness at source position')
#plt.ylabel('color measurement error')
#plt.title('012sims')
#plt.savefig('012ratio_maskingModified.pdf',dpi=300)
