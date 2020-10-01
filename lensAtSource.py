import matplotlib
matplotlib.use('Agg')
from astropy.io import fits
import os
import glob
import numpy
import math
import matplotlib.pyplot as plot

os.makedirs('lensResAtSource')
os.makedirs('lensResAtSource/sources')
os.makedirs('lensResAtSource/res')

NBach = 1
directoryNumbers = []
directoryNames = []
for i in range(NBach):
    directoryNumbers.append(i)
    j = str(10*i+1).zfill(3)
    k = str(10*i+10).zfill(3)
    directoryNames.append('%s..%s'%(j,k))

idAddress = {}
yer = []
names = []
for i in directoryNumbers:
    f = open('41-1-%s/progress/names.txt'%directoryNames[directoryNumbers[i]])
    lines = f.readlines()
    f.close()

    for line in lines:
        idAddress['%s'%line.split('\n')[0]] = i
        names.append('%s'%line.split('\n')[0])

    f = open('41-1-%s/output.txt'%directoryNames[directoryNumbers[i]])
    lines = f.readlines()
    f.close()

    for line in lines:
        yer.append(float(line.split(' ')[6]))

filters = ['g','r','i','z','y']
ratio = {}
newlines = []
newlines.append('name filter resAtsourceSB sourceSB res/source I(Re)\n')
for name in names:
    #if simNumber[name] == 0:
    #    f = open('suicide0vol2/progress/vol2outputs/config%s'%name,'r')
    #    lines = f.readlines()
    #    f.close()

    #if simNumber[name] == 1:
    #    f = open('suicide1vol2/progress/vol2outputs/config%s'%name,'r')
    #    lines = f.readlines()
    #    f.close()

    #if simNumber[name] == 2:
    #    f = open('suicide2vol2/progress/vol2outputs/config%s'%name,'r')
    #    lines = f.readlines()
    #    f.close()

    f = open('41-1-%s/progress/vol2outputs/config%s'%(directoryNames[idAddress[name]],name),'r')
    lines = f.readlines()
    f.close()

    re = float(lines[38].split(' ')[1])
    n = float(lines[40].split(' ')[1])
    mag = {}
    mag['g'] = lines[41].split(' ')[1].split('\n')[0]
    mag['r'] = lines[42].split(' ')[1].split('\n')[0]
    mag['i'] = lines[43].split(' ')[1].split('\n')[0]
    mag['z'] = lines[44].split(' ')[1].split('\n')[0]
    mag['y'] = lines[45].split(' ')[1].split('\n')[0]

    #step0 = fits.open('progress/0outputs/config%s_ML.fits'%name)
    sourcePath = glob.glob('41-1-%s/progress/vol2outputs/config%s_*'%(directoryNames[idAddress[name]],name))
    sourcePath = sourcePath[0].split('_')[-2]
    source = fits.open('41-1-%s/progress/%soutputs/config%s_%s_ML.fits'%(directoryNames[idAddress[name]],int(sourcePath)+1,name,sourcePath))
    step0 = fits.open('41-1-%s/progress/%soutputs/config%s_%s_ML.fits'%(directoryNames[idAddress[name]],int(sourcePath)+1,name,sourcePath))
    count = 1
    for f in filters:
        #os.popen('cp data/data/sim0_%s_%s_var.fits ./residuals0/'%(name,f))
        #os.popen('cp data/data/sim0_%s_%s_psf.fits ./residuals0/'%(name,f))
        #original = fits.open('data/data/sim0_%s_%s_sci.fits'%(name,f))
        original = fits.open('../rangavar/13tousand/3quarry/data/%s_%s_sci.fits'%(name,f))
        step0[count].data = original[0].data - step0[count].data
        step0[count].writeto('./lensResAtSource/sim0_%s_%s_sci.fits'%(name,f))
        source[count+5].writeto('./lensResAtSource/sources/source_%s_%s.fits'%(name,f))
        permSource = source[count+5].data

        iRange,jRange = numpy.shape(permSource)

        iMax = 0.
        for i in range(iRange):
            for j in range(jRange):
                if permSource[i][j] >= iMax:
                    iMax = permSource[i][j]

        m = float(mag['%s'%f])
        #iRe = numpy.exp(-1.) / ( (10.**(m-27.))**(1./2.5) * n * math.gamma(2.*n) * 2. * numpy.pi * re**2. )

        #this definition of Re relies on a non-lensed source, buthere source sizes are different, making it illogical to use the Re for achieving iRe. We rather use the sky variance files instead
        #iRe = iMax * numpy.exp( - ( (iMax * 2. * numpy.pi * (re**2.) * n * math.gamma(2.*n) )/( 10.**((27.-m)/2.5) ) )**(1./(2.*n)) )
        #print(iRe)

        #using sky variance to get a threshold pixel value for defining the source borders
        #and we want the source region from the g-band since in the other bands, it might never go above sky background 
        if f == 'g':
            var = fits.open('41-1-%s/data/data/sim0_%s_%s_var.fits' %(directoryNames[idAddress[name]],name,f))[0].data.copy()

            #calculating the mimimum of the variance file 
            xRange,yRange = numpy.shape(var)
            minimum = 1000.0
            for x_pixel in range(xRange):
                for y_pixel in range(yRange):
                    if var[x_pixel][y_pixel] < minimum:
                        minimum = var[x_pixel][y_pixel]

            sigmaSky = minimum**0.5
            cutoff = 3.*sigmaSky
            iRe = cutoff

        #bn = 2.*n - 1./3. + 4./(405.*n) + 46./(25515.*(n**2.)) + 131./(1148175.*(n**3.)) - 2194697./(30690717750.*(m**4.))
        #iRe = ( (10.**((27.-m)/2.5)) * (bn**(2.*n)) )/(numpy.exp(bn)*2.*numpy.pi*(re**2.)*n*math.gamma(2.*n))
        #print(iRe)

        resAtSource = numpy.zeros((iRange,jRange))
        sourceSurfaceB = 0
        resSurfaceB = 0
        for i in range(iRange):
            for j in range(jRange):
                if permSource[i][j] >= iRe:
                    resAtSource[i][j] = step0[count].data[i][j]
                    sourceSurfaceB = sourceSurfaceB + permSource[i][j]
                    resSurfaceB = resSurfaceB + resAtSource[i][j]

        step0[count].data = resAtSource
        step0[count].writeto('./lensResAtSource/res/res_%s_%s_sci.fits'%(name,f))

        print(name,resSurfaceB,sourceSurfaceB)
        newlines.append('%s %s %s %s %s %s\n' %(name,f,resSurfaceB,sourceSurfaceB,resSurfaceB/sourceSurfaceB,iRe) )


        ratio['%s.%s'%(name,f)] = resSurfaceB/sourceSurfaceB

        count = count + 1

f = open('lensResAtSource/resRatio.txt','w')
f.writelines(newlines)
f.close()

clines = []
colorError = {}
for i in range(NBach):
    f = open('41-1-%s/plots/colors.txt'%directoryNames[directoryNumbers[i]],'r')
    lines = f.readlines()
    f.close()

    for line in lines[1:]:
        clines.append(line)
        name = line.split(' ')[0]
        if name in names:
            error = line.split(' ')[3].split('\n')[0]
            colorError[name] = float(error)

gRatio = []
iRatio = []
cError = []
for name in names:
    gRatio.append(float(ratio['%s.%s'%(name,'g')]))
    iRatio.append(float(ratio['%s.%s'%(name,'i')]))
    cError.append(colorError[name])

#giLogRatii = []
#for i in range(20):
#    x = numpy.log10(abs(gRatio[i])*abs(iRatio[i]))
#    giLogRatio.append(x)

errorDependency = []
for i in range(len(names)):
    x = numpy.log10( (1.+iRatio[i])/(1.+gRatio[i]) )
    errorDependency.append(x)

plot.scatter(errorDependency,cError)
plot.ylabel('Error on G-I')
plot.xlabel('Error Depenedency on the Resiual (if it exists)')
plot.savefig('lensResAtSource/colorErrorVsResidual.pdf')
plot.close()

simColor = {}
trueColor = {}
colorError = {}
for line in clines:
    name = line.split(' ')[0]
    sim = float(line.split(' ')[1])
    true = float(line.split(' ')[2])
    error = true - sim
    simColor[name] = sim
    trueColor[name] = true
    colorError[name] = error

cError = []
for name in names:
    cError.append(colorError[name])

#plot.scatter(errorDependency,cError)
#plot.ylabel('Error on G minus I(color)')
#plot.xlabel('Error Depenedency on the Resiual (x)')
#plot.savefig('colorErrorVsResidualRealValue.pdf')
#plot.savefig('colorErrorVsResidualRealValue.png')
#plot.close()

#plot.scatter(errorDependency,cError)
#plot.ylabel('Error on G minus I(color)')
#plot.ylim(-.3,.3)
#plot.xlabel('Residuals Ratio error dependency(x)')
#plot.savefig('colorErrorVsResidualRealValueLimithed.pdf')
#plot.savefig('colorErrorVsResidualRealValueLimithed.png')
#plot.close()

plot.errorbar(errorDependency,cError,yerr=yer,marker='.',linestyle="None",markerfacecolor="red")
plot.ylim(-1.,1.)
plot.xlabel(r'$\tau$',fontsize=14)
plot.ylabel(r'$\delta C$',fontsize=14)
plot.savefig('lensResAtSource/residuals.eps',format='eps',dpi=1000)
