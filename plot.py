import matplotlib
matplotlib.use('Agg')
import h5py
import numpy
import matplotlib.pyplot as plt
import glob
import os
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 10
fig_size[1] = 13

os.makedirs('progress/vol2outputs')
os.makedirs('plots')

#all the names
f = open('progress/names.txt','r')
lines = f.readlines()
f.close()
names = []
for line in lines:
        names.append(line.split('\n')[0])

#going through ouputs adding the last fit of each object to the list
mags = {}
d = {}
calculatedNames = []
for i in range(9):
        tempNames = []
        outputList = glob.glob('progress/%soutputs/*_ML'%(11-i))
        for out in outputList:
                outName = out.split('/')[-1].split('_')[0].split('g')[1]
                if outName not in calculatedNames:
                        tempNames.append(outName)
                        calculatedNames.append(outName)
                        os.popen('cp progress/%soutputs/config%s_%s_ML progress/vol2outputs' %(11-i,outName,10-i) )
                        os.popen('cp progress/%soutputs/config%s_%s_ML progress/vol2outputs/config%s' %(11-i,outName,10-i,outName) )


        for tempName in tempNames:
                f = open('progress/%soutputs/config%s_%s_ML'%(11-i,tempName,10-i),'r')
                magLines = f.readlines()
                f.close()

                Gmag = magLines[41].split('\n')[0]
                Gmag = Gmag.split(' ')[1]
                Rmag = magLines[42].split('\n')[0]
                Rmag = Rmag.split(' ')[1]
                Imag = magLines[43].split('\n')[0]
                Imag = Imag.split(' ')[1]
                Zmag = magLines[44].split('\n')[0]
                Zmag = Zmag.split(' ')[1]
                Ymag = magLines[45].split('\n')[0]
                Ymag = Ymag.split(' ')[1]
                xai = magLines[-1].split('\n')[0]
                xai = xai.split(' ')[1]

                chain = h5py.File('progress/%soutputs/config%s_%s_chain.hdf5'%(11-i,tempName,10-i), 'r')
                Ger = chain['source1.mag_g'][:, -1].std()
                Rer = chain['source1.mag_r'][:, -1].std()
                Ier = chain['source1.mag_i'][:, -1].std()
                Zer = chain['source1.mag_z'][:, -1].std()
                Yer = chain['source1.mag_y'][:, -1].std()

                d['%s'%tempName] = '%s %s %s %s %s %s %s %s %s %s %s %s\n' %(tempName,Gmag,0.02,Rmag,0.03,Imag,0.04,Zmag,0.05,Ymag,0.06,xai)
                mags['%s'%tempName] = '%s %s %s %s %s %s %s %s %s %s %s %s\n'%(tempName,Gmag,Ger,Rmag,Rer,Imag,Ier,Zmag,Zer,Ymag,Yer,xai)
                calculatedNames.append('%s'%tempName)

#lastNames = []
#for name in names:
#if name not in calculatedNames:
#               lastNames.append(name)  

#for lastName in lastNames:
#       f = open('progress/4outputs/config%s_3_ML'%lastName)
#       magLines = f.readlines()
#       f.close()
#
#       Gmag = magLines[48].split('\n')[0]
#       Gmag = Gmag.split(' ')[1]
#       Rmag = magLines[49].split('\n')[0]
#       Rmag = Rmag.split(' ')[1]
#       Imag = magLines[50].split('\n')[0]
#       Imag = Imag.split(' ')[1]
#       Zmag = magLines[51].split('\n')[0]
#       Zmag = Zmag.split(' ')[1]
#       Ymag = magLines[52].split('\n')[0]
#       Ymag = Ymag.split(' ')[1]
#        xai = magLines[-1].split('\n')[0]
#        xai = xai.split(' ')[1]
#        
#        chain = h5py.File('progress/%soutputs/config%s_%s_chain.hdf5'%(10-i,tempName,9-i), 'r')
#        Ger = chain['light2.mag_g'][:, -1].std()
#        Rer = chain['light2.mag_r'][:, -1].std()
#        Ier = chain['light2.mag_i'][:, -1].std()
#        Zer = chain['light2.mag_z'][:, -1].std()
#        Yer = chain['light2.mag_y'][:, -1].std()
#       
#        mags['%s'%lastName] = '%s %s %s %s %s %s %s %s %s %s %s %s\n'%(lastName,Gmag,Ger,Rmag,Rer,Imag,Ier,Zmag,Zer,Ymag,Yer,xai)
#
#       calculatedNames.append('%s'%lastName)

newlines = []
newlines.append('ID Gmag Ger Rmag Rer Imag Ier Zmag Zer Ymag Yer Xai\n')
for name in names:
        #newlines.append(mags['%s'%name])
        newlines.append(d['%s'%name])

f = open('plots/mags.txt','w')
f.writelines(newlines)
f.close()

xai7500 = []
xai7500Id = []
xai25000 = []
xai25000Id = []
xaiRest = []
xaiRestId = []
for newline in newlines[1:]:
    xai = newline.split('\n')[0]
    xai = float(xai.split(' ')[-1])
    Id = newline.split(' ')[0]
    if xai >= -7500. :
        xai7500.append(newline)
        xai7500Id.append(Id)
    elif xai >= -25000. :
        xai25000.append(newline)
        xai25000Id.append(Id)
    else:
        xaiRest.append(newline)
        xaiRestId.append(Id)

realMagG = []
realMagR = []
realMagI = []
realMagZ = []
realMagY = []
specz = []
chain = h5py.File('data/sim0_parameters.hdf5','r')
for i in range(len(mags)):
    realMagG.append(chain['source_mag_g'][i])
    realMagR.append(chain['source_mag_r'][i])
    realMagI.append(chain['source_mag_i'][i])
    realMagZ.append(chain['source_mag_z'][i])
    realMagY.append(chain['source_mag_y'][i])
    specz.append(chain['z_source'][i])

i = 0
realMagLines =[]
for name in names:
    realMagLines.append('%s %s %s %s %s %s\n'%(name,realMagG[i],realMagR[i],realMagI[i],realMagZ[i],realMagY[i]))
    i += 1

magG = []
Ger = []
magR = []
Rer = []
magI = []
Ier = []
magZ = []
Zer = []
magY = []
Yer = []
Id = []
for newline in newlines[1:]:
    newline = newline.split('\n')[0]
    Id.append(newline.split(' ')[0])
    magG.append(float(newline.split(' ')[1]))
    magR.append(float(newline.split(' ')[3]))
    magI.append(float(newline.split(' ')[5]))
    magZ.append(float(newline.split(' ')[7]))
    magY.append(float(newline.split(' ')[9]))

    Ger.append(float(newline.split(' ')[2]))
    Rer.append(float(newline.split(' ')[4]))
    Ier.append(float(newline.split(' ')[6]))
    Zer.append(float(newline.split(' ')[8]))
    Yer.append(float(newline.split(' ')[10]))

plt.subplot(321)
plt.errorbar(realMagG,magG,yerr=Ger,marker='.',linestyle="None",markerfacecolor="red")
x = numpy.array(realMagG)
y = numpy.array(realMagG)
fit = numpy.polyfit(x,y,1)
fitFunction = numpy.poly1d(fit)
plt.plot(x,fitFunction(x),color="green")

plt.xlabel('$M_{g,true}$')
plt.ylabel('$M_{g,fit}$')
#plt.title('gBand Fitted Magnitude vs True Magnitude\nfit details:%1.3f   %1.3f'%(fit[0],fit[1]))
#plt.savefig('plots/gBandPlot.pdf')
#plt.close()

plt.subplot(322)
plt.errorbar(realMagR,magR,yerr=Rer,marker='.',linestyle="None",markerfacecolor="red")
x = numpy.array(realMagR)
y = numpy.array(realMagR)
fit = numpy.polyfit(x,y,1)
fitFunction = numpy.poly1d(fit)
plt.plot(x,fitFunction(x),color="green")
plt.xlabel('$M_{r,true}$')
plt.ylabel('$M_{r,fit}$')
#plt.savefig('plots/rBandPlot.pdf')
#plt.close()

plt.subplot(323)
plt.errorbar(realMagI,magI,yerr=Ier,marker='.',linestyle="None",markerfacecolor="red")
x = numpy.array(realMagI)
y = numpy.array(realMagI)
fit = numpy.polyfit(x,y,1)
fitFunction = numpy.poly1d(fit)
plt.plot(x,fitFunction(x),color="green")
plt.xlabel('$M_{i,true}$')
plt.ylabel('$M_{i,fit}$')
#plt.savefig('plots/iBandPlot.pdf')
#plt.close()

plt.subplot(324)
plt.errorbar(realMagZ,magZ,yerr=Zer,marker='.',linestyle="None",markerfacecolor="red")
x = numpy.array(realMagZ)
y = numpy.array(realMagZ)
fit = numpy.polyfit(x,y,1)
fitFunction = numpy.poly1d(fit)
plt.plot(x,fitFunction(x),color="green")
plt.xlabel('$M_{z,true}$')
plt.ylabel('$M_{z,fit}$')
#plt.savefig('plots/zBandPlot.pdf')
#plt.close()

plt.subplot(325)
plt.errorbar(realMagY,magY,yerr=Yer,marker='.',linestyle="None",markerfacecolor="red")
x = numpy.array(realMagY)
y = numpy.array(realMagY)
fit = numpy.polyfit(x,y,1)
fitFunction = numpy.poly1d(fit)
plt.plot(x,fitFunction(x),color="green")
plt.xlabel('$M_{y,true}$')
plt.ylabel('$M_{y,fit}$')
#plt.savefig('plots/yBandPlot.pdf')

realColor = []
simColor = []
for i in range(len(names)):
    simColor.append(float(magG[i]-magI[i]))
    realColor.append(float(realMagG[i]-realMagI[i]))

colorLines = []
colorLines.append('name simColor realColor colorError\n')
for i in range(len(names)):
    colorLine = '%s %s %s %s\n'%(names[i],simColor[i],realColor[i],abs(simColor[i]-realColor[i]))
    colorLines.append(colorLine)

f = open('plots/colors.txt','w')
f.writelines(colorLines)
f.close()

plt.subplot(326)
#plt.close()
plt.errorbar(realColor,simColor,yerr=Ier,marker='.',linestyle="None",markerfacecolor="red")
x = numpy.array(realColor)
y = numpy.array(realColor)
fit = numpy.polyfit(x,y,1)
fitFunction = numpy.poly1d(fit)
plt.plot(x,fitFunction(x),color="green")
plt.xlabel('True g-i colour')
plt.ylabel('Fit g-i colour')
plt.xlim(-.3,.5)
plt.ylim(-1.5,1.5)
#plt.title('FitColor vs TrueColor (g-i)\nfit details:%1.3f   %1.3f'%(fit[0],fit[1]))
#plt.savefig('plots/color.pdf')
#plt.savefig('plots/color.png')
plt.savefig('all.eps',format='eps',dpi=1000)

finalLines = []
finalDataLines = []
for i in range(len(realMagG)):
    line = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(names[i],realMagG[i],magG[i],Ger[i],realMagR[i],magR[i],Rer[i],realMagI[i],magI[i],Ier[i],realMagZ[i],magZ[i],Zer[i],realMagY[i],magY[i],Yer[i],realColor[i],simColor[i],Ier[i],specz[i])
    lineData = "%s %s %s %s %s %s %s %s %s %s %s %s\n" %(names[i],magG[i],Ger[i],magR[i],Rer[i],magI[i],Ier[i],magZ[i],Zer[i],magY[i],Yer[i],specz[i])
    finalLines.append(line)
    finalDataLines.append(lineData)

f = open('output.txt','w')
f.writelines(finalLines)
f.close()

f = open('data.txt','w')
f.writelines(finalDataLines)
f.close()

goodG = []
goodR = []
goodI = []
goodZ = []
goodY = []
goodGer = []
goodRer = []
goodIer = []
goodZer = []
goodYer = []
goodrealG = []
goodrealR = []
goodrealI = []
goodrealZ = []
goodrealY = []
for i in range(len(names)):
    if Id[i] in xai7500Id:
        goodI.append(magI[i])
        goodG.append(magG[i])
        goodIer.append(Ier[i])
        goodGer.append(Ger[i])
        goodrealG.append(realMagG[i])
        goodrealI.append(realMagI[i])

plt.close()
plt.errorbar(goodrealG,goodG,yerr=goodGer,marker='.',linestyle="None",markerfacecolor="red")
x = numpy.array(goodrealG)
y = numpy.array(goodrealG)
fit = numpy.polyfit(x,y,1)
fitFunction = numpy.poly1d(fit)
plt.plot(x,fitFunction(x),color="green")
plt.xlabel('gBand True Magnitude')
plt.ylabel('gBand Fit Magnitude')
plt.savefig('plots/gBandPlot(good).pdf')

plt.close()
plt.errorbar(goodrealI,goodI,yerr=goodIer,marker='.',linestyle="None",markerfacecolor="red")
x = numpy.array(goodrealI)
y = numpy.array(goodrealI)
fit = numpy.polyfit(x,y,1)
fitFunction = numpy.poly1d(fit)
plt.plot(x,fitFunction(x),color="green")
plt.xlabel('iBand True Magnitude')
plt.ylabel('iBand Fit Magnitude')
plt.savefig('plots/iBandPlot(good).pdf')

realColor = []
simColor = []
for i in range(len(xai7500Id)):
    simColor.append(float(goodG[i]-goodI[i]))
    realColor.append(float(goodrealG[i]-goodrealI[i]))

plt.close()
plt.errorbar(realColor,simColor,yerr=goodIer,marker='.',linestyle="None",markerfacecolor="red")
x = numpy.array(realColor)
y = numpy.array(realColor)
fit = numpy.polyfit(x,y,1)
fitFunction = numpy.poly1d(fit)
plt.plot(x,fitFunction(x),color="green")
plt.xlabel('True Color')
plt.ylabel('Fit Color')
plt.title('FitColor vs TrueColor (g-i)\nfit details:%1.3f   %1.3f'%(fit[0],fit[1]))
plt.savefig('plots/color(good).pdf')
