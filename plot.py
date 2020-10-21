import matplotlib
matplotlib.use('Agg')
import h5py
import numpy
import matplotlib.pyplot as plt
import glob
import os
from astropy.io import fits

os.makedirs('41-1-011..020/progress/vol2outputs')
os.makedirs('41-1-031..040/progress/vol2outputs')
directories = ['41-1-011..020','41-1-031..040']

names = []
names_address = {}
calculatedNames = []
mags = {}
d = {}
for directory in directories:
        #all the names
        f = open('%s/progress/names.txt'%directory,'r')
        lines = f.readlines()
        f.close()
        names_address[directory] = []
        for line in lines:
                names.append(line.split('\n')[0])
                names_address[directory].append(line.split('\n')[0])
            
        #going through ouputs adding the last fit of each object to the list
        for i in range(9):
                tempNames = []
                outputList = glob.glob('%s/progress/%soutputs/*_ML'%(directory,11-i))
                for out in outputList:
                        outName = out.split('/')[-1].split('_')[0].split('g')[1]
                        if outName not in calculatedNames:
                                tempNames.append(outName)
                                calculatedNames.append(outName)
                                os.popen('cp %s/progress/%soutputs/config%s_%s_ML %s/progress/vol2outputs' %(directory,11-i,outName,10-i,directory) )
                                os.popen('cp %s/progress/%soutputs/config%s_%s_ML %s/progress/vol2outputs/config%s' %(directory,11-i,outName,10-i,directory,outName) )


                for tempName in tempNames:
                        f = open('%s/progress/%soutputs/config%s_%s_ML'%(directory,11-i,tempName,10-i),'r')
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

                        chain = h5py.File('%s/progress/%soutputs/config%s_%s_chain.hdf5'%(directory,11-i,tempName,10-i), 'r')
                        Ger = chain['source1.mag_g'][:, -1].std()
                        Rer = chain['source1.mag_r'][:, -1].std()
                        Ier = chain['source1.mag_i'][:, -1].std()
                        Zer = chain['source1.mag_z'][:, -1].std()
                        Yer = chain['source1.mag_y'][:, -1].std()

                        d['%s'%tempName] = '%s %s %s %s %s %s %s %s %s %s %s %s\n' %(tempName,Gmag,0.02,Rmag,0.03,Imag,0.04,Zmag,0.05,Ymag,0.06,xai)
                        mags['%s'%tempName] = '%s %s %s %s %s %s %s %s %s %s %s %s\n'%(tempName,Gmag,Ger,Rmag,Rer,Imag,Ier,Zmag,Zer,Ymag,Yer,xai)
                        calculatedNames.append('%s'%tempName)

#saving the read magnitudes anad errors in a .txt file
newlines = []
newlines.append('ID Gmag Ger Rmag Rer Imag Ier Zmag Zer Ymag Yer Xai\n')
for name in names:
        newlines.append(mags['%s'%name])
        #newlines.append(d['%s'%name])

f = open('paper_plots/mags.txt','w')
f.writelines(newlines)
f.close()

#reading the real magnitudes and specz
real_data = {}
for directory in directories:
    chain = h5py.File('%s/data/sim0_parameters.hdf5'%directory,'r')
    for i in range(len(names_address[directory])):
        real_data[names_address[directory][i]] = []
        real_data[names_address[directory][i]].append(chain['source_mag_g'][i])
        real_data[names_address[directory][i]].append(chain['source_mag_r'][i])
        real_data[names_address[directory][i]].append(chain['source_mag_i'][i])
        real_data[names_address[directory][i]].append(chain['source_mag_z'][i])
        real_data[names_address[directory][i]].append(chain['source_mag_y'][i])
        real_data[names_address[directory][i]].append(chain['z_source'][i])

#calculating the g-i colours
for name in names:
    sim_color = float(mags[name].split(' ')[1])-float(mags[name].split(' ')[5])
    real_color = real_data[name][0]-real_data[name][2]
    color_error = (sim_color-real_color)
    real_data[name].append(sim_color)
    real_data[name].append(real_color)
    real_data[name].append(color_error)

#saving the colours
colorLines = []
colorLines.append('name sim_color real_color color_error\n')
for name in names:
    colorLine = '%s %s %s %s\n'%(name,real_data[name][-3],real_data[name][-2],real_data[name][-1])
    colorLines.append(colorLine)

f = open('paper_plots/colors.txt','w')
f.writelines(colorLines)
f.close()

#controlling the existing color errors
error_names = []
for name in names:
    for i in range(5):
        if float(mags[name].split(' ')[2+2*i])>1.0 or real_data[name][i]>26.0 or real_data[name][i]<20.0 or float(mags[name].split(' ')[1+2*i]) > 26.0 or float(mags[name].split(' ')[1+2*i]) < 20.0:
            error_names.append(name)
    if real_data[name][-2]<-0.5 or real_data[name][-2]>0.75 or real_data[name][-3]<-0.5 or real_data[name][-3]>0.75 or real_data[name][-1]>0.5:
        error_names.append(name)

#deleting the extremely wrong objects
error_names = numpy.unique(error_names)
for e_name in error_names:
    del(mags[e_name])
    del(real_data[e_name])

#preparing the plotters
realMagG = []
realMagR = []
realMagI = []
realMagZ = []
realMagY = []
magG = []
magR = []
magI = []
magZ = [] 
magY = []
Ger = []
Rer = []
Ier = []
Zer = []
Yer = []
specz = []
realColor = []
simColor = []
for name in real_data.keys():
    realMagG.append(real_data[name][0])
    realMagR.append(real_data[name][1])
    realMagI.append(real_data[name][2])
    realMagZ.append(real_data[name][3])
    realMagY.append(real_data[name][4])
    specz.append(real_data[name][5])
    simColor.append(real_data[name][6])
    realColor.append(real_data[name][7])

    magG.append(float(mags[name].split(' ')[1]))
    magR.append(float(mags[name].split(' ')[3]))
    magI.append(float(mags[name].split(' ')[5]))
    magZ.append(float(mags[name].split(' ')[7]))
    magY.append(float(mags[name].split(' ')[9]))

    Ger.append(float(mags[name].split(' ')[2]))
    Rer.append(float(mags[name].split(' ')[4]))
    Ier.append(float(mags[name].split(' ')[6]))
    Zer.append(float(mags[name].split(' ')[8]))
    Yer.append(float(mags[name].split(' ')[10]))

#real_mag = fit_mag line
x = numpy.arange(-30,30,1)
y = numpy.arange(-30,30,1)
fit = numpy.polyfit(x,y,1)
fitFunction = numpy.poly1d(fit)

#plotting the real vs fit mags and colours
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = '24',weight='bold')
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.left'] = plt.rcParams['ytick.labelleft'] = True
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.figure(figsize=(14, 21))
x_lim = 19.5
y_lim = 26.5
fontsize = 30
marker = 'x'
size = 12
plt.rcParams["errorbar.capsize"] = 0.
mfc = 'darkslategray'
mec = 'darkslategray'
color = 'olivedrab'
ecolor = 'lightslategray'
elinewidth = 0.75
alpha=0.01

#g-band
plt.subplot(321)
plt.plot(x,fitFunction(x),color=color, label = 'fit mag = true mag',alpha=alpha,linewidth=elinewidth+2,zorder=1) 
plt.errorbar(realMagG,magG,yerr=Ger,marker=marker,linestyle="None",markerfacecolor=mfc,label='simulated mock-sample',markersize=size,mec=mec,ecolor=ecolor,elinewidth=elinewidth,zorder=10)
plt.xlim(x_lim,y_lim)
plt.ylim(x_lim,y_lim)
plt.xlabel('$m_{g,true}$',fontsize=fontsize)
plt.ylabel('$m_{g,fit}$',fontsize=fontsize)
plt.legend(fontsize=22)

#r_band
plt.subplot(322)
plt.plot(x,fitFunction(x),color=color,alpha=alpha,linewidth=elinewidth+2,zorder=1) 
plt.errorbar(realMagR,magR,yerr=Rer,marker=marker,linestyle="None",markerfacecolor=mfc,markersize=size,mec=mec,ecolor=ecolor,elinewidth=elinewidth,zorder=10)
plt.xlim(x_lim,y_lim)
plt.ylim(x_lim,y_lim)
plt.xlabel('$m_{r,true}$',fontsize=fontsize)
plt.ylabel('$m_{r,fit}$',fontsize=fontsize)

#i-band
plt.subplot(323)
plt.plot(x,fitFunction(x),color=color,alpha=alpha,linewidth=elinewidth+2,zorder=1) 
plt.errorbar(realMagI,magI,yerr=Ier,marker=marker,linestyle="None",markerfacecolor=mfc,markersize=size,mec=mec,ecolor=ecolor,elinewidth=elinewidth,zorder=10)
plt.xlim(x_lim,y_lim)
plt.ylim(x_lim,y_lim)
plt.xlabel('$m_{i,true}$',fontsize=fontsize)
plt.ylabel('$m_{i,fit}$',fontsize=fontsize)

#y-band
plt.subplot(324)
plt.plot(x,fitFunction(x),color=color,alpha=alpha,linewidth=elinewidth+2,zorder=1) 
plt.errorbar(realMagZ,magZ,yerr=Zer,marker=marker,linestyle="None",markerfacecolor=mfc,markersize=size,mec=mec,ecolor=ecolor,elinewidth=elinewidth,zorder=10)
plt.xlim(x_lim,y_lim)
plt.ylim(x_lim,y_lim)
plt.xlabel('$m_{z,true}$',fontsize=fontsize)
plt.ylabel('$m_{z,fit}$',fontsize=fontsize)

#z-band
plt.subplot(325)
plt.plot(x,fitFunction(x),color=color,alpha=alpha,linewidth=elinewidth+2,zorder=1) 
plt.errorbar(realMagY,magY,yerr=Yer,marker=marker,linestyle="None",markerfacecolor=mfc,markersize=size,mec=mec,ecolor=ecolor,elinewidth=elinewidth,zorder=10)
plt.xlim(x_lim,y_lim)
plt.ylim(x_lim,y_lim)
plt.xlabel('$m_{y,true}$',fontsize=fontsize)
plt.ylabel('$m_{y,fit}$',fontsize=fontsize)

#colours
plt.subplot(326)
plt.plot(x,fitFunction(x),color=color,alpha=alpha,linewidth=elinewidth+2,zorder=1) 
plt.errorbar(realColor,simColor,yerr=Ier,marker=marker,linestyle="None",markerfacecolor=mfc,markersize=size,mec=mec,ecolor=ecolor,elinewidth=elinewidth,zorder=10)
plt.xlabel('true $g$-$i$ colour',fontsize=fontsize-8)
plt.ylabel('fit $g$-$i$ colour',fontsize=fontsize-8)
plt.xlim(-0.45,0.80)
plt.ylim(-0.45,0.80)
plt.xticks(numpy.arange(-0.25, 1.0, 0.5))
plt.yticks(numpy.arange(-0.25, 1.0, 0.25))
plt.tight_layout()
plt.savefig('paper_plots/all.eps',format='eps',dpi=1000)
plt.close()

#saving the data
finalLines = []
finalDataLines = []
for i in range(len(realMagG)):
    line = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(names[i],realMagG[i],magG[i],Ger[i],realMagR[i],magR[i],Rer[i],realMagI[i],magI[i],Ier[i],realMagZ[i],magZ[i],Zer[i],realMagY[i],magY[i],Yer[i],realColor[i],simColor[i],Ier[i],specz[i])
    lineData = "%s %s %s %s %s %s %s %s %s %s %s %s\n" %(names[i],magG[i],Ger[i],magR[i],Rer[i],magI[i],Ier[i],magZ[i],Zer[i],magY[i],Yer[i],specz[i])
    finalLines.append(line)
    finalDataLines.append(lineData)

f = open('paper_plots/output.txt','w')
f.writelines(finalLines)
f.close()

f = open('paper_plots/data.txt','w')
f.writelines(finalDataLines)
f.close()

################################
#for the lens to source SB ratio
os.makedirs('sources')

lines = []
lines.append('name ratio error\n')
for name in mags.keys():
    if name in names_address[directories[0]]: directory = directories[0]
    if name in names_address[directories[1]]: directory = directories[1]
    noSource = fits.open('../rangavar/13tousand/3quarry/data/%s_i_sci.fits'%name)
    imageFile = fits.open('%s/data/data/sim0_%s_i_sci.fits'%(directory,name))
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
    print(iMax)

    #lens again
    newnoSource = fits.open('../rangavar/13tousand/3quarry/data/%s_i_sci.fits'%name)
    lens = float(newnoSource[0].data[sourcex,sourcey])
    source = float(noSource[0].data[sourcex,sourcey])
    noSource.writeto('sources/%s_source.fits'%name)
    ratio = abs(float(lens/source))
    real_data[name].append(ratio)
    lines.append('%s %f %f\n'%(name,ratio,real_data[name][8]))
    count = count + 1

f = open('paper_plots/ratio.txt','w')
f.writelines(lines)
f.close()

#plotters
ratio = []
error = []
for name in mags.keys():
    if abs(real_data[name][-2]) > 3e-4:
        ratio.append(real_data[name][-1])
        error.append(abs(real_data[name][-2]))

#fit line
fitRatio = numpy.log10(ratio)
fitError = numpy.log10(error)
p = numpy.polyfit(fitRatio, fitError, 1)
fittedLine = numpy.poly1d(p)
x = numpy.arange(min(fitRatio-2),max(fitRatio+2),0.1)
y = fittedLine(x)

#plotting
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = '24',weight='bold')
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.left'] = plt.rcParams['ytick.labelleft'] = True
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['xtick.major.pad'] = 8
plt.figure(figsize=(10, 8))

plt.plot(10.**x, 10.**y, c='darkgoldenrod', linestyle='--', linewidth = 3, label='best-fit line',zorder=1)
plt.scatter(ratio,error,c='midnightblue',s=40,marker='x',label='simulated mock-images',zorder=10)
plt.xscale('log')
plt.yscale('log')
plt.ylim(3e-4,1)
plt.xlim(1e-3,10)
plt.legend(fontsize='22')
plt.xlabel('lens to source SB ratio ($\Lambda_{i}$)',fontsize=26)
plt.ylabel('$g$-$i$ error ($\delta C$)',fontsize=26)
plt.savefig('paper_plots/ratio.eps',format='eps', dpi=1000)
plt.close()
