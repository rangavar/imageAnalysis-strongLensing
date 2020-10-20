import matplotlib
matplotlib.use('Agg')
import h5py
import numpy
import matplotlib.pyplot as plt
import glob
import os

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
    color_error = abs(sim_color-real_color)
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
    if real_data[name][-1]<0.5 or real_data[name][-1]>0.75 or real_data[name][-2]<0.5 or real_data[name][-2]>0.75:
        real_data[name][-1]< 0.5

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
