from Training_halos import *
 
def exploratoryPlots():
    trainingHalo = TrainingHalo("Training_halos.csv")
    allDistances = numpy.array([])
    allAngles    = numpy.array([])
    allElliptAngle = numpy.array([])
    allElliptMag   = numpy.array([])
    allElliptTangential = numpy.array([]);
    allElliptCross = numpy.array([]);
    for i in range(1,100):
        # for the time being run the analysis only in cases when there is only 1 halo
        if trainingHalo.nhalo[i-1] == 1 :
            namefile = 'Train_Skies/Training_Sky%i.csv' % i
            skyRawData  = numpy.loadtxt(namefile, \
                             delimiter=',',unpack=True,usecols=(1,2,3,4),skiprows=1)
            exampleSky   = TrainingSky(skyRawData, i)
            exampleSky.computeHaloInformation(trainingHalo)
            allDistances   = numpy.concatenate((allDistances,   exampleSky.distance))
            allAngles      = numpy.concatenate((allAngles,      exampleSky.angle))
            allElliptMag   = numpy.concatenate((allElliptMag,   exampleSky.ellipmodulus))
            allElliptAngle = numpy.concatenate((allElliptAngle, exampleSky.elliptangle))
            allElliptTangential = numpy.concatenate((allElliptTangential, exampleSky.ellipticityNormal[0]))
            allElliptCross      = numpy.concatenate((allElliptCross, exampleSky.ellipticityNormal[1]))


# get the elliptical angle relative to the position of the halo
    relativeAngle =  numpy.array([putAngleInNormalRange(a) for a in allElliptAngle-allAngles])
#  bins the angle results in 2D
    ddist=750
    dangle=0.2
    dellipt=0.05
    distanceBins, angleBins,binCountAngle = collectStatistics.bin2D(allDistances, relativeAngle, ddist,dangle )
    distanceBins, ellipBins,binCountEllip = collectStatistics.bin2D(allDistances, allElliptMag, ddist, dellipt)    
    binCountAngle = collectStatistics.normaliseByFirstCoordinate(binCountAngle)
    binCountEllip = collectStatistics.normaliseByFirstCoordinate(binCountEllip)

    pdfFits = numpy.array( [collectStatistics.fitToCosine(angleBins, row) for row in binCountAngle] )
    sequenceOfInteRest=4
    (intercept, slope) = collectStatistics.fitToCosine(angleBins,binCountAngle[sequenceOfInteRest,:])
    fittedPDF = numpy.array([math.exp(intercept + slope * (cos(2*theta))) for theta in angleBins])

 # do all plotting           
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(angleBins, binCountAngle[sequenceOfInteRest,:])
    plt.hold(True)
    ax.plot(angleBins,fittedPDF)
    plt.hold(False)
    
    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(111)
    ax2.loglog(distanceBins, numpy.array([-fit[1] for fit in pdfFits]),'bo')
    plt.hold(True)
    ax2.loglog(distanceBins, 350**0.65 * distanceBins**(-0.65), 'ro')
    fig3 = plt.figure(3)
    ax3 = fig3.add_subplot(111)
    ax3.plot(distanceBins, numpy.array([fit[0] for fit in pdfFits])) 

    plt.show()
        
 