""" TrainingSets.py

load the training data and save it into a small 

"""
import numpy as numpy 
import matplotlib.pyplot as plt
import pickle
from scipy import stats

from findMaxima import *
from geometriclib import *

class TrainingHalo:
    """ TrainingHalo
    loads the information on the position of dark matter halos for training purposes
    """
    def __init__(self, namefile):
        self.nhalo  = numpy.loadtxt(namefile, \
                             delimiter=',',unpack=True,usecols=(1,),skiprows=1)

        xcoord1 ,ycoord1,xcoord2 ,ycoord2, xcoord3 ,ycoord3 = numpy.loadtxt(namefile, \
                             delimiter=',',unpack=True,usecols=(4,5,6,7,8,9),skiprows=1)

        self.coord1 = zip(xcoord1, ycoord1)
        self.coord2 = zip(xcoord2, ycoord2)
        self.coord3 = zip(xcoord3, ycoord3)

class TrainingSky:
    """TrainingSky
    This class contains the ellipticities (first and second) and the position (x and y) of each galaxy in a TrainingSky.
    number contains the index for the TrainingSky (counting from 1)
    """
    def __init__(self,(xcoord ,ycoord,ellipticity1 ,ellipticity2), numberOfSky):
    
        self.number = numberOfSky
        
        self.coord = zip(xcoord, ycoord)
        self.ellipticityCartesian = zip(ellipticity1,ellipticity2)

        self.ellipmodulus      = numpy.array([ modulus(c)for c in self.ellipticityCartesian] )
        self.elliptangle       = numpy.array([ angle(e)/2.0 for  e in self.ellipticityCartesian])
        self.xcoord = xcoord
        self.ycoord = ycoord
        
    def computeHaloInformation(self,trainingHalo):
        halocoord = trainingHalo.coord1[self.number-1]

        #compute the distance from the halo:
        dcoord = numpy.array( [diff(c, halocoord) for c in self.coord])
        self.distance = numpy.array([modulus(dc) for dc in dcoord ])
        #compute the angle to the center of the halo
        self.angle = numpy.array ([angle(dc) for dc in dcoord])
        #compute the components of the ellipticity perpendicular and at 45 degrees to the halo
        self.ellipticityNormal = numpy.array([rotateEllipticity(ellipse, a) for ellipse, a in zip (self.ellipticityCartesian, self.angle) ])        

    def calculateAngleAndDistance(self, x,y):
        halocoord = (x,y)

        #compute the distance from the halo:
        dx = self.xcoord - x
        dy = self.ycoord - y
        distance = numpy.sqrt( dx**2 + dy**2 )
        
        #compute the angle to the center of the halo
        angles = numpy.arctan2 (dy,dx)
        #compute the components of the ellipticity perpendicular and at 45 degrees to the halo
        relativeAngles = self.elliptangle - angles
        return zip ( relativeAngles, distance)
            
class collectStatistics:
    """collectStatistics
    given several Training_Sky files, generate the statistics on pdf for angles and ellipticity

    """
    @staticmethod
    def bin2D(x,y,dx,dy):
        """bin2D
        given a set of two coordinates x and y construct the bincount of when x is within dx of a certain value
        and when y is within dy of another one 
        """
        xmin = min(x)
        ymin = min(y)
        xmax = max(x)
        ymax = max(y)

        nx = int((xmax-xmin)/dx)+1
        ny = int((ymax-ymin)/dy)+1

        binCounts = numpy.zeros(nx*ny)
        binCounts = binCounts.reshape(nx,ny)

        for xit, yit in zip(x,y):
            indexX = (xit - xmin)/dx
            indexY = (yit - ymin)/dy

            binCounts[indexX, indexY] += 1

        xvalues = dx*numpy.array(range(nx))+xmin+dx/2.0
        yvalues = dy*numpy.array(range(ny))+ymin+dy/2.0    
        
        return xvalues, yvalues, binCounts

    @staticmethod    
    def normaliseByFirstCoordinate(binCounts2D):
        """ normaliseByFirstCoordinate
        takes a 2D binCount and normalises each row by their sum
        """
        row_sum = binCounts2D.sum(axis=1)
        return binCounts2D / row_sum[:,numpy.newaxis]

    @staticmethod   
    def fitToCosine(thetaRange, angleDependence ):
        """fitToSine fit a set of angle depenendent data to the formula: 
            intercept + amplitude (cos(2*theta))"""
        logAngleDependence = numpy.array([numpy.log(angle) for angle in angleDependence+1E-10])
        cosines = numpy.array([cos(2.0*theta)  for theta in thetaRange])
        amplitude , intercept, r_value, p_value, std_err = stats.linregress(cosines, logAngleDependence)
        

        return (intercept, amplitude)
        
    @staticmethod
    def fitAmplitudeOfExpCosRelationToPowerLaw(distanceRange, amplitudes):
        """fitAmplitudeOfExpCosRelationToPowerLaw 
        given a vector of distances and the amplitude B for the pdf at that distance to the halo
        it determines the coefficients of the relationship:
        B(r) = B0 r**alpha

        B(r) has the same meaning as the pdf fitted in fitToCosine, namely:
        pdf(theta) = exp(-B(r) * cos(2*theta) - A)
        """
        logdistance  = numpy.array([numpy.log(d) for d in distanceRange])
        logamplitude = numpy.array([numpy.log(a) for a in amplitudes])

        alpha , logB0, r_value, p_value, std_err = stats.linregress(logdistance, logamplitude)
        return (numpy.exp(logB0), alpha)
    
class modelpdf:
    """modelpdf
    get the probability density function for a given angle
    """

    def __init__(self, alpha=0., b0=0.):
        self.alpha = alpha
        self.b0 = b0

    def fitPDFToExpCosineModel(self, allDistances, relativeAngle, distBinWidth,angleBinWidth):
        """ fitPDFToExpCosineModel
        givena the value of the distances of all galaxies to dark matter
        halos allDistances and the relative angle relativeAngle between the 
        orientation of the ellipticity and the position of the galaxy, 
        determine the coefficients of the probability density function
        for the relative angle theta which is assumed to be of the form:

        pdf(theta) = exp(-B(r) * cos(2*theta) - A(r))

        with the B coefficient being of the form:
        B(r) = B0 r**alpha

        The B coefficient is determined numerically by fitting the relationship between the 
        frequency of occurance of relative angles as a function of distance. 

        The A coefficient is determined analytically by noting that:

        integral_-pi^pi (pdf(theta) dtheta ) =  1 = 2 pi * exp(B(r)) I0(A(r)) 
        where I0 is the modified bessel function of the first kind
        this means that the formula for A must be:

        A(r) = 1 / (2*pi * exp(B(r)))
        
        distBinWidth and angleBinWidth set the width of bins for distances and angles respectively 
        """

        distanceBins, angleBins,binCountAngle = collectStatistics.bin2D(allDistances, relativeAngle, distBinWidth,angleBinWidth )
        pdfFits = numpy.array( [collectStatistics.fitToCosine(angleBins, row) for row in binCountAngle] )
        bCoefficients = numpy.array([-fit[1] for fit in pdfFits])
        self.b0, self.alpha = collectStatistics.fitAmplitudeOfExpCosRelationToPowerLaw(distanceBins, bCoefficients)

    def likelyhood(self,theta, r):
        """likelyhood
        computes the likelyhood that an ellipticity with angle theta at distance r is recorded"""
        br = self.b0 * r**self.alpha
        expar = (1./(2.*numpy.pi * numpy.i0(br)))
        return numpy.exp(-(br*cos(2*theta))) * expar

    def preComputeLikelyhood(self, dtheta, dr, width=4200.):
        drmax = width * numpy.sqrt(2)
        rrange     = numpy.arange(0., drmax, dr)
        thetarange =  numpy.arange(-numpy.pi, numpy.pi, dtheta)
        self.preComputedLikelyhood = numpy.zeros([len(thetarange), len(rrange)])
        for ir in range(len(rrange)):
            for itheta in range(len(thetarange)):
                #taking the log of the likelyhood means I can simply sum rather than multiplying
                self.preComputedLikelyhood[itheta, ir] = numpy.log(self.likelyhood(itheta*dtheta+dtheta/2., ir*dr+dr/2.))               
        self.dtheta = dtheta
        self.dr = dr

    def rankLikelyhood(self, sky, resolution,xmin,xmax, ymin, ymax):
        xvalues = numpy.arange(xmin, xmax, resolution)
        yvalues = numpy.arange(ymin, ymax, resolution)
        maxsignal = 0.
        signal = 0.
        maxx = 0
        maxy = 0
        signals = numpy.zeros([len(yvalues), len(xvalues)])
        countx=0
        county=0
        for x in xvalues:
            for y in yvalues:
                allThetaAndD  = sky.calculateAngleAndDistance(x,y)
                signals[county,countx]    = \
                 sum(numpy.array([self.preComputedLikelyhood[int(theta/self.dtheta), int(r/self.dr)] for theta, r in allThetaAndD]))
                county+=1
            county = 0    
            countx += 1    
        return (xvalues, yvalues, signals)

def main():
    trainingHalo = TrainingHalo("Training_halos.csv")
    allDistances = numpy.array([])
    allAngles    = numpy.array([])
    allElliptAngle = numpy.array([])
    for i in range(1,150):
        # for the time being run the analysis only in cases when there is only 1 halo
        if trainingHalo.nhalo[i-1] == 1 :
            namefile = 'Train_Skies/Training_Sky%i.csv' % i
            skyRawData  = numpy.loadtxt(namefile, \
                             delimiter=',',unpack=True,usecols=(1,2,3,4),skiprows=1)
            exampleSky   = TrainingSky(skyRawData, i)
            exampleSky.computeHaloInformation(trainingHalo)
            allDistances   = numpy.concatenate((allDistances,   exampleSky.distance))
            allAngles      = numpy.concatenate((allAngles,      exampleSky.angle))
            allElliptAngle = numpy.concatenate((allElliptAngle, exampleSky.elliptangle))

    # get the elliptical angle relative to the position of the halo
    relativeAngle =  numpy.array([putAngleInNormalRange(a) for a in allElliptAngle-allAngles])
    model = modelpdf()
    model.fitPDFToExpCosineModel(allDistances,relativeAngle,1000,0.2)

    print "model parameters: " , model.alpha, model.b0

    #save the model object here:

    f = file("C:\Users\kirj\Documents\DarkHalo\modelFit.dat", 'w')

    #test it on the first sky:
    print "precompute likelyhood"
    model.preComputeLikelyhood(numpy.pi/200., 4200/200.)
    print "precomputing finished"
    exporter = pickle.Pickler(f)
    exporter.dump(model)

if __name__ == '__main__':
    main()
