import numpy as np
import matplotlib.pyplot as plt


def smoothen2D(sigma, coarsearray):

    ysize, xsize = coarsearray.shape
    
    windowWidth = 2*int(sigma)
    windowSize  = 2*windowWidth+1
    
    if windowSize >= ysize or windowSize >= xsize:
        print "error! use smaller sigma"
        print "the window size is too large for the number of points"
        return smoothen2D(sigma/2., coarsearray)

    window = np.zeros([windowSize, windowSize])
    normalise = 0.
    for i in range(windowSize):
        for j in range(windowSize):
            d = np.sqrt(float(i)**2+ float(j)**2)/sigma
            val = np.exp(-d)
            normalise += val
            window[i,j] = val
    window = window / normalise

    smootharray = np.zeros([ ysize -windowSize,xsize-windowSize])
    for j in range(ysize-windowSize):
        for i in range(xsize-windowSize):

            sliceCoarse = coarsearray[j:j+windowSize, i:i+windowSize]

            smootharray[j,i]  = sum(sum(sliceCoarse * window) )
    return smootharray, windowWidth

class findLocalMaxima:
    """ this class will identify the local maxima, 
    initialise with the xvalues, yvalues and the signals,
    then call get getTopNMaxima"""
    def __init__(self, xvalues, yvalues, signals):
        self.xvalues = np.array(xvalues)
        self.yvalues = np.array(yvalues)
        self.signals = np.array(signals)
        self.xsize = len(xvalues)
        self.ysize = len(yvalues)
        self.makemin = np.zeros([self.ysize, self.xsize])
        print self.xsize, self.ysize, self.makemin.shape,self.signals.shape, signals.shape

    def getTopNMaxima(self, n, ressofar = []):
        """ getTopNMaxima
        return the top n maxima"""    
        maxsignal = max(self.signals.reshape(self.signals.size))
        self.minsignal = min(self.signals.reshape(self.signals.size))
        
        maxfirstIndex = np.where(self.signals==maxsignal)[0][0]
        maxseconIndex = np.where(self.signals==maxsignal)[1][0]
        print maxfirstIndex, maxseconIndex

        res = ( self.xvalues[maxseconIndex],self.yvalues[maxfirstIndex])
        ressofar.append(res)
        n -= 1

        self.removeLocalMaximum(maxfirstIndex, maxseconIndex, maxsignal) 

        if n == 0 :
            return ressofar
        else:
            return self.getTopNMaxima(n,ressofar)

    def removeLocalMaximum(self, index1, index2, maxvalue):
        self.tagLocalMaximum(index1, index2, maxvalue)

        for i in range(self.ysize):
            for j in range(self.xsize):
                if self.makemin[i,j] == 1:
                    self.signals[i,j] = self.minsignal


    def tagLocalMaximum(self, indexY, indexX, maxvalue):

        #determine largest value and set others to min
        #also work out the average
        val = self.signals[indexY,indexX]
  
        # if the value in the neighbour is smaller to this box, exit otherwise tag as a minimum
        if maxvalue < val:
            return       
        self.makemin[indexY,indexX] = 1
        #otherwise remove othe maxims
        for up in range(indexY-1,indexY+2):
            for left  in range(indexX-1,indexX+2):
                nextx = min(max(left, 0), self.xsize-1) 
                nexty = min(max(up,   0), self.ysize-1)

                #if the neighbour is not already tagged, investigate it for taggin
                if self.makemin[nexty, nextx] == 0:
                    self.tagLocalMaximum( nexty,nextx, val)


if __name__ == '__main__':
    xvalues = np.arange(0.,100.)
    yvalues = np.arange(0, 50.)
    signals = np.array([[np.cos(x)+np.sin(y) +np.random.rand() for x in xvalues] for y in yvalues])
    smoothsig,windowWidth = smoothen2D(10.,signals)
    print smoothsig.shape
    xtmp = xvalues[windowWidth:len(xvalues)-windowWidth-1]
    ytmp =  yvalues[windowWidth:len(yvalues)-windowWidth-1]
    print xtmp.shape
    print ytmp.shape
    findmax = findLocalMaxima(xtmp, ytmp, smoothsig)


    print findmax.getTopNMaxima(1)

    plt.figure(1)
    plt.contour(xtmp, ytmp, findmax.signals)
    plt.figure(2)
    plt.contour(xvalues, yvalues, signals)
    plt.show()