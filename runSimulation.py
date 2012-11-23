
    for p in range(1, 110):
        i = p-1
        print p, " actual: ", trainingHalo.coord1[i], trainingHalo.coord2[i], trainingHalo.coord3[i]
        skyRawData  = numpy.loadtxt('Train_Skies/Training_Sky%i.csv' % p, \
                                 delimiter=',',unpack=True,usecols=(1,2,3,4),skiprows=1)
        sky   = TrainingSky(skyRawData, i )
        xvalues,yvalues, signals = model.rankLikelyhood(sky,100., 0., 4200., 0., 4200.)

        smoothsig,windowWidth = smoothen2D(3.,signals)     
        xtmp = xvalues[windowWidth:len(xvalues)-windowWidth-1]
        ytmp =  yvalues[windowWidth:len(yvalues)-windowWidth-1]
        findmax = findLocalMaxima(xtmp, ytmp, smoothsig)

        print "predicted: " , findmax.getTopNMaxima(3,[])


     #   plt.contour(xtmp,ytmp, smoothsig)
    #    plt.show()