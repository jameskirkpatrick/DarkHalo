import pickle as pkl
import csv as c
from TrainingSets import *
from findMaxima import *

if __name__ == '__main__':
    print "load data: "
    f = file("C:\Users\kirj\Documents\DarkHalo\modelFit.dat")
    modelImporter = pickle.Unpickler(f)
    model  = modelImporter.load()
    print "data loaded"

    countHalos = numpy.loadtxt('Training_haloCounts.csv' , \
                                 delimiter=',',unpack=True,usecols=(1,),skiprows=1)
    trainingHalo = TrainingHalo("Training_halos.csv")

    n_skies  =len(countHalos)
    position_halo = np.zeros([n_skies,2,3])
    p=1
    for nhalos in countHalos:
        i = p-1
        print p, " actual: ", trainingHalo.coord1[i], trainingHalo.coord2[i], trainingHalo.coord3[i]
        skyRawData  = numpy.loadtxt('Train_Skies/Training_Sky%i.csv' % p, \
                                 delimiter=',',unpack=True,usecols=(1,2,3,4),skiprows=1)
        sky   = TrainingSky(skyRawData, i )
        xvalues,yvalues, signals = model.rankLikelyhood(sky,100., 0., 4200., 0., 4200.)
        # if we are looking for multiple halos, smoothen the function
        if nhalos > 1:
            smoothsig,windowWidth = smoothen2D(1.,signals)     
            xtmp = xvalues[windowWidth:len(xvalues)-windowWidth-1]
            ytmp =  yvalues[windowWidth:len(yvalues)-windowWidth-1]
            findmax = findLocalMaxima(xtmp, ytmp, smoothsig)
        else:
            findmax =  findLocalMaxima(xvalues, yvalues, signals)   
        res =  findmax.getTopNMaxima(nhalos,[])
        #pad with 0s when there are less then 3
        while len(res) < 3:
            res.append([0.,0.])
        print "computed: " ,res
        res = np.array(res)
        for xy in range(2):
            for halo in range(3):
                position_halo[i][xy][halo] = res[halo][xy]

        p+=1

#output to cvs:
    c = c.writer(open("TrainingSets_pdfmethod.csv", "wb")) #Now write the array to a csv file
    c.writerow([str('SkyId'),str('pred_x1'),str( 'pred_y1'),str( 'pred_x2'),str( 'pred_y2'),str( 'pred_x3'),str(' pred_y3')])
    for k in xrange(n_skies):
        halostr=['Sky'+str(k+1)] #Create a string that will write to the file
                      #and give the first element the sky_id
        for n in xrange(3):
            halostr.append(position_halo[k,0,n]) #Assign each of the
                                             #halo x and y positions to the string
            halostr.append(position_halo[k,1,n])
        c.writerow(halostr) #Write the string to a csv
                        #file with the sky_id and the estimated positions