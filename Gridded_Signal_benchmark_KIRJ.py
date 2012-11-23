
import numpy as np
import csv as c

if __name__ == "__main__":
    n_skies = 3
    position_halo=np.zeros([n_skies,2,3],float) #Set up the array in which I will
                                                #assign my estimated positions
    
    nhalo=np.array([1,1,1])
            
    for k in xrange(n_skies):
        p=k+1
        #Read in the x,y,e1 and e2 positions of
        #each galaxy in the list for sky number k:
        x,y,e1,e2=np.loadtxt('Train_Skies/Training_Sky%i.csv'\
                             % p,delimiter=',',unpack=True,usecols=(1,2,3,4),skiprows=1)

        #So I need to grid the sky up. Here I set the parameters of the grid.
        nbin=15 #Number of bins in my grid
        image_size=4200.0 #Overall size of my image
        binwidth=float(image_size)/float(nbin) # The resulting width of each grid section
    
        average_tan_force=np.zeros([nbin,nbin],float) #Set up the signal array
                                                      #in which Im going to find
                                                      #the maximum of.
    
        for i in xrange(nbin):
            for j in xrange(nbin):
            
                x0=i*binwidth+binwidth/2. #I set the proposed x position of the halo
                y0=j*binwidth+binwidth/2. #I set the proposed y position of the halo
            
                angle_wrt_halo=np.arctan((y-y0)/(x-x0)) #I find the angle each
                                                    #galaxy is at with respects
                                                    #to the centre of the halo.               
                tangential_force=-(e1*np.cos(2.0*angle_wrt_halo)\
                               +e2*np.sin(2.0*angle_wrt_halo))
                               #Find out what the tangential force
                               #(or signal) is for each galaxy with
                               #respects to the halo centre, (x0,y0)
                tangential_force_in_bin=tangential_force[(x >= i*binwidth) & \
                                                     (x < (i+1)*binwidth) & \
                                                     (y >= j*binwidth) & \
                                                     (y < (j+1)*binwidth)]
                                    #Find out which galaxies lie within the gridded box


                if len(tangential_force_in_bin) > 0:
                    average_tan_force[i,j]=sum(tangential_force_in_bin)\
                        /len(tangential_force_in_bin) #Find the average signal per galaxy
                else:
                    average_tan_force[i,j]=0
                

        index=np.sort(average_tan_force,axis=None) #Sort the grid into the
                                                   #highest value bin first,
                                                   #which should be the centre
                                                   #of one of the halos
        index=index[::-1] #Reverse the array so the largest is first
        for n in xrange(int(nhalo[k])):
            position_halo[k,0,n]=np.where(average_tan_force\
                                          == index[n])[0][0]\
                                          *binwidth
                                          #For each halo in the sky find
                                          #the position and assign
            position_halo[k,1,n]=np.where(average_tan_force\
                                          == index[n])[1][0]\
                                          *binwidth
                                          #The three grid bins
                                          #with the highest signal should
                                          #contain the three halos.
    for i in range(n_skies):
      print position_halo[i][0][0], position_halo[i][1][0] 
  