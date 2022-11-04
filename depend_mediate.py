import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import pylab as py
from matplotlib import pyplot
import os, sys, re
from os.path import join
####################################  COOL STARS    ###########################################3
numf=int(5);     numm=int(22);      numw=int(586)
norm=np.zeros((numf))            
norm=[340.47576713,780.34789551,986.66688173,1050.56463272,692.7877672]
c=['m','b','g','y','r']
dirname="./Patrick_Harrington/mediate/"

filr=open("./prop.txt","a"); 
filr.close()
counter=121

param=np.zeros((4))
for fname in os.listdir(dirname):
    #print fname 
    fname1=join(dirname,fname)   
    print "reading the file: ",  fname1
    counter+=1
    
    param=np.array([float(i) for i in re.findall(r"[-+]?\d*\.\d+|\d+",str(fname))])
    fobj=open(fname1,'r')
    #line1=fobj.readline()    ###       Teff= 15000 log g= 2.0, Z= 1.0
    #dummy1,  Teff, dummy2, du2, logg, dummy3 ,FeH = line1.split()
    Teff=int(param[0]);  logg=float(param[1])
    print "The parameters:   Teff= ",  Teff," logg=",  logg,   "wave: ",  param[2],   param[3]
    #gra = logg.split(",")
    #print "logg: ",  gra[0], type(gra)
   
    mu=np.zeros((numm))
    mu=fobj.readline().split()
   # print mu 
   
    #### blink  
    blink=fobj.readline().split()
   
    wave=np.zeros((numw));  
    wave=fobj.readline().split()
    #print "wave: ",  wave
    blink=fobj.readline().split() 

    #####  blink    
     
    Ii=np.zeros((numw))  
    Ii=fobj.readline().split()
   # print "Ii: ",  Ii 
    
    blink=fobj.readline().split() 
    #input("Enter a number")
  
    #angle=np.zeros((numm))
    #angle=fobj.readline().split()
    #input ("Enter a number ")   
    '''
#### Here we read the transmission function for UBVRI filters#####################
    # numf= 5# number of filter:  U B  V  R  I 
    trans=np.zeros((5,1801,2));  
    trans[0,:,:]=np.loadtxt("./transmission/Bessel_U-1.txt")
    trans[1,:,:]=np.loadtxt("./transmission/Bessel_B-1.txt")
    trans[2,:,:]=np.loadtxt("./transmission/Bessel_V-1.txt")
    trans[3,:,:]=np.loadtxt("./transmission/Bessel_R-1.txt")
    trans[4,:,:]=np.loadtxt("./transmission/Bessel_I-1.txt")
    for i in range(1801):
        for j in range(5):
            trans[j,i,1]=float(abs(trans[j,i,1])/100.0) ### transmision function  normalized
            trans[j,i,0]=float(abs(trans[j,i,0])*10.0)### wave length in Angestrum   
######## averaging on transimssion for each filter over 'numw' wave length point########################     
    tranm=np.zeros((numf,numw))        
    for k in range(numf):    
        for i in range(numw):
            tranm[k,i]=0.0
            num=int(0.0) 
            if(i!=int(numw-1)):  up=float(wave[i+1])
            else:                up=float(wave[numw-1])+12.0
            for j in range(1801):
                wav=float(trans[k,j,0])
                if(wav>=float(wave[i]) and wav<=up):
                    tranm[k,i]+= float(trans[k,j,1]); num+=1;  
            if(num>0):     tranm[k,i]=float(tranm[k,i]/num/1.0)
            else:          tranm[k,i]=float(0.0)        
            #print "filter", k, "wavelenght: ", float(wave[i]), " ave_transmission: ",  tranm[k,i], "num: ", num,"Wave_band: ",  up, wave[i] 
            if(float(tranm[k,i])>1.0 or float(tranm[k,i])<0.0 ):
                print "Error ,  tranm",  tranm[k,i], "num: ", num
                input("Enter a number")  
        print "******************************************************"               
#########  saving tramsnission in a file and calculating normalization factor && plotting it   
    fii=open("./transmission/transmission_mediate.txt","w");  fii.close()           
    norm=np.zeros((numf))
    d=np.zeros((numw,numf+1))            
    for i in range(numw):            
        fii=open("./transmission/transmission_mediate.txt","a")
        d[i,0]=float(wave[i]); 
        if(i!=int(numw-1)):  up=float(wave[i+1])
        else:                up=float(wave[numw-1])+12.0  
        for j in range(numf): 
            norm[j]+= abs(tranm[j,i])*abs(up-float(wave[i]))
            d[i,j+1] =float(tranm[j,i])
        np.savetxt(fii,d[i,:].reshape((-1,6)),"%.5f   %.10f   %.10f   %.10f   %.10f   %.10f")
        fii.close()
    print "Normalization factor: ",  norm    
    plt.clf()
    plt.plot(d[:,0], d[:,1],'p-',lw=2.0)
    plt.plot(d[:,0], d[:,2],'b-',lw=2.0)
    plt.plot(d[:,0], d[:,3],'g-',lw=2.0)
    plt.plot(d[:,0], d[:,4],'y-',lw=2.0)
    plt.plot(d[:,0], d[:,5],'r-',lw=2.0)
    fig=plt.gcf()
    fig.savefig("./transmission/transmission_mediate.jpg",dpi=200)
    
    sys.exit()
    '''          
    
          
     
     
########## reading the rest of file#############################################


    tranq=np.zeros((numw,numf+1))    
    tranq=np.loadtxt('./transmission/transmission_mediate.txt')

    stokI=np.zeros((numf,numm))    
    stokQ=np.zeros((numf,numm))
    stok_I= np.zeros((numw,numm))
    stok_Q= np.zeros((numw,numm))
    pol= np.zeros((numw,numm))
    for i in range(numw):
        stok_I[i,:]= fobj.readline().split()
    blink=fobj.readline().split()     
    for i in range(numw):
        stok_Q[i,:]= fobj.readline().split()
    blink=fobj.readline().split() 
    for i in range(numw):    
        pol[i,:]   = fobj.readline().split()
    
        for j in range(numm):
            if(abs( float(stok_Q[i,j])/float(stok_I[i,j]) - float(pol[i,j]) ) >0.00001):
                print "pol: ", pol[i,j], "Q/I: ", float(stok_Q[i,j])/float(stok_I[i,j])
                input("Enetr a number ")  
        fil=open("./wave/W{0:f}.txt".format(float(wave[i])),"w");  fil.close(); tt=np.zeros((numm,3))
        for j in range(numm):
            fil=open("./wave/W{0:f}.txt".format(float(wave[i])),"a")
            tt[j,0]=float(mu[j]);  tt[j,1]=abs(float(stok_I[i,j]));  tt[j,2]=float(stok_Q[i,j])
            np.savetxt(fil,tt[j,:].reshape((1,3)),"%.10f   %.12f     %.12f")
            fil.close() 
    ################################################
        ''' 
        plt.clf() 
        plt.plot(1.0-tt[:,0], tt[:,1],'r-', lw=2.0)
        plt.plot(1.0-tt[:,0], tt[:,2],'b-', lw=2.0) 
        plt.plot(1.0-tt[:,0],np.divide(tt[:,2],tt[:,1]),'g-', lw=2.0)  
        plt.xlabel(r"$1-\mu$",fontsize=18,labelpad=0.1)
        plt.title(r"$wavelength={0:d}$".format(int(wave[i])),fontsize=18)
        fig=plt.gcf()    
        fig.savefig("./wave/IQ{0:d}.jpg".format(i),dpi=200)
        '''
    ################################################
        if(i!=int(numw-1)):  up=float(wave[i+1])
        else:                up=float(wave[numw-1])+12.0   
        for j in range(numf):
            stokI[j,:] += stok_I[i,:] * abs(tranq[i,j+1]) * abs(1000000000.0/float(wave[i])/float(wave[i])) * abs(up-float(wave[i]))
            stokQ[j,:] += stok_Q[i,:] * abs(tranq[i,j+1]) * abs(1000000000.0/float(wave[i])/float(wave[i])) * abs(up-float(wave[i]))

        
        
        
    ################################################    
    plt.clf() 
    fil=open("./filter2/{0:d}_{1:.1f}.txt".format(int(Teff),float(logg)),"w");  
    fil.close(); tt=np.zeros((6))
    for j in range(numm):
        tt[0]=float(mu[j]); 
        for i in range(numf):
            tt[i+1]=abs(float(stokI[i,j])*10000.0/norm[i]); # tt[j,2]=float(stokQ[i,j])/float(norm[i])
        fil=open("./filter2/{0:d}_{1:.1f}.txt".format(int(Teff),float(logg)),"a");  
        np.savetxt(fil,tt.reshape((1,6)),fmt="%.4f   %.9f    %.9f   %.9f    %.9f    %.9f")
        fil.close() 
        #plt.plot(1.0-tt[:,0], np.divide(abs(tt[:,2]),tt[:,1])*100.0,color=c[i], lw=2.0)
        #plt.ylabel(r"$Polarization [\%]$",fontsize=18,labelpad=0.1)
        #plt.xlabel(r"$1-\mu$",fontsize=18,labelpad=0.1)
    #fig=plt.gcf()    
    #fig.savefig("./filter/pol{0:d}_{1:.1f}.jpg".format(int(Teff),float(gra[0])),dpi=200)
    ################################################
    
    print "*******************************************"
    filr=open("./prop.txt","a"); 
    prop=np.zeros((3))
    prop[2]=numm;  prop[0]=Teff;  prop[1]=float(logg) 
    np.savetxt(filr,prop.reshape((1,3)),"%d    %.2f   %d")
    filr.close()
    
    print "*******************************************"
    #input ("Enter a number ")
    
    
    
