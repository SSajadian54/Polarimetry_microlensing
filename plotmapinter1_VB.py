import numpy as np
import matplotlib.pyplot as plt
import pylab as py
from matplotlib import rcParams
import Image
from numpy import ma
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
#import matplotlib.cm as cm
#import matplotlib.mlab as mlab
#matplotlib.rcParams['xtick.direction'] = 'out'
#matplotlib.rcParams['ytick.direction'] = 'out'



rhostar=float(0.0069717883314647172490752)
RA=float(np.pi/180.0)
n=36864
array=np.zeros((800,2))
array=np.loadtxt("caustic_inter1.txt") 
data=np.zeros((n,7))
data=np.loadtxt("Map_pol_inter1_VB.txt") 
###############################################################################

ny=int(192)
nx=int(192)

count=int(0)
mapp=np.zeros((ny*2-1,nx*2-1))
mapt=np.zeros((ny*2-1,nx*2-1))
mapm=np.zeros((ny*2-1,nx*2-1))
maxx=float(-10.0)
for i in range(nx):
    x= 0.0 +i*rhostar*0.75;
    for j in range(ny):
        y= 0.0 + j*rhostar*0.75;
        if(abs(x-data[count,0])<float(0.1*rhostar) and abs(y-data[count,1])<float(0.1*rhostar)):
            if(j==0 and i==0):
                mapp[ny-1,nx-1]=float(np.log10(data[count,3]))
                mapt[ny-1,nx-1]=float(abs(data[count,5]))
                mapm[ny-1,nx-1]=float(np.log10(data[count,2]))
            elif (j==0 and i!=0): 
                mapp[ny-1,nx-1+i]=float(np.log10(data[count,3]))
                mapp[ny-1,nx-1-i]=float(np.log10(data[count,3]))
                mapt[ny-1,nx-1+i]=float(abs(data[count,5]))
                mapt[ny-1,nx-1-i]=float(abs(data[count,5]))
                mapm[ny-1,nx-1+i]=float(np.log10(data[count,2]))
                mapm[ny-1,nx-1-i]=float(np.log10(data[count,2]))
            elif (j!=0 and i==0): 
                mapp[ny-1+j,nx-1]=float(np.log10(data[count,3]))
                mapp[ny-1-j,nx-1]=float(np.log10(data[count,3]))
                mapt[ny-1+j,nx-1]=float(abs(data[count,5]))
                mapt[ny-1-j,nx-1]=float(abs(data[count,5]))
                mapm[ny-1+j,nx-1]=float(np.log10(data[count,2]))
                mapm[ny-1-j,nx-1]=float(np.log10(data[count,2]))
            else:
                mapp[ny-1+j,nx-1+i]=float(np.log10(data[count,3]))
                mapp[ny-1-j,nx-1+i]=float(np.log10(data[count,3]))
                mapp[ny-1+j,nx-1-i]=float(np.log10(data[count,3]))
                mapp[ny-1-j,nx-1-i]=float(np.log10(data[count,3]))
                mapt[ny-1+j,nx-1+i]=float(abs(data[count,5]))
                mapt[ny-1-j,nx-1+i]=float(abs(data[count,5]))
                mapt[ny-1+j,nx-1-i]=float(abs(data[count,5]))
                mapt[ny-1-j,nx-1-i]=float(abs(data[count,5]))
                mapm[ny-1+j,nx-1+i]=float(np.log10(data[count,2]))
                mapm[ny-1-j,nx-1+i]=float(np.log10(data[count,2]))
                mapm[ny-1+j,nx-1-i]=float(np.log10(data[count,2]))
                mapm[ny-1-j,nx-1-i]=float(np.log10(data[count,2]))
            if(float(np.log10(data[count,3]))>maxx):
                maxx=np.log10(float(data[count,3]))
        else:
            print "ERORR:  counjt: ", count, "x : ", x, "y: ", y, "data[0]", data[count,0], "data[1]: ", data[count,3]
            input ("Enter a numer")
        count  +=1
        
print "maximum polarization amount: ", maxx
  
  
#for i in range(2*nx-1):
#    for j in range(2*ny-1):
#        if(float(mapp[j,i])<-1.0):  mapp[j,i]=-5.0  
#        if(float(mapm[j,i])<1.5):  mapm[j,i]=0.0  


          
###============================= MAP polarization ================================================
tt=int(8)
v=np.zeros((tt))
plt.clf()
plt.plot(array[:,0],array[:,1],'k-',lw=0.9)
fig3=Figure()
plt.imshow(mapp,cmap='viridis',extent=(-1.0,1.0,1.0,-1.0),interpolation='nearest',aspect='auto')
plt.clim()
minn=-5.0## np.min(mapp[:,:])
maxx=np.max(mapp[:,:])
step=float((maxx-minn)/(tt-1.0));
for m in range(tt):
    v[m]=round(float(minn+m*step),1)
cb=plt.colorbar(orientation='vertical',shrink=1.0,pad=0.05,ticks=v)
plt.clim(v[0]-0.0005*step,v[tt-1]+0.0005*step)

plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16, rotation=0)
cb.ax.tick_params(labelsize=16)
plt.xlabel(r"$x(R_{\rm{E}})$",fontsize=18)
plt.ylabel(r"$y(R_{\rm{E}})$",fontsize=18)
plt.title(r"$q=1.0,~d(R_{E})=1.0,~\rho_{\star}= 0.007$", fontsize=18)

py.xlim([-1.0,1.0])
py.ylim([-1.0,1.0])
fig3=plt.gcf()
fig3.savefig("polmapint1_VB.jpg",dpi=200)
print(">>>>>>>>>>>>>>>>>>>>>>>> Map polarization was made <<<<<<<<<<<<<<<")






###============================= MAP polarization ================================================
tt=int(8)
v=np.zeros((tt))
plt.clf()
plt.plot(array[:,0],array[:,1],'k-',lw=0.9)
plt.imshow(mapt,cmap='viridis',extent=(-1.0,1.0,1.0,-1.0),interpolation='nearest',aspect='auto')
plt.clim()
minn=np.min(mapt[:,:])
maxx=np.max(mapt[:,:])
step=float((maxx-minn)/(tt-1.0));
for m in range(tt):
    v[m]=round(float(minn+m*step),1)
cb=plt.colorbar(orientation='vertical',shrink=1.0,pad=0.05,ticks=v)
plt.clim(v[0]-0.0005*step,v[tt-1]+0.0005*step)

cb.ax.tick_params(labelsize=16)
plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16, rotation=0)
plt.xlabel(r"$x(R_{\rm{E}})$",fontsize=18)
plt.ylabel(r"$y(R_{\rm{E}})$",fontsize=18)
plt.title(r"$q=1.0,~d(R_{E})=1.0,~\rho_{\star}= 0.007$", fontsize=18)

py.xlim([-1.0,1.0])
py.ylim([-1.0,1.0])
fig3=plt.gcf()
fig3.savefig("angmapint1_VB.jpg",dpi=200)
print(">>>>>>>>>>>>>>>>>>>>>>>> Map angle pol was made <<<<<<<<<<<<<<<")








###============================= MAP magnification ================================================
tt=int(8)
v=np.zeros((tt))
plt.clf()
plt.plot(array[:,0],array[:,1],'k-',lw=0.9)
plt.imshow(mapm,cmap='viridis',extent=(-1.0,1.0,1.0,-1.0),interpolation='nearest',aspect='auto')
plt.clim()
minn=np.min(mapm[:,:])
maxx=np.max(mapm[:,:])
step=float((maxx-minn)/(tt-1.0));
for m in range(tt):
    v[m]=round(float(minn+m*step),1)
cb=plt.colorbar(orientation='vertical',shrink=1.0,pad=0.05,ticks=v)
plt.clim(v[0]-0.0005*step,v[tt-1]+0.0005*step)

plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16, rotation=0)
cb.ax.tick_params(labelsize=16)
plt.xlabel(r"$x(R_{\rm{E}})$",fontsize=18)
plt.ylabel(r"$y(R_{\rm{E}})$",fontsize=18)
plt.title(r"$q=1.0,~d(R_{E})=1.0,~\rho_{\star}= 0.007$",fontsize=18)
py.xlim([-1.0,1.0])
py.ylim([-1.0,1.0])
fig3=plt.gcf()
fig3.savefig("magmapint1_VB.jpg",dpi=200)
print(">>>>>>>>>>>>>>>>>>>>>>>> Map angle pol was made <<<<<<<<<<<<<<<")





