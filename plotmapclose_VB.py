import numpy as np
import matplotlib.pyplot as plt
import pylab as py
from matplotlib import rcParams
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



rhostar=float(0.0060625415798480344550803)
RA=float(np.pi/180.0)
array=np.zeros((800,2))
array=np.loadtxt("caustic_close1.txt")
n=int(59680)##int(129800) 
data=np.zeros((n,7))
data=np.loadtxt("Map_pol_close1.txt") 
###############################################################################

ny=int(373)
nx=int(160)

count=int(0)
mapp=np.zeros((ny*2-1,nx)) 
mapt=np.zeros((ny*2-1,nx))
mapm=np.zeros((ny*2-1,nx))
maxx=float(-10.0)
nmax=int(0)
for i in range(nx):
    x=-0.3+i*rhostar*0.62;
    for j in range(ny):
        y= 0.0 + j*rhostar*0.62;
        if(abs(x-data[count,0])<float(0.1*rhostar) and abs(y-data[count,1])<float(0.1*rhostar)):
            if(j==0):
                mapp[ny-1,i]=float(np.log10(data[count,3]))
                mapt[ny-1,i]=float(abs(data[count,5]))
                mapm[ny-1,i]=float(np.log10(data[count,2]))
            else:
                mapp[j+ny-1,i]=float(np.log10(data[count,3]))
                mapp[ny-1-j,i]=float(np.log10(data[count,3]))
                mapt[j+ny-1,i]=float(abs(data[count,5]))
                mapt[ny-1-j,i]=float(abs(data[count,5]))
                mapm[j+ny-1,i]=float(np.log10(data[count,2]))
                mapm[ny-1-j,i]=float(np.log10(data[count,2]))
            if(float(np.log10(data[count,3]))>maxx):
                maxx=np.log10(float(data[count,3]))
                nmax=count
        else:
            print "ERORR:  count: ", count, "x : ", x, "y: ", y, "data[0]", data[count,0], "data[1]: ", data[count,3]
            input ("Enter a number")
        count  +=1
print "maximum polarization amount: ", maxx, "polarization itself: ", pow(10.0,maxx)
print "Counter: ", nmax, "all: ",  data[nmax,:]
    



  
#for i in range(nx):
#    for j in range(2*ny-1):
#        if(float(mapp[j,i])<0.05):  mapp[j,i]=-5.0  
#        if(float(mapm[j,i])<1.5):  mapm[j,i]=0.0  


          
###============================= map polarization ================================================
tt=int(8)
v=np.zeros((tt))
plt.clf()
plt.plot(array[:200,0],array[:200,1],'k-',lw=0.9)
plt.plot(array[200:400,0],array[200:400,1],'k-',lw=0.9)
plt.plot(array[400:,0],array[400:,1],'k-',lw=0.9)
fig3=Figure()
#plt.figure(figsize=(5.0,5.0),facecolor='g')
plt.imshow(mapp,cmap='viridis',extent=(-0.3,0.3,1.4,-1.4),interpolation='nearest',aspect='auto')
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
plt.title(r"$q=0.6,~d(R_{E})=0.7,~\rho_{\star}=0.006$", fontsize=18)

py.xlim([-0.3,0.3])
py.ylim([-1.4,1.4])
fig3=plt.gcf()
fig3.savefig("polmapclose_VB.jpg",dpi=200)
print(">>>>>>>>>>>>>>>>>>>>>>>> Map polarization was made <<<<<<<<<<<<<<<")






###============================= map angle polarization ================================================
tt=int(8)
v=np.zeros((tt))
plt.clf()
plt.plot(array[:200,0],array[:200,1],'k-',lw=0.9)
plt.plot(array[200:400,0],array[200:400,1],'k-',lw=0.9)
plt.plot(array[400:,0],array[400:,1],'k-',lw=0.9)
plt.imshow(mapt,cmap='viridis',extent=(-0.3,0.3,1.4,-1.4),interpolation='nearest',aspect='auto')
plt.clim()
minn=np.min(mapt[:,:])
maxx=np.max(mapt[:,:])
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
plt.title(r"$q=0.6,~d(R_{E})=0.7,~\rho_{\star}=0.006$", fontsize=18)

py.xlim([-0.3,0.3])
py.ylim([-1.4,1.4])
fig3=plt.gcf()
fig3.savefig("angmapclose_VB.jpg",dpi=200)
print(">>>>>>>>>>>>>>>>>>>>>>>> Map angle pol was made <<<<<<<<<<<<<<<")








###============================= map magnification ================================================
tt=int(8)
v=np.zeros((tt))
plt.clf()
plt.plot(array[:200,0],array[:200,1],'k-',lw=0.9)
plt.plot(array[200:400,0],array[200:400,1],'k-',lw=0.9)
plt.plot(array[400:,0],array[400:,1],'k-',lw=0.9)
plt.imshow(mapm,cmap='viridis',extent=(-0.3,0.3,1.4,-1.4),interpolation='nearest',aspect='auto')
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
plt.title(r"$q=0.6,~d(R_{E})=0.7,~\rho_{\star}=0.006$", fontsize=18)

py.xlim([-0.3,0.3])
py.ylim([-1.4,1.4])
fig3=plt.gcf()
fig3.savefig("magmapclose_VB.jpg",dpi=200)
print(">>>>>>>>>>>>>>>>>>>>>>>> Map angle pol was made <<<<<<<<<<<<<<<")


plt.clf()
plt.plot(array[:,0],array[:,1],'ko',ms=0.5)
#plt.plot(0.18633709,0.93193389,'ro',ms=3.5)
plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16, rotation=0)
plt.xlabel(r"$x(R_{\rm{E}})$",fontsize=18)
plt.ylabel(r"$y(R_{\rm{E}})$",fontsize=18)
plt.title(r"$q=0.6,~d(R_{E})=0.7,~\rho_{\star}=0.006$", fontsize=18)
py.xlim([-0.3,0.3])
py.ylim([-1.4,1.4])
fig3=plt.gcf()
#fig3=plt.gcf().gca().add_artist(circle1)
fig3.savefig("test.jpg",dpi=200)




