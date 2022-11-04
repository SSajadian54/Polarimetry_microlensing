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




rhostar=float(0.00705783509461056172)
RA=float(np.pi/180.0)
n=21546
array=np.zeros((800,2))
array=np.loadtxt("caustic_wide1.txt") 
data=np.zeros((n,7))
data=np.loadtxt("Map_pol_wide_VB.txt") 


###############################################################################
count=int(0)
nx=int(378)
ny=int(57)
step=float(0.75)
mapp=np.zeros((ny*2-1,nx))
mapt=np.zeros((ny*2-1,nx))
mapm=np.zeros((ny*2-1,nx))
maxx=float(-10.0)
maxm=float(1.0); 
for i in range(nx):
    x=-0.8+i*rhostar*step;
    for j in range(ny):
        y= 0.0 + j*rhostar*step;
        if(abs(x-data[count,0])<float(step*rhostar*0.1) and abs(y-data[count,1])<float(step*rhostar*0.1)):
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
            if(float(data[count,2])>maxm):
                maxm=float(data[count,2]);  
        else:
            print "ERORR:  count: ", count, "x : ", x, "y: ", y, "data[0]", data[count,0], "data[1]: ", data[count,3]
            input ("Enter a number")
        count  +=1
print "maximum polarization amount: ", maxx
print "maximum magnification amount: ", maxm    

   

###============================= MAP polarization ================================================
tt=int(8)
v=np.zeros((tt))
plt.clf()
plt.plot(array[:400,0],array[:400,1],'k-',lw=0.9)
plt.plot(array[400:,0],array[400:,1],'k-',lw=0.9)
plt.imshow(mapp,cmap='viridis',extent=(-0.8,1.2,0.3,-0.3),interpolation='nearest',aspect='auto')
plt.clim()
minn=-5.0## np.min(mapp[:,:])
maxx=np.max(mapp[:,:])
step=float((maxx-minn)/(tt-1.0));
for m in range(tt):
    v[m]=round(float(minn+m*step),1)
cb=plt.colorbar(orientation='horizontal',shrink=1.0,pad=0.2,ticks=v)
plt.clim(v[0]-0.0005*step,v[tt-1]+0.0005*step)

plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16, rotation=0)
cb.ax.tick_params(labelsize=16)
plt.xlabel(r"$x(R_{\rm{E}})$",fontsize=18)
plt.ylabel(r"$y(R_{\rm{E}})$",fontsize=18)
plt.title(r"$q=0.7,~d(R_{E})=2,~\rho_{\star}=0.008$",fontsize=18)

py.xlim([-0.8,1.2])
py.ylim([-0.3,0.3])
fig3=plt.gcf()
fig3.savefig("polmapwide_VB.jpg",dpi=200)
print(">>>>>>>>>>>>>>>>>>>>>>>> Map polarization was made <<<<<<<<<<<<<<<")


###============================= MAP polarization ================================================
tt=int(8)
v=np.zeros((tt))
plt.clf()
plt.plot(array[:400,0],array[:400,1],'k-',lw=0.9)
plt.plot(array[400:,0],array[400:,1],'k-',lw=0.9)
plt.imshow(mapt,cmap='viridis',extent=(-0.8,1.2,0.3,-0.3),interpolation='nearest',aspect='auto')
plt.clim()
minn=np.min(mapt[:,:])
maxx=np.max(mapt[:,:])
step=float((maxx-minn)/(tt-1.0));
for m in range(tt):
    v[m]=round(float(minn+m*step),1)
cb=plt.colorbar(orientation='horizontal',shrink=1.0,pad=0.2,ticks=v)
plt.clim(v[0]-0.0005*step,v[tt-1]+0.0005*step)


plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16, rotation=0)
cb.ax.tick_params(labelsize=16)
plt.xlabel(r"$x(R_{\rm{E}})$",fontsize=18)
plt.ylabel(r"$y(R_{\rm{E}})$",fontsize=18)
plt.title(r"$q=0.7,~d(R_{E})=2,~\rho_{\star}=0.008$",fontsize=18)

py.xlim([-0.8,1.2])
py.ylim([-0.3,0.3])
fig3=plt.gcf()
fig3.savefig("angmapwide_VB.jpg",dpi=200)
print(">>>>>>>>>>>>>>>>>>>>>>>> Map angle pol was made <<<<<<<<<<<<<<<")


###============================= MAP magnification ================================================
tt=int(8)
v=np.zeros((tt))
plt.clf()
plt.plot(array[:400,0],array[:400,1],'k-',lw=0.9)
plt.plot(array[400:,0],array[400:,1],'k-',lw=0.9)
plt.imshow(mapm,cmap='viridis',extent=(-0.8,1.2,0.3,-0.3),interpolation='nearest',aspect='auto')
plt.clim()
minn=np.min(mapm[:,:])
maxx=np.max(mapm[:,:])
step=float((maxx-minn)/(tt-1.0));
for m in range(tt):
    v[m]=round(float(minn+m*step),1)
cb= plt.colorbar(orientation='horizontal',shrink=1.0,pad=0.2,ticks=v)
plt.clim(v[0]-0.0005*step,v[tt-1]+0.0005*step)


plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16, rotation=0)
cb.ax.tick_params(labelsize=16)
plt.xlabel(r"$x(R_{\rm{E}})$",fontsize=18)
plt.ylabel(r"$y(R_{\rm{E}})$",fontsize=18)
plt.title(r"$q=0.7,~d(R_{E})=2,~\rho_{\star}=0.008$",fontsize=18)

py.xlim([-0.8,1.2])
py.ylim([-0.3,0.3])
fig3=plt.gcf()
fig3.savefig("magmapwide_VB.jpg",dpi=200)
print(">>>>>>>>>>>>>>>>>>>>>>>> Map angle pol was made <<<<<<<<<<<<<<<")





