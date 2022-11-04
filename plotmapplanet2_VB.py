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




rhostar=float(0.0076681925206768678307534)##0.00766819252067686783) 
RA=float(np.pi/180.0)
n=int(121907)
nc=int(800)
array=np.zeros((nc,2))
array=np.loadtxt("caustic_planet.txt") 
data=np.zeros((n,7))
data=np.loadtxt("Map_pol_planet2_VBc.txt") 

###############################################################################

ny=int(101)
nx=int(1207)
count=int(0)
mapp=np.zeros((ny*2-1,nx))
mapt=np.zeros((ny*2-1,nx))
mapm=np.zeros((ny*2-1,nx))
maxx=float(-10.0)
for i in range(nx):
    x=-0.25+i*rhostar*0.032415;
    for j in range(ny):
        y= 0.0 + j*rhostar*0.032415;
        if(abs(x-data[count,0])<float(0.2*rhostar) and abs(y-data[count,1])<float(0.2*rhostar)):
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
        else:
            print ("ERORR:  count: ", count, "x : ", x, "y: ", y, "data[0]", float(data[count,0]), "data[1]: ", float(data[count,1]) )
            input("Enter a number")
        count  +=1
print "maximum polarization amount: ", maxx







    
dds=open("ring2.dat","w")    
dds.close()    
pp=np.zeros((3)) 
for i in range(nx):
    for j in range(2*ny-1):        
        if(float(mapp[j,i])>-0.2): 
            pp[0]=-0.25+i*rhostar*0.072415;
            pp[1]= 0.0 + (j-ny+1.0)*rhostar*0.072415;
            pp[2]= float(mapp[j,i])
            dds=open("ring2.dat","a")
            np.savetxt(dds,pp.reshape((1,3)),fmt="%e    %e   %e")
            dds.close()



        

###============================= MAP polarization ================================================
tt=int(8)
v=np.zeros((tt))
plt.clf()
plt.plot(array[:200,0],array[:200,1],'k-',lw=0.9)
plt.plot(array[200:400,0],array[200:400,1],'k-',lw=0.9)
plt.plot(array[400:,0],array[400:,1],'k-',lw=0.9)
plt.imshow(mapp,cmap='viridis',extent=(-0.25,0.05,0.025,-0.025),interpolation='nearest',aspect='auto')
plt.clim()
minn=np.min(mapp[:,:])
maxx=np.max(mapp[:,:])
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
plt.title(r"$q=0.0001,~d(R_{E})=0.9,~\rho_{\star}=0.0077$",fontsize=18)


py.xlim([-0.25,0.05])
py.ylim([-0.025,0.025])
fig3=plt.gcf()
fig3.savefig("polmapplanet2_VB.jpg",dpi=200)
print(">>>>>>>>>>>>>>>>>>>>>>>> Map polarization was made <<<<<<<<<<<<<<<")






###============================= MAP angle  ================================================
tt=int(8)
v=np.zeros((tt))
plt.clf()
plt.plot(array[:200,0],array[:200,1],'k-',lw=0.9)
plt.plot(array[200:400,0],array[200:400,1],'k-',lw=0.9)
plt.plot(array[400:,0],array[400:,1],'k-',lw=0.9)
plt.imshow(mapt,cmap='viridis',extent=(-0.25,0.05,0.3,-0.3),interpolation='nearest',aspect='auto')
plt.clim()
minn=np.min(mapt[:,:])
maxx=np.max(mapt[:,:])
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
plt.title(r"$q=0.0001,~d(R_{E})=0.9,~\rho_{\star}=0.0077$",fontsize=18)

py.xlim([-0.25,0.05])
py.ylim([-0.025,0.025])
fig3=plt.gcf()
fig3.savefig("angmapplanet2_VB.jpg",dpi=200)
print(">>>>>>>>>>>>>>>>>>>>>>>> Map angle pol was made <<<<<<<<<<<<<<<")








###============================= magnification map ================================================
tt=int(8)
v=np.zeros((tt))
plt.clf()
plt.plot(array[:200,0],array[:200,1],'k-',lw=0.9)
plt.plot(array[200:400,0],array[200:400,1],'k-',lw=0.9)
plt.plot(array[400:,0],array[400:,1],'k-',lw=0.9)
plt.imshow(mapm,cmap='viridis',extent=(-0.25,0.05,0.025,-0.025),interpolation='nearest',aspect='auto')
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
plt.title(r"$q=0.0001,~d(R_{E})=0.9,~\rho_{\star}=0.0077$",fontsize=18)

py.xlim([-0.25,0.05])
py.ylim([-0.025,0.025])
fig3=plt.gcf()
fig3.savefig("magmapplanet2_VB.jpg",dpi=200)
print(">>>>>>>>>>>>>>>>>>>>>>>> Magnification Map was made <<<<<<<<<<<<<<<")





