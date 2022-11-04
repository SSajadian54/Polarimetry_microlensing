/* Bename Izad Taala *///// 15/5/97
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "VBBinaryLensingLibrary.h"
using namespace std;
// we want to minimize the time needed for calculating a polarimetry light curve 


///=================== Constants  ===========================================///
const double RA=180.0/M_PI;
const double KP=3.08568025*pow(10.,19); // in meter.
const double G= 6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity= 3.0*pow(10.0,8.0);//velosity of light
const double Msun=1.98892*pow(10.,30); //in [kg].
const double AU=1.4960*pow(10.0,11.0);
const double year=364.5; ///days
const double Rsun= 6.957*pow(10.0,8.0); ///solar radius [meter]
const double sigma=5.1*pow(10.0,-31.0);///Rylieh cross section [m^2]
const int NM= 2000; ///No.  steps over x or y !!!
const int nl=5000; 
const double drop= 0.001; /// it indicates the outer radius
const double ml=0.02; ///0.028; ///0.14-1.6*pow(10.0,-5.0)*7000.0;
const double ql=4.2*pow(10.0,-4.0);
const double c2= 0.032; 
const double MH= 3.35*pow(10.0,-27); ///Hydogen malequle mass in Kg [H2]
const double n0= double(0.0000003/MH);// rho0/<M>, rho0=0.1 kg/m^3

struct source {
    double Ds, Rs, rho_star, Rh, beta, tau, Renv;
    double xc[nl], yc[nl], x,y,dxy,R, alfa; 
    double tetas, t0, t;
    double tmin, tmax,proj; 
    double Ts,I0, gamal,limb;///main-sequence stars 
    int type; 
};
struct Binarylens{ 
     double u0, tE, q,dis,Dl, Mt,RE;
     double xls, x,y;
};
struct polarimetry{
   double Ii, Iq, Iu, Ip,S0;
   double Sii, Sqi, Sui; 
   double Sim[nl], Sqm[nl], Sum[nl];
   double poli, teti;  
       
};
///=============================================================================
void parameters(source & s, Binarylens & l); 
void Intensity(source & s, polarimetry & p);  
double TETA(double, double, int);

///=============================================================================
///==============================================================//
///                                                              //
///                  Main program                                //
///                                                              //
///==============================================================//
///=============================================================================
int main() {


    VBBinaryLensing *VBBL=new VBBinaryLensing;
    VBBL->Tol=1.e-3; // Setting accuracy
    source s; 
    Binarylens l; 
    polarimetry p; 
    int j,i; 
    double ymax, magni, polar, Tet, y1, y2; 


    parameters(s,l);

    
///=============================================================================  
    FILE* magpat;
    magpat=fopen("./files/magpat2.txt","w");
    double dl=0.000001; ///in Einstei radius
    double xlmin=-0.4,  xlmax=+0.2,  ylmin=0.05,  ylmax=0.6; 
    int imax=int((xlmax-xlmin)/dl)+10; 
    int jmax=int((ylmax-ylmin)/dl)+10; 
    cout<<"imax: "<<imax<<"\t jmax:  "<<jmax<<endl;
 
    double **Mag; 
    Mag=new double *[imax];
    for ( int i=0; i<imax; i++ ) Mag[i] = new double[jmax];
    cout<<"2D array was made!!!"<<endl;
    for(i=0; i<imax; ++i){  l.x= double(xlmin+ i*dl); 
    for(j=0; j<jmax; ++j){  l.y= double(ylmin+ j*dl);  
    Mag[i][j]=VBBL->BinaryMag0(l.dis,l.q ,l.x ,l.y);
    fprintf(magpat,"%.5lf    %.5lf     %.10lf \n",l.x,l.y,Mag[i][j]);}
    if(i%300==0) cout<<"i: "<<i<<"\t x: "<<l.x<<"\t j: "<<j<<endl;   }
    fclose(magpat);
    delete VBBL; 
    cout<<"*********** magnification pattern was made ******************"<<endl;
    //int uue; cin>>uue; 



///=============================================================================
    FILE* poli2;
    poli2=fopen("./files/Intensitypat2.txt","w");
    s.dxy= double(dl*l.RE*AU/Rsun); 
    cout<<"No. devisionof Rs: "<<s.Rs/s.dxy<<endl;
    p.Sqi=p.Sui=p.Sii=0.0; j=0; 
    for(s.x=0.0; s.x<= s.Renv; s.x += s.dxy){
    ymax=sqrt(s.Renv*s.Renv-s.x*s.x); 
    for(s.y=0.0; s.y<=ymax;    s.y += s.dxy){
    ++j; }}
    cout<<"dxy: "<<s.dxy<<"\t total number of intensity calculations: "<<j<<endl; 
    double Ii[j+5]={0.0}; double Iq[j+5]={0.0};  double Iu[j+5]={0.0};     
    int JJ= j; j=0; 
    for(s.x=0.0; s.x<= s.Renv; s.x += s.dxy){
    ymax=sqrt(s.Renv*s.Renv-s.x*s.x); 
    for(s.y=0.0; s.y<=ymax;    s.y += s.dxy){
    s.R=sqrt(s.x*s.x+s.y*s.y); 
    s.alfa=TETA(s.x,s.y,0);
    Intensity(s,p);
    Ii[j]=p.Ii;   
    Iq[j]=p.Iq;
    Iu[j]=p.Iu; 
    p.Sqi += 4.0*p.Iq;
    p.Sui += 0.0*p.Iu;
    p.Sii += 4.0*p.Ii;
    p.poli=double(p.Ip*100.0/p.Ii);// in per cent 
    p.teti= TETA(p.Iq,p.Iu,1);
    fprintf(poli2,"%.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %g  %g  %g\n",s.x,s.y,s.alfa*RA,s.R,p.poli,p.teti,p.Iq,p.Iu,p.Ii); 
    ++j; 
    }}
    p.poli= double(sqrt(p.Sqi*p.Sqi+ p.Sui*p.Sui)*100.0/p.Sii);
    p.teti= TETA(p.Sqi,p.Sui,1); 
    cout<<"Total intrinsic polarization[%]: "<<p.poli<<endl;
    cout<<"Tet(degree): "<<p.teti<<endl;
    cout<<"*********** Intensity pattern was made **********************"<<endl;
    fclose(poli2);
    //int uuw; cin>>uuw; 




///=============================================================================
    FILE* poli;
    poli=fopen("./files/polmag2.txt","w");
    for (int k=0; k<nl; ++k) {
    p.Sim[i]=p.Sqm[i]=p.Sum[i]=0.0;
    int l=0; 
    for(s.x=0.0; s.x<= s.Renv; s.x += s.dxy){
    ymax=sqrt(s.Renv*s.Renv-s.x*s.x); 
    for(s.y=0.0; s.y<=ymax;    s.y += s.dxy){
  
///============================================
    y1= double(s.x*s.proj+s.xc[k]); 
    y2= double(s.y*s.proj+s.yc[k]);
    i=int((y1-xlmin)/dl);  j=int((y2-ylmin)/dl);
    if(i>=0 and i<imax and j>=0 and j<jmax){
    p.Sim[k] += Mag[i][j]* Ii[l];
    p.Sqm[k] += Mag[i][j]* Iq[l]; 
    p.Sum[k] += Mag[i][j]* Iu[l];}
    else{ cout<<"ERORR: i: "<<i<<"\t j: "<<j<<"\t y1:  "<<y1<<"\t y2: "<<y2<<endl;  int iie;  cin>>iie; }
///============================================
    y1= double(-s.x*s.proj+s.xc[k]); 
    y2= double(s.y*s.proj+s.yc[k]);
    i=int((y1-xlmin)/dl);  j=int((y2-ylmin)/dl);
    if(i>=0 and i<imax and j>=0 and j<jmax){
    p.Sim[k] += Mag[i][j]* Ii[l];
    p.Sqm[k] += Mag[i][j]* Iq[l]; 
    p.Sum[k] += -Mag[i][j]* Iu[l];}
    else{ cout<<"ERORR: i: "<<i<<"\t j: "<<j<<"\t y1:  "<<y1<<"\t y2: "<<y2<<endl;  int iie;  cin>>iie; }
///============================================
    y1= double(s.x*s.proj+s.xc[k]); 
    y2= double(-s.y*s.proj+s.yc[k]);
    i=int((y1-xlmin)/dl);  j=int((y2-ylmin)/dl);
    if(i>=0 and i<imax and j>=0 and j<jmax){
    p.Sim[k] += Mag[i][j]* Ii[l];
    p.Sqm[k] += Mag[i][j]* Iq[l]; 
    p.Sum[k] += -Mag[i][j]* Iu[l];}
    else{ cout<<"ERORR: i: "<<i<<"\t j: "<<j<<"\t y1:  "<<y1<<"\t y2: "<<y2<<endl;  int iie;  cin>>iie; }
///============================================
    y1= double(-s.x*s.proj+s.xc[k]); 
    y2= double(-s.y*s.proj+s.yc[k]);
    i=int((y1-xlmin)/dl);  j=int((y2-ylmin)/dl);
    if(i>=0 and i<imax and j>=0 and j<jmax){
    p.Sim[k] += Mag[i][j]* Ii[l];
    p.Sqm[k] += Mag[i][j]* Iq[l]; 
    p.Sum[k] += Mag[i][j]* Iu[l]; }
    else{ cout<<"ERORR: i: "<<i<<"\t j: "<<j<<"\t y1:  "<<y1<<"\t y2: "<<y2<<endl;  int iie;  cin>>iie; } ++l; }}
    if(l != JJ ) {cout<<"BIG ERROR l: "<<l<<"\t JJ: "<<JJ<<endl;  int uue;  cin>>uue; }
    magni= double(p.Sim[k]/p.Sii); 
    polar= double(sqrt(p.Sqm[k]*p.Sqm[k]+p.Sum[k]*p.Sum[k])*100.0/p.Sim[k]);///[%]
    Tet= TETA(p.Sqm[k],p.Sum[k],1);
    fprintf(poli,"%d  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %g  %g  %g  %g\n",
    i,s.xc[k],s.yc[k],magni,polar,Tet,p.Sim[k],p.Sqm[k],p.Sum[k],p.Sii);
    cout<<"================================================================"<<endl;
    cout<<"step: "<<i<<"\t xc: "<<s.xc[k]<<"\t yc: "<<s.yc[k]<<endl;
    cout<<"magnification: "<<magni<<"\t pol[%]: "<<polar<<"\t Tet(degree): "<<Tet<<endl;
    cout<<"Si: "<<p.Sim[k]<<"\t Sq: "<<p.Sqm[k]<<"\t Su: "<<p.Sum[k]<<endl;
    cout<<"================================================================="<<endl;} 
    fclose(poli); 
    return(0);
}
///=============================================================================
void parameters(source & s, Binarylens & l)
{
    l.u0=0.3; 
    l.tE=30.0;//day 
    l.q=1.0;   
    l.dis=1.0; 
    l.Dl= 6.5;  
    l.Mt=0.3* (1.0 +l.q); ///in solar mass
    s.Ds=8.0; 
    s.Rs=3.0;///in [solar unit]
    s.t0=0.0; 
    s.tetas=(0.0)/RA;///radian 
    s.tmin=-0.4*l.tE+s.t0; 
    s.tmax= 0.2*l.tE+s.t0;
    l.xls=double(l.Dl/s.Ds); 
    l.RE= sqrt(4.0*G*l.Mt*Msun*s.Ds*KP)*sqrt(l.xls*(1.0-l.xls))/velocity/AU; 
    s.rho_star= s.Rs*Rsun*l.Dl/(s.Ds*l.RE*AU); 
    s.proj=double(l.xls*Rsun/l.RE/AU); 
    s.I0=1.0;

    s.type=2; 
    s.Ts= 7700.0;///surface temperature 
    if(s.type==3){
    s.limb=0.7; 
    s.Rh= 1.1*s.Rs; ///[in solar unit] 
    s.beta= 2.0;///power of density profile
    s.tau= n0*sigma*s.Rh*Rsun/(s.beta-1.0); 
    s.Renv=s.Rh*pow(drop,-1.0/s.beta);}

    else {
    s.Rh=s.Rs; 
    s.Renv=s.Rs;
    s.limb= 0.64;///early-type stars 
    if(s.Ts>=6200.0)                      s.gamal=0.46;
     else if(s.Ts>=5700.0 && s.Ts<6200.0) s.gamal=0.53;
     else if(s.Ts>=5100.0 && s.Ts<5700.0) s.gamal=0.54;
     else if(s.Ts>=4800.0 && s.Ts<5100.0) s.gamal=0.59;
     else if(s.Ts>=4600.0 && s.Ts<4800.0) s.gamal=0.61;
     else if(s.Ts>=4400.0 && s.Ts<4600.0) s.gamal=0.63;
     else if(s.Ts<=4400.0 )               s.gamal=0.65;
     if(s.Ts>=7500.0) s.type=2;
     else s.type=1; }

    double dt= double(s.tmax-s.tmin)/nl/1.0;
    for(int i=0; i<nl; ++i ){
    s.t=s.tmin + dt*i;
    double pt= double(s.t-s.t0)/l.tE;
    s.xc[i]=pt*cos(s.tetas) - l.u0*sin(s.tetas);
    s.yc[i]=pt*sin(s.tetas) + l.u0*cos(s.tetas); }

    cout<<"Renv: "<<s.Renv<<"\t Rh: "<<s.Rh<<endl;
    cout<<"tau: "<<s.tau<<"\t rho_star: "<<s.rho_star<<endl;
   
}
///=============================================================================
double TETA(double x, double y,int flag) {

double alfa=0.0;

   if(y>=0.0 && x>0.0)       alfa=atan(fabs(y)/fabs(x));
   else if(y>0.0  && x<=0.0) alfa=atan(fabs(x)/fabs(y))+M_PI/2.0;
   else if(y<=0.0 && x<0.00) alfa=atan(fabs(y)/fabs(x))+M_PI;
   else if(y<0.0 &&  x>=0.0) alfa=atan(fabs(x)/fabs(y))+3.0*M_PI/2.0;
   if(x==0.0 && y==0.0) alfa=0.0;
   if(x==0.0 && y>0.0 ) alfa=M_PI/2.0;
   if(x==0.0 && y<0.0 ) alfa=3.0*M_PI/2.0;
   if(y==0.0 && x>0.0 ) alfa=0.0;
   if(y==0.0 && x<0.0 ) alfa=M_PI;///in radian
   if(flag==1) {
   alfa=alfa*0.5*RA;///in degree;فقط باید در ۱/۲ ضرب شود وهیچ تغییر دیگری نیاز ندارد
   if(alfa>360.0) alfa=alfa-360.0;
   if(alfa>360.0) alfa=alfa-360.0;}
   return(alfa);
}
///=============================================================================
void Intensity(source & s, polarimetry & p){

   p.Ii=0.0; 
   p.Iq=p.Iu=p.Ip=0.0;
   double mu; 
   if(s.R<s.Rs) mu=sqrt(fabs(1.0-s.R*s.R/s.Rs/s.Rs));
   else mu=0.0; 


if(s.type==1){/// main sequence stars     
    p.Ip=ql*s.I0*(1.0-mu*mu)/(mu+ml);
    p.Ii=s.I0*(1.0-s.gamal*(1.0-1.5*mu));
    p.Iq=-p.Ip*cos(2.0*s.alfa);
    p.Iu=-p.Ip*sin(2.0*s.alfa); }



else if(s.type==2){///early-type stars
   p.Ip=c2*s.I0*(1.0-mu);
   p.Ii=s.I0*(1.0-s.limb*(1.0-mu));
   p.Iq=-p.Ip*cos(2.0*s.alfa);
   p.Iu=-p.Ip*sin(2.0*s.alfa);}



else if(s.type==3){////red giant stars, Simmons et al. 2002   
   double z,zmin, J, K,r,density; 
   if(s.R> s.Rh)  zmin=0.0;
   if(s.R< s.Rh)  zmin= sqrt(fabs(s.Rh*s.Rh-s.R*s.R)); 
   if(s.R==s.Rh)  zmin=0.0; 
   z=zmin;   
   double dz=0.01*s.Rs;///[in solar unit]
   int flag=-1;  
   //double frac=0.0;  ///double first=n0;//fabs(n0*pow(s.Rh/sqrt(s.R*s.R+zmin*zmin),s.beta));
   int step=0;   



   do{
   flag=-1; 
   r= sqrt(fabs(s.R*s.R+ z*z) );

   if(r>=s.Rh) {
   ++step;    
   J= 0.5*s.I0*(1.0-sqrt(1.0-s.Rs*s.Rs/r/r));
   K= s.I0*(1.0-pow(1.0-s.Rs*s.Rs/r/r,1.5))/6.0; 
   density=fabs(n0*pow(s.Rh/r,s.beta)); 
   p.Ii+=  density*fabs( (3.0*J-K)+(3.0*K-J)*pow(z/r,2.0));/// in the unit I0
   p.Ip+=  density*fabs( (3.0*K-J)*pow(s.R/r,2.0) ); 
   if(double(density/n0)<=drop ) flag=1; 
   /*cout<<"********************************************************"<<endl;
   cout<<"step: "<<step<<endl; 
   cout<<"r/Rs: "<<r/s.Rs<<"\t density/n0: "<<density/n0<<endl; ///<"\t first: "<<first<<endl;  
   cout<<"z: "<<z<<"\t J: "<<J<<"\t K: "<<K<<endl;     int uue;  cin>>uue; */
   }//end of if
   if(r==0.0 || r<s.Rs || z>r || r<0.0){
   cout<<"ERROR:  r: "<<r<<"\t z: "<<z<<"\t R: "<<s.R<<"\t x: "<<s.x<<"\t Y: "<<s.y<<endl; int uue; cin>>uue; }
   z += dz; 
   }while(flag<0); 
   if(s.R>=s.Rs){p.Ii=2.0*p.Ii;    p.Ip=2.0*p.Ip;}
   p.Ii= s.I0 + 3.0*sigma*p.Ii*dz*Rsun/8.0;
   p.Iq=-p.Ip*cos(2.0*s.alfa)*dz*Rsun*3.0*sigma/8.0; 
   p.Iu=-p.Ip*sin(2.0*s.alfa)*dz*Rsun*3.0*sigma/8.0;



   /*cout<<"===========INTENSITY FUNCTION================================="<<endl;
   cout<<"zmin: "<<zmin<<"\t alfa: "<<s.alfa<<endl;
   cout<<"R: "<<s.R<<"\t x: "<<s.x<<"\t y: "<<s.y<<endl;
   cout<<"poli[%]: "<<sqrt(p.Iq*p.Iq+ p.Iu*p.Iu)*100.0/p.Ii<<"\t  No. step: "<<step<<endl;
   cout<<"p.Ii: "<<p.Ii<<"\t p.Iq: "<<p.Iq<<"\t Iu: "<<p.Iu<<endl; 
   cout<<"flag: "<<flag<<"\t break at: "<<density/n0<<"\t r/Renv: "<<r/s.Renv<<endl;
  // int yye; cin>>yye; 
   cout<<"=============================================================="<<endl;*/

}///end of red giant stars 

}
////============================================================================

