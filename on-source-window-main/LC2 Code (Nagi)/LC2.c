#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define day 86400.0		/*1 day in seconds*/
#define xmin 0.4		/*Minimum ionization radius*/	
#define dt 1.0			/*Time step (s)*/
#define Eni 3.89e10		/*Ni energy generation rate (erg/g/s)*/
#define Eco 6.8e9		/*Co energy generation rate (erg/g/s)*/
#define c 2.99e10		/*Speed of light (cm/s)*/
#define Msol 1.989e33		/*Solar mass (g)*/



double E,Mni,Ek,F,L,T0,T,Tc,Tion,Lion,Ith,IM,INi,Eth,Lsn;
double t,td,th,tmax,sigma,sigma1,p1,p2,p3,p4,p5,g,t0;
double xh,xi,dxi,Xni,Xco,dXni,dXco,logL,a,Z,A,Ag,Ep,s;
double kappa,M,v,R0,R,rho,ri,dr,Q,t0,tni,tco,tp,Em,Em1;

int data()			 /*Reading input file*/
{
FILE *f;

f=fopen("parameters.inp","rt");

   fscanf(f,"%lf",&R0);   	  /*Initial radius (cm)*/
   fscanf(f,"%lf",&M);    	  /*Ejecta mass (M_sun)*/
   fscanf(f,"%lf",&Tion); 	  /*Ionization/recombination temperature (K)*/
   fscanf(f,"%lf",&Mni);  	  /*Initial nickel mass (M_sun) */
   fscanf(f,"%lf",&Ek); 	  /*Initial kinetic energy (1e51 erg)*/ 
   fscanf(f,"%lf",&E); 		  /*Intital thermal energy (1e51 erg)*/    
   fscanf(f,"%lf",&a);            /*Exponential density profile exponent*/
   fscanf(f,"%lf",&s);		  /*Power-low density profile exponent*/
   fscanf(f,"%lf",&kappa);     	  /*Thomson scattering opacity (cm^2/g)*/
   fscanf(f,"%lf",&Ep);		  /*Initial magnetar rotational energy (erg)*/
   fscanf(f,"%lf",&tp);		  /*Characteristic time scale of magnetar spin-down (d)*/
   fscanf(f,"%lf",&Ag);		  /*Gamma-leak (d^2)*/
   fscanf(f,"%lf",&tmax);	  /*Final epoch (day)*/

M=Msol*M;	 		
Mni=Msol*Mni;
E=1e51*E;
Ek=1e51*Ek;
Ep=1e51*Ep;
tp=tp*86400.0;
Ag=Ag*pow(86400.0,2);		
T0=pow(E*M_PI/(4*pow(R0,3.0)*7.57e-15),0.25); 

if(s==3.0 || s==5.0){printf("Parameter error! The n = 3.0 and n = 5.0 are forbidden!\n"); exit(1);}

fclose(f);
}

double psi(double x)			/*Temperature profile assuming uniform density*/
{
double z;
if(x==1.0) z=0.0; else if(x>0) z=sin(M_PI*x)/(M_PI*x); else z=1.0;
return z;
}

double eta(double x, double a1)		/*Exponential density profile*/
{
double z;
if(x<xmin) z=1; else z=exp(-a1*(x-xmin));
return z;
}

double theta(double x, double a2)	/*Power-low density profile*/
{
double z;
if(x<xmin) z=1; else z=pow(x/xmin,-a2);
return z;
}

double temp(double T1,double y0)	/*Find ionization radius*/
{
double y,dy,zone;
y=y0;
dy=1e-9;
T=T1*pow(psi(1.0),0.25);

	while(y>=xmin && T<Tion)
	{
	   T=T1*pow(psi(y),0.25);
	   y=y-dy;
	} 

if(y<y0) zone=y+0.5*dy; else zone=y0;
return zone;	
}  

double IM_int(double b, double a1, double a2)
{
double sum=0.0,x=0.0,dx=1e-7;
while(x<b)
{
sum=sum+(eta(x,a1)*pow(x,a2)+eta(x+dx,a1)*pow(x+dx,a2))*dx*0.5; 	/*Trapesoid rule*/
x=x+dx;
}
return sum;
}

double Ith_int(double b)
{
double sum=0.0,x=0.0,dx=1e-7;

while(x<b)
{
sum=sum+(psi(x)*x*x+psi(x+dx)*pow(x+dx,2.0))*dx*0.5;
x=x+dx;
}
return sum;
}

double series(double x)			/*Taylor-series of sin(pi*x)/(pi*x)*/
{
double z;
z=1-pow(M_PI*x,2)/6.0+pow(M_PI*x,4)/120.0-pow(M_PI*x,6)/5040.0+pow(M_PI*x,8)/362880.0;
return z;
}

int main()
{
int j=0,i=0;
double dF1,dF2,dF3,dF4;
double opac,f,g1;

tni=8.8*day;				/*Ni decay time (s)*/
tco=111.3*day;				/*Co decay time (s)*/
Z=1.0;					/*Average atomic charge in ejecta*/
A=1.0;					/*Average atomic mass in ejecta*/

data();

  tmax=tmax*day; 			/*Final epoch (s)*/
  t=0.0; F=1.0;
  xh=1.0; xi=1.0; 
  Xni=1.0; Xco=0.0;
  Q=1.6e13*Z/A*pow(Z,1.333333);		/*Ionization energy release*/
  //Q = 3.22e13;
  
  Ith=Ith_int(1.0);
  IM=IM_int(1.0,a,2.0);
  
  
  if(s==0.0) f=IM_int(1.0,a,2.0);
  else f=(3.0*pow(xmin,s)-s*pow(xmin,3.0))/(3.0*(3.0-s));
  
  if(s==0.0) g1=IM_int(1.0,a,4.0);
  else g1=(5.0*pow(xmin,s)-s*pow(xmin,5.0))/(5.0*(5.0-s));
  
  v=sqrt(2.0*Ek*f/(g1*M));	
  
  rho=M/(4.0*M_PI*pow(R0,3)*f);
  
  td=3.0*kappa*rho*R0*R0/(pow(M_PI,2.0)*c);
  th=R0/v;
    
  Eth=4.0*M_PI*pow(R0,3)*7.57e-15*pow(T0,4.0)*Ith;
  t0=pow(2.0*th*td,0.5);	  
	   
   printf("#Supernova light curve model based on Arnett & Fu ApJ 340, 396 (1989)\n");
   printf("#\n");
   printf("#Inital model parameters\n");
   printf("#R0 = %lg cm\n",R0);
   printf("#Mej = %lg M_sol\n",M/Msol);
   printf("#MNi = %lg M_sol\n",Mni/Msol);
   printf("#Eth = %lg erg\n",Eth);
   printf("#Ekin = %lg erg\n",Ek);
   printf("#Trec = %lg K\n",Tion);
   printf("#kappa = %lg cm^2/g\n",kappa);
   printf("#a = %lg \n",a);
   printf("#s = %lg \n",s);   
   printf("#Ep = %lg erg\n",Ep);
   printf("#tp= %lg d\n",tp/day);
   printf("#Ag = %lg d^2\n",Ag/(day*day)); 
   printf("#\n");
   printf("#Calaculated physical properties \n");
   printf("#v = %lg km/s\n", v/1e5);
   printf("#rho0 = %lg g/cm^3\n", rho);
   printf("#td = %lg d\n",td/day);
   printf("#th = %lg d\n",th/day);
   printf("#t0 = %lg d\n",t0/day);
   printf("#\n");
   printf("#t/day    Lsn    logLsn     Ldiff      Lrec 	  xi \n");
       
  p1=Eni*Mni*tni/Eth;
  p2=tni/td;
  p3=tni/Eth;
  p4=tni/tco;
  p5=Eco/Eni; 
 
	while(t<=tmax)
	{ 
	  opac=1-exp(-Ag/(t*t));

  	  if(tp==0) Em=0;
  	  else Em=Ep/(tp*pow((1.0+t/tp),2.0));
  	  if(tp==0) Em1=0;
  	  else Em1=Ep/(tp*pow((1.0+(t+dt*0.5)/tp),2.0));
  	  
	  R=R0+v*t;	   
	  sigma=R/R0;			
	  sigma1=(R0+v*(t+dt*0.5))/R0;
	  Tc=T0*pow(F,0.25)/sigma;
	 
	 if(xi>=xmin && Tc>Tion) xi=temp(Tc,xh);
	 else xi=xmin;
	 
	  dxi=(series(xh)-series(xi))*xi;
	  	  	  	  
  	  g=Xni+p5*Xco;
	  dXni=-Xni*dt/tni;
	  dXco=(Xni-p4*Xco)*dt/tni;
	  ri=xi*R;
	  dr=dxi*R; if(dr>0.0) dr=0.0;
	    
	  if (s==0.0)Lion=-4.0*M_PI*rho*eta(xi,a)*pow(sigma,-3.0)*Q*ri*ri*dr/dt;
	  else Lion=-4.0*M_PI*rho*theta(xi,s)*pow(sigma,-3.0)*Q*ri*ri*dr/dt; 
	  if (Lion<0) Lion=0.0;	 
	  L=xi*F*Eth*opac/td; 
	   	  	     	  
	  Lsn=L+Lion;
	  logL=log10(Lsn);  
	  
	  if(j==0){
	  printf("%lf %lg %lg %lg %lg %lg \n",t/day, Lsn, logL, L, Lion, xi);
	  }
	  j++;
	  if(j==10000) j=0;
	  
	  /*Runge-Kutta method*/	     
	  dF1=sigma/pow(xi,3.0)*(p1*g-p2*F*xi-2.0*F*pow(xi,2.0)*dxi*tni/(sigma*dt)+p3*Em);
	  dF2=sigma1/pow(xi,3.0)*(p1*g-p2*(F+dF1*dt*0.5/tni)*xi-2.0*(F+dF1*dt*0.5/tni)*pow(xi,2)*dxi*tni/(sigma1*dt)+p3*Em1);
	  dF3=sigma1/pow(xi,3.0)*(p1*g-p2*(F+dF2*dt*0.5/tni)*xi-2.0*(F+dF2*dt*0.5/tni)*pow(xi,2)*dxi*tni/(sigma1*dt)+p3*Em1);
	  dF4=sigma1/pow(xi,3.0)*(p1*g-p2*(F+dF3*dt/tni)*xi-2.0*(F+dF3*dt/tni)*pow(xi,2)*dxi*tni/(sigma1*dt)+p3*Em1); 
	  F=F+(dF1+2.0*dF2+2.0*dF3+dF4)*dt/(6.0*tni);
	  
	  Xni=Xni+dXni;
	  Xco=Xco+dXco;
	  xh=xi;
	  t=t+dt;
	  
	} 	
}


