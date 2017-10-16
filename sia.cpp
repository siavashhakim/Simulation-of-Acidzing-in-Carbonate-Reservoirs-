#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;
// parameter
int main(){
const int t=100;
const int n=20;
const int m=9;
const int mm=0.3;
const int omega=0.4;
const int omega0=0.4;
const int beta=1;
const int alph=0.25;
const int a0=50;
const int C0=0.02;
const int rhos=2.71;
const int laength=10.4;
const int K0=7.896e-12;
const int viscoa=0.13;
const int densitya=1.7*7;
const int landax=0.5;
const int landaT=0.1;
const int landax=0.5;
const int Dm=3.6e-5;
const int r0=0.001;
const int Sh_infinit=3;
const int eps00=0.2;
const int u00=1.04;
const int alpha=1.37;
const int alpha0=0.5;
double viscokinetica;
double density;
double Sc;
double Sh0;
int e;
int mz;
int nz;
int tt;
int i;
int j;
int k;
int kk;
int iz;
double Rep0;
double dx;
double dyy;
double dt;
double dy;
double ka0;
double Da0;
double Nac;
double Pe;
double Rep;
double Sh;
// Condition of some parameters (constant)
viscokinetica= viscoa/density;
Rep0= 2*u00*r0/viscokinetica;
Sc=viscoa/(densitya*Dm);
Sh0=Sh_infinit+(0.7/(mm^0.5))*((Rep0)^0.5)*(Sc^(1/3));
ka0=(Sh0*Dm)/(2*r0);
Da0=ka0*a0*laength/u00;
Nac=alpha*C0/rhos;
Pe=laength*u00/Dm;
Rep=2*u00*r0*density/viscoa;
Sh=Sh_infinit+(0.7/(sqrt(mm)))*(sqrt(Rep0))*(Sc^(1/3));
// Matrix
float f1[n][m][t];
float f11[n][m][t];
float f111[n][m][t];
float f1111[n][m][t];
float f2[n][m][t];
float A[n*m][n*m][t];
float C[n*m][1][t];
float Cf[n][m][t];
float Kx[n][m][t];
float r[n][m][t];
float a[n][m][t];
float Da[n][m][t];
float eps[n][m][t];
float epss[n][m][t];
float P[n][m][t];
float PP[n][m][t];
float ue[n][m][t];
float uw[n][m][t];
float un[n][m][t];
float us[n][m][t];
float u[n][m][t];
float Dx[n][m][t];
float DT[n][m][t];
float ka[n][m][t];
float K1[n][m][t];
float K2[n][m][t];
float K3[n][m][t];
float K4[n][m][t];
float KK1[n][m][t];
float KK2[n][m][t];
float KK3[n][m][t];
float KK4[n][m][t];
// Some equations
dx=1/(n-1);
dyy=1/(m-1);
dt=1/(t-1);
dy=dyy
nz=n;
mz=m;

for (tt=0; tt<=t; tt++){
// Initial condition
for (j=1;j<=mz;j++){

for (e=1;e<=nz;e++){
eps[e][j][0]=0.4;
Cf[e][j][0]=0;
Kx[e][j][0]=1;
r[e][j][0]=1;
a[e][j][0]=1;
Da[e][j][0]=Da0;
}
}


// BC1 Injection
for (j=1;j<=mz;j++){
Cf[1][j][tt]=0;
}
Cf[1][(m-1)/2+1][tt]=0;
for (i=1; i<=n; i++){
for (j=1; j<=m; j++){
if (eps[i][j][tt]>1)
eps[i][j][tt]=1;
if (eps[i][j][tt]<0)
eps[i][j][tt]=0;
}
}
for (j=1; j<=m; j++){
if (epss[n][j][tt]>0.4)

cout<<"Breakthrough"<<tt<<endl; 

}
// pressure
for (j=2; j<=m-1; j++){
for (e=2; e<=n-1; e++){
k=(e-1)*m+j;
A[k][k][tt]=((2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e+1][j][tt]))))+(1/(alph^2))*(2*dx/(dyy*((1/Kx[e][j][tt])+(1/Kx[e][j+1][tt]))))+(1/(alph^2))*(2*dx/(dyy*((1/Kx[e][j][tt])+(1/Kx[e][j-1][tt]))))+(2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e-1][j][tt])))));
           
 A[k][k+m][tt]=-(2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e+1][j][tt]))));
        
A[k][k-m][tt]=-(2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e-1][j][tt]))));
        
A[k][k+1][tt]=-(1/(alph^2))*(2*dx/(dyy*((1/Kx[e][j][tt])+(1/Kx[e][j+1][tt]))));
        
A[k][k-1][tt]=-(1/(alph^2))*(2*dx/(dyy*((1/Kx[e][j][tt]+1/Kx[e][j-1][tt]))));
        
C[k][1][tt]=-1*Da[e][j][tt]*Nac*a[e][j][tt]*Cf[e][j][tt]*dx*dyy;
}
}
//e=1 , left Boundary Condition
for (j=2;j<=m-1; j++){
if (j==(m-1)/2+1) {
    
    e=1;
    k=(e-1)*m+j;
    A[k][k][tt]=-1;
    A[k][k+m][tt]=1;
    A[k][k-1][tt]=0;
    A[k][k+1][tt]=0;
	C[k][1][tt]=-(dx/Kx[e][j][tt]);
}
else{
e=1;
k=(e-1)*m+j;
A[k][k][tt]=1;
C[k][1][tt]=0;
}
// e=n, Right Boundary Condition
for (j=2; j<=m-1; j++){
    e=n;
    k=(e-1)*m+j;
    A[k][k][tt]=1;
    A[k][k-1][tt]=0;
    A[k][k-m][tt]=0;
    A[k][k+1][tt]=0;
    C[k][1][tt]=0;
}
//j=1 , Floor Boundary Condition
for (e=2; e<=n-1; e++){
    j=1;
    k=(e-1)*m+j;
    A[k][k][tt]=((2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e+1][j][tt]))))+(1/(alph^2))*(2*dx/(dyy*((1/Kx[e][j][tt])+(1/Kx[e][j+1][tt]))))+(2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e-1][j][tt])))));
    
A[k][k-m][tt]=-(2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e-1][j][tt]))));
    
A[k][k+m][tt]=-(2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e+1][j][tt]))));
    
A[k][k+1][tt]=-(1/(alph^2))*(2*dx/(dyy*((1/Kx[e][j][tt])+(1/Kx[e][j+1][tt]))));
   
C[k][1][tt]=-1*Da[e][j][tt]*Nac*a[e][j][tt]*Cf[e][j][tt]*dx*dyy;
}
// j=m , Ceil Boundary Condition
 
for (e=2; e<=n-1; e++){
    j=m;
    k=(e-1)*m+j;
    A[k][k][tt]=((2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e+1][j][tt]))))+(1/(alph^2))*(2*dx/(dyy*((1/Kx[e][j][tt])+(1/Kx[e][j-1][tt]))))+(2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e-1][j][tt])))));
    
A[k][k-1][tt]=-(1/(alph^2))*(2*dx/(dyy*((1/Kx[e][j][tt])+(1/Kx[e][j-1][tt]))));
    
A[k][k+m][tt]=-(2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e+1][j][tt]))));
   
A[k][k-m][tt]=-(2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e][j][tt]))));
    
C[k][1][tt]=-1*Da[e][j][tt]*Nac*a[e][j][tt]*Cf[e][j][tt]*dx*dyy;
}

// Left-floor BC
e=1; j=1;
k=(e-1)*m+j;
A[k][k][tt]=((2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e+1][j][tt]))))+(1/(alph^2))*(2*dx/(dyy*((1/Kx[e][j][tt])+(1/Kx[e][j+1][tt])))));
A[k][k+1][tt]=-(1/(alph^2))*(2*dx/(dyy*((1/Kx[e][j][tt])+(1/Kx[e][j+1][tt]))));
A[k][k+m][tt]=-(2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e+1][j][tt]))));
C[k][1][tt]=-1*Da[e][j][tt]*Nac*a[e][j][tt]*Cf[e][j][tt]*dx*dyy;

// Left-Ceil BC
e=1; j=m;
k=(e-1)*m+j;
A[k][k][tt]=((1/(alph^2))*(2*dx/(dyy*((1/Kx[e][j][tt])+(1/Kx[e][j-1][tt]))))+(2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e+1][j][tt])))));
A[k][k-1][tt]=-(1/(alph^2))*(2*dx/(dyy*((1/Kx[e][j][tt])+(1/Kx[e][j-1][tt]))));
A[k][k+m][tt]=-(2*dyy/(dx*((1/Kx[e][j][tt])+(1/Kx[e+1][j][tt]))));
C[k][1][tt]=-1*Da[e][j][tt]*Nac*a[e][j][tt]*Cf[e][j][tt]*dx*dyy;

// Right- Floor BC
e=n; j=1;
k=(e-1)*m+j;
A[k][k][tt]=1;
C[k][1][tt]=0;
 
// Right-Ceil BC
e=n; j=m;
k=(e-1)*m+j;
A[k][k][tt]=1;
C[k][1][tt]=0;

//////////////////////////////////////////////////////////////
//???? INV function
////////////////////////////////////////////////////////////
kk=1;
for (e=1;e<=n;e++){
    for (j=1;j<=m;j++){
        P[e][j][tt]=PP[kk][tt];
        kk=kk+1;
  }
}

P[1][1][tt]=(P[1][2][tt]+P[2][1][tt])/2;
P[1][m][tt]=(P[1][m-1][tt]+P[2][m][tt])/2;

// Velocity
for (j=1;j<=mz;j++){
   for (e=1;e<=nz-1;e++){
ue[e][j][tt]=(-2*(P[e+1][j][tt]-P[e][j][tt])/(dx*((1/Kx[e][j][tt])+(1/Kx[e+1][j][tt]))));
 
   }
}
for (j=1;j<=mz;j++){

ue[n][j][tt]=0;
}
 
for (j=1;j<=mz;j++){
   for (e=2;e<=nz-1;e++){
 
uw[e][j][tt]=(-2*(P[e][j][tt]-P[e-1][j][tt])/(dx*((1/Kx[e][j][tt])+(1/Kx[e-1][j][tt]))));
 
   }
}

for (j=1;j<=mz;j++){

uw[1][j][tt]=0;
}
 
for (j=1;j<=mz-1;j++){
   for (e=1;e<=nz;e++){
 
un[e][j][tt]=(-2*(P[e][j+1][tt]-P[e][j][tt])/(dx*((1/Kx[e][j][tt])+(1/Kx[e][j+1][tt]))));
 
  }
}
un[1][mz][tt]=0;
 
for (j=2;j<=mz;j++){

   for (e=1;e<=nz;e++){
 
us[e][j][tt]=(-2*(P[e][j][tt]-P[e][j-1][tt])/(dx*((1/Kx[e][j][tt])+(1/Kx[e][j-1][tt]))));
}
}
for (e=1;e<=nz;e++){
us[e][1][tt]=0;
}
for (j=1;j<=mz;j++){
    for (e=1;e<=nz;e++){
   u[e][j][tt]=sqrt((((ue[e][j][tt]+uw[e][j][tt])/2)^2)+(((un[e][j][tt]+us[e][j][tt])/2)^2)); 
  
 Dx[e][j][tt]=(alpha0*eps[e][j][tt]*Dm+2*landax*u00*u[e][j][tt]*r[e][j][tt]*r0)/Dm;
   
 DT[e][j][tt]=(alpha0*eps[e][j][tt]*Dm+2*landaT*u00*u[e][j][tt]*r[e][j][tt]*r0)/Dm;
   }
}
 
// Equations
for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

    ka[iz][j][tt]=(Sh*Dm)/(2*r0*r[iz][j][tt]);  
    
  
  }
}

for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

if (iz==1 && j!=1 && j!=mz)
     f1[iz][j][tt]=0;

else if (iz==n && j!=1 && j!=mz)
      
f1[iz][j][tt]=((-max(ue[iz][j][tt],0)*Cf[iz][j][tt]+max(uw[iz][j][tt],0)*Cf[iz-1][j][tt]-max(-uw[iz][j][tt],0)*Cf[iz][j][tt])/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*Cf[iz][j][tt]+max(-un[iz][j][tt],0)*Cf[iz][j+1][tt]+max(us[iz][j][tt],0)*Cf[iz][j-1][tt]-max(-us[iz][j][tt],0)*Cf[iz][j][tt])/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((-(Cf[iz][j][tt]-Cf[iz-1][j][tt])/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+((Cf[iz][j+1][tt]-Cf[iz][j][tt])/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt])))-(Cf[iz][j][tt]-Cf[iz][j-1][tt])/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*Cf[iz][j][tt]*(1+Nac*Cf[iz][j][tt]))/eps[iz][j][tt]; 

else if (j==1 && iz!=1 && iz!=n)
f1[iz][j][tt]=((-max(ue[iz][j][tt],0)*Cf[iz][j][tt]+max(-ue[iz][j][tt],0)*Cf[iz+1][j][tt]+max(uw[iz][j][tt],0)*Cf[iz-1][j][tt]-max(-uw[iz][j][tt],0)*Cf[iz][j][tt])/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*Cf[iz][j][tt]+max(-un[iz][j][tt],0)*Cf[iz][j+1][tt]-max(-us[iz][j][tt],0)*Cf[iz][j][tt])/(alph*dy*eps[iz][j][tt]))+(2/Pe)*(((Cf[iz+1][j][tt]-Cf[iz][j][tt])/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt])))-(Cf[iz][j][tt]-Cf[iz-1][j][tt])/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+((Cf[iz][j+1][tt]-Cf[iz][j][tt])/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*Cf[iz][j][tt]*(1+Nac*Cf[iz][j][tt]))/eps[iz][j][tt]; 

else if (j==mz && iz!=1 && iz!=n)
f1[iz][j][tt]=((-max(ue[iz][j][tt],0)*Cf[iz][j][tt]+max(-ue[iz][j][tt],0)*Cf[iz+1][j][tt]+max(uw[iz][j][tt],0)*Cf[iz-1][j][tt]-max(-uw[iz][j][tt],0)*Cf[iz][j][tt])/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*Cf[iz][j][tt]+max(us[iz][j][tt],0)*Cf[iz][j-1][tt]-max(-us[iz][j][tt],0)*Cf[iz][j][tt])/(alph*dy*eps[iz][j][tt]))+(2/Pe)*(((Cf[iz+1][j][tt]-Cf[iz][j][tt])/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt])))-(Cf[iz][j][tt]-Cf[iz-1][j][tt])/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(-(Cf[iz][j][tt]-Cf[iz][j-1][tt])/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*Cf[iz][j][tt]*(1+Nac*Cf[iz][j][tt]))/eps[iz][j][tt];

else if (iz==n && j==1)
f1[iz][j][tt]=((-max(ue[iz][j][tt],0)*Cf[iz][j][tt]+max(uw[iz][j][tt],0)*Cf[iz-1][j][tt]-max(-uw[iz][j][tt],0)*Cf[iz][j][tt])/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*Cf[iz][j][tt]+max(-un[iz][j][tt],0)*Cf[iz][j+1][tt]-max(-us[iz][j][tt],0)*Cf[iz][j][tt])/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((-(Cf[iz][j][tt]-Cf[iz-1][j][tt])/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+((Cf[iz][j+1][tt]-Cf[iz][j][tt])/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*Cf[iz][j][tt]*(1+Nac*Cf[iz][j][tt]))/eps[iz][j][tt];

else if (iz==n && j==mz)
f1[iz][j][tt]=((-max(ue[iz][j][tt],0)*Cf[iz][j][tt]+max(uw[iz][j][tt],0)*Cf[iz-1][j][tt]-max(-uw[iz][j][tt],0)*Cf[iz][j][tt])/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*Cf[iz][j][tt]+max(us[iz][j][tt],0)*Cf[iz][j-1][tt]-max(-us[iz][j][tt],0)*Cf[iz][j][tt])/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((-(Cf[iz][j][tt]-Cf[iz-1][j][tt])/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(-(Cf[iz][j][tt]-Cf[iz][j-1][tt])/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*Cf[iz][j][tt]*(1+Nac*Cf[iz][j][tt]))/eps[iz][j][tt];

else if (iz==1 && j==1)
f1[iz][j][tt]=((-max(ue[iz][j][tt],0)*Cf[iz][j][tt]+max(-ue[iz][j][tt],0)*Cf[iz+1][j][tt]-max(-uw[iz][j][tt],0)*Cf[iz][j][tt])/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*Cf[iz][j][tt]+max(-un[iz][j][tt],0)*Cf[iz][j+1][tt]-max(-us[iz][j][tt],0)*Cf[iz][j][tt])/(alph*dy*eps[iz][j][tt]))+(2/Pe)*(((Cf[iz+1][j][tt]-Cf[iz][j][tt])/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt]))))/dx+((Cf[iz][j+1][tt]-Cf[iz][j][tt])/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*Cf[iz][j][tt]*(1+Nac*Cf[iz][j][tt]))/eps[iz][j][tt]; 

else if (iz==1 && j==mz)
f1[iz][j][tt]=((-max(ue[iz][j][tt],0)*Cf[iz][j][tt]+max(-ue[iz][j][tt],0)*Cf[iz+1][j][tt]-max(-uw[iz][j][tt],0)*Cf[iz][j][tt])/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*Cf[iz][j][tt]+max(us[iz][j][tt],0)*Cf[iz][j-1][tt]-max(-us[iz][j][tt],0)*Cf[iz][j][tt])/(alph*dy*eps[iz][j][tt]))+(2/Pe)*(((Cf[iz+1][j][tt]-Cf[iz][j][tt])/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt]))))/dx+(-(Cf[iz][j][tt]-Cf[iz][j-1][tt])/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*Cf[iz][j][tt]*(1+Nac*Cf[iz][j][tt]))/eps[iz][j][tt];

else
f1[iz][j][tt]=((-max(ue[iz][j][tt],0)*Cf[iz][j][tt]+max(-ue[iz][j][tt],0)*Cf[iz+1][j][tt]+max(uw[iz][j][tt],0)*Cf[iz-1][j][tt]-max(-uw[iz][j][tt],0)*Cf[iz][j][tt])/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*Cf[iz][j][tt]+max(-un[iz][j][tt],0)*Cf[iz][j+1][tt]+max(us[iz][j][tt],0)*Cf[iz][j-1][tt]-max(-us[iz][j][tt],0)*Cf[iz][j][tt])/(alph*dy*eps[iz][j][tt]))+(2/Pe)*(((Cf[iz+1][j][tt]-Cf[iz][j][tt])/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt])))-(Cf[iz][j][tt]-Cf[iz-1][j][tt])/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+((Cf[iz][j+1][tt]-Cf[iz][j][tt])/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt])))-(Cf[iz][j][tt]-Cf[iz][j-1][tt])/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*Cf[iz][j][tt]*(1+Nac*Cf[iz][j][tt]))/eps[iz][j][tt]; 


}
}

// RK4
// K1
for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

K1[iz][j][tt]=dt*f1[iz][j][tt];
}
}

//K2
for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

if (iz==1 && j!=1 && j!=mz)
     f11[iz][j][tt]=0;

else if (iz==n && j!=1 && j!=mz)
      
f11[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+0.5*K1[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+0.5*K1[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ 0.5*K1[iz][j+1][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+0.5*K1[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((-((Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])-(Cf[iz-1][j][tt]+0.5*K1[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+0.5*K1[iz][j+1][tt])-(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt])))-((Cf[iz][j][tt]+0.5*K1[iz][j][tt])-(Cf[iz][j-1][tt]+0.5*K1[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+0.5*K1[iz][j][tt])))/eps[iz][j][tt]; 

else if (j==1 && iz!=1 && iz!=n)
f11[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ 0.5*K1[iz+1][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+ 0.5*K1[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ 0.5*K1[iz][j+1][tt]))-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ 0.5*K1[iz+1][j][tt])-(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt])))-((Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])-(Cf[iz-1][j][tt]+ 0.5*K1[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+ 0.5*K1[iz][j+1][tt])-(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])))/eps[iz][j][tt]; 

else if (j==mz && iz!=1 && iz!=n)
f11[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ 0.5*K1[iz+1][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+ 0.5*K1[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt] 0.5*K1[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+ 0.5*K1[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ 0.5*K1[iz+1][j][tt])-(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt])))-((Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])-(Cf[iz-1][j][tt]+ 0.5*K1[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(-((Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])-(Cf[iz][j-1][tt]+ 0.5*K1[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])))/eps[iz][j][tt];

else if (iz==n && j==1)
f11[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+ 0.5*K1[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ 0.5*K1[iz][j+1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((-((Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])-(Cf[iz-1][j][tt]+ 0.5*K1[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+ 0.5*K1[iz][j+1][tt])-(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])))/eps[iz][j][tt];

else if (iz==n && j==mz)
f11[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+ 0.5*K1[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+ 0.5*K1[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((-((Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])-(Cf[iz-1][j][tt]+ 0.5*K1[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(-((Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])-(Cf[iz][j-1][tt]+ 0.5*K1[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])))/eps[iz][j][tt];

else if (iz==1 && j==1)
f11[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ 0.5*K1[iz+1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ 0.5*K1[iz][j+1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ 0.5*K1[iz+1][j][tt])-(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+ 0.5*K1[iz][j+1][tt])-(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])))/eps[iz][j][tt]; 

else if (iz==1 && j==mz)
f11[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ 0.5*K1[iz+1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+ 0.5*K1[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ 0.5*K1[iz+1][j][tt])-(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt]))))/dx+(-((Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])-(Cf[iz][j-1][tt]+ 0.5*K1[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])))/eps[iz][j][tt];

else
f11[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ 0.5*K1[iz+1][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+ 0.5*K1[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ 0.5*K1[iz][j+1][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+ 0.5*K1[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ 0.5*K1[iz+1][j][tt])-(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt])))-((Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])-(Cf[iz-1][j][tt]+ 0.5*K1[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+ 0.5*K1[iz][j+1][tt])-(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt])))-((Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])-(Cf[iz][j-1][tt]+ 0.5*K1[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K1[iz][j][tt])))/eps[iz][j][tt]; 


}
}

for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

K2[iz][j][tt]=dt*f11[iz][j][tt];
}
}

//K3
for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

if (iz==1 && j!=1 && j!=mz)
     f111[iz][j][tt]=0;

else if (iz==n && j!=1 && j!=mz)
      
f111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+0.5*K2[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+0.5*K2[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ 0.5*K2[iz][j+1][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+0.5*K2[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((-((Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])-(Cf[iz-1][j][tt]+0.5*K2[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+0.5*K2[iz][j+1][tt])-(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt])))-((Cf[iz][j][tt]+0.5*K2[iz][j][tt])-(Cf[iz][j-1][tt]+0.5*K2[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+0.5*K2[iz][j][tt])))/eps[iz][j][tt]; 

else if (j==1 && iz!=1 && iz!=n)
f111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ 0.5*K2[iz+1][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+ 0.5*K2[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ 0.5*K2[iz][j+1][tt]))-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ 0.5*K2[iz+1][j][tt])-(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt])))-((Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])-(Cf[iz-1][j][tt]+ 0.5*K2[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+ 0.5*K2[iz][j+1][tt])-(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])))/eps[iz][j][tt]; 

else if (j==mz && iz!=1 && iz!=n)
f111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ 0.5*K2[iz+1][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+ 0.5*K2[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt] 0.5*K2[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+ 0.5*K2[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ 0.5*K2[iz+1][j][tt])-(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt])))-((Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])-(Cf[iz-1][j][tt]+ 0.5*K2[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(-((Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])-(Cf[iz][j-1][tt]+ 0.5*K2[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])))/eps[iz][j][tt];

else if (iz==n && j==1)
f111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+ 0.5*K2[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ 0.5*K2[iz][j+1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((-((Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])-(Cf[iz-1][j][tt]+ 0.5*K2[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+ 0.5*K2[iz][j+1][tt])-(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])))/eps[iz][j][tt];

else if (iz==n && j==mz)
f111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+ 0.5*K2[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+ 0.5*K2[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((-((Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])-(Cf[iz-1][j][tt]+ 0.5*K2[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(-((Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])-(Cf[iz][j-1][tt]+ 0.5*K2[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])))/eps[iz][j][tt];

else if (iz==1 && j==1)
f111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ 0.5*K2[iz+1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ 0.5*K2[iz][j+1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ 0.5*K2[iz+1][j][tt])-(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+ 0.5*K2[iz][j+1][tt])-(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])))/eps[iz][j][tt]; 

else if (iz==1 && j==mz)
f111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ 0.5*K2[iz+1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+ 0.5*K2[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ 0.5*K2[iz+1][j][tt])-(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt]))))/dx+(-((Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])-(Cf[iz][j-1][tt]+ 0.5*K2[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])))/eps[iz][j][tt];

else
f111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ 0.5*K2[iz+1][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+ 0.5*K2[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ 0.5*K2[iz][j+1][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+ 0.5*K2[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ 0.5*K2[iz+1][j][tt])-(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt])))-((Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])-(Cf[iz-1][j][tt]+ 0.5*K2[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+ 0.5*K2[iz][j+1][tt])-(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt])))-((Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])-(Cf[iz][j-1][tt]+ 0.5*K2[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])))/eps[iz][j][tt]; 


}
}

for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

K3[iz][j][tt]=dt*f111[iz][j][tt];
}
}

//K4
for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

if (iz==1 && j!=1 && j!=mz)
     f1111[iz][j][tt]=0;

else if (iz==n && j!=1 && j!=mz)
      
f1111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+ K3[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ K3[iz][j+1][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+ K3[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((-((Cf[iz][j][tt]+ K3[iz][j][tt])-(Cf[iz-1][j][tt]+0.5*K2[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+ K3[iz][j+1][tt])-(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt])))-((Cf[iz][j][tt]+ K3[iz][j][tt])-(Cf[iz][j-1][tt]+ K3[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ K3[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ K3[iz][j][tt])))/eps[iz][j][tt]; 

else if (j==1 && iz!=1 && iz!=n)
f1111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ K3[iz+1][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+ K3[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ K3[iz][j+1][tt]))-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ K3[iz+1][j][tt])-(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt])))-((Cf[iz][j][tt]+K3[iz][j][tt])-(Cf[iz-1][j][tt]+ K3[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+K3[iz][j+1][tt])-(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ K3[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+ K3[iz][j][tt])))/eps[iz][j][tt]; 

else if (j==mz && iz!=1 && iz!=n)
f1111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ K3[iz+1][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+K3[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt] K3[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+K3[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ K3[iz+1][j][tt])-(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt])))-((Cf[iz][j][tt]+ 0.5*K2[iz][j][tt])-(Cf[iz-1][j][tt]+K3[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(-((Cf[iz][j][tt]+K3[iz][j][tt])-(Cf[iz][j-1][tt]+K3[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ K3[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+K3[iz][j][tt])))/eps[iz][j][tt];

else if (iz==n && j==1)
f1111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+K3[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ K3[iz][j+1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((-((Cf[iz][j][tt]+ K3[iz][j][tt])-(Cf[iz-1][j][tt]+K3[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+K3[iz][j+1][tt])-(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ K3[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+K3[iz][j][tt])))/eps[iz][j][tt];

else if (iz==n && j==mz)
f1111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+K3[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+K3[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((-((Cf[iz][j][tt]+ K3[iz][j][tt])-(Cf[iz-1][j][tt]+K3[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(-((Cf[iz][j][tt]+K3[iz][j][tt])-(Cf[iz][j-1][tt]+K3[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ K3[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+K3[iz][j][tt])))/eps[iz][j][tt];

else if (iz==1 && j==1)
f1111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ K3[iz+1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ K3[iz][j+1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ K3[iz+1][j][tt])-(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+K3[iz][j+1][tt])-(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ K3[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+K3[iz][j][tt])))/eps[iz][j][tt]; 

else if (iz==1 && j==mz)
f1111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+K3[iz+1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+K3[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ K3[iz+1][j][tt])-(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt]))))/dx+(-((Cf[iz][j][tt]+K3[iz][j][tt])-(Cf[iz][j-1][tt]+K3[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ K3[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+K3[iz][j][tt])))/eps[iz][j][tt];

else
f1111[iz][j][tt]=((-max(ue[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(-ue[iz][j][tt],0)*(Cf[iz+1][j][tt]+ K3[iz+1][j][tt])+max(uw[iz][j][tt],0)*(Cf[iz-1][j][tt]+K3[iz-1][j][tt])-max(-uw[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dx*eps[iz][j][tt]))+((-max(un[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt])+max(-un[iz][j][tt],0)*(Cf[iz][j+1][tt]+ K3[iz][j+1][tt])+max(us[iz][j][tt],0)*(Cf[iz][j-1][tt]+K3[iz][j-1][tt])-max(-us[iz][j][tt],0)*(Cf[iz][j][tt]+ K3[iz][j][tt]))/(alph*dy*eps[iz][j][tt]))+(2/Pe)*((((Cf[iz+1][j][tt]+ K3[iz+1][j][tt])-(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz+1][j][tt])))-((Cf[iz][j][tt]+K3[iz][j][tt])-(Cf[iz-1][j][tt]+K3[iz-1][j][tt]))/(dx*((1/Dx[iz][j][tt])+(1/Dx[iz-1][j][tt]))))/dx+(((Cf[iz][j+1][tt]+K3[iz][j+1][tt])-(Cf[iz][j][tt]+ K3[iz][j][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j+1][tt])))-((Cf[iz][j][tt]+K3[iz][j][tt])-(Cf[iz][j-1][tt]+K3[iz][j-1][tt]))/(dy*((1/DT[iz][j][tt])+(1/DT[iz][j-1][tt]))))/((alph^2)*dy))/eps[iz][j][tt]-(Da[iz][j][tt]*a[iz][j][tt]*(Cf[iz][j][tt]+ K3[iz][j][tt])*(1+Nac*(Cf[iz][j][tt]+K3[iz][j][tt])))/eps[iz][j][tt]; 


}
}

for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

K4[iz][j][tt]=dt*f1111[iz][j][tt];
}
}
//
for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){
Cf[iz][j][tt+1]= Cf[iz][j][tt]+(K1[iz][j][tt]+2*K2[iz][j][tt]+2*K3[iz][j][tt]+K4[iz][j][tt])/6;
	}
}
//

// Porosity
for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

f2[iz][j][tt]=Da[iz][j][tt]*Nac*a[iz][j][tt]*Cf[iz][j][tt+1];
}
}
// KK1=KK2=KK3=KK4
for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

KK1[iz][j][tt]=dt* f2[iz][j][tt];
}
}

for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

eps[iz][j][tt+1]= eps[iz][j][tt]+ KK1[iz][j][tt];

}
}

for (j=1;j<=mz;j++){
    for (iz=1;iz<=nz;iz++){

// Permeability
Kx[iz][j][tt+1]=(((omega/omega0)^2)*(eps[iz][j][tt+1]/eps00(iz,j))*(((eps[iz][j][tt+1]*(1-eps00(iz,j)))/(eps00(iz,j)*(1-eps[iz][j][tt+1])))^(2*beta)))); 
// radius of pore   
r[iz][j][tt+1]=(((omega/omega0))*(((eps[iz][j][tt+1]*(1-eps00(iz,j)))/(eps00(iz,j)*(1-eps[iz][j][tt+1])))^(beta)))); 
   
// Interficial area
a[iz][j][tt+1]=(((omega/omega0))*(eps[iz][j][tt+1]/eps00(iz,j))*(((eps00(iz,j)*(1-eps[iz][j][tt+1]))/(eps[iz][j][tt+1]*(1-eps00(iz,j))))^beta))); 


ka[iz][j][tt+1]=(Sh*Dm)/(2*r0*r[iz][j][tt+1]); 
 
Da[iz][j][tt+1]=(ka[iz][j][tt+1]*a0*laength/u00); 
  
  }
}

}
}