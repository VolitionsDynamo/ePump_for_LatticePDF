#include <math.h>
#include <iostream>
using namespace std;
#include "random.h"
#include <stdlib.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long Random::idum = -1;

double Random::uniform(void)
/* 
   This is the random number generator >>ran1<< of the
   Numerical Recipes in C book. 
*/
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (idum <= 0 || !iy) {
     if (-(idum) < 1) idum=1;
     else idum = -(idum);
     for (j=NTAB+7;j>=0;j--) {
        k=(idum)/IQ;
        idum=IA*(idum-k*IQ)-IR*k;
        if (idum < 0) idum += IM;
        if (j < NTAB) iv[j] = idum;
     }
     iy = iv[0];
  }
  k=(idum)/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if (idum < 0) idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp; 
}


double Random::normal(const double &z0, const double &S)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(z) =  exp(-(z-z0)^2/(2*S^2))/sqrt(2*PI*S^2)
        distribution in (-inf,inf).
        Uses method from "Numerical Recipes". 
------------------------------------------------------------------*/
{
  static int stored=0;
  static double zExtra;
  double v1, v2, rsq; 

  if (!stored) {
      do {
         v1=2.0*uniform()-1.0;
         v2=2.0*uniform()-1.0;
         rsq=v1*v1+v2*v2;
      } while (rsq>=1.0||rsq==0.0);
      double fac = sqrt(-2.0*log(rsq)/rsq);
      zExtra=v1*fac;
      stored=1;
      return S*v2*fac+z0;
  } 
  stored=0;
  return S*zExtra+z0;
}


double gammaln(double xx) 
/*------------------------------------------------------------------
        The log of the gamma function.
        Uses method from "Numerical Recipes". 
------------------------------------------------------------------*/
{
double x,y,tmp,ser;
static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
int j;

y=x=xx;
tmp=x+5.5;
tmp -= (x+0.5)*log(tmp);
ser=1.000000000190015;
for (j=0;j<=5;j++) ser += cof[j]/++y;
return -tmp+log(2.5066282746310005*ser/x);
}


int Random::poisson(const double &mu)
/*------------------------------------------------------------------
        Produce a random integer r (>=0) with a Poisson distribution
        with mean value mu. 
             P(r) = mu^r * exp(-mu) / r!  .
        Uses method from "Numerical Recipes".    
------------------------------------------------------------------*/
{
  static double sq,alxm,g;
  static double oldm=-1.0;  
  double em,t,y;

  if (mu < 12.0) {
      if (mu != oldm) {
          oldm=mu;
          g=exp(-mu);
      }
      em=-1.0;
      t=1.0;
      do {
         ++em;
         t *= uniform();
      } while (t>g);
   } else {
      if (mu != oldm) {
        oldm=mu;
        sq=sqrt(2.0*mu);
        alxm=log(mu);
        g=mu*alxm-gammaln(mu+1.0);
      }
      do {
         do {
            y=tan(PI*uniform());
            em= sq*y+mu;
         } while (em < 0.0);
         em=floor(em);
         t=0.9*(1.0+y*y)*exp(em*alxm-gammaln(em+1.0)-g);
       } while (uniform()>t);
   }
   return (int) em;
}


#define MAXIT 100
#define ACC 1.0e-8

double Random::phiAB(const double &A, const double &B) 
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(Phi) =  (1 + A*cos(Phi) + B*sin(Phi)) /(2*PI)
        with A^2 + B^2 < 1.
        Uses Newton-Raphson as in "Numerical Recipes". 
------------------------------------------------------------------*/
{
  double phil, phih, dphi, dphiold;
  double phi, phi0;
  double sphi, cphi;
  double f, df, temp;
  int j;

  phil = 0.0;
  phih = dphi = dphiold = 2.0*PI;
  phi = phi0 = dphi*uniform();  
  sphi = sin(phi);
  cphi = cos(phi);
  f = A*sphi + B*(1.0-cphi);
  df = 1.0 + A*cphi + B*sphi;
  for (j=1;j<=MAXIT;j++) {
      if ((((phi-phih)*df-f)*((phi-phil)*df-f) >=0.0)
         || (fabs(2.0*f) > fabs(dphiold*df))) {
         dphiold=dphi;
         dphi= 0.5*(phih-phil);
         phi=phil+dphi;
         if (phil == phi) return phi;
      } else {
         dphiold=dphi;
         dphi=f/df;
         temp=phi;
         phi -= dphi;
         if (temp == phi) return phi;
      }
      if (fabs(dphi) < ACC) return phi;
      sphi = sin(phi);
      cphi = cos(phi);
      f = phi - phi0 + A*sphi + B*(1.0-cphi);
      df = 1.0 + A*cphi + B*sphi;
      if (f < 0.0)
         phil=phi;
      else
         phih=phi;
  }          
  cerr <<"\nNewton-Raphson didn't converge in phiAB!\n";
  exit(1);
  return 0.0;
}


double Random::oneMinusZOverZ(const double &zl, const double &zh)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(z) =  ((1-z)/z) / (log(zh/zl)-zh+zl) 
        distribution in (zl,zh).
        Uses Newton-Raphson as in "Numerical Recipes". 
------------------------------------------------------------------*/
{
  double zll, zhh, dz, dzold;
  double z;
  double y, r;
  double f, df, temp;
  int j;

  zll = zl;
  zhh = zh;
  dz = dzold = zh-zl;
  y = uniform();
  z = zh*pow(zl/zh,y);
  r = zh-log(zh) + y*(log(zh/zl)-zh+zl);
  f = log(z)-z+r;
  df = (1.0 - z)/z;
  for (j=1;j<=MAXIT;j++) {
      if ((((z-zhh)*df-f)*((z-zll)*df-f) >=0.0)
         || (fabs(2.0*f) > fabs(dzold*df))) {
         dzold=dz;
         dz= 0.5*(zhh-zll);
         z=zll+dz;
         if (zll == z) return z;
      } else {
         dzold=dz;
         dz=f/df;
         temp=z;
         z -= dz;
         if (temp == z) return z;
      }
      if (fabs(dz) < ACC) return z;
      f = log(z)-z+r;
      df = (1.0 - z)/z;
      if (f < 0.0)
         zll=z;
      else
         zhh=z;
  }          
  cerr <<"\nNewton-Raphson didn't converge in oneMinusZoverZ!\n";
  exit(1);
  return 0.0;
}


double Random::betaReg(const double &b)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(z) =  (1-z^2)/(1-(b*z)^2)/norm 
        distribution in (-1,1).
        Uses Newton-Raphson as in "Numerical Recipes". 
------------------------------------------------------------------*/
{
  double zl, zh, dz, dzold;
  double z;
  double y, bsq, r, g, h;
  double f, df, temp;
  int j;

  zl = -1.0;
  zh = 1.0;
  dz = dzold = 2.0;
  y = uniform();
  bsq = b*b;
  g = (1.0+bsq)/(2.0*b);
  h = 1.0-bsq;
  r = (0.5-y)*(g*log((1.0+b)/(1.0-b))-1.0);
  if (fabs(r/bsq)<1.0) z=-r/bsq;
  else  z=0.0;
  f = r+0.5*(g*log((1.0+b*z)/(1.0-b*z))-z*h/(1.0-bsq*z*z));
  df = bsq*(1.0-z*z)/(1.0-bsq*z*z)/(1.0-bsq*z*z);
  for (j=1;j<=MAXIT;j++) {
      if ((((z-zh)*df-f)*((z-zl)*df-f) >=0.0)
         || (fabs(2.0*f) > fabs(dzold*df))) {
         dzold=dz;
         dz= 0.5*(zh-zl);
         z=zl+dz;
         if (zl == z) return z;
      } else {
         dzold=dz;
         dz=f/df;
         temp=z;
         z -= dz;
         if (temp == z) return z;
      }
      if (fabs(dz) < ACC) return z;
      f = r+0.5*(g*log((1.0+b*z)/(1.0-b*z))-z*h/(1.0-bsq*z*z));
      df = bsq*(1.0-z*z)/(1.0-bsq*z*z)/(1.0-bsq*z*z);
      if (f < 0.0)
         zl=z;
      else
         zh=z;
  }          
  cerr <<"\nNewton-Raphson didn't converge in betaReg!\n";
  exit(1);
  return 0.0;
}


double Random::electronDist(const double &b)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(z) =  b/2*(1-z)^(b/2-1)*(1+3/8*b) - b*(1+z)/4. 
        distribution in (0,1).
        (The Fadin-Kuraev form of the electron PDF.)      
        Uses Newton-Raphson as in "Numerical Recipes". 
------------------------------------------------------------------*/
{
  double zl, zh, dz, dzold;
  double z;
  double y, r;
  double f, df, temp;
  int j;

  zl = 0.0;
  zh = dz = dzold = 1.0;
  y = uniform();
  r = 1.0 + 3.0*b/8.0;
  z = 1.0 - pow((1.0-y),2.0/b);
  f = r*y-b/4.0*(z+z*z/2.0) -y;
  df = r*b/2.0*pow((1.0-z),(b/2.0-1.0)) - b/4.0*(1.0+z);
  for (j=1;j<=MAXIT;j++) {
      if ((((z-zh)*df-f)*((z-zl)*df-f) >=0.0)
         || (fabs(2.0*f) > fabs(dzold*df))) {
         dzold=dz;
         dz= 0.5*(zh-zl);
         z=zl+dz;
         if (zl == z) return z;
      } else {
         dzold=dz;
         dz=f/df;
         temp=z;
         z -= dz;
         if (temp == z) return z;
      }
      if (fabs(dz) < ACC) return z;
      f = r*(1.0-pow((1.0-z),b/2.0))-b/4.0*(z+z*z/2.0) -y;
      df = r*b/2.0*pow((1.0-z),(b/2.0-1.0)) - b/4.0*(1.0+z);
      if (f < 0.0)
         zl=z;
      else
         zh=z;
  }          
  cerr <<"\nNewton-Raphson didn't converge in electronDist!\n";
  exit(1);
  return 0.0;
}

double Random::Pqtog(const double &zl, const double &zh)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
            P(z) =  (1+(1-z)^2)/z / [2*(zl-zh-log(zl/zh))+0.5*(zh^2-zl^2)]
        distribution in (zl,zh).
        (the quark to gluon splitting function)
        Uses Newton-Raphson as in "Numerical Recipes".
------------------------------------------------------------------*/
{
  double zll, zhh, dz, dzold;
  double z;
  double y, r;
  double f, df, temp;
  int j;

  zll = zl;
  zhh = zh;
  dz = dzold = zh-zl;
  y = uniform();
  z = zh*pow(zl/zh,y);
  r = -2.0*log(zh)+2.0*zh-0.5*zh*zh -y*(2.0*(log(zl/zh)-zl+zh)
	+0.5*(zl*zl-zh*zh));
  f = 2.0*(log(z)-z)+0.5*z*z+r;
  df = 2.0/z-2.0+z;
  for (j=1;j<=MAXIT;j++) {
      if ((((z-zhh)*df-f)*((z-zll)*df-f) >=0.0)
         || (fabs(2.0*f) > fabs(dzold*df))) {
         dzold=dz;
         dz= 0.5*(zhh-zll);
         z=zll+dz;
         if (zll == z) return z;
      } else {
         dzold=dz;
         dz=f/df;
         temp=z;
         z -= dz;
         if (temp == z) return z;
      }
      if (fabs(dz) < ACC) return z;
      f = 2.0*(log(z)-z)+0.5*z*z+r;
      df = 2.0/z-2.0+z;
      if (f < 0.0)
         zll=z;
      else
         zhh=z;
  }          
  cerr <<"\nNewton-Raphson didn't converge in Pqtog!\n";
  exit(1);
  return 0.0;
}


double Random::Pgtog(double zl, double zh)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
            P(z) =  ((1+(1-z)^3)/z + (1+z^3)/(1-z) )/ norm.
        distribution in (zl,zh).
        (the gluon to gluon splitting function)
------------------------------------------------------------------*/
{
    double zlm=1.0-zl;
    double zhm=1.0-zh;
    weight w(2);
    double norm0=2.0*log(zh/zl)-3.0*(zh-zl)+1.5*(zh*zh-zl*zl)
	-(zh*zh*zh-zl*zl*zl)/3.0;
    double norm1=2.0*log(zlm/zhm)-3.0*(zlm-zhm)+1.5*(zlm*zlm-zhm*zhm)
	-(zlm*zlm*zlm-zhm*zhm*zhm)/3.0;
    w[0]=norm0;
    w[1]=norm1;
    w.normalize();
    int wp=w.choose();
    switch(wp) {
        case 0:
              return ggfun(zl,zh,norm0);
              break;
        case 1:
              return 1.0-ggfun(zhm,zlm,norm1);
              break;
    }
    cerr <<"\nUnknown error in Pgtog!\n";
    exit(1);
    return 0.0;
}


double Random::ggfun(const double &zl, const double &zh, const double norm)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
            P(z) =  (1+(1-z)^3)/z / norm
        distribution in (zl,zh).
        (Used by the gluon to gluon splitting function)
        Uses Newton-Raphson as in "Numerical Recipes".
------------------------------------------------------------------*/
{
  double zll, zhh, dz, dzold;
  double z;
  double y, r;
  double f, df, temp;
  int j;

  zll = zl;
  zhh = zh;
  dz = dzold = zh-zl;
  y = uniform();
  z = zh*pow(zl/zh,y);
  r = -2.0*log(zh)+3.0*zh-1.5*zh*zh +zh*zh*zh/3.0+ y*norm;
  f = 2.0*log(z)-3.0*z+1.5*z*z-z*z*z/3.0+r;
  df = 2.0/z-3.0+3.0*z-z*z;
  for (j=1;j<=MAXIT;j++) {
      if ((((z-zhh)*df-f)*((z-zll)*df-f) >=0.0)
         || (fabs(2.0*f) > fabs(dzold*df))) {
         dzold=dz;
         dz= 0.5*(zhh-zll);
         z=zll+dz;
         if (zll == z) return z;
      } else {
         dzold=dz;
         dz=f/df;
         temp=z;
         z -= dz;
         if (temp == z) return z;
      }
      if (fabs(dz) < ACC) return z;
          f = 2.0*log(z)-3.0*z+1.5*z*z-z*z*z/3.0+r;
          df = 2.0/z-3.0+3.0*z-z*z;
      if (f < 0.0)
         zll=z;
      else
         zhh=z;
  }
  cerr <<"\nNewton-Raphson didn't converge in ggfun!\n";
  exit(1);
  return 0.0;
}


double Random::SudaTfixed(const double &AA, const double &zl, const double &zh)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
            P(z) =  (z^A-z)/z/norm
        distribution in (zl,zh).
        (The fixed alpha_s Sudakov factor for top decay.)
        The norm<0 for A>1 and norm>0 for A<1.
        The limit A-->1 is well-defined, but not implemented here.
        Uses Newton-Raphson as in "Numerical Recipes". 
------------------------------------------------------------------*/
{
  double zll, zhh, dz, dzold;
  double z;
  double y, r, A;
  double f, df, temp;
  int j;

  zll = zl;
  zhh = zh;
  A=AA;
  dz = dzold = zh-zl;
  y = uniform();
  z =1.0-sqrt(y*(1.0-zll)*(1.0-zll)+(1.0-y)*(1.0-zhh)*(1.0-zhh));
  r =pow(zhh,A)-A*zhh-y*(pow(zhh,A)-A*zhh-pow(zll,A)+A*zll);
  f = (pow(z,A)-A*z-r)/(1.0-A);
  df = A*(pow(z,A)-z)/z/(1.0-A);
  for (j=1;j<=MAXIT;j++) {
      if ((((z-zhh)*df-f)*((z-zll)*df-f) >=0.0)
         || (fabs(2.0*f) > fabs(dzold*df))) {
         dzold=dz;
         dz= 0.5*(zhh-zll);
         z=zll+dz;
         if (zll == z) return z;
      } else {
         dzold=dz;
         dz=f/df;
         temp=z;
         z -= dz;
         if (temp == z) return z;
      }
      if (fabs(dz) < ACC) return z;
      f = (pow(z,A)-A*z-r)/(1.0-A);
      df = A*(pow(z,A)-z)/z/(1.0-A);
      if (f < 0.0)
         zll=z;
      else
         zhh=z;
  }          
  cerr <<"\nNewton-Raphson didn't converge in SudaTfixed!\n";
  exit(1);
  return 0.0;
}

double Int1(double a, double b,const  double &d, const double  &b0, const
              double &asb);
double trapzd(double a, double b, int n,const  double &d, const double  &b0, 
              const double &asb);
double dInt(const double &y,const  double &d, const double  &b0, const
              double &asb);

const double EPSI=1.0e-6;
const int JMAX=20;

double Int1(double a, double b, const  double &d, const double  &b0, const
              double &asb)
{
  double s[JMAX], sold[JMAX], h[JMAX];
  int j, i;

  h[0]=1.0;
  s[0]=trapzd(a,b,0,d,b0,asb);
  for (j=1;j<JMAX;j++) {
     h[j]=0.25*h[j-1];
     sold[0] = s[0];
     s[0] = trapzd(a,b,j,d,b0,asb);
     for (i=1;i<=j;i++) {
         sold[i] = s[i];
         s[i] = (h[j-i]*s[i-1]-h[j]*sold[i-1])/(h[j-i]-h[j]);
     }
     if (fabs(s[j]-s[j-1])<EPSI*fabs(s[j])) return s[j];
  }
  cerr <<"Too many steps in routine qromb.\n";
  exit(1);

  return 0.0;
}

   double trapzd(double a, double b, int n, const  double &d, 
       const double  &b0, const double &asb) 
{
  double x,tnm,sum,del;
  static double s;
  int it,j;

  if (n==0) {
      return (s=0.5*(b-a)*(dInt(a,d,b0,asb)+dInt(b,d,b0,asb)));
  } else {
    for (it=1,j=1;j<n;j++) it <<=1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for(sum=0.0,j=1;j<=it;j++,x+=del) sum += dInt(x,d,b0,asb);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}



double dInt(const double &y, const  double &d, const double  &b0, const
              double &asb) {
             return pow(y,-d-1.0+16.0/3.0/b0/asb)*
                     (d+16.0/3.0/b0/asb*(y-1.0))*
                      exp((y-1.0)/asb-16.0/3.0/b0/asb*(y-1.0));}


double Random::SudaTop(const double &as, const double &r0, 
                       const double &zl)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
        for the emission of a gluon during T-->BW decay in the interval
        in (zl,1.0).
        (Uses running alpha_s.)
        The inputs are:
               as = alpha_s(Q_max^2)
               r0 = Q_min^2/m_t^2/(1-m_W^2/m_t^2)^2
               zl = lower limit (usually=r0*(1+m_W/m_t)^2)
               zh = 1.0 (fixed)
        Uses Newton-Raphson as in "Numerical Recipes". 
------------------------------------------------------------------*/
{
  double zll, zhh, dz, dzold;
  double z;
  double y, r, flavors;
  double f, df, temp;
  double d,b0,asb;
  int j;
  double F, aa, al;

  flavors=5.0;
  b0 = (33.0-2.0*flavors)/3.0;
  d = 16.0/3.0/b0*(-log(r0)-1.0+0.25);
  asb=as*b0/4.0/PI;
  zll = zl;
  zhh = 1.0;
  dz = dzold = zhh-zll;
  y = uniform();
  z = pow( (1.0-y)*pow(zll,asb*d)+y*pow(zhh,asb*d),1.0/asb/d);
  al = 1.0+asb*log(zl);
  r = (1.0-y)*(
      pow(al,d)*exp(16.0/3.0/b0*(log(zl)-log(al)/asb))
      +pow(al,d-16.0/3.0/b0/asb)*pow(zl,16.0/3.0/b0)
             *Int1(al,1.0,d,b0,asb) ) +y;
  aa=1.0+asb*log(z);
  F=pow(aa,d)*exp(16.0/3.0/b0*(log(z)-log(aa)/asb))
         +pow(aa,d-16.0/3.0/b0/asb)*pow(z,16.0/3.0/b0)
          *Int1(aa,1.0,d,b0,asb);
  f = F-r;
  df = asb/z*(d+16.0/3.0/b0*log(z))*(F-z)/aa;
  for (j=1;j<=MAXIT;j++) {
      if ((((z-zhh)*df-f)*((z-zll)*df-f) >=0.0)
         || (fabs(2.0*f) > fabs(dzold*df))) {
         dzold=dz;
         dz= 0.5*(zhh-zll);
         z=zll+dz;
         if (zll == z) return z;
      } else {
         dzold=dz;
         dz=f/df;
         temp=z;
         z -= dz;
         if (temp == z) return z;
      }
      if (fabs(dz) < ACC) return z;
      aa=1.0+asb*log(z);
      F=pow(aa,d)*exp(16.0/3.0/b0*(log(z)-log(aa)/asb))
            +pow(aa,d-16.0/3.0/b0/asb)*pow(z,16.0/3.0/b0)
             *Int1(aa,1.0,d,b0,asb);
      f = F-r;
      df = asb/z*(d+16.0/3.0/b0*log(z))*(F-z)/aa;
      if (f < 0.0)
         zll=z;
      else
         zhh=z;
  }          
  cerr <<"\nNewton-Raphson didn't converge in SudaTop!\n";
  exit(1);
  return 0.0;
}


double Random::SudaQuark(const double &as, const double &r0, 
                       const double &zl, const double &zh)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
        for the emission of a gluon off a quark in the interval
        in (zl,zh).
        (Uses running alpha_s.)
        The inputs are:
               as = alpha_s(Q_max^2)
               r0 = Q_min^2/Q_max^2
               zl = lower limit (usually=r0)
               zh = upper limit (usually=1.0)
        Uses Newton-Raphson as in "Numerical Recipes". 
------------------------------------------------------------------*/
{
  double zll, zhh, dz, dzold;
  double z;
  double y, r, d, b0, asb, flavors;
  double f, df, temp;
  int j;
  double F, aa, al, ah;

  flavors=5.0;
  b0 = (33.0-2.0*flavors)/3.0;
  d = 16.0/3.0/b0*(-log(r0)-1.0+0.25);
  asb=as*b0/4.0/PI;
  zll = zl;
  zhh = zh;
  dz = dzold = zh-zl;
  y = uniform();
  z = pow( (1.0-y)*pow(zl,asb*d)+y*pow(zh,asb*d),1.0/asb/d);
  al = 1.0+asb*log(zl);
  ah = 1.0+asb*log(zh);
  r = (1.0-y)*pow(al,d)*exp(16.0/3.0/b0*(log(zl)-log(al)/asb))
     +y*pow(ah,d)*exp(16.0/3.0/b0*(log(zh)-log(ah)/asb));
  aa=1.0+asb*log(z);
  F=pow(aa,d)*exp(16.0/3.0/b0*(log(z)-log(aa)/asb));
  f = F-r;
  df = asb/z*(d+16.0/3.0/b0*log(z))*F/aa;
  for (j=1;j<=MAXIT;j++) {
      if ((((z-zhh)*df-f)*((z-zll)*df-f) >=0.0)
         || (fabs(2.0*f) > fabs(dzold*df))) {
         dzold=dz;
         dz= 0.5*(zhh-zll);
         z=zll+dz;
         if (zll == z) return z;
      } else {
         dzold=dz;
         dz=f/df;
         temp=z;
         z -= dz;
         if (temp == z) return z;
      }
      if (fabs(dz) < ACC) return z;
      aa=1.0+asb*log(z);
      F=pow(aa,d)*exp(16.0/3.0/b0*(log(z)-log(aa)/asb));
      f = F-r;
      df = asb/z*(d+16.0/3.0/b0*log(z))*F/aa;
      if (f < 0.0)
         zll=z;
      else
         zhh=z;
  }          
  cerr <<"\nNewton-Raphson didn't converge in SudaQuark!\n";
  exit(1);
  return 0.0;
}


weight::weight(const int &size) 
{
   nwts=size;
   w = new double[nwts];
   n = new double[nwts];
   for (int i=0; i<nwts; i++) {
      w[i] = 0.0;
      n[i] = 0.0;
   }
}


weight::weight(const weight &wgt) 
{
   nwts=wgt.nwts;
   w=new double[nwts];
   n=new double[nwts];   
   for (int i=0;i<nwts;i++) {
      w[i]=wgt.w[i];
      n[i]=wgt.n[i];
   }
}


weight& weight::operator=(const weight &wgt) 
{
   if (nwts!=wgt.nwts) {
     delete [] w;
     delete [] n;
     nwts=wgt.nwts;
     w=new double[nwts];
     n=new double[nwts];
   }
   for (int i=0;i<nwts;i++) {
      w[i]=wgt.w[i];
      n[i]=wgt.n[i];
   }
   return *this;
}


double weight::normalize(void)
{
   int i;
   double norm = 0.0;
   for (i=0;i<nwts;i++)
      norm += w[i];
   w[0] /=norm;
   n[0] = w[0];
   for (i=1;i<nwts;i++) {
      w[i] /= norm;
      n[i] = n[i-1] + w[i];
   }
   return norm;
}

int weight::choose(void)
{ 
  Random r;
  double R = r.uniform();
  int idum = 0;
  int maxi = nwts-1;

  while (idum<maxi && n[idum]<R)
     idum++;
  return idum;
}
