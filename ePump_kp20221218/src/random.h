#ifndef RANDOMHDR
#define RANDOMHDR

#include <math.h>
#ifdef __GNUG__
#define PI M_PI
#else
#define PI _PI
#endif

double gammaln(double x);  // necessary for Poisson distribution.


class Random {
/*
   Produce real random numbers with various probability 
   distributions.
*/
private:
   static long idum;
public:
   void seed(const long iseed) {idum = iseed;}  // Initialize the seed
                                                // with a NEGATIVE integer.
                                                // Default is -1.

   double uniform(void);         //  Uniform distribution in (0,1).

   double normal(const double &z0,const double &S);
          //  normal distribution with mean z0 and deviation S.
          //  exp(-(z-z0)^2/(2*S^2))/sqrt(2*PI*S^2) 

   inline double plus(void);    //  (1+c)/2 distribution in (-1,1).

   inline double minus(void);   //  (1-c)/2 distribution in (-1,1).

   inline double plusSq(void);  //  (1+c)^2 * 3/8 distribution in (-1,1).

   inline double minusSq(void); //  (1-c)^2 * 3/8 distribution in (-1,1).

   inline double sinSq(void);   //  (1-c^2) * 3/4 distribution in (-1,1).

   double phiAB(const double &A, const double &B);   
                    //  (1 + A*cos(phi) + B*sin(phi)) / (2*PI)
                    //         distribution in (0,2*PI).

   inline double reciprocalZ(const double &zl, const double &zh=1.0);
                    //  (1/z) / log(zh/zl) distribution in (zl,zh).

   inline double linearZ(const double &z0);
                    //  2z/(1-z0^2) distribution in (z0,1).

   inline double powerZ(const double &n,const double &zl, 
                                 const double &zh=1.0);
                    // z^n / ((zh^(n+1)-zl^(n+1))/(n+1)) 
                    // distribution in (zl,zh).

   inline double lorentzian(const double &a, const double &g);
                    // g/PI / ( (x-a)^2 + g^2 ) distribution in (-inf,inf).

   inline double massReg(const double &b);
                    //  b / (1-b*z) / log((1+b)/(1-b)) distribution
                    //                                  in (-1,1).

   inline double eDecay(const double &z0, const double &t);
                    //  exp( (z0-z)/t ) /t distribution from (z0,inf).


   double oneMinusZOverZ(const double &zl, const double &zh=1.0);
                    //  ((1-z)/z) / (log(zh/zl)-zh+zl)) distribution 
                    //                                  in (zl,zh).

   inline double oneOverZZplusA(const double &zl, const double &a=1.0);
                    //  1/z/(z+a) * a / (log(1+a/zl)) distribution 
                    //                                  in (zl, inf).

   inline double oneOverZplusAZplusB(const double &a, const double &b);
                    // 1/(z+a)/(z+b) * (b-a)/log(b/a) distribution
                    //                                  in (0, inf).

   double betaReg(const double &b);
                    //  (1-z^2)/(1-(b*z)^2)^2/norm
                    //              distribution in (-1,1).

   double electronDist(const double &b);
                    //  returns the Fadin-Kuraev form of the electron parton 
                    //  distribution in (0,1).  
                    //  That is,
                    //  b/2*(1-z)^(b/2-1)*(1+3/8*b) - b*(1+z)/4. 

   double Pqtog(const double &zl, const double &zh=1.0);
                    // returns the quark to gluon splitting function
                    // on the interval (zl,zh).  That is,
                    // (1+(1-z)^2)/z / [2*(zl-zh-log(zl/zh))+0.5*(zh^2-zl^2)].
                    
   double Pgtog(double zl, double zh);
                    // returns the gluon to gluon splitting function
                    // on the interval (zl,zh).  That is,
                    // ((1+(1-z)^3)/z + (1+z^3)/(1-z) )/ norm.

   double Pgtoq(const double &zl, const double &zh=1.0);
                    // returns the gluon to quark splitting function
                    // on the interval (zl,zh).  That is,
                    // (z^2+(1-z)^2)/ norm.
                    
   double ggfun(const double &zl, const double &zh, const double norm);
		    // returns the function (1+(1-z)^3) / z / norm
		    // on the interval (zl,zh). Used by Pgtog.
                                        
   double SudaTfixed(const double &AA, const double &zl, const double &zh=1.0);
                    // returns the Sudakov factor for top decay (fixed alpha_s)
                    // on the interval (zl,zh).  That is,
                    // (z^A-z)/z/norm, 
                    // where norm>0 for A<1 and norm<0 for A>1.
                    // Although the distribution is well-defined for A-->1,
                    // that case does not work here.

   double SudaTop(const double &as, const double &r0,
                       const double &zl);
                    // returns the Sudakov factor for gluon emission
                    // during T-->BW decay (running alpha_s)
                    // on the interval (zl,1.0).  
                    // The inputs are :
                    //      as = alpha_s(Q_max^2)
                    //      r0 = Q_min^2/m_t^2/(1-m_W^2/m_t^2)^2     
                    // Usually zl=r0*(1+m_W/m_t)^2 
                    // The upper limit is fixed at zh=1.

   double SudaQuark(const double &as, const double &r0,
                       const double &zl, const double &zh=1.0);
                    // returns the Sudakov factor for gluon emission
                    // off a quark (running alpha_s)
                    // on the interval (zl,zh).  
                    // The inputs are :
                    //      as = alpha_s(Q_max^2)
                    //      r0 = Q_min^2/Q_max^2     
                    // Usually zl=r0 and zh=1.

    int poisson(const double &mu);
                    // Returns an random integer from a poisson
                    // distribution with mean value mu.
                    //     P(r) = mu^r * exp(-mu) / r! .
                    // with  r>=0.

};


inline double Random::plus(void) 
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(CosChi) = (1 + CosChi)/2
------------------------------------------------------------------*/
{
  return 2.0*sqrt(uniform()) - 1.0;
}
 

inline double Random::minus(void) 
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(CosChi) = (1 - CosChi)/2
------------------------------------------------------------------*/
{
  return 1.0 - 2.0*sqrt(uniform());
}


inline double Random::plusSq(void) 
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(CosChi) = (1 + CosChi)^2 * 3/8
------------------------------------------------------------------*/
{
  return 2.0*pow(uniform(),1.0/3.0) - 1.0;
}


inline double Random::minusSq(void) 
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(CosChi) = (1 - CosChi)^2 * 3/8
V------------------------------------------------------------------*/
{
  return 1.0 - 2.0*pow(uniform(),1.0/3.0);
}


inline double Random::sinSq(void) 
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(CosChi) =  SinChi^2 * 3/4
------------------------------------------------------------------*/
{
  double delta = acos(1.0 - 2.0*uniform());
  return  -2.0*cos((delta + PI)/3.0);
}

inline double Random::reciprocalZ(const double &zl, const double &zh) 
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(z) = (1/z) / log(zh/zl)
        on the interval (zl,zh).
------------------------------------------------------------------*/
{
  return zh*pow(zl/zh,uniform());
}

inline double Random::linearZ(const double &z0) 
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(z) = 2z/(1-z0^2) 
        on the interval (z0,1).
------------------------------------------------------------------*/
{
  return sqrt(z0*z0+uniform()*(1.0-z0*z0));
}

inline double Random::powerZ(const double &n, const double &zl,
               const double &zh) 
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(z) = z^n / ((zh^(n+1)-zl^(n+1))/(n+1)) 
         on the interval (zl,zh).
------------------------------------------------------------------*/
{
  double y=uniform();
  return pow( (1.0-y)*pow(zl,n+1.0) + y*pow(zh,n+1.0) , 1.0/(n+1.0) );
}

inline double Random::lorentzian(const double &a, const double &g)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(x) = g/PI / ( (x-a)^2 + g^2)
        on the interval (-inf,inf).
------------------------------------------------------------------*/
{
  return a + g/tan(PI*uniform());
}

inline double Random::massReg(const double &b)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(x) = b/(1-b*x)/log((1+b)/(1-b))
        on the interval (-1,1).
------------------------------------------------------------------*/
{
  double y=uniform();
  return (1.0-pow(1.0+b,1.0-y)*pow(1.0-b,y))/b;
}


inline double Random::eDecay(const double &z0, const double &t)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(x) = exp( (z0-z)/t ) / t
        on the interval (z0, inf).
------------------------------------------------------------------*/
{
  double y=uniform();
  return z0-t*log(y);
}


inline double Random::oneOverZZplusA(const double &zl, const double &a)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(z) = 1/z/(z+a) * a / log(1+a/zl)
        on the interval (zl, inf).
------------------------------------------------------------------*/
{
  double y=uniform();
  return a/(pow(1.0+a/zl,y)-1.0);
}


inline double Random::oneOverZplusAZplusB(const double &a, const double &b)
/*------------------------------------------------------------------
        Produce a random variable with the probability distribution
                P(z) = 1/(z+a)/(z+b) * (b-a)/log(b/a) 
        on the interval (0, inf).
------------------------------------------------------------------*/
{
  double y=uniform();
  double x=pow(b/a,y);  
  return a*(x-1.0)/(1-x*a/b);
}


class weight {
 private:
   int nwts;     //  number of possible outcomes

   double *w;    //  array of weights

   double *n;    //  noodle array

 public:
   weight(void) {nwts=0; w=0; n=0;}
   weight(const int &size) ;
   weight( const weight &wgt) ;
   ~weight(void) {delete [] w; delete [] n;}
   weight& operator=(const weight &wgt) ;
   double& operator[](int index) {return w[index];}
   int size(void) {return nwts;}
   double normalize(void) ;  // normalizes the weights, prepares the
                             // noodle, and returns the normalization 
                             // factor

   int choose(void) ;        // chooses one of the possible outcomes
                             // with the probability given by their
                             // weights
};


#endif
