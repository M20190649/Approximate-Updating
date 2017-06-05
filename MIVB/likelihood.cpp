// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp> // Making use of many of the distributions in Boost
#include "likelihood.h"
#include "evCopula.c"
#include "hfunc.h"
using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace boost::math;

#define UMAX  1-1e-10
#define UMIN  1e-10
#define XEPS 1e-4
#define XINFMAX DBL_MAX


///////////////////////////////////////////////////////
// New
double qgamma (double x, double shape, double rate, double a, double b) {
  boost::math::gamma_distribution<> aux (shape, 1.0/rate);
  return quantile(aux, x);
}

double pgamma (double x, double shape, double rate, double a, double b) {
  boost::math::gamma_distribution<> aux (shape, 1.0/rate);
  return cdf(aux, x);
}

double qnorm5(double x, double mean, double sd, double a, double b) {
  boost::math::normal_distribution<> aux (mean, sd);
  return quantile(aux, x);
}

double pnorm5(double x, double mean, double sd, double a, double b) {
  boost::math::normal_distribution<> aux (mean, sd);
  return cdf(aux, x);
}

double qt(double x, double df, double a, double b) {
  boost::math::students_t_distribution<> aux (df);
  return quantile(aux, x);
}

double dt(double x, double df, double a){
  boost::math::students_t_distribution<> aux (df);
  return pdf(aux, x);
}

double pt(double x, double df, double a){
  boost::math::students_t_distribution<> aux (df);
  return cdf(aux, x);
}

double pt(double x, double df, double a, double b){
  boost::math::students_t_distribution<> aux (df);
  return cdf(aux, x);
}



double StableGammaDivision(double x1, double x2)
{
  int i;
  double a1, a2, b1, b2, sum=1.0;
  a1 = fmod(max(x1,x2),1.0);
  a2 = max(x1,x2)-a1;
  b1 = fmod(min(x1,x2),1.0);
  b2 = min(x1,x2)-b1;
  if(a1==0.0 && b1==0.0)
  {
    for(i=1 ; i<(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=b2 ; i<(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
  }
  else if(a1>0.0 && b1==0.0)
  {
    for(i=1 ; i<(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=(int)b2 ; i<=(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum *= gamma(a1);
  }
  else if(a1==0.0 && b1>0.0)
  {
    for(i=1 ; i<=(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=((int)b2+1) ; i<(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum /= gamma(b1);
  }
  else if(a1>0.0 && b1>0.0)
  {
    for(i=1 ; i<=(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=((int)b2+1) ; i<=(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum *= gamma(a1)/gamma(b1);
  }
  if(x2 > x1) sum = 1.0/sum;
  return sum;
}


double log1mexp(double a)
{
  double result;
  if (a<log(2)) {
    result=log(-expm1(-a));
  }else{
    result=log1p(-exp(-a));
  }
  return result;
}

void archCDF(double* u, double* v, int* n, double* param, int* copula, double* out)
{
  int j;
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14;
  
  for(j=0;j<*n;j++)
  {
    if(u[j]>UMAX && v[j]>UMAX){ out[j]=1;}
    else if(u[j]>UMAX){ out[j]=v[j];}
    else if(v[j]>UMAX){ out[j]=u[j];}
    else if(u[j]<UMIN || v[j]<UMIN){ out[j]=0;}
    else
    {
      if(*copula==3)	//Clayton
      {
        t1 = pow(u[j],-param[0]);
        t2 = pow(v[j],-param[0]);
        t3 = t1+t2-1;
        out[j] = pow(t3,-1/param[0]);
      }
      else if(*copula==4)	//Gumbel
      {
        t1 = -log(u[j]);
        t2 = -log(v[j]);
        t3 = pow(t1,param[0]);
        t4 = pow(t2,param[0]);
        t5 = t3+t4;
        t6 = -pow(t5,1/param[0]);
        out[j] = exp(t6);
      }
      else if(*copula==5)	//Frank
      {
        if (param[0]>0) {
          t1=-log1p(exp(-param[0]) * expm1(param[0]-u[j]*param[0])/expm1(-param[0]));
          t2=-log1p(exp(-param[0]) * expm1(param[0]-v[j]*param[0])/expm1(-param[0]));
          out[j] = -log1mexp(t1+t2-log1mexp(param[0]))/param[0];
        } else {
          out[j] =-1/param[0] * log(1 + exp(-(-log((exp(-param[0] * u[j]) - 1)/(exp(-param[0]) - 1)) + -log((exp(-param[0] * v[j]) - 1)/(exp(-param[0]) - 1)))) * (exp(-param[0]) - 1));
        }
      }
      else if(*copula==6)	//Joe
      {
        t1 = 1-u[j];
        t2 = 1-v[j];
        t3 = pow(t1,param[0]);
        t4 = pow(t2,param[0]);
        t5 = t3*t4;
        out[j] = 1-pow(t3+t4-t5,1/param[0]);
      }
      else if(*copula==7)	//BB1
      {
        t1 = pow(u[j],-param[0]);
        t2 = pow(v[j],-param[0]);
        t3 = t1-1;
        t4 = t2-1;
        t5 = pow(t3,param[1]);
        t6 = pow(t4,param[1]);
        t7 = t5+t6;
        t8 = pow(t7,1/param[1]);
        out[j] = pow(1+t8,-1/param[0]);
      }
      else if(*copula==8)	//BB6
      {
        t1 = 1-u[j];
        t2 = 1-v[j];
        t3 = pow(t1,param[0]);
        t4 = pow(t2,param[0]);
        t5 = 1-t3;
        t6 = 1-t4;
        t7 = -log(t5);
        t8 = -log(t6);
        t9 = pow(t7,param[1]);
        t10 = pow(t8,param[1]);
        t11 = t9+t10;
        t12 = pow(t11,1/param[1]);
        t13 = exp(-t12);
        t14 = 1-t13;
        out[j] = 1-pow(t14,1/param[0]);
      }
      else if(*copula==9)	//BB7
      {
        t1 = 1-u[j];
        t2 = 1-v[j];
        t3 = pow(t1,param[0]);
        t4 = pow(t2,param[0]);
        t5 = 1-t3;
        t6 = 1-t4;
        t7 = pow(t5,-param[1]);
        t8 = pow(t6,-param[1]);
        t9 = t7+t8-1;
        t10 = pow(t9,-1/param[1]);
        t11 = 1-t10;
        t12 = pow(t11,1/param[0]);
        out[j] = 1-t12;
      }
      else if(*copula==10)    //BB8
      {
        double nu;
        t1 = param[1]*u[j];
        t2 = param[1]*v[j];
        t3 = 1-t1;
        t4 = 1-t2;
        t5 = pow(t3,param[0]);
        t6 = pow(t4,param[0]);
        t7 = 1-t5;
        t8 = 1-t6;
        nu = 1-param[1];
        nu = pow(nu,param[0]);
        nu = 1-nu;
        nu = 1/nu;
        t9 = 1-nu*t7*t8;
        t10 = pow(t9,1/param[0]);
        out[j] = 1/param[1]*(1-t10);
      }
      else if(*copula==41)
      {
        t1=qgamma(1.0-u[j],param[0],1,1,0);
        t2=qgamma(1.0-v[j],param[0],1,1,0);
        t3=pow(pow(t1,param[0])+pow(t2,param[0]),(1.0/param[0]));
        out[j]=1.0-pgamma(t3,param[0],1,1,0);
      }
    }
  }
  
}





void dbb1(double* u, double* v, int* n, double* param, double* out)
{
  int i;
  double th, de;
  double t1, t2, t3, t16, t17, t38, t39, t4, t5, t6, t7, t9, t10, t12, t13, t20, t24, t25, t27, t29, t32, t33, t34, t36, t43, t59;
  
  th = param[0];
  de = param[1];
  
  for(i=0;i<*n;i++)
  {
    t1 = pow(u[i],(-th));
    t2 = t1-1.0;
    t3 = pow(t2,de);
    t16 = 1./u[i];
    t17 = 1./t2;
    t38 = t1*t16;
    t39 = t38*t17;
    t4 = pow(v[i],(-th));
    t5 = t4-1.0;
    t6 = pow(t5,de);
    t7 = t3+t6;
    t9 = pow(t7,(1./de));
    t10 = 1.0+t9;
    t12 = pow(t10,(-1./th));
    t13 = t12*t9;
    t20 = 1./t10;
    t24 = t9*t9;
    t25 = t12*t24;
    t27 = 1./v[i];
    t29 = 1./t5;
    t32 = t7*t7;
    t33 = 1./t32;
    t34 = t10*t10;
    t36 = t33/t34;
    t43 = t4*th;
    t59 = t43*t27*t29;
    
    out[i] = t25*t6*t27*t4*t29*t36*t3*t39-t13*t6*t43*t27*t29*t33*t3*t38*t17*t20+
      t13*t3*t38*t17*t33*t20*t6*de*t59+t25*t3*t39*t36*t6*t59;
  }
  
}


void dbb6(double* u, double* v, int* n, double* param, double* out)
{
  int i;
  double th, de;
  double t1, t2, t3, t4, t5, t12, t16, t32, t38, t39, t40, t47, t50, t61, t90, t6, t7, t8, t9, t10, t11, t13, t14, t35, t36, t37, t42, t48, t53, t56, t57, t59, t78, t80, t87, t93;
  
  th = param[0];
  de = param[1];
  
  for(i=0;i<*n;i++)
  {
    t1 = 1.0-u[i];
    t2 = pow(t1,th);
    t3 = 1.0-t2;
    t4 = log(t3);
    t5 = pow(-t4,de);
    t12 = 1/de;
    t16 = 1/th;
    t32 = de-1.0;
    t38 = 2.0*de;
    t39 = -1.0+t38;
    t40 = pow(-t4,t39);
    t47 = 3.0*de-1.0;
    t50 = pow(-t4,t32);
    t61 = pow(-t4,t47);
    t90 = pow(-t4,t38);
    t6 = 1.0-v[i];
    t7 = pow(t6,th);
    t8 = 1.0-t7;
    t9 = log(t8);
    t10 = pow(-t9,de);
    t11 = t5+t10;
    t13 = pow(t11,t12);
    t14 = exp(-t13);
    t35 = pow(t11,-2.0*t32*t12);
    t36 = t35*th;
    t37 = exp(t13);
    t42 = pow(-t9,t39);
    t48 = pow(-t9,t47);
    t53 = t13*de;
    t56 = pow(-t9,t32);
    t57 = t37*t50*t56;
    t59 = t13*th;
    t78 = t37-1.0;
    t80 = pow(t78*t14,t16);
    t87 = t78*t78;
    t93 = pow(-t9,t38);
    
    out[i] = (2.0*t36*t37*t40*t42+t36*t37*t48*t50+t53*th*t57-t59*t57+
      t36*t37*t61*t56-2.0*t35*t40*t42-t35*t61*t56-t53*th*t50*t56+t59*t50*t56-
      t35*t48*t50) *t80*t7*t2/t3/t8/t87/(t90+2.0*t5*t10+t93)/t1/t6;
  }
  
}


void dbb7(double* u, double* v, int* n, double* param, double* out)
{
  int i;
  double th, de;
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t11, t12, t14, t15, t16, t18, t20, t23, t24, t25, t27, t30, t31, t32, t35, t37, t42, t54;
  
  th = param[0];
  de = param[1];
  
  for(i=0;i<*n;i++)
  {
    t1 = 1.0-u[i];
    t2 = pow(t1,th);
    t3 = 1.0-t2;
    t4 = pow(t3,-de);
    t5 = 1.0-v[i];
    t6 = pow(t5,th);
    t7 = 1.0-t6;
    t8 = pow(t7,-de);
    t9 = t4+t8-1.0;
    t11 = pow(t9,-1.0/de);
    t12 = 1.0-t11;
    t14 = pow(t12,1.0/th);
    t15 = t11*t11;
    t16 = t14*t15;
    t18 = 1./t5;
    t20 = 1./t7;
    t23 = t9*t9;
    t24 = 1./t23;
    t25 = t12*t12;
    t27 = t24/t25;
    t30 = t2/t1;
    t31 = 1./t3;
    t32 = t30*t31;
    t35 = t14*t11;
    t37 = t6*th;
    t42 = 1./t12;
    t54 = t37*t18*t20;
    
    out[i] = -t16*t8*t6*t18*t20*t27*t4*t32 + t35*t8*t37*t18*t20*t24*t4*t30*t31*t42+
      t35*t4*t30*t31*t24*t42*t8*de*t54+t16*t4*t32*t27*t8*t54;
  }
  
}


void dbb8(double* u, double* v, int* n, double* param, double* out)
{
  int i;
  double th, de;
  double t2, t3, t12, t16, t6, t7, t10, t11, t33, t38, t39, t49, t59, t69, t25, t26, t29, t44, t45, t50, t54, t62, t67;
  
  th = param[0];
  de = param[1];
  
  for(i=0;i<*n;i++)
  {
    t2 = 1.0-de*u[i];
    t3 = pow(t2,th);
    t10 = 1.0-de;
    t11 = pow(t10,th);
    t12 = 1.0-t11;
    t16 = 1/th;
    t33 = th*t3;
    t38 = 2.0*th;
    t39 = pow(t10,t38);
    t49 = pow(t2,t38);
    t59 = pow(t10,3.0*th);
    t69 = t12*t12;
    t6 = 1.0-de*v[i];
    t7 = pow(t6,th);
    t25 = t3*t7;
    t26 = t11-t7-t3+t25;
    t29 = pow(-t26/t12,t16);
    t44 = pow(t6,t38);
    t45 = t3*t44;
    t50 = t49*t7;
    t54 = t49*t44;
    t62 = -2.0*t25*t11+t25-t33*t7+3.0*t33*t7*t11-3.0*t33*t7*t39+t25*t39+
      2.0* t45*t11-t45*t39+2.0*t50*t11-t50*t39-2.0*t54*t11+t54*t39+t54-
      t50-t45+t33*t7*t59;
    t67 = t26*t26;
    out[i] = -de*t29*t62/t6/t2/t67/t69;
  }
  
}



void LL_mod2(double* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
  double* negv;
  double* negu;
  negv = (double *) malloc(*n*sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  double ntheta, nnu;
  double nfamily;
  ntheta = -*theta;
  nnu = -*nu;
  
  for(int i=0;i<*n;i++)
  {
    if(u[i]<UMIN) u[i]=UMIN;
    else if(u[i]>UMAX) u[i]=UMAX;
    if(v[i]<UMIN) v[i]=UMIN;
    else if(v[i]>UMAX) v[i]=UMAX;
  }
  if((*family)==43)
  {
    nfamily=3;
    if(*theta > 0){
      ntheta=2*(*theta)/(1-*theta);
      LL(&nfamily, n, u,  v, &ntheta, nu, loglik);
    }else{
      ntheta=-2*(*theta)/(1+*theta);
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      LL(&nfamily, n, negu,  v, &ntheta, &nnu, loglik);
    }
  }else if((*family)==44)
  {
    nfamily=4;
    if(*theta > 0){
      ntheta=1/(1-*theta);
      LL(&nfamily, n, u,  v, &ntheta, nu, loglik);
    }else{
      ntheta=1/(1+*theta);
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      LL(&nfamily, n, negu,  v, &ntheta, &nnu, loglik);
    }
  }else{
    
    if((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30) | (*family==61))	// 90? rotated copulas
    {
      nfamily = (*family)-20;
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      LL(&nfamily, n, negu,  v, &ntheta, &nnu, loglik);
    }
    else if((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40) | (*family==71))	// 270? rotated copulas
    {
      nfamily = (*family)-30;
      for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
      LL(&nfamily, n, u,  negv, &ntheta, &nnu, loglik);
    }
    else if((*family==124) | (*family==224))
    {
      nfamily = (*family)-20;
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      LL(&nfamily, n, v, negu, &ntheta, nu, loglik);
    }
    else if((*family==134) | (*family==234))
    {
      nfamily = (*family)-30;
      for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
      LL(&nfamily, n, negv, u, &ntheta, nu, loglik);
    }
    else {
      LL(family, n, u,  v, theta, nu, loglik);
    }
  }
  free(negv);
  free(negu);
}


//////////////////////////////////////////////////////////////
// Function to compute log-likelihood for bivariate copula
// Input:
// family    copula family (0=independent, 1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank)
// n         sample size
// u         first variable of data set
// v         second variable of data set
// theta     dependency parameter
// nu        degrees-of-freedom for students copula
// loglik    output
//////////////////////////////////////////////////////////////
void LL(double* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
  int j;
  double *dat, rho, ll=0.0, t1=0.0, t2=0.0, f;
  //Allocate memory:
  dat = Calloc(2,double);
  
  for(int i=0;i<*n;i++)
  {
    if(u[i]<UMIN) u[i]=UMIN;
    else if(u[i]>UMAX) u[i]=UMAX;
    if(v[i]<UMIN) v[i]=UMIN;
    else if(v[i]>UMAX) v[i]=UMAX;
  }
  
  //Compute log-likelihood:
  if(*family==0) //independent
    ll = 0;
  else if(*family==1) //Gaussian
  {
    rho=*theta;
    for(j=0;j<*n;j++)
    {
      dat[0]=u[j]; dat[1]=v[j];
      t1 = qnorm(dat[0],0.0,1.0,1,0); t2 = qnorm(dat[1],0.0,1.0,1,0);
      f = 1.0/sqrt(1.0-pow(rho,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*rho*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(rho,2.0))));
      if(log(f)>XINFMAX) ll += log(XINFMAX);
      else if(f < DBL_MIN) ll += log(DBL_MIN);
      else ll += log(f);
    }
  }
  else if(*family==2) //Student
  {
    rho=*theta;
    for(j=0;j<*n;j++)
    {
      dat[0] = u[j]; dat[1] = v[j];
      t1 = qt(dat[0],*nu,1,0); t2 = qt(dat[1],*nu,1,0);
      f = StableGammaDivision((*nu+2.0)/2.0,*nu/2.0)/(*nu*pi*sqrt(1.0-pow(rho,2.0))*dt(t1,*nu,0)*dt(t2,*nu,0))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho*t1*t2)/(*nu*(1.0-pow(rho,2.0))),-(*nu+2.0)/2.0);
      if(log(f)>XINFMAX) ll += log(XINFMAX);
      else if(f < DBL_MIN) ll += log(DBL_MIN);
      else ll += log(f);
    }
  }
  else if(*family==3) //Clayton
  {
    if(*theta == 0) ll = 0;
    else if(*theta < 1e-10) ll = 0;
    else
    {
      for(j=0;j<*n;j++)
      {
        dat[0] = u[j]; dat[1] = v[j];
        f=log1p(*theta)-(1.0+*theta)*log(dat[0]*dat[1])-(2.0+1.0/(*theta))*log(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0);
        if(f>XINFMAX) ll += log(XINFMAX);
        else if(f<log(DBL_MIN)) ll += log(DBL_MIN);
        else ll += f;
      }
    }
  }
  else if(*family==4) //Gumbel
  {
    for(j=0;j<*n;j++)
    {
      dat[0] = u[j]; dat[1] = v[j];
      t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
      f= -pow(t1,1.0/(*theta))+(2.0/(*theta)-2.0)*log(t1)+(*theta-1.0)*log(log(dat[0])*log(dat[1]))-log(dat[0]*dat[1])+log1p((*theta-1.0)*pow(t1,-1.0/(*theta)));
      
      if(f>XINFMAX) ll += log(XINFMAX);
      else if(f<log(DBL_MIN))ll += log(DBL_MIN);
      else ll += f;
    }
  }
  else if(*family==5) // Frank
  {
    for(j=0;j<*n;j++)
    {
      if (fabs(*theta) < 1e-10) {
        ll = 0;
      } else {
        dat[0] = u[j]; dat[1] = v[j];
        f = (*theta*(exp(*theta)-1.0)*exp(*theta*dat[1]+*theta*dat[0]+*theta))/pow(exp(*theta*dat[1]+*theta*dat[0])-exp(*theta*dat[1]+*theta)-exp(*theta*dat[0]+*theta)+exp(*theta),2.0);
        if(log(f)>XINFMAX) ll += log(XINFMAX);
        else if(f < DBL_MIN) ll += log(DBL_MIN);
        else ll += log(f);
      }
    }
  }
  else if(*family==6)	//Joe
  {
    for(j=0;j<*n;j++)
    {
      f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
      if(log(f)>XINFMAX) ll += log(XINFMAX);
      else if(f < DBL_MIN) ll += log(DBL_MIN);
      else ll += log(f);
    }
  }
  else if(*family==7)	//BB1
  {
    if(*theta == 0){
      for(j=0;j<*n;j++)
      {
        dat[0] = u[j]; dat[1] = v[j];
        t1 = pow(-log(dat[0]),*nu)+pow(-log(dat[1]),*nu);
        f= -pow(t1,1/(*nu))+(2/(*nu)-2)*log(t1)+(*nu-1)*log(log(dat[0])*log(dat[1]))-log(dat[0]*dat[1])+log(1+(*nu-1)*pow(t1,-1.0/(*nu)));
        if(f>XINFMAX) ll += log(XINFMAX);
        else if(f<log(DBL_MIN))ll += log(DBL_MIN);
        else ll += f;
      }
    }else{
      
      double *param, *fuc;
      param=Calloc(2,double);
      param[0]=*theta;
      param[1]=*nu;
      fuc = Calloc(*n,double);
      dbb1(u, v, n, param, fuc);
      for(j=0;j<*n;j++)
      {
        if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
        else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
        else ll += log(fuc[j]);
      }
      Free(fuc); Free(param);
    }
  }
  else if(*family==8)	//BB6
  {
    double *param, *fuc;
    param=Calloc(2,double);
    param[0]=*theta;
    param[1]=*nu;
    fuc = Calloc(*n,double);
    dbb6(u, v, n, param, fuc);
    for(j=0;j<*n;j++)
    {
      if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
      else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
      else ll += log(fuc[j]);
    }
    Free(fuc); Free(param);
  }
  else if(*family==9)	//BB7
  {
    if(*nu==0)
    {
      for(j=0;j<*n;j++)
      {
        f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
        if(log(f)>XINFMAX) ll += log(XINFMAX);
        else if(f < DBL_MIN) ll += log(DBL_MIN);
        else ll += log(f);
      }
    }
    else
    {
      double *param, *fuc;
      param=Calloc(2,double);
      param[0]=*theta;
      param[1]=*nu;
      fuc = Calloc(*n,double);
      dbb7(u, v, n, param, fuc);
      for(j=0;j<*n;j++)
      {
        if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
        else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
        else ll += log(fuc[j]);
      }
      Free(fuc); Free(param);
    }
  }
  else if(*family==10) //BB8
  {
    double *param, *fuc;
    param=Calloc(2,double);
    param[0]=*theta;
    param[1]=*nu;
    fuc = Calloc(*n,double);
    dbb8(u, v, n, param, fuc);
    for(j=0;j<*n;j++)
    {
      if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
      else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
      else ll += log(fuc[j]);
    }
    Free(fuc); Free(param);
    
  }
  else if(*family==13) //rotated Clayton (180?)
  {
    if(*theta == 0) ll = 0;
    else if(*theta < XEPS) ll = 0;
    else
    {
      for(j=0;j<*n;j++)
      {
        dat[0] = 1-u[j]; dat[1] = 1-v[j];
        f = (1.0+*theta)*pow(dat[0]*dat[1],-1.0-*theta)*pow(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0,-2.0-1.0/(*theta));
        f = max(f,0.0);
        if(log(f)>XINFMAX) ll += log(XINFMAX);
        else if(f < DBL_MIN) ll += log(DBL_MIN);
        else ll += log(f);
      }
    }
  }
  else if(*family==14) //rotated Gumbel (180?)
  {
    for(j=0;j<*n;j++)
    {
      dat[0] = 1-u[j]; dat[1] = 1-v[j];
      t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
      t2 = exp(-pow(t1,1.0/(*theta)));
      f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*theta))*pow(log(dat[0])*log(dat[1]),*theta-1.0)*(1.0+(*theta-1.0)*pow(t1,-1.0/(*theta)));
      if(log(f)>XINFMAX) ll += log(XINFMAX);
      else if(f < DBL_MIN) ll += log(DBL_MIN);
      else ll += log(f);
    }
  }
  else if(*family==16) //rotated Joe (180?)
  {
    for(j=0;j<*n;j++)
    {
      u[j]=1-u[j]; v[j]=1-v[j];
      f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
      if(log(f)>XINFMAX) ll += log(XINFMAX);
      else if(f < DBL_MIN) ll += log(DBL_MIN);
      else ll += log(f);
      u[j]=1-u[j]; v[j]=1-v[j];
    }
  }
  else if(*family==17) //rotated BB1
  {
    if(*theta == 0){
      for(j=0;j<*n;j++)
      {
        dat[0] = 1-u[j]; dat[1] = 1-v[j];
        t1 = pow(-log(dat[0]),*nu)+pow(-log(dat[1]),*nu);
        f= -pow(t1,1/(*nu))+(2/(*nu)-2)*log(t1)+(*nu-1)*log(log(dat[0])*log(dat[1]))-log(dat[0]*dat[1])+log(1+(*nu-1)*pow(t1,-1.0/(*nu)));
        if(f>XINFMAX) ll += log(XINFMAX);
        else if(f<log(DBL_MIN))ll += log(DBL_MIN);
        else ll += f;
      }
    }else{
      
      double *param, *fuc;
      param=Calloc(2,double);
      param[0]=*theta;
      param[1]=*nu;
      fuc = Calloc(*n,double);
      int k=1;
      for(j=0;j<*n;j++)
      {
        dat[0] = 1-u[j]; dat[1] = 1-v[j];
        
        dbb1(&dat[0], &dat[1], &k, param, &fuc[j]);
        
        if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
        else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
        else ll += log(fuc[j]);
      }
      Free(fuc); Free(param);
    }
  }
  else if(*family==18) //rotated BB6
  {
    double *param, *fuc;
    param=Calloc(2,double);
    param[0]=*theta;
    param[1]=*nu;
    fuc = Calloc(*n,double);
    int k=1;
    for(j=0;j<*n;j++)
    {
      dat[0] = 1-u[j]; dat[1] = 1-v[j];
      
      dbb6(&dat[0], &dat[1], &k, param, &fuc[j]);
      
      if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
      else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
      else ll += log(fuc[j]);
    }
    Free(fuc); Free(param);
  }
  else if(*family==19) //rotated BB7
  {
    if(*nu==0){
      for(j=0;j<*n;j++)
      {
        f = pow(pow(u[j],*theta)+pow(v[j],*theta)-pow(u[j],*theta)*pow(v[j],*theta),1/(*theta)-2)*pow(u[j],*theta-1)*pow(v[j],*theta-1)*(*theta-1+pow(u[j],*theta)+pow(v[j],*theta)-pow(u[j],*theta)*pow(v[j],*theta));
        if(log(f)>XINFMAX) ll += log(XINFMAX);
        else if(f < DBL_MIN) ll += log(DBL_MIN);
        else ll += log(f);
      }
    }else{
      double *param, *fuc;
      param=Calloc(2,double);
      param[0]=*theta;
      param[1]=*nu;
      fuc = Calloc(*n,double);
      int k=1;
      
      for(j=0;j<*n;j++)
      {
        dat[0] = 1-u[j]; dat[1] = 1-v[j];
        dbb7(&dat[0], &dat[1], &k, param, &fuc[j]);
        
        if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
        else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
        else ll += log(fuc[j]);
      }
      Free(fuc); Free(param);
    }
  }
  else if(*family==20) //rotated BB8
  {
    double *param, *fuc;
    param=Calloc(2,double);
    param[0]=*theta;
    param[1]=*nu;
    fuc = Calloc(*n,double);
    int k=1;
    
    for(j=0;j<*n;j++)
    {
      dat[0] = 1-u[j]; dat[1] = 1-v[j];
      dbb8(&dat[0], &dat[1], &k, param, &fuc[j]);
      
      if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
      else if(fuc[j]<DBL_MIN) ll += log(DBL_MIN);
      else ll += log(fuc[j]);
    }
    Free(fuc); Free(param);
  }
  else if(*family==41)		// New: 1-parametric asymmetric copula (from Harry Joe)
  {
    double tem1, tem2, con, sm, tem;
    
    for(j=0;j<*n;j++)
    {
      dat[0] = 1-u[j]; dat[1] = 1-v[j];
      tem1=qgamma(1.0-dat[0],*theta,1,1,0);
      tem2=qgamma(1.0-dat[1],*theta,1,1,0);
      con=gamma(1.0+(*theta))/(*theta);
      sm=pow(tem1,*theta)+pow(tem2,*theta);
      tem=pow(sm,(1.0/(*theta)));
      f=con*tem*exp(-tem+tem1+tem2)/sm;
      
      if(log(f)>XINFMAX) ll += log(XINFMAX);
      else if(f < DBL_MIN) ll += log(DBL_MIN);
      else ll += log(f);
    }
  }
  else if(*family==51)		// New: rotated 1-parametric asymmetric copula (from Harry Joe)
  {
    double tem1, tem2, con, sm, tem;
    
    for(j=0;j<*n;j++)
    {
      tem1=qgamma(1.0-u[j],*theta,1,1,0);
      tem2=qgamma(1.0-v[j],*theta,1,1,0);
      con=gamma(1.0+(*theta))/(*theta);
      sm=pow(tem1,*theta)+pow(tem2,*theta);
      tem=pow(sm,(1.0/(*theta)));
      f=con*tem*exp(-tem+tem1+tem2)/sm;
      
      if(log(f)>XINFMAX) ll += log(XINFMAX);
      else if(f < DBL_MIN) ll += log(DBL_MIN);
      else ll += log(f);
    }
  }
  else if(*family==104)		//New: Tawn
  {
    int T=1; //length of sample, different from T in TawnPDF
    double par3=1.0;
    for(j=0;j<*n;j++)
    {
      TawnPDF(&u[j], &v[j], &T, theta, nu, &par3, &f);
      if(log(f)>XINFMAX) ll += log(XINFMAX);
      else if(f < DBL_MIN) ll += log(DBL_MIN);
      else ll += log(f);
    }
  }
  else if(*family==114)		//New: rotated Tawn
  {
    int T=1; // length of sample, different from T in TawnPDF
    double par3=1.0;
    for(j=0;j<*n;j++)
    {
      dat[0] = 1-u[j]; dat[1] = 1-v[j];
      TawnPDF(&dat[0], &dat[1], &T, theta, nu, &par3, &f);
      if(log(f)>XINFMAX) ll += log(XINFMAX);
      else if(f < DBL_MIN) ll += log(DBL_MIN);
      else ll += log(f);
    }
  }
  else if(*family==204)		//New: Tawn2
  {
    int T=1; // length of sample, different from T in TawnPDF
    double par2=1.0;
    for(j=0;j<*n;j++)
    {
      TawnPDF(&u[j], &v[j], &T, theta, &par2, nu, &f);
      if(log(f)>XINFMAX) ll += log(XINFMAX);
      else if(f < DBL_MIN) ll += log(DBL_MIN);
      else ll += log(f);
    }
  }
  else if(*family==214)		//New: rotated Tawn2
  {
    int T=1; // length of sample, different from T in TawnPDF
    double par2=1.0;
    for(j=0;j<*n;j++)
    {
      dat[0] = 1-u[j]; dat[1] = 1-v[j];
      TawnPDF(&dat[0], &dat[1], &T, theta, &par2, nu, &f);
      if(log(f)>XINFMAX) ll += log(XINFMAX);
      else if(f < DBL_MIN) ll += log(DBL_MIN);
      else ll += log(f);
    }
  }
  else
  {
    Rprintf("%d\n",*family);
    Rcpp:: Rcout << "Error in LL: Unknown copula family" << std::endl;
  }
  //Free memory:
  Free(dat);
  //Write to output vector:
  *loglik = ll;
}

//////////////////////////////////////////////////////////////
// Function to compute likelihood for bivariate copula
// Input:
// family    copula family (0=independent, 1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank)
// n         sample size
// u         first variable of data set
// v         second variable of data set
// theta     dependency parameter
// nu        degrees-of-freedom for students copula
// coplik    output
//////////////////////////////////////////////////////////////
void copLik_mod(double* family, int* n, double* u, double* v, double* theta, double* nu, double* coplik)
{
  double* negv;
  double* negu;
  negv = (double *) malloc(*n*sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  double ntheta, nnu, nfamily;
  int i;
  ntheta = -*theta;
  nnu = -*nu;
  
  for(int i=0;i<*n;i++)
  {
    if(u[i]<UMIN) u[i]=UMIN;
    else if(u[i]>UMAX) u[i]=UMAX;
    if(v[i]<UMIN) v[i]=UMIN;
    else if(v[i]>UMAX) v[i]=UMAX;
  }
  
  if((*family)==43)
  {
    nfamily=3;
    if(*theta > 0){
      ntheta=2*(*theta)/(1-*theta);
      copLik(&nfamily, n, u,  v, &ntheta, nu, coplik);
    }else{
      ntheta=-2*(*theta)/(1+*theta);
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
      copLik(&nfamily, n, u,  negv, &ntheta, &nnu, coplik);
    }
  }else if((*family)==44)
  {
    nfamily=4;
    if(*theta > 0){
      ntheta=1/(1-*theta);
      copLik(&nfamily, n, u,  v, &ntheta, nu, coplik);
    }else{
      ntheta=1/(1+*theta);
      for (i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
      copLik(&nfamily, n, u,  negv, &ntheta, &nnu, coplik);
    }
  }else{
    
    if(((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30)) )	// 90? rotated copulas
    {
      nfamily = (*family)-20;
      for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
      copLik(&nfamily, n, u,  negv, &ntheta, &nnu, coplik);
    }
    else if(((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40)) )	// 270? rotated copulas
    {
      nfamily = (*family)-30;
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      copLik(&nfamily, n, negu,  v, &ntheta, &nnu, coplik);
    }
    else {
      copLik(family, n, u,  v, theta, nu, coplik);
    }
  }
  free(negv);
  free(negu);
}


void copLik(double* family, int* n, double* u, double* v, double* theta, double* nu, double* coplik)
{
  int j;
  double *dat, rho, lik=1.0, t1=0.0, t2=0.0, f;
  //Allocate memory:
  dat = Calloc(2,double);
  
  for(int i=0;i<*n;i++)
  {
    if(u[i]<UMIN) u[i]=UMIN;
    else if(u[i]>UMAX) u[i]=UMAX;
    if(v[i]<UMIN) v[i]=UMIN;
    else if(v[i]>UMAX) v[i]=UMAX;
  }
  
  //Compute likelihood:
  if(*family==0) //independent
    lik = 1.0;
  else if(*family==1) //Gaussian
  {
    rho=*theta;
    for(j=0;j<*n;j++)
    {
      dat[0]=u[j]; dat[1]=v[j];
      t1 = qnorm(dat[0],0.0,1.0,1,0); t2 = qnorm(dat[1],0.0,1.0,1,0);
      f = 1.0/sqrt(1.0-pow(rho,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*rho*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(rho,2.0))));
      lik *= f;
    }
  }
  else if(*family==2) //Student
  {
    rho=*theta;
    for(j=0;j<*n;j++)
    {
      dat[0] = u[j]; dat[1] = v[j];
      t1 = qt(dat[0],*nu,1,0); t2 = qt(dat[1],*nu,1,0);
      f = StableGammaDivision((*nu+2.0)/2.0,*nu/2.0)/(*nu*pi*sqrt(1.0-pow(rho,2.0))*dt(t1,*nu,0)*dt(t2,*nu,0))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho*t1*t2)/(*nu*(1.0-pow(rho,2.0))),-(*nu+2.0)/2.0);
      lik *= f;
    }
  }
  else if(*family==3) //Clayton
  {
    if(*theta == 0) lik = 1.0;
    if(*theta < XEPS) lik = 1.0;
    else
    {
      for(j=0;j<*n;j++)
      {
        dat[0] = u[j]; dat[1] = v[j];
        f = (1.0+*theta)*pow(dat[0]*dat[1],-1.0-*theta)*pow(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0,-2.0-1.0/(*theta));
        f = MAX(f,0);
        lik *= f;
      }
    }
  }
  else if(*family==4) //Gumbel
  {
    for(j=0;j<*n;j++)
    {
      dat[0] = u[j]; dat[1] = v[j];
      t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
      t2 = exp(-pow(t1,1.0/(*theta)));
      f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*theta))*pow(log(dat[0])*log(dat[1]),*theta-1.0)*(1.0+(*theta-1.0)*pow(t1,-1.0/(*theta)));
      lik *= f;
    }
  }
  else if(*family==5) // Frank
  {
    for(j=0;j<*n;j++)
    {
      dat[0] = u[j]; dat[1] = v[j];
      f = (*theta*(exp(*theta)-1.0)*exp(*theta*dat[1]+*theta*dat[0]+*theta))/pow(exp(*theta*dat[1]+*theta*dat[0])-exp(*theta*dat[1]+*theta)-exp(*theta*dat[0]+*theta)+exp(*theta),2.0);
      lik *= f;
    }
  }
  else if(*family==6)	//Joe
  {
    for(j=0;j<*n;j++)
    {
      f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
      lik *= f;
    }
  }
  else if(*family==7)	//BB1
  {
    if(*theta == 0){
      
      for(j=0;j<*n;j++)
      {
        dat[0] = u[j]; dat[1] = v[j];
        t1 = pow(-log(dat[0]),*nu)+pow(-log(dat[1]),*nu);
        t2 = exp(-pow(t1,1.0/(*nu)));
        f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*nu))*pow(log(dat[0])*log(dat[1]),*nu-1.0)*(1.0+(*nu-1.0)*pow(t1,-1.0/(*nu)));
        lik *= f;
      }
      
    }else{
      
      double *param, *fuc;
      param=Calloc(2,double);
      param[0]=*theta;
      param[1]=*nu;
      fuc = Calloc(*n,double);
      dbb1(u, v, n, param, fuc);
      for(j=0;j<*n;j++)
      {
        lik *= fuc[j];
      }
      Free(fuc); Free(param);
    }
  }
  else if(*family==8)	//BB6
  {
    double *param, *fuc;
    param=Calloc(2,double);
    param[0]=*theta;
    param[1]=*nu;
    fuc = Calloc(*n,double);
    dbb6(u, v, n, param, fuc);
    for(j=0;j<*n;j++)
    {
      lik *= fuc[j];
    }
    Free(fuc); Free(param);
  }
  else if(*family==9)	//BB7
  {
    if(*nu == 0){
      for(j=0;j<*n;j++)
      {
        f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
        lik *= f;
      }
    }else{
      
      double *param, *fuc;
      param=Calloc(2,double);
      param[0]=*theta;
      param[1]=*nu;
      fuc = Calloc(*n,double);
      dbb7(u, v, n, param, fuc);
      for(j=0;j<*n;j++)
      {
        lik *= fuc[j];
      }
      Free(fuc); Free(param);
    }
  }
  else if(*family==10)	//BB8
  {
    double *param, *fuc;
    param=Calloc(2,double);
    param[0]=*theta;
    param[1]=*nu;
    fuc = Calloc(*n,double);
    dbb8(u, v, n, param, fuc);
    for(j=0;j<*n;j++)
    {
      lik *= fuc[j];
    }
    Free(fuc); Free(param);
  }
  else if(*family==13) //rotated Clayton (180?)
  {
    if(*theta == 0) lik = 1.0;
    else
    {
      for(j=0;j<*n;j++)
      {
        dat[0] = 1-u[j]; dat[1] = 1-v[j];
        f = (1.0+*theta)*pow(dat[0]*dat[1],-1.0-*theta)*pow(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0,-2.0-1.0/(*theta));
        lik *= f;
      }
    }
  }
  else if(*family==14) //rotated Gumbel (180?)
  {
    for(j=0;j<*n;j++)
    {
      dat[0] = 1-u[j]; dat[1] = 1-v[j];
      t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
      t2 = exp(-pow(t1,1.0/(*theta)));
      f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*theta))*pow(log(dat[0])*log(dat[1]),*theta-1.0)*(1.0+(*theta-1.0)*pow(t1,-1.0/(*theta)));
      lik *= f;
    }
  }
  else if(*family==16) //rotated Joe (180?)
  {
    for(j=0;j<*n;j++)
    {
      dat[0]=1-u[j]; dat[1]=1-v[j];
      f = pow(pow(1-dat[0],*theta)+pow(1-dat[1],*theta)-pow(1-dat[0],*theta)*pow(1-dat[1],*theta),1/(*theta)-2)*pow(1-dat[0],*theta-1)*pow(1-dat[1],*theta-1)*(*theta-1+pow(1-dat[0],*theta)+pow(1-dat[1],*theta)-pow(1-dat[0],*theta)*pow(1-dat[1],*theta));
      lik *= f;
      
    }
  }
  else if(*family==17)	//rotated BB1
  {
    if(*theta == 0){
      
      for(j=0;j<*n;j++)
      {
        dat[0] = 1-u[j]; dat[1] = 1-v[j];
        t1 = pow(-log(dat[0]),*nu)+pow(-log(dat[1]),*nu);
        t2 = exp(-pow(t1,1.0/(*nu)));
        f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*nu))*pow(log(dat[0])*log(dat[1]),*nu-1.0)*(1.0+(*nu-1.0)*pow(t1,-1.0/(*nu)));
        lik *= f;
      }
      
    }else{
      
      double *param, *fuc;
      param=Calloc(2,double);
      param[0]=*theta;
      param[1]=*nu;
      fuc = Calloc(*n,double);
      int k=1;
      for(j=0;j<*n;j++)
      {
        dat[0] = 1-u[j]; dat[1] = 1-v[j];
        
        dbb1(&dat[0], &dat[1], &k, param, &fuc[j]);
        
        lik *= fuc[j];
      }
      Free(fuc); Free(param);
    }
  }
  else if(*family==18)	//rotated BB6
  {
    double *param, *fuc;
    param=Calloc(2,double);
    param[0]=*theta;
    param[1]=*nu;
    fuc = Calloc(*n,double);
    int k=1;
    for(j=0;j<*n;j++)
    {
      dat[0] = 1-u[j]; dat[1] = 1-v[j];
      
      dbb6(&dat[0], &dat[1], &k, param, &fuc[j]);
      
      lik *= fuc[j];
    }
    Free(fuc); Free(param);
  }
  else if(*family==19)	//rotated BB7
  {
    if(*nu == 0){
      for(j=0;j<*n;j++)
      {
        f = pow(pow(u[j],*theta)+pow(v[j],*theta)-pow(u[j],*theta)*pow(v[j],*theta),1/(*theta)-2)*pow(u[j],*theta-1)*pow(v[j],*theta-1)*(*theta-1+pow(u[j],*theta)+pow(v[j],*theta)-pow(u[j],*theta)*pow(v[j],*theta));
        lik *= f;
      }
    }else{
      
      double *param, *fuc;
      param=Calloc(2,double);
      param[0]=*theta;
      param[1]=*nu;
      fuc = Calloc(*n,double);
      int k=1;
      
      for(j=0;j<*n;j++)
      {
        dat[0] = 1-u[j]; dat[1] = 1-v[j];
        
        dbb7(&dat[0], &dat[1], &k, param, &fuc[j]);
        
        lik *= fuc[j];
      }
      Free(fuc); Free(param);
    }
  }
  else if(*family==20)	//rotated BB8
  {
    double *param, *fuc;
    param=Calloc(2,double);
    param[0]=*theta;
    param[1]=*nu;
    fuc = Calloc(*n,double);
    int k=1;
    
    for(j=0;j<*n;j++)
    {
      dat[0] = 1-u[j]; dat[1] = 1-v[j];
      
      dbb8(&dat[0], &dat[1], &k, param, &fuc[j]);
      
      lik *= fuc[j];
    }
    Free(fuc); Free(param);
  }
  
  else Rcpp::Rcout << "Error: Unknown Copula Family" << std::endl;
  //Free memory:
  Free(dat);
  //Write to output vector:
  *coplik = lik;
}









//// log likelihood for each observation --------
// unvectorized version
void LL_mod_seperate(double* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
  int kk=1;
  for(int i=0; i<(*n); i++){
    LL_mod2(family,&kk,&u[i],&v[i],theta,nu,&loglik[i]);
  };
}

// vectorized version
void LL_mod_seperate_vec(double* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
  int nn=1;
  for(int i=0; i<(*n); i++){
    LL_mod2(&family[i],&nn,&u[i],&v[i],&theta[i],&nu[i],&loglik[i]);
  };
}

//// density for each observation --------
// unvectorized version
void PDF_seperate(double* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
  int kk=1;
  for(int i=0; i<(*n); i++){
    LL_mod2(family,&kk,&u[i],&v[i],theta,nu,&loglik[i]);
    loglik[i] = exp(loglik[i]);
  };
}

// vectorized version
void PDF_seperate_vec(double* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
  int nn=1;
  for(int i=0; i<(*n); i++){
    LL_mod2(&family[i],&nn,&u[i],&v[i],&theta[i],&nu[i],&loglik[i]);
    loglik[i] = exp(loglik[i]);
  };
}

// h-func for BB1

void pcondbb1(double* u, double* v, int* n, double* param, double* out)
{
  int i;
  double th, de;
  double t1, t2, t3, t16, t17, t4, t5, t6, t7, t9, t10, t12, t13, t20;
  
  th = param[0];
  de = param[1];
  for(i=0;i<*n;i++)
  {
    t1 = pow(u[i],-th);
    t2 = t1-1.;
    t3 = pow(t2,de);
    t16 = 1./u[i];
    t17 = 1./t2;
    t4 = pow(v[i],-th);
    t5 = t4-1.;
    t6 = pow(t5,de);
    t7 = t3+t6;
    t9 = pow(t7,1/de);
    t10 = 1.0+t9;
    t12 = pow(t10,-1/th);
    t13 = t12*t9;
    t20 = 1./t10;
    out[i] = t13*t3*t1*t16*t17/t7*t20;
  }
  
}



void pcondbb6(double* u, double* v, int* n, double* param, double* out)
{
  int i;
  double th, de;
  double t1, t2, t3, t4, t5, t12, t16, t6, t7, t8, t9, t10, t11, t13, t14, t15, t17;
  
  th = param[0];
  de = param[1];
  
  for(i=0;i<*n;i++)
  {
    t1 = 1.0-u[i];
    t2 = pow(t1,th);
    t3 = 1.0-t2;
    t4 = log(t3);
    t5 = pow(-t4,de);
    t12 = 1/de;
    t16 = 1/th;
    t6 = 1.0-v[i];
    t7 = pow(t6,th);
    t8 = 1.0-t7;
    t9 = log(t8);
    t10 = pow(-t9,de);
    t11 = t5+t10;
    t13 = pow(t11,t12);
    t14 = exp(-t13);
    t15 = 1.0-t14;
    t17 = pow(t15,t16);
    
    out[i] = -t17*t13*t5*t2/t1/t3/t4/t11*t14/t15;
  }
  
}


void pcondbb7(double* u, double* v, int* n, double* param, double* out)
{
  int i;
  double th, de;
  double t1, t2, t3, t4, t6, t8, t9, t11, t12, t14;
  
  th = param[0];
  de = param[1];
  
  for(i=0;i<*n;i++)
  {
    t1 = 1.0-u[i];
    t2 = pow(t1,1.0*th);
    t3 = 1.0-t2;
    t4 = pow(t3,-1.0*de);
    t6 = pow(1.0-v[i],1.0*th);
    t8 = pow(1.0-t6,-1.0*de);
    t9 = t4+t8-1.0;
    t11 = pow(t9,-1.0/de);
    t12 = 1.0-t11;
    t14 = pow(t12,1.0/th);
    
    out[i] = t14*t11*t4*t2/t1/t3/t9/t12;
  }
  
}


void pcondbb8(double* u, double* v, int* n, double* param, double* out)
{
  int i;
  double th, de;
  double t2, t3, t12, t16, t6, t7, t8, t10, t11, t13, t15, t17;
  
  th = param[0];
  de = param[1];
  
  for(i=0;i<*n;i++)
  {
    t2 = 1.0-de*u[i];
    t3 = pow(t2,th);
    t10 = 1.0-de;
    t11 = pow(t10,th);
    t12 = 1.0-t11;
    t13 = 1/t12;
    t16 = 1/th;
    t6 = 1.0-de*v[i];
    t7 = pow(t6,th);
    t8 = 1.0-t7;
    t15 = 1.0-(1.0-t3)*t8*t13;
    t17 = pow(t15,t16);
    
    out[i] = t17*t3/t2*t8*t13/t15;
  }
  
}


// Since the h function is not symmetric in case of double Gumbel and double Clayton we have two implement both separately,
// i.e. Hfunc1 and Hfunc2
void  Hfunc1(double* family,int* n,double* u,double* v,double* theta,double* nu,double* out)
{
  double *negv, *negu;
  negv = (double *) malloc(*n* sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  double ntheta, nnu, nfamily;
  int j, T=1;
  ntheta = -*theta;
  nnu = -*nu;
  
  for(int i=0;i<*n;i++)
  {
    if(u[i]<UMIN) u[i]=UMIN;
    else if(u[i]>UMAX) u[i]=UMAX;
    if(v[i]<UMIN) v[i]=UMIN;
    else if(v[i]>UMAX) v[i]=UMAX;
  }
  
  if((*family)==43)
  {
    nfamily=3;
    if(*theta > 0){
      ntheta=2*(*theta)/(1-*theta);
      Hfunc(&nfamily, n, u, v, &ntheta, &nnu, out);
    }else{
      ntheta=-2*(*theta)/(1+*theta);
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hfunc(&nfamily, n, u, negv, &ntheta, &nnu, out);
    }
  }else if((*family)==44)
  {
    nfamily=4;
    if(*theta > 0){
      ntheta=1/(1-*theta);
      Hfunc (&nfamily, n, u, v, &ntheta, &nnu, out);
    }else{
      ntheta=1/(1+*theta);
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hfunc(&nfamily, n, u, negv, &ntheta, &nnu, out);
    }
  }else{
    
    if(((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30) | (*family==61) ))
    {
      nfamily=(*family)-20;
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hfunc(&nfamily, n, u, negv, &ntheta, &nnu, out);
    }
    else if(((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40) | (*family==71) ))
    {
      nfamily=(*family)-30;
      for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
      Hfunc(&nfamily, n, negu, v, &ntheta, &nnu, out);
      for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
    // u and v enter in wrong order from BiCopHfunc and have to be treated accordingly
    else if(*family==104)
    {
      double par3=1;
      dC_du(v,u,n,theta,nu,&par3,out);
    }
    else if(*family==114)
    {
      double par3=1;
      for(j=0;j<*n;j++)
      {
        negv[j]= 1-v[j];
        negu[j]= 1-u[j];
        dC_du(&negv[j],&negu[j],&T,theta,nu,&par3,&out[j]);
        out[j]= 1-out[j];
      }
    }
    else if(*family==124)
    {
      double par3=*nu;
      double par2=1;
      for(j=0;j<*n;j++)
      {
        negv[j]= 1-v[j];
        dC_du(&negv[j],&u[j],&T,&ntheta,&par2,&par3,&out[j]);
      }
      
    }
    else if(*family==134)
    {
      double par3=*nu;
      double par2=1;
      for(j=0;j<*n;j++)
      {
        negu[j]= 1-u[j];
        dC_du(&v[j],&negu[j],&T,&ntheta,&par2,&par3,&out[j]);
        out[j]=1-out[j];
      }
    }
    else if(*family==204)
    {
      double par3=*nu;
      double par2=1;
      dC_du(v,u,n,theta,&par2,&par3,out);
    }
    else if(*family==214)
    {
      double par3=*nu;
      double par2=1;
      for(j=0;j<*n;j++)
      {
        negv[j]= 1-v[j];
        negu[j]= 1-u[j];
        dC_du(&negv[j],&negu[j],&T,theta,&par2,&par3,&out[j]);
        out[j]= 1-out[j];
      }
    }
    else if(*family==224)
    {
      double par3=1;
      for(j=0;j<*n;j++)
      {
        negv[j]= 1-v[j];
        dC_du(&negv[j],&u[j],&T,&ntheta,nu,&par3,&out[j]);
      }
    }
    else if(*family==234)
    {
      double par3=1;
      for(j=0;j<*n;j++)
      {
        negu[j]= 1-u[j];
        dC_du(&v[j],&negu[j],&T,&ntheta,nu,&par3,&out[j]);
        out[j]=1-out[j];
        
      }
    }
    else {
      Hfunc (family, n, u, v, theta, nu, out);
    }
  }
  // ensure that results are in [0,1]
  for(int j=0; j <* n; ++j){out[j] = MIN(MAX(out[j], 0), 1);}
  free(negv);
  free(negu);
}

void  Hfunc2(double* family,int* n,double* v,double* u,double* theta,double* nu,double* out)
{
  double *negv, *negu;
  negv = (double *) malloc(*n * sizeof(double));
  negu = (double *) malloc(*n * sizeof(double));
  double ntheta, nnu, nfamily;
  ntheta = -*theta;
  nnu = -*nu;
  
  for(int i=0;i<*n;i++)
  {
    if(u[i]<UMIN) u[i]=UMIN;
    else if(u[i]>UMAX) u[i]=UMAX;
    if(v[i]<UMIN) v[i]=UMIN;
    else if(v[i]>UMAX) v[i]=UMAX;
  }
  
  if((*family)==43)
  {
    nfamily=3;
    if(*theta > 0){
      ntheta=2*(*theta)/(1-*theta);
      Hfunc (&nfamily, n, v, u, &ntheta, &nnu, out);
    }else{
      ntheta=-2*(*theta)/(1+*theta);
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hfunc(&nfamily, n, negv, u, &ntheta, &nnu, out);
      for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
  }
  else if((*family)==44)
  {
    nfamily=4;
    if(*theta > 0){
      ntheta=1/(1-*theta);
      Hfunc (&nfamily, n, v, u, &ntheta, &nnu, out);
    }else{
      ntheta=1/(1+*theta);
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hfunc(&nfamily, n, negv, u, &ntheta, &nnu, out);
      for (int i = 0; i < *n; i++) {out[i]=1-out[i];};		}
  }else{
    
    if(((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30) | (*family==61) ))
    {
      nfamily=(*family)-20;
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hfunc(&nfamily, n, negv, u, &ntheta, &nnu, out);
      for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
    else if(((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40) | (*family==71) ))
    {
      nfamily=(*family)-30;
      for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
      Hfunc(&nfamily, n, v, negu, &ntheta, &nnu, out);
    }
    else if((*family==104) | (*family==204) | (*family==114) | (*family==214))
    {
      // switch u and v and change type
      if((*family)/100 == 1) nfamily = (*family) + 100;
      if((*family)/100 == 2) nfamily = (*family) - 100;
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      Hfunc1(&nfamily, n, v, u, theta, nu, out);
      
    }
    else if((*family==124) | (*family==224) | (*family==134) | (*family==234))
    {
      // switch u and v and change type
      if((*family)/100 == 1) nfamily = (*family) + 100;
      if((*family)/100 == 2) nfamily = (*family) - 100;
      for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      Hfunc1(&nfamily, n, negv, negu, theta, nu, out);
      for (int i = 0; i < *n; i++) {out[i] = 1 - out[i];};
    }
    else
    {
      // switch u and v
      Hfunc(family, n, v, u, theta, nu, out);
    }
  }
  // ensure that results are in [0,1]
  for(int i=0; i < *n; ++i) {out[i] = MIN(MAX(out[i], UMIN), UMAX);}
  free(negv);
  free(negu);
}

// vectorized versions
void Hfunc1_vec(double* family,int* n,double* u,double* v,double* theta,double* nu,double* out)
{
  int nn=1;
  for(int i=0; i<(*n); i++){
    Hfunc1(&family[i], &nn, &u[i], &v[i], &theta[i], &nu[i], &out[i]);
  };
}

void Hfunc2_vec(double* family,int* n,double* u,double* v,double* theta,double* nu,double* out)
{
  int nn=1;
  for(int i=0; i<(*n); i++){
    Hfunc2(&family[i], &nn, &u[i], &v[i], &theta[i], &nu[i], &out[i]);
  };
}


//////////////////////////////////////////////////////////////
// Function to compute h-function for vine simulation and estimation
// Input:
// family   copula family (0=independent,  1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe, 7=BB1, 8=BB7)
// n        number of iterations
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// nu       degrees-of-freedom for the students copula
// out      output
//////////////////////////////////////////////////////////////
void Hfunc(double* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
  int j;
  double *h;
  h = Calloc(*n,double);
  double x;
  
  
  for(j=0;j<*n;j++)
  {
    if((v[j]==0) | ( u[j]==0)) h[j] = 0;
    else if (v[j]==1) h[j] = u[j];
    else
    {
      if(*family==0) //independent
      {
        h[j] = u[j];
      }
      else if(*family==1) //gaussian
      {
        x = (qnorm(u[j],0.0,1.0,1,0) - *theta*qnorm(v[j],0.0,1.0,1,0))/sqrt(1.0-pow(*theta,2.0));
        h[j] = pnorm(x,0.0,1.0,1,0);
      }
      else if(*family==2) //student
      {
        double t1, t2, mu, sigma2;
        t1 = qt(u[j],*nu,1,0); t2 = qt(v[j],*nu,1,0); mu = *theta*t2; sigma2 = ((*nu+t2*t2)*(1.0-*theta*(*theta)))/(*nu+1.0);
        h[j] = pt((t1-mu)/sqrt(sigma2),*nu+1.0,1,0);
      }
      else if(*family==3) //clayton
      {
        if(*theta == 0) h[j] = u[j] ;
        if(*theta < XEPS) h[j] = u[j] ;
        else
        {
          x = pow(u[j],-*theta)+pow(v[j],-*theta)-1.0 ;
          h[j] =   pow(v[j],-*theta-1.0)*pow(x,-1.0-1.0/(*theta));
          if(*theta < 0)
          {
            if(x < 0) h[j] = 0;
          }
        }
      }
      else if(*family==4) //gumbel
      {
        if(*theta == 1) h[j] = u[j] ;
        else
        {
          h[j] = -(exp(-pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)))*pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)-1.0)*pow(-log(v[j]),*theta))/(v[j]*log(v[j]));
        }
      }
      else if(*family==5) //frank
      {
        if(*theta==0) h[j]=u[j];
        else
        {
          h[j] = -(exp(*theta)*(exp(*theta*u[j])-1.0))/(exp(*theta*v[j]+*theta*u[j])-exp(*theta*v[j]+*theta)-exp(*theta*u[j]+*theta)+exp(*theta));
        }
      }
      else if(*family==6) //joe
      {
        if(*theta==1) h[j]=u[j];
        else
        {
          h[j] = pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
        }
      }
      else if(*family==7)	//BB1
      {
        double* param;
        param = Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        int T=1;
        if(*nu==1)
        {
          if(*theta==0) h[j]=u[j];
          else h[j]=pow(pow(u[j],-*theta)+pow(v[j],-*theta)-1,-1/(*theta)-1)*pow(v[j],-*theta-1);
        }
        else if(*theta==0)
        {
          h[j]=-(exp(-pow(pow(-log(v[j]),*nu)+pow(-log(u[j]),*nu),1.0/(*nu)))*pow(pow(-log(v[j]),*nu)+pow(-log(u[j]),*nu),1.0/(*nu)-1.0)*pow(-log(v[j]),*nu))/(v[j]*log(v[j]));
        }
        else
        {
          pcondbb1(&v[j],&u[j],&T,param,&h[j]);
        }
        Free(param);
      }
      else if(*family==8) //BB6
      {
        double* param;
        param = Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        int T=1;
        if(*theta==1)
        {
          if(*nu==1) h[j]=u[j];
          else h[j]=-(exp(-pow(pow(-log(v[j]),*nu)+pow(-log(u[j]),*nu),1.0/(*nu)))*pow(pow(-log(v[j]),*nu)+pow(-log(u[j]),*nu),1.0/(*nu)-1.0)*pow(-log(v[j]),*nu))/(v[j]*log(v[j]));
        }
        else if(*nu==1)
        {
          h[j]=pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
        }
        else
        {
          pcondbb6(&v[j],&u[j],&T,param,&h[j]);
        }
        Free(param);
      }
      else if(*family==9)	//BB7
      {
        double* param;
        param = Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        int T=1;
        if(*theta==1)
        {
          if(*nu==0) h[j]=u[j];
          else h[j]=pow(pow(u[j],-*nu)+pow(v[j],-*nu)-1,-1/(*nu)-1)*pow(v[j],-*nu-1);
        }
        else if(*nu==0)
        {
          h[j] = pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
        }
        else
        {
          pcondbb7(&v[j],&u[j],&T,param,&h[j]);
        }
        Free(param);
      }
      else if(*family==10) //BB8
      {
        double* param;
        param = Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        int T=1;
        if(*nu==0)
        {
          h[j]=u[j];
        }
        else if(*nu==1)
        {
          if(*theta==1) h[j]=u[j];
          else h[j]=pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
        }
        else
        {
          pcondbb8(&v[j],&u[j],&T,param,&h[j]);
        }
        Free(param);
      }
      else if(*family==13) //rotated clayton (180?)
      {
        if(*theta == 0) h[j] = u[j] ;
        if(*theta < XEPS) h[j] = u[j] ;
        else
        {
          u[j]=1-u[j];
          v[j]=1-v[j];
          x = pow(u[j],-*theta)+pow(v[j],-*theta)-1.0 ;
          h[j] =   pow(v[j],-*theta-1.0)*pow(x,-1.0-1.0/(*theta)); // pow(v[j],-*theta-1.0)*pow(pow(u[j],-*theta)+pow(v[j],-*theta)-1.0,-1.0-1.0/(*theta));
          h[j]= 1-h[j];
          u[j]=1-u[j];
          v[j]=1-v[j];
        }
      }
      else if(*family==14) //rotated gumbel (180?)
      {
        v[j]= 1-v[j];
        u[j]= 1-u[j];
        h[j]= -(exp(-pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)))*pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)-1.0)*pow(-log(v[j]),*theta))/(v[j]*	log(v[j]));
        h[j]= 1-h[j];
        u[j]=1-u[j];
        v[j]=1-v[j];
      }
      else if(*family==16)
      {
        v[j]= 1-v[j];
        u[j]= 1-u[j];
        h[j] = pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
        h[j]= 1-h[j];
        u[j]=1-u[j];
        v[j]=1-v[j];
      }
      else if(*family==17) //rotated BB1
      {
        double* param;
        param = Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        int T=1;
        if(*nu==1)
        {
          if(*theta==0) h[j]=u[j];
          else
          {
            h[j]=pow(pow(1-u[j],-*theta)+pow(1-v[j],-*theta)-1,-1/(*theta)-1)*pow(1-v[j],-*theta-1);
            h[j]= 1-h[j];
          }
        }
        else if(*theta==0)
        {
          h[j]=-(exp(-pow(pow(-log(1-v[j]),*nu)+pow(-log(1-u[j]),*nu),1.0/(*nu)))*pow(pow(-log(1-v[j]),*nu)+pow(-log(1-u[j]),*nu),1.0/(*nu)-1.0)*pow(-log(1-v[j]),*nu))/((1-v[j])*log(1-v[j]));
          h[j]= 1-h[j];
        }
        else
        {
          v[j]= 1-v[j];
          u[j]= 1-u[j];
          pcondbb1(&v[j],&u[j],&T,param,&h[j]);
          u[j]=1-u[j];
          v[j]=1-v[j];
          h[j]= 1-h[j];
        }
        Free(param);
      }
      else if(*family==18) //rotated BB6
      {
        double* param;
        param = Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        int T=1;
        if(*theta==1)
        {
          if(*nu==1) h[j]=u[j];
          else
          {
            h[j]=-(exp(-pow(pow(-log(1-v[j]),*nu)+pow(-log(1-u[j]),*nu),1.0/(*nu)))*pow(pow(-log(1-v[j]),*nu)+pow(-log(1-u[j]),*nu),1.0/(*nu)-1.0)*pow(-log(1-v[j]),*nu))/((1-v[j])*log(1-v[j]));
            h[j]= 1-h[j];
          }
        }
        else if(*nu==1)
        {
          h[j]=pow(pow(u[j],*theta) + pow(v[j],*theta) - pow(u[j],*theta)*pow(v[j],*theta),1.0/(*theta)-1) * pow(v[j],*theta-1.0)*(1-pow(u[j],*theta));
          h[j]= 1-h[j];
        }
        else
        {
          v[j]= 1-v[j];
          u[j]= 1-u[j];
          pcondbb6(&v[j],&u[j],&T,param,&h[j]);
          u[j]=1-u[j];
          v[j]=1-v[j];
          h[j]= 1-h[j];
        }
        Free(param);
      }
      else if(*family==19) //rotated BB7
      {
        double* param;
        param = Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        int T=1;
        if(*theta==1)
        {
          if(*nu==0) h[j]=u[j];
          else{
            h[j]=pow(pow(1-u[j],-*nu)+pow(1-v[j],-*nu)-1,-1/(*nu)-1)*pow(1-v[j],-*nu-1);
            h[j]= 1-h[j];
          }
        }
        else if(*nu==0)
        {
          h[j] = pow(pow(u[j],*theta) + pow(v[j],*theta) - pow(u[j],*theta)*pow(v[j],*theta),1.0/(*theta)-1) * pow(v[j],*theta-1.0)*(1-pow(u[j],*theta));
          h[j]= 1-h[j];
        }
        else
        {
          v[j]= 1-v[j];
          u[j]= 1-u[j];
          pcondbb7(&v[j],&u[j],&T,param,&h[j]);
          u[j]=1-u[j];
          v[j]=1-v[j];
          h[j]= 1-h[j];
        }
        Free(param);
      }
      else if(*family==20) //rotated BB8
      {
        double* param;
        param = Calloc(2,double);
        param[0]=*theta;
        param[1]=*nu;
        int T=1;
        if(*nu==0)
        {
          h[j]=u[j];
        }
        else if(*nu==1)
        {
          if(*theta==1) h[j]=u[j];
          else{
            h[j]=pow(pow(u[j],*theta) + pow(v[j],*theta) - pow(u[j],*theta)*pow(v[j],*theta),1.0/(*theta)-1) * pow(v[j],*theta-1.0)*(1-pow(u[j],*theta));
            h[j]= 1-h[j];
          }
        }
        else
        {
          v[j]= 1-v[j];
          u[j]= 1-u[j];
          pcondbb8(&v[j],&u[j],&T,param,&h[j]);
          u[j]=1-u[j];
          v[j]=1-v[j];
          h[j]= 1-h[j];
        }
        Free(param);
      }
      else if(*family==41)
      {
        double t1,t2,t3;
        t1=qgamma(1.0-u[j],*theta,1,1,0);
        t2=qgamma(1.0-v[j],*theta,1,1,0);
        t3=pow(pow(t1,*theta)+pow(t2,*theta),(1.0/(*theta)));
        h[j]=exp(-t3+t1);
      }
    }
    out[j] = MAX(MIN(h[j],UMAX),UMIN);
  }
  Free(h);
}

///////////////////////////////////////////////////////////////
void qcondgum(double* q, double* u, double* de, double* out)
{
  double a,p,g,gp,z1,z2,con,de1,dif;
  double mxdif;
  int iter;
  
  p = 1-*q;
  z1 = -log(*u);
  con=log(1.-p)-z1+(1.-*de)*log(z1); de1=*de-1.;
  a=pow(2.*pow(z1,*de),1./(*de));
  mxdif=1; iter=0;
  dif=.1;  // needed in case first step leads to NaN
  while ((mxdif > 1.e-6) && (iter < 20))
  {
    g=a+de1*log(a)+con;
    gp=1.+de1/a;
    dif=g/gp;
    a-=dif; iter++;
    int it = 0;
    while ((a <= z1) && (it < 20)) {
      dif /= 2.;
      a += dif;
      ++it;
    }
    mxdif=fabs(dif);
  }
  z2=pow(pow(a,*de)-pow(z1,*de),1./(*de));
  *out = exp(-z2);
}

void qcondjoe(double* q, double* u, double* de, double* out)
{ double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t13,t15,t16,t19,t23,t28,t31;
  double c21,pdf;
  int iter;
  double diff,v,de1,dtem,de1inv,tem;
  
  t1 = 1.0-*u;
  t2 = pow(t1,1.0*(*de));
  t7 = 1./(*de);
  t10 = t2*(*de);
  t11 = 1./t1;
  t19 = (*de)*(*de);
  de1=*de-1;  // may need better modification for large delta
  dtem=-de1/(1.+de1); de1inv=-1./de1;
  
  // v = 0.5 * (q+u); // starting guess
  
  // Use a better starting point based on reflected B4 copula
  // A good starting point is crucial when delta is large because
  //    C_{2|1} will be steep
  // C_{R,2|1}(v|u)=1-C_{2|1}(1-v|1-u),
  // C_{R,2|1}^{-1}(q|u)=1-C_{2|1}^{-1}(1-q|1-u)
  tem=pow(1.-*q,dtem)-1.;
  tem=tem*pow(1.-*u,-de1)+1.;
  v=pow(tem,de1inv); v=1.-v;
  diff=1; iter=0;
  while(fabs(diff)>1.e-6 && iter<20)
  { t3 = 1.-v;
    t4 = pow(t3,*de);
    t5 = t2*t4;
    t6 = t2+t4-t5;
    t8 = pow(t6,t7);
    t9 = t7*t8;
    t13 = t11*t4;
    t15 = -t10*t11+t10*t13;
    t16 = 1./t6;
    t23 = 1./t3;
    t28 = t6*t6;
    t31 = (-t4*(*de)*t23+t5*(*de)*t23)/t28*t15;
    c21 = -t9*t15*t16;
    pdf = -t8/t19*t31+t8*(*de)*t2*t13*t23*t16+t9*t31;
    iter++;
  
    diff=(c21-*q)/pdf;
    v-=diff;
    while(v<=0 || v>=1 || fabs(diff)>0.25 ) { diff/=2.; v+=diff; }
  }
  *out = v;
}




///////////////////////////////////////////////////////////////
// Function to compute inversion of H numerically through Bisection
// Input:
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// out      output
//////////////////////////////////////////////////////////////
void HNumInv(double* family, double* u, double* v, double* theta, double* nu, double* out)
{
  
  int br = 0, in = 1, it = 0;
  double ans = 0.0, tol = 1e-12, x0 = UMIN, x1 = UMAX, fl, fh, val;
  Hfunc1(family,&in,&x0,v,theta,nu,&fl);
  fl -= *u;
  Hfunc1(family,&in,&x1,v,theta,nu,&fh);
  fh -= *u;
  if (fabs(fl) <= tol) {
    ans = x0;
    br = 1;
  }
  if (fabs(fh) <= tol) {
    ans = x1;
    br = 1;
  }
  
  while (!br){
    ans = (x0 + x1) / 2.0;
    Hfunc1(family,&in,&ans,v,theta,nu,&val);
    val -= *u;
    
    //stop if values become too close (avoid infinite loop)
    if (fabs(val) <= tol) br = 1;
    if (fabs(x0-x1) <= tol) br = 1;
    
    if (val > 0.0) {
      x1 = ans;
      fh = val;
    } else {
      x0 = ans;
      fl = val;
    }
    
    //stop if too many iterations are required (avoid infinite loop)
    ++it;
    if (it > 50) br = 1;
  }
  
  *out = ans;
}

void HNumInv2(double* family, double* v, double* u, double* theta, double* nu, double* out)
{
  
  int br = 0, in = 1, it = 0;
  double ans = 0.0, tol = 1e-12, x0 = UMIN, x1 = UMAX, fl, fh, val;
  Hfunc2(family, &in, &x0, u, theta, nu, &fl);
  fl -= *v;
  Hfunc2(family, &in, &x1, u, theta, nu, &fh);
  fh -= *v;
  if (fabs(fl) <= tol) {
    ans = x0;
    br = 1;
  }
  if (fabs(fh) <= tol) {
    ans = x1;
    br = 1;
  }
  
  while (!br){
    ans = (x0 + x1) / 2.0;
    Hfunc2(family, &in, &ans, u, theta, nu, &val);
    val -= *v;
    
    //stop if values become too close (avoid infinite loop)
    if (fabs(val) <= tol) br = 1;
    if (fabs(x0-x1) <= tol) br = 1;
    
    if (val > 0.0) {
      x1 = ans;
      fh = val;
    } else {
      x0 = ans;
      fl = val;
    }
    
    //stop if too many iterations are required (avoid infinite loop)
    ++it;
    if (it > 50) br = 1;
  }
  *out = ans;
}

/////////////////////////////////////////////
// Function to invert h-function for vine simulation and estimation
/////////////////////////////////////////////
void Hinv1(double* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
  double *negv, *negu;
  negv = (double*) Calloc(*n,double);
  negu = (double*) Calloc(*n,double);
  double ntheta, nnu, nfamily;
  ntheta = -*theta;
  nnu = -*nu;
  
  for(int i=0;i<*n;i++)
  {
    if(u[i]<UMIN) u[i]=UMIN;
    else if(u[i]>UMAX) u[i]=UMAX;
    if(v[i]<UMIN) v[i]=UMIN;
    else if(v[i]>UMAX) v[i]=UMAX;
  }
  
  if((*family)==43)
  {
    nfamily=3;
    if(*theta > 0){
      ntheta=2*(*theta)/(1-*theta);
      Hinv(&nfamily, n, u, v, &ntheta, &nnu, out);
    }else{
      ntheta=-2*(*theta)/(1+*theta);
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hinv(&nfamily, n, u, negv, &ntheta, &nnu, out);
      for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
  }
  else if((*family)==44)
  {
    nfamily=4;
    if(*theta > 0){
      ntheta=1/(1-*theta);
      Hinv(&nfamily, n, u, v, &ntheta, &nnu, out);
    }else{
      ntheta=1/(1+*theta);
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hinv(&nfamily, n, u, negv, &ntheta, &nnu, out);
      for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
  }
  else if(((*family ==23) | (*family ==24) | (*family==26) | (*family ==27) | (*family ==28) | (*family==29) | (*family==30) | (*family==61) ))
  {
    nfamily=(*family)-20;
    for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
    Hinv(&nfamily,  n,  u,  negv,  &ntheta,  &nnu,  out);
  }
  else if(((*family==33) | (*family==34) | (*family==36) | (*family ==37) | (*family ==38) | (*family==39) | (*family==40) | (*family==71) ))
  {
    nfamily=(*family)-30;
    for (int i = 0; i < *n; i++) {negu[i]=1 - u[i];};
    Hinv(&nfamily,  n,  negu,  v,  &ntheta,  &nnu,  out);
    for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
  }
  else {
    Hinv(family,  n,  u,  v,  theta,  nu,  out);
  }
  Free(negv);
  Free(negu);
}

void Hinv2(double* family, int* n, double* v, double* u, double* theta, double* nu, double* out)
{
  double *negv, *negu;
  negv = (double *) malloc(*n*sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  double ntheta, nnu, nfamily;
  ntheta = -*theta;
  nnu = -*nu;
  
  for(int i=0;i<*n;i++)
  {
    if(u[i]<UMIN) u[i]=UMIN;
    else if(u[i]>UMAX) u[i]=UMAX;
    if(v[i]<UMIN) v[i]=UMIN;
    else if(v[i]>UMAX) v[i]=UMAX;
  }
  
  if((*family)==43)
  {
    nfamily=3;
    if(*theta > 0){
      ntheta=2*(*theta)/(1-*theta);
      Hinv(&nfamily, n, v, u, &ntheta, &nnu, out);
    }else{
      ntheta=-2*(*theta)/(1+*theta);
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hinv(&nfamily, n, negv, u, &ntheta, &nnu, out);
      for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
  }
  else if((*family)==44)
  {
    nfamily=4;
    if(*theta > 0){
      ntheta=1/(1-*theta);
      Hinv(&nfamily, n, v, u, &ntheta, &nnu, out);
    }else{
      ntheta=1/(1+*theta);
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hinv(&nfamily, n, negv, u, &ntheta, &nnu, out);
      for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
  }
  else if(((*family ==23) | (*family ==24) | (*family==26) | (*family ==27) | (*family ==28) | (*family==29) | (*family==30) | (*family==61) ))
  {
    nfamily = (*family)-20;
    for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
    Hinv(&nfamily,  n,  negv, u,  &ntheta,  &nnu,  out);
    for (int i = 0; i < *n; i++) {out[i] = 1 - out[i];};
  }
  else if(((*family==33) | (*family==34) | (*family==36) | (*family ==37) | (*family ==38) | (*family==39) | (*family==40) | (*family==71) ))
  {
    nfamily=(*family)-30;
    for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
    Hinv(&nfamily,  n,  v,  negu,  &ntheta,  &nnu,  out);
  }
  else if((*family==104) | (*family==204) | (*family==114) | (*family==214))
  {
    // change type
    //if((*family)/100 == 1) nfamily = (*family) + 100;
    //if((*family)/100 == 2) nfamily = (*family) - 100;
    for (int i = 0; i < *n; ++i) {
      HNumInv2(family, &v[i], &u[i], theta, nu, &out[i]);
    }
  }
  else if((*family==124) | (*family==224))
  {
    // change type
    //if((*family)/100 == 1) nfamily = (*family) + 100;
    //if((*family)/100 == 2) nfamily = (*family) - 100;
    //nfamily = nfamily + 10;
    for (int i = 0; i < *n; ++i) {
      HNumInv2(family, &v[i], &u[i], theta, nu, &out[i]);
    }
  }
  else if((*family==134) | (*family==234))
  {
    // change type
    //if((*family)/100 == 1) nfamily = (*family) + 100;
    //if((*family)/100 == 2) nfamily = (*family) - 100;
    //nfamily = nfamily - 10;
    for (int i = 0; i < *n; ++i) {
      HNumInv2(family, &v[i], &u[i], theta, nu, &out[i]);
    }
  }
  else {
    Hinv(family,  n,  v,  u,  theta,  nu,  out);
  }
  free(negv);
  free(negu);
}

// vectorized versions
void Hinv1_vec(double* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
  int nn=1;
  for(int i=0; i<(*n); i++){
    Hinv1(&family[i], &nn, &u[i], &v[i], &theta[i], &nu[i], &out[i]);
  };
}

void Hinv2_vec(double* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
  int nn=1;
  for(int i=0; i<(*n); i++){
    Hinv2(&family[i], &nn, &u[i], &v[i], &theta[i], &nu[i], &out[i]);
  };
}




//////////////////////////////////////////////////////////////
// Function to invert h-function for vine simulation and estimation
// Input:
// family   copula family (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe)
// n        number of iterations
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// nu       degrees-of-freedom for the students copula
//////////////////////////////////////////////////////////////
void Hinv(double* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
  int j;
  double *hinv;
  hinv = Calloc(*n,double);
  
  for(int i=0;i<*n;i++)
  {
    if(u[i]<UMIN) u[i]=UMIN;
    else if(u[i]>UMAX) u[i]=UMAX;
    if(v[i]<UMIN) v[i]=UMIN;
    else if(v[i]>UMAX) v[i]=UMAX;
  }
  
  for(j=0;j<*n;j++)
  {
    if(*family==0)
    {
      hinv[j]=u[j];
    }
    else if(*family==1) //gaussian
    {
      hinv[j] = pnorm(qnorm(u[j],0.0,1.0,1,0)*sqrt(1.0-pow(*theta,2.0))+*theta*qnorm(v[j],0.0,1.0,1,0),0.0,1.0,1,0);
    }
    else if(*family==2) //student
    {
      double temp1, temp2, mu, var;
      temp1 = qt(u[j],*nu+1.0,1,0); temp2 = qt(v[j],*nu,1,0); mu = *theta*temp2; var=((*nu+(temp2*temp2))*(1.0-(*theta*(*theta))))/(*nu+1.0);
      hinv[j] = pt((sqrt(var)*temp1)+mu,*nu,1,0);
    }
    else if(*family==3) //clayton
    {
      if(*theta < XEPS) {
        hinv[j]=u[j];
      } else if (*theta < 75) {
        hinv[j] = pow(pow(u[j]*pow(v[j],*theta+1.0),-*theta/(*theta+1.0))+1.0-pow(v[j],-*theta),-1.0/(*theta));
      } else {
        double nu=0.0;
        HNumInv(family,&u[j],&v[j],theta,&nu,&hinv[j]);
      }
    }
    else if(*family==4) //gumbel - must turn to numerical inversion
    {
      qcondgum(&u[j],&v[j],theta,&hinv[j]);
    }
    else if(*family==5) //frank - numerical inversion
    {
      //hinv[j] = -1/(*theta)*log(1-(1-exp(-*theta)) / ((1/u[j]-1)*exp(-*theta*v[j])+1));
      double nu=0.0;
      HNumInv(family,&u[j],&v[j],theta,&nu,&hinv[j]);
    }
    else if(*family==6) //joe - numerical inversion
    {
      if(*theta<40)
      {
        qcondjoe(&u[j],&v[j],theta,&hinv[j]);
      }
      else
      {
        double nu=0.0;
        HNumInv(family,&u[j],&v[j],theta,&nu,&hinv[j]);
      }
    }
    else if(*family==7) //BB1
    {
      HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
    }
    else if(*family==8) //BB6
    {
      HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
    }
    else if(*family==9) //BB7
    {
      HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
    }
    else if(*family==10) //BB8
    {
      HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
    }
    else if(*family==13)
    {
      u[j]=1-u[j];
      v[j]=1-v[j];
      hinv[j] = pow(pow(u[j]*pow(v[j],*theta+1.0),-*theta/(*theta+1.0))+1.0-pow(v[j],-*theta),-1.0/(*theta));
      hinv[j]=1-hinv[j];
      u[j]=1-u[j];
      v[j]=1-v[j];
    }
    else if(*family==14) //rotated gumbel (180?) - must turn to numerical inversion
    {
      u[j]=1-u[j];
      v[j]=1-v[j];
      qcondgum(&u[j],&v[j],theta,&hinv[j]);
      hinv[j]=1-hinv[j];
      u[j]=1-u[j];
      v[j]=1-v[j];
    }
    else if(*family==16) //rotated joe (180?) - must turn to numerical inversion
    {
      u[j]=1-u[j];
      v[j]=1-v[j];
      if(*theta<40)
      {
        qcondjoe(&u[j],&v[j],theta,&hinv[j]);
      }
      else
      {
        double jj=6;
        double nu=0.0;
        HNumInv(&jj,&u[j],&v[j],theta,&nu,&hinv[j]);
      }
      hinv[j]=1-hinv[j];
      u[j]=1-u[j];
      v[j]=1-v[j];
    }
    else if(*family==17) //rotated BB1 (180?) - must turn to numerical inversion
    {
      double jj=7;
      u[j]=1-u[j];
      v[j]=1-v[j];
      HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
      hinv[j]=1-hinv[j];
      u[j]=1-u[j];
      v[j]=1-v[j];
    }
    else if(*family==18) //rotated BB6 (180?) - must turn to numerical inversion
    {
      double jj=8;
      u[j]=1-u[j];
      v[j]=1-v[j];
      HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
      hinv[j]=1-hinv[j];
      u[j]=1-u[j];
      v[j]=1-v[j];
    }
    else if(*family==19) //rotated BB7 (180?) - must turn to numerical inversion
    {
      double jj=9;
      u[j]=1-u[j];
      v[j]=1-v[j];
      HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
      hinv[j]=1-hinv[j];
      u[j]=1-u[j];
      v[j]=1-v[j];
    }
    else if(*family==20) //rotated BB8 (180?) - must turn to numerical inversion
    {
      double jj=10;
      u[j]=1-u[j];
      v[j]=1-v[j];
      HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
      hinv[j]=1-hinv[j];
      u[j]=1-u[j];
      v[j]=1-v[j];
    }
    else if(*family==41)	// 1-parametric asymmetric copula (Harry Joe)
    {
      double de=*theta;
      double tem1,y;
      tem1=qgamma(1.0-v[j],de,1,1,0);
      y=pow(tem1-log(u[j]),de) - pow(tem1,de);
      y=pow(y,(1.0/de));
      hinv[j]=1.0-pgamma(y,de,1,1,0);
    }
    else if(*family==51)	// rotated (180) 1-parametric asymmetric copula (Harry Joe)
    {
      double de=*theta;
      double tem1,y;
      u[j]=1-u[j];
      v[j]=1-v[j];
      tem1=qgamma(1.0-v[j],de,1,1,0);
      y=pow(tem1-log(u[j]),de) - pow(tem1,de);
      y=pow(y,(1.0/de));
      hinv[j]=1.0-pgamma(y,de,1,1,0);
      hinv[j]=1-hinv[j];
      u[j]=1-u[j];
      v[j]=1-v[j];
    }
    else if(((*family)/100 == 1) | ((*family)/100 == 2)) //Tawn
    {
      HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
    }
    
    out[j] = MAX(MIN(hinv[j],UMAX),UMIN);
  }
  Free(hinv);
}


