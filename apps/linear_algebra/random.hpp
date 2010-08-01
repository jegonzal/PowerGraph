#ifndef MYMYRANDOM_H
#define MYMYRANDOM_H

//#include <iostream>
#include <cmath>

/*
#define FALSE 0
#define TRUE 1
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define IM1 2147483563
#define IM2 2147483399
#define AM1 (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
#define PI 3.141592654
#define E 2.71828182845904523536
#define IB1 1
#define IB2 2
#define IB5 16
#define IB18 131072
#define MASK (IB1+IB2+IB5)
#define NITER 4

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j >= 0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j >= 0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM1*iy) > RNMX) return RNMX;
	else return temp;
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
float ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i <= 54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k <= 4;k++)
			for (i=1;i <= 55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
void psdes(unsigned long *lword, unsigned long *irword)
{
	unsigned long i,ia,ib,iswap,itmph=0,itmpl=0;
	static unsigned long c1[NITER]={
		0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
	static unsigned long c2[NITER]={
		0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};

	for (i=0;i < NITER;i++) {
		ia=(iswap=(*irword)) ^ c1[i];
		itmpl = ia & 0xffff;
		itmph = ia >> 16;
		ib=itmpl*itmpl+ ~(itmph*itmph);
		*irword=(*lword) ^ (((ia = (ib >> 16) |
			((ib & 0xffff) << 16)) ^ c2[i])+itmpl*itmph);
		*lword=iswap;
	}
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */
/*
float ran4(long *idum)
{
	void psdes(unsigned long *lword, unsigned long *irword);
	unsigned long irword,itemp,lword;
	static long idums = 0;
#if defined(vax) || defined(_vax_) || defined(__vax__) || defined(VAX)
	static unsigned long jflone = 0x00004080;
	static unsigned long jflmsk = 0xffff007f;
#else
	static unsigned long jflone = 0x3f800000;
	static unsigned long jflmsk = 0x007fffff;
#endif

	if (*idum < 0) {
		idums = -(*idum);
		*idum=1;
	}
	irword=(*idum);
	lword=idums;
	psdes(&lword,&irword);
	itemp=jflone | (jflmsk & irword);
	++(*idum);
	return (*(float *)&itemp)-1.0;
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
float expdev(long *idum)
{
	float ran1(long *idum);
	float dum;

	do
		dum=ran1(idum);
	while (dum == 0.0);
	return -log(dum);
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
float gasdev(long *idum)
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

float gamdev(int ia)
{
	int j;
	float am,e,s,v1,v2,x,y;

	assert(ia>=1);
        if (ia < 6) {
		x=1.0;
		for (j=1;j <= ia;j++) x *= drand48();
		x = -log(x);
	} else {
		do {
			do {
				do {
					v1=drand48();
					v2=2.0*drand48()-1.0;
				} while (v1*v1+v2*v2 > 1.0);
				y=v2/v1;
				am=ia-1;
				s=sqrt(2.0*am+1.0);
				x=s*y+am;
			} while (x <= 0.0);
			e=(1.0+y*y)*exp(am*log(x/am)-s*y);
		} while (drand48() > e);
	}
	return x;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */
/*
float gammln(float xx)
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
	for (j=0;j <= 5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
float poidev(float xm, long *idum)
{
	float gammln(float xx);
	float ran1(long *idum);
	static float sq,alxm,g,oldm=(-1.0);
	float em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran1(idum);
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*ran1(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran1(idum) > t);
	}
	return em;
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
float bnldev(float pp, int n, long *idum)
{
	float gammln(float xx);
	float ran1(long *idum);
	int j;
	static int nold=(-1);
	float am,em,g,angle,p,bnl,sq,t,y;
	static float pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j <= n;j++)
			if (ran1(idum) < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j <= n;j++) {
			t *= ran1(idum);
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran1(idum);
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
				-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran1(idum) > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
int irbit1(unsigned long *iseed)
{
	unsigned long newbit;

	newbit = (*iseed & IB18) >> 17
		^ (*iseed & IB5) >> 4
		^ (*iseed & IB2) >> 1
		^ (*iseed & IB1);
	*iseed=(*iseed << 1) | newbit;
	return (int) newbit;
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
int irbit2(unsigned long *iseed)
{
	if (*iseed & IB18) {
		*iseed=((*iseed ^ MASK) << 1) | IB1;
		return 1;
	} else {
		*iseed <<= 1;
		return 0;
	}
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
float factrl(int n)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	static int ntop=4;
	static float a[33]={1.0,1.0,2.0,6.0,24.0};
	int j;

	if (n < 0) nrerror("Negative factorial in routine factrl");
	if (n > 32) return exp(gammln(n+1.0));
	while (ntop < n) {
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
float factln(int n)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	static float a[101];

	if (n < 0) nrerror("Negative factorial in routine factln");
	if (n <= 1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	else return gammln(n+1.0);
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
float bico(int n, int k)
{
	float factln(int n);

	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}*/
/* (C) Copr. 1986-92 Numerical Recipes Software 0NL-Z%. */

/*
float gs(float alpha, long *idum)
{
  float ran1(long *idum);
  float b,p,u,u1,x;
  int found;

  found=FALSE;
step1:
  u=ran1(idum);
    b=(E+alpha)/E;
  p=b*u;
  if (p > 1)
    goto step3;
step2:
  x=exp(log(p)/alpha);
  if (u1 > exp(-x))
    goto step1;
  else
    return x;
step3:
  x=-log((b-p)/alpha);
  u1=ran1(idum);
  if ( u1 > exp((alpha-1)*log(x)) )
    goto step1;
  else
    return x;
}*/
/* Algoritm GS from W J Kennedy & J E Gentle, Statistical Computing,
   New York, NY, and Basel: Dekker 1980, p. 213.  Original reference is
   J R Ahrens and U Dieter, Computer methods for sampling from gamma, beta,
   Poisson, and binomial distributions, Computing 12 (1974), 223-246 */

/*
float gf1(float alpha, long *idum)
{
  float ran1(long *idum);
  float y,z,w,x,u;
  int found;

  found=FALSE;
  while (found==FALSE)
    {
      y=gamdev(1);
      z=y/exp(y+1);
      w=exp((alpha-1)*log(z));
      x=alpha*y;
      u=ran1(idum);
      if (u <= w)
        found=TRUE;
    }
  return x;
}*/
/* Algoritm GF1 from W J Kennedy & J E Gentle, Statistical Computing, 
   New York, NY, and Basel: Dekker 1980, p. 214.  Original reference is G S 
   Fishman, Sampling from the gamma distribution on a computer, CACM
   (Communications of the Association for Computing Machinery) 19 (1976),
   407-409 */

/*
float g0(float alpha, long *idum)
{
  float ran1(long *idum);
  void nrerror(char error_text[]);
  double mu,v,sigma2,sigma,w,d,b,u0,logu0;
  float u,s,t,x;
  int found;

  found=FALSE;
step0:
  if (alpha <= 2.5327805161251) nrerror("Error in routine G0");
step1:
  mu=alpha-1;
  v=sqrt(alpha);
  sigma2=alpha+1.632993161855*v;
  sigma=sqrt(sigma2);
  w=sigma2/mu;
  d=2.44948974278318*sigma;
  b=mu+d;
step2:
  u=ran1(idum);
  u0=0.009572265238289;
  if (u < u0)
    goto step8;
step3:
  s=gasdev(idum);
  x=mu+sigma*s;
  if((x < 0)||(x > b))
    goto step2;
step4:
  u=ran1(idum);
  t=s*s/2;
step5:
  if (u <= 1 - t*((1-2*s/v)*w -1))
    return x;
  else
    goto step7;
step6:
  if (u <= 1 - t*(w-1))
    return x;
step7:
  if (log(u) <= mu*(1 + log(x/mu)) - x + t)
    return x;
  else
    goto step2;
step8:
  s=gamdev(1);
  x=b*(1 + s/d);
step9:
  u=ran1(idum);
  logu0 = mu*(2 + log(x/mu) - x/b) + 3.7203284924588 - b - log(sigma*d/b);
  if (log(u) > logu0)
    goto step2;
  else
    return x;
}*/
/* Algoritm G0 from W J Kennedy & J E Gentle, Statistical Computing,
   New York, NY, and Basel: Dekker 1980, p. 215.  Original reference is
   J H Ahrens and U Dieter, Computer methods for sampling from gamma, beta,
   Poisson and binomial distributions, Computing 12 (1974), 223-246. */
  
/*
float gammadev(float alpha, long *idum)
{
  float ran1(long *idum);
  float gs(float alpha, long *idum);
  float gf1(float alpha, long *idum);
  float x;
  int ia;
  
  ia=int(alpha);
  if (alpha==ia)
    x=gamdev(ia);
  else
    {
      if (alpha<1)
        x=gs(alpha,idum);
      else
        if (alpha<3)
          x=gf1(alpha,idum);
        else
          x=g0(alpha,idum);
    }
  return x;
}*/

/*
void nrerror(char error_text[])*/
/* Numerical Recipes standard error handler */
/*
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}*/


/*                             DISCLAIMER OF WARRANTY                     */
/*          THE PROGRAMS ACCESSED  BY THIS ROUTINE  (AND ON THE ORIGINAL  */
/*     DISKETTE)  ARE PROVIDED   "AS IS"   WITHOUT WARRANTY OF ANY KIND.  */
/*     WE MAKE NO WARRANTIES, EXPRESS OR IMPLIED,  THAT THEY ARE FREE OF  */
/*     ERROR,   OR  ARE  CONSISTENT  WITH  ANY  PARTICULAR  STANDARD  OF  */
/*     MERCHANTABILITY, OR THAT THEY WILL MEET YOUR REQUIREMENTS FOR ANY  */
/*     PARTICULAR APPLICATION.  THEY SHOULD NOT BE RELIED ON FOR SOLVING  */
/*     A PROBLEM  WHOSE INCORRECT SOLUTION  COULD  RESULT IN INJURY TO A  */
/*     PERSON OR LOSS OF PROPERTY.  IF YOU DO USE THEM IN SUCH A MANNER,  */
/*     IT IS  AT YOUR OWN RISK.  THE AUTHORS AND PUBLISHER  DISCLAIM ALL  */
/*     LIABILITY   FOR  DIRECT,   INDIRECT,   OR  CONSEQUENTIAL  DAMAGES  */
/*     RESULTING FROM YOUR USE OF THE PROGRAMS.                           */



#endif //MYMYRANDOM_H
