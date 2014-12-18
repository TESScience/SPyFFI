/* test_poidev_v1.c
 * Alan M. Levine
 * December 16, 2014
 * Heritage: cr_img4_z1.c, poidev2t.c
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "twister.h"
#include "seed_tw_ran.h"

/******************************************************************************/
// ZKBT says this function was copied and pasted from gasdev2t.c
double gasdev()
{
	static int iset=0;
	static double gset;
	double fac,r,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ranMT()-1.0;
			v2=2.0*ranMT()-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

/******************************************************************************/
double gammln(double xx)
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
  int j;

  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}

/******************************************************************************/
// From poidev2t.c

double poidev(double xm)
{
  static double sq, alxm, g, oldm=(-1.0);
  double em, t, y;

  if (xm < 12.0) {
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em = -1;
    t=1.0;
    do {
      em += 1.0;
      t *= ranMT();
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
	y=tan(M_PI*ranMT());
	em=sq*y+xm;
      } while (em < 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
    } while (ranMT() > t);
  }
  return em;
}

/******************************************************************************/
// Tasks that we don't want to do unnecessarily when doing multiple images

void setup()
{
  unsigned long seed;

  seed = get_tw_seed();
  seedMT(seed);
}

/******************************************************************************/

int main(int argc, char *argv[])
{
  double x, y;
  int i, j;

  setup();

  x = 0.01;
  for(j=0;j<20;++j) {
    printf("\n");
    for(i=0;i<20;i++) {
      y = poidev(x);
      printf("%f %f\n",x,y);
    }
    x *= 2.0;
  }
  
}

/******************************************************************************/
