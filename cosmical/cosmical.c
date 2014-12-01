/* cosmical[full]|[small].c
 * 
 * Alan M. Levine
 * July 14, 2014
 * Heritage: cr_img3.c
 *  
 *  on August 21, 2014, Zach Berta-Thompson modified to allow wrapping with Python
 *   (also, renamed to prevent confusion an honor its developer)
 */

///////////////////////////////////////////////////////////////
/* ZKBT says: change NX and NY for bigger or smaller images! */

//#define NX 32
//#define NY 32
///////////////////////////////////////////////////////////////

 
/* cosmic ray interaction in a CCD simulation code 
 *
 * assume track is straight and infinitely thin, and that the average energy
 * deposition is that of a minimum-ionizing particle, about 80 e's/micron.
 *
 * The equation of the line is:
 *
 * (x - x0)/a = (y - y0)/b = (z - Z0)/c
 *
 * where x0, y0, z0, a, b, and c are given having been obtained from a 
 * routine that calls a random number generator.
 *
 * The CCD is NX by NY pixels with each pixel dx x dy x dz in size.
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "twister.h"
#include "seed_tw_ran.h"

#define DR (M_PI/180.0)
#define NELEC_MICRON 80.0
#define SQRT_FANO_SI 0.34 // square root of the Fano factor for silicon (0.115)


double origin[3], dx[3];   // CCD geometric parameters
int npix[3];   // Set to NX, NY, 1 in initialization


/******************************************************************************/
/* exposure time - should include frame store time (which is approximately the
 * same as the exposure time).
 *
 * exp_time_low_row = exposure time in imaging region (no time in frame store)
 * exp_time_high_row = exposure time in imaging region + max in frame store
 *                     (more or less twice the imaging time)
 *
 * The frame store region has pixels that are the same in size as those in the 
 * primary image pixels.
 */

void usage(char *prog)
{
  fprintf(stderr,"Usage: %s <exp_time_low_row_(s)> <exp_time_high_row_(s)>\n", prog);
  fprintf(stderr,"   <CR_flux N/cm^2/s> > output_file\n");
  exit(-1);
}

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
/* Find intersection point of a line and plane in the case there is one such 
 * point exactly - the number of such points
 * may be, in general, 0, 1, or infinity.
 * This routine only does planes normal to a coordinate axis
 * v0, a, specify the line
 * iplax, plval specify the plane
 */

void line_plane(double v0[3], double a[3], int iplax, double plval,
		int *nint, double pt[3])
{
  int iax;

  if (a[iplax] == 0.0) {
    // The line lies precisely parallel to the plane.
    // The no. of intersection points is zero or infinity.
    if (plval == v0[iplax])
      *nint = 2;   // A two here represents an infinite number of points.
                   /* Ignore these cases for now.
                    * The right thing to do may be to get intersection points
                    * with edges of the rectangle (not to done in this function).
		    */
    else
      *nint = 0;
  }
  else {
    // The line intersects the plane at precisely one point.
    *nint = 1;
    for(iax=0;iax<3;++iax) {
      if (iax == iplax)
	pt[iax] = plval;
      else
	pt[iax] = v0[iax] + ((a[iax]/a[iplax])*(plval - v0[iplax]));
    }
  }
  
}

/******************************************************************************/
void line_box(double v0[3], double a[3], double corn[2][3], int *npts,
	      double pts[2][3])
{

  int i, iax, iax2, iax3, ip, nint[2][3], np;
  double pt[3];

  np = 0;

  for(iax=0;iax<3;++iax) {
    for(ip=0;ip<2;++ip) {

      line_plane(v0,a,iax,corn[ip][iax],&nint[ip][iax],pt);
      iax2 = (iax + 1) % 3;
      iax3 = (iax2 + 1)% 3;
      if (nint[ip][iax] == 1) {

	if ( (pt[iax2] >= corn[0][iax2]) && (pt[iax2] < corn[1][iax2]) &&
	     (pt[iax3] >= corn[0][iax3]) && (pt[iax3] < corn[1][iax3]) ) {
	  pts[np][iax] = pt[iax];
	  pts[np][iax2] = pt[iax2];
	  pts[np][iax3] = pt[iax3];
	  ++np;
	  if (np >= 2)
	    break;
	}

      }

    }
    if (np >= 2)
      break;
  }

  *npts = np;
}

/******************************************************************************/

void points_to_range(double p[2][3], double xr[2], double yr[2])
{
  if (p[0][0] <= p[1][0]) {
    xr[0] = p[0][0];
    xr[1] = p[1][0];
  }
  else {
    xr[0] = p[1][0];
    xr[1] = p[0][0];
  }
  if (p[0][1] <= p[1][1]) {
    yr[0] = p[0][1];
    yr[1] = p[1][1];
  }
  else {
    yr[0] = p[1][1];
    yr[1] = p[0][1];
  }
}

/******************************************************************************/

double path_len(double pts[2][3])
{
  double del, lensq, plen;
  int i;

  lensq = 0.0;
  for(i=0;i<3;++i) {
    del = pts[1][i] - pts[0][i];
    lensq += del*del;
  }

  plen = sqrt(lensq);
  return(plen);
}

/******************************************************************************/
// Find intersection points of ray with CCD surfaces. 
// Use these to do a range in jx or jy (whichever is smaller).
// Find intersection points of the ray with the outer surfaces of the row
// or column of pixels.
// Determine which pixels in the row or column contain charge and how much.

void do_cosmic_ray(double v0[3], double a[3], double *image, long NX, long NY)
{
  double corn[2][3], ptsccd[2][3], ptscol[2][3], ptspix[2][3];
  double xr[2], yr[2], xrc[2], yrc[2], plen, nelec, nelecmn;
  int i, j, iax, npts, ixa, ixb, iya, iyb, jya, jyb;

  for(iax=0;iax<3;++iax) {
    corn[0][iax] = origin[iax];
    corn[1][iax] = origin[iax] + (npix[iax]*dx[iax]);
    // fprintf(stdout,"corners: %f %f\n",corn[0][iax],corn[1][iax]);
  }

  line_box(v0,a,corn,&npts,ptsccd);  // entire CCD
  if (npts != 2)
    return;
  for(iax=0;iax<3;++iax) {
    // fprintf(stdout,"ccd pts: %f %f\n",ptsccd[0][iax],ptsccd[1][iax]);
  }

  points_to_range(ptsccd,xr,yr);

  ixa = floor((xr[0] - origin[0])/dx[0]);
  if (ixa < 0) ixa = 0;
  ixb = ceil((xr[1] - origin[0])/dx[0]);
  if (ixb >= NX) ixb = NX - 1;
  iya = floor((yr[0] - origin[1])/dx[1]);
  if (iya < 0) iya = 0;
  iyb = ceil((yr[1] - origin[1])/dx[1]);
  if (iyb >= NY) iyb = NY - 1;
  // fprintf(stdout,"  %f %f %d %d\n",xr[0],xr[1],ixa,ixb);
  // fprintf(stdout,"  %f %f %d %d\n",yr[0],yr[1],iya,iyb);
  
  // Do each column in the range of intersection
  for(i=ixa;i<=ixb;++i) {
    corn[0][0] = origin[0] + (i*dx[0]);
    corn[1][0] = corn[0][0] + dx[0];
    corn[0][1] = origin[1];
    corn[1][1] = corn[0][1] + (npix[1]*dx[1]);
    for(iax=0;iax<3;++iax) {
      // fprintf(stdout,"%d corners: %f %f\n",i,corn[0][iax],corn[1][iax]);
    }

    line_box(v0,a,corn,&npts,ptscol);  // column i
    if (npts != 2) {
      continue;
    }

    points_to_range(ptscol,xrc,yrc);
    jya = floor((yrc[0] - origin[1])/dx[1]);
    if (jya < 0) jya = 0;
    jyb = ceil((yrc[1] - origin[1])/dx[1]);
    if (jyb >= NY) jyb = NY - 1;

    // Do each pixel in the range of intersection
    for(j=jya;j<=jyb;++j) {
      // fprintf(stdout,"i, j = %d %d\n",i,j);

      corn[0][1] = origin[0] + (j*dx[1]);
      corn[1][1] = corn[0][1] + dx[1];

      line_box(v0,a,corn,&npts,ptspix);  // pixel i, j
      if (npts != 2) {
	continue;
      }

      plen = path_len(ptspix);
      nelecmn = NELEC_MICRON*plen;
      nelec = 0.0;
      if (nelecmn > 0.0)
	nelec = nelecmn + gasdev()*SQRT_FANO_SI*sqrt(nelecmn);
      *(image + i*NY + j) += nelec;
    }

  }

}

/******************************************************************************/

void print_image(double *image, long NX, long NY)
{
  int i, j, ne;

  for(j=(NY-1);j>=0;--j) {
    for(i=0;i<NX;++i) {
      ne = *(image + i*NY + j);
      fprintf(stdout," %4d",ne);
    }
    fprintf(stdout,"\n");
  }

}

/******************************************************************************/
/* expos_fac = factor (>= 1.0) that gives ratio of effective exposure times to
 *             cosmic rays from the top of the imaging array to the bottom,
 *             due to residence time in the frame store region
 */

void get_ran_cr(double origin[3], double ccddim[3], double v0[3], double a[3],
		double ratefac)
{
  double cth, sth, ph, x, y, c1, c2, xc;
  int i;

  c1 = 2.0/(1.0 + ratefac);
  c2 = 2.0*(1.0 - c1);

  for(i=0;i<3;++i) {
    x = ranMT();
    if (i != 1)
      y = x;
    else {
      if ( (x >= 0.0) && ( x < 1.0) ) {
	xc = x/c1;
	y = 2.0*xc/(1.0 + sqrt(1.0 + (2.0*c2*xc/c1)));
      }
      else {
	fprintf(stderr,"get_ran_cr(): error x = %lf\n",x);
	y = x;
      }
    }
    v0[i] = origin[i] + (y*ccddim[i]);
  }

  ph = 2.0*M_PI*ranMT();
  cth = -1.0 + 2.0*ranMT();
  if (fabs(cth) <= 1.0)
    sth = sqrt(1.0 - cth*cth);
  else {
    cth = 1.0;
    sth = 0.0;
  }

  a[0] = sth*cos(ph);
  a[1] = sth*sin(ph);
  a[2] = cth;

}


/******************************************************************************/

void print_cr(FILE *fp, double v0[3], double a[3])
{
  int i;

  for(i=0;i<3;++i) {
    fprintf(fp,"  %f",v0[i]);
  }
  fprintf(fp,"\n");
  for(i=0;i<3;++i) {
    fprintf(fp,"  %f",a[i]);
  }
  fprintf(fp,"\n\n");

}

/******************************************************************************/
// ZKBT says: cosmic generation moved from main() to a separate function,
//    which can be called directly from Python and returns a 1D image array
double * cosmical(double exptm1, double exptm2, long ncr, long NX, long NY)
{
  clock_t t;
  t = clock();
  double v0[3], a[3], ccddim[3], th, ph, exptm, rfac;
  //double image[NX*NY];
  double* image = malloc(NX*NY*sizeof(double));
  int i, j, icr;
  unsigned long seed;
  
  rfac = exptm2/exptm1;
  exptm = 0.5*(exptm1 + exptm2);

  for(i=0;i<NX;++i) {
    for(j=0;j<NY;++j) {
      image[i*NY + j] = 0.0;
    }
  }

  origin[0] = 0.0;
  origin[1] = 0.0;
  origin[2] = 0.0;

  dx[0] = 15.0;  // microns
  dx[1] = 15.0;
  dx[2] = 100.0;

  ccddim[0] = NX*dx[0];
  ccddim[1] = NY*dx[1];
  ccddim[2] = dx[2];

  npix[0] = NX;
  npix[1] = NY;
  npix[2] = 1;

  // Invent a CR for testing

  v0[0] = 107.0;
  v0[1] = 99.9;
  v0[2] = 20.2;
  th = 65.3*DR;
  ph = 31.4*DR;
  a[0] = sin(th)*cos(ph);
  a[1] = sin(th)*sin(ph);
  a[2] = cos(th);

  seed = get_tw_seed();
  fprintf(stderr," using CosmicAl to generate simulated cosmic rays for a %ldx%ld pixel image, assuming:\n",NX,NY);
  fprintf(stderr,"     exposure times spanning %f to %f s\n", exptm1, exptm2);
  fprintf(stderr,"     a total of %ld events \n", ncr);
  fprintf(stderr,"     deposition of %f electrons/micron\n", NELEC_MICRON);
  fprintf(stderr,"     Fano factor of %f\n", SQRT_FANO_SI*SQRT_FANO_SI);
  fprintf(stderr,"     and a seed of %lu \n",seed);
  seedMT(seed);

  //fprintf(stderr,"ncrmn, ncr = %d, %d\n\n",ncrmn,ncr);


  for (icr=0;icr<ncr;++icr) {
    get_ran_cr(origin,ccddim,v0,a,rfac);
    //fprintf(stderr," cosmic ray #%d\n",icr+1);
    //print_cr(stderr,v0,a);
    do_cosmic_ray(v0,a,image,NX,NY);
  }
  t = clock() - t;
  fprintf(stderr," added %ld cosmic rays in %f seconds \n",ncr,((float)t)/CLOCKS_PER_SEC);
  //print_image(image,NX,NY);
  return image;
}


int main(int argc, char *argv[])
{

  double exptm1, exptm2;
  long ncr;
  if (argc != 4)
    usage(argv[0]);
  exptm1 = atof(argv[1]);
  exptm2 = atof(argv[2]);
  ncr = atol(argv[3]);

  long NX, NY;
  NX = 16;
  NY = 16;
  double *test = cosmical(exptm1,exptm2,ncr,NX,NY);
  //for(int i=0;i<NX*NX;i++)
  //	printf("*(test + %d) = %f ", i, *(test + i));
  

}

/******************************************************************************/
