/* cosmical[full]|[small].c
 */
/* cr_img5_z1.c
 * Alan M. Levine
 * December 17, 2014
 * Heritage: cr_img4_z1.c
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
 *
 * Oct. 31, 2014 - Begin putting in code to simulate energy straggling
 * more accurately.  The energy-loss code was tested in "test_landau2.c".
 *
 * August 21, 2014 - Zach Berta-Thompson modified to allow wrapping with Python.
 *
 * Dec. 15, 2014 - Merge cr_img4.c and cosmical.c (latter is from cr_img3.c).
 *
 * Dec. 17, 2014 - Add code to optionally convolve image with a 3x3 kernel
 *                 that grossly represents the effect of charge diffusion.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "twister.h"
#include "seed_tw_ran.h"
#include "fmemopen.h"
#include "st_dat.h"

#define DR (M_PI/180.0)
#define NELEC_MICRON 80.0
#define SQRT_FANO_SI 0.34 // square root of the Fano factor for silicon (0.115)

#define NLMAX 6000
#define M_e_KEV 511.0
#define N_AVO (6.02e23)  // Avogadro's number
#define ELEC_RADIUS (2.8e-13) // classical electronic radius in cm
#define ELEC_RADIUS_SQ (7.8e-26) // square of classical electron radius in cm^2
#define RHO_SI 2.32  // density of silicon
#define N_ION 3.68  // eV per electron-hole pair
#define M_PROTON_MEV 938.0

double lamb[NLMAX], phi[NLMAX], phicum[NLMAX];
int nlamb;
double eprot, pgamma, bfac, beta;
double xlen, xi, xlam, delta, capi, lnepp;

double origin[3], dx[3], ccddim[3];   // CCD geometric parameters
// The 'difker' array contains the values for rssq = 5.0 microns.
double difker[3][3] = {
  {0.0034, 0.0516, 0.0034},
  {0.0516, 0.7798, 0.0516},
  {0.0034, 0.0516, 0.0034}
};
int npix[3];   // Set to NX, NY, 1 in initialization
int ncr;
double ncrmn, rfac;
double *image, *imaged, *temp;

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
  fprintf(stderr,"Usage: %s <exp_time_low_row_(s)> <exp_time_high_row_(s)>\n",prog);
  fprintf(stderr,"   <CR_flux N/cm^2/s> <no_of_images_to_sim> <idodif=0 or 1>\n");
  fprintf(stderr,"    > output_file\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"     'idodif' = 0  => do NOT do diffusion\n");
  fprintf(stderr,"              = 1  => do diffusion\n");
  fprintf(stderr,"\n");
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
// For reading the file st.dat (in tess2/cr3 on Oct. 30, 2014)

void read_phi_lambda(FILE *fp)
{
  int i, j;
  double pl0;

  for(i=0;i<NLMAX;++i) {
    if (fscanf(fp,"%d %lf %lf %lf %lf",&j,&lamb[i],&pl0,&phi[i],&phicum[i]) != 5) {
      break;
    }
  }
  nlamb = i;
}

/******************************************************************************/

double get_lambda_random()
{
  double x, xlam, cumdif, lamdif;
  int i;

  x = ranMT();
  xlam = 0.0;
  for(i=0;i<(nlamb-1);++i) {
    if ( (x >= phicum[i]) && (x <= phicum[i+1]) ) {
      cumdif = phicum[i+1] - phicum[i];
      lamdif = lamb[i+1] - lamb[i];
      if (cumdif > 0.0) {
        xlam = lamb[i] + (lamdif*(x - phicum[i])/cumdif);
      }
      else {
        xlam = lamb[i];
      }
      break;
    }
  }
  // fprintf(stderr,"x,i,xlam = %f %d %f\n",x,i,xlam);
  return(xlam);
}

/******************************************************************************/
//  0.0 < beta < 1.0
//  capi \sim Z*13.5 eV (= 189 eV for silicon)

double get_ln_epsprime(double beta, double capi)
{
  double lnepp, bsq;

  bsq = beta*beta;
  lnepp = log( (1.0-bsq)*capi*capi/(2.0*M_e_KEV*bsq) ) + bsq;

  return(lnepp);
}

/******************************************************************************/
// xlen = path length in microns
// This routine assumes the particle is traveling in pure silicon.
// xi is in keV

double get_xi(double xlen, double beta)
{
  double xi, xcm;

  xcm = xlen/10000.0;
  xi = 2.0*M_PI*M_e_KEV*xcm*ELEC_RADIUS_SQ*(RHO_SI*N_AVO*0.5)/(beta*beta);

  return(xi);
}

/******************************************************************************/
// return value is in kev

double get_energy_loss(double beta, double xi, double xlam, double lnepp)
{
  double delta;

  delta = xi*(xlam + log(xi) - lnepp + 0.423);
  return(delta);
}

/******************************************************************************/
// delta is the energy loss in keV

double get_num_elect(double delta)
{
  double nelec;

  nelec = 1000.0*delta/N_ION;
  return(nelec);
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

  int iax, iax2, iax3, ip, nint[2][3], np;
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

void do_cosmic_ray(double v0[3], double a[3], long NX, long NY)
{
  double corn[2][3], ptsccd[2][3], ptscol[2][3], ptspix[2][3];
  double xr[2], yr[2], xrc[2], yrc[2], plen, nelec;
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

      /*  Inaccurate old ionization simulation
          nelecmn = NELEC_MICRON*plen;
          nelec = 0.0;
          if (nelecmn > 0.0)
          nelec = nelecmn + gasdev()*SQRT_FANO_SI*sqrt(nelecmn);
      */

      // New ionization simulation ala Landau, Bichsel, et al.
      // Neglect fluctuations in conversion of energy loss (delta)
      // to no. of electron-hole pairs (nelec).
      xlen = plen;
      xi = get_xi(xlen,beta);
      xlam = get_lambda_random();
      delta = get_energy_loss(beta,xi,xlam,lnepp);
      nelec = get_num_elect(delta);
      /* fprintf(stdout,"xlen,xi,xlam,delta,nelec = %f %f %f %f %f\n\n",
         xlen,xi,xlam,delta,nelec); */

      *(image + i*NY + j) += nelec;
    }
  }
}

/******************************************************************************/

void print_image(long NX, long NY, FILE *fp)
{
  int i, j, ne;

  for(j=(NY-1);j>=0;--j) {
    for(i=0;i<NX;++i) {
      ne = *(image + i*NY + j);
      fprintf(fp," %4d",ne);
    }
    fprintf(fp,"\n");
  }

}

/******************************************************************************/
/* expos_fac = factor (>= 1.0) that gives ratio of effective exposure times to
 *             cosmic rays from the top of the imaging array to the bottom,
 *             due to residence time in the frame store region
 */

void get_ran_cr(double origin[3], double v0[3], double a[3], double ratefac)
{
  double cth, sth, ph, x, y, c1, c2, xc;
  int i;
  //fprintf(stdout," %lf\n",ratefac);

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
// Convolve the 2-d image array w/ the 3x3 kernel difker.

void do_diffusion(long NX, long NY)
{
  int i, j, il, im, jl, jm, ip, jp;
  double ximg;

  // Swap 'image' and 'imaged' so the final results are in 'image'.

  temp = image;
  image = imaged;
  imaged = temp;

  for(j=0;j<NY;++j) {
    jl = j - 1;
    if (jl < 0 ) jl = 0;
    jm = j + 1;
    if (jm == NY) jm = NY - 1;

    for(i=0;i<NX;++i) {
      il = i - 1;
      if (il < 0 ) il = 0;
      im = i + 1;
      if (im == NX) im = NX - 1;

      for(jp=jl;jp<=jm;++jp) {
          for(ip=il;ip<=im;++ip) {
              ximg = imaged[ip*NY + jp];
              image[i*NY + j] += difker[ip-i+1][jp-j+1]*ximg;
          }
      }
    }
  }
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
// Tasks that we don't want to do unnecessarily when doing multiple images

void cosmical_setup(double crfl, double exptm1, double exptm2, long NX, long NY, int idodif)
{
  double exptm, tot;
  int ip, jp;
  unsigned long seed;
  FILE *fpin = fmemopen(phi_lambda_data, strlen(phi_lambda_data), "r");
  read_phi_lambda(fpin);
  fclose(fpin);

  eprot = 8000.0;  // MeV
  pgamma = eprot/M_PROTON_MEV;
  bfac = 1.0 - (1.0/(pgamma*pgamma));
  if (bfac >= 0.0)
    beta = sqrt(bfac);
  else
    beta = 0.0;
  capi = 13.5*14.0/1000.0;
  lnepp = get_ln_epsprime(beta,capi);

  dx[0] = 15.0;  // microns
  dx[1] = 15.0;
  dx[2] = 100.0;

  rfac = exptm2/exptm1;
  exptm = 0.5*(exptm1 + exptm2);
  ncrmn = exptm*crfl*(NX*dx[0]/10000.0)*(NY*dx[1]/10000.0);
  // fprintf(stderr,"setup: ncrmn = %f\n",ncrmn);

  image = malloc(NX*NY*sizeof(double));
  if (idodif == 1)
     imaged = malloc(NX*NY*sizeof(double));

  origin[0] = 0.0;
  origin[1] = 0.0;
  origin[2] = 0.0;

  ccddim[0] = NX*dx[0];
  ccddim[1] = NY*dx[1];
  ccddim[2] = dx[2];

  npix[0] = NX;
  npix[1] = NY;
  npix[2] = 1;

  if(idodif == 1) {
    tot = 0.0;
    for(jp=0;jp<3;++jp) {
      for(ip=0;ip<3;++ip) {
    tot += difker[ip][jp];
      }
    }
    for(jp=0;jp<3;++jp) {
      for(ip=0;ip<3;++ip) {
    difker[ip][jp] /= tot;
    fprintf(stderr,"difker[%d][%d] = %f\n",ip,jp,difker[ip][jp]);
      }
    }
  }

  // Invent a CR for testing
  /*
  v0[0] = 107.0;
  v0[1] = 99.9;
  v0[2] = 20.2;
  th = 65.3*DR;
  ph = 31.4*DR;
  a[0] = sin(th)*cos(ph);
  a[1] = sin(th)*sin(ph);
  a[2] = cos(th);
  */

  seed = get_tw_seed();
  //fprintf(stderr," using CosmicAl to generate simulated cosmic rays for a %ldx%ld pixel image, assuming:\n", NX,NY);
  //fprintf(stderr,"     exposure times spanning %f to %f s\n", exptm1, exptm2);
  //fprintf(stderr,"     a mean of %f events \n", ncrmn);
  //fprintf(stderr,"     Fano factor of %f\n", SQRT_FANO_SI*SQRT_FANO_SI);
  //fprintf(stderr,"     and a seed of %lu \n",seed);
  seedMT(seed);
}


/******************************************************************************/
// ZKBT says: cosmic generation moved from main() to a separate function,
//    which can be called directly from Python and returns a 1D image array

double * cosmical(double crfl, double exptm1, double exptm2, long NX, long NY, int idodif)
{
  cosmical_setup(crfl, exptm1, exptm2, NX, NY, idodif);

  double v0[3], a[3], rfac;
  int i, j, icr;
  rfac = exptm2/exptm1;

  for(i=0;i<NX;++i) {
    for(j=0;j<NY;++j) {
      image[i*NY + j] = 0.0;
      if (idodif == 1)
        imaged[i*NY + j] = 0.0;
      }
    }

    // ncr = ncrmn + gasdev()*sqrt(ncrmn);
    ncr = poidev(ncrmn);
    // ncr = 3;  // testing only
    // fprintf(stderr,"ncrmn, ncr = %f, %d\n\n",ncrmn,ncr);

    for (icr=0;icr<ncr;++icr) {
      get_ran_cr(origin,v0,a,rfac);
      //fprintf(stderr," cosmic ray #%d\n",icr+1);
      //print_cr(stderr,v0,a);
      do_cosmic_ray(v0,a,NX,NY);
      //print_image(NX,NY,stderr);
    }

    if (idodif == 1)
      do_diffusion(NX,NY);
    free(imaged);
    free(temp);
    return image;
  }


/******************************************************************************/

int main(int argc, char *argv[])
{
  char fileout[64];
  double exptm1, exptm2, crfl;
  int iimg, nimg, idodif;
  long NX, NY;
  //FILE *fp;
  clock_t t;

  t = clock();

  NX = 32;
  NY = 256;

  if (argc != 6)
    usage(argv[0]);
  exptm1 = atof(argv[1]);
  exptm2 = atof(argv[2]);
  crfl = atof(argv[3]);
  nimg = atoi(argv[4]);
  idodif = atoi(argv[5]);  // = 0 or 1 => no or yes, to diffusion modelling

  cosmical_setup(crfl,exptm1,exptm2,NX,NY, idodif);

  // loop here to end to do multiple images (nimg) with one setup
  for(iimg=0;iimg<nimg;++iimg) {
    cosmical(crfl,exptm1,exptm2,NX,NY,idodif);

    /* for(i=0;i<NX*NX;i++)
       printf("*(image + %d) = %f ", i, *(image + i));
    */
    sprintf(fileout,"image_%05d.dat",iimg);
    // 'fp" could be stdout instead of a pointer to a specific file
    //fp = fopen(fileout,"w");
    //print_image(NX,NY,fp);
    //fclose(fp);
    t = clock() - t;
  }
}

/******************************************************************************/
