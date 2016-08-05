/* seed_tw_ran.c
 * Alan M. Levine
 * July 8, 2014
 * Heritage: ran_mt.c
 *
 * Use time function to obtain a seed for twister.c 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

unsigned long get_tw_seed()
{
  time_t tr;

  tr = clock();
  if ( (tr % 2) == 0)
    ++tr;

  //fprintf(stderr,"time tr = %lu \n",tr);

  return( (unsigned long) tr);
}
