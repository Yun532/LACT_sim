/* ============================================================================

Copyright (C) 2021  Konrad Bernloehr

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

============================================================================ */

/** @file select_iact.c
 *  @short A program reading simulated CORSIKA data written through the
 *      IACT interface and, if it contains extra information on particles
 *      emitting Cherenkov light, reduce to light from selected particles.
 *
 *  Relevant environment variables:
 *     PRINT_TEL_VERBOSE
 *     MAX_PRINT_ARRAY
 *
 *  @author Konrad Bernloehr
 *
 *  @date   2021
 */

/** @defgroup select_iact_c The select_iact program */
/** @{ */

#include "initial.h"      /* This file includes others as required. */
#include "io_basic.h"     /* This file includes others as required. */
#include "mc_tel.h"
#include "fileopen.h"
// #include "eventio_version.h"
#include <sys/time.h>
#include <strings.h>

int array_select_mc_phot (IO_BUFFER *iobuf);
int tel_select_mc_phot (IO_BUFFER *iobuf);
int tel_select_mc_phot3d (IO_BUFFER *iobuf);
void add_selector (double m1, double m2, double E1, double E2, int c);
int select_bunches(struct bunch *bunches, int *nbunches, double *photons);
int select_bunches3d(struct bunch3d *bunches, int *nbunches, double *photons);
void ioerrorcheck (void);

#ifndef MAXTEL
# ifdef MAXIMUM_TELESCOPES
#  define MAXTEL MAXIMUM_TELESCOPES
# else
#  define MAXTEL 5
# endif
#endif

struct bunch *tel_bunches[MAXTEL];
struct bunch3d *tel_bunches3d[MAXTEL];
int max_bunches[MAXTEL];
int max_bunches3d[MAXTEL];
int tel_nbunches[MAXTEL];
int tel_nbunches3d[MAXTEL];
double tel_photons[MAXTEL];
double tel_photons3d[MAXTEL];

struct selector
{
   double min_mass, max_mass;
   double min_energy, max_energy;
   int charge;
};

typedef struct selector Selector;

Selector *selectors = NULL;
size_t nselect = 0;
static int verbose = 0;

void add_selector (double m1, double m2, double e1, double e2, int c)
{
   if ( selectors == NULL )
   {
      selectors = (Selector *) malloc(sizeof(Selector));
      if ( selectors == NULL )
      {
         fprintf(stderr,"Allocation failed.\n");
         exit(1);
      }
   }
   else
   {
      Selector *new_selectors = (Selector *) realloc(selectors,sizeof(Selector)*(nselect+1));
      if ( new_selectors == NULL )
      {
         fprintf(stderr,"Allocation failed.\n");
         exit(1);
      }
      selectors = new_selectors;
   }

   selectors[nselect].min_mass = m1;
   selectors[nselect].max_mass = m2;
   selectors[nselect].min_energy = e1;
   selectors[nselect].max_energy = e2;
   selectors[nselect].charge = c;

   nselect++;
}

struct bunch *sel_bunch = NULL;
struct bunch3d *sel_bunch3d = NULL;
int sel_max = 0;
int sel_max3d = 0;

int select_bunches(struct bunch *b, int *nb, double *photons)
{
   int ib, nsel = 0;
   double phot_sel = 0.;
   
   if ( b == NULL || nb == NULL || photons == NULL || selectors == NULL )
      return -1;
   if ( (*nb) <= 0 )
      return 0;

   if ( sel_max < (*nb) || sel_bunch == NULL )
   {
      struct bunch *sbn = NULL;
      if ( sel_bunch == NULL )
         sbn = (struct bunch *) malloc(sizeof(struct bunch)*(*nb));
      else
         sbn = (struct bunch *) realloc(sel_bunch,sizeof(struct bunch)*(*nb));
      if ( sbn == NULL )
      {
         fprintf(stderr,"Bunches allocation error.\n");
         exit(1);
      }
      sel_bunch = sbn;
      sel_max = (*nb);
   }

   for ( ib=0; ib+1<(*nb); ib++ )
   {
      int selected = 0;
      if ( b[ib].lambda < 9000. && b[ib+1].lambda >= 9000. ) /* Normal photon bunch plus emitting particle ? */
      {
         size_t is;
         for ( is=0; is<nselect; is++ )
         {
            int sm = 0, se = 0, sc = 0;
            /* Either the mass of the emitting particle is in the requested range or we don't care. */
            double mass = b[ib+1].cx; /* [GeV/c^2] */
            if ( (mass >= selectors[is].min_mass && mass <= selectors[is].max_mass) ||
                 (selectors[is].min_mass == 0. && selectors[is].max_mass == 0.) )
               sm = 1;
            /* Either the charge matches or we don't care */
            int charge = Nint(b[ib+1].cy);
            if ( charge == selectors[is].charge || selectors[is].charge == 0 )
               sc = 1;
            /* Either the energy or the emitting particle is in the requested range or we don't care. */
            double energy = b[ib+1].photons; /* [GeV] */
            if ( (energy >= selectors[is].min_energy && energy <= selectors[is].max_energy) ||
                 (selectors[is].min_energy == 0. && selectors[is].max_energy == 0.) )
               se = 1;
            if ( sm && se && sc ) /* All criteria passed */
            {
#ifdef DEBUG_SELECT
printf("Bunch at %f m, %f m in direction %f,%f, arrival time %f ns, "
       "emission at %f m height, with %f photons of wavelength %f nm.\n",
       b[ib].x*0.01, b[ib].y*0.01, b[ib].cx, b[ib].cy, b[ib].ctime, b[ib].zem*0.01, b[ib].photons, b[ib].lambda);
printf("emitted by particle of mass %5.3f MeV/c^2, charge %1.0f, energy %7.5f GeV at time %f ns\n",
       b[ib+1].cx*1000., b[ib+1].cy, b[ib+1].photons, b[ib+1].zem);
printf("selected by selector for mass range %f to %f GeV, energy range %f to %f GeV, charge %d.\n",
       selectors[is].min_mass, selectors[is].max_mass, 
       selectors[is].min_energy, selectors[is].max_energy,
       selectors[is].charge);
#endif
               /* Raw copy of photon bunch plus the extra data on the emitting particle */
               memcpy(sel_bunch+nsel,b+ib,sizeof(struct bunch)*2);
               selected = 1;
               break;
            }
         }
         if ( selected )
         {
            phot_sel += b[ib].photons;
            nsel += 2; /* We copied photon bunch and emitting particle info */
         }
         ib++;
      }
   }

   if ( verbose )
    printf("Remaining: %d of %d photon bunches with %4.2f of %4.2f photons.\n",
      nsel, *nb, phot_sel, *photons);
   
   int nb_orig = (*nb);
   *nb = nsel;
   *photons = phot_sel;
   if ( nsel != 0 && nsel != nb_orig ) /* copy needed */
   {
      memcpy(b,sel_bunch,sizeof(struct bunch)*nsel);
      if ( verbose )
         printf("Original photon bunch list replaced with list of %d selected bunches.\n", nsel);
   }

   return nsel;
}

int select_bunches3d(struct bunch3d *b, int *nb, double *photons)
{
   int ib, nsel = 0;
   double phot_sel = 0.;
   
   if ( b == NULL || nb == NULL || photons == NULL || selectors == NULL )
      return -1;
   if ( (*nb) <= 0 )
      return 0;

   if ( sel_max < (*nb) || sel_bunch == NULL )
   {
      struct bunch3d *sbn = NULL;
      if ( sel_bunch3d == NULL )
         sbn = (struct bunch3d *) malloc(sizeof(struct bunch3d)*(*nb));
      else
         sbn = (struct bunch3d *) realloc(sel_bunch,sizeof(struct bunch3d)*(*nb));
      if ( sbn == NULL )
      {
         fprintf(stderr,"Bunches allocation error.\n");
         exit(1);
      }
      sel_bunch3d = sbn;
      sel_max = (*nb);
   }

   for ( ib=0; ib+1<(*nb); ib++ )
   {
      int selected = 0;
      if ( b[ib].lambda < 9000. && b[ib+1].lambda >= 9000. ) /* Normal photon bunch plus emitting particle ? */
      {
         size_t is;
         for ( is=0; is<nselect; is++ )
         {
            int sm = 0, se = 0, sc = 0;
            /* Either the mass of the emitting particle is in the requested range or we don't care. */
            double mass = b[ib+1].cx; /* [GeV/c^2] */
            if ( (mass >= selectors[is].min_mass && mass <= selectors[is].max_mass) ||
                 (selectors[is].min_mass == 0. && selectors[is].max_mass == 0.) )
               sm = 1;
            /* Either the charge matches or we don't care */
            int charge = Nint(b[ib+1].cy);
            if ( charge == selectors[is].charge || selectors[is].charge == 0 )
               sc = 1;
            /* Either the energy or the emitting particle is in the requested range or we don't care. */
            double energy = b[ib+1].photons; /* [GeV] */
            if ( (energy >= selectors[is].min_energy && energy <= selectors[is].max_energy) ||
                 (selectors[is].min_energy == 0. && selectors[is].max_energy == 0.) )
               se = 1;
            if ( sm && se && sc ) /* All criteria passed */
            {
#ifdef DEBUG_SELECT
printf("3D bunch at %f m, %f m, %f m in direction %f,%f,%f, arrival time %f ns, "
       "emission at %f m distance, with %f photons of wavelength %f nm.\n",
       b[ib].x*0.01, b[ib].y*0.01, b[ib].z*0.01, b[ib].cx, b[ib].cy, b[ib].cz, 
       b[ib].ctime, b[ib].dist*0.01, b[ib].photons, b[ib].lambda);
printf("emitted by particle of mass %5.3f MeV/c^2, charge %1.0f, energy %7.5f GeV at time %f ns\n",
       b[ib+1].cx*1000., b[ib+1].cy, b[ib+1].photons, b[ib+1].dist);
printf("selected by selector for mass range %f to %f GeV, energy range %f to %f GeV, charge %d.\n",
       selectors[is].min_mass, selectors[is].max_mass, 
       selectors[is].min_energy, selectors[is].max_energy,
       selectors[is].charge);
#endif
               /* Raw copy of photon bunch plus the extra data on the emitting particle */
               memcpy(sel_bunch3d+nsel,b+ib,sizeof(struct bunch3d)*2);
               selected = 1;
               break;
            }
         }
         if ( selected )
         {
            phot_sel += b[ib].photons;
            nsel += 2; /* We copied photon bunch and emitting particle info */
         }
         ib++;
      }
   }

   if ( verbose )
    printf("Remaining: %d of %d 3D photon bunches with %4.2f of %4.2f photons.\n",
      nsel, *nb, phot_sel, *photons);
   
   int nb_orig = (*nb);
   *nb = nsel;
   *photons = phot_sel;
   if ( nsel != 0 && nsel != nb_orig ) /* copy needed */
   {
      memcpy(b,sel_bunch,sizeof(struct bunch3d)*nsel);
      if ( verbose )
         printf("Original 3D photon bunch list replaced with list of %d selected 3D bunches.\n", nsel);
   }

   return nsel;
}

int tel_select_mc_phot (IO_BUFFER *iobuf)
{
   int rc;
   int nbunches=0;
   double photons=0.;
   int itel=0, iarray=0;

   rc = read_tel_photons(iobuf,0,&iarray,&itel,&photons,NULL,&nbunches);
   if ( rc != -10 )
   {
      fprintf(stderr,"\nNot proper MC photons data block.\n");
      skip_subitem(iobuf);
      return -3;
   }
   if ( iarray == 999 && itel == 999 )
   {
      printf("Got an unexpected particles data block.\n");
      skip_subitem(iobuf);
      return 0;
   }
   if ( verbose )
    printf("Got %d bunches with %f photons for telescope %d in array %d.\n",
      nbunches,photons,itel,iarray);
   if ( itel < 0 || itel >= MAXTEL )
   {
      printf("Telescope %d is outside valid range.\n", itel);
      skip_subitem(iobuf);
      return -3;
   }
   if ( tel_bunches[itel] == NULL )
   {
      tel_bunches[itel] = (struct bunch *) malloc(sizeof(struct bunch)*nbunches);
      max_bunches[itel] = nbunches;
   }
   else if ( nbunches > max_bunches[itel] )
   {
      tel_bunches[itel] = (struct bunch *) realloc(tel_bunches[itel],sizeof(struct bunch)*nbunches);
      max_bunches[itel] = nbunches;
   }
   if ( tel_bunches[itel] == NULL )
   {
      fflush(stdout);
      fprintf(stderr,"(Re-) Allocation of memory for photon bunches failed.\n");
      exit(1);
   }

   rc = read_tel_photons(iobuf,nbunches,&iarray,&itel,&tel_photons[itel],
      tel_bunches[itel],&tel_nbunches[itel]);

   if ( rc != 0 )
   {
      if ( itel >= 0 && itel < MAXTEL )
      {
         tel_photons[itel] = 0.;
         tel_nbunches[itel] = 0;
      }
      return rc;
   }

   select_bunches(tel_bunches[itel],&tel_nbunches[itel],&tel_photons[itel]);

   if ( tel_nbunches[itel] > 0 )
   {
      if ( verbose )
         printf("Writing photon bunch list of %d bunches for telescope no. %d in array %d.\n", 
            tel_nbunches[itel], itel, iarray);
      write_tel_photons(iobuf, iarray, itel, tel_photons[itel], tel_bunches[itel],
         tel_nbunches[itel], 0, NULL);
      tel_nbunches[itel] = 0;
      tel_photons[itel] = 0.;
   }

   return 0;
}

int tel_select_mc_phot3d (IO_BUFFER *iobuf)
{
   int rc;
   int nbunches=0;
   double photons=0.;
   int itel=0, iarray=0;

   rc = read_tel_photons3d(iobuf,0,&iarray,&itel,&photons,NULL,&nbunches);
   if ( rc != -10 )
   {
      fprintf(stderr,"\nNot proper MC 3D photons data block.\n");
      skip_subitem(iobuf);
      return -3;
   }
   if ( iarray == 999 && itel == 999 )
   {
      printf("Got an unexpected particles data block.\n");
      skip_subitem(iobuf);
      return 0;
   }
   if ( verbose )
    printf("Got %d 3D bunches with %f photons for telescope %d in array %d.\n",
      nbunches,photons,itel,iarray);
   if ( itel < 0 || itel >= MAXTEL )
   {
      printf("Telescope %d is outside valid range.\n", itel);
      skip_subitem(iobuf);
      return -3;
   }
   if ( tel_bunches3d[itel] == NULL )
   {
      tel_bunches3d[itel] = (struct bunch3d *) malloc(sizeof(struct bunch3d)*nbunches);
      max_bunches3d[itel] = nbunches;
   }
   else if ( nbunches > max_bunches3d[itel] )
   {
      tel_bunches3d[itel] = (struct bunch3d *) realloc(tel_bunches3d[itel],sizeof(struct bunch3d)*nbunches);
      max_bunches3d[itel] = nbunches;
   }
   if ( tel_bunches3d[itel] == NULL )
   {
      fflush(stdout);
      fprintf(stderr,"(Re-) Allocation of memory for photon 3D bunches failed.\n");
      exit(1);
   }

   rc = read_tel_photons3d(iobuf,nbunches,&iarray,&itel,&tel_photons3d[itel],
      tel_bunches3d[itel],&tel_nbunches3d[itel]);

   if ( rc != 0 )
   {
      if ( itel >= 0 && itel < MAXTEL )
      {
         tel_photons3d[itel] = 0.;
         tel_nbunches3d[itel] = 0;
      }
      return rc;
   }

   select_bunches3d(tel_bunches3d[itel],&tel_nbunches3d[itel],&tel_photons3d[itel]);

   if ( tel_nbunches3d[itel] > 0 )
   {
      if ( verbose )
         printf("Writing 3D photon bunch list of %d bunches for telescope no. %d in array %d.\n", 
            tel_nbunches3d[itel], itel, iarray);
      write_tel_photons3d(iobuf, iarray, itel, tel_photons3d[itel], tel_bunches3d[itel],
         tel_nbunches3d[itel], 0, NULL);
      tel_nbunches3d[itel] = 0;
      tel_photons3d[itel] = 0.;
   }

   return 0;
}

/* --------------------- array_select_mc_phot ---------------------- */
/**
 *  Select Monte Carlo photons.
*/  
 
int array_select_mc_phot (IO_BUFFER *iobuf)
{
   int type;
   int rc;
   IO_ITEM_HEADER item_header;
   int nbunches=0;
   double photons=0.;
   int itel=0, iarray=0;
   
   if ( (rc = begin_read_tel_array(iobuf, &item_header, &iarray)) < 0 )
      return rc;
   if ( verbose )
      printf("\nMC photon or photo-electron data for array %d\n", iarray);

   while ( (type = next_subitem_type(iobuf)) > 0 )
   {
      nbunches=0;
      photons=0.;
      
      switch (type)
      {
         case IO_TYPE_MC_PHOTONS:
            rc = read_tel_photons(iobuf,0,&iarray,&itel,&photons,NULL,&nbunches);
            if ( rc != -10 )
            {
               fprintf(stderr,"\nNot a proper MC photons data block.\n");
               skip_subitem(iobuf);
               continue;
            }
            if ( verbose )
             printf("Got %d bunches with %f photons for telescope %d in array %d.\n",
               nbunches,photons,itel,iarray);
            if ( itel < 0 || itel >= MAXTEL )
            {
               printf("Telescope %d is outside valid range.\n", itel);
               skip_subitem(iobuf);
               continue;
            }
            if ( tel_bunches[itel] == NULL )
            {
               tel_bunches[itel] = (struct bunch *) malloc(sizeof(struct bunch)*nbunches);
               max_bunches[itel] = nbunches;
            }
            else if ( nbunches > max_bunches[itel] )
            {
               tel_bunches[itel] = (struct bunch *) realloc(tel_bunches[itel],sizeof(struct bunch)*nbunches);
               max_bunches[itel] = nbunches;
            }
            if ( tel_bunches[itel] == NULL )
            {
               fflush(stdout);
               fprintf(stderr,"(Re-) Allocation of memory for photon bunches failed.\n");
               exit(1);
            }
            rc = read_tel_photons(iobuf,nbunches,&iarray,&itel,&tel_photons[itel],
               tel_bunches[itel],&tel_nbunches[itel]);
            if ( rc < 0 )
            {
               if ( itel >= 0 && itel < MAXTEL )
               {
                  tel_photons[itel] = 0.;
                  tel_nbunches[itel] = 0;
               }
               get_item_end(iobuf,&item_header);
               return rc;
            }
            select_bunches(tel_bunches[itel],&tel_nbunches[itel],&tel_photons[itel]);
            break;
         case IO_TYPE_MC_PHOTONS3D:
            rc = read_tel_photons3d(iobuf,0,&iarray,&itel,&photons,NULL,&nbunches);
            if ( rc != -10 )
            {
               fprintf(stderr,"\nNot a proper MC 3D photons data block.\n");
               skip_subitem(iobuf);
               continue;
            }
            if ( verbose )
             printf("Got %d 3D bunches with %f photons for telescope %d in array %d.\n",
               nbunches,photons,itel,iarray);
            if ( itel < 0 || itel >= MAXTEL )
            {
               printf("Telescope %d is outside valid range.\n", itel);
               skip_subitem(iobuf);
               continue;
            }
            if ( tel_bunches3d[itel] == NULL )
            {
               tel_bunches3d[itel] = (struct bunch3d *) malloc(sizeof(struct bunch3d)*nbunches);
               max_bunches3d[itel] = nbunches;
            }
            else if ( nbunches > max_bunches3d[itel] )
            {
               tel_bunches3d[itel] = (struct bunch3d *) realloc(tel_bunches3d[itel],sizeof(struct bunch3d)*nbunches);
               max_bunches3d[itel] = nbunches;
            }
            if ( tel_bunches[itel] == NULL )
            {
               fflush(stdout);
               fprintf(stderr,"(Re-) Allocation of memory for photon bunches failed.\n");
               exit(1);
            }
            rc = read_tel_photons3d(iobuf,nbunches,&iarray,&itel,&tel_photons3d[itel],
               tel_bunches3d[itel],&tel_nbunches3d[itel]);
            if ( rc < 0 )
            {
               if ( itel >= 0 && itel < MAXTEL )
               {
                  tel_photons3d[itel] = 0.;
                  tel_nbunches3d[itel] = 0;
               }
               get_item_end(iobuf,&item_header);
               return rc;
            }
            select_bunches3d(tel_bunches3d[itel],&tel_nbunches3d[itel],&tel_photons3d[itel]);
            break;
         case IO_TYPE_MC_PE:
            fflush(stdout);
            fprintf(stderr,"Not handling photo-electrons here.\n");
            skip_subitem(iobuf);
            break;
         default:
            fflush(stdout);
            fprintf(stderr,
               "Fix me: unexpected item type %d in my_hess_mc_phot()\n",type);
            skip_subitem(iobuf);
            break;
      }
   }

   end_read_tel_array(iobuf, &item_header);
   

   /* By default, a full array is in one block. */
   begin_write_tel_array(iobuf,&item_header,iarray);

   for ( itel=0; itel<MAXTEL; itel++ )
   {
      if ( tel_nbunches[itel] > 0 )
      {
         if ( verbose )
            printf("Writing photon bunch list of %d bunches for telescope no. %d in array %d.\n", 
               tel_nbunches[itel], itel, iarray);
         write_tel_photons(iobuf, iarray, itel, tel_photons[itel], tel_bunches[itel],
            tel_nbunches[itel], 0, NULL);
         tel_nbunches[itel] = 0;
      }
      if ( tel_nbunches3d[itel] > 0 )
      {
         if ( verbose )
            printf("Writing 3D photon bunch list of %d bunches for telescope no. %d in array %d.\n", 
               tel_nbunches3d[itel], itel, iarray);
         write_tel_photons3d(iobuf, iarray, itel, tel_photons3d[itel], tel_bunches3d[itel],
            tel_nbunches3d[itel], 0, NULL);
         tel_nbunches3d[itel] = 0;
      }
   }

   if ( end_write_tel_array(iobuf,&item_header) != 0 )
      ioerrorcheck();
      
   return 0;
}

void ioerrorcheck (void)
{
   static size_t ioerrcnt = 0;
   
   if ( ++ioerrcnt > 100 )
   {
      fflush(stdout);
      fprintf(stderr,"\n\nFatal error:\n"
          "Too many errors on the IACT output file were detected.\n"
          "It makes little sense to continue with CORSIKA when no output\n"
          "can be written or the output file is badly corrupted.\n"
          "This situation could arise when your disk is full or the output is\n"
          "piped into a detector simulation program and that happened to fail.\n");
      fflush(NULL);
      exit(9);
   }
}

void syntax(void);

void syntax()
{
   fprintf(stderr,"select_iact: Select photons by emitting particle in CORSIKA IACT data.\n");
   fprintf(stderr,"Syntax: select_iact [ options ] { - | filename ... }\n");
   fprintf(stderr,"Options:\n");
   fprintf(stderr,"   -s selection   (add selector for emitting particle)\n");
   fprintf(stderr,"   -V             (more verbose output)\n");
   fprintf(stderr,"   -q             (less verbose output)\n");
   fprintf(stderr,"   -o output      (set output file name)\n");
   fprintf(stderr,"There can be multiple '-s' selections, with logical OR applied.\n");
   fprintf(stderr,"Selections can be by name of the particle type (group), e.g.\n");
   fprintf(stderr,"   electron/positron/e+-/mu-/mu+/proton/nucleus/...\n");
   fprintf(stderr,"but for a generic selection it should provide five numbers:\n");
   fprintf(stderr,"   minimum mass [GeV/c**2], maximum mass [GeV/c**2],\n");
   fprintf(stderr,"   minumum energy [GeV], maximum energy [GeV], charge number.\n");
   fprintf(stderr,"If minimum and maximum mass and/or energy are both zero or the\n");
   fprintf(stderr,"the charge number is zero, that aspect of the selection gets ignored.\n");
   fprintf(stderr,"For example, '-s 0.13,0.14,0,0,0' selects muons of either charge sign\n");
   fprintf(stderr,"and of any energy, while '-s 0.13,0.14,2,10,1' selects only mu+\n");
   fprintf(stderr,"with energies between 2 and 10 GeV.\n");
   fprintf(stderr,"\nNote that the selection process can only work if the CORSIKA binary\n");
   fprintf(stderr,"was compiled with IACTEXT defined, and then run with 'IACT STORE-EMITTER'\n");
   fprintf(stderr,"enabled in the CORSIKA inputs cards, in order to write the extra\n");
   fprintf(stderr,"data needed for the selection.\n");
   exit(1);
}

/* -------------------- main program ---------------------- */
/** 
 *  @short Main program 
 *
 *  Main program function of select_iact program.
 */

int main (int argc, char **argv)
{
   IO_BUFFER *iobuf = NULL;
   IO_ITEM_HEADER item_header;
   const char *input_fname = NULL;
   const char *output_fname = NULL;
   int iarg;
   int itel;

   if ( (iobuf = allocate_io_buffer(5000000L)) == NULL )
   {
      Error("Cannot allocate I/O buffer");
      exit(1);
   }
   if ( iobuf->max_length < 1000000000L )
      iobuf->max_length = 1000000000L; /* Set a maximum of 1 GB */

   if ( argc < 2 )
      syntax();

   /* ================= Big loop over all input files =============== */

   for ( iarg=1; iarg<argc; iarg++ )
   {
      if ( strcmp(argv[iarg],"-v") == 0 )
      {
	 verbose++;
         continue;
      }
      else if ( strcmp(argv[iarg],"-V") == 0 )
      {
         setenv("PRINT_TEL_VERBOSE","0",1);
	 verbose++;
         continue;
      }
      else if ( strcmp(argv[iarg],"-q") == 0 )
      {
         setenv("PRINT_TEL_VERBOSE","0",1);
	 verbose--;
         continue;
      }
      else if ( (strcmp(argv[iarg],"-o") == 0 || strcmp(argv[iarg],"--output-file") == 0)
             && iarg+1 < argc )
      {
         output_fname = argv[++iarg];
         iobuf->output_file = fileopen(output_fname,WRITE_BINARY);
         if ( iobuf->output_file == NULL )
         {
            perror(output_fname);
            exit(1);
         }
         printf("\nOutput file '%s' has been opened.\n", output_fname);
         continue;
      }
      else if ( (strcmp(argv[iarg],"-s") == 0 || strcmp(argv[iarg],"--select") == 0)
             && iarg+1 < argc )
      {
         iarg++;
         char *s = argv[iarg];
         printf("Selector requested: %s\n", s);
         if ( strcmp(s,"electron") == 0 || strcmp(s,"e-") == 0 )
            add_selector(0.510e-3,0.512e-3,0.,0.,-1);
         else if ( strcmp(s,"positron") == 0 || strcmp(s,"e+") == 0 )
            add_selector(0.510e-3,0.512e-3,0.,0.,1);
         else if ( strcmp(s,"e+-") == 0 )
            add_selector(0.510e-3,0.512e-3,0.,0.,0);
         else if ( strcmp(s,"mu-") == 0 || strcmp(s,"muon-") == 0 )
            add_selector(105e-3,107e-3,0.,0.,-1);
         else if ( strcmp(s,"mu+") == 0 || strcmp(s,"muon+") == 0 )
            add_selector(105e-3,107e-3,0.,0.,1);
         else if ( strcmp(s,"mu+-") == 0 || strcmp(s,"muon+-") == 0 || strcmp(s,"muon") == 0)
            add_selector(105e-3,107e-3,0.,0.,0);
         else if ( strcmp(s,"pion") == 0 )
            add_selector(130e-3,140e-3,0.,0.,0);
         else if ( strcmp(s,"proton") == 0 )
            add_selector(0.930,0.940,0.,0.,1);
         else if ( strcmp(s,"nucleus") == 0 )
            add_selector(3.7,400.,0.,0.,0); /* Might include some other short-lived particles */
         else
         {
            double m1=0., m2=0., e1=0., e2=0.;
            int c=0;
            if ( sscanf(s,"%lf,%lf,%lf,%lf,%d", &m1, &m2, &e1, &e2, &c) < 2 )
            {
               fprintf(stderr,"Invalid selector: %s\n", s);
               exit(1);
            }
            printf("Add selector for %f <= M <= %f GeV/c**2, %f < E < %f GeV, charge %d\n",
               m1, m2, e1, e2, c);
            add_selector(m1, m2, e1, e2, c);
         }
         continue;
      }
      else if ( strcmp(argv[iarg],"-h") == 0 || strcmp(argv[iarg],"--help") == 0 )
         syntax();
      else if ( argv[iarg][0] == '-' && argv[iarg][1] != '\0' )
         syntax();

      input_fname = argv[iarg];
      if ( strcmp(input_fname ,"-") == 0 )
      {
         iobuf->input_file = stdin;
      }
      else 
      {
         if ( (iobuf->input_file = fileopen(input_fname,READ_BINARY)) == NULL )
         {
            perror(input_fname);
            fprintf(stderr, "Cannot open input file.");
            break;
         }
      }

      fflush(stdout);
      printf("\nInput file '%s' has been opened.\n",input_fname);

      for (;;) /* Loop over all data in the input file */
      {
         /* Find and read the next block of data. */
         /* In case of problems with the data, just give up. */
         if ( find_io_block(iobuf,&item_header) != 0 )
            break;

         if ( read_io_block(iobuf,&item_header) != 0 )
            break;

 
         /* What did we actually get? */
         switch ( (int) item_header.type )
         {
            /* =================================================== */
            /* That is data for a whole array packed into a common envelope */
            case IO_TYPE_MC_TELARRAY:
               array_select_mc_phot(iobuf);
               break;

            /* =================================================== */
            /* With extended output option activated, the particles
               arriving at ground level would be stored as seemingly
               stray photon bunch block. */
            /* Treat bunches in 3D format just the same as in the traditional format. */
            case IO_TYPE_MC_PHOTONS:
            case IO_TYPE_MC_PHOTONS3D:
               if ( item_header.ident == 999999 ) /* Particle data block not subject to selection process */
               {
                  if ( write_io_block(iobuf) != 0 )
                  {
                     fflush(NULL);
                     fprintf(stderr,"Error writing output data block of type %d to file '%s'.\n",
                        (int) item_header.type, output_fname);
                     exit(2);
                  }
                  continue;
               }
               if ( item_header.type == IO_TYPE_MC_PHOTONS )
                  tel_select_mc_phot(iobuf);
               else if ( item_header.type == IO_TYPE_MC_PHOTONS3D )
                  tel_select_mc_phot3d(iobuf);

               break;

            /* =================================================== */
            case IO_TYPE_MC_TELARRAY_HEAD:
               /* Reset photon bunch contents for all telescopes */
               for ( itel=0; itel<MAXTEL; itel++ )
               {
                  tel_nbunches[itel] = 0;
                  tel_nbunches3d[itel] = 0;
                  tel_photons[itel] = 0.;
                  tel_photons3d[itel] = 0.;
               }
               /* Apart from that write the data block as-is to output file. */
               if ( write_io_block(iobuf) != 0 )
               {
                  fflush(NULL);
                  fprintf(stderr,"Error writing output data block of type %d to file '%s'.\n",
                     (int) item_header.type, output_fname);
                  exit(2);
               }
               break;

            default:
               /* All other data blocks are written as-is, without any further action. */
               if ( write_io_block(iobuf) != 0 )
               {
                  fflush(NULL);
                  fprintf(stderr,"Error writing output data block of type %d to file '%s'.\n",
                     (int) item_header.type, output_fname);
                  exit(2);
               }

         } /* switch for data types */

      } /* Loop over all data in a file */
    
      /* ================ Done with this input data file ============== */

      if ( iobuf->input_file != NULL && iobuf->input_file != stdin )
         fileclose(iobuf->input_file);
      iobuf->input_file = NULL;
      reset_io_block(iobuf);

   } /* Loop over all input files */

   if ( iobuf != NULL )
   {
      free_io_buffer(iobuf);
      iobuf = NULL;
   }

   return 0;
}

/** @} */
