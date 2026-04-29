/* ============================================================================

Copyright (C) 2001, 2003, 2006, 2007, 2010, 2014, 2018  Konrad Bernloehr

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

/** @file tohbook.c
 *  @short Convert my histograms to HBOOK (PAW) histograms.
 *
 *  @author Konrad Bernloehr
 *  @date   2001 to 2018
 */

#include "initial.h"
#include "histogram.h"
#include "tohbook.h"

#ifdef USING_SCCS
static char SccsId[] = "@(#)tohbook.c   1.5  10/18/95";
#endif

/* ------------------ convert_histograms_to_hbook ------------------- */
/*
@@@convert_histograms_to_hbook()
@@    Initialize HBOOK and write all histograms to an HBOOK file.
@*
@*    Arguments:  fname  --  Name of HBOOK output file.
@*
@*    Return value:  (none)
@*
 */

void convert_histograms_to_hbook (const char *fname)
{
   HISTOGRAM *thisto;
   int ihisto;
   int id=0;

   BEGIN_HBOOK();

   ihisto = 1;
   for ( thisto = get_first_histogram();
         thisto != (HISTOGRAM *) NULL;
         thisto = thisto->next )
      if ( histogram_to_hbook(ihisto++,thisto) < 0 )
         ihisto--;

   SAVE_HBOOK(fname,id);

fprintf(stderr,
  "%d Histogramme wurden in HBOOK-Datei '%s' geschrieben.\n",
  ihisto-1,fname);
}

/* -------------------- histogram_to_hbook ------------------ */
/*
@@@histogram_to_hbook()
@@    Create an HBOOK histogram and fill it with the contents of
@@    the given histogram, if it contains any entries.
@@    If the histogram has an ID number, it is booked with this
@@    Id. Otherwise, 90000 + a sequential number is used.
@*
@*    Arguments:  ihisto  --  Histogram sequential number
@*                histo   --  Histogram pointer.
@*
@*    Return value:  0 (ok),  -1 (invalid histogram)
@*
 */

int histogram_to_hbook (int ihisto, HISTOGRAM *histo)
{
   const char *title;
   static char string[82];
   float *data;
   float x, y, z;
   double xlow, xhigh, ylow, yhigh;
   int i, j, hnum;
   static double rzone[3] = { -0.5, 1.5, 0.5 };

   if ( histo == (HISTOGRAM *) NULL )
      return -1;
   if ( histo->entries == 0 || histo->nbins <= 0 )
      return -1;
   if ( histo->type == 'F' || histo->type == 'D' )
   {
      if ( histo->extension == (struct Histogram_Extension *) NULL )
         return -1;
      if ( (histo->type == 'F' && (histo->extension)->fdata == (float *) NULL) ||
           (histo->type == 'D' && (histo->extension)->ddata == (double *) NULL) )
         return -1;
   }
   else if ( histo->counts == (unsigned long *) NULL )
      return -1;

   /* If the histogram has a (so far) unique identification number, */
   /* use that number as the HBOOK histogram number. Otherwise look */
   /* for a unused number. */
   hnum = 90000 + ihisto;
   if ( histo->ident > 0 )
      if ( !HEXIST(histo->ident) )
         hnum = histo->ident;
   for (i=hnum; i<100000; i++)
      if ( !HEXIST(hnum) )
         break;

   if ( histo->title != (char *) NULL )
   {
      strncpy(string,histo->title,80);
      string[80] = '\0';
      strcat(string,"$");
      title = string;
   }
   else
      title = "UNSPECIFIED HISTOGRAM$";
   if ( histo->type == 'I' || histo->type == 'i')
   {
      xlow  = (double) histo->specific.integer.lower_limit;
      xhigh = (double) histo->specific.integer.upper_limit;
   }
   else
   {
      xlow  = (double) histo->specific.real.lower_limit;
      xhigh = (double) histo->specific.real.upper_limit;
   }
   if ( histo->nbins_2d > 0 )  /* It's a 2-D histogram */
   {
     if ( histo->type == 'I' || histo->type == 'i' )
      {
         ylow  = (double) histo->specific_2d.integer.lower_limit;
         yhigh = (double) histo->specific_2d.integer.upper_limit;
      }
      else
      {
         ylow  = (double) histo->specific_2d.real.lower_limit;
         yhigh = (double) histo->specific_2d.real.upper_limit;
      }
      HBOOK2(hnum,title,histo->nbins,xlow,xhigh,histo->nbins_2d,ylow,yhigh,1000000.);
      if ( histo->type == 'F' || histo->type == 'D' )
      {
         int ixzone, iyzone, iq;
         double co;
         for (ixzone=0; ixzone<3; ixzone++)
            for (iyzone=0; iyzone<3; iyzone++)
               if ((iq=ixzone+3*iyzone)<8)
                  if ( (co=(histo->extension)->content_outside[iq]) != 0. )
                  {
                     if ( (co > 0. || co < 0. || co == 0.) &&
                          co < 1e35 && co > -1e35 ) /* Check for NaN or Infty */
                        HFILL(hnum,xlow+rzone[ixzone]*(xhigh-xlow),
                           ylow+rzone[iyzone]*(yhigh-ylow),(float)co);
                     else
                        fprintf(stderr,
                            "Bad value outside boundaries in histogram %ld\n",
                            histo->ident);
                  }
      }
      else
      {
         if ( histo->underflow > 0 )
            HFILL(hnum,xlow-(xhigh-xlow),0.5*(ylow+yhigh),
               (float)histo->underflow);
         if ( histo->overflow > 0 )
            HFILL(hnum,xhigh+(xhigh-xlow),0.5*(ylow+yhigh),
               (float)histo->overflow);
         if ( histo->underflow_2d > 0 )
            HFILL(hnum,0.5*(xlow+xhigh),ylow-(yhigh-ylow),
               (float)histo->underflow_2d);
         if ( histo->overflow_2d > 0 )
            HFILL(hnum,0.5*(xlow+xhigh),yhigh+(yhigh-ylow),
               (float)histo->overflow_2d);
      }
      if ( histo->type == 'F' )
      {
         REGISTER struct Histogram_Extension *he = histo->extension;
         REGISTER int nb = histo->nbins;
         REGISTER double w;
         int bad = 0;
         for (i=0; i<histo->nbins; i++)
            for (j=0; j<histo->nbins_2d; j++) /* Check for NaN or +-Infinity */
               if ( !((w=he->fdata[j*nb+i]) > 0. || w < 0. || w == 0.) ||
                     w > 1e35 || w < -1e35 )
               {
                  bad++;
                  he->fdata[j*nb+i] = 0.;
               }
         if ( bad > 0. )
            fprintf(stderr,"%d bad entr%s in histogram %ld set to zero\n",
                bad,(bad==1)?"y":"ies",histo->ident);
         HPAK(hnum,(histo->extension)->fdata);
      }
      else
      {
         REGISTER struct Histogram_Extension *he = histo->extension;
         REGISTER int nb = histo->nbins;
         REGISTER double w = 0.;
         int bad = 0;
         double yscale = (yhigh-ylow)/histo->nbins_2d;
         for (i=0; i<histo->nbins; i++)
         {
            x = (float) (xlow + (i+0.5)*(xhigh-xlow)/histo->nbins);
            for (j=0; j<histo->nbins_2d; j++)
            {
               y = (float) (ylow + (j+0.5)*yscale);
               if ( histo->type == 'D' )
                  z = (float) (w=he->ddata[j*nb+i]);
               else
                  z = (float) histo->counts[j*nb+i];
               if ( !(w > 0. || w < 0. || w == 0.)  ||
                     w > 1e35 || w < -1e35) /* Check for NaN or +-Infinity */
               {
                  bad++;
                  z = 0.;
               }
               if ( z != (float) 0. )
                  HFILL(hnum,x,y,z);
            }
         }
         if ( bad > 0. )
            fprintf(stderr,"%d bad entr%s in histogram %ld set to zero\n",
                bad,(bad==1)?"y":"ies",histo->ident);
      }
   }
   else  /* It's a 1-D histogram */
   {
      if ( histo->type == 'F' )
         data = (histo->extension)->fdata;
      else
      {
         if ( (data = (float *) malloc(histo->nbins*sizeof(float))) ==
             (float *) NULL )
            return -2;
         if ( histo->type == 'D' )
            for ( i=0; i<histo->nbins; i++ )
               data[i] = (float) (histo->extension)->ddata[i];
         else
            for ( i=0; i<histo->nbins; i++ )
               data[i] = (float) histo->counts[i];
      }
      HBOOK1(hnum,title,histo->nbins,xlow,xhigh,0.);
      HPAK(hnum,data);
      if ( histo->type == 'F' || histo->type == 'D' )
      {
         HFILL(hnum,xlow-(xhigh-xlow),0.,
            (float)(histo->extension)->content_outside[0]);
         HFILL(hnum,xlow+(xhigh-xlow),0.,
            (float)(histo->extension)->content_outside[1]);
      }
      else
      {
         if ( histo->underflow > 0 )
            HFILL(hnum,xlow-(xhigh-xlow),0.,(float)histo->underflow);
         if ( histo->overflow > 0 )
            HFILL(hnum,xhigh+(xhigh-xlow),0.,(float)histo->overflow);
      }
      if ( histo->type != 'F' )
         free(data);
   }

   return 0;
}
