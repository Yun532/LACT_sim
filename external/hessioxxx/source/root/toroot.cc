/* ============================================================================

Copyright (C) 2002, 2003, 2005, 2006, 2007, 2009, 2010, 2011, 2018 Konrad Bernloehr (Konrad Bernl&ouml;hr)

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

/** @file toroot.cc
 *  @short Functions for conversion of eventio histograms to ROOT format.
 *
 *  @author  Konrad Bernloehr
 *  @date    2002 to 2018
 */

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <string>
#include <sstream>

extern "C" {
#include "initial.h"
#include "histogram.h"
}
#include "toroot.hh"

using std::string;

string num2str (int i);
string num2str (double d);
template <class T> string num2str (T num);

/// Convert an int to a string using the STL.

string num2str (int i)
{
   std::ostringstream ost;
   ost << i;
   return ost.str();
}

/// Convert a double to a string using the STL.

string num2str (double d)
{
   std::ostringstream ost;
   ost << d;
   return ost.str();
}

/**
 *  Convert various sorts of numbers to a string
 */
 
template <class T> string num2str (T num)
{
   std::ostringstream ost;
   ost << num;
   return ost.str();
}

/* ------------------ convert_histograms_to_root ------------------- */
/**
 *  @short Open a ROOT file for output, convert all histograms known and write to file.
 *
 *  @param  fname  Name of ROOT output file.
 *
 */

void convert_histograms_to_root (const char *fname)
{
   HISTOGRAM *thisto;
   int ihisto;

   TFile hfile(fname,"RECREATE");

   ihisto = 1;
   for ( thisto = get_first_histogram();
         thisto != (HISTOGRAM *) NULL;
         thisto = thisto->next )
      if ( histogram_to_root(ihisto++,thisto) < 0 )
         ihisto--;

  // Save all objects in this file
  hfile.Write();

  // Close the file. Note that this is automatically done when you leave
  // the application.
  hfile.Close();

fprintf(stderr,
  "%d Histograms written to ROOT file '%s.\n",
  ihisto-1,fname);
}

/* -------------------- histogram_to_root ------------------ */
/**
 *  @short Create a ROOT histogram from the eventio histogram.
 *
 *  Create a ROOT histogram and fill it with the contents of
 *  the given histogram, if it contains any entries.
 *  If the histogram has an ID number, it is booked with this
 *  Id. Otherwise, 90000 + a sequential number is used.
 *
 *  @param  ihisto  Histogram sequential number
 *  @param  histo   Histogram pointer
 *
 *  @return  0 (ok),  -1 (invalid histogram)
 *
 */

int histogram_to_root (int ihisto, HISTOGRAM *histo)
{
   const char *title;
   static char tstring[82];
   double x, y;
   double xlow, xhigh, ylow, yhigh;
   int i, j, hnum;
   static double rzone[3] = { -0.5, 1.5, 0.5 };
   TH1F *th1f = 0;
   TH2F *th2f = 0;
   //TH1D *th1d = 0; // For now use float type histograms.
   //TH2D *th2d = 0;

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
   /* use that number as the ROOT histogram number. Otherwise look */
   /* for a unused number. */
   hnum = 90000 + ihisto;
   if ( histo->ident > 0 )
      hnum = histo->ident;
   else // take care of unnumbered histograms
   {
      for (i=hnum; i<100000; i++)
         if ( get_histogram_by_ident(hnum) == 0 )
            break;
      hnum = 1;
   }

   if ( histo->title != (char *) NULL )
   {
      strncpy(tstring,histo->title,80);
      tstring[80] = '\0';
      strcat(tstring,"$");
      title = tstring;
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
      th2f = new TH2F(num2str(hnum).c_str(),title,histo->nbins,xlow,xhigh,histo->nbins_2d,ylow,yhigh);
      if ( histo->type == 'F' || histo->type == 'D' )
      {
         int ixzone, iyzone, iq;
         double co;
         for (ixzone=0; ixzone<3; ixzone++)
            for (iyzone=0; iyzone<3; iyzone++)
               if ((iq=ixzone+3*iyzone)<8)
                  if ( (co=(histo->extension)->content_outside[iq]) != 0. )
                  {
                     if ( (co == co) &&
                          co < 1e35 && co > -1e35 ) /* Check for NaN or Infty */
                        th2f->Fill(xlow+rzone[ixzone]*(xhigh-xlow),
                           ylow+rzone[iyzone]*(yhigh-ylow),co);
                     else
                        fprintf(stderr,
                            "Bad value outside boundaries in histogram %ld\n",
                            histo->ident);
                  }
      }
      else
      {
         if ( histo->underflow > 0 )
            th2f->Fill(xlow-(xhigh-xlow),0.5*(ylow+yhigh),
               histo->underflow);
         if ( histo->overflow > 0 )
            th2f->Fill(xhigh+(xhigh-xlow),0.5*(ylow+yhigh),
               histo->overflow);
         if ( histo->underflow_2d > 0 )
            th2f->Fill(0.5*(xlow+xhigh),ylow-(yhigh-ylow),
               histo->underflow_2d);
         if ( histo->overflow_2d > 0 )
            th2f->Fill(0.5*(xlow+xhigh),yhigh+(yhigh-ylow),
               histo->overflow_2d);
      }

      {
         struct Histogram_Extension *he = histo->extension;
         int nb = histo->nbins;
         double w = 0.;
         int bad = 0;
         double yscale = (yhigh-ylow)/histo->nbins_2d;
         for (i=0; i<histo->nbins; i++)
         {
            x = (xlow + (i+0.5)*(xhigh-xlow)/histo->nbins);
            for (j=0; j<histo->nbins_2d; j++)
            {
               y = (ylow + (j+0.5)*yscale);
               if ( histo->type == 'D' )
                  w = he->ddata[j*nb+i];
               else if ( histo->type == 'F' )
                  w = he->fdata[j*nb+i];
               else
                  w = histo->counts[j*nb+i];
               if ( !(w == w)  ||
                     w > 1e35 || w < -1e35) /* Check for NaN or +-Infinity */
               {
                  bad++;
                  w = 0.;
               }
               if ( w != 0. )
                  th2f->Fill(x,y,w);
            }
         }
         if ( bad > 0 )
            fprintf(stderr,"%d bad entr%s in histogram %ld set to zero\n",
                bad,(bad==1)?"y":"ies",histo->ident);
      }
   }
   else  /* It's a 1-D histogram */
   {
      th1f = new TH1F(num2str(hnum).c_str(),title,histo->nbins,xlow,xhigh);
      double w = 0.;
      int bad = 0;
      for (i=0; i<histo->nbins; i++)
      {
         x = (xlow + (i+0.5)*(xhigh-xlow)/histo->nbins);
         if ( histo->type == 'F' )
            w = (histo->extension)->fdata[i];
         else if ( histo->type == 'D' )
            w = (histo->extension)->ddata[i];
         else
            w =  histo->counts[i];
         if ( !(w == w)  ||
               w > 1e35 || w < -1e35) /* Check for NaN or +-Infinity */
         {
            bad++;
            w = 0.;
         }
         th1f->Fill(x,w);
      }
      if ( histo->type == 'F' || histo->type == 'D' )
      {
         th1f->Fill(xlow-(xhigh-xlow),
            (histo->extension)->content_outside[0]);
         th1f->Fill(xhigh+(xhigh-xlow),
            (histo->extension)->content_outside[1]);
      }
      else
      {
         if ( histo->underflow > 0 )
            th1f->Fill(xlow-(xhigh-xlow),histo->underflow);
         if ( histo->overflow > 0 )
            th1f->Fill(xhigh+(xhigh-xlow),histo->overflow);
      }
   }

   return 0;
}
