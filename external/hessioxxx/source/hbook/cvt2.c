/* ============================================================================

Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009, 2010, 2011, 2014  Konrad Bernloehr

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

/** @file cvt2.c
 *  @short Utility program for converting histograms to HBOOK format.
@verbatim
Syntax:  hdata2hbook [ input_file [ output_file ] ]
    or:  hdata2hbook -a input_files ... -o output_file
@endverbatim
 *  The program was originally called @c cvt2.
 *  The default input file name is 'testpattern.hdata',
 *  the default output file name is 'testpattern.hbook' or
 *  the input file name with extension '.hbook' (instead of '.hdata').
 *  The histograms may be within multiple I/O blocks of the
 *  input file. Only non-empty histograms are written to output.
 *
 *  With the '-a' option, all identical histograms in the input
 *  files will be added up before writing them to output.
 *
 *  @author Konrad Bernloehr
 *  @date   2001 to 2014
 */

/** @defgroup cvt2_c The hdata2hbook program (cvt2) */
/** @{ */

#include "initial.h"
#include "histogram.h"
#include "io_basic.h"
#include "tohbook.h"
#include "io_histogram.h"
#include "fileopen.h"

/** Main program. */

int main (int argc, char **argv)
{
   const char *fname1, *fname2, *prgm = argv[0];
   char tmpname[1024];
   int add_flag = 0, add_flag2 = 0, verbose = 0, hdata_fmt = 0;
   int list_flag = 0;

   while ( argc >= 2 )
   {
      if ( strcmp(argv[1],"-V") == 0 || strcmp(argv[1],"--verbose") == 0 )
      {
         argc--;
         argv++;
         verbose++;
      }
      else if ( strcmp(argv[1],"-VV") == 0 )
      {
         argc--;
         argv++;
         verbose+=2;
      }
      else if ( strcmp(argv[1],"-L") == 0 )
      {
         argc--;
         argv++;
         list_flag = 1;
      }
      else if ( strcmp(argv[1],"-e") == 0 )
      {
         argc--;
         argv++;
         hdata_fmt = 1;
      }
      else
         break;
   }

   if ( argc >= 2 )
   {
      if ( strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0 )
      {
         fprintf(stderr,"Syntax: %s [-V | -VV ] [ -L ] [ input_file [ output_file ] ]\n"
                        "or:     %s [-V | -VV ] [ -L ] [-e] -a input_files ... -o output_file\n"
                        "or:     %s [-V | -VV ] [ -L ] [-e] -A -o output_file\n",
                  prgm, prgm, prgm);
         fprintf(stderr,"The '-V'/'-VV' results in more verbose screen output.\n");
         fprintf(stderr,"The '-L' flag will list histograms as seen on input.\n");
         fprintf(stderr,"With the '-a' option, input histograms are added up.\n");
         fprintf(stderr,"The '-A' is similar but input file names are read from stdin.\n");
         fprintf(stderr,"With the '-e' option, the output is again in EventIO format.\n");
         fprintf(stderr,"Otherwise the output is in HBOOK format.\n");
         exit(0);
      }
      if ( strcmp(argv[1],"-a") == 0 )
      {
         add_flag = 1;
         argc--;
         argv++;
      }
      else if ( strcmp(argv[1],"-A") == 0 )
      {
         add_flag2 = 1;
         argc--;
         argv++;
      }
   }

   if ( argc < 2 )
      fname1 = "test.hdata";
   else
      fname1 = argv[1];
   if ( hdata_fmt )
      fname2 = "test.hdat-out.gz";
   else
      fname2 = "test.hbook";

   if ( argc < 3 && strcmp(fname2,"-") != 0 )
   {
      size_t l = strlen(fname1);
      strncpy(tmpname,fname1,sizeof(tmpname)-7);
      if ( l>9 && strcmp(tmpname+l-9,".hdata.gz") == 0 )
         strcpy(tmpname+l-9,".hbook");
      else if ( l>10 && strcmp(tmpname+l-10,".hdata.bz2") == 0 )
         strcpy(tmpname+l-9,".hbook");
      else if ( l>6 && strcmp(tmpname+l-6,".hdata") == 0 )
      {
         if ( hdata_fmt )
            strcpy(tmpname+l-6,".hdat-out.gz");
         else
            strcpy(tmpname+l-6,".hbook");
      }
      else
      {
         if ( hdata_fmt )
            strcat(tmpname,".hdat-out.gz");
         else
            strcat(tmpname,".hbook");
      }
      fname2 = tmpname;
   }
   else
      fname2 = argv[2];

   if ( add_flag )
   {
      int iarg;
      for (iarg=1; iarg<argc; iarg++)
      {
         if ( strcmp(argv[iarg],"-o") == 0 )
         {
            if ( iarg+1 < argc )
               fname2 = argv[iarg+1];
            break;
         }
         read_histogram_file(argv[iarg],1+16*list_flag);
      }
   }
   else if ( add_flag2 )
   {
      char line[2048];
      while ( fgets(line,sizeof(line)-1,stdin) != NULL )
      {
         char *b = line;
         size_t l = strlen(line);
         if ( l > 0 )
         {
            if ( line[l-1] == '\n' )
            {
               line[l-1] = '\0';
               l--;
            }
            if ( l == 0 )
               continue;
            while ( l>0 && (line[l-1] == ' ' || line[l-1] == '\t') )
            {
               line[l-1] = '\0';
               l--;
            }
         }
         if ( l == 0 )
            continue;
         while (*b == ' ' || *b == '\t')
            b++;
         read_histogram_file(b,1+16*list_flag);
      }
   }
   else
      read_histogram_file(fname1,0+16*list_flag);

   sort_histograms();
   if ( get_first_histogram() == 0 )
   {
      fprintf(stderr,"No histograms available for conversion.\n");
      exit(1);
   }

   if ( verbose )
   {
      HISTOGRAM *h = get_first_histogram();
      while ( h != NULL )
      {
         printf("Histogram of type %c, ID=%ld, title=\"%s\" is %cD: ",
            h->type, h->ident, (h->title!=NULL)?h->title:"(none)",
            (h->nbins_2d > 0) ? '2' : '1');
         if ( h->nbins_2d > 0)
            printf("%d * %d bins",h->nbins,h->nbins_2d);
         else
            printf("%d bins",h->nbins);
         printf(", %ld entries.\n",h->entries);
         h = h->next;
      }

      if ( verbose >= 2 )
         display_all_histograms(); 
   }
   
   if ( hdata_fmt && (add_flag || add_flag2) )
   {
      write_all_histograms(fname2);
      return 0;
   }

   BEGIN_HBOOK();

   /* If the user forgot the '-a' flag he could overwrite his input files. */
   /* Take care for the most obvious case. */
   if ( strstr(fname2,".hdata") != NULL )
   {
      char tmpnam2[1032];
      strcpy(tmpnam2,fname2);
      strcat(tmpnam2,".hbook");
      strcpy(tmpname,tmpnam2);
      fname2 = tmpname;
   }

   convert_histograms_to_hbook(fname2);

   exit(0);
   return(0);
}

/** @} */
