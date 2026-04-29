/* ============================================================================

Copyright (C) 2002, 2003, 2004, 2006, 2007, 2009, 2010, 2011, 2016, 2018, 2022 Konrad Bernloehr (Konrad Bernl&ouml;hr)

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

/** @file cvt3.cc
 *  @short Conversion of eventio histograms to ROOT format.
 *
@verbatim
   Syntax:  hdata2root [ input_file [ output_file ] ]
       or:  hdata2root -a input_files ... -o output_file
@endverbatim
 *    The program was originally called @c cvt3.
 *    The default input file name is 'testpattern.hdata',
 *    the default output file name is 'testpattern.root' or
 *    the input file name with extension '.root' (instead of '.hdata').
 *    The histograms may be within multiple I/O blocks of the
 *    input file. Only non-empty histograms are written to output.
 *    Take care not to replace any ROOT data format you wanted to keep.
 *
 *    With the '-a' option, all identical histograms in the input
 *    files will be added up before writing them to output.
 *
 *  @author  Konrad Bernloehr
 *  @date    2002 to 2022
 */
 
/** @defgroup cvt3_cc The hdata2root program (cvt3) */
/** @{ */

#include "initial.h"
#include "histogram.h"
#include "io_basic.h"
#include "io_histogram.h"
#include "warning.h"
#include "fileopen.h"
#include "straux.h"

#include "toroot.hh"

#include <vector>

using std::vector;

int read_file (IO_BUFFER *iobuf, const char *fname, int add_flag, int list_flag);

int read_file (IO_BUFFER *iobuf, const char *fname, int add_flag, int list_flag)
{
   FILE *hdata_file;
   IO_ITEM_HEADER item_header;
   int rc, nblocks, nhist=0;

   if ( strcmp(fname,"-") == 0 )
   {
      Warning("Reading data from standard input");
      hdata_file = stdin;
   }
   else if ( (hdata_file = fileopen(fname,READ_BINARY)) == (FILE *) NULL )
   {
      printf("File '%s' not opened\n",fname);
      return -1;
   }
   iobuf->input_file = hdata_file;
   iobuf->max_length = 1000000000;

   nblocks = 0;
   while ( (rc = find_io_block(iobuf,&item_header)) >= 0 )
   {
      int n;
      if ( item_header.type != 100 )
      {
         char message[1024];
         (void) skip_io_block(iobuf,&item_header);
         (void) sprintf(message,"Data in input file are not histograms but type %ld.\n",
              item_header.type);
         Warning(message);
         continue;
      }
      if ( read_io_block(iobuf,&item_header) < 0 )
      {
         Warning("Input data read error.");
         return -1;
      }
      nblocks++;

      if ( list_flag )
         print_histograms(iobuf);

      if ( (n=read_histograms((HISTOGRAM **) NULL, add_flag?-1:0, iobuf )) < 0 )
      {
         Warning("There are problems with the input histograms");
      }
      nhist += n;
   }

   if ( rc == -1 )
      Warning("Input data error. Stop.");
   else if ( nblocks != 1 )
   {
      char message[1024];
      sprintf(message,"End of input data after %d histogram blocks.",nblocks);
      Warning(message);
   }
   else
      printf("Read %d histograms from %s\n",nhist,fname);
   
   if ( hdata_file != stdin )
      fclose(hdata_file);
   
   return 0;
}

int main (int argc, char **argv)
{
   const char *fname1, *fname2;
   char tmpname[1024];
   IO_BUFFER *iobuf;
   int add_flag = 0, add_flag2 = 0, verbose = 0;
   int list_flag = 0;
   vector<int> selected_histograms;

   while ( argc >= 2 )
   {
      if ( argv[1][0] != '-' )
         break;

      if ( strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0 )
      {
         fprintf(stderr,"Conversion of eventio histograms to ROOT format.\n\n");

         fprintf(stderr,"Syntax: %s [-V | -VV ] [ -L ] [ input_file [ output_file ] ]\n"
                        "or:     %s [-V | -VV ] -a input_files ... -o output_file\n"
                        "or:     %s [-V | -VV ] -A -o output_file\n",
                  argv[0], argv[0], argv[0]);
         fprintf(stderr,"The '-V'/'-VV' results in more verbose screen output.\n");
         fprintf(stderr,"The '-L' flag will list histograms as seen on input.\n");
         fprintf(stderr,"With the '-a' option, input histograms are added up.\n");
         fprintf(stderr,"The '-A' is similar but input file names are read from stdin.\n");
         fprintf(stderr,"Other options: --only n1,n2,...\n");
         exit(0);
      }

      if ( strcmp(argv[1],"-V") == 0 || strcmp(argv[1],"--verbose") == 0 )
      {
         argc--;
         argv++;
         verbose++;
         continue;
      }
      else if ( strcmp(argv[1],"-VV") == 0 )
      {
         argc--;
         argv++;
         verbose+=2;
         continue;
      }
      else if ( strcmp(argv[1],"-L") == 0 )
      {
         argc--;
         argv++;
         list_flag = 1;
         continue;
      }
      else if ( strcmp(argv[1],"-a") == 0 )
      {
         add_flag = 1;
         argc--;
         argv++;
         continue;
      }
      else if ( strcmp(argv[1],"-A") == 0 )
      {
         add_flag2 = 1;
         argc--;
         argv++;
         continue;
      }
      else if ( strcmp(argv[1],"--only") == 0 && argc > 2 )
      {
         char word[20];
         int ipos = 0;
         while ( getword(argv[2],&ipos,word,sizeof(word)-1,',','\n') > 0 )
         {
            int k = atoi(word);
            if ( k > 0 )
               selected_histograms.push_back(k);
         }
         argc -= 2;
         argv += 2;
         continue;
      }
      
      fprintf(stderr,"Invalid option '%s'.\n", argv[1]);
      exit(1);
   }

   if ( argc < 2 )
      fname1 = "test.hdata";
   else
      fname1 = argv[1];
   fname2 = "test.root";

   if ( argc < 3 && strcmp(fname2,"-") != 0 )
   {
      size_t l = strlen(fname1);
      strncpy(tmpname,fname1,sizeof(tmpname)-6);
      if ( l>9 && strcmp(tmpname+l-9,".hdata.gz") == 0 )
         strcpy(tmpname+l-9,".root");
      else if ( l>10 && strcmp(tmpname+l-10,".hdata.bz2") == 0 )
         strcpy(tmpname+l-9,".root");
      else if ( l>6 && strcmp(tmpname+l-6,".hdata") == 0 )
         strcpy(tmpname+l-6,".root");
      else
         strcat(tmpname,".root");
      fname2 = tmpname;
   }
   else
      fname2 = argv[2];

   if ( (iobuf = allocate_io_buffer(0)) == (IO_BUFFER *) NULL )
   {
      fprintf(stderr,"No I/O buffer\n");
      exit(1);
   }

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
         read_file(iobuf, argv[iarg], 1, list_flag);
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
         read_file(iobuf, b, 1, list_flag);
      }
   }
   else
      read_file(iobuf, fname1, 0, list_flag);

   sort_histograms();
   if ( get_first_histogram() == 0 )
   {
      fprintf(stderr,"No histograms available for conversion.\n");
      exit(1);
   }
   
   if ( selected_histograms.size() > 0 )
   {
      HISTOGRAM *h = get_first_histogram();
      size_t n = selected_histograms.size();
      while ( h != NULL )
      {
         bool selected = false;
         for ( size_t i = 0; i<n; i++ )
         {
            if ( h->ident == selected_histograms[i] )
            {
               selected = true;
               break;
            }
         }
         HISTOGRAM *hn = h->next;
         if ( !selected )
            free_histogram(h);
         h = hn;
         continue;
      }
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

   /* If the user forgot the '-a' flag he could overwrite his input files. */
   /* Take care for the most obvious case. */
   if ( strstr(fname2,".hdata") != NULL )
   {
      char tmpnam2[1032];
      strcpy(tmpnam2,fname2);
      strcat(tmpnam2,".root");
      strcpy(tmpname,tmpnam2);
      fname2 = tmpname;
   }

   convert_histograms_to_root(fname2);

   exit(0);
   return 0;
}

/** @} */
