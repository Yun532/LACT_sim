/* ============================================================================

Copyright (C) 2003, 2009, 2010, 2014, 2018, 2019, 2021, 2022  Konrad Bernloehr

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

/** @file extract_hess.c
 *  @short Extract simulated calibration type event data originally encapsulated.
 *
 *  @author Konrad Bernloehr
 *  @date   2003 to 2022
 */

/** @defgroup extract_hess_c The extract_hess program */
/** @{ */

#include "initial.h"      /* This file includes others as required. */
#include "io_basic.h"     /* This file includes others as required. */
#include "mc_tel.h"
#include "io_history.h"
#include "io_hess.h"
#include "histogram.h"
#include "io_histogram.h"
#include "fileopen.h"

#include <signal.h>

static int interrupted;

/** Show program syntax */

static void syntax (char *program);

static void syntax (char *program)
{
   printf("Extract simulated calibration type event data originally encapsulated.\n\n");
   
   printf("Syntax: %s [ options ] [ - | input_fname ... ]\n",program);
   printf("Options:\n");
   printf("   -o fname     (Set output file name.)\n");
   printf("   -t type      (Extract calibration data of given type.)\n");
   printf("      Types:  0 Pedestals with closed lid (default).\n");
   printf("              1 Pedestals with open lid.\n");
   printf("              2 LED events (separate LED per pixel).\n");
   printf("            >=3 Flatfield 'laser' events of some amplitude level.\n");
   printf("   --dark       (Equivalent to type 0.)\n");
   printf("   --pedestal   (Equivalent to type 1.)\n");
   printf("   --led        (Equivalent to type 2.)\n");
   printf("   --laser n    (Equivalent to type 2+n. Same: --flatfield n)\n");
   printf("   -b nb        (Change maximum size of I/O buffers.)\n");
   exit(1);
}

/* -------------------- main program ---------------------- */
/** 
 *  @short Main program 
 *
 *  Main program function of extract_hess.c program.
 */

int main (int argc, char **argv)
{
   IO_BUFFER *iobuf, *iobuf2;
   IO_ITEM_HEADER item_header, item_header2;
   const char *input_fname = NULL, *output_fname = "iact.simhess.extract";
   // int verbose = 0;
   char *program = argv[0];
   int iarg;
   int rc;
   int type_selected = 0;

   push_command_history(argc,argv);

   if ( argc < 2 )
      input_fname = "iact.simhess";

   if ( (iobuf = allocate_io_buffer(1000000L)) == NULL )
   {
      Error("Cannot allocate I/O buffer");
      exit(1);
   }
   if ( iobuf->max_length < 200000000 )
      iobuf->max_length = 200000000;
   if ( (iobuf2 = allocate_io_buffer(1000000L)) == NULL )
   {
      Error("Cannot allocate I/O buffer 2");
      exit(1);
   }
   if ( iobuf2->max_length < 200000000 )
      iobuf2->max_length = 200000000;

   for (iarg=1; iarg<argc && input_fname == NULL; iarg++)
   {
      if ( strcmp(argv[iarg],"-o") == 0 && iarg+1 < argc )
        output_fname = argv[++iarg];
      else if ( strncmp(argv[iarg],"-o",2) == 0 && strlen(argv[iarg]) > 2 )
        output_fname = argv[iarg]+2;
      // else if ( strcmp(argv[iarg],"-v") == 0 )
      //         verbose = 1;
      else if ( strcmp(argv[iarg],"-t") == 0 && iarg+1 < argc )
         type_selected = atoi(argv[++iarg]);
      else if ( strcmp(argv[iarg],"--dark") == 0 )
         type_selected = 0;
      else if ( strcmp(argv[iarg],"--pedestal") == 0 )
         type_selected = 1;
      else if ( strcmp(argv[iarg],"--led") == 0 )
         type_selected = 2;
      else if ( strcmp(argv[iarg],"--laser") == 0 && iarg+1 < argc )
         type_selected = 2 + atoi(argv[++iarg]);
      else if ( strcmp(argv[iarg],"--flatfield") == 0 && iarg+1 < argc )
         type_selected = 2 + atoi(argv[++iarg]);
      else if ( strcmp(argv[iarg],"-b") == 0 && iarg+1 < argc )
         iobuf->max_length = iobuf2->max_length = atol(argv[++iarg]);
      else if ( argv[iarg][0] == '-' && argv[iarg][1] != '\0' )
        syntax(program);
      else
         input_fname = argv[iarg];
    }
    
    fprintf(stderr,"Extracting calibration data of type %d to %s\n",type_selected,output_fname);

    if ( input_fname == 0 )
    {
       Error("No input file.\n");
       syntax(program);
    }
    if ( strcmp(input_fname ,"-") == 0 )
      iobuf->input_file = stdin;
    else if ( (iobuf->input_file = fileopen(input_fname,READ_BINARY)) == NULL )
    {
      perror(input_fname);
      Error("Cannot open input file.");
      exit(1);
    }
    printf("\nInput file '%s' has been opened.\n",input_fname);
    input_fname = NULL;

    if ( (iobuf2->output_file = fileopen(output_fname,WRITE_BINARY)) == NULL )
    {
      perror(output_fname);
      Error("Cannot open output file.");
      exit(1);
    }
    printf("\nOutput file '%s' has been opened.\n",output_fname);

    /* Save the command line history */
    write_history(0,iobuf2);

    for (;;) /* Loop over all data in the input file */
    {
      if ( interrupted )
         break;

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
         case IO_TYPE_SIMTEL_RUNHEADER:    /* 2000 */
         case IO_TYPE_SIMTEL_MCRUNHEADER:  /* 2001 */
         case IO_TYPE_MC_INPUTCFG:       /* 1212 */
         case IO_TYPE_MC_ATMPROF:        /* 1216 */
         case IO_TYPE_SIMTEL_CAMSETTINGS:  /* 2002 */
         case IO_TYPE_SIMTEL_CAMORGAN:     /* 2003 */
         case IO_TYPE_SIMTEL_PIXELSET:     /* 2004 */
         case IO_TYPE_SIMTEL_PIXELDISABLE: /* 2005 */
         case IO_TYPE_SIMTEL_CAMSOFTSET:   /* 2006 */
         case IO_TYPE_SIMTEL_POINTINGCOR:  /* 2007 */
         case IO_TYPE_SIMTEL_TRACKSET:     /* 2008 */
         case IO_TYPE_SIMTEL_MC_PE_SUM:    /* 2026 */
         case IO_TYPE_SIMTEL_TEL_MONI:     /* 2022 */
         case IO_TYPE_SIMTEL_LASCAL:       /* 2023 */
         case IO_TYPE_SIMTEL_RUNSTAT:      /* 2024 */
         case IO_TYPE_SIMTEL_MC_RUNSTAT:   /* 2025 */
         case IO_TYPE_HISTORY:           /* 70 */
         case IO_TYPE_METAPARAM:         /* 75 */
            /* Copy it to output. */
            reset_io_block(iobuf2);
            get_item_begin(iobuf,&item_header);
            copy_item_to_io_block(iobuf2,iobuf,&item_header);
            get_item_end(iobuf,&item_header);
            write_io_block(iobuf2);
            break;

         /* =================================================== */
         case IO_TYPE_SIMTEL_CALIBEVENT: /* 2028 */
            /* Unwrap it or get rid of it. */
            reset_io_block(iobuf2);
            if ( (rc = get_item_begin(iobuf,&item_header)) < 0 )
               continue;
            if ( item_header.ident == type_selected )
            {
               item_header2.type = IO_TYPE_SIMTEL_EVENT;
               if ( (rc = get_item_begin(iobuf,&item_header2)) < 0 )
                  continue;
               copy_item_to_io_block(iobuf2,iobuf,&item_header2);
               get_item_end(iobuf,&item_header2);
               write_io_block(iobuf2);
            }
            get_item_end(iobuf,&item_header);
            break;

         /* =================================================== */
         case IO_TYPE_SIMTEL_CALIB_PE: /* 2034 */
            /* Unwrap it or get rid of it. */
            reset_io_block(iobuf2);
            if ( (rc = get_item_begin(iobuf,&item_header)) < 0 )
               continue;
            if ( item_header.ident == type_selected )
            {
               item_header2.type = IO_TYPE_MC_TELARRAY;
               if ( (rc = get_item_begin(iobuf,&item_header2)) < 0 )
                  continue;
               copy_item_to_io_block(iobuf2,iobuf,&item_header2);
               get_item_end(iobuf,&item_header2);
               write_io_block(iobuf2);
            }
            get_item_end(iobuf,&item_header);
            break;

         /* =================================================== */
         case IO_TYPE_MC_TELARRAY:    /* 1204 */
         case IO_TYPE_SIMTEL_EVENT:     /* 2010 */
         case IO_TYPE_SIMTEL_MC_SHOWER: /* 2020 */
         case IO_TYPE_SIMTEL_MC_EVENT:  /* 2021 */
         case 100: /* histograms */
            /* Get rid of it; not belonging to internal calibration event. */
            break;
      }
   }

   fileclose(iobuf2->output_file);

   return 0;
}

/** @} */
