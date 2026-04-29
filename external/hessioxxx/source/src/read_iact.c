/* ============================================================================

Copyright (C) 2018, 2019, 2021, 2022, 2023  Konrad Bernloehr

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

/** @file read_iact.c
 *  @short A program reading simulated CORSIKA data written through the
 *      IACT interface and shows the contents as readable text.
 *
 *  Relevant environment variables:
 *     PRINT_TEL_VERBOSE
 *     MAX_PRINT_ARRAY
 *
 *  @author Konrad Bernloehr
 *
 *  @date   2018 to 2023
 */

/** @defgroup read_iact_c The read_iact program */
/** @{ */

#include "initial.h"      /* This file includes others as required. */
#include "io_basic.h"     /* This file includes others as required. */
#include "io_history.h"
#include "mc_tel.h"
#include "fileopen.h"
// #include "eventio_version.h"
#include <sys/time.h>
#include <strings.h>

int my_print_simtel_mc_phot (IO_BUFFER *iobuf);

/* --------------------- print_simtel_mc_phot ---------------------- */
/**
 *  Print Monte Carlo photons and photo-electrons.
*/  
 
int my_print_simtel_mc_phot (IO_BUFFER *iobuf)
{
   int iarray=0, type;
   int rc;
   IO_ITEM_HEADER item_header;
   if ( (rc = begin_read_tel_array(iobuf, &item_header, &iarray)) < 0 )
      return rc;
   printf("\nMC photon or photo-electron data for array %d\n", iarray);
   while ( (type = next_subitem_type(iobuf)) > 0 )
   {
      switch (type)
      {
         /* The usual photon bunch data block (normal or compact) */
         case IO_TYPE_MC_PHOTONS:
            rc = print_tel_photons(iobuf);
            if ( rc < 0 )
            {
               get_item_end(iobuf,&item_header);
               return rc;
            }
            break;
         /* The 3D bunch format also capable to handle horizontal and upward photons etc. */
         case IO_TYPE_MC_PHOTONS3D:
            rc = print_tel_photons3d(iobuf);
            if ( rc < 0 )
            {
               get_item_end(iobuf,&item_header);
               return rc;
            }
            break;
         /* Photo-electrons in the sensor */
         case IO_TYPE_MC_PE:
            rc = print_photo_electrons(iobuf);
            if ( rc < 0 )
            {
               get_item_end(iobuf,&item_header);
               return rc;
            }
            break;
         /* Anything else cannot be handled here */
         default:
            fprintf(stderr,
               "Fix me: unexpected item type %d in print_simtel_mc_phot()\n",type);
            skip_subitem(iobuf);
      }
   }

   return end_read_tel_array(iobuf, &item_header);
}

void syntax(void);

void syntax()
{
   fprintf(stderr,"read_iact: Read and show CORSIKA IACT data.\n");
   fprintf(stderr,"Syntax: read_iact [ -v ] [ -n max_values ] { - | filename ... }\n");
   fprintf(stderr,"Options:  -v            Verbose output\n"
                  "          -n max_values Change how many items (for example \n"
                  "                        bunches per telescope) are shown.\n");
   fprintf(stderr,"Also note environment variables MAX_PRINT_ARRAY (-n),\n"
                  "PRINT_TEL_VERBOSE (-v), and MAX_UNKNOWN.\n");
   exit(1);
}

void show_header (IO_ITEM_HEADER *item_header);

void show_header (IO_ITEM_HEADER *item_header)
{
   printf("\nItem type %lu, version %u, ident %ld, flag=%d, ext=%d, length %ju\n",
      item_header->type, item_header->version, item_header->ident,
      item_header->user_flag, item_header->use_extension, item_header->length);
}


/* -------------------- main program ---------------------- */
/** 
 *  @short Main program 
 *
 *  Main program function of read_hess.c program.
 */

int main (int argc, char **argv)
{
   IO_BUFFER *iobuf = NULL;
   IO_ITEM_HEADER item_header;
   const char *input_fname = NULL;
   int showdata = 1;
   int iarg;
   int n_unknown = 0, max_unknown=100;
   int verbose = 0;

   if ( (iobuf = allocate_io_buffer(5000000L)) == NULL )
   {
      Error("Cannot allocate I/O buffer");
      exit(1);
   }
   if ( iobuf->max_length < 1000000000L )
      iobuf->max_length = 1000000000L; /* Set a maximum of 1 GB */

   if ( argc < 2 )
      syntax();

   if ( getenv("MAX_UNKNOWN") != NULL )
      max_unknown = atoi(getenv("MAX_UNKNOWN"));

   /* ================= Big loop over all input files =============== */

   for ( iarg=1; iarg<argc; iarg++ )
   {
      if ( strcmp(argv[iarg],"-v") == 0 )
      {
         char t[20];
         verbose++;
         sprintf(t,"%d",verbose);
         setenv("PRINT_TEL_VERBOSE",t,1);
         continue;
      }
      else if ( strcmp(argv[iarg],"-n") == 0 )
      {
         setenv("MAX_PRINT_ARRAY",argv[iarg+1],1);
         iarg++;
         continue;
      }
      else if ( strcmp(argv[iarg],"-h") == 0 || strcmp(argv[iarg],"--help") == 0 )
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

         if ( verbose >= 2 )
            show_header(&item_header);
 
         /* What did we actually get? */
         switch ( (int) item_header.type )
         {
            /* =================================================== */
            case IO_TYPE_MC_TELARRAY:
               if ( showdata )
                  my_print_simtel_mc_phot(iobuf);
               break;

            /* =================================================== */
            case IO_TYPE_MC_TELARRAY_HEAD:
               if ( showdata )
                  printf("\nStart of split photon bunch blocks (id=%ld)\n",
                     item_header.ident);
               break;

            /* =================================================== */
            /* With extended output option activated, the particles
               arriving at ground level would be stored as seemingly
               stray photon bunch block. */
            case IO_TYPE_MC_PHOTONS:
               if ( showdata )
               {
                  printf("\nParticles arriving at observation level:\n");
                  print_tel_photons(iobuf);
               }
               break;

            /* =================================================== */
            /* In some future version, particles might come as 3D. */
            case IO_TYPE_MC_PHOTONS3D:
               if ( showdata )
               {
                  printf("\nParticles arriving at observation level (3D):\n");
                  print_tel_photons3d(iobuf);
               }
               break;

            /* =================================================== */
            case IO_TYPE_MC_TELARRAY_END:
               if ( showdata )
                  printf("\nEnd of split photon bunch blocks (id=%ld)\n",
                     item_header.ident);
               break;

            /* =================================================== */
            case IO_TYPE_MC_LONGI:
               if ( showdata )
                  print_shower_longitudinal(iobuf);
               break;

            /* =================================================== */
            case IO_TYPE_MC_RUNH:
            case IO_TYPE_MC_EVTH:
            case IO_TYPE_MC_EVTE:
            case IO_TYPE_MC_RUNE:
               if ( showdata )
                  print_tel_block(iobuf);
               break;

            /* =================================================== */
            case IO_TYPE_MC_TELPOS:
               if ( showdata )
                  print_tel_pos(iobuf);
               break;

            /* =================================================== */
            case IO_TYPE_MC_TELOFF:
               if ( showdata )
                  print_tel_offset(iobuf);
               break;

            /* =================================================== */
            case IO_TYPE_MC_EXTRA_PARAM:
               if ( showdata )
                  print_shower_extra_parameters(iobuf);
               break;

            /* =================================================== */
            case IO_TYPE_MC_ATMPROF: /* Atmospheric density profile, tabulated and/or 5-layer. */
               if ( showdata )
                  print_atmprof(iobuf);
               break;

            /* =================================================== */
            case IO_TYPE_HISTORY: /* How program was run and how it was configured. */
               if ( showdata )
                  list_history(iobuf,NULL);
               break;

            /* =================================================== */
            case IO_TYPE_METAPARAM: /* Any notable parameters passed? */
               if ( showdata )
                  print_metaparam(iobuf);
               break;

            /* =================================================== */
     	    /* CORSIKA inputs */
	    case IO_TYPE_MC_INPUTCFG: /* CORSIKA inputs as seen by the IACT module */
            {
               struct linked_string corsika_inputs;
               corsika_inputs.text = NULL;
               corsika_inputs.next = NULL;
               read_input_lines(iobuf,&corsika_inputs);
               if ( corsika_inputs.text != NULL )
               {
                  struct linked_string *xl = NULL, *xln = NULL;
                  printf("\nCORSIKA was run with the following input lines:\n");
                  for (xl = &corsika_inputs; xl!=NULL; xl=xln)
                  {
                     printf("   %s\n",xl->text);
                     free(xl->text);
                     xl->text = NULL;
                     xln = xl->next;
                     xl->next = NULL;
                     if ( xl != &corsika_inputs )
                        free(xl);
                  }
               }
            }
	       break;

            default:
               fflush(stdout);
               if ( n_unknown < max_unknown )
                  fprintf(stderr,"Ignoring unknown data block type %ld\n",item_header.type);
               else if ( n_unknown == max_unknown )
                  fprintf(stderr,"More unknown data blocks ignored ...\n");
               n_unknown++;
               break;

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
