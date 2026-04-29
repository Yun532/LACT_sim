/* ============================================================================

Copyright (C) 2010, 2014, 2018  Konrad Bernloehr

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

/**
 *  @file fcat.c
 *  @short Trivial test and utility program for the fileopen/fileclose functions. 
 *
 *  @author Konrad Bernloehr
 *  @date   2010 to 2018
 */

/** @defgroup fact_c The fcat program */
/** @{ */

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

#include "fileopen.h"

#define BSIZE 8192

void syntax(void);

void syntax()
{
   fprintf(stderr,"fcat: read and write (concatenate) files, using the fileopen() functionalities.\n"
        "Syntax: fcat [ -o new-output-file | -a append-to-file ] [ input file(s) ... ]\n"
        "By default reading from standard input and writing to standard input (like 'cat').\n");
   exit(1);
}

int main (int argc, char **argv)
{
   FILE *f = NULL, *output = stdout;
   int iarg;
   size_t n;
   char buffer[BSIZE];
   
   if ( argc > 1 && strcmp(argv[1],"--help") == 0 )
      syntax();
   
   if ( argc > 2 && (strcmp(argv[1],"-o") == 0 || strcmp(argv[1],"-a") == 0) )
   {
      if ( strcmp(argv[1],"-o") == 0 )
         output = fileopen(argv[2],"w");
      else
         output = fileopen(argv[2],"a");
      argc -= 2;
      argv += 2;
      if ( output == NULL )
      {
         perror(argv[0]);
         exit(1);
      }
   }

   for (iarg=1; iarg<argc; iarg++)
   {
      if ( (f = fileopen(argv[iarg],"r")) == NULL )
      {
         if ( errno != 0 )
            perror(argv[iarg]);
         fflush(0);
         fprintf(stderr,"(An error opening the file)\n");
         exit(1);
      }
      while ( (n=fread(buffer,1,BSIZE,f)) == BSIZE )
         fwrite(buffer,1,n,output);
      if ( n > 0 )
         fwrite(buffer,1,n,output);
      if ( ferror(f) )
      {
         if ( errno != 0 )
            perror(argv[iarg]);
         fflush(0);
         fprintf(stderr,"(An error reading from the file)\n");
         if ( output != stdout )
            fileclose(output);
         exit(1);
      }
      if ( fileclose(f) != 0 )
      {
         if ( errno != 0 )
            perror(argv[iarg]);
         fflush(0);
         fprintf(stderr,"(An error closing the file)\n");
         if ( output != stdout )
            fileclose(output);
         exit(1);
      }
   }

   if ( output != stdout )
      fileclose(output);

   return 0;
}

/** @} */
