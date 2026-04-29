/* ============================================================================

Copyright (C) 1992, 2001, 2006, 2007, 2010, 2014, 2020  Konrad Bernloehr

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
   @file tohbook.h
   @brief Macros and function declarations to call CERN Library HBOOK functions.

   @author  Konrad Bernloehr
   @date    1992 to 2020
*/

#ifndef __TOHBOOK_LOADED

#define __TOHBOOK_LOADED 1

/* +++ SccsId: "@(#)tohbook.h	1.1  5/21/92" +++ */

#ifndef __INITIAL_LOADED
# include "initial.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

void beginhbook_(void);
#define BEGIN_HBOOK() beginhbook_()

void savehbook_(const char *, int *, size_t);
#define SAVE_HBOOK(TITLE,ID) savehbook_(TITLE,&ID,strlen(TITLE))

void hbook1_(int *, const char *, int *, float*, float *, 
  float *, size_t);
#define HBOOK1(ID,CHTITLE,NX,XMI,XMA,VMX) \
  do { int _arg_ID=ID, _arg_NX=NX; \
       float _arg_XMI=XMI, _arg_XMA=XMA, _arg_VMX=VMX; \
       hbook1_(&_arg_ID,CHTITLE,&_arg_NX,&_arg_XMI,&_arg_XMA,&_arg_VMX, \
       strlen(CHTITLE)); } while(0)
void hbook2_(int *, const char *, int *,
    float *, float *, int *, float *, float *, float *, size_t);
#define HBOOK2(ID,CHTITLE,NX,XMI,XMA,NY,YMI,YMA,VMX) \
  do { int _arg_ID=ID, _arg_NX=NX, _arg_NY=NY; \
       float _arg_XMI=XMI, _arg_XMA=XMA, \
         _arg_YMI=YMI, _arg_YMA=YMA, _arg_VMX=VMX; \
       hbook2_(&_arg_ID,CHTITLE,&_arg_NX,&_arg_XMI,&_arg_XMA, \
         &_arg_NY,&_arg_YMI,&_arg_YMA,&_arg_VMX, \
       strlen(CHTITLE)); } while(0)
void hpak_(int *, float *);
#define HPAK(ID,DATA) \
  do { int _arg_ID=ID; hpak_(&_arg_ID,DATA); } while(0)
void hfill_(int *, float *, float *, float *);
#define HFILL(ID,X,Y,WEIGHT) \
  do { int _arg_ID=ID; float _arg_X=X, _arg_Y=Y, _arg_WEIGHT=WEIGHT; \
     hfill_(&_arg_ID,&_arg_X,&_arg_Y,&_arg_WEIGHT); } while(0)

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"

int hexist_(int *);
static int hexist_call (int id);
static int hexist_call (int id) { return hexist_(&id); }
#define HEXIST(ID) \
   hexist_call(ID)

#pragma GCC diagnostic pop

/* Function prototypes of C language interface functions */
void convert_histograms_to_hbook (const char *fname);
int histogram_to_hbook (int ihisto, HISTOGRAM *histo);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
