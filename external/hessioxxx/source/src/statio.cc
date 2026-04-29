/* ============================================================================

Copyright (C) 2006, 2010, 2012, 2014, 2019 Konrad Bernloehr (Konrad Bernl&ouml;hr)

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

/** @file statio.cc
 *  @short A program for statistics of eventio data blocks by block type.
 *
@verbatim
Show statistics of EventIO blocks in given files.
Syntax: statio [ -v ] [ -t ] filename [ ... ]
Options:
  -v  Verbose output
  -t  Show total statistics
@endverbatim
 *
 *  @author Konrad Bernloehr
 *  @date 2006 to 2019
 *
 */

/** @defgroup statio_cc The statio program */
/** @{ */

#define __STDC_LIMIT_MACROS 1
#include <map>
using std::map;
#include <iostream>
#include <string>
#include "fileopen.h"
#include "EventIO.hh"

using namespace eventio;
using std::string;
using std::cout;

struct iostats
{
   size_t count;
   uint64_t bytes;
   unsigned version_low;
   unsigned version_high;
   int depth_min;
   int depth_max;
   
   iostats() {count=0; bytes=0; version_low=0; version_high=0; depth_min=0; depth_max=0; }
   iostats(size_t c, uint64_t b, unsigned v, int dp) { count=c; bytes=b;
      version_low = version_high = v; depth_min=depth_max = dp; }
   ~iostats() {}
};

void syntax(const char *prg)
{
   fprintf(stderr,"Show statistics of EventIO blocks in given files.\n");
   fprintf(stderr,"Syntax: %s [ -v ] [ -s ] [ -t ] filename [ ... ]\n", prg); 
   fprintf(stderr,"Options:\n");
   fprintf(stderr,"  -v  Verbose output\n");
   fprintf(stderr,"  -s  Include sub-structure items in statistics\n");
   fprintf(stderr,"  -t  Show total statistics\n");
   exit(1);
}

void add_stats ( map<int,iostats> &sm, int it, uint64_t ln, unsigned iv, int dp)
{
   map<int,iostats>::iterator ti = sm.find(it);
   if ( ti != sm.end() )
   {
      ti->second.count++;
      ti->second.bytes += ln;
      if ( iv < ti->second.version_low )
          ti->second.version_low = iv;
      if ( iv > ti->second.version_high )
          ti->second.version_high = iv;
      if ( ti->second.depth_min > dp )
         ti->second.depth_min = dp;
      if ( ti->second.depth_max < dp )
         ti->second.depth_max = dp;
   }
   else
   {
      sm[it] = iostats(1,ln,iv,dp); /* Add a new type */
   }
}

void collect_substats ( map<int,iostats> &sm, EventIO::Item& item )
{
   int it = item.Type(), its;
   unsigned iv = item.Version();
   size_t ln = item.size() + 12 + (item.Depth()==0?4:0) + (item.IsExtended()?4:0);
   int dp = item.Depth();
   
   add_stats(sm, it, ln, iv, dp);

   if ( item.IsSearchable() )
   {
      while ( item.Status() == 0 && (its=item.NextSubItemType()) > 0 && its < 1000000 )
      {
         EventIO::Item subitem(item,"get",its);
         if ( subitem.Status() == 0 && subitem.Type() > 0 )
         {
            collect_substats(sm,subitem);
         }
         else
            break;
      }
   }
}

int main (int argc, char **argv)
{
   EventIO iobuf;
   bool verbose = false;
   bool totals = false;
   char *prg = argv[0];
   size_t nf = 0;
   bool substat = false;
   size_t max_size=4000000000; /* Only relevant with '-s' */

   map<int,iostats> sm;
   
#ifdef ALWAYS_WITH_REGISTRY
   /* Use default registry of known types with any compiler. */
   /* Requires linking against the library (not available with */
   /* bare-bone eventio, as for CORSIKA IACT/ATMO extension). */
   /* Not necessary when compiled with GCC. */
   set_ev_reg_std();
#endif

   for (int iarg=1; iarg<argc; iarg++)
   {
      if ( string(argv[iarg]) == "-v" )
      {
         verbose = true;
         continue;
      } 
      if ( string(argv[iarg]) == "-t" )
      {
         totals = true;
         continue;
      } 
      if ( string(argv[iarg]) == "-s" )
      {
         substat = true;
         iobuf.SetMaxSize(max_size);
         continue;
      } 
   
      if ( string(argv[iarg]) == "--help" )
      {
         syntax(prg);
      } 

      FILE *f = 0;
      if ( string(argv[iarg]) == "-" )
         f = stdin;
      else
         f = fileopen(argv[iarg],READ_BINARY);
      if ( f == 0 )
      {
         perror(argv[iarg]);
         continue;
      }

      nf++;

      iobuf.OpenInput(f);

      for (;;)
      {
         if ( iobuf.Find() < 0 )
            break;
         unsigned long it = iobuf.ItemType();
         unsigned iv = iobuf.ItemVersion();
         size_t ln = iobuf.size();
         if ( substat )
         {
            if ( iobuf.Read() < 0 )
               break;
            // add_stats(sm, it, ln, iv, 0);
            EventIO::Item top_item(iobuf,"get");
            collect_substats(sm, top_item);
         }
         else
         {
            if ( iobuf.Skip() < 0 )
               break;
            add_stats(sm, it, ln, iv, 0);
         }
      }

      iobuf.CloseInput();
   }

   if ( nf == 0 )
   {
      syntax(prg);
   }
   
   size_t total_c = 0;
   uint64_t total_b = 0;
   if ( sm.size() > 0 && ! verbose )
   {
      cout << "Type\tBlocks\tBytes";
      if ( substat )
      {
         cout << "\t      Depth  Version(s)\tName\n";
      }
      else
         cout << "\t    Version(s)\tName\n";
   }
   
   for (map<int,iostats>::iterator ti = sm.begin(); 
         ti != sm.end(); ++ti)
   {
      const char *name = eventio_registered_typename(ti->first);

      if ( verbose )
      {
         const char *desc = eventio_registered_description(ti->first);
         cout << "Type " << ti->first << ": " 
              << ti->second.count << " blocks with " 
              << ti->second.bytes << " bytes";
         if ( substat )
         {
            if ( ti->second.depth_min == ti->second.depth_max )
               cout << " at depth " << ti->second.depth_min;
            else
               cout << " at depths " << ti->second.depth_min << " to " << ti->second.depth_max;
         }
         if ( ti->second.version_low == ti->second.version_high )
            cout << " (version " << ti->second.version_low << ")";
         else
            cout << " (versions " << ti->second.version_low 
                 << " to " << ti->second.version_high << ")";
         if ( name != NULL && *name != '\0' )
            cout << "\t[" << name << "] " << desc << "\n";
         else
            cout << "\n";
      }
      else
      {
         cout << ti->first << "\t" << ti->second.count << "\t" 
              << ti->second.bytes << "\t"
              << (ti->second.bytes < 10000000 ? "\t" : "");
         if ( substat )
         {
            if ( ti->second.depth_min == ti->second.depth_max )
               cout << ti->second.depth_min << "\t";
            else
               cout << ti->second.depth_min << "-" << ti->second.depth_max << "\t";
         }
         cout << ti->second.version_low;
         if ( ti->second.version_low != ti->second.version_high )
            cout << "-" << ti->second.version_high;
         if ( name != NULL && *name != '\0' )
            cout << "\t[" << name << "]\n";
         else
            cout << "\n";
      }

      if ( ti->second.depth_min == 0 ) /* No double counting for total */
      {
         total_c += ti->second.count;
         total_b += ti->second.bytes;
      }
   }

   if ( totals )
   {
      if ( verbose )
      {
         cout << "Total:" << "\t" << total_c << " blocks containing " << total_b << " bytes";
         if ( total_b >= (1UL << 30) )
            cout << " (" << total_b/(double)(1UL << 30) << " GiB)";
         cout << "\n";
      }
      else
         cout << "Total:" << "\t" << total_c << "\t" << total_b << "\n";
   }

   return 0;
}

/** @} */
