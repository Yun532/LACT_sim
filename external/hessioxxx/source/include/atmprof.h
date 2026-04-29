/* ============================================================================

Copyright (C) 2008, 2010, 2019  Konrad Bernloehr

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

/** @file atmprof.h
 *  @short Function prototypes for atmprof.c.
 *
 *  @author  Konrad Bernloehr 
 *  @date    2008 to 2019
 */

#ifdef __cplusplus
extern "C" {
#endif

int init_atmprof (int atmosphere);
int init_atmprof_s (AtmProf *aprof);

double rhofx (double height);
double thickx (double height);
double refidx (double height);
double heighx (double thick);

#ifdef __cplusplus
}
#endif
