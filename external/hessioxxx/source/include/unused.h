/* ============================================================================

Copyright (C) 2023  Konrad Bernloehr

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

/** @file unused.h
 *  @short Pre-processor macro definitions used to tell the 
 *         compiler/user/documentation that a function parameter 
 *         may be or definitely is unused. No point in warning about it.
 *         Compiler-dependent.
 *
 *  @author  Konrad Bernloehr 
 *  @date    2023
 */

#ifndef UNUSED_H__LOADED
#define UNUSED_H__LOADED

/* Unused parameters in function declarations:
   For function parameters known to be not used, the variable name is recommended
   to be omitted in C++ - but this might not be supported by every compiler.
   In C that would not work anyway, thus making use of attributes where possible.
   For documentation, syntax checking, and debugging we might still
   want to have the variable around and marked up as unused.
*/

/* Compiler-dependent mark-up of unused parameters; no separation of type and name. */
#ifndef UNUSED_PAR
# if defined(__GNUC__) 
#  define UNUSED_PAR(x) x ## _unused __attribute__((unused)) 
# elif defined(__LCLINT__) 
#  define UNUSED_PAR(x) /*@unused@*/ x 
# else 
#  define UNUSED_PAR(x) x 
# endif
#endif

/* Compiler-dependent mark-up of unused parameters; for C and C++. */
#ifndef UNUSED_PAR2
# if defined(__GNUC__) 
#  define UNUSED_PAR2(x,y) x y ## _unused __attribute__((unused)) 
# elif defined(__LCLINT__) 
#  define UNUSED_PAR2(x, y) x /*@unused@*/ y
# elif defined (__cplusplus)
#  define UNUSED_PAR2(x,y) x 
# else
#  define UNUSED_PAR2(x,y) x y
# endif
#endif

/* Compiler-dependent mark-up of unused parameters; for C++ only. */
#ifndef UNUSED_PAR2XX
# if defined(__LCLINT__) 
#  define UNUSED_PAR2XX(x, y) x /*@unused@*/ y
# elif defined (__cplusplus)
#  define UNUSED_PAR2XX(x,y) x 
# else
#  define UNUSED_PAR2XX(x,y) x y
# endif
#endif

/* In C++17 we also have the [[maybe_unused]] attribute but that is not the same. */

#endif
