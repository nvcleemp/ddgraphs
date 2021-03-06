/* util.h
 * =========================================================================
 * This file is part of the ddgraphs project
 *
 * Copyright (C) 2011 - 2013 Universiteit Gent
 *
 * Author: Nicolas Van Cleemput
 * In collaboration with Gunnar Brinkmann
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * A copy of the GNU General Public License can be found in the file
 * LICENSE.txt provided with this distribution. This license can also
 * be found on the GNU website at http://www.gnu.org/licenses/gpl.html.
 *
 * If you did not receive a copy of the GNU General Public License along
 * with this program, contact the lead developer, or write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */

#ifndef _UTIL_H //if not defined
#define _UTIL_H

#define HALFFLOOR(n) ((n)%2==0 ? (n)/2 : ((n)-1)/2)
#define MAX(a, b) (((a) >= (b)) ? (a) : (b))
#define MIN(a, b) (((a) <= (b)) ? (a) : (b))

#if WORDSIZE == 64
#define BIT(i) (1ULL << (i))
#else //i.e. WORDSIZE == 32
#define BIT(i) (1U << (i))
#endif

#define ERRORMSG(msg) { fprintf(stderr, "%s:%u %s\n", __FILE__, __LINE__, msg); fflush(stderr); exit(1); }

/******************Debugging macros**********************/

//#define _DEBUG

#ifdef _DEBUG

#define _DEBUGMESSAGES
#define _DEBUGDUMPS
#define _DEBUGASSERTS
#define _DEBUGINTERMEDIATE
#define _DEBUGMETHODS
#define _DEBUGTRACING

#endif

//=============== MESSAGE MACRO'S ===============

#ifdef _DEBUGMESSAGES

#define DEBUGMSG(msg) { fprintf(stderr, "%s:%u %s\n", __FILE__, __LINE__, msg); fflush(stderr); }
#define DEBUGCONDITIONALMSG(condition, msg) if(condition){ fprintf(stderr, "%s:%u %s\n", __FILE__, __LINE__, msg); fflush(stderr);}

#else

#define DEBUGMSG(msg)
#define DEBUGCONDITIONALMSG(condition, msg)

#endif

//=============== TRACING MACRO'S ===============

#ifdef _DEBUGTRACING

#define DEBUGTRACE_ENTER { fprintf(stderr, "%s:%u Entering %s\n", __FILE__, __LINE__, __func__); fflush(stderr); }
#define DEBUGTRACE_EXIT { fprintf(stderr, "%s:%u Exiting %s\n", __FILE__, __LINE__, __func__); fflush(stderr); }

#else

#define DEBUGTRACE_ENTER
#define DEBUGTRACE_EXIT

#endif

//=============== DUMP MACRO'S ===============

#ifdef _DEBUGDUMPS

#define DEBUGDUMP(var, format) { fprintf(stderr, "%s:%u %s=" format "\n", __FILE__, __LINE__, #var, var); fflush(stderr); }

#define DEBUGARRAYDUMP(var, size, format) { \
                                            if(size > 0) {\
                                                fprintf(stderr, "%s:%u %s= [" format, __FILE__, __LINE__, #var, var[0]);\
                                                if(size > 0) {\
                                                    int debugarraydumpcounter;\
                                                    for(debugarraydumpcounter=1; debugarraydumpcounter<size-1; debugarraydumpcounter++){ \
                                                        fprintf(stderr, ", " format, var[debugarraydumpcounter]);\
                                                    }\
                                                    fprintf(stderr, ", " format "]", var[size-1]);\
                                                 }\
                                                 fprintf(stderr,"\n"); fflush(stderr);\
                                            } else {\
                                                fprintf(stderr, "%s:%u %s= []\n", __FILE__, __LINE__, #var);\
                                            }\
                                          }

#define DEBUG2DARRAYDUMP(var, size1, size2, format) { \
                                            fprintf(stderr, "%s:%u %s=\n", __FILE__, __LINE__, #var);\
                                            int debug2darraydumpcounter1,debug2darraydumpcounter2;\
                                            for(debug2darraydumpcounter1=0; debug2darraydumpcounter1<size1; debug2darraydumpcounter1++){ \
                                              fprintf(stderr, "[" format, var[debug2darraydumpcounter1][0]);\
                                              for(debug2darraydumpcounter2=1; debug2darraydumpcounter2<size2-1; debug2darraydumpcounter2++){ \
                                                 fprintf(stderr, ", " format, var[debug2darraydumpcounter1][debug2darraydumpcounter2]);\
                                              }\
                                              fprintf(stderr, ", " format "]\n", var[debug2darraydumpcounter1][size2-1]);\
                                            }\
                                            fflush(stderr);\
                                          }

#else

#define DEBUGDUMP(var, format)
#define DEBUGARRAYDUMP(var, size, format)
#define DEBUG2DARRAYDUMP(var, size1, size2, format)

#endif

//=============== ASSERT MACRO'S ===============

#ifdef _DEBUGASSERTS

#define DEBUGASSERT(assertion) if(!(assertion)) {fprintf(stderr, "%s:%u Assertion failed: %s\n", __FILE__, __LINE__, #assertion); fflush(stderr); exit(1);}

#define DEBUGASSERTMSG(assertion, msg) if(!(assertion)) {fprintf(stderr, "%s:%u Assertion failed: %s\n", __FILE__, __LINE__, #assertion);\
                                                         fprintf(stderr, "%s:%u %s\n", __FILE__, __LINE__, msg); exit(1);}

#define DEBUGASSERTMSGDUMP(assertion, msg, var, format) if(!(assertion)) {fprintf(stderr, "%s:%u Assertion failed: %s\n", __FILE__, __LINE__, #assertion);\
                                                         fprintf(stderr, "%s:%u %s: " format "\n", __FILE__, __LINE__, msg, var); exit(1);}

#else

#define DEBUGASSERT(assertion)
#define DEBUGASSERTMSG(assertion, msg)
#define DEBUGASSERTMSGDUMP(assertion, msg, var, format)

#endif

//=============== PROFILING MACRO'S ===============

#ifdef _PROFILING

#define PROFILINGINCREMENT(counter) counter++;

#else

#define PROFILINGINCREMENT(counter)

#endif

//================== Timing =======================

#define time_factor sysconf(_SC_CLK_TCK)

//=============== Common defines ==================

//defined by nauty
//typedef int boolean;

//#define TRUE 1;
//#define FALSE 0;

#endif // end if not defined, and end the header file
