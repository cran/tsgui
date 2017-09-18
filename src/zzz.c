/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2017 -- 2017 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/


#include <R_ext/Rdynload.h> 
// #include "tsgui.h"
#define NULL 0

  
#define CALLDEF_DO(name, n) {#name, (DL_FUNC) &name, n}
static R_CallMethodDef callMethods[]  = {
  // in die respectiven C-Dateien muss tsgui.h eingebunden sein
  //  CALLDEF_DO(tsgui, 13),
  //  CALLDEF_DO(),
  {NULL, NULL, 0}
};




#define CALLABLE(FCTN)  R_RegisterCCallable("tsgui", #FCTN, (DL_FUNC)  FCTN)
void R_init_tsgui(DllInfo  *dll) {
  R_registerRoutines(dll, NULL, // .C
		     NULL, // callMethods,
		     NULL, // .Fortran
		     NULL); // ext
  R_useDynamicSymbols(dll, FALSE);
}



void R_unload_tsgui(DllInfo *info) {
  /* Release resources. */
}

