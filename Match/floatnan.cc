//   Match 2.0.1  -- aligns data series using dynamic programming
//   Copyright (C) 2001-2003  Lorraine E. Lisiecki and Philip A. Lisiecki 
//
//This program is free software; you can redistribute it and/or
//modify it under the terms of the GNU General Public License
//as published by the Free Software Foundation; version 2
//of the License.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
//02111-1307, USA.
//
//The original authors reserve the right to license this program or
//modified versions of this program under other licenses at our
//discretion.
//
//If you publish results generated by this software please cite
//Lisiecki, L. E. and P. A. Lisiecki, The application of dynamic
//  programming to the correlation of paleoclimate records, 
//  submitted to Paleoceanography, 2001
//
//Any questions regarding this license or the operation of this
//software may be directed to Lorraine Lisiecki
//<zogalum@alum.mit.edu>.

#include <iostream>
#include <stdio.h>
#include "floatnan.hh"

namespace std {}
using namespace std;


floatnan::floatnan(string s, int *position) {
  int dummy_pos=0;
  if(!position) position=&dummy_pos;
  int n, j;
  float b;
  n=sscanf(s.c_str()+*position, "%f%n", &b, &j);
  if(n==1) {
    *position+=j;
    value=b;
    nan=false;
  } else {
    char dummy;
    n=sscanf(s.c_str()+*position, "%*1[Nn]%*1[Aa]%1[Nn]%n", &dummy, &j);
    if(n!=1)
      throw parse_error();
    *position+=j;
    value=0.0;
    nan=true;
  }
}

ostream &operator<<(ostream &os, const floatnan &fn) {
  if(fn.nan)
    os<<"NaN";
  else
    os<<fn.value;
  return os;
}

