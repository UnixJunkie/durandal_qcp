/*
 *  **************************************************************************
 *  Copyright 2009 Shuai Cheng Li and Yen Kaow Ng
 *  **************************************************************************
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  **************************************************************************/

#include <cmath>

#include "Stru.h"

Stru::Stru(SimPDB* pdb) {

  mPDB = pdb;
  mCAlpha = pdb->mCAlpha;
  mSIG = NULL;
}

void Stru::init_lower_bound_carmsd(int len) {

  if (mSIG == NULL) { // init only once
    mSIG = new float [len];
    for (int i = 0 ; i < len; ++i) {
      mSIG[i] = dist(mCAlpha[3*i], mCAlpha[3*i+1], mCAlpha[3*i+2]);
    }
  }
}

Stru::~Stru() {

  delete mPDB;
  if (mSIG != NULL) {
    delete [] mSIG;
  }
}

float Stru::dist(float x, float y, float z, float *zz) {
  float xd=x-zz[0];
  float yd=y-zz[1];
  float zd=z-zz[2];
  return sqrt(xd*xd+yd*yd+zd*zd);
}

float Stru::dist(float x, float y, float z) {
  return sqrt(x*x+y*y+z*z);
}
