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

#ifndef STRU_H
#define STRU_H

#include "SimpPDB.h"

class Stru
{

 public:

  float* mCAlpha;
  float* mSIG; //signature
  SimPDB* mPDB;

  Stru(SimPDB* pdb);
  ~Stru();
  float dist(float x, float y, float z, float *zz);
  float dist(float x, float y, float z);
  void init_lower_bound_carmsd(int len);
};

#endif /* STRU_H */
