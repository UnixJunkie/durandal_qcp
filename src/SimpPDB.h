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

#ifndef _SIMP_PDB
#define _SIMP_PDB

#define LONGEST_CHAIN 4000
#include <iostream>

using namespace std;
class SimPDB
{
 public:
  const char* mProteinFileName;
  int mNumResidue;
  //double mSquaredSum;
  float * mCAlpha;
  void read(int expected_count, bool use_URMSD);
  SimPDB(const char* aProteinFileName, bool use_URMSD);
  SimPDB(const char* aProteinFileName, int len, bool use_URMSD);
  ~SimPDB();
};

#endif /* _SIMP_PDB */
