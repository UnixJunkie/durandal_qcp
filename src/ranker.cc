/**
 *  Copyright (C) 2011, Zhang Initiative Research Unit,
 *  Advance Science Institute, Riken
 *  2-1 Hirosawa, Wako, Saitama 351-0198, Japan
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 4 of the License, or
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
 *****************************************************************************/

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "SimpPDB.h"
#include "Stru.h"
#include "rmsd.h"
#include "qcprot.h"
#include "Singleton.h"

int _mLen        = 0;
float* _previous = NULL;
double** _mPOS1  = NULL;
double** _mPOS2  = NULL;
RMSD_method _method = KABSCH;

// WARNING: this function is not thread safe, not good for OpenMP
float ca_rmsd(float* coor1, float* coor2, size_t idx1, size_t idx2,
              bool crazy) {

  if (coor1 == coor2) {
    return 0.0;
  }

  double rmsd;
  int k3;

  if (coor1 == _previous) {
    for (int k = 0; k < _mLen; ++k) {
      k3 = k*3;
      _mPOS2[0][k] = coor2[k3];
      _mPOS2[1][k] = coor2[k3 + 1];
      _mPOS2[2][k] = coor2[k3 + 2];
    }
  } else {
    for (int k = 0; k < _mLen; ++k) {
      k3 = k*3;
      _mPOS1[0][k] = coor1[k3];
      _mPOS1[1][k] = coor1[k3 + 1];
      _mPOS1[2][k] = coor1[k3 + 2];
      _mPOS2[0][k] = coor2[k3];
      _mPOS2[1][k] = coor2[k3 + 1];
      _mPOS2[2][k] = coor2[k3 + 2];
    }
  }

  switch (_method) {
  case THEOBALD:
    rmsd = rmsd_without_rotation_matrix(_mPOS1, _mPOS2, _mLen,
                                        crazy or (coor1 != _previous),
                                        idx1, idx2);
    break;
  case KABSCH:
    fast_rmsd(_mPOS1, _mPOS2, _mLen, &rmsd);
    break;
  }
  _previous = coor1;

  return (float)rmsd;
}

int main (int argc, char** argv) {

  bool crazy = false;

  if (argc < 2) {
    cout << "Output on stdout the CARMSD of each PDB in the list to the\n"
         << "first one in the list (for example, you may want to put as "
         << "first PDB the one for the known structure).\n"
         << "---" << endl;
    //              0        1
    cout << "usage: ./ranker PDB_list [--theobald] [--crazy n]"    << endl;
    cout << "       --theobald : use THEOBALD instead of KABSCH"   << endl;
    cout << "       --crazy n  : n times ca_rmsd method benchmark" << endl;
    return 1;
  }

  ifstream input_stream(argv[1]);
  string current_line;

  if (contains(argc, argv, "--theobald")) {
    _method = THEOBALD;
  }

  crazy = contains(argc, argv, "--crazy");

  if (not input_stream.is_open()) {
    cout << "error: can't read file: " << string(argv[1]) << endl;
    exit(1);
  }

  bool initialized  = false;
  SimPDB* reference = NULL;
  size_t line_number = 0;

  while (getline(input_stream, current_line)) {

    if (not initialized) {
      reference = new SimPDB(current_line.c_str(), false);
      _mLen  = reference->mNumResidue;
      _mPOS1 = MatInit(3, _mLen);
      _mPOS2 = MatInit(3, _mLen);
      initialized = true;
    } else {
      SimPDB* sim = new SimPDB(current_line.c_str(), _mLen, false);
      cout << current_line << ":"
           << ca_rmsd(reference->mCAlpha, sim->mCAlpha, 0, line_number, crazy)
           << '\n';
      if (crazy) {
        int n = atoi(get_option(argc, argv, "--crazy").c_str());
        for (int i = 0 ; i < n ; ++i) {
          ca_rmsd(reference->mCAlpha, sim->mCAlpha, 0, line_number, crazy);
        }
        exit(0);
      }
      delete sim;
    }
    ++line_number;
  }

  delete reference;

  MatDestroy(&_mPOS1);
  MatDestroy(&_mPOS2);

  return 0;
}
