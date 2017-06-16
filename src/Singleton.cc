/**
 *  Copyright (C) 2010, Zhang Initiative Research Unit,
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

#include <iostream>

#include "Singleton.h"

Single::Single() {
}

// <never change their implementation> !!! ------------------------------------
Single::~Single() {}

Single::Single(const Single&) {}

Single& Single::operator=(const Single&) {
  return Single::instance();
}
// </never change their implementation> !!!

Strategy cstring_to_strategy(const char* s) {

  if (string(s) == "seq") {
    cout << "using sequential reference picking" << endl;
    return SEQUENTIAL;
  }
  if (string(s) == "maxdiff") {
    cout << "using max_diff reference picking" << endl;
    return MAX_DIFF;
  }
  if (string(s) == "rand") {
    cout << "using random reference picking" << endl;
    return RANDOM;
  }
  if (string(s) == "maxent") {
    cout << "using max entropy reference picking" << endl;
    return MAX_ENTROPY;
  }
  cout << "error: unsupported strategy: " << s << endl;
  exit(1);
}

// return a random number in [lower_bound; upper_bound[
size_t rand_between(size_t lower_bound, size_t upper_bound) {
  // cf. "man rand" from  Linux Programmer's Manual for details
  return (lower_bound + (size_t) (upper_bound * (rand() / (RAND_MAX + 1.0))));
}

// 'man times' for details (Linux Programmer's Manual (2))
clock_t get_user_plus_system_times() {
  struct tms t;
  times(&t);
  return (t.tms_utime + t.tms_stime);
}

// remove i from v
void remove(vector<int>& v, int i) {

  for (vector<int>::iterator it = v.begin() ; it != v.end() ; ++it) {
    if (*it == i) {
      v.erase(it);
      break;
    }
  }
}

// return median value of a sorted vector
float median(vector<float>& v) {

  float result = 0.0;
  size_t end   = v.size();

  if (end > 0) {
    if (end % 2 == 0) {
      result = (v[(int)((end-1) / 2)]  + v[(int)(end / 2)]) / 2.0;
    } else {
      result = v[(int)(end / 2)];
    }
  }

  return result;
}

// return average value of a vector
float average(vector<float>& v) {

  float result = 0.0;
  size_t end   = v.size();

  if (end > 0) {
    for (size_t i = 0 ; i < end ; ++i) {
      result += v[i];
    }
    result /= (float)end;
  }

  return result;
}

bool contains(int argc, char** argv, string param) {

  for (int i = 0 ; i < argc ; ++i) {

    if (param == string(argv[i])) {
      return true;
    }
  }
  return false;
}

bool contains(int argc, char** argv, string short_opt, string long_opt) {
  return (contains(argc, argv, short_opt) or
          contains(argc, argv, long_opt));
}

string get_option(int argc, char** argv, string option_flag) {

  for (int i = 0 ; i < argc ; ++i) {

    if (option_flag == string(argv[i])) {
      return string(argv[i + 1]);
    }
  }
  return string("");
}

string get_option(int argc, char** argv, string short_opt, string long_opt) {

  string res = get_option(argc, argv, short_opt);

  if (res == "") {
    res = get_option(argc, argv, long_opt);
  }

  return res;
}
