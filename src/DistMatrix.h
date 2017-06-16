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
 ******************************************************************************
 *
 * A distance matrix containing ranges of distances.
 * Only its upper half (excluding diagonal of zeroes) is filled.
 * This matrix is diagonally symmetric.
 */

#ifndef DIST_MATRIX_H
#define DIST_MATRIX_H

#include <deque>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "DistRange.h"
#include "Singleton.h"
#include "Stru.h"

using namespace std;

enum NCN_OR_NFN {
  NCN,
  NFN
};

class DistMatrix {

 private:

  DistRange _DIST_RANGE_ZERO;
  int _nb_measured;
  int _nb_guessed;
  int _nb_undecided;
  int _nb_pairs;
  clock_t _start;
  size_t _tabu_FIFO_max_size;
  DistRange** _matrix;
  float _clustering_radius;
  float _energy_present;
  size_t _nb_PDBs;
  size_t _previous_i_used;
  size_t _nb_energies;
  int _last_brute_delta;
  int _last_smart_delta;
  int _nb_ca_rmsd_calls;
  int _smart_delta_nb_undecided;
  int _brute_delta_nb_undecided;
  int _smart_delta_t;
  int _brute_delta_t;
  int _nb_references_used;
  RMSD_method _method;
  vector<Stru*>* _read_pdbs;
  vector<int> _all_PDBs_index;
  vector<int> _nb_undecided_neighbors;
  vector<int> _references;
  deque<int> _tabu_FIFO;
  vector< vector<float> > _index_to_energy;
  vector<string> _index_to_pdbname;

  // <for RMSD computations, from Calibur>
  int _mLen;       // nb residues
  double **_mPOS1;
  double **_mPOS2;
  // </for RMSD computations, from Calibur>

  void set(size_t i, size_t j, float exact_distance);
  size_t compare_to_reference();
  void propagate(size_t current_reference);
  void fast_propagate(size_t current_reference);
  void brute_finish();
  void brute_fix_ranges();
  float ca_rmsd(size_t i, size_t j);
  size_t get_next_reference();
  void output_entropy();
  float get_or_measure(size_t i, size_t j);
  float lower_bound_carmsd(int i, int j);
  float upper_bound_carmsd(int i, int j);
  float get_or_compute(size_t i, size_t j);
  float superimposeAndReplace(size_t reference, size_t moving);
  void update_if_necessary(size_t b, size_t c,
                           DistRange& old_range,
                           DistRange& new_range);
  size_t get_highest_entropy_index();
  float sample(vector<float>& sampled_distances);
  void find_neighbors(float dist_to_current,
                      size_t current,
                      vector< pair<float, size_t> >& dists_to_ref,
                      vector<size_t>& neighbors);
  void single_linkage_agglomerative(const int nb_clusters_to_find, 
                                    map<size_t, vector<size_t> >&
                                    clusters_found);

 public:

  DistMatrix(const char* input_filename, float clustering_radius,
             clock_t start);

  DistMatrix(const vector<string>& PDB_files,
             const int nb_clusters_to_find,
             map<size_t, vector<size_t> >& clusters_found);

  virtual ~DistMatrix();

  size_t get_nb_PDBs();
  DistRange& get(size_t i, size_t j);
  void brute_init();
  void smart_init();
  void csv_dump();
  void txt_dump();
  void output_snapshot();
  void dump_neighbors();
  vector<int> get_all_PDBs_index();
  void get_biggest_cluster(const vector<int>& remaining_pdbs,
                           vector<int>& biggest_cluster_found,
                           vector< vector<int> >& pole_position_clusters);
  string& resolve(size_t i);
  float find_threshold(float percent_structures, float percent_dists);
  float find_clustering_distance(float cutoff_ratio, float bin_width,
                                 float* histo_median);
  float compute_D(vector<int>& cluster);
  float compute_D_ziru(vector<int>& cluster);
  float compute_packing(vector<int>& cluster);
  float rand_sample_avg_rmsd();
  void set_clustering_distance(float d);
  void only_rank();
  void sample_and_output(int nb_samples);
  void output_avg_RMSD_to_all();
  void estimate_packing(vector<int>& cluster,
                        double* min_pack, double* max_pack);
  void display_energy(ostream& out, vector<int>& cluster);
  float resolve_energy(size_t pdb_index, size_t energy_index);
  bool has_energies();
  size_t get_nb_energies();

  friend class Tests;
};

#endif /* DIST_MATRIX_H */
