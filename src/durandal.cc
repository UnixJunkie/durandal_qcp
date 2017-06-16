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

#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "DistMatrix.h"
#include "Singleton.h"

using namespace std;

int usage() {
  cout << "usage: ./durandal -i pdb_list -d clustering_distance -o file\n"
       << "                  {-b|-s|-p}\n"
       << "                  [{-v|--csv|--txt|--entropy}]\n"
       << "  --semi-auto dist_percent : choose cutoff to accept\n"
       << "                             approximately only\n"
       << "                             dist_percent smallest distances\n"
       << "                             ([0.03 .. 0.05] is reasonable)\n"
       << "  -b|--brute   : brute mode\n"
       << "  -m max_clusters : number of clusters to output\n"
       << "                    (default is 3, use -1 for all)\n"
       << "  -o|--output output_file : save clusters to file\n"
       << "  -s|--smart   : smart mode\n"
       << "  -v|--verbose : verbose clusters listing\n"
       << "                 (default is index in input file only)\n"
       << "  -h|--help\n"
       << "\n"
       << "  FOR DEVELOPER ONLY:\n"
       << "  --stable     : output all pole position clusters\n"
       << "  --compute-D  : output D of found clusters\n"
       << "  --csv        : matrix csv dump\n"
       << "  -e|--entropy entropy_file : output 'entropy clock'\n"
       << "                              to entropy_file\n"
       << "  --sample nb_samples : random sampling of pairwise CARMSDs\n"
       << "  -p|--parano  : check consistency of matrices\n"
       << "  --kabsch     : use KABSCH instead of THEOBALD for RMSD\n"
       << "  --txt        : matrix txt dump\n"
       << "################################################################"
       << "####\n"
       << "IF YOU USE THIS SOFTWARE, PLEASE CITE THE CORRESPONDING "
       << "PUBLICATION:\n"
       << "author  = {Francois Berenger, Rojan Shrestha, Yong Zhou,\n"
       << "           David Simoncini and Kam Y. J. Zhang}\n"
       << "title   = {Durandal: fast exact clustering of protein decoys}\n"
       << "journal = {Journal of Computational Chemistry}\n"
       << "################################################################"
       << "####\n";
  exit(1);
}

// output cluster, remove its members from remaining_pdbs, unless this is
// a pole_position cluster
void print_cluster(ostream& out,
                   vector<int>& cluster,
                   vector<int>& remaining_pdbs,
                   DistMatrix& dm,
                   bool verbose,
                   bool pole_position,
                   bool compute_D) {

  int elt = cluster[0];
  bool has_energy = dm.has_energies();

  if (has_energy) {
    dm.display_energy(out, cluster);
  }
  if (compute_D) {
    out << "D:    " << dm.compute_D(cluster) << endl;
  }
  out << "members(" << cluster.size() << "):";
  if (verbose) {
    out << '\n' << dm.resolve(elt);
    if (has_energy) {
      for (size_t energy_index = 0 ; energy_index < dm.get_nb_energies();
           ++energy_index) {
        out << ':' << dm.resolve_energy(elt, energy_index);
      }
    }
  } else {
    out << ' ' << elt;
  }
  if (not pole_position) {
    remove(remaining_pdbs, elt);
  }
  for (size_t i = 1 ; i < cluster.size() ; ++i) {
    elt = cluster[i];
    if (verbose) {
      out << '\n' << dm.resolve(elt);
      if (has_energy) {
        for (size_t energy_index = 0 ; energy_index < dm.get_nb_energies();
             ++energy_index) {
          out << ':' << dm.resolve_energy(elt, energy_index);
        }
      }
    } else {
      out << ' ' << elt;
    }
    if (not pole_position) {
      remove(remaining_pdbs, elt);
    }
  }
  out << '\n';
}

void output_cluster(ostream& out,
                    vector<int>& cluster,
                    vector<int>& remaining_pdbs,
                    vector< vector<int> >& pole_position_clusters,
                    DistMatrix& dm,
                    bool verbose ,
                    bool stable,
                    bool compute_D) {

  if (pole_position_clusters.size() > 0) {
    if (stable) {
      out << "<pole position clusters("
          << pole_position_clusters.size() << ")>:\n";
      for (size_t i = 0 ; i < pole_position_clusters.size() ; ++i) {
        print_cluster(out, pole_position_clusters[i], remaining_pdbs,
                      dm, verbose, true, compute_D);
      }
      out << "</pole position clusters>\n";
    } else {
      out << "pole position centers("
          << pole_position_clusters.size() << "):";
      if (verbose) {
        out << '\n';
        for (size_t i = 0 ; i < pole_position_clusters.size() ; ++i) {
          out << dm.resolve(pole_position_clusters[i][0]) << '\n';
        }
      } else {
        for (size_t i = 0 ; i < pole_position_clusters.size() ; ++i) {
          out << ' ' << pole_position_clusters[i][0];
        }
        out << '\n';
      }
    }
  }
  print_cluster(out, cluster, remaining_pdbs, dm, verbose, false, compute_D);
}

int main(int argc, char** argv) {
  // <KEEP THIS AS FIRST THING IN MAIN>
  clock_t start = get_user_plus_system_times();
  // </KEEP THIS AS FIRST THING IN MAIN>

  if (argc == 1 or contains(argc, argv, "-h", "--help")) {
    usage();
  }

  // seed Random Number Generator
  srand(time(NULL)); // comment   here for repeatability
  // srand(0);       // uncomment here for repeatability

  Single& single = Single::instance();

  bool matrix_csv_dump     = contains(argc, argv, "--csv");
  bool paranoiac_mode      = contains(argc, argv, "-p", "--parano");
  bool matrix_txt_dump     = contains(argc, argv, "--txt");
  bool stable              = contains(argc, argv, "--stable");
  bool auto_threshold      = contains(argc, argv, "--semi-auto");
  bool histo               = false;
  bool sample              = contains(argc, argv, "--sample");
  single._compute_D        = contains(argc, argv, "--compute-D");
  single._only_rank        = false;
  single._smart_mode       = contains(argc, argv, "-s", "--smart");
  single._verbose          = contains(argc, argv, "-v", "--verbose");
  single._brute_mode       = contains(argc, argv, "-b", "--brute");
  single._output_entropy   = contains(argc, argv, "-e", "--entropy");
  single._auto_switch      = true;
  single._smart_fix        = false;
  single._brute_fix        = false;
  single._no_lower_bound   = false;
  single._no_upper_bound   = false;
  single._output_to_file   = contains(argc, argv, "-o", "--output");
  single._use_URMSD        = false;
  single._snapshot_entropy = false;
  single._use_theobald     = not contains(argc, argv, "--kabsch");

  string input_file = get_option(argc, argv, "-i");
  if (input_file == "") {
    cout << "error: -i pdb_list is mandatory" << endl;
    usage();
  }
  string distance = get_option(argc, argv, "-d");
  if (auto_threshold) {
    distance = string("10"); // will be overwritten later
  }
  if (distance == "" and not (histo or sample)) {
    cout << "error: -d clustering_distance is mandatory" << endl;
    usage();
  }
  if (single._output_entropy) {
    string entropy_file = get_option(argc, argv, "-e", "--entropy");
    if (entropy_file == "") {
      cout << "error: -e|--entropy needs a filename parameter" << endl;
      usage();
    }
    single._entropy_file = fopen(entropy_file.c_str(), "w");
    if (single._entropy_file == NULL) {
      cout << "error: cannot write to: " << entropy_file << endl;
      usage();
    }
  }
  string strategy = get_option(argc, argv, "--strat");
  if (strategy != "") {
    single._reference_choosing_strategy = \
      cstring_to_strategy(strategy.c_str());
  } else {
    single._reference_choosing_strategy = \
      cstring_to_strategy("seq");
  }
  string output_file = get_option(argc, argv, "-o", "--output");
  if (not (histo or sample)) {
    if (output_file == "") {
      cout << "error: -o is mandatory" << endl;
      usage();
    }
  }
  int max_output = numeric_limits<int>::max();
  if (contains(argc, argv, "-m")) {
    int new_max_output = atoi(get_option(argc, argv, "-m").c_str());
    if (new_max_output != -1) {
      max_output = new_max_output; // limit set by user
    }
  } else {
    max_output = 3; // default output limit
  }
  if (not single._brute_mode and not single._smart_mode and
      not paranoiac_mode     and not (histo or sample)) {
    cout << "error: -b, -s or --parano is mandatory" << endl;
    return usage();
  }
  if (single._brute_mode and single._smart_mode) {
    cout << "error: -b XOR -s" << endl;
    return usage();
  }
  float clustering_distance = atof(distance.c_str());
  if (paranoiac_mode) {
    DistMatrix brute_matrix(input_file.c_str(), clustering_distance, start);
    brute_matrix.brute_init();
    DistMatrix smart_matrix(input_file.c_str(), clustering_distance, start);
    smart_matrix.smart_init();
    bool alright = true;

    for (size_t i = 0 ; i < brute_matrix.get_nb_PDBs() ; ++i) {
      for (size_t j = i + 1 ; j < brute_matrix.get_nb_PDBs() ; ++j) {

        float real_dist = brute_matrix.get(i, j)._mini;
        DistRange range = smart_matrix.get(i, j);
        float mini = range._mini;
        float maxi = range._maxi;
        float diff;

        if (mini > real_dist) {
          diff = mini - real_dist;
          cout << "mini > real by: " << diff;
          if (diff < 0.001) {
            cout << " SAFE TO IGNORE: so small";
          }
          cout << '\n';
          alright = false;
        }
        if (maxi < real_dist) {
          diff = real_dist - maxi;
          cout << "maxi < real by: " << diff;
          if (diff < 0.001) {
            cout << " SAFE TO IGNORE: so small";
          }
          cout << '\n';
          alright = false;
        }
        //cout << mini << ' ' << real_dist << ' ' << maxi << endl;
      }
    }
    if (alright) {
      cout << "Distance and range matrices are coherent\n";
    }
    return 0;
  }

  DistMatrix dm(input_file.c_str(), clustering_distance, start);

  if (auto_threshold) {

    string cutoff_str = get_option(argc, argv, "--semi-auto");
    float cutoff      = atof(cutoff_str.c_str());
    float threshold, median;

    cout << "using semi automatic threshold finding (" << cutoff
         << ')' << endl;
    threshold = dm.find_clustering_distance(cutoff, 0.0, &median);
    cout << "threshold found: " << threshold << endl;
    dm.set_clustering_distance(threshold);
  }
  if (sample) {


    string nb_samples_str = get_option(argc, argv, "--sample");
    int nb_samples        = atoi(nb_samples_str.c_str());

    dm.sample_and_output(nb_samples);
    return 0;
  }
  // will output the biggest clusters, in the order they are found
  // if two clusters have same size, only the first one is taken into account
  if (single._brute_mode) {
    dm.brute_init();
  } else if (single._smart_mode) {
    dm.smart_init();
  }
  if (matrix_csv_dump or matrix_txt_dump) {
    if (matrix_csv_dump) {
      dm.csv_dump();
    } else {
      dm.txt_dump();
    }
    return 0;
  }

  ofstream out (output_file.c_str());
  if (not out.is_open()) {
    cout << "error: cannot write to: " << output_file << endl;
    usage();
  }
  // find and print all clusters
  vector<int> remaining_pdbs = dm.get_all_PDBs_index();
  for (int printed = 0 ;
       remaining_pdbs.size() > 0 and printed < max_output ;
       ++printed) {

    vector<int> biggest_cluster;
    vector< vector<int> > pole_position_clusters;

    // find biggest
    dm.get_biggest_cluster(remaining_pdbs,
                           biggest_cluster,
                           pole_position_clusters);
    // output and remove from remaining
    output_cluster(out, biggest_cluster, remaining_pdbs,
                   pole_position_clusters, dm, single._verbose, stable,
                   single._compute_D);
  }
  out.close();
  if (single._output_entropy) {
    fclose(single._entropy_file);
  }

  return 0;
}
