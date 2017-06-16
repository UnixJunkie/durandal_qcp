#!/bin/bash

set -x

# tests using KABSCH
rm -f out current.out
./durandal --kabsch -i ./very_few_pdbs -d 2.0 --brute -o out -m -1 > /dev/null
egrep "^members|^pole" out > current.out
diff reference.out current.out
./durandal --kabsch -i ./very_few_pdbs -d 2.0 --smart -o out -m -1 > /dev/null
egrep "^members|^pole" out > current.out
diff reference.out current.out
./durandal --kabsch -i ./very_few_pdbs -d 2.0 --parano -o out -m -1 | \
egrep -v "SAFE TO IGNORE"
# tests using THEOBALD
rm -f out current.out
./durandal -i ./very_few_pdbs -d 2.0 --brute -o out -m -1 > \
/dev/null
egrep "^members|^pole" out > current.out
diff reference.out current.out
./durandal -i ./very_few_pdbs -d 2.0 --smart -o out -m -1 > \
/dev/null
egrep "^members|^pole" out > current.out
diff reference.out current.out
./durandal -i ./very_few_pdbs -d 2.0 --parano -o out -m -1 | \
egrep -v "SAFE TO IGNORE"
# test clusters' energy characteristics
./durandal -i ./very_few_pdbs_e -d 2.0 -s -o out -m -1 -v
diff out cluster_energies_reference.out
# test ranker
rm -f very_few_pdbs_rmsd_current
./ranker ./very_few_pdbs > very_few_pdbs_rmsd_current
diff very_few_pdbs_rmsd_current very_few_pdbs_rmsd_reference
rm -f very_few_pdbs_rmsd_current
./ranker ./very_few_pdbs --theobald > very_few_pdbs_rmsd_current
diff very_few_pdbs_rmsd_current very_few_pdbs_rmsd_reference
