all:    tags
	scons -j5 opt=true

dryrun:
	scons -j5 opt=true -n

fast:
	scons -j5

tags:
	find . -name "*.cc" >  src_files
	find . -name "*.h"  >> src_files
	cat src_files | xargs etags -a -o src/TAGS
	rm -f src_files

test: all
	./test.sh

debug:
	scons -j5 debug=true

clean:
	scons -c
	rm -rf current.out TAGS build very_few_pdbs_rmsd_current out
