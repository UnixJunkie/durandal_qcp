import os, sys, commands

debug_mode   = True
profile_mode = False

for key, value in ARGLIST:
    if key == 'help':
        print "# production compile:"
        print "scons -j5"
        sys.exit(0)
    if key == 'debug':
        debug_mode = True
    if key == 'opt':
        debug_mode = False
    if key == 'prof':
        profile_mode = True

parse_flags_content = (" -D_USE_FAST_RMSD_ -D_SHOW_PERCENTAGE_COMPLETE_ " +
                       "-D_LARGE_DECOY_SET_ ")

link_flags = " "
common_flags = " -W -Wall "
if debug_mode:
    parse_flags_content = " -g -O0 " + common_flags
if not debug_mode:
    parse_flags_content = " -O3 -DNDEBUG " + common_flags
if profile_mode:
    parse_flags_content = " -pg -g -O3 -DNDEBUG " + common_flags
    link_flags = " -pg "

# dont' pollute current dir with object files
VariantDir('build', 'src', duplicate=0)

env = Environment (ENV = {'PATH' : os.environ['PATH'], # used by colorgcc
                          'TERM' : os.environ['TERM'],
                          'HOME' : os.environ['HOME']},
                   CXX = "g++",
                   CCFLAGS = parse_flags_content,
                   LINKFLAGS = link_flags)

durandal_prog = env.Program(
    'durandal',
    ['build/durandal.cc','build/DistMatrix.cc','build/DistRange.cc',
     'build/Stru.cc','build/SimpPDB.cc','build/rmsd.cc','build/qcprot.cc',
     'build/Singleton.cc'])

ranker_prog = env.Program('ranker', ['build/ranker.cc','build/Stru.cc',
                                     'build/SimpPDB.cc','build/rmsd.cc',
                                     'build/qcprot.cc','build/Singleton.cc'])

carmsd_prog = env.Program('carmsd', ['build/carmsd.cc','build/Stru.cc',
                                     'build/SimpPDB.cc','build/rmsd.cc',
                                     'build/qcprot.cc'])
