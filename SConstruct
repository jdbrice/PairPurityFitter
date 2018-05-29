#!python
import subprocess
import os
import platform

# SConscript('color_SConscript')
# Import( 'env' )

SConscript('color_SConscript')
Import( 'env' )

env.Append(CPPPATH    = [ "modules/" ])
env.Append(CXXFLAGS   = ['-std=c++11'])
JDB_LIB = os.environ.get("JDB_LIB", "" )
env.Append(LIBPATH    = [ JDB_LIB + "/lib/" ] )
env.Append(LIBS       = [ "libRooBarbCore.a", "libRooBarbConfig.a", "libRooBarbTasks.a", "libRooBarbRootAna.a", "libRooBarbUnitTest.a", "libRooBarbExtra.a" ] )

debug = ARGUMENTS.get('debug', 0)
if int(debug):
	env.Append(CXXFLAGS = '-g')

# # ROOT
ROOTCFLAGS      = subprocess.check_output( ['root-config',  '--cflags'] ).rstrip().split( " " )
ROOTLIBS        = subprocess.check_output( ["root-config",  "--libs"] )
ROOT_SYS        = os.environ.get( "ROOTSYS", "" )


env.Append(CPPPATH   = [ ROOT_SYS + "/include/", JDB_LIB + "/include" ])
env.Append(CXXFLAGS  = ROOTCFLAGS )
env.Append(LINKFLAGS = ROOTCFLAGS )

env.Append(CXXFLAGS 		= "-DJDB_LOG_LEVEL=0" )

env[ "_LIBFLAGS" ] = env[ "_LIBFLAGS" ] + " " + ROOTLIBS + " "
# REMOVE "-pthread" from clang link options in OS X
# assuming darwin=clang but must be a better way...
if "Darwin" in platform.platform() :
	env[ "LINKFLAGS" ].remove( "-pthread" )

env.Program( target="bin/purity.app", source=[ "Engine.cpp"] )
