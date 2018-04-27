

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>

#include "PairPurityFitter.h"


#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	Logger::setGlobalLogLevel( "none" );
	TaskFactory::registerTaskRunner<PairPurityFitter>( "PairPurityFitter" );
	TaskEngine engine( argc, argv, "PairPurityFitter" );
	return 0;
}
