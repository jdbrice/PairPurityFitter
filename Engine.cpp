

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>

#include "PairPurityFitter.h"
#include "RatioMaker.h"
#include "CutOptimizationMaker.h"

#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	Logger::setGlobalLogLevel( "none" );
	TaskFactory::registerTaskRunner<PairPurityFitter>( "PairPurityFitter" );
	TaskFactory::registerTaskRunner<RatioMaker>( "RatioMaker" );
	TaskFactory::registerTaskRunner<CutOptimizationMaker>( "CutOptimizationMaker" );

	TaskEngine engine( argc, argv, "PairPurityFitter" );
	return 0;
}
