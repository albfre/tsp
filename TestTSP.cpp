#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <iterator>
#include <fstream>
#include <vector>

#include "TravelingSalesmanSolver.h"

using namespace std;
using namespace TravelingSalesmanSolver;

double computeDistance( const vector< double >& point1, const vector< double >& point2 )
{
  assert( point1.size() == point2.size() );
  double dist = 0.0;
  for ( size_t i = 0; i < point1.size(); ++i ) {
    double diff = point1[ i ] - point2[ i ];
    dist += diff * diff;
  }
  return sqrt( dist );
}

void testTSP( size_t numOfPoints )
{
  vector< vector< double > > points( numOfPoints, vector< double >( 2 ) );
  for ( size_t i = 0; i < numOfPoints; ++i ) {
    points[ i ][ 0 ] = double( rand() ) / RAND_MAX;
    points[ i ][ 1 ] = double( rand() ) / RAND_MAX;
  }
  vector< vector< double > > distances( numOfPoints, vector< double >( numOfPoints ) );
  for ( size_t i = 0; i < numOfPoints; ++i ) {
    for ( size_t j = i + 1; j < numOfPoints; ++j ) {
      distances[ i ][ j ] = computeDistance( points[ i ], points[ j ] );
      distances[ j ][ i ] = distances[ i ][ j ];
    }
  }

  cout << "Running test with " << numOfPoints << " points." << endl;
  double start( clock() );
  vector< size_t > path = computePath( distances );
  cout << "CPU seconds to run test: " << setprecision( 4 ) << ( clock() - start ) / CLOCKS_PER_SEC << endl;;
}

int main( int argc, const char* argv[] )
{
  bool standardTest = argc == 1;

  if ( standardTest ) {
    testTSP( 50 );
  }
  else {
    size_t numOfPoints = atoi( argv[ 1 ] );
    testTSP( numOfPoints );
  }
}
