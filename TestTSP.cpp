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

#include "drill.h"
#include "TravelingSalespersonProblemSolver.h"

using namespace std;
using namespace TravelingSalespersonProblemSolver;

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

void testTSPRegular( size_t numOfPoints )
{
  size_t numOfPointsBefore = numOfPoints;
  for ( size_t iii = 0; iii < 1; ++iii ) {
    numOfPoints = numOfPointsBefore;
  vector< vector< double > > points;
  size_t rows = size_t( sqrt( double( numOfPoints ) ) );
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < rows; ++j ) {
      points.push_back( vector< double >( 2 ) );
      points.back()[ 0 ] = 0.2 * ( j + 0.5 * ( i % 2 ) );
      points.back()[ 1 ] = sqrt( 0.2 * 0.2 - 0.1 * 0.1 ) * i;
    }
  }
  for ( size_t i = 0; i < numOfPoints / 2; ++i ) {
    size_t ind = rand() % points.size();
    points.erase( points.begin() + ind );
  }
  numOfPoints = points.size();
  cerr << "numPoints: " << numOfPoints << endl;


  vector< vector< double > > distances( numOfPoints + 1, vector< double >( numOfPoints + 1, 1e5 ) );
  for ( size_t i = 0; i < numOfPoints; ++i ) {
    for ( size_t j = i + 1; j < numOfPoints; ++j ) {
      distances[ i ][ j ] = computeDistance( points[ i ], points[ j ] );
      distances[ j ][ i ] = distances[ i ][ j ];
    }
    distances[ i ][ i ] = 0.0;
  }
  distances.back().back() = 0.0;
  distances.back()[ 0 ] = 1e-6;
  distances[ 0 ].back() = 1e-6;

  cout << "Running test on regular instance with " << numOfPoints << " points." << endl;
  double start( clock() );
  vector< size_t > path = computeTour( distances );
  size_t ii = find( path.begin(), path.end(), numOfPoints ) - path.begin();
  cerr << "First point: " << path[ ( ii + 1 ) % path.size() ] << " " << path[ ( ii + path.size() - 1 ) % path.size() ] << endl;
  cout << "CPU seconds to run test: " << setprecision( 4 ) << ( clock() - start ) / CLOCKS_PER_SEC << endl;;
  }
}

void testTSPDrill( int arg )
{
  vector< vector< double > > points = getDrill();
  size_t numOfPoints = points.size();
  vector< vector< double > > distances( numOfPoints, vector< double >( numOfPoints ) );
  for ( size_t i = 0; i < numOfPoints; ++i ) {
    for ( size_t j = i + 1; j < numOfPoints; ++j ) {
      distances[ i ][ j ] = computeDistance( points[ i ], points[ j ] );
      distances[ j ][ i ] = distances[ i ][ j ];
    }
    distances[ i ][ i ] = 0.0;
  }

  cout << "Running test on random instance with " << numOfPoints << " points." << endl;
  double start( clock() );
  vector< size_t > path = computeTour( distances );
  cout << "CPU seconds to run test: " << setprecision( 4 ) << ( clock() - start ) / CLOCKS_PER_SEC << endl;;
}

void testTSPRandom( size_t numOfPoints )
{
  for ( size_t i = 0; i < numOfPoints; ++i ) {
    rand();
  }
  for ( size_t iii = 0; iii < 1; ++iii ) {
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
    distances[ i ][ i ] = 0.0;
  }

  cout << "Running test on random instance with " << numOfPoints << " points." << endl;
  double start( clock() );
  vector< size_t > path = computeTour( distances );
  cout << "CPU seconds to run test: " << setprecision( 4 ) << ( clock() - start ) / CLOCKS_PER_SEC << endl;;
  }
}

int main( int argc, const char* argv[] )
{
  srand( 1729 );
  bool standardTest = argc == 1;
  if ( standardTest ) {
    testTSPRandom( 50 );
  }
  else if ( argc == 2 ) {
    size_t numOfPoints = atoi( argv[ 1 ] );
    testTSPRandom( numOfPoints );
  }
  else {
    if ( atoi( argv[ 1 ] ) == 1 ) {
      testTSPDrill( atoi( argv[ 2 ] ) );
    }
    else {
      size_t numOfPoints = atoi( argv[ 2 ] );
      testTSPRegular( numOfPoints );
    }
  }
}

