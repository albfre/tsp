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

  MatrixDistances distanceMat( points );
  distanceMat.setMatrix( distances );

  cout << "Running test on regular instance with " << numOfPoints << " points." << endl;
  double start( clock() );
  vector< size_t > path = computeTour( distanceMat );
  size_t ii = find( path.begin(), path.end(), numOfPoints ) - path.begin();
  cerr << "First point: " << path[ ( ii + 1 ) % path.size() ] << " " << path[ ( ii + path.size() - 1 ) % path.size() ] << endl;
  cout << "CPU seconds to run test: " << setprecision( 4 ) << ( clock() - start ) / CLOCKS_PER_SEC << endl;;
  }
}

void trim( string& s )
{
  s.erase( s.find_last_not_of( " \n\r\t" ) + 1 );
}

vector< vector< double > > readTSP( string fileName )
{
  ifstream file( "lib/" + fileName + ".tsp" );
  assert( file.is_open() );
  vector< vector< double > > points;
  points.reserve( 1000 );
  string line;
  while ( !file.eof() ) {
    getline( file, line );
    trim( line );
    if ( line == "NODE_COORD_SECTION" ) {
      getline( file, line );
      trim( line );
      while ( !file.eof() && line != "EOF" ) {
        istringstream iss( line );
        string token;
        assert( iss );
        iss >> token;
        assert( iss );
        iss >> token;
        points.push_back( vector< double >( 2 ) );
        points.back()[ 0 ] = atof( token.c_str() );
        assert( iss );
        iss >> token;
        points.back()[ 1 ] = atof( token.c_str() );
        getline( file, line );
        trim( line );
      }
    }
  }
  assert( points.size() > 0 );
  return points;
}

void testTSPLib( string name )
{
  vector< vector< double > > points = readTSP( name );
  size_t numOfPoints = points.size();
  vector< vector< double > > emptyPoints;
  bool useFullMatrix = numOfPoints < 2000;
  MatrixRoundedDistances mDistances( useFullMatrix ? points : emptyPoints );
  OnTheFlyRoundedDistances sDistances( useFullMatrix ? emptyPoints : points );
  VDistances* distances;
  if ( useFullMatrix ) {
    distances = &mDistances;
  }
  else {
    distances = &sDistances;
  }

  cout << "Running test on " + name + " with " << numOfPoints << " points." << endl;
  double start( clock() );
  vector< size_t > path = computeTour( *distances );
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
  MatrixDistances distances( points );

  cout << "Running test on random instance with " << numOfPoints << " points." << endl;
  double start( clock() );
  vector< size_t > path = computeTour( distances );
  cout << "CPU seconds to run test: " << setprecision( 4 ) << ( clock() - start ) / CLOCKS_PER_SEC << endl;;
  }
}

void testTSPRandomClusters( size_t numOfPoints )
{
  for ( size_t i = 0; i < numOfPoints; ++i ) {
    rand();
  }
  for ( size_t iii = 0; iii < 1; ++iii ) {
  vector< vector< double > > points;
  points.reserve( numOfPoints );
  size_t numOfClusters = min( size_t( ( rand() % 30 ) + 5 ), size_t( numOfPoints / 2 ) );
  for ( size_t c = 0; c < numOfClusters; ++c ) {
    vector< double > cluster( 2 );
    cluster[ 0 ] = 10.0 * ( double( rand() ) / RAND_MAX );
    cluster[ 1 ] = 10.0 * ( double( rand() ) / RAND_MAX );
    for ( size_t i = 0; i * numOfClusters < numOfPoints; ++i ) {
      points.push_back( vector< double >( 2 ) );
      points.back()[ 0 ] = cluster[ 0 ] + double( rand() ) / RAND_MAX;
      points.back()[ 1 ] = cluster[ 1 ] + double( rand() ) / RAND_MAX;
    }
  }
  numOfPoints = points.size();
  MatrixDistances distances( points );

  cout << "Running test on random instance with " << numOfPoints << " points and " << numOfClusters << " clusters." << endl;
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
    switch ( atoi( argv[ 1 ] ) ) {
      case 0: testTSPRegular( atoi( argv[ 2 ] ) ); break;
      case 1: testTSPLib( "pcb1173" ); break;
      case 2: testTSPRandomClusters( atoi( argv[ 2 ] ) ); break;
      default: testTSPLib( string( argv[ 2 ] ) );
    }
  }
}
