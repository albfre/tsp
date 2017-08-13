#include <assert.h>
#include <cmath>
#include "Distances.h"

using namespace std;
using namespace TravelingSalespersonProblemSolver;

namespace {
  inline double distanceRound( double d ) { return double( long( d + 0.5 ) ); }
}

// VDistances
VDistances::VDistances( const vector< vector< double > >& points,
                        function< double ( double ) > rounding ) :
  rounding_( rounding )
{
  for ( size_t i = 1; i < points.size(); ++i ) {
    assert( points[ i ].size() == points[ 0 ].size() );
  }
}

VDistances::~VDistances() {}

bool VDistances::empty() const { return size() == 0; }

inline double VDistances::computeDistance_( const double* point1, const double* point2, size_t pointDimension ) const
{
  double dist = 0.0;
  for ( size_t i = 0; i < pointDimension; ++i ) {
    const double diff = point1[ i ] - point2[ i ];
    dist += diff * diff;
  }
  return rounding_( sqrt( dist ) );
}

inline double VDistances::computeDistance_( const vector< double >& point1, const vector< double >& point2 ) const
{
  return computeDistance_( &point1[ 0 ], &point2[ 0 ], point1.size() );
}


// MatrixDistances
MatrixDistances::MatrixDistances( const vector< vector< double > >& points,
                                  function< double ( double ) > rounding ) :
  VDistances( points, rounding ),
  size_( points.size() ),
  distances_( size_ * size_, 0.0 )
{
  for ( size_t i = 0; i < size_; ++i ) {
    const size_t iSize = i * size_;
    for ( size_t j = 0; j < size_; ++j ) {
      distances_[ iSize + j ] = computeDistance_( points[ i ], points[ j ] );
    }
  }
}
MatrixDistances::~MatrixDistances() {}
inline double MatrixDistances::operator()( size_t i, size_t j ) const { return distances_[ i * size_ + j ]; }
inline size_t MatrixDistances::size() const { return size_; }
void MatrixDistances::setMatrix( vector< vector< double > >& distances )
{
  for ( size_t i = 0; i < distances.size(); ++i ) {
    assert( distances[ i ].size() == distances.size() );
    assert( distances[ i ][ i ] == 0.0 );
    for ( size_t j = i + 1; j < distances[ i ].size(); ++j ) {
      assert( fabs( distances[ i ][ j ] - distances[ j ][ i ] ) < 1e-9 );
      assert( distances[ i ][ j ] >= 0.0 );
    }
  }
  size_ = distances.size();
  distances_.resize( size_ * size_ );
  for ( size_t i = 0; i < size_; ++i ) {
    const size_t iSize = i * size_;
    for ( size_t j = 0; j < size_; ++j ) {
      distances_[ iSize + j ] = distances[ i ][ j ];
    }
  }
}


// MatrixRoundedDistances
MatrixRoundedDistances::MatrixRoundedDistances( const vector< vector< double > >& points ) :
  MatrixDistances( points, [] ( double d ) { return distanceRound( d ); } )
{}
MatrixRoundedDistances::~MatrixRoundedDistances() {}


// OnTheFlyDistances
OnTheFlyDistances::OnTheFlyDistances( const vector< vector< double > >& points,
                                      function< double ( double ) > rounding ) :
  VDistances( points, rounding ),
  size_( points.size() ),
  pointDimension_( points.empty() ? 0 : points[ 0 ].size() ),
  points_( size_ * pointDimension_, 0.0 )
{
  for ( size_t i = 0; i < points.size(); ++i ) {
    assert( points[ i ].size() == pointDimension_ );
    for ( size_t j = 0; j < pointDimension_; ++j ) {
      points_[ i * pointDimension_ + j ] = points[ i ][ j ];
    }
  }
}
OnTheFlyDistances::~OnTheFlyDistances() {}
inline double OnTheFlyDistances::operator()( size_t i, size_t j ) const { return computeDistance_( &points_[ i * pointDimension_ ], &points_[ j * pointDimension_ ], pointDimension_ ); }
inline size_t OnTheFlyDistances::size() const { return size_; }


// OnTheFlyRoundedDistances
OnTheFlyRoundedDistances::OnTheFlyRoundedDistances( const vector< vector< double > >& points ) :
  OnTheFlyDistances( points, [] ( double d ) { return distanceRound( d ); } )
{}

OnTheFlyRoundedDistances::~OnTheFlyRoundedDistances() {}
