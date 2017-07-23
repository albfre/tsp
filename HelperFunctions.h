#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

namespace TravelingSalespersonProblemSolver {
  const double tolerance = 1e-9;

  void updatePosition( const std::vector< size_t >& tour,
                       std::vector< size_t >& position )
  {
    position.resize( tour.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }
  }

  bool inBetterTour( size_t t1, size_t t2, const std::vector< size_t >& betterTour )
  {
    for ( size_t i = 0; i + 1 < betterTour.size(); ++i ) {
      if ( ( betterTour[ i ] == t1 && betterTour[ i + 1 ] == t2 ) ||
           ( betterTour[ i ] == t2 && betterTour[ i + 1 ] == t1 ) ) {
        return true;
      }
    }
    return betterTour.size() > 1 &&
           ( ( betterTour.front() == t1 && betterTour.back() == t2 ) ||
             ( betterTour.front() == t2 && betterTour.back() == t1 ) );
  }

  void printTour_( std::string label, const std::vector< size_t >& tour )
  {
    std::cerr << label << ": ";
    for ( size_t i = 0; i < tour.size(); ++i ) {
      std::cerr << tour[ i ] << " ";
    }
    std::cerr << std::endl;
  }

  double getLength_( const std::vector< size_t >& tour, const VDistances& distances )
  {
    double distance = distances( tour.back(), tour.front() );
    for ( size_t i = 0; i + 1 < tour.size(); ++i ) {
      distance += distances( tour[ i ], tour[ i + 1 ] );
    }
    return distance;
  }

  bool isTour_( const std::vector< size_t >& tour, const std::vector< size_t >& position )
  {
    assert( tour.size() == position.size() );
    for ( size_t i = 0; i < position.size(); ++i ) {
      if ( position[ i ] >= tour.size() || tour[ position[ i ] ] != i ) {
        std::cerr << "position[ i ]: " << position[ i ] << " tour[ position[ i ] ]: " << tour[ position[ i ] ] << " i: " << i << std::endl;
        return false;
      }
    }
    return true;
  }

  bool isTour_( const std::vector< size_t >& tour, const VDistances& distances )
  {
    assert( tour.size() == distances.size() );
    std::vector< size_t > position( tour.size() );
    updatePosition( tour, position );
    return isTour_( tour, position );
  }

  size_t previous_( size_t node, const std::vector< size_t >& tour, const std::vector< size_t >& position )
  {
    return position[ node ] > 0 ? tour[ position[ node ] - 1 ] : tour.back();
  }

  size_t next_( size_t node, const std::vector< size_t >& tour, const std::vector< size_t >& position )
  {
    return position[ node ] + 1 < tour.size() ? tour[ position[ node ] + 1 ] : tour.front();
  }

  bool between_( size_t a, size_t b, size_t c, const std::vector< size_t >& position )
  {
    return position[ a ] <= position[ c ] ? position[ a ] >= position[ b ] || position[ b ] >= position[ c ]
                                          : position[ b ] <= position[ a ] && position[ b ] >= position[ c ];
  }

  double getDistance_( const VDistances& distances,
                       const std::vector< double >& lagrangeMultipliers,
                       size_t i,
                       size_t j )
  {
    return distances( i, j ) + lagrangeMultipliers[ i ] + lagrangeMultipliers[ j ];
  }

}

#endif // HELPER_FUNCTIONS_H
