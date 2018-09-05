#ifndef PERFORM_HELSGAUN_MOVE_H
#define PERFORM_HELSGAUN_MOVE_H

#include "HelperFunctions.h"

// For a description of the following, see K. Helsgaun, General k-opt submoves for the Lin-Kernighan TSP heuristic.  Mathematical Programming Computation, 2009
namespace TravelingSalespersonProblemSolver {
  // Sets up the vectors p, q, and incl needed to perform moves
  void getPQI_( std::vector< size_t >& p,
                std::vector< size_t >& q,
                std::vector< size_t >& incl,
                const std::vector< size_t >& ts,
                const std::vector< size_t >& tour,
                const std::vector< size_t >& position )
  {
    assert( !ts.empty() );
    assert( ts.size() % 2 == 0 );
    std::vector< size_t > pHalf;
    pHalf.reserve( ts.size() / 2 );
    for ( size_t i = 0; i < ts.size(); i += 2 ) {
      pHalf.push_back( ts[ i ] == next_( ts[ i + 1 ], tour, position ) ? i + 1 : i );
    }
    std::sort( pHalf.begin(), pHalf.end(), [&] ( const size_t i, const size_t j ) {
      return position[ ts[ i ] ] < position[ ts[ j ] ] || ( position[ ts[ i ] ] == position[ ts[ j ] ] && i < j );
    } );

    p.clear();
    p.reserve( ts.size() );
    for ( const auto ph : pHalf ) {
      p.push_back( ph );
      p.push_back( ph % 2 == 0 ? ph + 1 : ph - 1 );
    }

    q = std::vector< size_t >( p.size() );
    for ( size_t i = 0; i < p.size(); ++i ) {
      q[ p[ i ] ] = i;
    }

    incl = std::vector< size_t >( p.size() );
    incl.front() = ts.size() - 1;
    incl.back() = 0;
    for ( size_t i = 1; i + 1 < ts.size(); i += 2 ) {
      incl[ i ] = i + 1;
      incl[ i + 1 ] = i;
    }
  }

  // Perform the k-opt move defined by the edge switches ts using the helper vectors incl, p, and q
  void performKOptMove_( const std::vector< size_t >& ts,
                         const std::vector< size_t >& incl,
                         const std::vector< size_t >& p,
                         const std::vector< size_t >& q,
                         std::vector< size_t >& tour,
                         std::vector< size_t >& position )
  {
    const std::vector< size_t > tourCopy( tour );
    size_t index = 0;
    size_t currentIndex = 0;
    for ( size_t step = 0; step < ts.size() / 2; ++step ) {
      currentIndex = incl[ currentIndex ];
      size_t currentNode = ts[ currentIndex ];
      size_t i = q[ currentIndex ];
      currentIndex = i % 2 == 0 ? p[ i == 0 ? p.size() - 1 : i - 1 ]
                                : p[ i + 1 == p.size() ? 0 : i + 1 ];
      size_t nextNode = ts[ currentIndex ];
      bool increasing = i % 2 == 1;
      if ( increasing ) {
        for ( ; currentNode != nextNode; currentNode = next_( currentNode, tourCopy, position ), ++index ) {
          tour[ index ] = currentNode;
        }
      }
      else {
        for ( ; currentNode != nextNode; currentNode = previous_( currentNode, tourCopy, position ), ++index ) {
          tour[ index ] = currentNode;
        }
      }
      assert( index < tour.size() );
      tour[ index ] = nextNode;
      ++index;
    }
    updatePosition( tour, position );
  }

  void performKOptMove_( const std::vector< size_t >& ts,
                         const std::vector< size_t >& incl,
                         std::vector< size_t >& tour,
                         std::vector< size_t >& position )
  {
    std::vector< size_t > p, q, inclTmp;
    getPQI_( p, q, inclTmp, ts, tour, position );
    performKOptMove_( ts, incl, p, q, tour, position );
  }

  void performKOptMove_( const std::vector< size_t >& ts,
                         std::vector< size_t >& tour,
                         std::vector< size_t >& position )
  {
    std::vector< size_t > p, q, incl;
    getPQI_( p, q, incl, ts, tour, position );
    performKOptMove_( ts, incl, p, q, tour, position );
  }

  // Indicates whether the path switches in ts define a legal tour
  bool makesTour_( const std::vector< size_t >& ts,
                   const std::vector< size_t >& tour,
                   const std::vector< size_t >& position )
  {
    std::vector< size_t > p, q, incl;
    getPQI_( p, q, incl, ts, tour, position );

    // For a description of the following no-body loop, see K. Helsgaun, General k-opt submoves for the Lin-Kernighan TSP heuristic.  Mathematical Programming Computation, 2009
    size_t count = 1;
    for ( size_t i = ts.size(); ( i = ( q[ incl[ p[ i - 1 ] ] ] + 1 ) ^ 1 ) != 0; ++count );
    return 2 * count == ts.size();
  }
}

#endif // PERFORM_HELSGAUN_MOVE_H
