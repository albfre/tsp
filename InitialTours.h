#ifndef INITIAL_TOURS_H
#define INITIAL_TOURS_H

#include "HelperFunctions.h"

namespace TravelingSalespersonProblemSolver {

  std::vector< size_t > getRandomTour_( const VDistances& distances, std::function<int()> random )
  {
    std::vector< size_t > tour( distances.size() );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      tour[ i ] = i;
    }

    for ( size_t i = 0; i < distances.size(); ++i ) {
      size_t ind1 = random() % distances.size();
      size_t ind2 = random() % distances.size();
      std::swap( tour[ ind1 ], tour[ ind2 ] );
    }
    return tour;
  }

  std::vector< size_t > getNearestNeighborTour_( const VDistances& distances, std::function<int()> random  )
  {
    std::vector< size_t > tour;
    tour.reserve( distances.size() );
    size_t startNode = random() % distances.size();
    tour.push_back( startNode );
    std::vector< bool > usedNodes( distances.size(), false );
    usedNodes[  startNode ] = true;
    while ( tour.size() < distances.size() ) {
      size_t currentNode = tour.back();
      size_t minUnusedIndex = static_cast< size_t >( -1 );
      double minUnusedDistance = std::numeric_limits< double >::max();
      for ( size_t i = 0; i < distances.size(); ++i ) {
        if ( !usedNodes[ i ] ) {
          double distance = distances( currentNode, i );
          if ( distance < minUnusedDistance ) {
            minUnusedIndex = i;
            minUnusedDistance = distance;
          }
        }
      }
      assert( minUnusedIndex != static_cast< size_t >( -1 ) );
      tour.push_back( minUnusedIndex );
      usedNodes[ minUnusedIndex ] = true;
    }
    assert( isTour_( tour, distances ) );
    return tour;
  }

  std::vector< size_t > getGreedyTour_( const VDistances& distances,
                                                    const std::vector< double >& lagrangeMultipliers )
  {
    // The greedy heuristic of matroid theory adds the shortest edge that neither
    // makes the degree of a vertex greater than 2 nor creates a cycle shorter than N.
    assert( distances.size() == lagrangeMultipliers.size() );
    std::vector< size_t > degree( distances.size() );

    std::vector< size_t > fragmentIndices( distances.size() );
    iota( fragmentIndices.begin(), fragmentIndices.end(), 0 );
    std::vector< std::vector< size_t > > fragments( distances.size() );
    for ( const auto& i : fragmentIndices ) {
      fragments[ i ] = std::vector< size_t >( 1, i  );
    }
    std::vector< std::pair< double, std::pair< size_t, size_t > > > edges;
    edges.reserve( distances.size() * ( distances.size() - 1 ) / 2 );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      for ( size_t j = i + 1; j < distances.size(); ++j ) {
        edges.push_back( std::make_pair( getDistance_( distances, lagrangeMultipliers, i, j ), std::make_pair( i, j ) ) );
      }
    }
    std::sort( edges.rbegin(), edges.rend() ); // sort from rbegin to rend => smallest element last
    while ( !edges.empty() ) {
      auto edgeNodes = edges.back().second;
      edges.pop_back();
      const size_t v1 = edgeNodes.first;
      const size_t v2 = edgeNodes.second;
      if ( degree[ v1 ] < 2 && degree[ v2 ] < 2 && fragmentIndices[ v1 ] != fragmentIndices[ v2 ] ) {
        assert( fragmentIndices[ v1 ] != static_cast< size_t >( -1 ) );
        assert( fragmentIndices[ v2 ] != static_cast< size_t >( -1 ) );
        std::vector< size_t >& f1 = fragments[ fragmentIndices[ v1 ] ];
        std::vector< size_t >& f2 = fragments[ fragmentIndices[ v2 ] ];
        assert( f1.front() == v1 || f1.back() == v1 );
        assert( f2.front() == v2 || f2.back() == v2 );

        if ( f1.front() == v1 ) {
          reverse( f1.begin(), f1.end() );
        }
        if ( f2.back() == v2 ) {
          reverse( f2.begin(), f2.end() );
        }
        f1.insert( f1.end(), f2.begin(), f2.end() );
        fragmentIndices[ f2.back() ] = fragmentIndices[ v1 ];
        f2.clear();
        ++degree[ v1 ];
        ++degree[ v2 ];
        if ( degree[ v1 ] >= 2 ) {
          fragmentIndices[ v1 ] = static_cast< size_t >( -1 );
        }
        if ( degree[ v2 ] >= 2 ) {
          fragmentIndices[ v2 ] = static_cast< size_t >( -1 );
        }
      }
    }
    std::vector< size_t > tour;
    for ( size_t i = 0; i < fragments.size(); ++i ) {
      if ( !fragments[ i ].empty() ) {
        tour.swap( fragments[ i ] );
        break;
      }
    }

    assert( isTour_( tour, distances ) );
    return tour;
  }

  std::vector< size_t > getGreedyTour_( const VDistances& distances )
  {
    std::vector< double > lagrangeMultipliers( distances.size(), 0.0 );
    return getGreedyTour_( distances, lagrangeMultipliers );
  }

  std::vector< size_t > getHelsgaunInitialTour_( const std::vector< std::vector< size_t > >& nearestNeighbors,
                                                 const std::vector< std::vector< size_t > >& helsgaunNeighbors,
                                                 const std::vector< std::vector< double > >& helsgaunDistances,
                                                 const std::vector< size_t >& bestTour,
                                                 std::function<int()> random )
  {
    assert( helsgaunNeighbors.size() == helsgaunDistances.size() );
    assert( bestTour.size() == nearestNeighbors.size() );
    std::vector< size_t > tour( 1, random() % nearestNeighbors.size() );
    size_t current = tour.front();
    std::vector< bool > added( nearestNeighbors.size(), false );
    added[ current ] = true;

    std::vector< size_t > bestPosition( bestTour.size() );
    for ( size_t i = 0; i < bestTour.size(); ++i ) {
      bestPosition[ bestTour[ i ] ] = i;
    }

    while ( tour.size() < bestTour.size() ) {
      current = tour.back();
      bool found = false;
      size_t bestNext = static_cast< size_t >( -1 );
      for ( size_t i = 0; i < helsgaunNeighbors[ current ].size(); ++i ) {
        size_t next = helsgaunNeighbors[ current ][ i ];
        if ( !added[ next ] &&
             helsgaunDistances[ current ][ i ] < tolerance &&
             ( bestTour[ ( bestPosition[ current ] + 1 ) % bestTour.size() ] == next ||
               bestTour[ ( bestPosition[ current ] + bestTour.size() - 1 ) % bestTour.size() ] == next ) ) {
          bestNext = next;
          added[ next ] = true;
          found = true;
          break;
        }
      }
      if ( !found ) {
        for ( size_t i = 0; i < nearestNeighbors[ current ].size(); ++i ) {
          size_t next = nearestNeighbors[ current ][ i ];
          if ( !added[ next ] ) {
            bestNext = next;
            added[ next ] = true;
            found = true;
            break;
          }
        }
      }
      if ( !found ) {
        for ( size_t i = 0; i < added.size(); ++i ) {
          if ( !added[ i ] ) {
            bestNext = i;
            added[ i ] = true;
            break;
          }
        }
      }
      assert( bestNext != static_cast< size_t >( -1 ) );
      tour.push_back( bestNext );
    }

    std::vector< size_t > position( tour.size() );
    updatePosition( tour, position );
    assert( isTour_( tour, position ) );
    return tour;
  }

  // Computes the optimal tour using brute force
  std::vector< size_t > getBruteForceTour_( const VDistances& distances )
  {
    assert( !distances.empty() );
    assert( distances.size() < 10 );
    // It would be possible to speed up the brute force by only considering tours
    // of one orientation, but since we use it only for small instances, this improvement
    // is unnecessary.
    std::vector< size_t > tour( distances.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      tour[ i ] = i;
    }
    std::vector< size_t > bestTour( tour );
    double bestDistance = getLength_( tour, distances );
    do {
      double distance = getLength_( tour, distances );
      if ( distance < bestDistance ) {
        bestDistance = distance;
        bestTour = tour;
      }
     } while ( std::next_permutation( tour.begin() + 1, tour.end() ) ); // tour.begin() + 1 => fixed first node
    return bestTour;
  }

}

#endif // INITIAL_TOURS_H
