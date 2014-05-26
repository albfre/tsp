/* SYSTEM INCLUDES */
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <limits>
#include <set>

/* HEADER */
#include "TravelingSalespersonProblemSolver.h"

using namespace std;

namespace {
  void assertIsTour_( const vector< size_t >& tour,
                      const vector< vector< double > >& distances )
  {
    assert( tour.size() == distances.size() );
    vector< size_t > tourCopy( tour );
    sort( tourCopy.begin(), tourCopy.end() );
    for ( size_t i = 0; i < tourCopy.size(); ++i ) {
      if ( tourCopy[ i ] != i ) {
        cerr << "sort( tourCopy )[ i ]: " << tourCopy[ i ] << ", i: " << i << endl;
      }
      assert( tourCopy[ i ] == i );
    }
  }

  vector< size_t > getRandomtour_( const vector< vector< double > >& distances )
  {
    vector< size_t > tour( distances.size() );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      tour[ i ] = i;
    }

    for ( size_t i = 0; i < distances.size(); ++i ) {
      size_t ind1 = rand() % distances.size();
      size_t ind2 = rand() % distances.size();
      swap( tour[ ind1 ], tour[ ind2 ] );
    }
    return tour;
  }

  /*
  void crossingReverse( vector< size_t >& tour, const size_t i, const size_t j )
  {
    assert( i > j );
    vector< size_t >::iterator iIt = tour.begin() + i;
    vector< size_t >::iterator jIt = tour.begin() + j;
    while ( iIt != tour.end() && jIt != tour.begin() ) {
      --jIt;
      iter_swap( iIt, jIt );
      ++iIt;
    }
    if ( iIt == tour.end() ) {
      reverse( tour.begin(), jIt );
    }
    else {
      reverse( iIt, tour.end() );
    }
  }

  void reverseShortest( vector< size_t >& tour, const size_t i, const size_t j )
  {
    if ( j - i < tour.size() + i - j ) {
      reverse( tour.begin() + i, tour.begin() + j );
    }
    else {
      vector< size_t >::iterator iIt = tour.begin() + i + 1;
      vector< size_t >::iterator jIt = tour.begin() + j;
      while ( iIt-- != tour.begin() && jIt != tour.end() ) {
        iter_swap( iIt, jIt );
        ++jIt;
      }
      if ( jIt == tour.end() ) {
        reverse( tour.begin(), iIt );
      }
      else {
        reverse( jIt, tour.end() );
      }
    }
  }
  */

  vector< size_t > getNearestNeighbortour_( const vector< vector< double > >& distances )
  {
    vector< size_t > tour;
    tour.reserve( distances.size() );
    size_t startNode = 0; //rand() % distances.size();
    tour.push_back( startNode );
    set< size_t > usedNodes;
    usedNodes.insert( startNode );
    while ( tour.size() < distances.size() ) {
      size_t currentNode = tour.back();
      size_t minUnusedIndex = (size_t)-1;
      double minUnusedDistance = numeric_limits< double >::max();
      for ( size_t i = 0; i < distances.size(); ++i ) {
        if ( usedNodes.find( i ) == usedNodes.end() ) {
          double distance = distances[ currentNode ][ i ];
          if ( distance < minUnusedDistance ) {
            minUnusedIndex = i;
            minUnusedDistance = distance;
          }
        }
      }
      assert( minUnusedIndex != (size_t)-1 );
      tour.push_back( minUnusedIndex );
      usedNodes.insert( minUnusedIndex );
    }
    assertIsTour_( tour, distances );
    return tour;
  }

  vector< size_t > getGreedyTour_( const vector< vector< double > >& distances )
  {
    // The greedy heuristic of matroid theory adds the shortest edge that neither
    // makes the degree of a vertex greater than 2 nor creates a cycle shorter than N.
    vector< size_t > degree( distances.size() );

    vector< size_t > fragmentIndices( distances.size() );
    vector< vector< size_t > > fragments( distances.size() );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      fragments[ i ] = vector< size_t >( 1, i );
      fragmentIndices[ i ] = i;
    }
    vector< pair< double, pair< size_t, size_t > > > edges;
    edges.reserve( distances.size() * ( distances.size() - 1 ) / 2 );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      for ( size_t j = i + 1; j < distances[ i ].size(); ++j ) {
        edges.push_back( make_pair( distances[ i ][ j ], make_pair( i, j ) ) );
      }
    }
    vector< pair< size_t, size_t > > addedEdges;
    sort( edges.rbegin(), edges.rend() ); // sort from rbegin to rend => smallest element last
    while ( edges.size() > 0 ) {
      pair< size_t, size_t > edgeNodes = edges.back().second;
      edges.pop_back();
      size_t v1 = edgeNodes.first;
      size_t v2 = edgeNodes.second;
      if ( degree[ v1 ] < 2 && degree[ v2 ] < 2 && fragmentIndices[ v1 ] != fragmentIndices[ v2 ] ) {
        assert( fragmentIndices[ v1 ] != (size_t)-1 );
        assert( fragmentIndices[ v2 ] != (size_t)-1 );
        vector< size_t >& f1 = fragments[ fragmentIndices[ v1 ] ];
        vector< size_t >& f2 = fragments[ fragmentIndices[ v2 ] ];
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
          fragmentIndices[ v1 ] = (size_t)-1;
        }
        if ( degree[ v2 ] >= 2 ) {
          fragmentIndices[ v2 ] = (size_t)-1;
        }
      }
    }
    vector< size_t > tour;
    for ( size_t i = 0; i < fragments.size(); ++i ) {
      if ( fragments[ i ].size() > 0 ) {
        tour.swap( fragments[ i ] );
        break;
      }
    }

    assertIsTour_( tour, distances );
    return tour;
  }

  vector< vector< size_t > > computeNearestNeighbors_( const vector< vector< double > >& distances,
                                                       size_t numberOfNeighbors )
  {
    numberOfNeighbors = min( numberOfNeighbors, distances.size() - 1 );
    vector< vector< size_t > > nearestNeighbors( distances.size(), vector< size_t >( numberOfNeighbors ) );
    vector< pair< double, size_t > > tmpNeighbors( distances.size() );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      for ( size_t j = 0; j < distances[ i ].size(); ++j ) {
        tmpNeighbors[ j ] = make_pair( distances[ i ][ j ], j );
      }
      sort( tmpNeighbors.begin(), tmpNeighbors.end() );
      for ( size_t j = 0; j < numberOfNeighbors; ++j ) {
        // Take neighbor j + 1 in order to avoid adding self as neighbor
        nearestNeighbors[ i ][ j ] = tmpNeighbors[ j + 1 ].second;
      }
    }
    for ( size_t i = 0; i < nearestNeighbors.size(); ++i ) {
      assert( find( nearestNeighbors[ i ].begin(), nearestNeighbors[ i ].end(), i ) == nearestNeighbors[ i ].end() );
    }
    return nearestNeighbors;
  }

  double compute1TreeLength_( vector< size_t >& vertexDegrees,
                              const vector< vector< double > >& distances,
                              const vector< double >& lagrangeMultipliers )
  {
    vertexDegrees.assign( distances.size(), 0 );

    // 1. Compute length of the minimum spanning tree excluding one vertex, using Prim's algorithm.

    // Select the vertex with minimum maximal distances to a neighbor in
    // order to allow for problems with artificial nodes (such as those used
    // to construct Hamiltonian tours) with epsilon distance to all neighbors.
    size_t minimaxVertex = 0;
    double minimaxValue = *max_element( distances[ 0 ].begin(), distances[ 0 ].end() );
    for ( size_t i = 1; i < distances.size(); ++i ) {
      double value = *max_element( distances[ i ].begin(), distances[ i ].end() );
      if ( value < minimaxValue ) {
        minimaxValue = value;
        minimaxVertex = i;
      }
    }

    vector< size_t > unusedVertices;
    unusedVertices.reserve( distances.size() - 1 );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      if ( i != minimaxVertex ) {
        unusedVertices.push_back( i );
      }
    }

    // Popping one element means that it is added to the tree; add the root to the tree.
    size_t rootVertex = unusedVertices.back();
    unusedVertices.pop_back();

    // For each unused vertex i, closestTreeNode[ i ] points to the vertex in the tree which is closest to i.
    vector< size_t > closestTreeNode( distances.size(), rootVertex );

    double treeLength = 0.0;
    while ( unusedVertices.size() > 0 ) {
      size_t fromIndex = 0;
      vector< size_t >::iterator minIt = unusedVertices.begin();
      double minDistance = numeric_limits< double >::max();
      for ( vector< size_t >::iterator it = unusedVertices.begin(); it != unusedVertices.end(); ++it ) {
        size_t mstIndex = closestTreeNode[ *it ];
        double distance = distances[ mstIndex ][ *it ] + lagrangeMultipliers[ mstIndex ] + lagrangeMultipliers[ *it ];
        if ( distance < minDistance ) {
          minDistance = distance;
          fromIndex = mstIndex;
          minIt = it;
        }
      }
      size_t indexOfNewTreeNode = *minIt;
      treeLength += minDistance;
      ++vertexDegrees[ fromIndex ];
      ++vertexDegrees[ indexOfNewTreeNode ];
      *minIt = unusedVertices.back();
      unusedVertices.pop_back();

      for ( vector< size_t >::iterator it = unusedVertices.begin(); it != unusedVertices.end(); ++it ) {
        size_t mstIndex = closestTreeNode[ *it ];
        double oldDistance = distances[ mstIndex ][ *it ] + lagrangeMultipliers[ mstIndex ] + lagrangeMultipliers[ *it ];
        double newDistance = distances[ indexOfNewTreeNode ][ *it ] + lagrangeMultipliers[ indexOfNewTreeNode ] + lagrangeMultipliers[ *it ];
        if ( newDistance < oldDistance ) {
          closestTreeNode[ *it ] = indexOfNewTreeNode;
        }
      }
    }

    // 2. Add the two shortest edges connecting to the first vertex
    double minElement = numeric_limits< double >::max();
    double secondMinElement = numeric_limits< double >::max();
    for ( size_t i = 0; i < distances[ minimaxVertex ].size(); ++i ) {
      double value = distances[ minimaxVertex ][ i ] + lagrangeMultipliers[ minimaxVertex ] + lagrangeMultipliers[ i ];
      if ( value < secondMinElement ) {
        secondMinElement = value;
        if ( value < minElement ) {
          secondMinElement = minElement;
          minElement = value;
        }
      }
    }
    vertexDegrees[ minimaxVertex ] = 2;

    treeLength += minElement + secondMinElement;
    return treeLength - 2.0 * accumulate( lagrangeMultipliers.begin(), lagrangeMultipliers.end(), 0.0 );
  }

  double getHeldKarpLowerBound_( const vector< vector< double > >& distances )
  {
    double bestLength = numeric_limits< double >::min();
    double treeLength = numeric_limits< double >::min();
    vector< double > lagrangeMultipliers( distances.size() );
    double lambda = 0.1;
    for ( size_t i = 0; i < 50; ++i ) {
      vector< size_t > vertexDegrees;
      treeLength = compute1TreeLength_( vertexDegrees, distances, lagrangeMultipliers );
      bestLength = max( bestLength, treeLength );
      double denominator = 0.0;
      for ( size_t j = 0; j < lagrangeMultipliers.size(); ++j ) {
        double d = double( vertexDegrees[ j ] ) - 2.0;
        denominator += d * d;
      }
      if ( denominator == 0.0 ) {
        return bestLength;
      }
      double t = 2.0 * lambda * treeLength / denominator;

      for ( size_t j = 0; j < lagrangeMultipliers.size(); ++j ) {
        lagrangeMultipliers[ j ] += t * ( int( vertexDegrees[ j ] ) - 2 );
      }
      lambda = 1.0 / ( 20.0 + 10 * i );
    }

    return bestLength;
  }

  double getLength_( vector< size_t > tour,
                     const vector< vector< double > >& distances )
  {
    double distance = distances[ tour.back() ][ tour.front() ];
    for ( size_t i = 0; i + 1 < tour.size(); ++i ) {
      distance += distances[ tour[ i ] ][ tour[ i + 1 ] ];
    }
    return distance;
  }

  void update3OptIntervals_( vector< size_t >& tour,
                             size_t i1Begin, size_t i1End, bool reverseI1,
                             size_t i2Begin, size_t i2End, bool reverseI2,
                             size_t i3Begin, size_t i3End, bool reverseI3,
                             const vector< vector< double > >& distances )
  {
    vector< size_t > tourCopy( tour );
    size_t tourIdx = 0;
    for ( size_t i = 0; i < 3; ++i ) {
      size_t iBegin, iEnd;
      bool reverse;
      switch ( i ) {
        case 0: iBegin = i1Begin; iEnd = i1End; reverse = reverseI1; break;
        case 1: iBegin = i2Begin; iEnd = i2End; reverse = reverseI2; break;
        case 2: iBegin = i3Begin; iEnd = i3End; reverse = reverseI3; break;
        default: assert( false ); iBegin = 0; iEnd = 0; reverse = false;
      }
      if ( reverse ) {
        if ( iBegin < iEnd ) {
          iBegin += tour.size();
        }
        for ( size_t idx = iBegin; idx > iEnd; --idx, ++tourIdx ) {
          tour[ tourIdx ] = tourCopy[ idx % tour.size() ];
        }
      }
      else {
        if ( iEnd < iBegin ) {
          iEnd += tour.size();
        }
        for ( size_t idx = iBegin; idx < iEnd; ++idx, ++tourIdx ) {
          tour[ tourIdx ] = tourCopy[ idx % tour.size() ];
        }
      }
    }
    if ( tourIdx != tour.size() ) {
      cerr << tourIdx << " " << tour.size() << endl;
    }
    assert( tourIdx == tour.size() );
    assert( getLength_( tour, distances ) < getLength_( tourCopy, distances ) );
  }

  bool update3Opt_( size_t i,
                    size_t j,
                    size_t k,
                    vector< size_t >& tour,
                    const vector< vector< double > >& distances )
  {
    // Because all combinations are tried, the order of i, j, and k doesn't matter.
    vector< size_t > vertices( 3 );
    vertices[ 0 ] = i; vertices[ 1 ] = j; vertices[ 2 ] = k;
    sort( vertices.begin(), vertices.end() );
    i = vertices[ 0 ]; j = vertices[ 1 ]; k = vertices[ 2 ];
    assert( i < j && j < k );
    const size_t tourI = tour[ i ];
    const size_t iMinus1 = ( i + tour.size() - 1 ) % tour.size();
    const size_t tourIminus1 = tour[ iMinus1 ];
    const size_t tourJ = tour[ j ];
    const size_t jMinus1 = ( j + tour.size() - 1 ) % tour.size();
    const size_t tourJminus1 = tour[ jMinus1 ];
    const size_t tourK = tour[ k ];
    const size_t kMinus1 = ( k + tour.size() - 1 ) % tour.size();
    const size_t tourKminus1 = tour[ kMinus1 ];
    const double eps = 1e-9;
    const double removedDistance = distances[ tourIminus1 ][ tourI ] +
                                   distances[ tourJminus1 ][ tourJ ] +
                                   distances[ tourKminus1 ][ tourK ] - eps; // subtract a little something to avoid numerical errors
    vector< double > newDistances( 4 );
    newDistances[ 0 ] = distances[ tourI ][ tourJ ] + distances[ tourK ][ tourJminus1 ] + distances[ tourIminus1 ][ tourKminus1 ];
    newDistances[ 1 ] = distances[ tourJ ][ tourK ] + distances[ tourI ][ tourKminus1 ] + distances[ tourJminus1 ][ tourIminus1 ];
    newDistances[ 2 ] = distances[ tourK ][ tourI ] + distances[ tourJ ][ tourIminus1 ] + distances[ tourKminus1 ][ tourJminus1 ];
    newDistances[ 3 ] = distances[ tourI ][ tourKminus1 ] + distances[ tourJ ][ tourIminus1 ] + distances[ tourK ][ tourJminus1 ];
    size_t minIndex = min_element( newDistances.begin(), newDistances.end() ) - newDistances.begin();
    if ( newDistances[ minIndex ] < removedDistance ) {
      switch ( minIndex ) {
        case 0: reverse( tour.begin() + i, tour.begin() + j ); reverse( tour.begin() + i, tour.begin() + k ); break;
        case 1: reverse( tour.begin() + i, tour.begin() + j ); reverse( tour.begin() + j, tour.begin() + k ); break;
//        case 2: reverse( tour.begin() + i, tour.begin() + k ); crossingReverse( tour, k, i ); break;
//        case 0: update3OptIntervals_( tour, i, j, false, k, i, false, kMinus1, jMinus1, true, distances ); break;
//        case 1: update3OptIntervals_( tour, i, j, false, iMinus1, kMinus1, true, j, k, false, distances ); break;
        case 2: update3OptIntervals_( tour, i, j, false, kMinus1, jMinus1, true, iMinus1, kMinus1, true, distances ); break;
        case 3: update3OptIntervals_( tour, i, j, false, k, i, false, j, k, false, distances ); break;
        default: assert( false );
      }
      return true;
    }
    return false;
  }

  /*
  double getGain_( size_t i,
                   size_t j,
                   size_t k,
                   vector< size_t >& tour,
                   const vector< vector< double > >& distances )
  {
    vector< size_t > vertices( 3 );
    vertices[ 0 ] = i; vertices[ 1 ] = j; vertices[ 2 ] = k;
    sort( vertices.begin(), vertices.end() );
    i = vertices[ 0 ]; j = vertices[ 1 ]; k = vertices[ 2 ];
    assert( i < j && j < k );
    const size_t tourI = tour[ i ];
    const size_t iMinus1 = ( i + tour.size() - 1 ) % tour.size();
    const size_t tourIminus1 = tour[ iMinus1 ];
    const size_t tourJ = tour[ j ];
    const size_t jMinus1 = ( j + tour.size() - 1 ) % tour.size();
    const size_t tourJminus1 = tour[ jMinus1 ];
    const size_t tourK = tour[ k ];
    const size_t kMinus1 = ( k + tour.size() - 1 ) % tour.size();
    const size_t tourKminus1 = tour[ kMinus1 ];
    const double eps = 1e-9;
    const double removedDistance = distances[ tourIminus1 ][ tourI ] +
                                   distances[ tourJminus1 ][ tourJ ] +
                                   distances[ tourKminus1 ][ tourK ] - eps; // subtract a little something to avoid numerical errors
    vector< double > newDistances( 4 );
    newDistances[ 0 ] = distances[ tourI ][ tourJ ] + distances[ tourK ][ tourJminus1 ] + distances[ tourIminus1 ][ tourKminus1 ];
    newDistances[ 1 ] = distances[ tourJ ][ tourK ] + distances[ tourI ][ tourKminus1 ] + distances[ tourJminus1 ][ tourIminus1 ];
    newDistances[ 2 ] = distances[ tourK ][ tourI ] + distances[ tourJ ][ tourIminus1 ] + distances[ tourKminus1 ][ tourJminus1 ];
    newDistances[ 3 ] = distances[ tourI ][ tourKminus1 ] + distances[ tourJ ][ tourIminus1 ] + distances[ tourK ][ tourJminus1 ];
    size_t minIndex = min_element( newDistances.begin(), newDistances.end() ) - newDistances.begin();
    return removedDistance - newDistances[ minIndex ];
  }

  void compute3OptTour2_( vector< size_t >& tour,
                         const vector< vector< double > >& distances,
                         const vector< vector< size_t > >& nearestNeighbors )
  {
    vector< size_t > position( tour.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < tour.size(); ++i ) {
        for ( vector< size_t >::const_iterator jIt = nearestNeighbors[ i ].begin(); jIt != nearestNeighbors[ i ].end(); ++jIt ) {
          size_t bestI = 0;
          size_t bestJ = 0;
          size_t bestK = 0;
          double bestGain = 0.0;
          size_t numGains = 0;

          size_t indexOfIIntour = position[ i ];
          size_t indexOfJIntour = position[ *jIt ];
          assert( indexOfJIntour != indexOfIIntour );
          for ( vector< size_t >::const_iterator kIt = nearestNeighbors[ *jIt ].begin(); kIt != nearestNeighbors[ *jIt ].end(); ++kIt ) {
            size_t indexOfKIntour = position[ *kIt ];
            assert( indexOfKIntour != indexOfJIntour );
            if ( indexOfKIntour == indexOfIIntour ) {
              continue;
            }
            double gain = getGain_( indexOfIIntour, indexOfJIntour, indexOfKIntour, tour, distances ) ;
            if ( gain > 0 ) {
              ++numGains;
            }
            if ( gain > bestGain ) {
              bestGain = gain;
              bestI = indexOfIIntour;
              bestJ = indexOfJIntour;
              bestK = indexOfKIntour;
              changed = true;
            }
            if ( numGains >= 8 ) {
              break;
            }
          }
          if ( bestGain > 0.0 ) {
            assert( update3Opt_( bestI, bestJ, bestK, tour, distances ) );
            for ( size_t pi = 0; pi < tour.size(); ++pi ) {
              position[ tour[ pi ] ] = pi;
            }
            break;
          }
        }
      }
    }
  }
  */

  void compute3OptTourOld_( vector< size_t >& tour,
                            const vector< vector< double > >& distances,
                            const vector< vector< size_t > >& nearestNeighbors )
  {
    vector< size_t > position( tour.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < tour.size(); ++i ) {
        for ( vector< size_t >::const_iterator jIt = nearestNeighbors[ i ].begin(); jIt != nearestNeighbors[ i ].end(); ++jIt ) {
          size_t indexOfIIntour = position[ i ];
          size_t indexOfJIntour = position[ *jIt ];
          assert( indexOfJIntour != indexOfIIntour );
          for ( vector< size_t >::const_iterator kIt = nearestNeighbors[ *jIt ].begin(); kIt != nearestNeighbors[ *jIt ].end(); ++kIt ) {
            size_t indexOfKIntour = position[ *kIt ];
            assert( indexOfKIntour != indexOfJIntour );
            if ( indexOfKIntour == indexOfIIntour ) {
              continue;
            }
            if ( update3Opt_( indexOfIIntour, indexOfJIntour, indexOfKIntour, tour, distances ) ) {
              changed = true;
              for ( size_t i = 0; i < tour.size(); ++i ) {
                position[ tour[ i ] ] = i;
              }
              break;
            }
          }
        }
      }
    }
  }

  bool update5Opt_( size_t i,
                    size_t j,
                    size_t k,
                    size_t l,
                    size_t m,
                    vector< size_t >& tour,
                    const vector< vector< double > >& distances )
  {
    vector< size_t > vertices( 5 );
    vertices[ 0 ] = i; vertices[ 1 ] = j; vertices[ 2 ] = k; vertices[ 3 ] = l; vertices[ 4 ] = m;
    sort( vertices.begin(), vertices.end() );
    i = vertices[ 0 ]; j = vertices[ 1 ]; k = vertices[ 2 ]; l = vertices[ 3 ]; m = vertices[ 4 ];
    assert( i < j && j < k && k < l && l < m);
    const size_t tourI = tour[ i ];
    const size_t iMinus1 = ( i + tour.size() - 1 ) % tour.size();
    const size_t tourIminus1 = tour[ iMinus1 ];
    const size_t tourJ = tour[ j ];
    const size_t jMinus1 = ( j + tour.size() - 1 ) % tour.size();
    const size_t tourJminus1 = tour[ jMinus1 ];
    const size_t tourK = tour[ k ];
    const size_t kMinus1 = ( k + tour.size() - 1 ) % tour.size();
    const size_t tourKminus1 = tour[ kMinus1 ];
    const size_t tourL = tour[ l ];
    const size_t lMinus1 = ( l + tour.size() - 1 ) % tour.size();
    const size_t tourLminus1 = tour[ lMinus1 ];
    const size_t tourM = tour[ m ];
    const size_t mMinus1 = ( m + tour.size() - 1 ) % tour.size();
    const size_t tourMminus1 = tour[ mMinus1 ];
    const double eps = 1e-9;
    const double removedDistance = distances[ tourIminus1 ][ tourI ] +
                                   distances[ tourJminus1 ][ tourJ ] +
                                   distances[ tourKminus1 ][ tourK ] +
                                   distances[ tourLminus1 ][ tourL ] +
                                   distances[ tourMminus1 ][ tourM ] - eps; // subtract a little something to avoid numerical errors
    double newDistance = distances[ tourI ][ tourKminus1 ] +
                         distances[ tourJ ][ tourLminus1 ] +
                         distances[ tourK ][ tourMminus1 ] +
                         distances[ tourL ][ tourIminus1 ] +
                         distances[ tourM ][ tourJminus1 ];
    if ( newDistance < removedDistance ) {
      vector< size_t > tourCopy( tour );
      size_t tourIndex = i;
      copy( tourCopy.begin() + l, tourCopy.begin() + m, tour.begin() + tourIndex );
      tourIndex += m - l;
      copy( tourCopy.begin() + k, tourCopy.begin() + l, tour.begin() + tourIndex );
      tourIndex += l - k;
      copy( tourCopy.begin() + j, tourCopy.begin() + k, tour.begin() + tourIndex );
      tourIndex += k - j;
      copy( tourCopy.begin() + i, tourCopy.begin() + j, tour.begin() + tourIndex );
      return true;
    }
    return false;
  }

  void compute5OptTour_( vector< size_t >& tour,
                         const vector< vector< double > >& distances,
                         const vector< vector< size_t > >& nearestNeighbors )
  {
    vector< size_t > position( tour.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }
    bool changed = true;
    size_t maxNeighbors = 10;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < tour.size(); ++i ) {
        for ( vector< size_t >::const_iterator jIt = nearestNeighbors[ i ].begin(); jIt != nearestNeighbors[ i ].end(); ++jIt ) {
          if ( jIt - nearestNeighbors[ i ].begin() > int( maxNeighbors ) ) {
            break;
          }
          for ( vector< size_t >::const_iterator kIt = nearestNeighbors[ *jIt ].begin(); kIt != nearestNeighbors[ *jIt ].end(); ++kIt ) {
            if ( kIt - nearestNeighbors[ *jIt ].begin() > int( maxNeighbors ) ) {
              break;
            }
            for ( vector< size_t >::const_iterator lIt = nearestNeighbors[ *kIt ].begin(); lIt != nearestNeighbors[ *kIt ].end(); ++lIt ) {
              if ( lIt - nearestNeighbors[ *kIt ].begin() > int( maxNeighbors ) ) {
                break;
              }
              for ( vector< size_t >::const_iterator mIt = nearestNeighbors[ *lIt ].begin(); mIt != nearestNeighbors[ *lIt ].end(); ++mIt ) {
                if ( mIt - nearestNeighbors[ *lIt ].begin() > int( maxNeighbors ) ) {
                  break;
                }
                size_t indexOfIIntour = position[ i ];
                size_t indexOfJIntour = position[ *jIt ];
                size_t indexOfKIntour = position[ *kIt ];
                size_t indexOfLIntour = position[ *lIt ];
                size_t indexOfMIntour = position[ *mIt ];
                assert( indexOfJIntour != indexOfIIntour );
                assert( indexOfKIntour != indexOfJIntour );
                assert( indexOfLIntour != indexOfKIntour );
                assert( indexOfMIntour != indexOfLIntour );
                if ( indexOfKIntour == indexOfIIntour || indexOfLIntour == indexOfIIntour || indexOfLIntour == indexOfJIntour ||
                     indexOfMIntour == indexOfIIntour || indexOfMIntour == indexOfJIntour || indexOfMIntour == indexOfKIntour ) {
                  continue;
                }
                if ( update5Opt_( indexOfIIntour, indexOfJIntour, indexOfKIntour, indexOfLIntour, indexOfMIntour, tour, distances ) ) {
                  changed = true;
                  for ( size_t i = 0; i < tour.size(); ++i ) {
                    position[ tour[ i ] ] = i;
                  }
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  bool contains_( const vector< pair< size_t, size_t > >& xs, const pair< size_t, size_t >& x )
  {
    for ( vector< pair< size_t, size_t > >::const_iterator it = xs.begin(); it != xs.end(); ++it ) {
      if ( it->first == x.first && it->second == x.second ) {
        return true;
      }
      if ( it->second == x.first && it->first == x.second ) {
        return true;
      }
    }
    return false;
  }

  void make2OptMove_( size_t it0, size_t it1, size_t it2, size_t it3, vector< size_t >& tour )
  {
    cerr << it0 << " " << it1 << " " << it2 << " " << it3 << endl;
    assert( fabs( fabs( float( it0 ) - float( it1 ) ) - 1.0 ) < 1e-6 || fabs( fabs( float( it0 ) - float( it1 ) ) + 1.0 - float( tour.size() ) ) < 1e-6  );
    assert( fabs( fabs( float( it2 ) - float( it3 ) ) - 1.0 ) < 1e-6 || fabs( fabs( float( it2 ) - float( it3 ) ) + 1.0 - float( tour.size() ) ) < 1e-6  );
    size_t ind1 = max( it0, it1 );
    size_t ind2 = max( it2, it3 );
    if ( ind1 > ind2 ) {
      swap( ind1, ind2 );
    }
    reverse( tour.begin() + ind1, tour.begin() + ind2 );
  }

  void computeLinKernighantour_( vector< size_t >& tour,
                                 const vector< vector< double > >& distances )
  {
    bool changed = true;
    vector< pair< size_t, size_t > > x;
    vector< pair< size_t, size_t > > y;
    vector< size_t > t;
    const double eps = 1e-9;
    while ( changed ) {
      changed = false;
restartLinKernighan:
      // Step 2. Let ind = 0. Choose t_0
      for ( size_t i = 1; i < tour.size(); ++i ) {
        vector< pair< size_t, size_t > > T;
        for ( size_t j = 0; j < tour.size(); ++j ) {
          T.push_back( make_pair( tour[ j ], tour[ ( j + 1 ) % tour.size() ] ) );
        }
        size_t ind = 0;
        x.clear();
        y.clear();
        t.clear();
        t.push_back( tour[ i == 0 ? tour.size() - 1 : i - 1 ] );
        t.push_back( tour[ i ] );
        // Step 3. Choose x_0 = ( t_0, t_1 ) in T
        x.push_back( make_pair( t[ 0 ], t[ 1 ] ) );
        for ( size_t j = i + 2; j < tour.size(); ++j ) {
          // Step 4. Choose y_0 = ( t_1, t_2 ) not in T such that G > 0
          if ( i < 2 && j + 1 == tour.size() ) {
            continue;
          }
          if ( contains_( T, make_pair( t[ 1 ], tour[ j ] ) ) ) {
            // y is in T
            continue;
          }
          double G0 = distances[ t[ 0 ] ][ t[ 1 ] ] - distances[ t[ 1 ] ][ tour[ j ] ];
          if ( G0 <= eps ) {
            continue;
          }
          double G = G0;

          // Found y not in T with positive gain
          t.push_back( tour[ j ] );
          y.push_back( make_pair( t[ 1 ], t[ 2 ] ) ); // y_0

          size_t maxKForIndEquals2 = 0;
          bool nextYExists = true;
          bool nextXExists = true;
          while ( nextYExists ) {
            // Step 5. Let ind = ind + 1
            ++ind;
            // Step 6. Choose x_i = (t_(2i), t_(2i+1)) in T such that
            // (a) if t_(2i+1) is joined to t_0, the resulting configuration is a tour T'
            // (b) x_i != y_s for all s < i
            size_t tNext = tour[ ( ( find( tour.begin(), tour.end(), t.back() ) - tour.begin() ) + tour.size() - 1 ) % tour.size() ]; // element in tour previous to t_(2i)
            assert( nextXExists ); // condition (b): x_i != y_s for all s < i
            x.push_back( make_pair( t.back(), tNext ) ); // Add x_i = (t_(2i), t_(2i+1))
            t.push_back( tNext ); // Add t_(2i+1)

            cerr << ind << ". x size: " << x.size() << ", G: " << G << ", return: " << G + distances[ x.back().first ][ x.back().second ] - distances[ t.back() ][ t.front() ] << endl;
            if ( G + distances[ x.back().first ][ x.back().second ] - distances[ t.back() ][ t.front() ] > eps ) {
              y.push_back( make_pair( t.back(), t.front() ) );
              changed = true;
              // Take tour
              vector< size_t > tourCopy( tour );
              assert( t.size() % 2 == 0 );
              assert( x.size() == y.size() );
              cerr << "Take tour" << endl;
              cerr << "x.size() " << x.size() << endl;
              cerr << "G: " << G << " dist: " << distances[ t.back() ][ t.front() ] << endl;
              cerr << "tour: ";
              for ( size_t k = 0; k < tour.size(); ++k ) {
                cerr << tour[ k ] << " ";
              }
              cerr << endl;
              cerr << "x: ";
              for ( size_t k = 0; k < x.size(); ++k ) {
                cerr << "(" << x[ k ].first << ", " << x[ k ].second << "), ";
              }
              cerr << endl;
              cerr << "y: ";
              for ( size_t k = 0; k < y.size(); ++k ) {
                cerr << "(" << y[ k ].first << ", " << y[ k ].second << "), ";
              }
              cerr << endl;
              for ( size_t k = 1; k < x.size(); ++k ) {
                size_t t0 = find( tour.begin(), tour.end(), x[ 0 ].first ) - tour.begin();
                size_t t1 = find( tour.begin(), tour.end(), x[ k - 1 ].second ) - tour.begin();
                size_t t2 = find( tour.begin(), tour.end(), x[ k ].first ) - tour.begin();
                size_t t3 = find( tour.begin(), tour.end(), x[ k ].second ) - tour.begin();
                cerr << "2opt: (" << x[ 0 ].first << ", " << x[ k - 1].second << "), (" << x[ k ].first << ", " << x[ k ].second << ")" << endl;
                make2OptMove_( t0, t1, t2, t3, tour );
              }
              assert( tour != tourCopy );

              if ( getLength_( tour, distances ) >= getLength_( tourCopy, distances ) ) {
                cerr << setprecision( 28 ) << getLength_( tour, distances ) << " " << getLength_( tourCopy, distances ) << endl;
              }
              assert( getLength_( tour, distances ) < getLength_( tourCopy, distances ) );
              goto restartLinKernighan;
            }

step7:
            nextYExists = false;
            // Step 7. Select y_i = (t_(2i+1), t_(2i+2)) not in T such that
            // (a) G_i > 0
            // (b) y_i != x_s for all s <= i
            // (c) x_(i+1) exists
            // If such y_i exists, go to step 5
            size_t startK = y.size() == 1 ? maxKForIndEquals2 : 0;
            for ( size_t k = startK; k < tour.size(); ++k ) {
              if ( y.size() == 1 ) {
                ++maxKForIndEquals2;
              }
              if ( k == i || k == j || t.back() == tour[ k ] ) {
                continue;
              }
              pair< size_t, size_t > yCandidate;
              yCandidate = make_pair( t.back(), tour[ k ] );

              // y is not in T
              if ( contains_( T, yCandidate ) ) {
                continue;
              }
              // (a) G_i > 0
              double gain = distances[ x.back().first ][ x.back().second ] - distances[ t.back() ][ tour[ k ] ];
              if ( G + gain <= eps ) {
                continue;
              }
              // (b) y_i != x_s for all s <= i
              if ( contains_( x, yCandidate ) ) {
                continue;
              }
              // (c) x_(i+1) exists
              size_t tNextCandidate = tour[ ( ( find( tour.begin(), tour.end(), tour[ k ] ) - tour.begin() ) + tour.size() - 1 ) % tour.size() ]; // element in tour previous to t_(2i+2) candidate
              pair< size_t, size_t > xCandidate;
              xCandidate = make_pair( tour[ k ], tNextCandidate );
              nextXExists = !contains_( x, xCandidate ) &&
                            !contains_( y, xCandidate ) &&
                            tNextCandidate != t.back();
              if ( !nextXExists ) {
                continue;
              }
              G += gain;
              y.push_back( yCandidate );
              t.push_back( tour[ k ] );

              // Found y, goto step 5
              nextYExists = true;
              break;
            }

            if ( !nextYExists && maxKForIndEquals2 < tour.size() ) {
              // Step 8. If there is an untried alternative for y_1, let i = 1 and go to Step 7
              y.resize( 1 );
              x.resize( 2 );
              t.resize( 4 );
              ind = 1;
              G = G0;
              goto step7;
            }
          }
        }
      }
    }
  }

  size_t previous_( size_t node, const vector< size_t >& tour, const vector< size_t >& position )
  {
    return tour[ ( position[ node ] + tour.size() - 1 ) % tour.size() ];
  }

  size_t next_( size_t node, const vector< size_t >& tour, const vector< size_t >& position )
  {
    return tour[ ( position[ node ] + 1 ) % tour.size() ];
  }

  bool between_( size_t a, size_t b, size_t c, const vector< size_t >& position )
  {
    return ( position[ a ] < position[ c ] && position[ c ] < position[ b ] ) ||
           ( position[ b ] < position[ a ] && position[ a ] < position[ c ] ) ||
           ( position[ c ] < position[ b ] && position[ b ] < position[ a ] );
  }

  void flip_( vector< size_t >& tour, vector< size_t >& position, size_t t1, size_t t2, size_t t3, size_t t4 )
  {

    if ( !( t2 == next_( t1, tour, position ) && t3 == next_( t4, tour, position ) ) ) {
      cerr << "tpos: " << position[ t2 ] << " " << position[ t1 ] << " " << position[ t3 ] << " " << position[ t4 ] << endl;
    }
    assert( t2 == next_( t1, tour, position ) && t3 == next_( t4, tour, position ) );
    size_t pt2 = position[ t2 ];
    size_t pt3 = position[ t3 ];
    if ( pt2 < pt3 ) {
      reverse( tour.begin() + pt2, tour.begin() + pt3 );
    }
    else {
      vector< size_t >::iterator iIt = tour.begin() + pt2;
      vector< size_t >::iterator jIt = tour.begin() + pt3;
      while ( iIt != tour.end() && jIt != tour.begin() ) {
        --jIt;
        iter_swap( iIt, jIt );
        ++iIt;
      }
      if ( iIt == tour.end() ) {
        reverse( tour.begin(), jIt );
      }
      else {
        reverse( iIt, tour.end() );
      }
    }
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }
  }

  void performMove_( const vector< size_t >& ts, vector< size_t >& tour, vector< size_t >& position )
  {
    assert( ts.size() == 4 || ts.size() == 6 );
    if ( ts.size() == 4 ) {
      size_t t1 = ts[ 0 ];
      size_t t2 = ts[ 1 ];
      size_t t3 = ts[ 2 ];
      size_t t4 = ts[ 3 ];
      if ( t1 == next_( t2, tour, position ) ) {
        swap( t1, t2 );
      }
      if ( t4 == next_( t3, tour, position ) ) {
        swap( t3, t4 );
      }
      flip_( tour, position, t1, t2, t3, t4 );
    }
    else {
      size_t t1 = ts[ 0 ];
      size_t t2 = ts[ 1 ];
      size_t t3 = ts[ 2 ];
      size_t t4 = ts[ 3 ];
      size_t t5 = ts[ 4 ];
      size_t t6 = ts[ 5 ];

//      cerr << "ts: " << t1 << " " << t2 << " " << t3 << " " << t4 << " " << t5 << " " << t6 << endl;
//      cerr << "pts: " << position[t1] << " " << position[t2] << " " << position[t3] << " " << position[t4] << " " << position[t5] << " " << position[t6] << endl;

//      cerr << "tour before: ";
//      for ( size_t i = 0; i< tour.size(); ++i ) {
//        cerr << tour[ i ] << " ";
//      }
//      cerr << endl;
      if ( t1 == next_( t2, tour, position ) ) {
        swap( t1, t2 );
      }
      if ( t4 == next_( t3, tour, position ) ) {
        swap( t3, t4 );
      }
      flip_( tour, position, t1, t2, t3, t4 );
//      cerr << "tour intermediate: ";
//      for ( size_t i = 0; i< tour.size(); ++i ) {
//        cerr << tour[ i ] << " ";
//      }
//      cerr << endl;

      t1 = ts[ 0 ];
      t4 = ts[ 3 ];
      t5 = ts[ 4 ];
      t6 = ts[ 5 ];

      t2 = ts[ 1 ];
      t3 = ts[ 2 ];

//      cerr << "pts2: " << position[t1] << " " << position[t2] << " " << position[t3] << " " << position[t4] << " " << position[t5] << " " << position[t6] << endl;

      if ( t1 == next_( t4, tour, position ) ) {
        swap( t1, t4 );
      }
      if ( t6 == next_( t5, tour, position ) ) {
        swap( t5, t6 );
      }
      flip_( tour, position, t1, t4, t5, t6 );
    }
  }

  void compute3OptTour_( vector< size_t >& tour,
                         const vector< vector< double > >& distances,
                         const vector< vector< size_t > >& nearestNeighbors )
  {
    const double eps = 1e-9;
    vector< size_t > position( tour.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }

    bool changed = true;
    vector< size_t > bestTs;
    while ( changed ) {
      changed = false;
      for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
        bestTs.assign( 1, t1 );
        double G = 0.0;
        for ( size_t t2choice = 0; t2choice < 2; ++t2choice  ) {
          size_t t2 = t2choice == 0 ? previous_( t1, tour, position ) : next_( t1, tour, position );
          for ( size_t t3index = 0; t3index < nearestNeighbors[ t2 ].size(); ++t3index ) {
            size_t t3 = nearestNeighbors[ t2 ][ t3index ];
            if ( t3 == previous_( t2, tour, position ) || t3 == next_( t2, tour, position ) ) {
              continue;
            }
            double g1 = distances[ t1 ][ t2 ] - distances[ t2 ][ t3 ] - eps;
            if ( g1 <= 0.0 ) {
              continue;
            }
            // First choice of t4
            size_t t4 = t2choice == 0 ? next_( t3, tour, position ) : previous_( t3, tour, position );
            if ( t4 == previous_( t2, tour, position ) || t4 == next_( t2, tour, position ) ) {
              continue;
            }
            {
              // Test for improving 2-opt move
              double gain = g1 + distances[ t3 ][ t4 ] - distances[ t4 ][ t1 ];
              if ( gain > G ) {
                G = gain;
                bestTs.resize( 1 );
                bestTs.push_back( t2 );
                bestTs.push_back( t3 );
                bestTs.push_back( t4 );
              }
            }
            for ( size_t t5index = 0; t5index < nearestNeighbors[ t4 ].size(); ++t5index ) {
              size_t t5 = nearestNeighbors[ t4 ][ t5index ];
              if ( t5 == previous_( t4, tour, position ) || t5 == next_( t4, tour, position ) ) {
                continue;
              }
              double g2 = distances[ t3 ][ t4 ] - distances[ t4 ][ t5 ];
              if ( g1 + g2 <= 0.0 ) {
                continue;
              }

              // Select t6 such that a valid tour is created
//              bool betweenT3T1 = between_( t3, t1, t5, tour, position );
//                cerr << "type 1 " << t5 << " " << next_( t5, tour, position ) << " " << previous_( t5, tour, position ) << endl;
              size_t t6 = between_( t2, t4, t5, position ) ? next_( t5, tour, position ) : previous_( t5, tour, position );
              if ( t6 == t1 ) {
                continue;
              }
              double g3 = distances[ t5 ][ t6 ] - distances[ t6 ][ t1 ];
              double gain = g1 + g2 + g3;
              if ( gain > G ) {
                G = gain;
//                cerr << "set gain1 " << position[ t1 ] << " " << position[ t2] << endl;
                bestTs.resize( 1 );
                bestTs.push_back( t2 );
                bestTs.push_back( t3 );
                bestTs.push_back( t4 );
                bestTs.push_back( t5 );
                bestTs.push_back( t6 );
              }
            }

            // Second choice of t4
            continue;

            t4 = t2choice == 0 ? previous_( t3, tour, position ) : next_( t3, tour, position );
            for ( size_t t5index = 0; t5index < nearestNeighbors[ t4 ].size(); ++t5index ) {
              size_t t5 = nearestNeighbors[ t4 ][ t5index ];
              if ( t5 == previous_( t4, tour, position ) || t5 == next_( t4, tour, position ) ) {
                continue;
              }
              if ( ( t2choice == 0 && !between_( t3, t2, t5, position ) ) ||
                   ( t2choice == 1 && !between_( t2, t3, t5, position ) ) ) {
                continue;
              }
              double g2 = distances[ t3 ][ t4 ] - distances[ t4 ][ t5 ];
              if ( g1 + g2 <= 0.0 ) {
                continue;
              }
              for ( size_t t6choice = 0; t6choice < 2; ++t6choice ) {
                size_t t6 = t6choice == 0 ? next_( t5, tour, position ) : previous_( t5, tour, position );
                if ( t6 == t3 || t6 == t2 ) {
                  continue;
                }
                if ( t6 == previous_( t1, tour, position ) || t6 == next_( t1, tour, position ) ||
                     t6 == previous_( t2, tour, position ) || t6 == next_( t2, tour, position ) ||
                     t6 == previous_( t4, tour, position ) || t6 == next_( t4, tour, position ) ) {
                  continue;
                }
/*                cerr << "type 2" << endl;
               cerr << "t2choice " << t2choice << endl;
              cerr << "t1-5: " << t1 << " " << t2 << " " << t3 << " " << t4 << " " << t5 << " " << t6 << endl;
              cerr << "p t1-5: " << position[ t1 ] << " " << position[ t2 ] << " " << position[ t3 ] << " " << position[ t4 ] << " " << position[ t5 ] << " " << position[ t6 ] << endl;
              */
                double g3 = distances[ t5 ][ t6 ] - distances[ t6 ][ t1 ];
                double gain = g1 + g2 + g3;
                if ( gain > G ) {
                  G = gain;
//                  cerr << "set gain2 " << position[t1] << " " << position[t2] << endl;
                  bestTs.resize( 1 );
                  bestTs.push_back( t2 );
                  bestTs.push_back( t5 );
                  bestTs.push_back( t6 );
                  bestTs.push_back( t3 );
                  bestTs.push_back( t4 );
                }
              }
            }
          }
        }
        if ( G > 0.0 ) {
          changed = true;
          vector< size_t > tourCopy( tour );
          performMove_( bestTs, tour, position );
          if ( getLength_( tour, distances ) >= getLength_( tourCopy, distances ) ) {
            cerr << getLength_( tour, distances ) << " " << getLength_( tourCopy, distances ) << endl;
          }
          assert( getLength_( tour, distances ) < getLength_( tourCopy, distances ) );
        }
      }
    }
  }

  void compute2OptTour_( vector< size_t >& tour,
                         const vector< vector< double > >& distances,
                         const vector< vector< size_t > >& nearestNeighbors )
  {
    const double eps = 1e-9;
    vector< size_t > position( tour.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }

    vector< size_t > bestTs;
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
        bestTs.assign( 4, (size_t)-1 );
        double maxGain = 0.0;
        for ( size_t t2choice = 0; t2choice < 2; ++t2choice  ) {
          size_t t2 = t2choice == 0 ? previous_( t1, tour, position ) : next_( t1, tour, position );
          for ( size_t t3index = 0; t3index < nearestNeighbors[ t2 ].size(); ++t3index ) {
            size_t t3 = nearestNeighbors[ t2 ][ t3index ];
            if ( t3 == previous_( t2, tour, position ) || t3 == next_( t2, tour, position ) ) {
              continue;
            }
            if ( distances[ t2 ][ t3 ] >= distances[ t1 ][ t2 ] ) {
              break;
            }
            size_t t4 = t2choice == 0 ? next_( t3, tour, position ) : previous_( t3, tour, position );
            double gain = distances[ t1 ][ t2 ] + distances[ t3 ][ t4 ] - ( distances[ t1 ][ t4 ] + distances[ t2 ][ t3 ] ) - eps;
            if ( gain > maxGain ) {
              maxGain = gain;
              bestTs[ 0 ] = t1;
              bestTs[ 1 ] = t2;
              bestTs[ 2 ] = t3;
              bestTs[ 3 ] = t4;
            }
          }
        }
        if ( maxGain > 0.0 ) {
          changed = true;
          performMove_( bestTs, tour, position );
          for ( size_t i = 0; i < tour.size(); ++i ) {
            position[ tour[ i ] ] = i;
          }
        }
      }
    }
  }

  bool update2Opt_( size_t i,
                    size_t j,
                    vector< size_t >& tour,
                    const vector< vector< double > >& distances )
  {
    assert( i < j );
    const double eps = 1e-9;
    const size_t t0 = tour[ i == 0 ? tour.size() - 1 : i - 1 ];
    const size_t t1 = tour[ i ];
    const size_t t2 = tour[ j == 0 ? tour.size() - 1 : j - 1 ];
    const size_t t3 = tour[ j % tour.size() ];
    if ( distances[ t0 ][ t2 ] + distances[ t1 ][ t3 ] <
         distances[ t0 ][ t1 ] + distances[ t2 ][ t3 ] - eps ) {
      reverse( tour.begin() + i, tour.begin() + j );
      return true;
    }
    return false;
  }

  double getGain_( size_t i,
                   size_t j,
                   vector< size_t >& tour,
                   const vector< vector< double > >& distances )
  {
    assert( i < j );
    const double eps = 1e-9;
    const size_t t0 = tour[ i == 0 ? tour.size() - 1 : i - 1 ];
    const size_t t1 = tour[ i ];
    const size_t t2 = tour[ j == 0 ? tour.size() - 1 : j - 1 ];
    const size_t t3 = tour[ j % tour.size() ];
    return ( distances[ t0 ][ t1 ] + distances[ t2 ][ t3 ] - ( distances[ t0 ][ t2 ] + distances[ t1 ][ t3 ] ) ) - eps;
  }

  void compute2OptTour_( vector< size_t >& tour,
                         const vector< vector< double > >& distances )
  {
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < tour.size(); ++i ) {
        size_t bestI = (size_t)-1;
        size_t bestJ = (size_t)-1;
        double bestGain = 0.0;
        // Cutting the edges of two neighboring nodes does not lead to a new tour. Hence start from i + 2.
        for ( size_t j = i + 2; j < tour.size(); ++j ) {
          double gain = getGain_( i, j, tour, distances ) ;
          if ( gain > bestGain ) {
            bestI = i;
            bestJ = j;
            bestGain = gain;
            changed = true;
          }
        }
        if ( bestGain > 0.0 ) {
          assert( update2Opt_( bestI, bestJ, tour, distances ) );
        }
      }
    }
    assertIsTour_( tour, distances );
  }

  void doubleBridgeSwap_( vector< size_t >& tour, vector< size_t >& position, size_t t1, size_t t2, size_t t3, size_t t4, size_t t5, size_t t6, size_t t7, size_t t8 )
  {
    assert( t2 == next_( t1, tour, position ) && t4 == next_( t3, tour, position ) && t6 == next_( t5, tour, position ) && t8 == next_( t7, tour, position ) );
    vector< size_t > tourCopy( tour );
    size_t i = 0;
    for ( size_t t = t4; t != t8; t = next_( t, tourCopy, position ), ++i ) {
      tour[ i ] = t;
    }
    for ( size_t t = t6; t != t4; t = next_( t, tourCopy, position ), ++i ) {
      tour[ i ] = t;
    }
    for ( size_t t = t2; t != t6; t = next_( t, tourCopy, position ), ++i ) {
      tour[ i ] = t;
    }
    for ( size_t t = t8; t != t2; t = next_( t, tourCopy, position ), ++i ) {
      tour[ i ] = t;
    }
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }
  }

  void computeDoubleBridgeTour_( vector< size_t >& tour,
                                 const vector< vector< double > >& distances,
                                 const vector< vector< size_t > >& nearestNeighbors )
  {
    double eps = 1e-9;
    vector< size_t > position( tour.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }
    vector< size_t > bestTs;
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
        size_t t2 = next_( t1, tour, position );
        double maxGain = 0.0;

        for ( size_t t3index = 0; t3index < nearestNeighbors[ t2 ].size(); ++t3index ) {
          size_t t3 = nearestNeighbors[ t2 ][ t3index ];
          size_t t4 = next_( t3, tour, position );
          if ( t3 == t1 ||
               t4 == t1 ||
               t3 == next_( t2, tour, position ) ||
               t4 == previous_( t1, tour, position ) ) {
            continue;
          }
          double gainFirstBridge = distances[ t1 ][ t2 ] + distances[ t3 ][ t4 ] - distances[ t2 ][ t3 ] - distances[ t1 ][ t4 ] - eps;
          if ( gainFirstBridge <= 0.0 ) {
            continue;
          }

          for ( size_t t5 = next_( t2, tour, position ); between_( t2, t3, t5, position ); t5 = next_( t5, tour, position  ) ) {
            size_t t6 = next_( t5, tour, position );
            if ( t6 == t3 ) {
              continue;
            }

            for ( size_t t7 = next_( t4, tour, position ); between_( t4, t1, t7, position ); t7 = next_( t7, tour, position  ) ) {
              size_t t8 = next_( t7, tour, position );
              if ( t8 == t1 ) {
                continue;
              }

              double gainSecondBridge = distances[ t5 ][ t6 ] + distances[ t7 ][ t8 ] - distances[ t6 ][ t7 ] - distances[ t5 ][ t8 ];
              if ( gainFirstBridge + gainSecondBridge > maxGain ) {
                maxGain = gainFirstBridge + gainSecondBridge;
                bestTs.clear();
                bestTs.push_back( t1 );
                bestTs.push_back( t2 );
                bestTs.push_back( t3 );
                bestTs.push_back( t4 );
                bestTs.push_back( t5 );
                bestTs.push_back( t6 );
                bestTs.push_back( t7 );
                bestTs.push_back( t8 );
              }
            }
          }
        }
        if ( maxGain > 0.0 ) {
          doubleBridgeSwap_( tour, position, bestTs[ 0 ], bestTs[ 1 ], bestTs[ 2 ], bestTs[ 3 ], bestTs[ 4 ], bestTs[ 5 ], bestTs[ 6 ], bestTs[ 7 ] );
          changed = true;
        }
      }
    }
  }

} // anonymous namespace


namespace TravelingSalespersonProblemSolver {
  vector< size_t > computeTour( const vector< vector< double > >& distances )
  {
    assert( distances.size() > 0 );
    srand( 1729 );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      assert( distances.size() == distances[ i ].size() );
      for ( size_t j = 0; j < distances.size(); ++j ) {
        if ( i != j ) {
          assert( distances[ i ][ j ] > 0.0 );
        }
        else {
          assert( distances[ i ][ i ] == 0.0 );
        }
      }
    }
    const vector< size_t > tourRand = getRandomtour_( distances );
    const vector< size_t > tourNN = getNearestNeighbortour_( distances );
    const vector< size_t > tourGreedy = getGreedyTour_( distances );
    vector< size_t > tour( tourGreedy );
    cerr << setprecision( 8 );

    vector< vector< size_t > > nearestNeighbors;
    {
      double start( clock() );
      nearestNeighbors = computeNearestNeighbors_( distances, 20 );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "Time to compute " << nearestNeighbors.front().size() << " nearest neighbors: " << time << endl;
    }

    cerr << "Initial distance: " << getLength_( tour, distances ) << endl;
    cerr << "Nearest neighbor distance: " << getLength_( tourNN, distances ) << endl;
    cerr << "Greedy distance: " << getLength_( tourGreedy, distances ) << endl;

    if ( true ) {
      cerr << "1-tree distance: ";
      double start( clock() );
      double lowerBound = getHeldKarpLowerBound_( distances );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << lowerBound << ", time: " << time << endl;
    }

    if ( true ) {
      double start( clock() );
      compute2OptTour_( tour, distances, nearestNeighbors );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "2-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
      assertIsTour_( tour, distances );
    }

    if ( true ) {
      tour = tourGreedy;
      double start( clock() );
      compute3OptTour_( tour, distances, nearestNeighbors );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "3-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
      assertIsTour_( tour, distances );
    }

    if ( true ) {
      double start( clock() );
      computeDoubleBridgeTour_( tour, distances, nearestNeighbors );
      compute3OptTour_( tour, distances, nearestNeighbors );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "4-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
      assertIsTour_( tour, distances );
    }

    if ( false ) {
      computeLinKernighantour_( tour, distances );
      compute5OptTour_( tour, distances, nearestNeighbors );
    }

    return tour;
  }
}
