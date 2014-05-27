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

// For profiling
#if 0
#define INLINE_ATTRIBUTE __attribute__ ((noinline))
#else
#define INLINE_ATTRIBUTE
#endif

using namespace std;

namespace TravelingSalespersonProblemSolver {
  namespace {
  void INLINE_ATTRIBUTE assertIsTour_( const vector< size_t >& tour,
                                       const vector< vector< double > >& distances )
  {
    assert( tour.size() == distances.size() );
    set< size_t > tourSet( tour.begin(), tour.end() );
    assert( tourSet.size() == tour.size() );
    assert( *tourSet.begin() == 0 );
    assert( *tourSet.rbegin() + 1 == tour.size() );
  }

  vector< size_t > INLINE_ATTRIBUTE getRandomTour_( const vector< vector< double > >& distances )
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

  vector< size_t > INLINE_ATTRIBUTE getNearestNeighborTour_( const vector< vector< double > >& distances )
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

  vector< size_t > INLINE_ATTRIBUTE getGreedyTour_( const vector< vector< double > >& distances )
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

  vector< vector< size_t > > INLINE_ATTRIBUTE computeNearestNeighbors_( const vector< vector< double > >& distances,
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

  double INLINE_ATTRIBUTE compute1TreeLength_( vector< size_t >& vertexDegrees,
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

  double INLINE_ATTRIBUTE getHeldKarpLowerBound_( const vector< vector< double > >& distances )
  {
    double maxLength = numeric_limits< double >::min();
    vector< double > lagrangeMultipliers( distances.size() );
    double lambda = 0.1;
    for ( size_t i = 0; i < 50; ++i ) {
      vector< size_t > vertexDegrees;
      double treeLength = compute1TreeLength_( vertexDegrees, distances, lagrangeMultipliers );
      maxLength = max( maxLength, treeLength );
      double denominator = 0.0;
      for ( size_t j = 0; j < lagrangeMultipliers.size(); ++j ) {
        double d = double( vertexDegrees[ j ] ) - 2.0;
        denominator += d * d;
      }
      if ( denominator == 0.0 ) {
        return maxLength;
      }
      double t = 2.0 * lambda * treeLength / denominator;

      for ( size_t j = 0; j < lagrangeMultipliers.size(); ++j ) {
        lagrangeMultipliers[ j ] += t * ( int( vertexDegrees[ j ] ) - 2 );
      }
      lambda = 1.0 / ( 20.0 + 10 * i );
    }

    return maxLength;
  }

  double INLINE_ATTRIBUTE getLength_( const vector< size_t >& tour,
                                      const vector< vector< double > >& distances )
  {
    double distance = distances[ tour.back() ][ tour.front() ];
    for ( size_t i = 0; i + 1 < tour.size(); ++i ) {
      distance += distances[ tour[ i ] ][ tour[ i + 1 ] ];
    }
    return distance;
  }

  bool INLINE_ATTRIBUTE update5Opt_( size_t i,
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

  void INLINE_ATTRIBUTE compute5OptTour_( vector< size_t >& tour,
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

  size_t INLINE_ATTRIBUTE previous_( size_t node, const vector< size_t >& tour, const vector< size_t >& position )
  {
    return position[ node ] > 0 ? tour[ position[ node ] - 1 ] : tour.back();
  }

  size_t INLINE_ATTRIBUTE next_( size_t node, const vector< size_t >& tour, const vector< size_t >& position )
  {
    return position[ node ] + 1 < tour.size() ? tour[ position[ node ] + 1 ] : tour.front();
  }

  bool INLINE_ATTRIBUTE between_( size_t a, size_t b, size_t c, const vector< size_t >& position )
  {
    return position[ a ] <= position[ c ] ? position[ a ] >= position[ b ] || position[ b ] >= position[ c ]
                                          : position[ b ] <= position[ a ] && position[ b ] >= position[ c ];
  }

  void INLINE_ATTRIBUTE flip_( vector< size_t >& tour, vector< size_t >& position, size_t t1, size_t t2, size_t t3, size_t t4 )
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

  void INLINE_ATTRIBUTE performMove_( const vector< size_t >& ts, vector< size_t >& tour, vector< size_t >& position )
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

      if ( t1 == next_( t2, tour, position ) ) {
        swap( t1, t2 );
      }
      if ( t4 == next_( t3, tour, position ) ) {
        swap( t3, t4 );
      }
      flip_( tour, position, t1, t2, t3, t4 );

      t1 = ts[ 0 ];
      t4 = ts[ 3 ];
      t5 = ts[ 4 ];
      t6 = ts[ 5 ];

      t2 = ts[ 1 ];
      t3 = ts[ 2 ];

      if ( t1 == next_( t4, tour, position ) ) {
        swap( t1, t4 );
      }
      if ( t6 == next_( t5, tour, position ) ) {
        swap( t5, t6 );
      }
      flip_( tour, position, t1, t4, t5, t6 );
    }
  }

  void INLINE_ATTRIBUTE compute3OptTour_( vector< size_t >& tour,
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
        double G = 0.0;
        for ( size_t t2choice = 0; t2choice < 2; ++t2choice  ) {
          size_t t2 = t2choice == 0 ? previous_( t1, tour, position ) : next_( t1, tour, position );
          for ( size_t t3index = 0; t3index < nearestNeighbors[ t2 ].size(); ++t3index ) {
            size_t t3 = nearestNeighbors[ t2 ][ t3index ];
            if ( t3 == previous_( t2, tour, position ) || t3 == next_( t2, tour, position ) ) {
              continue;
            }
            double g1 = distances[ t1 ][ t2 ] - distances[ t2 ][ t3 ];
            if ( g1 <= eps ) {
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
                bestTs = { t1, t2, t3, t4 };
              }
            }
            for ( size_t t5index = 0; t5index < nearestNeighbors[ t4 ].size(); ++t5index ) {
              size_t t5 = nearestNeighbors[ t4 ][ t5index ];
              if ( t5 == previous_( t4, tour, position ) || t5 == next_( t4, tour, position ) ) {
                continue;
              }
              double g2 = distances[ t3 ][ t4 ] - distances[ t4 ][ t5 ];
              if ( g1 + g2 <= eps ) {
                continue;
              }

              // Select t6 such that a valid tour is created
              size_t t6 = between_( t2, t4, t5, position ) ? next_( t5, tour, position ) : previous_( t5, tour, position );
              if ( t6 == t1 ) {
                continue;
              }
              double g3 = distances[ t5 ][ t6 ] - distances[ t6 ][ t1 ];
              double gain = g1 + g2 + g3;
              if ( gain > G ) {
                G = gain;
                bestTs = { t1, t2, t3, t4, t5, t6 };
              }
            }

            // Second choice of t4
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
              if ( g1 + g2 <= eps ) {
                continue;
              }
              // Only consider one choice of t6. The other choice is possible, but clutters the code and doesn't lead to a significant improvement.
              size_t t6choice = t2choice;
              size_t t6 = t6choice == 0 ? next_( t5, tour, position ) : previous_( t5, tour, position );
              if ( t6 == t3 || t6 == t2 || t6 == t1 ) {
                continue;
              }
              double g3 = distances[ t5 ][ t6 ] - distances[ t6 ][ t1 ];
              double gain = g1 + g2 + g3;
              if ( gain > G ) {
                G = gain;
                bestTs = { t3, t4, t5, t6, t1, t2 };
              }
            }
          }
        }
        if ( G > eps ) {
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

  void INLINE_ATTRIBUTE compute2OptTour_( vector< size_t >& tour,
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
            double gain = distances[ t1 ][ t2 ] + distances[ t3 ][ t4 ] - ( distances[ t1 ][ t4 ] + distances[ t2 ][ t3 ] );
            if ( gain > maxGain ) {
              maxGain = gain;
              bestTs[ 0 ] = t1;
              bestTs[ 1 ] = t2;
              bestTs[ 2 ] = t3;
              bestTs[ 3 ] = t4;
            }
          }
        }
        if ( maxGain > eps ) {
          changed = true;
          performMove_( bestTs, tour, position );
        }
      }
    }
  }

  void INLINE_ATTRIBUTE doubleBridgeSwap_( vector< size_t >& tour, vector< size_t >& position, size_t t1, size_t t2, size_t t3, size_t t4, size_t t5, size_t t6, size_t t7, size_t t8 )
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

  void INLINE_ATTRIBUTE computeDoubleBridgeTour_( vector< size_t >& tour,
                                                  const vector< vector< double > >& distances,
                                                  const vector< vector< size_t > >& nearestNeighbors )
  {
    if ( tour.size() < 5 ) {
      return;
    }
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
        size_t t2 = next_( t1, tour, position );
        double maxGain = 0.0;

        for ( size_t t3index = 0; t3index < nearestNeighbors[ t2 ].size(); ++t3index ) {
          size_t t3 = nearestNeighbors[ t2 ][ t3index ];
          size_t t4 = next_( t3, tour, position );
          if ( t3 == t1 || t4 == t1 ) {
            continue;
          }
          double gainFirstBridge = distances[ t1 ][ t2 ] + distances[ t3 ][ t4 ] - distances[ t2 ][ t3 ] - distances[ t1 ][ t4 ];
          if ( gainFirstBridge <= 0.0 ) {
            continue;
          }

          for ( size_t t5 = t2; t5 != t3; t5 = next_( t5, tour, position  ) ) {
            size_t t6 = next_( t5, tour, position );
            for ( size_t t7index = 0; t7index < nearestNeighbors[ t6 ].size(); ++t7index ) {
              size_t t7 = nearestNeighbors[ t6 ][ t7index ];
              if ( !between_( t4, t1, t7, position ) || t7 == t1 ) {
                continue;
              }
              size_t t8 = next_( t7, tour, position );

              double gainSecondBridge = distances[ t5 ][ t6 ] + distances[ t7 ][ t8 ] - distances[ t6 ][ t7 ] - distances[ t5 ][ t8 ];
              if ( gainFirstBridge + gainSecondBridge > maxGain ) {
                maxGain = gainFirstBridge + gainSecondBridge;
                bestTs = { t1, t2, t3, t4, t5, t6, t7, t8 };
              }
            }
          }
        }
        if ( maxGain > eps ) {
          doubleBridgeSwap_( tour, position, bestTs[ 0 ], bestTs[ 1 ], bestTs[ 2 ], bestTs[ 3 ], bestTs[ 4 ], bestTs[ 5 ], bestTs[ 6 ], bestTs[ 7 ] );
          changed = true;
        }
      }
    }
  }
  } // anonymous namespace

  vector< size_t > INLINE_ATTRIBUTE computeTour( const vector< vector< double > >& distances )
  {
    assert( distances.size() > 0 );
    srand( 1729 );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      assert( distances.size() == distances[ i ].size() );
      for ( size_t j = 0; j < distances.size(); ++j ) {
        assert( fabs( distances[ i ][ j ] - distances[ j ][ i ] ) < 1e-9 );
        if ( i != j ) {
          assert( distances[ i ][ j ] > 0.0 );
        }
        else {
          assert( distances[ i ][ i ] == 0.0 );
        }
      }
    }
    const vector< size_t > tourRand = getRandomTour_( distances );
    const vector< size_t > tourNN = getNearestNeighborTour_( distances );
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
      tour = tourGreedy;
      double start( clock() );
      compute3OptTour_( tour, distances, nearestNeighbors );
      computeDoubleBridgeTour_( tour, distances, nearestNeighbors );
      compute3OptTour_( tour, distances, nearestNeighbors );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "4-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
      assertIsTour_( tour, distances );
    }

    if ( false ) {
      compute5OptTour_( tour, distances, nearestNeighbors );
    }

    return tour;
  }
} // TravelingSalespersonProblemSolver
