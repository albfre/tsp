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
  void assertIsPath_( const vector< size_t >& path,
                      const vector< vector< double > >& distances )
  {
    assert( path.size() == distances.size() );
    vector< size_t > pathCopy( path );
    sort( pathCopy.begin(), pathCopy.end() );
    for ( size_t i = 0; i < pathCopy.size(); ++i ) {
      if ( pathCopy[ i ] != i ) {
        cerr << "sort( pathCopy )[ i ]: " << pathCopy[ i ] << ", i: " << i << endl;
      }
      assert( pathCopy[ i ] == i );
    }
  }

  vector< size_t > getRandomPath_( const vector< vector< double > >& distances )
  {
    vector< size_t > path( distances.size() );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      path[ i ] = i;
    }

    for ( size_t i = 0; i < distances.size(); ++i ) {
      size_t ind1 = rand() % distances.size();
      size_t ind2 = rand() % distances.size();
      swap( path[ ind1 ], path[ ind2 ] );
    }
    return path;
  }

  /*
  void crossingReverse( vector< size_t >& path, const size_t i, const size_t j )
  {
    assert( i > j );
    vector< size_t >::iterator iIt = path.begin() + i;
    vector< size_t >::iterator jIt = path.begin() + j;
    while ( iIt != path.end() && jIt != path.begin() ) {
      --jIt;
      iter_swap( iIt, jIt );
      ++iIt;
    }
    if ( iIt == path.end() ) {
      reverse( path.begin(), jIt );
    }
    else {
      reverse( iIt, path.end() );
    }
  }

  void reverseShortest( vector< size_t >& path, const size_t i, const size_t j )
  {
    if ( j - i < path.size() + i - j ) {
      reverse( path.begin() + i, path.begin() + j );
    }
    else {
      vector< size_t >::iterator iIt = path.begin() + i + 1;
      vector< size_t >::iterator jIt = path.begin() + j;
      while ( iIt-- != path.begin() && jIt != path.end() ) {
        iter_swap( iIt, jIt );
        ++jIt;
      }
      if ( jIt == path.end() ) {
        reverse( path.begin(), iIt );
      }
      else {
        reverse( jIt, path.end() );
      }
    }
  }
  */

  vector< size_t > getNearestNeighborPath_( const vector< vector< double > >& distances )
  {
    vector< size_t > path;
    path.reserve( distances.size() );
    size_t startNode = 0; //rand() % distances.size();
    path.push_back( startNode );
    set< size_t > usedNodes;
    usedNodes.insert( startNode );
    while ( path.size() < distances.size() ) {
      size_t currentNode = path.back();
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
      path.push_back( minUnusedIndex );
      usedNodes.insert( minUnusedIndex );
    }
    assertIsPath_( path, distances );
    return path;
  }

  vector< size_t > getGreedyPath_( const vector< vector< double > >& distances )
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
    vector< size_t > path;
    for ( size_t i = 0; i < fragments.size(); ++i ) {
      if ( fragments[ i ].size() > 0 ) {
        path.swap( fragments[ i ] );
        break;
      }
    }

    assertIsPath_( path, distances );
    return path;
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
    // to construct Hamiltonian paths) with epsilon distance to all neighbors.
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
    for ( size_t i = 0; i + 2 < distances.size(); ++i ) {
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
    for ( size_t i = 0; i < distances[ 0 ].size(); ++i ) {
      double value = distances[ 0 ][ i ] + lagrangeMultipliers[ 0 ] + lagrangeMultipliers[ i ];
      if ( value < secondMinElement ) {
        secondMinElement = value;
        if ( value < minElement ) {
          secondMinElement = minElement;
          minElement = value;
        }
      }
    }
    vertexDegrees[ 0 ] = 2;

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
      double t = 2.0 * lambda * treeLength / denominator;

      for ( size_t j = 0; j < lagrangeMultipliers.size(); ++j ) {
        lagrangeMultipliers[ j ] += t * ( int( vertexDegrees[ j ] ) - 2 );
      }
      lambda = 1.0 / ( 20.0 + 10 * i );
    }

    return bestLength;
  }

  double getLength_( vector< size_t > path,
                     const vector< vector< double > >& distances )
  {
    double distance = distances[ path.back() ][ path.front() ];
    for ( size_t i = 0; i + 1 < path.size(); ++i ) {
      distance += distances[ path[ i ] ][ path[ i + 1 ] ];
    }
    return distance;
  }

  void update3OptIntervals_( vector< size_t >& path,
                             size_t i1Begin, size_t i1End, bool reverseI1,
                             size_t i2Begin, size_t i2End, bool reverseI2,
                             size_t i3Begin, size_t i3End, bool reverseI3,
                             const vector< vector< double > >& distances )
  {
    vector< size_t > pathCopy( path );
    size_t pathIdx = 0;
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
          iBegin += path.size();
        }
        for ( size_t idx = iBegin; idx > iEnd; --idx, ++pathIdx ) {
          path[ pathIdx ] = pathCopy[ idx % path.size() ];
        }
      }
      else {
        if ( iEnd < iBegin ) {
          iEnd += path.size();
        }
        for ( size_t idx = iBegin; idx < iEnd; ++idx, ++pathIdx ) {
          path[ pathIdx ] = pathCopy[ idx % path.size() ];
        }
      }
    }
    if ( pathIdx != path.size() ) {
      cerr << pathIdx << " " << path.size() << endl;
    }
    assert( pathIdx == path.size() );
    assert( getLength_( path, distances ) < getLength_( pathCopy, distances ) );
  }

  bool update3Opt_( size_t i,
                    size_t j,
                    size_t k,
                    vector< size_t >& path,
                    const vector< vector< double > >& distances )
  {
    vector< size_t > vertices( 3 );
    vertices[ 0 ] = i; vertices[ 1 ] = j; vertices[ 2 ] = k;
    sort( vertices.begin(), vertices.end() );
    i = vertices[ 0 ]; j = vertices[ 1 ]; k = vertices[ 2 ];
    assert( i < j && j < k );
    const size_t pathI = path[ i ];
    const size_t iMinus1 = ( i + path.size() - 1 ) % path.size();
    const size_t pathIminus1 = path[ iMinus1 ];
    const size_t pathJ = path[ j ];
    const size_t jMinus1 = ( j + path.size() - 1 ) % path.size();
    const size_t pathJminus1 = path[ jMinus1 ];
    const size_t pathK = path[ k ];
    const size_t kMinus1 = ( k + path.size() - 1 ) % path.size();
    const size_t pathKminus1 = path[ kMinus1 ];
    const double eps = 1e-9;
    const double removedDistance = distances[ pathIminus1 ][ pathI ] +
                                   distances[ pathJminus1 ][ pathJ ] +
                                   distances[ pathKminus1 ][ pathK ] - eps; // subtract a little something to avoid numerical errors
    vector< double > newDistances( 4 );
    newDistances[ 0 ] = distances[ pathI ][ pathJ ] + distances[ pathK ][ pathJminus1 ] + distances[ pathIminus1 ][ pathKminus1 ];
    newDistances[ 1 ] = distances[ pathJ ][ pathK ] + distances[ pathI ][ pathKminus1 ] + distances[ pathJminus1 ][ pathIminus1 ];
    newDistances[ 2 ] = distances[ pathK ][ pathI ] + distances[ pathJ ][ pathIminus1 ] + distances[ pathKminus1 ][ pathJminus1 ];
    newDistances[ 3 ] = distances[ pathI ][ pathKminus1 ] + distances[ pathJ ][ pathIminus1 ] + distances[ pathK ][ pathJminus1 ];
    size_t minIndex = min_element( newDistances.begin(), newDistances.end() ) - newDistances.begin();
    if ( newDistances[ minIndex ] < removedDistance ) {
      switch ( minIndex ) {
        case 0: reverse( path.begin() + i, path.begin() + j ); reverse( path.begin() + i, path.begin() + k ); break;
        case 1: reverse( path.begin() + i, path.begin() + j ); reverse( path.begin() + j, path.begin() + k ); break;
//        case 2: reverse( path.begin() + i, path.begin() + k ); crossingReverse( path, k, i ); break;
//        case 0: update3OptIntervals_( path, i, j, false, k, i, false, kMinus1, jMinus1, true, distances ); break;
//        case 1: update3OptIntervals_( path, i, j, false, iMinus1, kMinus1, true, j, k, false, distances ); break;
        case 2: update3OptIntervals_( path, i, j, false, kMinus1, jMinus1, true, iMinus1, kMinus1, true, distances ); break;
        case 3: update3OptIntervals_( path, i, j, false, k, i, false, j, k, false, distances ); break;
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
                   vector< size_t >& path,
                   const vector< vector< double > >& distances )
  {
    vector< size_t > vertices( 3 );
    vertices[ 0 ] = i; vertices[ 1 ] = j; vertices[ 2 ] = k;
    sort( vertices.begin(), vertices.end() );
    i = vertices[ 0 ]; j = vertices[ 1 ]; k = vertices[ 2 ];
    assert( i < j && j < k );
    const size_t pathI = path[ i ];
    const size_t iMinus1 = ( i + path.size() - 1 ) % path.size();
    const size_t pathIminus1 = path[ iMinus1 ];
    const size_t pathJ = path[ j ];
    const size_t jMinus1 = ( j + path.size() - 1 ) % path.size();
    const size_t pathJminus1 = path[ jMinus1 ];
    const size_t pathK = path[ k ];
    const size_t kMinus1 = ( k + path.size() - 1 ) % path.size();
    const size_t pathKminus1 = path[ kMinus1 ];
    const double eps = 1e-9;
    const double removedDistance = distances[ pathIminus1 ][ pathI ] +
                                   distances[ pathJminus1 ][ pathJ ] +
                                   distances[ pathKminus1 ][ pathK ] - eps; // subtract a little something to avoid numerical errors
    vector< double > newDistances( 4 );
    newDistances[ 0 ] = distances[ pathI ][ pathJ ] + distances[ pathK ][ pathJminus1 ] + distances[ pathIminus1 ][ pathKminus1 ];
    newDistances[ 1 ] = distances[ pathJ ][ pathK ] + distances[ pathI ][ pathKminus1 ] + distances[ pathJminus1 ][ pathIminus1 ];
    newDistances[ 2 ] = distances[ pathK ][ pathI ] + distances[ pathJ ][ pathIminus1 ] + distances[ pathKminus1 ][ pathJminus1 ];
    newDistances[ 3 ] = distances[ pathI ][ pathKminus1 ] + distances[ pathJ ][ pathIminus1 ] + distances[ pathK ][ pathJminus1 ];
    size_t minIndex = min_element( newDistances.begin(), newDistances.end() ) - newDistances.begin();
    return removedDistance - newDistances[ minIndex ];
  }

  void compute3OptPath2_( vector< size_t >& path,
                         const vector< vector< double > >& distances,
                         const vector< vector< size_t > >& nearestNeighbors )
  {
    vector< size_t > position( path.size() );
    for ( size_t i = 0; i < path.size(); ++i ) {
      position[ path[ i ] ] = i;
    }
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) {
        for ( vector< size_t >::const_iterator jIt = nearestNeighbors[ i ].begin(); jIt != nearestNeighbors[ i ].end(); ++jIt ) {
          size_t bestI = 0;
          size_t bestJ = 0;
          size_t bestK = 0;
          double bestGain = 0.0;
          size_t numGains = 0;

          size_t indexOfIInPath = position[ i ];
          size_t indexOfJInPath = position[ *jIt ];
          assert( indexOfJInPath != indexOfIInPath );
          for ( vector< size_t >::const_iterator kIt = nearestNeighbors[ *jIt ].begin(); kIt != nearestNeighbors[ *jIt ].end(); ++kIt ) {
            size_t indexOfKInPath = position[ *kIt ];
            assert( indexOfKInPath != indexOfJInPath );
            if ( indexOfKInPath == indexOfIInPath ) {
              continue;
            }
            double gain = getGain_( indexOfIInPath, indexOfJInPath, indexOfKInPath, path, distances ) ;
            if ( gain > 0 ) {
              ++numGains;
            }
            if ( gain > bestGain ) {
              bestGain = gain;
              bestI = indexOfIInPath;
              bestJ = indexOfJInPath;
              bestK = indexOfKInPath;
              changed = true;
            }
            if ( numGains >= 8 ) {
              break;
            }
          }
          if ( bestGain > 0.0 ) {
            assert( update3Opt_( bestI, bestJ, bestK, path, distances ) );
            for ( size_t pi = 0; pi < path.size(); ++pi ) {
              position[ path[ pi ] ] = pi;
            }
            break;
          }
        }
      }
    }
  }
  */

  void compute3OptPath_( vector< size_t >& path,
                         const vector< vector< double > >& distances,
                         const vector< vector< size_t > >& nearestNeighbors )
  {
    vector< size_t > position( path.size() );
    for ( size_t i = 0; i < path.size(); ++i ) {
      position[ path[ i ] ] = i;
    }
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) {
        for ( vector< size_t >::const_iterator jIt = nearestNeighbors[ i ].begin(); jIt != nearestNeighbors[ i ].end(); ++jIt ) {
          size_t indexOfIInPath = position[ i ];
          size_t indexOfJInPath = position[ *jIt ];
          assert( indexOfJInPath != indexOfIInPath );
          for ( vector< size_t >::const_iterator kIt = nearestNeighbors[ *jIt ].begin(); kIt != nearestNeighbors[ *jIt ].end(); ++kIt ) {
            size_t indexOfKInPath = position[ *kIt ];
            assert( indexOfKInPath != indexOfJInPath );
            if ( indexOfKInPath == indexOfIInPath ) {
              continue;
            }
            if ( update3Opt_( indexOfIInPath, indexOfJInPath, indexOfKInPath, path, distances ) ) {
              changed = true;
              for ( size_t i = 0; i < path.size(); ++i ) {
                position[ path[ i ] ] = i;
              }
              break;
            }
          }
        }
      }
    }
  }

  bool update4Opt_( size_t i,
                    size_t j,
                    size_t k,
                    size_t l,
                    vector< size_t >& path,
                    const vector< vector< double > >& distances )
  {
    vector< size_t > vertices( 4 );
    vertices[ 0 ] = i; vertices[ 1 ] = j; vertices[ 2 ] = k; vertices[ 3 ] = l;
    sort( vertices.begin(), vertices.end() );
    i = vertices[ 0 ]; j = vertices[ 1 ]; k = vertices[ 2 ]; l = vertices[ 3 ];
    assert( i < j && j < k );
    const size_t pathI = path[ i ];
    const size_t iMinus1 = ( i + path.size() - 1 ) % path.size();
    const size_t pathIminus1 = path[ iMinus1 ];
    const size_t pathJ = path[ j ];
    const size_t jMinus1 = ( j + path.size() - 1 ) % path.size();
    const size_t pathJminus1 = path[ jMinus1 ];
    const size_t pathK = path[ k ];
    const size_t kMinus1 = ( k + path.size() - 1 ) % path.size();
    const size_t pathKminus1 = path[ kMinus1 ];
    const size_t pathL = path[ l ];
    const size_t lMinus1 = ( l + path.size() - 1 ) % path.size();
    const size_t pathLminus1 = path[ lMinus1 ];
    const double eps = 1e-9;
    const double removedDistance = distances[ pathIminus1 ][ pathI ] +
                                   distances[ pathJminus1 ][ pathJ ] +
                                   distances[ pathKminus1 ][ pathK ] +
                                   distances[ pathLminus1 ][ pathL ] - eps; // subtract a little something to avoid numerical errors
    double newDistance = distances[ pathI ][ pathKminus1 ] + distances[ pathJ ][ pathLminus1 ] + distances[ pathK ][ pathIminus1 ] + distances[ pathL ][ pathJminus1 ];
    if ( newDistance < removedDistance ) {
      vector< size_t > pathCopy( path );
      // {i, j-1}
      copy( pathCopy.begin() + i, pathCopy.begin() + j, path.begin() );
      size_t pathIndex = j - i;

      // {l, i-1}
      for ( size_t ind = l; ind < i + path.size(); ++ind, ++pathIndex ) {
        path[ pathIndex ] = pathCopy[ ind % path.size() ];
      }

      // {k, l-1}
      copy( pathCopy.begin() + k, pathCopy.begin() + l, path.begin() + pathIndex );
      pathIndex += l - k;

      // {j, k-1}
      copy( pathCopy.begin() + j, pathCopy.begin() + k, path.begin() + pathIndex );
      pathIndex += k - j;
      if ( pathIndex != path.size() ) {
        cerr << "i: " << i << ", j: " << j << ", k: " << k << ", l: " << l << ", pathIndex: " << pathIndex << ", path.size(): " << path.size() << endl;
      }

      assert( pathIndex == path.size() );
      return true;
    }
    return false;
  }

  void compute4OptPath_( vector< size_t >& path,
                         const vector< vector< double > >& distances,
                         const vector< vector< size_t > >& nearestNeighbors )
  {
    vector< size_t > position( path.size() );
    for ( size_t i = 0; i < path.size(); ++i ) {
      position[ path[ i ] ] = i;
    }
    bool changed = true;
    size_t maxNeighbors = 10;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) {
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
              size_t indexOfIInPath = position[ i ];
              size_t indexOfJInPath = position[ *jIt ];
              size_t indexOfKInPath = position[ *kIt ];
              size_t indexOfLInPath = position[ *lIt ];
              assert( indexOfJInPath != indexOfIInPath );
              assert( indexOfKInPath != indexOfJInPath );
              assert( indexOfLInPath != indexOfKInPath );
              if ( indexOfKInPath == indexOfIInPath || indexOfLInPath == indexOfIInPath || indexOfLInPath == indexOfJInPath ) {
                continue;
              }
              if ( update4Opt_( indexOfIInPath, indexOfJInPath, indexOfKInPath, indexOfLInPath, path, distances ) ) {
                changed = true;
                for ( size_t i = 0; i < path.size(); ++i ) {
                  position[ path[ i ] ] = i;
                }
                break;
              }
            }
          }
        }
      }
    }
  }


  bool update2Opt_( size_t i,
                    size_t j,
                    vector< size_t >& path,
                    const vector< vector< double > >& distances )
  {
    assert( i < j );
    const double eps = 1e-9;
    const size_t t0 = path[ i == 0 ? path.size() - 1 : i - 1 ];
    const size_t t1 = path[ i ];
    const size_t t2 = path[ j == 0 ? path.size() - 1 : j - 1 ];
    const size_t t3 = path[ j % path.size() ];
    if ( distances[ t0 ][ t2 ] + distances[ t1 ][ t3 ] <
         distances[ t0 ][ t1 ] + distances[ t2 ][ t3 ] - eps ) {
      reverse( path.begin() + i, path.begin() + j );
      return true;
    }
    return false;
  }

  double getGain_( size_t i,
                   size_t j,
                   vector< size_t >& path,
                   const vector< vector< double > >& distances )
  {
    assert( i < j );
    const double eps = 1e-9;
    const size_t t0 = path[ i == 0 ? path.size() - 1 : i - 1 ];
    const size_t t1 = path[ i ];
    const size_t t2 = path[ j == 0 ? path.size() - 1 : j - 1 ];
    const size_t t3 = path[ j % path.size() ];
    return ( distances[ t0 ][ t1 ] + distances[ t2 ][ t3 ] - ( distances[ t0 ][ t2 ] + distances[ t1 ][ t3 ] ) ) - eps;
  }

  void compute2OptPath_( vector< size_t >& path,
                         const vector< vector< double > >& distances )
  {
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) {
        size_t bestI = (size_t)-1;
        size_t bestJ = (size_t)-1;
        double bestGain = 0.0;
        // Cutting the edges of two neighboring nodes does not lead to a new path. Hence start from i + 2.
        for ( size_t j = i + 2; j < path.size(); ++j ) {
          double gain = getGain_( i, j, path, distances ) ;
          if ( gain > bestGain ) {
            bestI = i;
            bestJ = j;
            bestGain = gain;
            changed = true;
          }
        }
        if ( bestGain > 0.0 ) {
          assert( update2Opt_( bestI, bestJ, path, distances ) );
        }
      }
    }
    assertIsPath_( path, distances );
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

  void make2OptMove_( size_t it0, size_t it1, size_t it2, size_t it3, vector< size_t >& path )
  {
    cerr << it0 << " " << it1 << " " << it2 << " " << it3 << endl;
    assert( fabs( fabs( float( it0 ) - float( it1 ) ) - 1.0 ) < 1e-6 || fabs( fabs( float( it0 ) - float( it1 ) ) + 1.0 - float( path.size() ) ) < 1e-6  );
    assert( fabs( fabs( float( it2 ) - float( it3 ) ) - 1.0 ) < 1e-6 || fabs( fabs( float( it2 ) - float( it3 ) ) + 1.0 - float( path.size() ) ) < 1e-6  );
    size_t ind1 = max( it0, it1 );
    size_t ind2 = max( it2, it3 );
    if ( ind1 > ind2 ) {
      swap( ind1, ind2 );
    }
    reverse( path.begin() + ind1, path.begin() + ind2 );
  }

  void computeLinKernighanPath_( vector< size_t >& path,
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
      for ( size_t i = 1; i < path.size(); ++i ) {
        vector< pair< size_t, size_t > > T;
        for ( size_t j = 0; j < path.size(); ++j ) {
          T.push_back( make_pair( path[ j ], path[ ( j + 1 ) % path.size() ] ) );
        }
        size_t ind = 0;
        x.clear();
        y.clear();
        t.clear();
        t.push_back( path[ i == 0 ? path.size() - 1 : i - 1 ] );
        t.push_back( path[ i ] );
        // Step 3. Choose x_0 = ( t_0, t_1 ) in T
        x.push_back( make_pair( t[ 0 ], t[ 1 ] ) );
        for ( size_t j = i + 2; j < path.size(); ++j ) {
          // Step 4. Choose y_0 = ( t_1, t_2 ) not in T such that G > 0
          if ( i < 2 && j + 1 == path.size() ) {
            continue;
          }
          if ( contains_( T, make_pair( t[ 1 ], path[ j ] ) ) ) {
            // y is in T
            continue;
          }
          double G0 = distances[ t[ 0 ] ][ t[ 1 ] ] - distances[ t[ 1 ] ][ path[ j ] ];
          if ( G0 <= eps ) {
            continue;
          }
          double G = G0;

          // Found y not in T with positive gain
          t.push_back( path[ j ] );
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
            size_t tNext = path[ ( ( find( path.begin(), path.end(), t.back() ) - path.begin() ) + path.size() - 1 ) % path.size() ]; // element in path previous to t_(2i)
            assert( nextXExists ); // condition (b): x_i != y_s for all s < i
            x.push_back( make_pair( t.back(), tNext ) ); // Add x_i = (t_(2i), t_(2i+1))
            t.push_back( tNext ); // Add t_(2i+1)

            cerr << ind << ". x size: " << x.size() << ", G: " << G << ", return: " << G + distances[ x.back().first ][ x.back().second ] - distances[ t.back() ][ t.front() ] << endl;
            if ( G + distances[ x.back().first ][ x.back().second ] - distances[ t.back() ][ t.front() ] > eps ) {
              y.push_back( make_pair( t.back(), t.front() ) );
              changed = true;
              // Take tour
              vector< size_t > pathCopy( path );
              assert( t.size() % 2 == 0 );
              assert( x.size() == y.size() );
              cerr << "Take tour" << endl;
              cerr << "x.size() " << x.size() << endl;
              cerr << "G: " << G << " dist: " << distances[ t.back() ][ t.front() ] << endl;
              cerr << "path: ";
              for ( size_t k = 0; k < path.size(); ++k ) {
                cerr << path[ k ] << " ";
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
                size_t t0 = find( path.begin(), path.end(), x[ 0 ].first ) - path.begin();
                size_t t1 = find( path.begin(), path.end(), x[ k - 1 ].second ) - path.begin();
                size_t t2 = find( path.begin(), path.end(), x[ k ].first ) - path.begin();
                size_t t3 = find( path.begin(), path.end(), x[ k ].second ) - path.begin();
                cerr << "2opt: (" << x[ 0 ].first << ", " << x[ k - 1].second << "), (" << x[ k ].first << ", " << x[ k ].second << ")" << endl;
                make2OptMove_( t0, t1, t2, t3, path );
              }
              assert( path != pathCopy );

              if ( getLength_( path, distances ) >= getLength_( pathCopy, distances ) ) {
                cerr << setprecision( 28 ) << getLength_( path, distances ) << " " << getLength_( pathCopy, distances ) << endl;
              }
              assert( getLength_( path, distances ) < getLength_( pathCopy, distances ) );
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
            for ( size_t k = startK; k < path.size(); ++k ) {
              if ( y.size() == 1 ) {
                ++maxKForIndEquals2;
              }
              if ( k == i || k == j || t.back() == path[ k ] ) {
                continue;
              }
              pair< size_t, size_t > yCandidate;
              yCandidate = make_pair( t.back(), path[ k ] );

              // y is not in T
              if ( contains_( T, yCandidate ) ) {
                continue;
              }
              // (a) G_i > 0
              double gain = distances[ x.back().first ][ x.back().second ] - distances[ t.back() ][ path[ k ] ];
              if ( G + gain <= eps ) {
                continue;
              }
              // (b) y_i != x_s for all s <= i
              if ( contains_( x, yCandidate ) ) {
                continue;
              }
              // (c) x_(i+1) exists
              size_t tNextCandidate = path[ ( ( find( path.begin(), path.end(), path[ k ] ) - path.begin() ) + path.size() - 1 ) % path.size() ]; // element in path previous to t_(2i+2) candidate
              pair< size_t, size_t > xCandidate;
              xCandidate = make_pair( path[ k ], tNextCandidate );
              nextXExists = !contains_( x, xCandidate ) &&
                            !contains_( y, xCandidate ) &&
                            tNextCandidate != t.back();
              if ( !nextXExists ) {
                continue;
              }
              G += gain;
              y.push_back( yCandidate );
              t.push_back( path[ k ] );

              // Found y, goto step 5
              nextYExists = true;
              break;
            }

            if ( !nextYExists && maxKForIndEquals2 < path.size() ) {
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
} // anonymous namespace

namespace TravelingSalespersonProblemSolver {
  vector< size_t > computePath( const vector< vector< double > >& distances )
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
    vector< size_t > pathRand = getRandomPath_( distances );
    vector< size_t > pathNN = getNearestNeighborPath_( distances );
    vector< size_t > pathGreedy = getGreedyPath_( distances );
    vector< size_t > path( pathGreedy );
    vector< vector< size_t > > nearestNeighbors = computeNearestNeighbors_( distances, 20 );

    cerr << setprecision( 7 );
    cerr << "Initial distance: " << getLength_( path, distances ) << endl;
    cerr << "Nearest neighbor distance: " << getLength_( pathNN, distances ) << endl;
    cerr << "Greedy distance: " << getLength_( pathGreedy, distances ) << endl;

    if ( true ) {
      cerr << "1-tree distance: ";
      double start( clock() );
      double lowerBound = getHeldKarpLowerBound_( distances );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << lowerBound << ", time: " << time << endl;
    }

    if ( true ) {
      double start( clock() );
      compute2OptPath_( path, distances );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "2-opt path distance: " << getLength_( path, distances ) << ", time: " << time << endl;
      assertIsPath_( path, distances );
    }

    if ( true ) {
      double start( clock() );
      compute3OptPath_( path, distances, nearestNeighbors );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "3-opt path distance: " << getLength_( path, distances ) << ", time: " << time << endl;
      assertIsPath_( path, distances );
    }

    if ( true ) {
      double start( clock() );
//      computeLinKernighanPath_( path, distances );
      compute4OptPath_( path, distances, nearestNeighbors );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "4-opt path distance: " << getLength_( path, distances ) << ", time: " << time << endl;
      assertIsPath_( path, distances );
    }

    if ( false ) {
      computeLinKernighanPath_( path, distances );
    }

    return path;
  }
}
