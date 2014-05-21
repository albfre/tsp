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
#include "TravelingSalesmanSolver.h"

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

  vector< vector< size_t > > computeNearestNeighbors_( const vector< vector< double > >& distances,
                                                       size_t numberOfNeighbors )
  {
    assert( numberOfNeighbors < distances.size() );
    vector< vector< size_t > > nearestNeighbors( distances.size(), vector< size_t >( numberOfNeighbors ) );
    vector< pair< double, size_t > > tmpNeighbors( distances.size() );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      for ( size_t j = 0; j < distances[ i ].size(); ++j ) {
        tmpNeighbors[ j ] = make_pair( distances[ i ][ j ], j );
      }
      sort( tmpNeighbors.begin(), tmpNeighbors.end() );
      // Start from 1 to avoid adding self as neighbor
      for ( size_t j = 0; j < numberOfNeighbors; ++j ) {
        nearestNeighbors[ i ][ j ] = tmpNeighbors[ j + 1 ].second;
      }
    }
    for ( size_t i = 0; i < nearestNeighbors.size(); ++i ) {
      assert( find( nearestNeighbors[ i ].begin(), nearestNeighbors[ i ].end(), i ) == nearestNeighbors[ i ].end() );
    }
    return nearestNeighbors;
  }

  struct Vertex {
    Vertex( size_t index ) : index_( index ), degree_( 0 ) {}
    size_t getIndex() const { return index_; }
    const vector< size_t >& getChildren() const { return children_; }
    void addChild( size_t index ) { children_.push_back( index ); }
    void increaseDegree( size_t increment = 1) { degree_ += increment; }
    size_t getDegree() const { return degree_; }
    private:
      size_t index_;
      size_t degree_;
      vector< size_t > children_;
  };

  double get1Tree_( vector< Vertex >& nodes,
                    const vector< vector< double > >& distances,
                    const vector< double >& lagrangeMultipliers )
  {
    // Compute minimum spanning tree of the vertices excluding the first
    nodes.clear();
    for ( size_t i = 0; i < distances.size(); ++i ) {
      nodes.push_back( Vertex( i ) );
    }

    vector< size_t > minimumSpanningTree;
    minimumSpanningTree.reserve( distances.size() );
    minimumSpanningTree.push_back( 1 );
    vector< size_t > unusedVertices;
    for ( size_t i = 2; i < distances.size(); ++i ) {
      unusedVertices.push_back( i );
    }

    double length = 0.0;
    while ( minimumSpanningTree.size() + 1 < distances.size() ) {
      size_t fromIndex = 0;
      vector< size_t >::iterator minIt = unusedVertices.begin();
      double minDistance = numeric_limits< double >::max();
      for ( size_t i = 0; i < minimumSpanningTree.size(); ++i ) {
        size_t mstIndex = minimumSpanningTree[ i ];
        for ( vector< size_t >::iterator it = unusedVertices.begin(); it != unusedVertices.end(); ++it ) {
          double distance = distances[ mstIndex ][ *it ] + lagrangeMultipliers[ mstIndex ] + lagrangeMultipliers[ *it ];
          if ( distance < minDistance ) {
            minDistance = distance;
            fromIndex = mstIndex;
            minIt = it;
          }
        }
      }
      minimumSpanningTree.push_back( *minIt );
      length += minDistance;
      nodes[ fromIndex ].addChild( *minIt );
      nodes[ fromIndex ].increaseDegree();
      nodes[ *minIt ].increaseDegree();
      unusedVertices.erase( minIt );
    }

    double minElement = numeric_limits< double >::max();
    double secondMinElement = numeric_limits< double >::max();
    size_t minIndex = (size_t)-1;
    size_t secondMinIndex = (size_t)-1;
    for ( size_t i = 0; i < distances[ 0 ].size(); ++i ) {
      double value = distances[ 0 ][ i ] + lagrangeMultipliers[ 0 ] + lagrangeMultipliers[ i ];
      if ( value < secondMinElement ) {
        secondMinElement = value;
        secondMinIndex = i;
        if ( value < minElement ) {
          secondMinElement = minElement;
          secondMinIndex = minIndex;
          minElement = value;
          minIndex = i;
        }
      }
    }
    nodes[ 0 ].addChild( minIndex );
    nodes[ 0 ].addChild( secondMinIndex );
    nodes[ 0 ].increaseDegree( 2 );

    length += minElement + secondMinElement;
    return length - 2.0 * accumulate( lagrangeMultipliers.begin(), lagrangeMultipliers.end(), 0.0 );
  }

  double getHeldKarpLowerBound_( const vector< vector< double > >& distances )
  {
    double bestLength = numeric_limits< double >::min();
    double length = numeric_limits< double >::min();
    vector< double > lagrangeMultipliers( distances.size() );
    double delta = 3e-3;
    for ( size_t i = 0; i < 10; ++i ) {
      vector< Vertex > nodes;
      length = get1Tree_( nodes, distances, lagrangeMultipliers );
      bestLength = max( bestLength, length );
      for ( size_t j = 0; j < lagrangeMultipliers.size(); ++j ) {
        lagrangeMultipliers[ j ] += ( int( nodes[ j ].getDegree() ) - 2 ) * delta;
      }
      delta *= 0.95;
    }

    return bestLength;
  }

  double getLength_( vector< size_t > path,
                     const vector< vector< double > >& distances )
  {
    double distance = 0.0;
    for ( size_t i = 0; i + 1 < path.size(); ++i ) {
      distance += distances[ path[ i ] ][ path[ i + 1 ] ];
    }
    distance += distances[ path.back() ][ path.front() ];
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
      size_t iBegin;
      size_t iEnd;
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
    vector< size_t > nums( 3 );
    nums[ 0 ] = i; nums[ 1 ] = j; nums[ 2 ] = k;
    sort( nums.begin(), nums.end() );
    i = nums[ 0 ]; j = nums[ 1 ]; k = nums[ 2 ];
    assert( i < j && j < k );
    const size_t pathI = path[ i ];
    const size_t iMinus1 = i == 0 ? path.size() - 1 : i - 1;
    const size_t pathIminus1 = path[ iMinus1 ];
    const size_t pathJ = path[ j ];
    const size_t jMinus1 = j == 0 ? path.size() - 1 : j - 1;
    const size_t pathJminus1 = path[ jMinus1 ];
    const size_t pathK = path[ k ];
    const size_t kMinus1 = k == 0 ? path.size() - 1 : k - 1;
    const size_t pathKminus1 = path[ kMinus1 ];
    const double eps = 1e-9;
    const double removedDistance = distances[ pathIminus1 ][ pathI ] +
                                   distances[ pathJminus1 ][ pathJ ] +
                                   distances[ pathKminus1 ][ pathK ] - eps; // subtract a little something to avoid numerical errors
    size_t bestIndex = (size_t)-1;
    double bestDistance = numeric_limits< double >::max();
    for ( size_t idx = 0; idx < 3; ++idx ) {
      double newDistance = numeric_limits< double >::max();
      switch ( idx ) {
        case 0: newDistance = distances[ pathJminus1 ][ pathK ] + distances[ pathIminus1 ][ pathKminus1 ] + distances[ pathJ ][ pathI ]; break;
        case 1: newDistance = distances[ pathJminus1 ][ pathIminus1 ] + distances[ pathK ][ pathJ ] + distances[ pathKminus1 ][ pathI ]; break;
        case 2: newDistance = distances[ pathJminus1 ][ pathKminus1 ] + distances[ pathJ ][ pathIminus1 ] + distances[ pathK ][ pathI ]; break;
        default: assert( false );
      }
      if ( newDistance < removedDistance && newDistance < bestDistance ) {
        bestIndex = idx;
        bestDistance = newDistance;
      }
    }
    switch( bestIndex ) {
      case 0: update3OptIntervals_( path, i, j, false, k, i, false, k - 1, j - 1, true, distances ); return true;
      case 1: update3OptIntervals_( path, i, j, false, i + path.size() - 1, k - 1, true, j, k, false, distances ); return true;
      case 2: update3OptIntervals_( path, i, j, false, k - 1, j - 1, true, i + path.size() - 1, k - 1, true, distances ); return true;
      default: return false;
    }
  }

  void compute3OptPathRandom_( vector< size_t >& path,
                               const vector< vector< double > >& distances )
  {
    bool changed = true;
    size_t outerIter = 0;
    while ( changed && outerIter < 100 ) {
      changed = false;
      ++outerIter;
      vector< size_t > randNums( 3 );
      for ( size_t iter = 0; iter < 10000; ++iter ) {
        randNums[ 0 ] = rand() % distances.size();
        randNums[ 1 ] = rand() % distances.size();
        randNums[ 2 ] = rand() % distances.size();
        sort( randNums.begin(), randNums.end() );
        if ( randNums[ 0 ] == randNums[ 1 ] || randNums[ 1 ] == randNums[ 2 ] ) {
          continue;
        }
        if ( update3Opt_( randNums[ 0 ], randNums[ 1 ], randNums[ 2 ], path, distances ) ) {
          changed = true;
        }
      }
    }
  }

  void compute3OptPath_( vector< size_t >& path,
                         const vector< vector< double > >& distances )
  {
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) {
        for ( size_t j = i + 1; j < path.size(); ++j ) {
          for ( size_t k = j + 1; k < path.size(); ++k ) {
            if ( update3Opt_( i, j, k, path, distances ) ) {
              changed = true;
              break;
            }
          }
        }
      }
    }
  }

  void compute3OptPath_( vector< size_t >& path,
                         const vector< vector< double > >& distances,
                         const vector< vector< size_t > >& nearestNeighbors )
  {
restart3opt:
    for ( size_t i = 0; i < path.size(); ++i ) {
      size_t indexOfIInPath = find( path.begin(), path.end(), i ) - path.begin();
      for ( vector< size_t >::const_iterator jIt = nearestNeighbors[ i ].begin(); jIt != nearestNeighbors[ i ].end(); ++jIt ) {
        size_t indexOfJInPath = find( path.begin(), path.end(), *jIt ) - path.begin();
        assert( indexOfJInPath != indexOfIInPath );
        for ( vector< size_t >::const_iterator kIt = nearestNeighbors[ *jIt ].begin(); kIt != nearestNeighbors[ *jIt ].end(); ++kIt ) {
          size_t indexOfKInPath = find( path.begin(), path.end(), *kIt ) - path.begin();
          assert( indexOfKInPath != indexOfJInPath );
          if ( indexOfKInPath == indexOfIInPath ) {
            continue;
          }
          if ( update3Opt_( indexOfIInPath, indexOfJInPath, indexOfKInPath, path, distances ) ) {
            goto restart3opt;
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
    if ( i < 2 && j + 1 == path.size() ) {
      return false;
    }
    const double eps = 1e-9;
    const size_t pathI = path[ i ];
    const size_t pathJ = path[ j ];
    const size_t pathJplus1 = path[ j + 1 == path.size() ? 0 : j + 1 ];
    const size_t pathIminus1 = path[ i == 0 ? path.size() - 1 : i - 1 ];
    if ( distances[ pathIminus1 ][ pathJ ] + distances[ pathI ][ pathJplus1 ] <
         distances[ pathIminus1 ][ pathI ] + distances[ pathJ ][ pathJplus1 ] - eps ) {
      reverse( path.begin() + i, path.begin() + j + 1 );
      return true;
    }
    return false;
  }

  void compute2OptPathRandom_( vector< size_t >& path,
                               const vector< vector< double > >& distances )
  {
    bool changed = true;
    size_t outerIter = 0;
    while ( changed && outerIter < 100 ) {
      changed = false;
      ++outerIter;
      vector< size_t > randNums( 2 );
      for ( size_t iter = 0; iter < 10000; ++iter ) {
        randNums[ 0 ] = rand() % distances.size();
        randNums[ 1 ] = rand() % distances.size();
        if ( randNums[ 0 ] == randNums[ 1 ] ) {
          continue;
        }
        sort( randNums.begin(), randNums.end() );
        if ( update2Opt_( randNums[ 0 ], randNums[ 1 ], path, distances ) ) {
          changed = true;
        }
      }
    }
  }

  void compute2OptPath_( vector< size_t >& path,
                         const vector< vector< double > >& distances )
  {
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) {
        for ( size_t j = i + 1; j < path.size(); ++j ) {
          if ( update2Opt_( i, j, path, distances ) ) {
            changed = true;
            break;
          }
        }
      }
    }
  }

  /*
  void computeLinKernighanPath_( vector< size_t >& path,
                                 const vector< vector< double > >& distances )
  {
    bool changed = true;
    vector< pair< size_t, size_t > > x;
    vector< pair< size_t, size_t > > y;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) { // Step 2
        x.clear();
        y.clear();
        x.push_back( make_pair( path[ i == 0 ? path.size() - 1 : i - 1 ], path[ i ] ) ); // Step 3
        size_t ind = 0;
        double G = 0.0;
        for ( size_t j = i + 2; j < path.size(); ++j ) {
          if ( i < 2 && j + 1 == path.size() ) {
            continue;
          }
          double gain = distances[ x[ ind ].first ][ x[ ind ].second ] - distances[ x[ ind ].second ][ path[ j ] ];
          if ( gain < 0.0 ) {
            continue;
          }
          G = gain;
          y.push_back( make_pair( x[ ind ].second, path[ j ] ) ); // Step 4
          ++ind; // Step 5

          x.push_back( make_pair( y[ ind - 1 ].second, path[ j - 1 ] ) )
          const double deltaDistance =
            distances[ t[ 1 ] ][ t[ 2 ] ] + distances[ t[ 3 ] ][ t[ 1 ] ] -
            distances[ t[ 0 ] ][ t[ 1 ] ] - distances[ t[ 2 ] ][ t[ 3 ] ];
          if ( deltaDistance < 0.0 ) {
            reverse( path.begin() + i, path.begin() + j + 1 );
            changed = true;
            break;
          }

          for ( size_t k = 0; k < path.size(); ++k ) {

          }
        }
        if ( !changed ) { // Step 12
        }

          const size_t pathI = path[ i ];
          const size_t pathJ = path[ j ];
          const size_t pathJplus1 = path[ j + 1 == path.size() ? 0 : j + 1 ];
          const size_t pathIminus1 = path[ i == 0 ? path.size() - 1 : i - 1 ];
          const double deltaDistance =
            distances[ pathIminus1 ][ pathJ ] + distances[ pathI ][ pathJplus1 ] -
            distances[ pathIminus1 ][ pathI ] - distances[ pathJ ][ pathJplus1 ];
          if ( deltaDistance < 0.0 ) {
            reverse( path.begin() + i, path.begin() + j + 1 );
            changed = true;
            break;

          }
        }
      }
    }
  }
  */
} // anonymous namespace

namespace TravelingSalesmanSolver {
  vector< size_t > computePath( const vector< vector< double > >& distances )
  {
    assert( distances.size() > 0 );
    srand( 1729 );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      assert( distances.size() == distances[ i ].size() );
    }
    vector< size_t > pathRand = getRandomPath_( distances );
    vector< size_t > pathNN = getNearestNeighborPath_( distances );
    vector< size_t > path( pathNN );
    vector< vector< size_t > > nearestNeighbors = computeNearestNeighbors_( distances, 10 );

    cerr << "Initial distance: " << getLength_( path, distances ) << endl;
    cerr << "Nearest neighbor distance: " << getLength_( pathNN, distances ) << endl;
    if ( true ) {
      cerr << "1-tree distance: ";
      double start( clock() );
      double lowerBound = getHeldKarpLowerBound_( distances );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << lowerBound << ", time: " << setprecision( 4 ) << time << endl;
    }

    if ( true ) {
      double start( clock() );
      compute2OptPath_( path, distances );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "2-opt path distance: " << getLength_( path, distances ) << ", time: " << setprecision( 4 ) << time << endl;
      assertIsPath_( path, distances );
    }

    if ( true ) {
      double start( clock() );
      compute3OptPath_( path, distances, nearestNeighbors );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "3-opt path distance: " << getLength_( path, distances ) << ", time: " << setprecision( 4 ) << time << endl;
      assertIsPath_( path, distances );
    }

    if ( false ) {
      compute2OptPathRandom_( path, distances );
      compute3OptPathRandom_( path, distances );
    }

    return path;
  }
}
