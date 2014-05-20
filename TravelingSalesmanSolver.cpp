#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <limits>
//#include <numeric>
//#include <stdexcept>
#include <set>

/* HEADER */
#include "TravelingSalesmanSolver.h"

using namespace std;

namespace TravelingSalesmanSolver {
  void assertIsPath( const vector< size_t >& path, const vector< vector< double > >& distances )
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

  size_t nSteps;
  void reverse( vector< size_t >& path, const size_t i, const size_t j )
  {
    if ( j - i < path.size() - j + i ) {
      reverse( path.begin() + i, path.begin() + j + 1 );
      nSteps += ( j - i ) / 2;
    }
    else {
      size_t nnn = 0;
      vector< size_t >::iterator iIt = path.begin() + i + 1;
      vector< size_t >::iterator jIt = path.begin() + j;
      for ( ; iIt != path.begin() && jIt != path.end(); ++jIt ) {
        --iIt;
        swap( *iIt, *jIt );
        ++nSteps;
        ++nnn;
      }
      if ( jIt == path.end() ) {
        jIt = path.begin();
      }
      else if ( iIt == path.begin() ) {
        iIt = path.end() - 1;
      }
      for ( ; iIt > jIt; --iIt, ++jIt ) {
        swap( *iIt, *jIt );
        ++nSteps;
        ++nnn;
      }
      if ( fabs( int(nnn) - int( path.size() - j + i ) / 2 ) > 5 ) {
        cerr << "n: " << nnn << " path.size() " << path.size() << " i: " << i << ", j: " << j << ", " <<  int( path.size() - j + i ) / 2  << endl;
        assert( false );
      }
    }
  }

  double oneTreeLength( const vector< vector< double > >& distances )
  {
    // Compute minimum spanning tree of the vertices excluding the first
    vector< size_t > minimumSpanningTree;
    minimumSpanningTree.reserve( distances.size() );
    minimumSpanningTree.push_back( 1 );

    set< size_t > unusedVertices;
    for ( size_t i = 2; i < distances.size(); ++i ) {
      unusedVertices.insert( i );
    }

    double length = 0.0;
    while ( minimumSpanningTree.size() + 1 < distances.size() ) {
      size_t minIndex = 0;
      double minDistance = numeric_limits< double >::max();
      for ( size_t i = 0; i < minimumSpanningTree.size(); ++i ) {
        size_t mstIndex = minimumSpanningTree[ i ];
        for ( set< size_t >::const_iterator it = unusedVertices.begin(); it != unusedVertices.end(); ++it ) {
          double distance = distances[ mstIndex ][ *it ];
          if ( distance < minDistance ) {
            minIndex = *it;
            minDistance = distance;
          }
        }
      }
      minimumSpanningTree.push_back( minIndex );
      unusedVertices.erase( minIndex );
      length += minDistance;
    }
    double minElement = numeric_limits< double >::max();
    double secondMinElement = numeric_limits< double >::max();
    for ( size_t i = 0; i < distances.front().size(); ++i ) {
      double value = distances.front()[ i ];
      if ( value < secondMinElement ) {
        secondMinElement = value;
        if ( value < minElement ) {
          secondMinElement = minElement;
          minElement = value;
        }
      }
    }
    length += minElement + secondMinElement;
    return length;
  }

  double getLength( vector< size_t > path, const vector< vector< double > >& distances )
  {
    double distance = 0.0;
    for ( size_t i = 0; i + 1 < path.size(); ++i ) {
      distance += distances[ path[ i ] ][ path[ i + 1 ] ];
    }
    distance += distances[ path.back() ][ path.front() ];
    return distance;
  }

  struct Vertex {
    Vertex( size_t index ) : index_( index ) {}
    size_t getIndex() const { return index_; }
    const vector< size_t >& getChildren() const { return children_; }
    void addChild( size_t index ) { children_.push_back( index ); }
    private:
      size_t index_;
      vector< size_t > children_;
  };

  vector< size_t > computeMinimumSpanningTreePath( const vector< vector< double > >& distances )
  {
    // Compute minimum spanning tree
    vector< Vertex > nodes;
    nodes.reserve( distances.size() );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      nodes.push_back( Vertex( i ) );
    }

    vector< size_t > minimumSpanningTree;
    minimumSpanningTree.reserve( distances.size() );
    minimumSpanningTree.push_back( 0 );

    set< size_t > unusedVertices;
    for ( size_t i = 1; i < distances.size(); ++i ) {
      unusedVertices.insert( i );
    }

    double length = 0.0;
    while ( minimumSpanningTree.size() < distances.size() ) {
      size_t fromIndex = 0;
      size_t minIndex = 0;
      double minDistance = numeric_limits< double >::max();
      for ( size_t i = 0; i < minimumSpanningTree.size(); ++i ) {
        size_t mstIndex = minimumSpanningTree[ i ];
        for ( set< size_t >::const_iterator it = unusedVertices.begin(); it != unusedVertices.end(); ++it ) {
          double distance = distances[ mstIndex ][ *it ];
          if ( distance < minDistance ) {
            minIndex = *it;
            minDistance = distance;
            fromIndex = mstIndex;
          }
        }
      }
      minimumSpanningTree.push_back( minIndex );
      unusedVertices.erase( minIndex );
      length += minDistance;
      nodes[ fromIndex ].addChild( minIndex );
    }

    set< size_t > testSet( minimumSpanningTree.begin(), minimumSpanningTree.end() );
    assert( testSet.size() == distances.size() );

    vector< size_t > path;
    path.reserve( distances.size() );
    vector< Vertex* > stack;
    stack.push_back( &nodes.front() );
    while ( stack.size() > 0 ) {
      Vertex* v = stack.back();
      stack.pop_back();
      path.push_back( v->getIndex() );
      const vector< size_t >& children = v->getChildren();
      for ( size_t i = children.size(); i > 0; --i ) {
        stack.push_back( &nodes[ children[ i - 1 ] ] );
      }
    }

    return path;
  }

  vector< size_t > constructNearestNeighborPath( const vector< vector< double > >& distances )
  {
    vector< size_t > path;
    path.reserve( distances.size() );
    path.push_back( 0 );
    set< size_t > usedNodes;
    usedNodes.insert( 0 );
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
    assertIsPath( path, distances );
    return path;
  }

  bool update3Opt( const size_t i, const size_t j, const size_t k, vector< size_t >& path, const vector< vector< double > >& distances )
  {
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
    // subtract a little something to avoid numerical errors
    const double eps = 1e-9;
    const double removedDistance = distances[ pathIminus1 ][ pathI ] + distances[ pathJminus1 ][ pathJ ] + distances[ pathKminus1 ][ pathK ] - eps;

    if ( distances[ pathJminus1 ][ pathK ] + distances[ pathIminus1 ][ pathKminus1 ] + distances[ pathJ ][ pathI ] < removedDistance ) {
      vector< size_t > pathCopy( path );
      copy( pathCopy.begin() + i, pathCopy.begin() + j, path.begin() );
      size_t pathIdx = j - i;
      for ( size_t idx = k; idx < i + path.size(); ++idx, ++pathIdx ) {
        path[ pathIdx ] = pathCopy[ idx % path.size() ];
      }
      for ( size_t idx = kMinus1; idx >= j; --idx, ++pathIdx ) {
        path[ pathIdx ] = pathCopy[ idx ];
      }
      assert( pathIdx == path.size() );
      assert( getLength( path, distances ) < getLength( pathCopy, distances ) );
      return true;
    }

    if ( distances[ pathJminus1 ][ pathIminus1 ] + distances[ pathK ][ pathJ ] + distances[ pathKminus1 ][ pathI ] < removedDistance ) {
      vector< size_t > pathCopy( path );
      copy( pathCopy.begin() + i, pathCopy.begin() + j, path.begin() );
      size_t pathIdx = j - i;
      for ( size_t idx = i + path.size() - 1; idx >= k; --idx, ++pathIdx ) {
        path[ pathIdx ] = pathCopy[ idx % path.size() ];
      }
      copy( pathCopy.begin() + j, pathCopy.begin() + k, path.begin() + pathIdx );
      pathIdx += k - j;
      assert( pathIdx == path.size() );
      assert( getLength( path, distances ) < getLength( pathCopy, distances ) );
      return true;
    }

    if ( distances[ pathJminus1 ][ pathKminus1 ] + distances[ pathJ ][ pathIminus1 ] + distances[ pathK ][ pathI ] < removedDistance ) {
      vector< size_t > pathCopy( path );
      copy( pathCopy.begin() + i, pathCopy.begin() + j, path.begin() );
      size_t pathIdx = j - i;
      for ( size_t idx = kMinus1; idx >= j; --idx, ++pathIdx ) {
        path[ pathIdx ] = pathCopy[ idx ];
      }
      for ( size_t idx = i + path.size() - 1; idx >= k; --idx, ++pathIdx ) {
        path[ pathIdx ] = pathCopy[ idx % path.size() ];
      }
      assert( pathIdx == path.size() );
      assert( getLength( path, distances ) < getLength( pathCopy, distances ) );
      return true;
    }

    return false;
  }

  void compute3OptPathRandom( vector< size_t >& path, const vector< vector< double > >& distances )
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
        if ( update3Opt( randNums[ 0 ], randNums[ 1 ], randNums[ 2 ], path, distances ) ) {
          changed = true;
        }
      }
    }
  }

  void compute3OptPath( vector< size_t >& path, const vector< vector< double > >& distances )
  {
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) {
        for ( size_t j = i + 1; j < path.size(); ++j ) {
          for ( size_t k = j + 1; k < path.size(); ++k ) {
            if ( update3Opt( i, j, k, path, distances ) ) {
              changed = true;
              break;
            }
          }
        }
      }
    }
  }

  bool update2Opt( size_t i, size_t j, vector< size_t >& path, const vector< vector< double > >& distances )
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

  void compute2OptPathRandom( vector< size_t >& path, const vector< vector< double > >& distances )
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
        if ( update2Opt( randNums[ 0 ], randNums[ 1 ], path, distances ) ) {
          changed = true;
        }
      }
    }
  }

  void compute2OptPath( vector< size_t >& path, const vector< vector< double > >& distances )
  {
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) {
        for ( size_t j = i + 1; j < path.size(); ++j ) {
          if ( update2Opt( i, j, path, distances ) ) {
            changed = true;
            break;
          }
        }
      }
    }
  }


  /*
  void computeLinKernighanPath( vector< size_t >& path, const vector< vector< double > >& distances )
  {
    vector< size_t > t;
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) {
        t.clear();
        t.push_back( path[ i == 0 ? path.size() - 1 : i - 1 ] );
        t.push_back( path[ i ] );
        size_t ind = 1;
        for ( size_t j = i + 2; j < path.size(); ++j ) {
          if ( i < 2 && j + 1 == path.size() ) {
            continue;
          }
          if ( distances[ t[ 0 ] ][ t[ 1 ] ] < distances[ t[ 1 ] ][ path[ j ] ] ) {
            continue;
          }
          t.push_back( path[ j ] );
          ++ind;

          t.push_back( path[ j == 0 ? path.size() - 1 : j - 1 ] );
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

  vector< size_t > computePath( const vector< vector< double > >& distances )
  {
    nSteps = 0;
    assert( distances.size() > 0 );
    assert( distances.size() < RAND_MAX );
    srand( 1729 );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      assert( distances.size() == distances[ i ].size() );
    }
    vector< size_t > path( distances.size() );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      path[ i ] = i;
    }

    for ( size_t i = 0; i < distances.size(); ++i ) {
      size_t ind1 = rand() % distances.size();
      size_t ind2 = rand() % distances.size();
      swap( path[ ind1 ], path[ ind2 ] );
    }
    vector< size_t > pathNN = constructNearestNeighborPath( distances );
    vector< size_t > path2( pathNN );
    vector< size_t > path3( path );

    cerr << "Initial distance: " << getLength( path, distances ) << endl;
    cerr << "Nearest neighbor distance: " << getLength( pathNN, distances ) << endl;
    if ( true ) {
      cerr << "1-tree distance: ";
      double start( clock() );
      double oneTreeL = oneTreeLength( distances );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << oneTreeL << ", time: " << setprecision( 4 ) << time << endl;
    }

    if ( false ) {
      double start( clock() );
      vector< size_t > mstPath = computeMinimumSpanningTreePath( distances );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "mst path distance: " << getLength( mstPath, distances ) << ", time: " << setprecision( 4 ) << time << endl;
      assertIsPath( mstPath, distances );
    }

    if ( true ) {
      double start( clock() );
      compute2OptPath( path2, distances );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "2-opt path distance: " << getLength( path2, distances ) << ", time: " << setprecision( 4 ) << time << endl;
      assertIsPath( path2, distances );
    }

    if ( true ) {
      double start( clock() );
      compute3OptPath( path2, distances );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << "3-opt path distance: " << getLength( path2, distances ) << ", time: " << setprecision( 4 ) << time << endl;
      assertIsPath( path2, distances );
    }

    assertIsPath( path, distances );

    return path;
  }
}
