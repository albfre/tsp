#include <algorithm>
#include <assert.h>
//#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
//#include <numeric>
//#include <stdexcept>
#include <set>

/* HEADER */
#include "TravelingSalesmanSolver.h"

using namespace std;

namespace TravelingSalesmanSolver {
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

  void compute2OptPath( vector< size_t >& path, const vector< vector< double > >& distances )
  {
    assert( path.size() == distances.size() );
    double bestDistance = getLength( path, distances );
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) {
        for ( size_t j = i + 1; j < path.size(); ++j ) {
          swap( path[ i ], path[ j ] );
          double distance = getLength( path, distances );
          if ( distance > bestDistance ) {
            swap( path[ i ], path[ j ] );
          }
          else {
            bestDistance = distance;
            changed = true;
            break;
          }
        }
      }
    }
  }

  void compute3OptPath( vector< size_t >& path, const vector< vector< double > >& distances )
  {
    assert( path.size() == distances.size() );
    double bestDistance = getLength( path, distances );
    for ( size_t i = 0; i < path.size(); ++i ) {
      for ( size_t j = i + 1; j < path.size(); ++j ) {
        for ( size_t k = j + 1; k < path.size(); ++k ) {
          swap( path[ i ], path[ j ] );
          swap( path[ j ], path[ k ] );
          double distance = getLength( path, distances );
          if ( distance > bestDistance ) {
            swap( path[ i ], path[ j ] );
          }
          else {
            bestDistance = distance;
          }
        }
      }
    }
  }

  void compute2Reverse( vector< size_t >& path, const vector< vector< double > >& distances )
  {
    double bestDistance = getLength( path, distances );
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) {
        for ( size_t j = i + 1; j < path.size(); ++j ) {
          reverse( path.begin() + i, path.begin() + j + 1 );
          double distance = getLength( path, distances );
          if ( distance < bestDistance ) {
            bestDistance = distance;
            changed = true;
            break;
          }
          else {
            reverse( path.begin() + i, path.begin() + j + 1 );
          }
        }
      }
    }
  }

  void compute2ReverseRandom( vector< size_t >& path, const vector< vector< double > >& distances )
  {
    double bestDistance = getLength( path, distances );
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
        reverse( path.begin() + randNums[ 0 ], path.begin() + randNums[ 1 ] + 1 );
        double distance = getLength( path, distances );
        if ( distance < bestDistance ) {
          bestDistance = distance;
          changed = true;
        }
        else {
          reverse( path.begin() + randNums[ 0 ], path.begin() + randNums[ 1 ] + 1 );
        }
      }
    }
  }

  void compute3ReverseRandom( vector< size_t >& path, const vector< vector< double > >& distances )
  {
    double bestDistance = getLength( path, distances );
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
        if ( randNums[ 0 ] == randNums[ 1 ] || randNums[ 1 ] == randNums[ 2 ] || randNums[ 0 ] + 1 == randNums[ 1 ] ) {
          continue;
        }
        reverse( path.begin() + randNums[ 0 ], path.begin() + randNums[ 1 ] );
        double distance1 = getLength( path, distances );
        if ( distance1 < bestDistance ) {
          bestDistance = distance1;
          changed = true;
          continue;
        }
        reverse( path.begin() + randNums[ 1 ], path.begin() + randNums[ 2 ] + 1 );
        double distance2 = getLength( path, distances );
        if ( distance2 < bestDistance ) {
          bestDistance = distance2;
          changed = true;
          continue;
        }
        reverse( path.begin() + randNums[ 0 ], path.begin() + randNums[ 1 ] );
        double distance3 = getLength( path, distances );
        if ( distance3 < bestDistance ) {
          bestDistance = distance3;
          changed = true;
          continue;
        }
        reverse( path.begin() + randNums[ 1 ], path.begin() + randNums[ 2 ] + 1 );
      }
    }
  }

  void compute3Reverse( vector< size_t >& path, const vector< vector< double > >& distances )
  {
    double bestDistance = getLength( path, distances );
    bool changed = true;
    while ( changed ) {
      changed = false;
      for ( size_t i = 0; i < path.size(); ++i ) {
        size_t j = i + 2;
        for ( j = i + 2; j < path.size(); ++j ) {
          size_t k = j + 1;
          for ( k = j + 1; k < path.size(); ++k ) {
            reverse( path.begin() + i, path.begin() + j );
            double distance = getLength( path, distances );
            if ( distance < bestDistance ) {
              bestDistance = distance;
              changed = true;
              break;
            }
            reverse( path.begin() + j, path.begin() + k + 1 );
            distance = getLength( path, distances );
            if ( distance < bestDistance ) {
              bestDistance = distance;
              changed = true;
              break;
            }
            reverse( path.begin() + i, path.begin() + j );
            distance = getLength( path, distances );
            if ( distance < bestDistance ) {
              bestDistance = distance;
              changed = true;
              break;
            }
            reverse( path.begin() + j, path.begin() + k + 1 );
          }
          if ( k < path.size() ) {
            break;
          }
        }
      }
    }
  }

  void assertIsPath( const vector< size_t >& path, const vector< vector< double > >& distances )
  {
    assert( path.size() == distances.size() );
    vector< size_t > pathCopy( path );
    sort( pathCopy.begin(), pathCopy.end() );
    for ( size_t i = 0; i < pathCopy.size(); ++i ) {
      assert( pathCopy[ i ] == i );
    }
  }

  vector< size_t > computePath( const vector< vector< double > >& distances )
  {
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
    vector< size_t > path2( path );
    vector< size_t > path3( path );

    cerr << "1-tree distance: ";
    cerr << oneTreeLength( distances ) << endl;
    cerr << "Initial distance: " << getLength( path, distances ) << endl;
//    compute2OptPath( path, distances );
//    double twoOptDistance = getLength( path, distances );
//    cerr << "2-opt distance: " << twoOptDistance << endl;
    vector< size_t > mstPath = computeMinimumSpanningTreePath( distances );
    cerr << "mst path distance : " << getLength( mstPath, distances ) << endl;
    assertIsPath( mstPath, distances );
//    compute2OptPath( mstPath, distances );
//    cerr << "2-opt mst path distance : " << getLength( mstPath, distances ) << endl;
    compute2Reverse( path2, distances );
    cerr << "2-reverse path distance : " << getLength( path2, distances ) << endl;
    assertIsPath( path2, distances );

//    compute2OptPath( path2, distances );
//    cerr << "2-opt 2-reverse path distance : " << getLength( path2, distances ) << endl;

    compute3ReverseRandom( path3, distances );
    cerr << "3-reverse path distance : " << getLength( path3, distances ) << endl;
    assertIsPath( path3, distances );

    compute2OptPath( path3, distances );
    cerr << "2-opt 3-reverse path distance : " << getLength( path3, distances ) << endl;

    assertIsPath( path, distances );

    return path;
  }
}
