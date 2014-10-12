/* SYSTEM INCLUDES */
#include <algorithm>
#include <assert.h>
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
using namespace TravelingSalespersonProblemSolver;

namespace {
const double tolerance = 1e-9;
const bool useInfeasibleMoves = true;
const size_t maxGainMoves = 1000;
vector< size_t > better;
bool inBetterTour( size_t t1, size_t t2 )
{
  for ( size_t i = 0; i + 1 < better.size(); ++i ) {
    if ( ( better[ i ] == t1 && better[ i + 1 ] == t2 ) ||
         ( better[ i ] == t2 && better[ i + 1 ] == t1 ) ) {
      return true;
    }
  }
  return better.size() > 1 &&
         ( ( better.front() == t1 && better.back() == t2 ) ||
           ( better.front() == t2 && better.back() == t1 ) );
}

// Updates the nearest neighbor matrix on the basis of two previous tours
void updateNearest( const vector< size_t >& tour,
                    vector< size_t >& tour1,
                    vector< size_t >& tour2,
                    vector< vector< size_t > >& nearest )
{
  tour2 = tour1;
  tour1 = tour;

  set< vector< size_t > > e1;
  set< vector< size_t > > e2;
  for ( size_t i = 0; i < tour1.size(); ++i ) {
    vector< size_t > e( { tour1[ i ], tour1[ i + 1 == tour1.size() ? 0 : i + 1 ] } );
    sort( e.begin(), e.end() );
    e1.insert( e );

    e = { tour2[ i ], tour2[ i + 1 == tour2.size() ? 0 : i + 1 ] };
    sort( e.begin(), e.end() );
    e2.insert( e );
  }

  /*
  for ( const auto& v : e1 ) {
    size_t i1 = v[ 0 ];
    size_t i2 = v[ 1 ];
    if ( find( nearest[ i1 ].begin(), nearest[ i1 ].end(), i2 ) == nearest[ i1 ].end() ) {
      nearest[ i1 ].insert( nearest[ i1 ].end(), i2 );
    }
    if ( find( nearest[ i2 ].begin(), nearest[ i2 ].end(), i1 ) == nearest[ i2 ].end() ) {
      nearest[ i2 ].insert( nearest[ i2 ].end(), i1 );
    }
  }
  */

  vector< vector< size_t > > v( tour.size() );
  vector< vector< size_t > >::iterator it = set_intersection( e1.begin(), e1.end(), e2.begin(), e2.end(), v.begin() );
  v.resize( it - v.begin() );

  for ( size_t i = 0; i < v.size(); ++i ) {
    size_t i1 = v[ i ][ 0 ];
    size_t i2 = v[ i ][ 1 ];
    nearest[ i1 ].insert( nearest[ i1 ].begin(), i2 );
    size_t size1 = remove_if( nearest[ i1 ].begin() + 1, nearest[ i1 ].end(), [&] ( size_t index ) { return index == i2; } ) - nearest[ i1 ].begin();
    nearest[ i1 ].resize( size1 );

    nearest[ i2 ].insert( nearest[ i2 ].begin(), i1 );
    size_t size2 = remove_if( nearest[ i2 ].begin() + 1, nearest[ i2 ].end(), [&] ( size_t index ) { return index == i1; } ) - nearest[ i2 ].begin();
    nearest[ i2 ].resize( size2 );
  }
}

double INLINE_ATTRIBUTE getLength_( const vector< size_t >& tour,
                                    const VDistances& distances )
{
  double distance = distances( tour.back(), tour.front() );
  for ( size_t i = 0; i + 1 < tour.size(); ++i ) {
    distance += distances( tour[ i ], tour[ i + 1 ] );
  }
  return distance;
}

void printTour_( string label, const vector< size_t >& tour )
{
  cerr << label << ": ";
  for ( size_t i = 0; i < tour.size(); ++i ) {
    cerr << tour[ i ] << " ";
  }
  cerr << endl;
}

bool INLINE_ATTRIBUTE isTour_( const vector< size_t >& tour,
                               const vector< size_t >& position )
{
  assert( tour.size() == position.size() );
  for ( size_t i = 0; i < position.size(); ++i ) {
    if ( position[ i ] >= position.size() || tour[ position[ i ] ] != i ) {
      cerr << "position[ i ]: " << position[ i ] << " tour[ position[ i ] ]: " << tour[ position[ i ] ] << " i: " << i << endl;
      return false;
    }
    assert( position[ i ] < tour.size() );
  }
  return true;
}

bool INLINE_ATTRIBUTE isTour_( const vector< size_t >& tour,
                               const VDistances& distances )
{
  assert( tour.size() == distances.size() );
  vector< size_t > position( tour.size() );
  for ( size_t i = 0; i < tour.size(); ++i ) {
    position[ tour[ i ] ] = i;
  }
  return isTour_( tour, position );
}

size_t previous_( size_t node, const vector< size_t >& tour, const vector< size_t >& position )
{
  return position[ node ] > 0 ? tour[ position[ node ] - 1 ] : tour.back();
}

size_t next_( size_t node, const vector< size_t >& tour, const vector< size_t >& position )
{
  return position[ node ] + 1 < tour.size() ? tour[ position[ node ] + 1 ] : tour.front();
}

bool between_( size_t a, size_t b, size_t c, const vector< size_t >& position )
{
  return position[ a ] <= position[ c ] ? position[ a ] >= position[ b ] || position[ b ] >= position[ c ]
                                        : position[ b ] <= position[ a ] && position[ b ] >= position[ c ];
}

double getDistance_( const VDistances& distances,
                     const vector< double >& lagrangeMultipliers,
                     size_t i,
                     size_t j )
{
  return distances( i, j ) + lagrangeMultipliers[ i ] + lagrangeMultipliers[ j ];
}


vector< size_t > INLINE_ATTRIBUTE getRandomTour_( const VDistances& distances )
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

vector< size_t > INLINE_ATTRIBUTE getNearestNeighborTour_( const VDistances& distances )
{
  vector< size_t > tour;
  tour.reserve( distances.size() );
  size_t startNode = rand() % distances.size();
  tour.push_back( startNode );
  set< size_t > usedNodes;
  usedNodes.insert( startNode );
  while ( tour.size() < distances.size() ) {
    size_t currentNode = tour.back();
    size_t minUnusedIndex = (size_t)-1;
    double minUnusedDistance = numeric_limits< double >::max();
    for ( size_t i = 0; i < distances.size(); ++i ) {
      if ( usedNodes.find( i ) == usedNodes.end() ) {
        double distance = distances( currentNode, i );
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
  assert( isTour_( tour, distances ) );
  return tour;
}

vector< size_t > INLINE_ATTRIBUTE getGreedyTour_( const VDistances& distances,
                                                  const vector< double >& lagrangeMultipliers )
{
  // The greedy heuristic of matroid theory adds the shortest edge that neither
  // makes the degree of a vertex greater than 2 nor creates a cycle shorter than N.
  assert( distances.size() == lagrangeMultipliers.size() );
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
    for ( size_t j = i + 1; j < distances.size(); ++j ) {
      edges.push_back( make_pair( getDistance_( distances, lagrangeMultipliers, i, j ), make_pair( i, j ) ) );
    }
  }
  sort( edges.rbegin(), edges.rend() ); // sort from rbegin to rend => smallest element last
  while ( !edges.empty() ) {
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
    if ( !fragments[ i ].empty() ) {
      tour.swap( fragments[ i ] );
      break;
    }
  }

  assert( isTour_( tour, distances ) );
  return tour;
}

vector< size_t > INLINE_ATTRIBUTE getGreedyTour_( const VDistances& distances )
{
  vector< double > lagrangeMultipliers( distances.size(), 0.0 );
  return getGreedyTour_( distances, lagrangeMultipliers );
}

vector< size_t > INLINE_ATTRIBUTE getHelsgaunInitialTour_( const vector< vector< size_t > >& nearestNeighbors,
                                                           const vector< vector< size_t > >& helsgaunNeighbors,
                                                           const vector< vector< double > >& helsgaunDistances,
                                                           const vector< size_t >& bestTour )
{
  assert( helsgaunNeighbors.size() == helsgaunDistances.size() );
  assert( bestTour.size() == nearestNeighbors.size() );
  vector< size_t > tour( 1, rand() % nearestNeighbors.size() );
  size_t current = tour.front();
  vector< bool > added( nearestNeighbors.size(), false );
  added[ current ] = true;

  vector< size_t > bestPosition( bestTour.size() );
  for ( size_t i = 0; i < bestTour.size(); ++i ) {
    bestPosition[ bestTour[ i ] ] = i;
  }

  while ( tour.size() < bestTour.size() ) {
    current = tour.back();
    bool found = false;
    size_t bestNext = (size_t)-1;
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
/*    if ( !found ) {
      for ( size_t i = 0; i < helsgaunNeighbors[ current ].size(); ++i ) {
        size_t next = helsgaunNeighbors[ current ][ i ];
        if ( !added[ next ] &&
             helsgaunDistances[ current ][ i ] < tolerance ) {
          bestNext = next;
          added[ next ] = true;
          found = true;
          break;
        }
      }
    }
    */
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
    assert( bestNext != (size_t)-1 );
    tour.push_back( bestNext );
  }

  vector< size_t > position( tour.size() );
  for ( size_t i = 0; i < tour.size(); ++i ) {
    position[ tour[ i ] ] = i;
  }
  assert( isTour_( tour, position ) );
  return tour;
}

// Computes the optimal tour using brute force
vector< size_t > getBruteForceTour_( const VDistances& distances )
{
  assert( !distances.empty() );
  assert( distances.size() < 10 );
  // It would be possible to speed up the brute force by only considering tours
  // of one orientation, but since we use it only for small instances, this improvement
  // is unnecessary.
  vector< size_t > tour( distances.size() );
  for ( size_t i = 0; i < tour.size(); ++i ) {
    tour[ i ] = i;
  }
  vector< size_t > bestTour( tour );
  double bestDistance = getLength_( tour, distances );
  do {
    double distance = getLength_( tour, distances );
    if ( distance < bestDistance ) {
      bestDistance = distance;
      bestTour = tour;
    }
   } while ( next_permutation( tour.begin() + 1, tour.end() ) ); // tour.begin() + 1 => fixed first node
  return bestTour;
}

vector< vector< size_t > > INLINE_ATTRIBUTE computeNearestNeighbors_( const VDistances& distances,
                                                                      size_t numberOfNeighbors )
{
  numberOfNeighbors = min( numberOfNeighbors, distances.size() - 1 );
  assert( numberOfNeighbors > 0 );
  vector< vector< size_t > > nearestNeighbors( distances.size(), vector< size_t >( numberOfNeighbors ) );
  set< pair< double, size_t > > tmpNeighbors;
  for ( size_t i = 0; i < distances.size(); ++i ) {
    tmpNeighbors.clear();
    double worstNearNeighbor = numeric_limits< double >::max();
    for ( size_t j = 0; j < distances.size(); ++j ) {
      if ( j != i ) {
        if ( tmpNeighbors.size() < numberOfNeighbors || distances( i, j ) < worstNearNeighbor ) {
          if ( tmpNeighbors.size() >= numberOfNeighbors ) {
            set< pair< double, size_t > >::iterator itLast = tmpNeighbors.end();
            --itLast;
            tmpNeighbors.erase( itLast );
          }
          tmpNeighbors.insert( make_pair( distances( i, j ), j ) );
          worstNearNeighbor = tmpNeighbors.rbegin()->first;
        }
      }
    }
    assert( tmpNeighbors.size() == numberOfNeighbors );
    set< pair< double, size_t > >::const_iterator it = tmpNeighbors.begin();
    for ( size_t j = 0; j < numberOfNeighbors; ++j, ++it ) {
      nearestNeighbors[ i ][ j ] = it->second;
    }
  }
  return nearestNeighbors;
}

struct Vertex {
  Vertex( size_t idx ) : nodeIndex( idx ), cost( numeric_limits< double >::max() ) {}
  size_t nodeIndex;
  size_t parentIndexInVertexList;
  double cost;
};

double INLINE_ATTRIBUTE compute1Tree_( vector< Vertex >& vertices,
                                       vector< pair< size_t, size_t > >& edges,
                                       vector< size_t >& vertexDegrees,
                                       const VDistances& distances,
                                       const vector< double >& lagrangeMultipliers )
{
  vertices.clear();
  vertices.reserve( distances.size() );
  edges.clear();
  edges.reserve( distances.size() );
  vertexDegrees.assign( distances.size(), 0 );

  // 1. Compute length of the minimum spanning tree excluding one vertex, using Prim's algorithm.

  // Select the vertex with minimum maximal distance to a neighbor in
  // order to allow for problems with artificial nodes (such as those used
  // to construct Hamiltonian tours) with toleranceilon distance to all neighbors.
  size_t minimaxVertex = 0;
  double minimaxValue = numeric_limits< double >::max();
  for ( size_t i = 0; i + 1 < distances.size(); ++i ) {
    double value = 0.0;
    for ( size_t j = i + 1; j < distances.size(); ++j ) {
      value = max( value, distances( i, j ) );
    }
    if ( value < minimaxValue  ) {
      minimaxValue = value;
      minimaxVertex = i;
    }
  }

  vertices.push_back( Vertex( minimaxVertex ) );
  size_t rootVertex = ( minimaxVertex + 1 ) % distances.size();
  vertices.push_back( Vertex( rootVertex ) );

  vector< Vertex > unusedVertices;
  unusedVertices.reserve( distances.size() - 2 );
  for ( size_t i = 0; i < distances.size(); ++i ) {
    if ( i != minimaxVertex && i != rootVertex ) {
      unusedVertices.push_back( Vertex( i ) );
      unusedVertices.back().parentIndexInVertexList = 1;
      unusedVertices.back().cost = getDistance_( distances, lagrangeMultipliers, i, vertices[ 1 ].nodeIndex );
    }
  }

  double treeLength = 0.0;
  while ( !unusedVertices.empty() ) {
    vector< Vertex >::iterator closestUnusedVertexIt = unusedVertices.begin();
    double minDistance = numeric_limits< double >::max();
    for ( vector< Vertex >::iterator it = unusedVertices.begin(); it != unusedVertices.end(); ++it ) {
      if ( it->cost < minDistance ) {
        minDistance = it->cost;
        closestUnusedVertexIt = it;
      }
    }
    treeLength += minDistance;
    size_t indexOfNewTreeNode = closestUnusedVertexIt->nodeIndex;
    size_t indexOfClosestTreeNode = vertices[ closestUnusedVertexIt->parentIndexInVertexList ].nodeIndex;
    ++vertexDegrees[ indexOfNewTreeNode ];
    ++vertexDegrees[ indexOfClosestTreeNode ];
    vertices.push_back( *closestUnusedVertexIt );
    edges.push_back( make_pair( indexOfClosestTreeNode, indexOfNewTreeNode ) );
    *closestUnusedVertexIt = unusedVertices.back();
    unusedVertices.pop_back();

    for ( vector< Vertex >::iterator it = unusedVertices.begin(); it != unusedVertices.end(); ++it ) {
      double newDistance = getDistance_( distances, lagrangeMultipliers, indexOfNewTreeNode, it->nodeIndex );
      if ( newDistance < it->cost ) {
        it->cost = newDistance;
        it->parentIndexInVertexList = vertices.size() - 1;
      }
    }
  }

  // 2. Add the two shortest edges connecting to the first vertex
  double minElement = numeric_limits< double >::max();
  double secondMinElement = numeric_limits< double >::max();
  size_t minElementIndex = 0;
  size_t secondMinElementIndex = 0;
  for ( size_t i = 0; i < distances.size(); ++i ) {
    if ( i == minimaxVertex ) {
      continue;
    }
    double value = getDistance_( distances, lagrangeMultipliers, minimaxVertex, i );
    if ( value < secondMinElement ) {
      secondMinElement = value;
      secondMinElementIndex = i;
      if ( value < minElement ) {
        secondMinElement = minElement;
        secondMinElementIndex = minElementIndex;
        minElement = value;
        minElementIndex = i;
      }
    }
  }
  vertexDegrees[ minimaxVertex ] = 2;
  ++vertexDegrees[ secondMinElementIndex ];
  ++vertexDegrees[ minElementIndex ];
  edges.push_back( make_pair( minimaxVertex, secondMinElementIndex ) );
  edges.push_back( make_pair( minimaxVertex, minElementIndex ) );

  treeLength += minElement + secondMinElement;
  return treeLength - 2.0 * accumulate( lagrangeMultipliers.begin(), lagrangeMultipliers.end(), 0.0 );
}

double INLINE_ATTRIBUTE getHeldKarpLowerBound_( const VDistances& distances, vector< double >& lagrangeMultipliers )
{
  double maxLength = numeric_limits< double >::lowest();
  lagrangeMultipliers.assign( distances.size(), 0.0 );
  vector< double > localLagrangeMultipliers( lagrangeMultipliers );
  vector< Vertex > vertices;
  vector< size_t > vertexDegrees;
  vector< pair< size_t, size_t > > edges;
  size_t maxIter = 50;
  if ( distances.size() > 50000 ) {
    maxIter = 1;
  }
  double lambda = 0.1;
  for ( size_t i = 0; i < maxIter; ++i ) {
    vector< size_t > vertexDegrees;
    double treeLength = compute1Tree_( vertices, edges, vertexDegrees, distances, localLagrangeMultipliers );
    if ( treeLength > maxLength ) {
      maxLength = treeLength;
      lagrangeMultipliers = localLagrangeMultipliers;
    }

    double denominator = 0.0;
    for ( size_t j = 0; j < localLagrangeMultipliers.size(); ++j ) {
      double d = double( vertexDegrees[ j ] ) - 2.0;
      denominator += d * d;
    }
    if ( denominator == 0.0 ) {
      return maxLength;
    }
    double t = 2.0 * lambda * treeLength / denominator;

    for ( size_t j = 0; j < localLagrangeMultipliers.size(); ++j ) {
      localLagrangeMultipliers[ j ] += t * ( int( vertexDegrees[ j ] ) - 2 );
    }
    lambda = 1.0 / ( 20.0 + 10 * i );
  }

    return maxLength;
  }

  double INLINE_ATTRIBUTE getHeldKarpLowerBound_( const VDistances& distances )
  {
    vector< double > lagrangeMultipliers( distances.size() );
    return getHeldKarpLowerBound_( distances, lagrangeMultipliers );
  }

  vector< vector< size_t > > INLINE_ATTRIBUTE computeHelsgaunNeighbors_( const VDistances& distances,
                                                                         vector< vector< double > >& alphaDistances,
                                                                         size_t numberOfNeighbors )
  {
    numberOfNeighbors = min( numberOfNeighbors, distances.size() - 1 );
    assert( numberOfNeighbors > 0 );

    vector< Vertex > vertices;
    vector< size_t > vertexDegrees;
    vector< pair< size_t, size_t > > edges;
    vector< double > lagrangeMultipliers( distances.size(), 0.0 );
    getHeldKarpLowerBound_( distances, lagrangeMultipliers );

    compute1Tree_( vertices,
                   edges,
                   vertexDegrees,
                   distances,
                   lagrangeMultipliers );

    vector< vector< size_t > > nearestNeighbors( distances.size(), vector< size_t >( numberOfNeighbors ) );
    alphaDistances.resize( distances.size(), vector< double >( numberOfNeighbors ) );
    set< pair< double, size_t > > tmpNeighbors;
    vector< double > a( distances.size(), 0.0 );
    vector< double > b( distances.size() );
    vector< size_t > mark( distances.size(), 0 );
    size_t frontInd = vertices.front().nodeIndex;
    const double maxEdgeWeight = max( getDistance_( distances, lagrangeMultipliers, edges.end()[ -2 ].first, edges.end()[ -2 ].second ),
                                      getDistance_( distances, lagrangeMultipliers, edges.end()[ -1 ].first, edges.end()[ -1 ].second ) );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      const size_t iInd = vertices[ i ].nodeIndex;
      if ( i == 0 ) {
        for ( size_t j = 1; j < distances.size(); ++j ) {
          size_t jInd = vertices[ j ].nodeIndex;
          a[ jInd ] = getDistance_( distances, lagrangeMultipliers, jInd, iInd ) - maxEdgeWeight;
        }
      }
      else {
        b[ i ] = numeric_limits< double >::lowest();
        size_t j = 0;
        for ( size_t k = i; k != 1; k = j ) {
          j = vertices[ k ].parentIndexInVertexList;
          b[ j ] = max( b[ k ], getDistance_( distances, lagrangeMultipliers, vertices[ j ].nodeIndex, vertices[ k ].nodeIndex ) );
          mark[ j ] = i;
        }
        fill( a.begin(), a.end(), 0.0 );
        for ( j = 1; j < vertices.size(); ++j ) {
          if ( j != i ) {
            size_t jInd = vertices[ j ].nodeIndex;
            if ( mark[ j ] != i ) {
              size_t jParent = vertices[ j ].parentIndexInVertexList;
              b[ j ] = max( b[ jParent ], getDistance_( distances, lagrangeMultipliers, jInd, vertices[ jParent ].nodeIndex ) );
            }
            a[ jInd ] = getDistance_( distances, lagrangeMultipliers, iInd, jInd ) - b[ j ];
          }
        }
        a[ frontInd ] = getDistance_( distances, lagrangeMultipliers, iInd, frontInd ) - maxEdgeWeight;
      }
      tmpNeighbors.clear();
      double worstNearNeighbor = numeric_limits< double >::max();
      for ( size_t j = 0; j < distances.size(); ++j ) {
        if ( j != i ) {
          size_t jInd = vertices[ j ].nodeIndex;
          if ( tmpNeighbors.size() < numberOfNeighbors || a[ jInd ] < worstNearNeighbor ) {
            if ( tmpNeighbors.size() >= numberOfNeighbors ) {
              set< pair< double, size_t > >::iterator itLast = tmpNeighbors.end();
              --itLast;
              tmpNeighbors.erase( itLast );
            }
            tmpNeighbors.insert( make_pair( a[ jInd ], jInd ) );
            worstNearNeighbor = tmpNeighbors.rbegin()->first;
          }
        }
      }
      assert( tmpNeighbors.size() == numberOfNeighbors );
      set< pair< double, size_t > >::const_iterator it = tmpNeighbors.begin();
      for ( size_t j = 0; j < numberOfNeighbors; ++j, ++it ) {
        alphaDistances[ iInd ][ j ] = it->first;
        nearestNeighbors[ iInd ][ j ] = it->second;
      }
    }
    return nearestNeighbors;
  }

  void INLINE_ATTRIBUTE getPQI_( vector< size_t >& p,
                                 vector< size_t >& q,
                                 vector< size_t >& incl,
                                 const vector< size_t >& ts,
                                 const vector< size_t >& tour,
                                 const vector< size_t >& position )
  {
    assert( !ts.empty() );
    assert( ts.size() % 2 == 0 );
    vector< size_t > pHalf;
    pHalf.reserve( ts.size() / 2 );
    for ( size_t i = 0; i < ts.size(); i += 2 ) {
      pHalf.push_back( ts[ i ] == next_( ts[ i + 1 ], tour, position ) ? i + 1 : i );
    }
    sort( pHalf.begin(), pHalf.end(), [&] ( const size_t i, const size_t j ) {
      return position[ ts[ i ] ] < position[ ts[ j ] ] || ( position[ ts[ i ] ] == position[ ts[ j ] ] && i < j );
    } );

    p.clear();
    p.reserve( ts.size() );
    for ( size_t i = 0; i < pHalf.size(); ++i ) {
      p.push_back( pHalf[ i ] );
      p.push_back( pHalf[ i ] % 2 == 0 ? pHalf[ i ] + 1 : pHalf[ i ] - 1 );
    }

    q = vector< size_t >( p.size() );
    for ( size_t i = 0; i < p.size(); ++i ) {
      q[ p[ i ] ] = i;
    }

    incl = vector< size_t >( p.size() );
    incl.front() = ts.size() - 1;
    incl.back() = 0;
    for ( size_t i = 1; i + 1 < ts.size(); i += 2 ) {
      incl[ i ] = i + 1;
      incl[ i + 1 ] = i;
    }
  }

  void INLINE_ATTRIBUTE performKOptMove_( const vector< size_t >& ts,
                                          const vector< size_t >& incl,
                                          const vector< size_t >& p,
                                          const vector< size_t >& q,
                                          vector< size_t >& tour,
                                          vector< size_t >& position )
  {
    vector< size_t > tourCopy( tour );
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
      for ( ; currentNode != nextNode; currentNode = increasing ? next_( currentNode, tourCopy, position ) : previous_( currentNode, tourCopy, position ), ++index ) {
        tour[ index ] = currentNode;
      }
      tour[ index ] = nextNode;
      ++index;
    }
    position.assign( position.size(), (size_t)-1 );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }
  }

  void INLINE_ATTRIBUTE performKOptMove_( const vector< size_t >& ts,
                                          const vector< size_t >& incl,
                                          vector< size_t >& tour,
                                          vector< size_t >& position )
  {
    vector< size_t > p, q, inclTmp;
    getPQI_( p, q, inclTmp, ts, tour, position );
    performKOptMove_( ts, incl, p, q, tour, position );
  }

  void INLINE_ATTRIBUTE performKOptMove_( const vector< size_t >& ts,
                                          vector< size_t >& tour,
                                          vector< size_t >& position )
  {
    vector< size_t > p, q, incl;
    getPQI_( p, q, incl, ts, tour, position );
    performKOptMove_( ts, incl, p, q, tour, position );
  }

  bool INLINE_ATTRIBUTE makesTour_( const vector< size_t >& ts,
                                    const vector< size_t >& tour,
                                    const vector< size_t >& position )
  {
    vector< size_t > p, q, incl;
    getPQI_( p, q, incl, ts, tour, position );

    size_t count = 1;
    for ( size_t i = ts.size(); ( i = ( q[ incl[ p[ i - 1 ] ] ] + 1 ) ^ 1 ) != 0; ++count );
    return 2 * count == ts.size();
    /*
       vector< size_t > visited;
       visited.reserve( ts.size() );
       size_t currentIndex = 0;
       for ( size_t step = 0; step < ts.size() / 2; ++step ) {
       currentIndex = incl[ currentIndex ];
       visited.push_back( ts[ currentIndex ] );
       size_t i = q[ currentIndex ];
       currentIndex = i % 2 == 0 ? p[ i == 0 ? p.size() - 1 : i - 1 ]
       : p[ i + 1 == p.size() ? 0 : i + 1 ];
       visited.push_back( ts[ currentIndex ] );
       }
       sort( visited.begin(), visited.end() );
       vector< size_t > tsCopy( ts );
       sort( tsCopy.begin(), tsCopy.end() );
       return visited == tsCopy;
       */
  }

  bool INLINE_ATTRIBUTE twoOptInnerLoop_( vector< size_t >& tour,
                                          vector< bool >& dontLook,
                                          const VDistances& distances,
                                          const vector< vector< size_t > >& nearestNeighbors )
  {
    vector< size_t > position( tour.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }

    bool changed = false;
    vector< size_t > bestTs;
    for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
      if ( !dontLook.empty() && dontLook[ t1 ] ) {
        continue;
      }
      bool found = false;
      double maxGain = 0.0;
      for ( size_t t2choice = 0; t2choice < 2; ++t2choice  ) {
        size_t t2 = t2choice == 0 ? previous_( t1, tour, position ) : next_( t1, tour, position );
        for ( size_t t3index = 0; t3index < nearestNeighbors[ t2 ].size(); ++t3index ) {
          size_t t3 = nearestNeighbors[ t2 ][ t3index ];
          if ( t3 == previous_( t2, tour, position ) || t3 == next_( t2, tour, position ) ) {
            continue;
          }
          if ( distances( t2, t3 ) >= distances( t1, t2 ) ) {
            continue;
          }
          size_t t4 = t2choice == 0 ? next_( t3, tour, position ) : previous_( t3, tour, position );
          double gain = distances( t1, t2 ) + distances( t3, t4 ) - distances( t1, t4 ) - distances( t2, t3 );
          if ( gain > maxGain ) {
            maxGain = gain;
            bestTs = { t1, t2, t3, t4 };
          }
        }
      }
      if ( maxGain > tolerance ) {
        if ( !dontLook.empty() ) {
          for ( size_t i = 0; i < bestTs.size(); ++i ) {
            dontLook[ bestTs[ i ] ] = false;
          }
        }
        performKOptMove_( bestTs, tour, position );
        assert( isTour_( tour, position ) );
        changed = true;
        found = true;
      }
      if ( !dontLook.empty() && !found ) {
        dontLook[ t1 ] = true;
      }
    }
    return changed;
  }

  bool INLINE_ATTRIBUTE improveTour2Opt_( vector< size_t >& tour,
                                          const VDistances& distances,
                                          const vector< vector< size_t > >& nearestNeighbors )
  {
    bool change = false;
    vector< bool > dontLook;
    while ( twoOptInnerLoop_( tour, dontLook, distances, nearestNeighbors ) ) {
      change = true;
    }
    return change;
  }

  bool INLINE_ATTRIBUTE threeOptInnerLoop_( vector< size_t >& tour,
                                            vector< bool >& dontLook,
                                            const VDistances& distances,
                                            const vector< vector< size_t > >& nearestNeighbors )
  {
    bool changed = false;
    vector< size_t > position( tour.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }
    vector< size_t > bestTs;
    for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
      if ( !dontLook.empty() && dontLook[ t1 ] ) {
        continue;
      }
      bool found = false;
      double G = 0.0;
      for ( size_t t2choice = 0; t2choice < 2; ++t2choice  ) {
        size_t t2 = t2choice == 0 ? previous_( t1, tour, position ) : next_( t1, tour, position );
        for ( size_t t3index = 0; t3index < nearestNeighbors[ t2 ].size(); ++t3index ) {
          size_t t3 = nearestNeighbors[ t2 ][ t3index ];
          if ( t3 == previous_( t2, tour, position ) || t3 == next_( t2, tour, position ) ) {
            continue;
          }
          double g1 = distances( t1, t2 ) - distances( t2, t3 );
          if ( g1 <= tolerance ) {
            continue;
          }
          // First choice of t4
          size_t t4 = t2choice == 0 ? next_( t3, tour, position ) : previous_( t3, tour, position );
          if ( t4 == previous_( t2, tour, position ) || t4 == next_( t2, tour, position ) ) {
            continue;
          }
          {
            // Test for improving 2-opt move
            double gain = g1 + distances( t3, t4 ) - distances( t4, t1 );
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
            double g2 = distances( t3, t4 ) - distances( t4, t5 );
            if ( g1 + g2 <= tolerance ) {
              continue;
            }

            // Select t6 such that a valid tour is created
            size_t t6 = between_( t2, t4, t5, position ) ? next_( t5, tour, position ) : previous_( t5, tour, position );
            if ( t6 == t1 ) {
              continue;
            }
            double g3 = distances( t5, t6 ) - distances( t6, t1 );
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
            double g2 = distances( t3, t4 ) - distances( t4, t5 );
            if ( g1 + g2 <= tolerance ) {
              ;
            }
            // Only consider one choice of t6. The other choice is possible, but clutters the code and doesn't lead to a significant improvement.
            size_t t6choice = t2choice;
            size_t t6 = t6choice == 0 ? next_( t5, tour, position ) : previous_( t5, tour, position );
            if ( t6 == t3 || t6 == t2 || t6 == t1 ) {
              continue;
            }
            double g3 = distances( t5, t6 ) - distances( t6, t1 );
            double gain = g1 + g2 + g3;
            if ( gain > G ) {
              G = gain;
              bestTs = { t3, t4, t5, t6, t1, t2 };
            }
          }
        }
      }
      if ( G > tolerance ) {
        if ( !dontLook.empty() ) {
          for ( size_t i = 0; i < bestTs.size(); ++i ) {
            dontLook[ bestTs[ i ] ] = false;
          }
        }
        performKOptMove_( bestTs, tour, position );
        assert( isTour_( tour, position ) );
        found = true;
        changed = true;
      }
      if ( !dontLook.empty() && !found ) {
        dontLook[ t1 ] = true;
      }
    }
    return changed;
  }

  bool INLINE_ATTRIBUTE improveTour3Opt_( vector< size_t >& tour,
                                          const VDistances& distances,
                                          const vector< vector< size_t > >& nearestNeighbors )
  {
    bool change = false;
    vector< bool > dontLook;
    while ( threeOptInnerLoop_( tour, dontLook, distances, nearestNeighbors ) ) {
      change = true;
    }
    return change;
  }

  bool INLINE_ATTRIBUTE improveTour23InnerLoop_( vector< size_t >& tour,
                                                 vector< bool >& dontLook,
                                                 vector< bool >& otherDontLook,
                                                 const VDistances& distances,
                                                 const vector< vector< size_t > >& nearestNeighbors )
  {
    // Performs infeasible 2-opt moves followed by a 2- or 3-opt move.
    // Includes double bridge moves.
    vector< size_t > position( tour.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }
    bool changed = false;
    vector< size_t > bestTs;
    for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
      if ( !dontLook.empty() && dontLook[ t1 ] ) {
        continue;
      }
      bool found = false;
      double maxGain = 0.0;
      vector< size_t > incl;
      for ( size_t t2choice = 0; t2choice < 1; ++t2choice ) {
        size_t t2 = t2choice == 0 ? next_( t1, tour, position ) : previous_( t1, tour, position );

        for ( size_t t3index = 0; t3index < nearestNeighbors[ t2 ].size(); ++t3index ) {
          size_t t3 = nearestNeighbors[ t2 ][ t3index ];
          size_t t4 = t2choice == 0 ? next_( t3, tour, position ) : previous_( t3, tour, position );
          if ( t3 == t1 || t4 == t1 ) {
            continue;
          }
          double gainFirstBridge = distances( t1, t2 ) + distances( t3, t4 ) - distances( t2, t3 ) - distances( t1, t4 );
          if ( gainFirstBridge <= tolerance ) {
            continue;
          }

          for ( size_t t5 = t2; t5 != t3; t5 = t2choice == 0 ? next_( t5, tour, position ) : previous_( t5, tour, position ) ) {
            for ( size_t t6choice = 0; t6choice < 1; ++t6choice ) {
              size_t t6 = t6choice == 0 ? next_( t5, tour, position ) : previous_( t5, tour, position );
              if ( t5 == t2 && t6 == t1 ) {
                continue;
              }
              for ( size_t t7index = 0; t7index < nearestNeighbors[ t6 ].size(); ++t7index ) {
                size_t t7 = nearestNeighbors[ t6 ][ t7index ];
                for ( size_t t8choice = 0; t8choice < 2; ++t8choice ) {
                  size_t t8 = t8choice == 0 ? next_( t7, tour, position ) : previous_( t7, tour, position );
                  if ( !between_( t4, t1, t7, position ) || !between_( t4, t1, t8, position ) ) {
                    continue;
                  }
                  double gain2 =
                    t8choice == t6choice ? distances( t5, t6 ) + distances( t7, t8 ) - distances( t6, t7 ) - distances( t5, t8 )
                                         : distances( t5, t6 ) + distances( t7, t8 ) - distances( t5, t7 ) - distances( t6, t8 );
                  if ( gainFirstBridge + gain2 > maxGain ) {
                    maxGain = gainFirstBridge + gain2;
                    bestTs = { t1, t2, t3, t4, t5, t6, t7, t8 };
                    if ( t8choice == t6choice ) {
                      incl = { 3, 2, 1, 0, 7, 6, 5, 4 };
                    }
                    else {
                      incl = { 3, 2, 1, 0, 6, 7, 4, 5 };
                    }
                  }

                  for ( size_t t9index = 0; t9index < nearestNeighbors[ t8 ].size(); ++t9index ) {
                    size_t t9 = nearestNeighbors[ t8 ][ t9index ];
                    size_t t10 = t8choice == 0 || between_( t2, t3, t9, position ) ? previous_( t9, tour, position ) : next_( t9, tour, position );
                    if ( //t9 == t5 || t9 == t6 || t9 == next_( t5, tour, position ) || t9 == previous_( t5, tour, position ) ||
                         ( t9 == t1 && t10 == t2 ) ||
                         ( t9 == t2 && t10 == t1 ) ||
                         ( t9 == t3 && t10 == t4 ) ||
                         ( t9 == t4 && t10 == t3 ) ||
                         ( t9 == t5 && t10 == t6 ) ||
                         ( t9 == t6 && t10 == t5 ) ||
                         ( t9 == t7 && t10 == t8 ) ||
                         ( t9 == t8 && t10 == t7 ) ) {
                      continue;
                    }

                    double gain3 =
                      distances( t5, t6 ) + distances( t7, t8 ) + distances( t9, t10 ) - distances( t6, t7 ) - distances( t8, t9 ) - distances( t10, t5 );
                    if ( gainFirstBridge + gain3 > maxGain ) {
                      maxGain = gainFirstBridge + gain3;
                      bestTs = { t1, t2, t3, t4, t5, t6, t7, t8, t9, t10 };
                      incl = { 3, 2, 1, 0, 9, 6, 5, 8, 7, 4 };
                    }
                  }
                }
              }
            }
          }
        }
      }
      if ( maxGain > tolerance ) {
        if ( !dontLook.empty() ) {
          for ( size_t i = 0; i < bestTs.size(); ++i ) {
            dontLook[ bestTs[ i ] ] = false;
          }
        }
        if ( !otherDontLook.empty() ) {
          for ( size_t i = 0; i < bestTs.size(); ++i ) {
            otherDontLook[ bestTs[ i ] ] = false;
          }
        }
        performKOptMove_( bestTs, incl, tour, position );
        assert( isTour_( tour, position ) );
        changed = true;
        found = true;
      }
      if ( !found && !dontLook.empty() ) {
        dontLook[ t1 ] = true;
      }
    }
    return changed;
  }

  bool INLINE_ATTRIBUTE improveTour23_( vector< size_t >& tour,
                                        const VDistances& distances,
                                        const vector< vector< size_t > >& nearestNeighbors )
  {
    if ( tour.size() < 5 ) {
      return false;
    }
    bool change = false;
    vector< bool > dontLook( tour.size(), false );
    vector< bool > otherDontLook;
    while ( improveTour23InnerLoop_( tour, dontLook, otherDontLook, distances, nearestNeighbors ) ) {
      change = true;
    }
    return change;
  }

  double INLINE_ATTRIBUTE bestKOptMove_( size_t depth,
                                         size_t k,
                                         vector< size_t >& ts,
                                         double G0,
                                         vector< pair< size_t, size_t > >& added,
                                         vector< pair< size_t, size_t > >& removed,
                                         vector< size_t >& bestTs,
                                         double& bestG,
                                         vector< size_t >& bestInfeasibleTs,
                                         double& bestInfeasibleG,
                                         const vector< size_t >& tour,
                                         const vector< size_t >& position,
                                         const VDistances& distances,
                                         const vector< vector< size_t > >& nearestNeighbors )
  {
    assert( ts.size() > 1 );
    const size_t tb = ts.back();
    vector< pair< size_t, size_t > > tcTdPairs;
    for ( size_t tcIndex = 0; tcIndex < nearestNeighbors[ tb ].size(); ++tcIndex ) {
      size_t tc = nearestNeighbors[ tb ][ tcIndex ];
      double G1 = G0 - distances( tb, tc );
      if ( G1 <= tolerance ) {
        continue;
      }
      if ( tc == next_( tb, tour, position ) ||
           tc == previous_( tb, tour, position ) ) {
        // The added edge should not belong to T
        continue;
      }
      for ( size_t tdChoice = 0; tdChoice < 2; ++tdChoice ) {
        size_t td = tdChoice == 0 ? previous_( tc, tour, position ) : next_( tc, tour, position );
        if ( depth + 1 == k ) {
          double G2 = G1 + distances( tc, td );
          if ( G2 - distances( ts.front(), td ) <= tolerance && ( G2 <= bestG && G2 <= bestInfeasibleG ) ) {
            continue;
          }
          if ( find( added.begin(), added.end(), make_pair( tc, td ) ) != added.end() ) {
            continue;
          }
        }
        if ( find( removed.begin(), removed.end(), make_pair( tc, td ) ) != removed.end() ) {
          continue;
        }
        tcTdPairs.push_back( make_pair( tc, td ) );
      }
    }
/*    sort( tcTdPairs.begin(), tcTdPairs.end(), [&] ( const pair< size_t, size_t >& p1,
                                                    const pair< size_t, size_t >& p2 ) {
      const size_t tc1 = p1.first;
      const size_t td1 = p1.second;
      const size_t tc2 = p2.first;
      const size_t td2 = p2.second;
      return distances( tc1, td1 ) - distances( tb, tc1 )
             < distances( tc2, td2 ) - distances( tb, tc2 );
    } );
    */
    reverse( tcTdPairs.begin(), tcTdPairs.end() );

    while ( !tcTdPairs.empty() ) {
      size_t tc = tcTdPairs.back().first;
      size_t td = tcTdPairs.back().second;
      tcTdPairs.pop_back();

      ts.push_back( tc );
      ts.push_back( td );
      double G1 = G0 - distances( tb, tc );
      double G2 = G1 + distances( tc, td );
      double gain = G2 - distances( ts.front(), ts.back() );
      if ( gain > tolerance && makesTour_( ts, tour, position ) ) {
        return gain;
      }

      if ( depth + 1 < k ) {
        added.push_back( make_pair( tb, tc ) );
        added.push_back( make_pair( tc, tb ) );
        removed.push_back( make_pair( tc, td ) );
        removed.push_back( make_pair( td, tc ) );
        gain = bestKOptMove_( depth + 1,
                              k,
                              ts,
                              G2,
                              added,
                              removed,
                              bestTs,
                              bestG,
                              bestInfeasibleTs,
                              bestInfeasibleG,
                              tour,
                              position,
                              distances,
                              nearestNeighbors );
        if ( gain > tolerance ) {
          return gain;
        }
        removed.resize( removed.size() - 2 );
        added.resize( added.size() - 2 );
      }
      else if ( G2 > bestG && makesTour_( ts, tour, position ) ) {
        bestG = G2;
        bestTs = ts;
      }
      else if ( G2 > bestInfeasibleG ) {
        bestInfeasibleG = G2;
        bestInfeasibleTs = ts;
      }

      ts.resize( ts.size() - 2 );
    }
    return -1.0;
  }

  bool INLINE_ATTRIBUTE kOptOuterLoop_( size_t k,
                                        vector< size_t >& tour,
                                        vector< bool >& dontLook,
                                        const VDistances& distances,
                                        const vector< vector< size_t > >& nearestNeighbors,
                                        bool linKernighan )
  {
    vector< size_t > position( tour.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }
    bool anyChange = false;
    for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
      if ( dontLook[ t1 ] ) {
        continue;
      }
      bool found = false;
      for ( size_t t2choice = 0; t2choice < 2; ++t2choice ) {
        size_t t2 = t2choice == 0 ? previous_( t1, tour, position ) : next_( t1, tour, position );
        if ( inBetterTour( t1, t2 ) ) {
          continue;
        }

        vector< size_t > ts( { t1, t2 } );
        double G0 = distances( t1, t2 );
        vector< pair< size_t, size_t > > added;
        vector< pair< size_t, size_t > > removed( { make_pair( t1, t2 ), make_pair( t2, t1 ) } );
        vector< size_t > bestTs;
        double bestG = tolerance;
        vector< size_t > bestInfeasibleTs;
        double bestInfeasibleG = tolerance;
        bool testChange = false;
        vector< size_t > tourCopy( tour );
        vector< size_t > positionCopy( position );
        vector< size_t > tsHistory;
        size_t lkDepth = 0;
        do {
          ++lkDepth;
          bestG = tolerance;
          bestInfeasibleG = useInfeasibleMoves ? tolerance : numeric_limits< double >::max();
          double gain = bestKOptMove_( 1,
                                       k,
                                       ts,
                                       G0,
                                       added,
                                       removed,
                                       bestTs,
                                       bestG,
                                       bestInfeasibleTs,
                                       bestInfeasibleG,
                                       tour,
                                       position,
                                       distances,
                                       nearestNeighbors );
          if ( gain > tolerance ) {
            performKOptMove_( ts, tour, position );
            assert( isTour_( tour, position ) );
            tsHistory.insert( tsHistory.end(), ts.begin(), ts.end() );
            for ( size_t i = 0; i < tsHistory.size(); ++i ) {
              dontLook[ tsHistory[ i ] ] = false;
            }
            found = true;
            anyChange = true;
            testChange = false;
            break;
          }
          if ( !linKernighan ) {
            break;
          }
          if ( !useInfeasibleMoves ) {
            bestInfeasibleG = numeric_limits< double >::lowest();
          }
          if ( bestG > tolerance ) {
            performKOptMove_( bestTs, tour, position );
            testChange = true;
            G0 = bestG;
            ts = { bestTs.front(), bestTs.back() };
            tsHistory.insert( tsHistory.end(), bestTs.begin(), bestTs.end() );
            removed.clear();
            added.clear();
            removed.push_back( make_pair( bestTs.front(), bestTs.back() ) );
            removed.push_back( make_pair( bestTs.back(), bestTs.front() ) );
            for ( size_t i = 1; i + 1 < tsHistory.size(); i += 2 ) {
              added.push_back( make_pair( tsHistory[ i ], tsHistory[ i + 1 ] ) );
              added.push_back( make_pair( tsHistory[ i + 1 ], tsHistory[ i ] ) );
            }
          }
          else if ( bestInfeasibleG > tolerance ) {
            G0 = bestInfeasibleG;
            ts = bestInfeasibleTs;
            removed.clear();
            added.clear();
            for ( size_t i = 0; i + 1 < ts.size(); i += 2 ) {
              removed.push_back( make_pair( ts[ i ], ts[ i + 1 ] ) );
              removed.push_back( make_pair( ts[ i + 1 ], ts[ i ] ) );
            }
            for ( size_t i = 1; i + 1 < ts.size(); i += 2 ) {
              added.push_back( make_pair( ts[ i ], ts[ i + 1 ] ) );
              added.push_back( make_pair( ts[ i + 1 ], ts[ i ] ) );
            }
          }
        } while ( ( bestG > tolerance && lkDepth < maxGainMoves ) || ( bestInfeasibleG > tolerance && ts.size() < 2 * k + 1 ) );
        if ( testChange ) {
          tour = tourCopy;
          position = positionCopy;
        }
      }
      if ( !found ) {
        dontLook[ t1 ] = true;
      }
    }
    return anyChange;
  }


  bool INLINE_ATTRIBUTE improveTourKOpt_( size_t k,
                                          bool linKernighan,
                                          vector< size_t >& tour,
                                          const VDistances& distances,
                                          const vector< vector< size_t > >& nearestNeighbors )
  {
    vector< bool > dontLook( tour.size(), false );
    bool change = false;
    bool change4 = true;
    while ( change4 ) {
      change4 = false;
      while ( kOptOuterLoop_( k, tour, dontLook, distances, nearestNeighbors, linKernighan ) ) {
        change = true;
      }

      vector< bool > dontLook4( tour.size(), false );
      while ( improveTour23InnerLoop_( tour, dontLook4, dontLook, distances, nearestNeighbors ) ) {
        change4 = true;
        change = true;
      }
    }
    return change;
  }

  void kSwapKick( size_t k,
                  vector< size_t >& tour,
                  vector< bool >& dontLook )
  {
    vector< size_t > position( tour.size() );
    for ( size_t i = 0; i < tour.size(); ++i ) {
      position[ tour[ i ] ] = i;
    }

    vector< size_t > swapEdges;
    swapEdges.reserve( k );
    while ( swapEdges.size() < k && swapEdges.size() < tour.size() / 2 ) {
      size_t i = rand() % tour.size();
      if ( find( swapEdges.begin(), swapEdges.end(), i ) == swapEdges.end() ) {
        swapEdges.push_back( i );
      }
    }
    sort( swapEdges.begin(), swapEdges.end(), [&] ( const size_t ta, const size_t tb ) { return position[ ta ] < position[ tb ]; } );
    vector< size_t > ts;
    ts.reserve( 2 * swapEdges.size() );
    for ( size_t i = 0; i < swapEdges.size(); ++i ) {
      ts.push_back( swapEdges[ i ] );
      ts.push_back( next_( ts.back(), tour, position ) );
    }
    vector< size_t > tourCopy( tour );
    size_t index = 0;
    for ( size_t i = 0; i + 1 < ts.size(); i += 2 ) {
      for ( size_t t = ts[ i ]; t != ( i == 0 ? ts.back() : ts[ i - 1 ] ); t = previous_( t, tourCopy, position ), ++index ) {
        tour[ index ] = t;
      }
      tour[ index ] = i == 0 ? ts.back() : ts[ i - 1 ];
      ++index;
    }
    assert( index == tour.size() );

    for ( size_t i = 0; i < ts.size(); ++i ) {
      dontLook[ ts[ i ] ] = false;
    }
  }

  bool INLINE_ATTRIBUTE improveTourIterated_( function< bool ( vector< size_t >&,
                                                               vector< bool >&,
                                                               const VDistances&,
                                                               const vector< vector< size_t > >& ) > improveTour,
                                              size_t iterations,
                                              bool useGain23,
                                              vector< size_t >& tour,
                                              const VDistances& distances,
                                              const vector< vector< size_t > >& nearestNeighbors )
  {
    assert( tour.size() > 8 );
    bool change = false;
    vector< bool > dontLook( tour.size(), false );
    vector< size_t > bestTour( tour );
    double bestLength = getLength_( tour, distances );
    vector< size_t > tour1( tour ), tour2( tour );
    vector< vector< size_t > > nn( nearestNeighbors );
    for ( size_t i = 0; i < iterations; ++i ) {
      bool change4 = true;
      while ( change4 ) {
        change4 = false;
        while ( improveTour( tour, dontLook, distances, nn ) ) {
          change = true;
        }

        if ( useGain23 ) {
          vector< bool > dontLook4( tour.size(), false );
          while ( improveTour23InnerLoop_( tour, dontLook4, dontLook, distances, nn ) ) {
            change4 = true;
            change = true;
          }
        }
      }
      double length = getLength_( tour, distances );
      if ( length < bestLength ) {
        bestTour = tour;
        bestLength = length;
        better = bestTour;
        updateNearest( tour, tour1, tour2, nn );
      }
      if ( i + 1 == iterations ) {
        break;
      }
      tour = bestTour;
      kSwapKick( 6, tour, dontLook );
    }
    tour = bestTour;
    return change;
  }

  bool INLINE_ATTRIBUTE improveTourIterated2Opt_( size_t iterations,
                                                  bool useGain23,
                                                  vector< size_t >& tour,
                                                  const VDistances& distances,
                                                  const vector< vector< size_t > >& nearestNeighbors )
  {
    return improveTourIterated_( &twoOptInnerLoop_, iterations, useGain23, tour, distances, nearestNeighbors );
  }

  bool INLINE_ATTRIBUTE improveTourIterated3Opt_( size_t iterations,
                                                  bool useGain23,
                                                  vector< size_t >& tour,
                                                  const VDistances& distances,
                                                  const vector< vector< size_t > >& nearestNeighbors )
  {
    return improveTourIterated_( &threeOptInnerLoop_, iterations, useGain23, tour, distances, nearestNeighbors );
  }

  bool INLINE_ATTRIBUTE improveTourIteratedLKH_( size_t k,
                                                 size_t iterations,
                                                 bool useGain23,
                                                 bool linKernighan,
                                                 vector< size_t >& tour,
                                                 const VDistances& distances,
                                                 const vector< vector< size_t > >& nearestNeighbors )
  {
    auto improveTour = [ k, linKernighan ] ( vector< size_t >& tour,
                                             vector< bool >& dontLook,
                                             const VDistances& distances,
                                             const vector< vector< size_t > >& nearestNeighbors ) {
      return kOptOuterLoop_( k, tour, dontLook, distances, nearestNeighbors, linKernighan );
    };
    return improveTourIterated_( improveTour, iterations, useGain23, tour, distances, nearestNeighbors );
  }

  double inline distanceRound( double d ) { return double( long( d + 0.5 ) ); }
} // anonymous namespace

// VDistances
VDistances::VDistances( const vector< vector< double > >& points,
                        function< double ( double ) > rounding ) :
  rounding_( rounding )
{
  for ( size_t i = 1; i < points.size(); ++i ) {
    assert( points[ i ].size() == points[ 0 ].size() );
  }
}
bool VDistances::empty() const { return size() == 0; }
double VDistances::computeDistance_( const vector< double >& point1, const vector< double >& point2 ) const
{
  double dist = 0.0;
  for ( size_t i = 0; i < point1.size(); ++i ) {
    double diff = point1[ i ] - point2[ i ];
    dist += diff * diff;
  }
  return rounding_( sqrt( dist ) );
}

// MatrixDistances
MatrixDistances::MatrixDistances( const vector< vector< double > >& points,
                                  function< double ( double ) > rounding ) :
  VDistances( points, rounding )
{
  distances_.assign( points.size(), vector< double >( points.size(), 0.0 ) );
  for ( size_t i = 0; i < points.size(); ++i ) {
    for ( size_t j = i + 1; j < points.size(); ++j ) {
      distances_[ i ][ j ] = distances_[ j ][ i ] = computeDistance_( points[ i ], points[ j ] );
    }
  }
}
double MatrixDistances::operator()( size_t i, size_t j ) const { return distances_[ i ][ j ]; }
size_t MatrixDistances::size() const { return distances_.size(); }
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
  distances_ = distances;
}

// MatrixRoundedDistances
MatrixRoundedDistances::MatrixRoundedDistances( const vector< vector< double > >& points ) :
  MatrixDistances( points, [] ( double d ) { return distanceRound( d ); } )
{}

// OnTheFlyDistances
OnTheFlyDistances::OnTheFlyDistances( const vector< vector< double > >& points,
                                  function< double ( double ) > rounding ) :
  VDistances( points, rounding ),
  points_( points )
{
  for ( size_t i = 0; i < points_.size(); ++i ) {
    assert( points_[ i ].size() == points_[ 0 ].size() );
  }
}
double OnTheFlyDistances::operator()( size_t i, size_t j ) const { return computeDistance_( points_[ i ], points_[ j ] ); }
size_t OnTheFlyDistances::size() const { return points_.size(); }

// OnTheFlyRoundedDistances
OnTheFlyRoundedDistances::OnTheFlyRoundedDistances( const vector< vector< double > >& points ) :
  OnTheFlyDistances( points, [] ( double d ) { return distanceRound( d ); } )
{}

vector< size_t > INLINE_ATTRIBUTE TravelingSalespersonProblemSolver::computeTour( const VDistances& distances )
{
  {
    double start( clock() );
    cerr << "Starting sanity check ... ";
    assert( !distances.empty() );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      for ( size_t j = 0; j < distances.size(); ++j ) {
        assert( fabs( distances( i, j ) - distances( j, i ) ) < 1e-9 );
        if ( i != j ) {
          assert( distances( i, j ) > 0.0 );
        }
        else {
          if ( distances( i, i ) != 0.0 ) {
            cerr << "i " << i << " " << distances( i, i ) << endl;
          }
          assert( distances( i, i ) == 0.0 );
        }
      }
    }
    double timeSanity( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "DONE: " << timeSanity << endl;
  }
  double start( clock() );
  const vector< size_t > tourRand = getRandomTour_( distances );
  double timeRand( ( clock() - start ) / CLOCKS_PER_SEC );
  start = clock();
  const vector< size_t > tourNN = getNearestNeighborTour_( distances );
  double timeNN( ( clock() - start ) / CLOCKS_PER_SEC );
  cerr << setprecision( 8 );

  vector< vector< size_t > > nearestNeighbors5( distances.size() );
  vector< vector< size_t > > nearestNeighbors10( distances.size() );
  vector< vector< size_t > > nearestNeighbors20( distances.size() );
  vector< vector< size_t > > nearestNeighbors30( distances.size() );
  vector< vector< size_t > > helsgaun10( distances.size() );
  vector< vector< size_t > > helsgaun5( distances.size() );
  vector< vector< double > > helsgaunDistances10;
  {
    double start( clock() );
    nearestNeighbors30 = computeNearestNeighbors_( distances, 30 );
    double timeNN( ( clock() - start ) / CLOCKS_PER_SEC );
    start = clock();
    helsgaun10 = computeHelsgaunNeighbors_( distances, helsgaunDistances10, 10 );
    double timeHelsgaun( ( clock() - start ) / CLOCKS_PER_SEC );
    start = clock();
    for ( size_t i = 0; i < distances.size(); ++i ) {
      nearestNeighbors5[ i ] = vector< size_t >( nearestNeighbors30[ i ].begin(), nearestNeighbors30[ i ].begin() + min( size_t( 5 ), nearestNeighbors30[ i ].size() ) );
      nearestNeighbors10[ i ] = vector< size_t >( nearestNeighbors30[ i ].begin(), nearestNeighbors30[ i ].begin() + min( size_t( 10 ), nearestNeighbors30[ i ].size() ) );
      nearestNeighbors20[ i ] = vector< size_t >( nearestNeighbors30[ i ].begin(), nearestNeighbors30[ i ].begin() + min( size_t( 20 ), nearestNeighbors30[ i ].size() ) );
      helsgaun5[ i ] = vector< size_t >( helsgaun10[ i ].begin(), helsgaun10[ i ].begin() + min( size_t( 5 ), helsgaun10[ i ].size() ) );
    }
    auto addVector = [&] ( vector< vector< size_t > >& nn, const vector< vector< size_t > >& hn ) {
      assert( nn.size() == hn.size() );
      for ( size_t i = 0; i < nn.size(); ++i ) {
        nn[ i ].insert( nn[ i ].begin(), hn[ i ].begin(), hn[ i ].end() );
        vector< size_t >::iterator it = remove_if( nn[ i ].begin() + hn[ i ].size(), nn[ i ].end(),
                                                   [&] ( size_t index ) { return find( hn[ i ].begin(), hn[ i ].end(), index ) != hn[ i ].end(); } );
        nn[ i ].resize( it - nn[ i ].begin() );
      }
    };
    addVector( nearestNeighbors30, helsgaun10 );
    addVector( nearestNeighbors20, helsgaun5 );
    addVector( nearestNeighbors10, helsgaun5 );
    addVector( nearestNeighbors5, helsgaun5 );

    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "Time to compute " << nearestNeighbors30.front().size() << " nearest neighbors: " << timeNN << endl;
    cerr << "Time to compute " << helsgaun10.front().size() << " Helsgaun neighbors: " << timeHelsgaun << endl;
    cerr << "Time to compute rest " << time << endl;
  }

  start = clock();
  const vector< size_t > tourGreedy = getGreedyTour_( distances );
  double timeGreedy( ( clock() - start ) / CLOCKS_PER_SEC );
  vector< size_t > tour( tourGreedy );

  cerr << "TimeRand: " << timeRand << ", timeNN " << timeNN << endl;
  cerr << "Greedy distance: " << getLength_( tourGreedy, distances ) << ", time: " << timeGreedy << endl;
  if ( distances.size() < 10 ) {
    tour = getBruteForceTour_( distances );
  }

  if ( true ) {
    cerr << "1-tree distance: ";
    double start( clock() );
    double lowerBound = getHeldKarpLowerBound_( distances );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << lowerBound << ", time: " << time << endl;
  }

  if ( true ) {
    double start( clock() );
    improveTour2Opt_( tour, distances, nearestNeighbors30 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "2-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
    assert( isTour_( tour, distances ) );
  }

  if ( true ) {
    tour = tourGreedy;
    double start( clock() );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "3-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
    assert( isTour_( tour, distances ) );
  }

  if ( true ) {
    tour = tourGreedy;
    double start( clock() );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    improveTour23_( tour, distances, nearestNeighbors20 );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    improveTour23_( tour, distances, nearestNeighbors20 );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "4-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  if ( true ) {
    tour = tourGreedy;
    double start( clock() );
    improveTourKOpt_( 5, false, tour, distances, nearestNeighbors20 );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    improveTour23_( tour, distances, nearestNeighbors20 );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "5-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  if ( false ) {
    tour = tourGreedy;
    double start( clock() );
    improveTourIterated2Opt_( 1000, false, tour, distances, nearestNeighbors5 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "I2    tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  if ( true ) {
    tour = tourGreedy;
    double start( clock() );
    improveTourIterated3Opt_( 1000, false, tour, distances, nearestNeighbors30 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "I3    tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  for ( size_t k = 5; k < 6; ++k ) {
    if ( true ) {
      tour = tourGreedy;
      double start( clock() );
      vector< size_t > tour1( tour );
      vector< size_t > tour2( tour );
      vector< vector< size_t > > nn( nearestNeighbors5 );
      vector< size_t > bestTour( tour );
      vector< vector< double > > points;
      for ( size_t i = 0; i < 10; ++i ) {
        tour = getHelsgaunInitialTour_( nn, helsgaun10, helsgaunDistances10, bestTour );
        bool useGain23 = true;
        while ( improveTourKOpt_( k, useGain23, tour, distances, nn ) ) {
          updateNearest( tour, tour1, tour2, nn );
        }
        if ( getLength_( tour, distances ) < getLength_( bestTour, distances ) ) {
          bestTour = tour;
          better = bestTour;
        }
      }
      tour = bestTour;

      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << k << "-LK  tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
    }
  }

  for ( size_t k = 5; k < 6; ++k ) {
    if ( true ) {
      tour = tourGreedy;
      double start( clock() );
      bool useGain23 = false;
      bool linKernighan = true;
      improveTourIteratedLKH_( k, 100, useGain23, linKernighan, tour, distances, nearestNeighbors5 );

      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << k << "-ILK tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
    }
  }

  if ( false ) {
    printTour_( "a", tourGreedy );
  }

  return tour;
}
