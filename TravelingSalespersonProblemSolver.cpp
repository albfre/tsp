/* SYSTEM INCLUDES */
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <limits>
#include <set>

/* HEADER */
#include "TravelingSalespersonProblemSolver.h"

// For profiling
#if 1
#define INLINE_ATTRIBUTE __attribute__ ((noinline))
#else
#define INLINE_ATTRIBUTE
#endif

using namespace std;

namespace {
const double tolerance = 1e-9;

vector< size_t > tour1, tour2;
void updateNearest( const vector< size_t >& tour, vector< vector< size_t > >& nearest )
{
  tour2 = tour1;
  tour1 = tour;

  set< vector< size_t > > e1;
  set< vector< size_t > > e2;
  for ( size_t i = 0; i + 1 < tour1.size(); ++i ) {
    vector< size_t > e( { tour1[ i ], tour1[ i + 1 ] } );
    sort( e.begin(), e.end() );
    e1.insert( e );

    e = { tour2[ i ], tour2[ i + 1 ] };
    sort( e.begin(), e.end() );
    e2.insert( e );
  }
  vector< size_t > e( { tour1.front(), tour1.back() } );
  sort( e.begin(), e.end() );
  e1.insert( e );

  e = { tour2.front(), tour2.back() };
  sort( e.begin(), e.end() );
  e2.insert( e );


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

  vector< vector< size_t > > v( tour.size() );
  vector< vector< size_t > >::iterator it = set_intersection( e1.begin(), e1.end(), e2.begin(), e2.end(), v.begin() );
  v.resize( it - v.begin() );

  for ( size_t i = 0; i < v.size(); ++i ) {
    size_t i1 = v[ i ][ 0 ];
    size_t i2 = v[ i ][ 1 ];
    nearest[ i1 ].insert( nearest[ i1 ].begin(), i2 );
    if ( find( nearest[ i1 ].begin(), nearest[ i1 ].end(), i2 ) != nearest[ i1 ].end() ) {
      size_t ii = 1;
      while ( ii < nearest[ i1 ].size() && nearest[ i1 ][ ii ] != i2 ) {
        ++ii;
      }
      while ( ii + 1 < nearest[ i1 ].size() ) {
        nearest[ i1 ][ ii ] = nearest[ i1 ][ ii + 1 ];
        ++ii;
      }
      nearest[ i1 ].resize( nearest[ i1 ].size() - 1 );
    }
    nearest[ i2 ].insert( nearest[ i2 ].begin(), i1 );
    if ( find( nearest[ i2 ].begin(), nearest[ i2 ].end(), i1 ) != nearest[ i2 ].end() ) {
      size_t ii = 1;
      while ( ii < nearest[ i2 ].size() && nearest[ i2 ][ ii ] != i1 ) {
        ++ii;
      }
      while ( ii + 1 < nearest[ i2 ].size() ) {
        nearest[ i2 ][ ii ] = nearest[ i2 ][ ii + 1 ];
        ++ii;
      }
      nearest[ i2 ].resize( nearest[ i2 ].size() - 1 );
    }
  }
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
                               const vector< vector< double > >& distances )
{
  assert( tour.size() == distances.size() );
  vector< size_t > position( tour.size() );
  for ( size_t i = 0; i < tour.size(); ++i ) {
    position[ tour[ i ] ] = i;
  }
  return isTour_( tour, position );
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

vector< size_t > INLINE_ATTRIBUTE getNearestNeighborTour_( const vector< vector< double > >& distances )
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
  assert( isTour_( tour, distances ) );
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

vector< size_t > INLINE_ATTRIBUTE getHelsgaunInitialTour_( const vector< vector< size_t > >& nearestNeighbors,
                                                           const vector< vector< double > >& distances,
                                                           const vector< size_t >& bestTour )
{
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
    for ( size_t i = 0; i < nearestNeighbors[ current ].size(); ++i ) {
      size_t next = nearestNeighbors[ current ][ i ];
      if ( !added[ next ] &&
           distances[ current ][ next ] < tolerance &&
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
vector< size_t > getBruteForceTour_( const vector< vector< double > >& distances )
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

vector< vector< size_t > > INLINE_ATTRIBUTE computeNearestNeighbors_( const vector< vector< double > >& distances,
                                                                      size_t numberOfNeighbors )
{
  numberOfNeighbors = min( numberOfNeighbors, distances.size() - 1 );
  vector< vector< size_t > > nearestNeighbors( distances.size(), vector< size_t >( numberOfNeighbors ) );
  set< pair< double, size_t > > tmpNeighbors;
  for ( size_t i = 0; i < distances.size(); ++i ) {
    tmpNeighbors.clear();
    for ( size_t j = 0; j < distances[ i ].size(); ++j ) {
      if( j != i ) {
        if ( tmpNeighbors.size() < numberOfNeighbors || distances[ i ][ j ] < tmpNeighbors.rbegin()->first ) {
          if ( tmpNeighbors.size() >= numberOfNeighbors ) {
            set< pair< double, size_t > >::iterator itLast = tmpNeighbors.end();
            --itLast;
            tmpNeighbors.erase( itLast );
          }
          tmpNeighbors.insert( make_pair( distances[ i ][ j ], j ) );
        }
      }
    }
    set< pair< double, size_t > >::const_iterator it = tmpNeighbors.begin();
    for ( size_t j = 0; j < numberOfNeighbors; ++j, ++it ) {
      // Take neighbor j + 1 in order to avoid adding self as neighbor
      nearestNeighbors[ i ][ j ] = it->second;
    }
  }
  return nearestNeighbors;
}

double getDistance_( const vector< vector< double > >& distances,
                     const vector< double >& lagrangeMultipliers,
                     size_t i,
                     size_t j )
{
  return distances[ i ][ j ] + lagrangeMultipliers[ i ] + lagrangeMultipliers[ j ];
}

struct Vertex {
  Vertex( size_t idx ) : nodeIndex( idx ) {}
  size_t nodeIndex;
  size_t parentIndexInVertexList;
};

double INLINE_ATTRIBUTE compute1Tree_( vector< Vertex >& vertices,
                                       vector< pair< size_t, size_t > >& edges,
                                       vector< size_t >& vertexDegrees,
                                       const vector< vector< double > >& distances,
                                       const vector< double >& lagrangeMultipliers )
{
  vertices.clear();
  vertices.reserve( distances.size() );
  edges.clear();
  edges.reserve( distances.size() );
  vertexDegrees.assign( distances.size(), 0 );

  // 1. Compute length of the minimum spanning tree excluding one vertex, using Prim's algorithm.

  // Select the vertex with minimum maximal distances to a neighbor in
  // order to allow for problems with artificial nodes (such as those used
  // to construct Hamiltonian tours) with toleranceilon distance to all neighbors.
  size_t minimaxVertex = 0;
  double minimaxValue = *max_element( distances[ 0 ].begin(), distances[ 0 ].end() );
  for ( size_t i = 1; i < distances.size(); ++i ) {
    double value = *max_element( distances[ i ].begin(), distances[ i ].end() );
    if ( value < minimaxValue ) {
      minimaxValue = value;
      minimaxVertex = i;
    }
  }

  vertices.push_back( Vertex( minimaxVertex ) );

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
  vertices.push_back( Vertex( rootVertex ) );

  // For each unused vertex i, closestTreeNode[ i ] points to the vertex in the tree which is closest to i.
  vector< size_t > closestTreeNode( distances.size(), 1 ); // 1 is index of rootVertex in vertex list

  double treeLength = 0.0;
  while ( !unusedVertices.empty() ) {
    size_t indexOfClosestTreeNode = 0;
    vector< size_t >::iterator closestUnusedVertexIt = unusedVertices.begin();
    double minDistance = numeric_limits< double >::max();
    for ( vector< size_t >::iterator it = unusedVertices.begin(); it != unusedVertices.end(); ++it ) {
      size_t indexOfTreeNode = vertices[ closestTreeNode[ *it ] ].nodeIndex;
      double distance = getDistance_( distances, lagrangeMultipliers, indexOfTreeNode, *it );
      if ( distance < minDistance ) {
        minDistance = distance;
        indexOfClosestTreeNode = indexOfTreeNode;
        closestUnusedVertexIt = it;
      }
    }
    treeLength += minDistance;
    size_t indexOfNewTreeNode = *closestUnusedVertexIt;
    ++vertexDegrees[ indexOfClosestTreeNode ];
    ++vertexDegrees[ indexOfNewTreeNode ];
    vertices.push_back( Vertex( indexOfNewTreeNode ) );
    vertices.back().parentIndexInVertexList = closestTreeNode[ *closestUnusedVertexIt ];
    edges.push_back( make_pair( indexOfClosestTreeNode, indexOfNewTreeNode ) );
    *closestUnusedVertexIt = unusedVertices.back();
    unusedVertices.pop_back();

    for ( vector< size_t >::iterator it = unusedVertices.begin(); it != unusedVertices.end(); ++it ) {
      size_t indexOfTreeNode = vertices[ closestTreeNode[ *it ] ].nodeIndex;
      double oldDistance = getDistance_( distances, lagrangeMultipliers, indexOfTreeNode, *it );
      double newDistance = getDistance_( distances, lagrangeMultipliers, indexOfNewTreeNode, *it );
      if ( newDistance < oldDistance ) {
        closestTreeNode[ *it ] = vertices.size() - 1;
      }
    }
  }

  // 2. Add the two shortest edges connecting to the first vertex
  double minElement = numeric_limits< double >::max();
  double secondMinElement = numeric_limits< double >::max();
  size_t minElementIndex = 0;
  size_t secondMinElementIndex = 0;
  for ( size_t i = 0; i < distances[ minimaxVertex ].size(); ++i ) {
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

double INLINE_ATTRIBUTE getHeldKarpLowerBound_( const vector< vector< double > >& distances, vector< double >& lagrangeMultipliers )
{
  double maxLength = numeric_limits< double >::lowest();
  lagrangeMultipliers.assign( distances.size(), 0.0 );
  vector< Vertex > vertices;
  vector< pair< size_t, size_t > > edges;
  vector< size_t > vertexDegrees;
  double lambda = 0.1;
  for ( size_t i = 0; i < 50; ++i ) {
    vector< size_t > vertexDegrees;
    double treeLength = compute1Tree_( vertices, edges, vertexDegrees, distances, lagrangeMultipliers );
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

double INLINE_ATTRIBUTE getHeldKarpLowerBound_( const vector< vector< double > >& distances )
{
  vector< double > lagrangeMultipliers( distances.size() );
  return getHeldKarpLowerBound_( distances, lagrangeMultipliers );
}

vector< vector< size_t > > INLINE_ATTRIBUTE computeHelsgaunNeighbors_( const vector< vector< double > >& distances,
                                                                       vector< vector< double > >& alphaDistances,
                                                                       size_t numberOfNeighbors )
{
  vector< Vertex > vertices;
  vector< pair< size_t, size_t > > edges;
  vector< size_t > vertexDegrees;
  vector< double > lagrangeMultipliers( distances.size(), 0.0 );
  getHeldKarpLowerBound_( distances, lagrangeMultipliers );
  compute1Tree_( vertices,
                 edges,
                 vertexDegrees,
                 distances,
                 lagrangeMultipliers );

  // insert (i,j) in T, this creates a cycle containing (i,j) in the spanning tree part of T. Then T(i,j) is obtained by removing the longest of the other edges in this cycle.
  vector< vector< double > > beta( distances.size(), vector< double >( distances.size() ) );
  for ( size_t i = 1; i < vertices.size(); ++i ) {
    size_t iInd = vertices[ i ].nodeIndex;
    beta[ iInd ][ iInd ] = numeric_limits< double >::lowest();
    for ( size_t j = i + 1; j < vertices.size(); ++j ) {
      size_t jInd = vertices[ j ].nodeIndex;
      size_t jParentInd = vertices[ vertices[ j ].parentIndexInVertexList ].nodeIndex;
      beta[ iInd ][ jInd ] = beta[ jInd ][ iInd ] = max( beta[ iInd ][ jParentInd ], getDistance_( distances, lagrangeMultipliers, jInd, jParentInd ) );
    }
  }
  alphaDistances.assign( distances.size(), vector< double >( distances.size() ) );
  for ( size_t i = 0; i < beta.size(); ++i ) {
    for ( size_t j = 0; j < beta[ i ].size(); ++j ) {
      alphaDistances[ i ][ j ] = getDistance_( distances, lagrangeMultipliers, i, j ) - beta[ i ][ j ];
    }
  }

  // if i or j has minimaxNode as end node, then T(i,j) is obtained from T by replacing the longest of the two edges of T incident to minimaxNode with (i,j)
  double maxEdgeWeight = max( getDistance_( distances, lagrangeMultipliers, edges[ edges.size() - 2 ].first, edges[ edges.size() - 2 ].second ),
                              getDistance_( distances, lagrangeMultipliers, edges[ edges.size() - 1 ].first, edges[ edges.size() - 1 ].second ) );
  for ( size_t i = 0; i < distances.size(); ++i ) {
    alphaDistances[ vertices.front().nodeIndex ][ i ] = alphaDistances[ i ][ vertices.front().nodeIndex ] = getDistance_( distances, lagrangeMultipliers, i, vertices.front().nodeIndex ) - maxEdgeWeight;
  }
  alphaDistances[ vertices.front().nodeIndex ][ vertices.front().nodeIndex ] = 0.0;

  // if (i,j) belongs to T, then T(i,j) is equal to T
  for ( size_t i = 0; i < edges.size(); ++i ) {
    alphaDistances[ edges[ i ].first ][ edges[ i ].second ] = alphaDistances[ edges[ i ].second ][ edges[ i ].first ] = 0.0;
  }

  return computeNearestNeighbors_( alphaDistances, numberOfNeighbors );
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

void INLINE_ATTRIBUTE flip_( vector< size_t >& tour,
                             vector< size_t >& position,
                             size_t t1,
                             size_t t2,
                             size_t t3,
                             size_t t4 )
{

  if ( !( t2 == next_( t1, tour, position ) && t3 == next_( t4, tour, position ) ) ) {
    cerr << "t: " << t1 << " " << t2 << " " << t3 << " " << t4 << endl;
    cerr << "tpos: " << position[ t2 ] << " " << position[ t1 ] << " " << position[ t3 ] << " " << position[ t4 ] << endl;
  }
  assert( t2 == next_( t1, tour, position ) );
  assert( t3 == next_( t4, tour, position ) );
  size_t pt2 = position[ t2 ];
  size_t pt3 = position[ t3 ];

  if ( pt3 < pt2 ) {
    swap( pt2, pt3 );
  }
  if ( pt3 - pt2 <= tour.size() - pt3 + pt2 ) {
    reverse( tour.begin() + pt2, tour.begin() + pt3 );
  }
  else {
    vector< size_t >::iterator iIt = tour.begin() + pt3;
    vector< size_t >::iterator jIt = tour.begin() + pt2;
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

void INLINE_ATTRIBUTE performMove_( const vector< size_t >& ts,
                                    vector< size_t >& tour,
                                    vector< size_t >& position )
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
                                        vector< size_t >& tour,
                                        vector< size_t >& position )
{
  vector< size_t > p, q, inclTmp;
  getPQI_( p, q, inclTmp, ts, tour, position );

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
                                        vector< size_t >& tour,
                                        vector< size_t >& position )
{
  vector< size_t > p, q, incl;
  getPQI_( p, q, incl, ts, tour, position );
  performKOptMove_( ts, incl, tour, position );
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

bool INLINE_ATTRIBUTE twoMakesTour_( const vector< size_t >& ts,
                                     const vector< size_t >& ss,
                                    const vector< size_t >& tour,
                                    const vector< size_t >& position )
{
  assert( ts.size() > 1 );
  vector< size_t > p, q, incl, incl2, ts2;

  getPQI_( p, q, incl, ts, tour, position );

  ts2 = ss;
  ts2.insert( ts2.end(), ts.begin(), ts.end() );

  getPQI_( p, q, incl2, ts2, tour, position );
  for ( size_t i = 0; i < incl.size(); ++i ) {
    incl[ i ] += ss.size();
  }
  if ( ss.size() > 0 ) {
    assert( ss.size() == 4 || ss.size() == 6 );
    if ( ss.size() == 4 ) {
      incl2 = { 3, 2, 1, 0 };
    }
    else {
      incl2 = { 5, 2, 1, 4, 3, 0 };
    }
    incl2.insert( incl2.end(), incl.begin(), incl.end() );
    incl = incl2;
  }


  size_t count = 1;
  for ( size_t i = ts2.size(); ( i = ( q[ incl[ p[ i - 1 ] ] ] + 1 ) ^ 1 ) != 0; ++count );
  return 2 * count == ts2.size();
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

bool INLINE_ATTRIBUTE improveTour2Opt_( vector< size_t >& tour,
                                        const vector< vector< double > >& distances,
                                        const vector< vector< size_t > >& nearestNeighbors )
{
  vector< size_t > position( tour.size() );
  for ( size_t i = 0; i < tour.size(); ++i ) {
    position[ tour[ i ] ] = i;
  }

  bool anyChange = false;
  vector< size_t > bestTs;
  bool changed = true;
  while ( changed ) {
    changed = false;
    for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
      double maxGain = 0.0;
      for ( size_t t2choice = 0; t2choice < 2; ++t2choice  ) {
        size_t t2 = t2choice == 0 ? previous_( t1, tour, position ) : next_( t1, tour, position );
        for ( size_t t3index = 0; t3index < nearestNeighbors[ t2 ].size(); ++t3index ) {
          size_t t3 = nearestNeighbors[ t2 ][ t3index ];
          if ( t3 == previous_( t2, tour, position ) || t3 == next_( t2, tour, position ) ) {
            continue;
          }
          if ( distances[ t2 ][ t3 ] >= distances[ t1 ][ t2 ] ) {
            continue;
          }
          size_t t4 = t2choice == 0 ? next_( t3, tour, position ) : previous_( t3, tour, position );
          double gain = distances[ t1 ][ t2 ] + distances[ t3 ][ t4 ] - distances[ t1 ][ t4 ] - distances[ t2 ][ t3 ];
          if ( gain > maxGain ) {
            maxGain = gain;
            bestTs = { t1, t2, t3, t4 };
          }
        }
      }
      if ( maxGain > tolerance ) {
        performMove_( bestTs, tour, position );
        assert( isTour_( tour, position ) );
        anyChange = true;
        changed = true;
      }
    }
  }
  return anyChange;
}

bool INLINE_ATTRIBUTE threeOptInnerLoop_( vector< size_t >& tour,
                                          vector< bool >& dontLook,
                                          const vector< vector< double > >& distances,
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
        double g1 = distances[ t1 ][ t2 ] - distances[ t2 ][ t3 ];
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
          if ( g1 + g2 <= tolerance ) {
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
          if ( g1 + g2 <= tolerance ) {
            ;
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
    if ( G > tolerance ) {
      if ( !dontLook.empty() ) {
        for ( size_t i = 0; i < bestTs.size(); ++i ) {
          dontLook[ bestTs[ i ] ] = false;
        }
      }
      performMove_( bestTs, tour, position );
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
                                        const vector< vector< double > >& distances,
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
                                               const vector< vector< double > >& distances,
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
        double gainFirstBridge = distances[ t1 ][ t2 ] + distances[ t3 ][ t4 ] - distances[ t2 ][ t3 ] - distances[ t1 ][ t4 ];
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
                  t8choice == t6choice ? distances[ t5 ][ t6 ] + distances[ t7 ][ t8 ] - distances[ t6 ][ t7 ] - distances[ t5 ][ t8 ]
                                       : distances[ t5 ][ t6 ] + distances[ t7 ][ t8 ] - distances[ t5 ][ t7 ] - distances[ t6 ][ t8 ];
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
                    distances[ t5 ][ t6 ] + distances[ t7 ][ t8 ] + distances[ t9 ][ t10 ] - distances[ t6 ][ t7 ] - distances[ t8 ][ t9 ] - distances[ t10 ][ t5 ];
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
      /*
      cerr << "23-ts:" << endl;
      for ( size_t i = 0; i < bestTs.size(); ++i ) {
        cerr << bestTs[i] << " ";
      }
      cerr << endl;
      cerr << "incl ts:" << endl;
      for ( size_t i = 0; i < incl.size(); ++i ) {
        cerr << incl[ i ]<< " ";
      }
      cerr << endl;
      */
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
                                      const vector< vector< double > >& distances,
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

bool INLINE_ATTRIBUTE linKernighanInnerLoop_( vector< size_t >& tour,
                                              vector< size_t >& position,
                                              vector< vector< size_t > >& added,
                                              vector< vector< size_t > >& removed,
                                              vector< size_t >& tcUntestedInStep2,
                                              size_t ta,
                                              size_t tb,
                                              double G,
                                              const size_t t1,
                                              const double lengthBefore,
                                              const vector< vector< double > >& distances,
                                              const vector< vector< size_t > >& nearestNeighbors )
{
  assert( added.size() == 2 || added.size() == 4 );
  assert( removed.size() == 4 || removed.size() == 6 );
  bool positiveGain = true;
  size_t numberOfEdgesToRemove = 2;
  while ( positiveGain ) {
    positiveGain = false;
    vector< size_t > mutableNearestNeighborList( nearestNeighbors[ tb ] );
    vector< size_t >& tcUntested = numberOfEdgesToRemove == 2 ? tcUntestedInStep2 : mutableNearestNeighborList;
    assert( tb == next_( t1, tour, position ) || tb == previous_( t1, tour, position ) );
    const bool tBchoice = tb == previous_( t1, tour, position );

    while ( !tcUntested.empty() ) {
      // Select the most promising untested candidate for tc
      size_t tc = 0, td = 0;
      double bestDiff = numeric_limits< double >::lowest();
      for ( size_t tcIndex = 0; tcIndex < tcUntested.size(); ++tcIndex ) {
        size_t testTc = tcUntested[ tcIndex ];
        size_t testTd = tBchoice ? next_( testTc, tour, position ) : previous_( testTc, tour, position );
        double diff = distances[ testTc ][ testTd ] - distances[ tb ][ testTc ];
        if ( diff > bestDiff ) {
          bestDiff = diff;
          tc = testTc;
          td = testTd;
        }
      }
      vector< size_t >::iterator it = find( tcUntested.begin(), tcUntested.end(), tc );
      assert( it != tcUntested.end() );
      *it = tcUntested.back();
      tcUntested.pop_back();

      double gn = distances[ ta ][ tb ] - distances[ tb ][ tc ];
      if ( G + gn <= tolerance ) {
        continue;
      }

      vector< size_t > bc( { tb, tc } );
      vector< size_t > cd( { tc, td } );
      if ( find( removed.begin(), removed.end(), bc ) != removed.end() ) {
        continue;
      }
      if ( find( added.begin(), added.end(), cd ) != added.end() ) {
        continue;
      }

      G += gn;
      performMove_( { t1, tb, tc, td }, tour, position );
      added.insert( added.end(), { { tb, tc }, { tc, tb } } );
      removed.insert( removed.end(), { { tc, td }, { td, tc } } );

      // If connecting back to t1 leads to an improvement, then take it!
      double gain = G + distances[ tc ][ td ] - distances[ td ][ t1 ];
      if ( gain > tolerance ) {
        added.insert( added.end(), { { t1, td }, { td, t1 } } );
        assert( isTour_( tour, position ) );
        return true;
      }

      ta = tc;
      tb = td;
      ++numberOfEdgesToRemove;
      positiveGain = true;
      break;
    }
  }
  return false;
}

bool INLINE_ATTRIBUTE linKernighanOuterLoop_( vector< size_t >& tour,
                                              vector< bool >& dontLook,
                                              const vector< vector< double > >& distances,
                                              const vector< vector< size_t > >& nearestNeighbors )
{
  vector< size_t > position( tour.size() );
  for ( size_t i = 0; i < tour.size(); ++i ) {
    position[ tour[ i ] ] = i;
  }
  const vector< size_t > tourCopy( tour );
  const vector< size_t > positionCopy( position );
  const double lengthBefore = getLength_( tour, distances );
  for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
    if ( dontLook[ t1 ] ) {
      continue;
    }
    double G = 0.0;
    for ( size_t t2choice = 0; t2choice < 2; ++t2choice ) {
      size_t t2 = t2choice == 0 ? previous_( t1, tour, position ) : next_( t1, tour, position );
      for ( size_t t3index = 0; t3index < nearestNeighbors[ t2 ].size(); ++t3index ) {
        size_t t3 = nearestNeighbors[ t2 ][ t3index ];
        if ( t3 == previous_( t2, tour, position ) || t3 == next_( t2, tour, position ) ) {
          continue;
        }
        double g1 = distances[ t1 ][ t2 ] - distances[ t2 ][ t3 ];
        if ( g1 <= tolerance ) {
          continue;
        }

        // First choice of t4
        {
          size_t t4 = t2choice == 0 ? next_( t3, tour, position ) : previous_( t3, tour, position );
          performMove_( { t1, t2, t3, t4 }, tour, position );
          // Test for improving 2-opt move
          double gain = g1 + distances[ t3 ][ t4 ] - distances[ t4 ][ t1 ];
          if ( gain > tolerance ) {
            assert( getLength_( tour, distances ) < lengthBefore );
            assert( isTour_( tour, position ) );
            return true;
          }

          G = g1;

          // Save the state after x2 to enable backtracking
          const double G2 = G;
          const vector< size_t > tour2( tour );
          const vector< size_t > position2( position );
          const size_t ta2 = t3;
          const size_t tb2 = t4;

          vector< size_t > tcUntestedInStep2 = nearestNeighbors[ tb2 ];

          // While loop corresponds to backtracking over x2 (second removed edge)
          while ( !tcUntestedInStep2.empty() ) {
            // Reset state to that after x0
            vector< vector< size_t > > added( { { t2, t3 }, { t3, t2 } } );
            vector< vector< size_t > > removed( { { t1, t2 }, { t2, t1 }, { t3, t4 }, { t4, t3 } } );
            tour = tour2;
            position = position2;
            if ( linKernighanInnerLoop_( tour,
                                         position,
                                         added,
                                         removed,
                                         tcUntestedInStep2,
                                         ta2,
                                         tb2,
                                         G2,
                                         t1,
                                         lengthBefore,
                                         distances,
                                         nearestNeighbors ) ) {
              for ( size_t i = 0; i < added.size(); i += 2 ) {
                dontLook[ added[ i ][ 0 ] ] = false;
                dontLook[ added[ i ][ 1 ] ] = false;
              }
              return true;
            }
          }
          // Found no gaining move. Reset tour
          tour = tourCopy;
          position = positionCopy;
        }

        if ( true ) {
          // Second choice of t4
          size_t t4 = t2choice == 0 ? previous_( t3, tour, position ) : next_( t3, tour, position );
          vector< size_t > t5Untested = nearestNeighbors[ t4 ];
          while ( !t5Untested.empty() ) {
            size_t t5 = 0;
            size_t t6 = 0;
            double bestDiff = numeric_limits< double >::lowest();
            for ( size_t t5index = 0; t5index < t5Untested.size(); ++t5index ) {
              size_t t5tmp = t5Untested[ t5index ];

              // Only consider one choice of t6tmp. The other choice is possible, but clutters the code and doesn't lead to a significant improvement.
              size_t t6choice = t2choice;
              size_t t6tmp = t6choice == 0 ? next_( t5tmp, tour, position ) : previous_( t5tmp, tour, position );
              double diff = distances[ t5tmp ][ t6tmp ] - distances[ t4 ][ t5tmp ];
              if ( diff > bestDiff ) {
                bestDiff = diff;
                t5 = t5tmp;
                t6 = t6tmp;
              }
            }
            vector< size_t >::iterator it = find( t5Untested.begin(), t5Untested.end(), t5 );
            assert( it != t5Untested.end() );
            *it = t5Untested.back();
            t5Untested.pop_back();

            if ( t5 == t3 || t5 == t2 ) {
              continue;
            }
            if ( t6 == t3 || t6 == t2 || t6 == t1 ) {
              continue;
            }
            if ( ( t2choice == 0 && !between_( t3, t2, t5, position ) ) ||
                ( t2choice == 1 && !between_( t2, t3, t5, position ) ) ) {
              continue;
            }

            double g2 = distances[ t3 ][ t4 ] - distances[ t4 ][ t5 ];
            if ( g1 + g2 <= tolerance ) {
              continue;
            }
            G = g1 + g2;

            performMove_( { t3, t4, t5, t6, t1, t2 }, tour, position );

            // Save the state after x2 to enable backtracking
            const double G2 = G;
            const vector< size_t > tour2( tour );
            const vector< size_t > position2( position );
            const size_t ta2 = t5;
            const size_t tb2 = t6;

            vector< size_t > tcUntestedInStep2 = nearestNeighbors[ tb2 ];

            // While loop corresponds to backtracking over x2 (second removed edge)
            while ( !tcUntestedInStep2.empty() ) {
              // Reset state to that after x0
              vector< vector< size_t > > added( { { t2, t3 }, { t3, t2 }, { t4, t5 }, { t5, t4 } } );
              vector< vector< size_t > > removed( { { t1, t2 }, { t2, t1 }, { t3, t4 }, { t4, t3 }, { t5, t6 }, { t6, t5 } } );
              tour = tour2;
              position = position2;
              if ( linKernighanInnerLoop_( tour,
                                           position,
                                           added,
                                           removed,
                                           tcUntestedInStep2,
                                           ta2,
                                           tb2,
                                           G2,
                                           t1,
                                           lengthBefore,
                                           distances,
                                           nearestNeighbors ) ) {
                for ( size_t i = 0; i < added.size(); i += 2 ) {
                  dontLook[ added[ i ][ 0 ] ] = false;
                  dontLook[ added[ i ][ 1 ] ] = false;
                }
                return true;
              }
            }
            // Found no gaining move. Reset tour
            tour = tourCopy;
            position = positionCopy;
          }
        }
      }
    }
    dontLook[ t1 ] = true;
  }
  return false;
}

bool INLINE_ATTRIBUTE improveTourLinKernighan_( vector< size_t >& tour,
                                                const vector< vector< double > >& distances,
                                                const vector< vector< size_t > >& nearestNeighbors )
{
  bool change = false;
  vector< bool > dontLook( tour.size(), false );
  while ( linKernighanOuterLoop_( tour, dontLook, distances, nearestNeighbors ) ) {
    change = true;
  }
  return change;
}


/*
bool INLINE_ATTRIBUTE removed_( const vector< size_t > ts, size_t ta, size_t tb )
{
  for ( size_t i = 0; i < ts.size(); i += 2 ) {
    if ( ( ta == ts[ i ] && tb == ts[ i + 1 ] ) ||
         ( ta == ts[ i + 1 ] && tb == ts[ i ] ) ) {
      return true;
    }
  }
  return false;
}
*/
/*
  bool INLINE_ATTRIBUTE added_( const vector< size_t > ts, size_t ta, size_t tb ) {
  for ( size_t i = 1; i + 1 < ts.size(); i += 2 ) {
    if ( ( ta == ts[ i ] && tb == ts[ i + 1 ] ) ||
         ( ta == ts[ i + 1 ] && tb == ts[ i ] ) ) {
      return true;
    }
  }
  return false;
}
*/

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
                                       const vector< vector< double > >& distances,
                                       const vector< vector< size_t > >& nearestNeighbors )
{
  assert( ts.size() > 1 );
  const size_t tb = ts.back();
  vector< pair< size_t, size_t > > tcTdPairs;
  for ( size_t tcIndex = 0; tcIndex < nearestNeighbors[ tb ].size(); ++tcIndex ) {
    size_t tc = nearestNeighbors[ tb ][ tcIndex ];
    double G1 = G0 - distances[ tb ][ tc ];
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
        double G2 = G1 + distances[ tc ][ td ];
        if ( G2 - distances[ ts.front() ][ td ] <= tolerance && ( G2 <= bestG && G2 <= bestInfeasibleG ) ) {
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
  sort( tcTdPairs.begin(), tcTdPairs.end(), [&] ( const pair< size_t, size_t >& p1,
                                                  const pair< size_t, size_t >& p2 ) {
    const size_t tc1 = p1.first;
    const size_t td1 = p1.second;
    const size_t tc2 = p2.first;
    const size_t td2 = p2.second;
    return distances[ tc1 ][ td1 ] - distances[ tb ][ tc1 ]
           < distances[ tc2 ][ td2 ] - distances[ tb ][ tc2 ];
  } );

  while ( !tcTdPairs.empty() ) {
    size_t tc = tcTdPairs.back().first;
    size_t td = tcTdPairs.back().second;
    tcTdPairs.pop_back();

    ts.push_back( tc );
    ts.push_back( td );
    double G1 = G0 - distances[ tb ][ tc ];
    double G2 = G1 + distances[ tc ][ td ];
    double gain = G2 - distances[ ts.front() ][ ts.back() ];
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

double INLINE_ATTRIBUTE twoBestKOptMove_( size_t depth,
                                       size_t k,
                                       const vector< size_t >& ss,
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
                                       const vector< vector< double > >& distances,
                                       const vector< vector< size_t > >& nearestNeighbors )
{
  assert( ts.size() > 1 );
  const size_t tb = ts.back();
  vector< pair< size_t, size_t > > tcTdPairs;
  for ( size_t tcIndex = 0; tcIndex < nearestNeighbors[ tb ].size(); ++tcIndex ) {
    size_t tc = nearestNeighbors[ tb ][ tcIndex ];
    double G1 = G0 - distances[ tb ][ tc ];
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
        double G2 = G1 + distances[ tc ][ td ];
        if ( G2 - distances[ ts.front() ][ td ] <= tolerance && ( G2 <= bestG && G2 <= bestInfeasibleG ) ) {
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
  sort( tcTdPairs.begin(), tcTdPairs.end(), [&] ( const pair< size_t, size_t >& p1,
                                                  const pair< size_t, size_t >& p2 ) {
    const size_t tc1 = p1.first;
    const size_t td1 = p1.second;
    const size_t tc2 = p2.first;
    const size_t td2 = p2.second;
    return distances[ tc1 ][ td1 ] - distances[ tb ][ tc1 ]
           < distances[ tc2 ][ td2 ] - distances[ tb ][ tc2 ];
  } );

  while ( !tcTdPairs.empty() ) {
    size_t tc = tcTdPairs.back().first;
    size_t td = tcTdPairs.back().second;
    tcTdPairs.pop_back();

    ts.push_back( tc );
    ts.push_back( td );
    double G1 = G0 - distances[ tb ][ tc ];
    double G2 = G1 + distances[ tc ][ td ];
    double gain = G2 - distances[ ts.front() ][ ts.back() ];
    if ( gain > tolerance && twoMakesTour_( ts, ss, tour, position ) ) {
      return gain;
    }

    if ( depth + 1 < k ) {
      added.push_back( make_pair( tb, tc ) );
      added.push_back( make_pair( tc, tb ) );
      removed.push_back( make_pair( tc, td ) );
      removed.push_back( make_pair( td, tc ) );
      gain = twoBestKOptMove_( depth + 1,
                            k,
                            ss,
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
    else if ( G2 > bestG && twoMakesTour_( ts, ss, tour, position ) ) {
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
                                      const vector< vector< double > >& distances,
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

      vector< size_t > ts( { t1, t2 } );
      double G0 = distances[ t1 ][ t2 ];
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
        bestInfeasibleG = tolerance;
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
      } while ( ( bestG > tolerance && lkDepth < 30 ) || ( bestInfeasibleG > tolerance && ts.size() < 2 * k + 1 ) );
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

bool INLINE_ATTRIBUTE twokOptOuterLoop_( size_t k,
                                      vector< size_t >& tour,
                                      vector< bool >& dontLook,
                                      vector< bool >& otherDontLook,
                                      const vector< vector< double > >& distances,
                                      const vector< vector< size_t > >& nearestNeighbors,
                                      bool linKernighan )
{
  vector< size_t > position( tour.size() );
  for ( size_t i = 0; i < tour.size(); ++i ) {
    position[ tour[ i ] ] = i;
  }
  bool anyChange = false;

  for ( size_t s1 = 0; s1 < tour.size(); ++s1 ) {
    if ( !dontLook.empty() && dontLook[ s1 ] ) {
      continue;
    }
    bool found = false;
apa:
    for ( size_t s2choice = 0; s2choice < 2; ++s2choice ) {
      size_t s2 = s2choice == 0 ? next_( s1, tour, position ) : previous_( s1, tour, position );

      for ( size_t s3index = 0; s3index < nearestNeighbors[ s2 ].size(); ++s3index ) {
        size_t s3 = nearestNeighbors[ s2 ][ s3index ];
        size_t s4 = s2choice == 0 ? next_( s3, tour, position ) : previous_( s3, tour, position );
        if ( s3 == s1 || s4 == s1 ) {
          continue;
        }
        double gainFirstBridge = distances[ s1 ][ s2 ] + distances[ s3 ][ s4 ] - distances[ s2 ][ s3 ] - distances[ s1 ][ s4 ];
        if ( gainFirstBridge <= tolerance ) {
          continue;
        }

        vector< size_t > ss( { s1, s2, s3, s4 } );
  for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
    for ( size_t t2choice = 0; t2choice < 1; ++t2choice ) {
      size_t t2 = t2choice == 0 ? previous_( t1, tour, position ) : next_( t1, tour, position );

      vector< size_t > ts( { t1, t2 } );
      double G0 = gainFirstBridge + distances[ t1 ][ t2 ];
      vector< pair< size_t, size_t > > added( { make_pair( s2, s3 ), make_pair( s3, s2 ), make_pair( s1, s4 ), make_pair( s4, s1 ) } );
      vector< pair< size_t, size_t > > removed( { make_pair( s1, s2 ), make_pair( s2, s1 ), make_pair( s3, s4 ), make_pair( s4, s3 ) } );
      if ( find( removed.begin(), removed.end(), make_pair( t1, t2 ) ) != removed.end() ) {
        continue;
      }
      removed.push_back( make_pair( t1, t2 ) );
      removed.push_back( make_pair( t2, t1 ) );
      vector< size_t > bestTs;
      double bestG = tolerance;
      vector< size_t > bestInfeasibleTs;
      double bestInfeasibleG = tolerance;
      bool testChange = false;
      vector< size_t > tourCopy( tour );
      vector< size_t > positionCopy( position );
      vector< size_t > tsHistory;
      size_t lkDepth = 0;
      vector< size_t > ss2;
      do {
        ++lkDepth;
        bestG = tolerance;
        bestInfeasibleG = tolerance;
        double gain = twoBestKOptMove_( 1,
                                     k,
                                     lkDepth == 1 ? ss : ss2,
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
        cerr << "gain " << lkDepth << ". " << s1 << " " << s2 << " " << s3 << " " << s4 << endl;
          vector< size_t > p, q, incl, incl2, ts2;
          if ( lkDepth < 2 ) {
            getPQI_( p, q, incl, ts, tour, position );
            for ( size_t i = 0; i < incl.size(); ++i ) {
              incl[ i ] += ss.size();
            }
            incl2 = { 3, 2, 1, 0 };
            incl2.insert( incl2.end(), incl.begin(), incl.end() );
            incl = incl2;
            ts2 = ss;
            ts2.insert( ts2.end(), ts.begin(), ts.end() );
            performKOptMove_( ts2, incl, tour, position );
          }
          else {
            performKOptMove_( ts, tour, position );
          }

          if ( !isTour_( tour, position ) ) {
            cerr << "ts: ";
            for ( size_t j = 0; j < ts2.size(); ++j ) {
              cerr << ts2[ j ] << " ";
            }
            cerr << endl;

            cerr << "incl: ";
            for ( size_t j = 0; j < incl.size(); ++j ) {
              cerr << incl[ j ] << " ";
            }
            cerr << endl;
          }

          assert( isTour_( tour, position ) );
          tsHistory.insert( tsHistory.end(), ts.begin(), ts.end() );
          for ( size_t i = 0; i < tsHistory.size(); ++i ) {
            dontLook[ tsHistory[ i ] ] = false;
            otherDontLook[ tsHistory[ i ] ] = false;
          }
          found = true;
          anyChange = true;
          testChange = false;
          goto apa;
          break;
        }
        if ( !linKernighan ) {
          break;
        }
        bestInfeasibleG = 0.0;
        if ( bestG > tolerance ) {
          if ( lkDepth < 2 ) {
            vector< size_t > p, q, ts2, incl, incl2;
            getPQI_( p, q, incl, bestTs, tour, position );
            for ( size_t i = 0; i < incl.size(); ++i ) {
              incl[ i ] += ss.size();
            }
            incl2 = { 3, 2, 1, 0 };
            incl2.insert( incl2.end(), incl.begin(), incl.end() );
            incl = incl2;
            ts2 = ss;
            ts2.insert( ts2.end(), bestTs.begin(), bestTs.end() );
            tsHistory.insert( tsHistory.end(), ts2.begin(), ts2.end() );

            performKOptMove_( ts2, incl, tour, position );
          }
          else {
            performKOptMove_( bestTs, tour, position );
            tsHistory.insert( tsHistory.end(), bestTs.begin(), bestTs.end() );
          }
          testChange = true;
          G0 = bestG;
          ts = { bestTs.front(), bestTs.back() };
          added = { make_pair( s2, s3 ), make_pair( s3, s2 ), make_pair( s1, s4 ), make_pair( s4, s1 ) };
          removed = { make_pair( s1, s2 ), make_pair( s2, s1 ), make_pair( s3, s4 ), make_pair( s4, s3 ) };
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
          added = { make_pair( s2, s3 ), make_pair( s3, s2 ), make_pair( s1, s4 ), make_pair( s4, s1 ) };
          removed = { make_pair( s1, s2 ), make_pair( s2, s1 ), make_pair( s3, s4 ), make_pair( s4, s3 ) };
          for ( size_t i = 0; i + 1 < ts.size(); i += 2 ) {
            removed.push_back( make_pair( ts[ i ], ts[ i + 1 ] ) );
            removed.push_back( make_pair( ts[ i + 1 ], ts[ i ] ) );
          }
          for ( size_t i = 1; i + 1 < ts.size(); i += 2 ) {
            added.push_back( make_pair( ts[ i ], ts[ i + 1 ] ) );
            added.push_back( make_pair( ts[ i + 1 ], ts[ i ] ) );
          }
        }
      } while ( lkDepth < 35 && ( bestG > tolerance || ( bestInfeasibleG > tolerance && ts.size() < 2 * k + 1 ) ) );
      if ( testChange ) {
        tour = tourCopy;
        position = positionCopy;
      }
    }
  }
      }
    }
    if ( !found ) {
      dontLook[ s1 ] = true;
    }
  }
  return anyChange;
}

bool INLINE_ATTRIBUTE threekOptOuterLoop_( size_t k,
                                      vector< size_t >& tour,
                                      vector< bool >& dontLook,
                                      vector< bool >& otherDontLook,
                                      const vector< vector< double > >& distances,
                                      const vector< vector< size_t > >& nearestNeighbors,
                                      bool linKernighan )
{
  vector< size_t > position( tour.size() );
  for ( size_t i = 0; i < tour.size(); ++i ) {
    position[ tour[ i ] ] = i;
  }
  bool anyChange = false;

  for ( size_t s1 = 0; s1 < tour.size(); ++s1 ) {
    if ( !dontLook.empty() && dontLook[ s1 ] ) {
      continue;
    }
    bool found = false;
apa:
    for ( size_t s2choice = 0; s2choice < 2; ++s2choice ) {
      size_t s2 = s2choice == 0 ? next_( s1, tour, position ) : previous_( s1, tour, position );

      for ( size_t s3index = 0; s3index < nearestNeighbors[ s2 ].size(); ++s3index ) {
        size_t s3 = nearestNeighbors[ s2 ][ s3index ];
        double g1 = distances[ s1 ][ s2 ] - distances[ s2 ][ s3 ];
        if ( g1 <= tolerance ) {
          continue;
        }
        for ( size_t s4choice = 0; s4choice < 2; ++s4choice ) {
          size_t s4 = s4choice == 0 ? next_( s3, tour, position ) : previous_( s3, tour, position );
          if ( s3 == s1 || s4 == s1 ) {
            continue;
          }
          for ( size_t s5index = 0; s5index < nearestNeighbors[ s4 ].size(); ++s5index ) {
            size_t s5 = nearestNeighbors[ s4 ][ s5index ];
            double g2 = distances[ s3 ][ s4 ] - distances[ s4 ][ s5 ];
            if ( g1 + g2 <= tolerance ) {
              continue;
            }
            for ( size_t s6choice = 0; s6choice < 2; ++s6choice ) {
              size_t s6 = s6choice == 0 ? next_( s5, tour, position ) : previous_( s5, tour, position );
              if ( ( s5 == s1 && s6 == s2 ) ||
                   ( s6 == s1 && s5 == s2 ) ||
                   ( s5 == s3 && s6 == s4 ) ||
                   ( s6 == s3 && s5 == s4 ) ) {
                continue;
              }
              if ( s6 == s1 ) {
                continue;
              }
              double g3 = distances[ s5 ][ s6 ] - distances[ s6 ][ s1 ];
              double gg = g1 + g2 + g3;
              if ( gg <= tolerance ) {
                continue;
              }
              vector< size_t > ss( { s1, s2, s3, s4, s5, s6 } );
              if ( makesTour_( ss, tour, position ) ) {
                continue;
              }




  for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
    for ( size_t t2choice = 0; t2choice < 1; ++t2choice ) {
      size_t t2 = t2choice == 0 ? previous_( t1, tour, position ) : next_( t1, tour, position );

      vector< size_t > ts( { t1, t2 } );
      double G0 = gg + distances[ t1 ][ t2 ];
      vector< pair< size_t, size_t > > added( { make_pair( s2, s3 ), make_pair( s3, s2 ), make_pair( s1, s4 ), make_pair( s4, s1 ) } );
      vector< pair< size_t, size_t > > removed( { make_pair( s1, s2 ), make_pair( s2, s1 ), make_pair( s3, s4 ), make_pair( s4, s3 ) } );
      added = { make_pair( s2, s3 ), make_pair( s3, s2 ), make_pair( s4, s5 ), make_pair( s5, s4 ), make_pair( s6, s1 ), make_pair( s1, s6 ) };
      removed = { make_pair( s1, s2 ), make_pair( s2, s1 ), make_pair( s3, s4 ), make_pair( s4, s3 ), make_pair( s5, s6 ), make_pair( s6, s5 )  };
      if ( find( removed.begin(), removed.end(), make_pair( t1, t2 ) ) != removed.end() ) {
        continue;
      }
      removed.push_back( make_pair( t1, t2 ) );
      removed.push_back( make_pair( t2, t1 ) );
      vector< size_t > bestTs;
      double bestG = tolerance;
      vector< size_t > bestInfeasibleTs;
      double bestInfeasibleG = tolerance;
      bool testChange = false;
      vector< size_t > tourCopy( tour );
      vector< size_t > positionCopy( position );
      vector< size_t > tsHistory;
      size_t lkDepth = 0;
      vector< size_t > ss2;
      do {
        ++lkDepth;
        bestG = tolerance;
        bestInfeasibleG = tolerance;
        double gain = twoBestKOptMove_( 1,
                                     k,
                                     lkDepth == 1 ? ss : ss2,
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
          cerr << "gain " << lkDepth << ". " <<s1 << " " << s2 << " " << s3 << " " << s4 << " " << s5 << " " << s6 << endl;
          vector< size_t > p, q, incl, incl2, ts2;
          if ( lkDepth < 2 ) {
            getPQI_( p, q, incl, ts, tour, position );
            for ( size_t i = 0; i < incl.size(); ++i ) {
              incl[ i ] += ss.size();
            }
            incl2 = { 5, 2, 1, 4, 3, 0 };
            incl2.insert( incl2.end(), incl.begin(), incl.end() );
            incl = incl2;
            ts2 = ss;
            ts2.insert( ts2.end(), ts.begin(), ts.end() );
            performKOptMove_( ts2, incl, tour, position );
          }
          else {
            performKOptMove_( ts, tour, position );
          }

          if ( !isTour_( tour, position ) ) {
            cerr << "ts: ";
            for ( size_t j = 0; j < ts2.size(); ++j ) {
              cerr << ts2[ j ] << " ";
            }
            cerr << endl;

            cerr << "incl: ";
            for ( size_t j = 0; j < incl.size(); ++j ) {
              cerr << incl[ j ] << " ";
            }
            cerr << endl;
          }

          assert( isTour_( tour, position ) );
          tsHistory.insert( tsHistory.end(), ts.begin(), ts.end() );
          for ( size_t i = 0; i < tsHistory.size(); ++i ) {
            dontLook[ tsHistory[ i ] ] = false;
            otherDontLook[ tsHistory[ i ] ] = false;
          }
          found = true;
          anyChange = true;
          testChange = false;
          goto apa;
          break;
        }
        if ( !linKernighan ) {
          break;
        }
        bestInfeasibleG = 0.0;
        if ( bestG > tolerance ) {
          if ( lkDepth < 2 ) {
            vector< size_t > p, q, ts2, incl, incl2;
            getPQI_( p, q, incl, bestTs, tour, position );
            for ( size_t i = 0; i < incl.size(); ++i ) {
              incl[ i ] += ss.size();
            }
            incl2 = { 5, 2, 1, 4, 3, 0 };
            incl2.insert( incl2.end(), incl.begin(), incl.end() );
            incl = incl2;
            ts2 = ss;
            ts2.insert( ts2.end(), bestTs.begin(), bestTs.end() );
            tsHistory.insert( tsHistory.end(), ts2.begin(), ts2.end() );

            performKOptMove_( ts2, incl, tour, position );
          }
          else {
            performKOptMove_( bestTs, tour, position );
            tsHistory.insert( tsHistory.end(), bestTs.begin(), bestTs.end() );
          }
          testChange = true;
          G0 = bestG;
          ts = { bestTs.front(), bestTs.back() };
          added = { make_pair( s2, s3 ), make_pair( s3, s2 ), make_pair( s4, s5 ), make_pair( s5, s4 ), make_pair( s6, s1 ), make_pair( s1, s6 ) };
          removed = { make_pair( s1, s2 ), make_pair( s2, s1 ), make_pair( s3, s4 ), make_pair( s4, s3 ), make_pair( s5, s6 ), make_pair( s6, s5 )  };
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
          added = { make_pair( s2, s3 ), make_pair( s3, s2 ), make_pair( s4, s5 ), make_pair( s5, s4 ), make_pair( s6, s1 ), make_pair( s1, s6 ) };
          removed = { make_pair( s1, s2 ), make_pair( s2, s1 ), make_pair( s3, s4 ), make_pair( s4, s3 ), make_pair( s5, s6 ), make_pair( s6, s5 )  };
          for ( size_t i = 0; i + 1 < ts.size(); i += 2 ) {
            removed.push_back( make_pair( ts[ i ], ts[ i + 1 ] ) );
            removed.push_back( make_pair( ts[ i + 1 ], ts[ i ] ) );
          }
          for ( size_t i = 1; i + 1 < ts.size(); i += 2 ) {
            added.push_back( make_pair( ts[ i ], ts[ i + 1 ] ) );
            added.push_back( make_pair( ts[ i + 1 ], ts[ i ] ) );
          }
        }
      } while ( lkDepth < 35 && ( bestG > tolerance || ( bestInfeasibleG > tolerance && ts.size() < 2 * k + 1 ) ) );
      if ( testChange ) {
        tour = tourCopy;
        position = positionCopy;
      }
    }
  }
            }
          }
        }
      }
    }
    if ( !found ) {
      dontLook[ s1 ] = true;
    }
  }
  return anyChange;
}
bool INLINE_ATTRIBUTE improveTourKOpt_( size_t k,
                                        bool linKernighan,
                                        vector< size_t >& tour,
                                        const vector< vector< double > >& distances,
                                        const vector< vector< size_t > >& nearestNeighbors )
{
  vector< bool > dontLook( tour.size(), false );
  bool change = false;
  bool change4 = true;
  vector< vector< size_t > > nnn( nearestNeighbors );
  while ( change4 ) {
    change4 = false;
    while ( kOptOuterLoop_( k, tour, dontLook, distances, nnn, linKernighan ) ) {
      change = true;
    }

    vector< bool > dontLook4( tour.size(), false );
   if ( 
       twokOptOuterLoop_( 2, tour, dontLook4, dontLook, distances, nnn, linKernighan  ) ) {
//       threekOptOuterLoop_( 2, tour, dontLook4, dontLook, distances, nnn, linKernighan  ) ) {
//        improveTour23InnerLoop_( tour, dontLook4, dontLook, distances, nnn ) ) {
      cerr << "update" << endl;
      change4 = true;
      change = true;
    }
    updateNearest( tour, nnn );
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
                                                             const vector< vector< double > >&,
                                                             const vector< vector< size_t > >& ) > improveTour,
                                            size_t iterations,
                                            vector< size_t >& tour,
                                            const vector< vector< double > >& distances,
                                            const vector< vector< size_t > >& nearestNeighbors )
{
  assert( tour.size() > 8 );
  bool change = false;
  vector< bool > dontLook( tour.size(), false );
  vector< size_t > bestTour( tour );
  double bestLength = getLength_( tour, distances );
  for ( size_t i = 0; i < iterations; ++i ) {
    bool change4 = true;
    while ( change4 ) {
      change4 = false;
      while ( improveTour( tour, dontLook, distances, nearestNeighbors ) ) {
        change = true;
      }

      vector< bool > dontLook4( tour.size(), false );
      while ( improveTour23InnerLoop_( tour, dontLook4, dontLook, distances, nearestNeighbors ) ) {
        change4 = true;
        change = true;
      }
    }
    double length = getLength_( tour, distances );
    if ( length < bestLength ) {
      bestTour = tour;
      bestLength = length;
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

bool INLINE_ATTRIBUTE improveTourIteratedLinKernighan_( size_t iterations,
                                                        vector< size_t >& tour,
                                                        const vector< vector< double > >& distances,
                                                        const vector< vector< size_t > >& nearestNeighbors )
{
  return improveTourIterated_( &linKernighanOuterLoop_, iterations, tour, distances, nearestNeighbors );
}

bool INLINE_ATTRIBUTE improveTourIterated3Opt_( size_t iterations,
                                                vector< size_t >& tour,
                                                const vector< vector< double > >& distances,
                                                const vector< vector< size_t > >& nearestNeighbors )
{
  return improveTourIterated_( &threeOptInnerLoop_, iterations, tour, distances, nearestNeighbors );
}

bool INLINE_ATTRIBUTE improveTourIteratedLKH_( size_t k,
                                               size_t iterations,
                                               bool linKernighan,
                                               vector< size_t >& tour,
                                               const vector< vector< double > >& distances,
                                               const vector< vector< size_t > >& nearestNeighbors )
{
  auto improveTour = [ k, linKernighan ] ( vector< size_t >& tour,
                                           vector< bool >& dontLook,
                                           const vector< vector< double > >& distances,
                                           const vector< vector< size_t > >& nearestNeighbors ) {
    return kOptOuterLoop_( k, tour, dontLook, distances, nearestNeighbors, linKernighan );
  };
  return improveTourIterated_( improveTour, iterations, tour, distances, nearestNeighbors );
}

} // anonymous namespace

vector< size_t > INLINE_ATTRIBUTE TravelingSalespersonProblemSolver::computeTour( const vector< vector< double > >& distances )
{
  assert( !distances.empty() );
  for ( size_t i = 0; i < distances.size(); ++i ) {
    assert( distances.size() == distances[ i ].size() );
    for ( size_t j = 0; j < distances.size(); ++j ) {
      assert( fabs( distances[ i ][ j ] - distances[ j ][ i ] ) < 1e-9 );
      if ( i != j ) {
        assert( distances[ i ][ j ] > 0.0 );
      }
      else {
        if ( distances[ i][i ] != 0.0 ) {
          cerr << "i " << i << " " << distances[ i ][ i ] << endl;
        }
        assert( distances[ i ][ i ] == 0.0 );
      }
    }
  }
  const vector< size_t > tourRand = getRandomTour_( distances );
  const vector< size_t > tourNN = getNearestNeighborTour_( distances );
  const vector< size_t > tourGreedy = getGreedyTour_( distances );
  vector< size_t > bestTour( distances.size(), 0 );
  vector< size_t > tour( tourGreedy );
  cerr << setprecision( 8 );

  vector< vector< size_t > > nearestNeighbors( distances.size() );
  vector< vector< size_t > > nearestNeighbors5( distances.size() );
  vector< vector< size_t > > nearestNeighbors10( distances.size() );
  vector< vector< size_t > > nearestNeighbors30( distances.size() );
  vector< vector< size_t > > helsgaun10( distances.size() );
  vector< vector< size_t > > helsgaun5( distances.size() );
  vector< vector< double > > helsgaunDistances;
  {
    double start( clock() );
    nearestNeighbors30 = computeNearestNeighbors_( distances, 30 );
    helsgaun10 = computeHelsgaunNeighbors_( distances, helsgaunDistances, 10 );
    if ( false ) {
      vector< size_t > helsgainInitialTour = getHelsgaunInitialTour_( nearestNeighbors30, helsgaunDistances, tourGreedy );
    }
    for ( size_t i = 0; i < distances.size(); ++i ) {
      nearestNeighbors5[ i ] = vector< size_t >( nearestNeighbors30[ i ].begin(), nearestNeighbors30[ i ].begin() + min( size_t( 5 ), nearestNeighbors30.size() ) );
      nearestNeighbors10[ i ] = vector< size_t >( nearestNeighbors30[ i ].begin(), nearestNeighbors30[ i ].begin() + min( size_t( 10 ), nearestNeighbors30.size() ) );
      nearestNeighbors[ i ] = vector< size_t >( nearestNeighbors30[ i ].begin(), nearestNeighbors30[ i ].begin() + min( size_t( 20 ), nearestNeighbors30.size() ) );
      helsgaun5[ i ] = vector< size_t >( helsgaun10[ i ].begin(), helsgaun10[ i ].begin() + min( size_t( 5 ), helsgaun10[ i ].size() ) );
    }
    for ( size_t i = 0; i < nearestNeighbors30.size(); ++i ) {
      for ( size_t ji = helsgaun10[ i ].size(); ji > 0; --ji ) {
        size_t j = ji - 1;
        if ( find( nearestNeighbors30[ i ].begin(), nearestNeighbors30[ i ].end(), helsgaun10[ i ][ j ] ) == nearestNeighbors30[ i ].end() ) {
          nearestNeighbors30[ i ].insert( nearestNeighbors30[ i ].begin(), helsgaun10[ i ][ j ] );
        }
      }
      for ( size_t ji = helsgaun5[ i ].size(); ji > 0; --ji ) {
        size_t j = ji - 1;
        if ( find( nearestNeighbors10[ i ].begin(), nearestNeighbors10[ i ].end(), helsgaun5[ i ][ j ] ) == nearestNeighbors10[ i ].end() ) {
          nearestNeighbors10[ i ].insert( nearestNeighbors10[ i ].begin(), helsgaun5[ i ][ j ] );
        }
        if ( find( nearestNeighbors5[ i ].begin(), nearestNeighbors5[ i ].end(), helsgaun5[ i ][ j ] ) == nearestNeighbors5[ i ].end() ) {
          nearestNeighbors5[ i ].insert( nearestNeighbors5[ i ].begin(), helsgaun5[ i ][ j ] );
        }
      }
    }
    /*    for ( size_t i = 0; i < nearestNeighbors5.size(); ++i ) {
          for ( size_t j = 0; j < nearestNeighbors5[ i ].size(); ++j ) {
          cerr << nearestNeighbors5[ i ][ j ] << " ";
          }
          cerr << endl;
          }
          */
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "Time to compute " << nearestNeighbors.front().size() << " nearest neighbors: " << time << endl;
  }
  tour1 = tour2 = tourGreedy;

  cerr << "Greedy distance: " << getLength_( tourGreedy, distances ) << endl;
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
    improveTour23_( tour, distances, nearestNeighbors );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    improveTour23_( tour, distances, nearestNeighbors );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "4-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  if ( false ) {
    tour = tourGreedy;
    double start( clock() );
    improveTourKOpt_( 5, false, tour, distances, nearestNeighbors );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    improveTour23_( tour, distances, nearestNeighbors );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "5-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  if ( true ) {
    tour = tourGreedy;
    double start( clock() );
    improveTourLinKernighan_( tour, distances, nearestNeighbors10 );
    bool threeOpt = improveTour3Opt_( tour, distances, nearestNeighbors30 );
    bool doubleBridge = improveTour23_( tour, distances, nearestNeighbors10 );
    if ( threeOpt || doubleBridge ) {
      improveTourLinKernighan_( tour, distances, nearestNeighbors10 );
      improveTour3Opt_( tour, distances, nearestNeighbors30 );
      improveTour23_( tour, distances, nearestNeighbors10 );
      improveTour3Opt_( tour, distances, nearestNeighbors30 );
    }
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "LK    tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  if ( false ) {
    tour = tourGreedy;
    double start( clock() );
    improveTourIteratedLinKernighan_( 10, tour, distances, nearestNeighbors10 );
    bool threeOpt = improveTour3Opt_( tour, distances, nearestNeighbors30 );
    bool doubleBridge = improveTour23_( tour, distances, nearestNeighbors10 );
    if ( threeOpt || doubleBridge ) {
      improveTourLinKernighan_( tour, distances, nearestNeighbors10 );
      improveTour3Opt_( tour, distances, nearestNeighbors30 );
      improveTour23_( tour, distances, nearestNeighbors10 );
      improveTour3Opt_( tour, distances, nearestNeighbors30 );
    }
    improveTourIterated3Opt_( 1000, tour, distances, nearestNeighbors10 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "ILK   tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  for ( size_t k = 5; k < 6; ++k ) {
    if ( true ) {
      tour = tourGreedy;
      double start( clock() );
      vector< vector< size_t > > nn( nearestNeighbors5 );
      improveTourKOpt_( k, true, tour, distances, nn );

      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << k << "-LK  tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
    }
  }

  for ( size_t k = 5; k < 6; ++k ) {
    if ( true ) {
      tour = tourGreedy;
      double start( clock() );
      improveTourIteratedLKH_( k, 50, true, tour, distances, nearestNeighbors5 );

      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << k << "-ILK tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
    }
  }

  if ( false ) {
    printTour_( "a", tourGreedy );
  }

  return tour;
}
