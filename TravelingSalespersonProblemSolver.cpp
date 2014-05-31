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

namespace {
void printTour_( size_t n, const vector< size_t >& tour )
{
  cerr << n << ". ";
  for ( size_t i = 0; i < tour.size(); ++i ) {
    cerr << tour[ i ] << " ";
  }
  cerr << endl;
}

void INLINE_ATTRIBUTE assertIsTour_( const vector< size_t >& tour,
                                     const vector< size_t >& position )
{
  assert( tour.size() == position.size() );
  for ( size_t i = 0; i < position.size(); ++i ) {
    assert( tour[ position[ i ] ] == i );
  }
}

void INLINE_ATTRIBUTE assertIsTour_( const vector< size_t >& tour,
                                     const vector< vector< double > >& distances )
{
  assert( tour.size() == distances.size() );
  vector< size_t > position( tour.size() );
  for ( size_t i = 0; i < tour.size(); ++i ) {
    position[ tour[ i ] ] = i;
  }
  assertIsTour_( tour, position );
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
  double maxLength = numeric_limits< double >::lowest();
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

bool INLINE_ATTRIBUTE improveTour2Opt_( vector< size_t >& tour,
                                        const vector< vector< double > >& distances,
                                        const vector< vector< size_t > >& nearestNeighbors )
{
  const double eps = 1e-9;
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
        double lengthBefore = getLength_( tour, distances );
        performMove_( bestTs, tour, position );
        assert( getLength_( tour, distances ) < lengthBefore );
        assertIsTour_( tour, position );
        anyChange = true;
        changed = true;
      }
    }
  }
  return anyChange;
}

bool INLINE_ATTRIBUTE improveTour3Opt_( vector< size_t >& tour,
                                        const vector< vector< double > >& distances,
                                        const vector< vector< size_t > >& nearestNeighbors )
{
  const double eps = 1e-9;
  vector< size_t > position( tour.size() );
  for ( size_t i = 0; i < tour.size(); ++i ) {
    position[ tour[ i ] ] = i;
  }

  bool anyChange = false;
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
        double lengthBefore = getLength_( tour, distances );
        performMove_( bestTs, tour, position );
        assert( getLength_( tour, distances ) < lengthBefore );
        assertIsTour_( tour, position );
        changed = true;
        anyChange = true;
      }
    }
  }
  return anyChange;
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

bool INLINE_ATTRIBUTE improveTourDoubleBridge_( vector< size_t >& tour,
                                                const vector< vector< double > >& distances,
                                                const vector< vector< size_t > >& nearestNeighbors )
{
  if ( tour.size() < 5 ) {
    return false;
  }
  const double eps = 1e-9;
  vector< size_t > position( tour.size() );
  for ( size_t i = 0; i < tour.size(); ++i ) {
    position[ tour[ i ] ] = i;
  }
  vector< size_t > bestTs;
  bool anyChange = false;
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
        double lengthBefore = getLength_( tour, distances );
        doubleBridgeSwap_( tour, position, bestTs[ 0 ], bestTs[ 1 ], bestTs[ 2 ], bestTs[ 3 ], bestTs[ 4 ], bestTs[ 5 ], bestTs[ 6 ], bestTs[ 7 ] );
        assert( getLength_( tour, distances ) < lengthBefore );
        assertIsTour_( tour, position );
        anyChange = true;
        changed = true;
      }
    }
  }
  return anyChange;
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
  const double eps = 1e-9;
  bool positiveGain = true;
  size_t numberOfEdgesToRemove = 2;
  while ( positiveGain ) {
    positiveGain = false;
    vector< size_t > mutableNearestNeighborList( nearestNeighbors[ tb ] );
    vector< size_t >& tcUntested = numberOfEdgesToRemove == 2 ? tcUntestedInStep2 : mutableNearestNeighborList;
    assert( tb == next_( t1, tour, position ) || tb == previous_( t1, tour, position ) );
    const bool tBchoice = tb == previous_( t1, tour, position );

    while ( tcUntested.size() > 0 ) {
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
      if ( G + gn <= eps ) {
        continue;
      }

      vector< size_t > bc( { tb, tc } );
      vector< size_t > cd( { tc, td } );
      if ( find( removed.begin(), removed.end(), cd ) != removed.end() ) {
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
      if ( gain > eps ) {
        added.insert( added.end(), { { t1, td }, { td, t1 } } );
        assert( getLength_( tour, distances ) < lengthBefore );
        assertIsTour_( tour, position );
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
  const double eps = 1e-9;
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
        if ( g1 <= eps ) {
          continue;
        }

        // First choice of t4
        {
          size_t t4 = t2choice == 0 ? next_( t3, tour, position ) : previous_( t3, tour, position );
          if ( t4 == previous_( t2, tour, position ) || t4 == next_( t2, tour, position ) ) {
            continue;
          }

          performMove_( { t1, t2, t3, t4 }, tour, position );
          // Test for improving 2-opt move
          double gain = g1 + distances[ t3 ][ t4 ] - distances[ t4 ][ t1 ];
          if ( gain > eps ) {
            assert( getLength_( tour, distances ) < lengthBefore );
            assertIsTour_( tour, position );
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
          while ( tcUntestedInStep2.size() > 0 ) {
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
          if ( t3 == previous_( t2, tour, position ) || t3 == next_( t2, tour, position ) ) {
            continue;
          }
          size_t t4 = t2choice == 0 ? previous_( t3, tour, position ) : next_( t3, tour, position );
          vector< size_t > t5Untested = nearestNeighbors[ t4 ];
          while ( t5Untested.size() > 0 ) {
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
            if ( g1 + g2 <= eps ) {
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
            while ( tcUntestedInStep2.size() > 0 ) {
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

} // anonymous namespace

vector< size_t > INLINE_ATTRIBUTE TravelingSalespersonProblemSolver::computeTour( const vector< vector< double > >& distances )
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
  vector< vector< size_t > > nearestNeighbors5;
  vector< vector< size_t > > nearestNeighbors10;
  vector< vector< size_t > > nearestNeighbors30;
  {
    double start( clock() );
    nearestNeighbors = computeNearestNeighbors_( distances, 20 );
    nearestNeighbors30 = computeNearestNeighbors_( distances, 30 );
    nearestNeighbors10 = computeNearestNeighbors_( distances, 10 );
    nearestNeighbors5 = computeNearestNeighbors_( distances, 5 );
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
    improveTour2Opt_( tour, distances, nearestNeighbors );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "2-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
    assertIsTour_( tour, distances );
  }

  if ( true ) {
    tour = tourGreedy;
    double start( clock() );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "3-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
    assertIsTour_( tour, distances );
  }

  if ( true ) {
    tour = tourGreedy;
    double start( clock() );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    improveTourDoubleBridge_( tour, distances, nearestNeighbors );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "4-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  if ( true ) {
    tour = tourGreedy;
    double start( clock() );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    improveTourDoubleBridge_( tour, distances, nearestNeighbors10 );
    improveTourLinKernighan_( tour, distances, nearestNeighbors10 );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    improveTourDoubleBridge_( tour, distances, nearestNeighbors10 );
    improveTour3Opt_( tour, distances, nearestNeighbors30 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "V-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  if ( false ) {
    printTour_( 1, tourGreedy );
  }

  return tour;
}
