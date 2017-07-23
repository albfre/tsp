#ifndef NEIGHBOR_ALGORITHMS_H
#define NEIGHBOR_ALGORITHMS_H

namespace TravelingSalespersonProblemSolver {

  // Updates the nearest neighbor matrix on the basis of two previous tours
  void updateNearest( const std::vector< size_t >& tour,
                      std::vector< size_t >& tour1,
                      std::vector< size_t >& tour2,
                      std::vector< std::vector< size_t > >& nearest )
  {
    tour2 = tour1;
    tour1 = tour;

    std::set< std::vector< size_t > > e1;
    std::set< std::vector< size_t > > e2;
    for ( size_t i = 0; i < tour1.size(); ++i ) {
      std::vector< size_t > e( { tour1[ i ], tour1[ i + 1 == tour1.size() ? 0 : i + 1 ] } );
      std::sort( e.begin(), e.end() );
      e1.insert( e );

      e = { tour2[ i ], tour2[ i + 1 == tour2.size() ? 0 : i + 1 ] };
      std::sort( e.begin(), e.end() );
      e2.insert( e );
    }

    std::vector< std::vector< size_t > > v( tour.size() );
    v.erase( std::set_intersection( e1.begin(), e1.end(), e2.begin(), e2.end(), v.begin() ), v.end() );

    for ( const auto& vi : v ) {
      const size_t i1 = vi[ 0 ];
      const size_t i2 = vi[ 1 ];
      nearest[ i1 ].insert( nearest[ i1 ].begin(), i2 );
      nearest[ i1 ].erase( remove_if( nearest[ i1 ].begin() + 1, nearest[ i1 ].end(), [&] ( size_t index ) { return index == i2; } ), nearest[ i1 ].end() );

      nearest[ i2 ].insert( nearest[ i2 ].begin(), i1 );
      nearest[ i2 ].erase( remove_if( nearest[ i2 ].begin() + 1, nearest[ i2 ].end(), [&] ( size_t index ) { return index == i1; } ), nearest[ i2 ].end() );
    }
  }

  std::vector< std::vector< size_t > > computeNearestNeighbors_( const VDistances& distances,
                                                                 size_t numberOfNeighbors )
  {
    assert( !distances.empty() );
    numberOfNeighbors = std::min( numberOfNeighbors, distances.size() - 1 );
    assert( numberOfNeighbors > 0 );
    std::vector< std::vector< size_t > > nearestNeighbors( distances.size(), std::vector< size_t >( numberOfNeighbors ) );
    std::set< std::pair< double, size_t > > tmpNeighbors;
    for ( size_t i = 0; i < distances.size(); ++i ) {
      tmpNeighbors.clear();
      double worstNearNeighbor = std::numeric_limits< double >::max();
      for ( size_t j = 0; j < distances.size(); ++j ) {
        if ( j != i ) {
          if ( tmpNeighbors.size() < numberOfNeighbors || distances( i, j ) < worstNearNeighbor ) {
            if ( tmpNeighbors.size() >= numberOfNeighbors ) {
              auto itLast = tmpNeighbors.end();
              --itLast;
              tmpNeighbors.erase( itLast );
            }
            tmpNeighbors.insert( std::make_pair( distances( i, j ), j ) );
            worstNearNeighbor = tmpNeighbors.rbegin()->first;
          }
        }
      }
      assert( tmpNeighbors.size() == numberOfNeighbors );
      auto it = tmpNeighbors.begin();
      for ( size_t j = 0; j < numberOfNeighbors; ++j, ++it ) {
        nearestNeighbors[ i ][ j ] = it->second;
      }
    }
    return nearestNeighbors;
  }

  struct Vertex {
    Vertex( size_t idx ) : nodeIndex( idx ), cost( std::numeric_limits< double >::max() ) {}
    size_t nodeIndex;
    size_t parentIndexInVertexList;
    double cost;
  };

  double compute1Tree_( std::vector< Vertex >& vertices,
                        std::vector< std::pair< size_t, size_t > >& edges,
                        std::vector< size_t >& vertexDegrees,
                        const VDistances& distances,
                        const std::vector< double >& lagrangeMultipliers )
  {
    assert( distances.size() > 1 );
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
    double minimaxValue = std::numeric_limits< double >::max();
    for ( size_t i = 0; i + 1 < distances.size(); ++i ) {
      double value = 0.0;
      for ( size_t j = i + 1; j < distances.size(); ++j ) {
        value = std::max( value, distances( i, j ) );
      }
      if ( value < minimaxValue  ) {
        minimaxValue = value;
        minimaxVertex = i;
      }
    }

    vertices.push_back( Vertex( minimaxVertex ) );
    const size_t rootVertex = ( minimaxVertex + 1 ) % distances.size();
    vertices.push_back( Vertex( rootVertex ) );

    std::vector< Vertex > unusedVertices;
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
      Vertex& closestUnusedVertex = *min_element( unusedVertices.begin(), unusedVertices.end(), [] ( const Vertex& v1, const Vertex& v2 ) { return v1.cost < v2.cost; } );
      treeLength += closestUnusedVertex.cost;
      const size_t indexOfNewTreeNode = closestUnusedVertex.nodeIndex;
      const size_t indexOfClosestTreeNode = vertices[ closestUnusedVertex.parentIndexInVertexList ].nodeIndex;
      ++vertexDegrees[ indexOfNewTreeNode ];
      ++vertexDegrees[ indexOfClosestTreeNode ];
      vertices.push_back( closestUnusedVertex );
      edges.push_back( std::make_pair( indexOfClosestTreeNode, indexOfNewTreeNode ) );
      closestUnusedVertex = unusedVertices.back();
      unusedVertices.pop_back();

      for ( Vertex& v : unusedVertices ) {
        const double newDistance = getDistance_( distances, lagrangeMultipliers, indexOfNewTreeNode, v.nodeIndex );
        if ( newDistance < v.cost ) {
          v.cost = newDistance;
          v.parentIndexInVertexList = vertices.size() - 1;
        }
      }
    }

    // 2. Add the two shortest edges connecting to the first vertex
    double minElement = std::numeric_limits< double >::max();
    double secondMinElement = std::numeric_limits< double >::max();
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
    edges.push_back( std::make_pair( minimaxVertex, secondMinElementIndex ) );
    edges.push_back( std::make_pair( minimaxVertex, minElementIndex ) );

    treeLength += minElement + secondMinElement;
    return treeLength - 2.0 * accumulate( lagrangeMultipliers.begin(), lagrangeMultipliers.end(), 0.0 );
  }

  double getHeldKarpLowerBound_( const VDistances& distances, std::vector< double >& lagrangeMultipliers )
  {
    double maxLength = std::numeric_limits< double >::lowest();
    lagrangeMultipliers.assign( distances.size(), 0.0 );
    std::vector< double > localLagrangeMultipliers( lagrangeMultipliers );
    std::vector< Vertex > vertices;
    std::vector< size_t > vertexDegrees;
    std::vector< std::pair< size_t, size_t > > edges;
    size_t maxIter = 50;
    if ( distances.size() > 50000 ) {
      maxIter = 1;
    }
    double lambda = 0.1;
    for ( size_t i = 0; i < maxIter; ++i ) {
      std::vector< size_t > vertexDegrees;
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

  double getHeldKarpLowerBound_( const VDistances& distances )
  {
    std::vector< double > lagrangeMultipliers( distances.size() );
    return getHeldKarpLowerBound_( distances, lagrangeMultipliers );
  }

  std::vector< std::vector< size_t > > computeHelsgaunNeighbors_( const VDistances& distances,
                                                                  std::vector< std::vector< double > >& alphaDistances,
                                                                  size_t numberOfNeighbors )
  {
    assert( distances.size() > 1 );
    numberOfNeighbors = std::min( numberOfNeighbors, distances.size() - 1 );
    assert( numberOfNeighbors > 0 );

    std::vector< Vertex > vertices;
    std::vector< size_t > vertexDegrees;
    std::vector< std::pair< size_t, size_t > > edges;
    std::vector< double > lagrangeMultipliers( distances.size(), 0.0 );
    getHeldKarpLowerBound_( distances, lagrangeMultipliers );

    compute1Tree_( vertices,
                   edges,
                   vertexDegrees,
                   distances,
                   lagrangeMultipliers );

    std::vector< std::vector< size_t > > nearestNeighbors( distances.size(), std::vector< size_t >( numberOfNeighbors ) );
    alphaDistances.resize( distances.size(), std::vector< double >( numberOfNeighbors ) );
    std::set< std::pair< double, size_t > > tmpNeighbors;
    std::vector< double > a( distances.size(), 0.0 );
    std::vector< double > b( distances.size() );
    std::vector< size_t > mark( distances.size(), 0 );
    size_t frontInd = vertices.front().nodeIndex;
    const double maxEdgeWeight = std::max( getDistance_( distances, lagrangeMultipliers, edges.end()[ -2 ].first, edges.end()[ -2 ].second ),
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
        b[ i ] = std::numeric_limits< double >::lowest();
        size_t j = 0;
        for ( size_t k = i; k != 1; k = j ) {
          j = vertices[ k ].parentIndexInVertexList;
          b[ j ] = std::max( b[ k ], getDistance_( distances, lagrangeMultipliers, vertices[ j ].nodeIndex, vertices[ k ].nodeIndex ) );
          mark[ j ] = i;
        }
        fill( a.begin(), a.end(), 0.0 );
        for ( j = 1; j < vertices.size(); ++j ) {
          if ( j != i ) {
            size_t jInd = vertices[ j ].nodeIndex;
            if ( mark[ j ] != i ) {
              size_t jParent = vertices[ j ].parentIndexInVertexList;
              b[ j ] = std::max( b[ jParent ], getDistance_( distances, lagrangeMultipliers, jInd, vertices[ jParent ].nodeIndex ) );
            }
            a[ jInd ] = getDistance_( distances, lagrangeMultipliers, iInd, jInd ) - b[ j ];
          }
        }
        a[ frontInd ] = getDistance_( distances, lagrangeMultipliers, iInd, frontInd ) - maxEdgeWeight;
      }
      tmpNeighbors.clear();
      double worstNearNeighbor = std::numeric_limits< double >::max();
      for ( size_t j = 0; j < distances.size(); ++j ) {
        if ( j != i ) {
          size_t jInd = vertices[ j ].nodeIndex;
          if ( tmpNeighbors.size() < numberOfNeighbors || a[ jInd ] < worstNearNeighbor ) {
            if ( tmpNeighbors.size() >= numberOfNeighbors ) {
              auto itLast = tmpNeighbors.end();
              --itLast;
              tmpNeighbors.erase( itLast );
            }
            tmpNeighbors.insert( std::make_pair( a[ jInd ], jInd ) );
            worstNearNeighbor = tmpNeighbors.rbegin()->first;
          }
        }
      }
      assert( tmpNeighbors.size() == numberOfNeighbors );
      auto it = tmpNeighbors.cbegin();
      for ( size_t j = 0; j < numberOfNeighbors; ++j, ++it ) {
        alphaDistances[ iInd ][ j ] = it->first;
        nearestNeighbors[ iInd ][ j ] = it->second;
      }
    }
    return nearestNeighbors;
  }

}

#endif // NEIGHBOR_ALGORITHMS_H
