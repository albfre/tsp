/* SYSTEM INCLUDES */
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <ctime>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <set>
#include <random>
#include <functional>

/* HEADER */
#include "TravelingSalespersonProblemSolver.h"
#include "PerformHelsgaunMove.h"
#include "NeighborAlgorithms.h"
#include "InitialTours.h"
#include "HelperFunctions.h"

// For profiling
#if 0
#define INLINE_ATTRIBUTE __attribute__ ((noinline))
#else
#define INLINE_ATTRIBUTE
#endif

using namespace std;
using namespace TravelingSalespersonProblemSolver;

namespace {
  constexpr auto maxNumOfInfeasibleMoves = static_cast< size_t >( 1 );
  constexpr auto useInfeasibleMoves = maxNumOfInfeasibleMoves > 0;
  constexpr auto maxGainMoves = static_cast< size_t >( 1000 );
  constexpr auto takeGreedyMove = true;

  // The basic 2-opt and 3-opt algorithms, which consider all 2/3-opt moves and take the best ones
  bool INLINE_ATTRIBUTE twoOrThreeOptImpl_( const bool only2Opt,
                                            vector< size_t >& tour,
                                            vector< size_t >& position,
                                            vector< bool >& dontLook,
                                            const VDistances& distances,
                                            const vector< vector< size_t > >& nearestNeighbors )
  {
    updatePosition( tour, position );
    auto changed = false;
    vector< size_t > bestTs;

    // For each node, consider all 3-opt moves and take the best one if it improves the tour
    for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
      if ( !dontLook.empty() && dontLook[ t1 ] ) {
        continue;
      }
      auto found = false;
      auto maxGain = 0.0;
      for ( const auto t2choice : { false, true } ) {
        const auto t2 = t2choice ? previous_( t1, tour, position ) : next_( t1, tour, position );
        for ( const auto& t3 : nearestNeighbors[ t2 ] ) {
          if ( isTourNeighbor_( t3, t2, tour, position ) ) {
            continue;
          }
          const auto g1 = distances( t1, t2 ) - distances( t2, t3 );
          if ( g1 <= tolerance ) {
            continue;
          }

          for ( const auto t4choice : only2Opt ? vector< bool >{ !t2choice } : vector< bool >{ false, true } ) {
            const auto t4 = t4choice ? previous_( t3, tour, position ) : next_( t3, tour, position );
            const auto g2 = g1 + distances( t3, t4 );
            const auto equalChoice = t4choice == t2choice;
            if ( !equalChoice ) {
              // Test for improving 2-opt move
              const auto gain = g2 - distances( t4, t1 );
              if ( gain > maxGain ) {
                maxGain = gain;
                bestTs = { t1, t2, t3, t4 };
                if ( takeGreedyMove ) {
                  goto foundGainfulMove;
                }
              }
            }
            if ( only2Opt ) {
              continue;
            }
            for ( const auto& t5 : nearestNeighbors[ t4 ] ) {
              if ( isTourNeighbor_( t5, t4, tour, position ) ) {
                continue;
              }
              const auto g3 = g2 - distances( t4, t5 );
              if ( g3 <= tolerance ||
                   (equalChoice &&
                     ( ( t4choice && !between_( t2, t5, t4, position ) ) ||
                       ( !t4choice && !between_( t2, t4, t5, position ) ) ) ) ) {
                continue;
              }
              for ( const auto t6choice : equalChoice ? vector< bool >{ false, true } : vector< bool >{ false } ) {
                const auto t6 = ( equalChoice && t6choice ) || ( !equalChoice && !between_( t2, t4, t5, position ) )
                  ? previous_( t5, tour, position ) : next_( t5, tour, position );
                if ( t6 == t1 ) {
                  continue;
                }
                const auto g4 = g3 + distances( t5, t6 );
                const auto gain = g4 - distances( t6, t1 );
                if ( gain > maxGain ) {
                  maxGain = gain;
                  bestTs = { t1, t2, t3, t4, t5, t6 };
                  if ( takeGreedyMove ) {
                    goto foundGainfulMove;
                  }
                }
              }
            }
          }
        }
      }
foundGainfulMove:
      if ( maxGain > tolerance ) {
        if ( !dontLook.empty() ) {
          for ( const auto& i : bestTs ) {
            dontLook[ i ] = false;
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

  // The 2-3-opt algorithm, which performs (possibly infeasible) 2-opt moves followed by 2- or 3-opt moves
  // This includes double bridge moves, which are not considered in 5-opt
  bool INLINE_ATTRIBUTE twoThreeOptImpl_( vector< size_t >& tour,
                                          vector< size_t >& position,
                                          vector< bool >& dontLook,
                                          vector< bool >& otherDontLook,
                                          const VDistances& distances,
                                          const vector< vector< size_t > >& nearestNeighbors )
  {
    updatePosition( tour, position );
    auto changed = false;
    vector< size_t > bestTs;
    for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
      if ( !dontLook.empty() && dontLook[ t1 ] ) {
        continue;
      }
      auto found = false;
      auto maxGain = 0.0;
      vector< size_t > incl;
      for ( const auto t2choice : { true } ) {
        const auto t2 = t2choice ? next_( t1, tour, position ) : previous_( t1, tour, position );
        for ( const auto& t3 : nearestNeighbors[ t2 ] ) {
          const auto t4 = t2choice ? previous_( t3, tour, position ) : next_( t3, tour, position );
          if ( t3 == t1 || t4 == t1 ) {
            continue;
          }
          const auto gainFirstBridge = distances( t1, t2 ) + distances( t3, t4 ) - distances( t2, t3 ) - distances( t1, t4 );
          if ( gainFirstBridge <= tolerance ) {
            continue;
          }

          for ( size_t t5 = t2; t5 != t3; t5 = t2choice == 0 ? next_( t5, tour, position ) : previous_( t5, tour, position ) ) {
            for ( const auto t6choice : { false, true } ) {
              auto t6 = t6choice ? next_( t5, tour, position ) : previous_( t5, tour, position );
              if ( t5 == t2 && t6 == t1 ) {
                continue;
              }
              for ( const auto& t7 : nearestNeighbors[ t6 ] ) {
                for ( const auto t8choice : { false, true } ) {
                  const auto t8 = t8choice ? next_( t7, tour, position ) : previous_( t7, tour, position );
                  if ( !between_( t4, t1, t7, position ) || !between_( t4, t1, t8, position ) ) {
                    continue;
                  }
                  const auto gain2 =
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
                    if ( takeGreedyMove ) {
                      goto foundGainfulMove;
                    }
                  }

                  for ( const auto& t9 : nearestNeighbors[ t8 ] ) {
                    const auto t10 = t8choice || between_( t2, t3, t9, position ) ? previous_( t9, tour, position ) : next_( t9, tour, position );
                    if ( ( t9 == t1 && t10 == t2 ) ||
                         ( t9 == t2 && t10 == t1 ) ||
                         ( t9 == t3 && t10 == t4 ) ||
                         ( t9 == t4 && t10 == t3 ) ||
                         ( t9 == t5 && t10 == t6 ) ||
                         ( t9 == t6 && t10 == t5 ) ||
                         ( t9 == t7 && t10 == t8 ) ||
                         ( t9 == t8 && t10 == t7 ) ) {
                      continue;
                    }

                    const auto gain3 =
                      distances( t5, t6 ) + distances( t7, t8 ) + distances( t9, t10 ) - distances( t6, t7 ) - distances( t8, t9 ) - distances( t10, t5 );
                    if ( gainFirstBridge + gain3 > maxGain ) {
                      maxGain = gainFirstBridge + gain3;
                      bestTs = { t1, t2, t3, t4, t5, t6, t7, t8, t9, t10 };
                      incl = { 3, 2, 1, 0, 9, 6, 5, 8, 7, 4 };
                      if ( takeGreedyMove ) {
                        goto foundGainfulMove;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
foundGainfulMove:
      if ( maxGain > tolerance ) {
        if ( !dontLook.empty() ) {
          for ( const auto& bi : bestTs ) {
            dontLook[ bi ] = false;
          }
        }
        if ( !otherDontLook.empty() ) {
          for ( const auto& bi : bestTs ) {
            otherDontLook[ bi ] = false;
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

  // Find the best feasible and infeasible k-opt moves
  double INLINE_ATTRIBUTE getBestKOptMove_( const size_t currentDepth,
                                            const size_t k,
                                            const double G0,
                                            vector< size_t >& ts,
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
    const auto t1 = ts.front();
    const auto t2 = ts.back();
    for ( const auto& t3 : nearestNeighbors[ t2 ] ) {
      if ( isTourNeighbor_( t3, t2, tour, position ) ) {
        continue;
      }
      const auto G1 = G0 - distances( t2, t3 );
      if ( G1 <= tolerance ) {
        continue;
      }
      if ( find( added.begin(), added.end(), make_pair( t2, t3 ) ) != added.end() ) {
        continue;
      }
      for ( const auto& t4 : { previous_( t3, tour, position ), next_( t3, tour, position ) } ) {
        const auto t3t4Pair = make_pair( t3, t4 );
        if ( currentDepth + 1 == k ) {
          const auto G2 = G1 + distances( t3, t4 );
          if ( G2 - distances( t1, t4 ) <= tolerance && ( G2 <= bestG && G2 <= bestInfeasibleG ) ) {
            continue;
          }
          if ( find( added.begin(), added.end(), t3t4Pair ) != added.end() ) {
            continue;
          }
        }
        if ( find( removed.begin(), removed.end(), t3t4Pair ) != removed.end() ) {
          continue;
        }

        ts.push_back( t3 );
        ts.push_back( t4 );
        const auto G1 = G0 - distances( t2, t3 );
        const auto G2 = G1 + distances( t3, t4 );
        const auto gain1 = G2 - distances( t1, t4 );
        if ( gain1 > tolerance && makesTour_( ts, tour, position ) ) {
          return gain1;
        }

        if ( currentDepth + 1 < k ) {
          added.push_back( { t2, t3 } );
          added.push_back( { t3, t2 } );
          removed.push_back( { t3, t4 } );
          removed.push_back( { t4, t3 } );
          const auto gain2 = getBestKOptMove_( currentDepth + 1,
                                               k,
                                               G2,
                                               ts,
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
          if ( gain2 > tolerance ) {
            return gain2;
          }
          removed.resize( removed.size() - 2 );
          added.resize( added.size() - 2 );
        }
        else if ( G2 > bestG && gain1 <= tolerance && makesTour_( ts, tour, position ) ) {
          // if gain1 > tolerance && makesTour, this conditional will not be reached
          bestG = G2;
          bestTs = ts;
        }
        else if ( G2 > bestInfeasibleG ) {
          bestInfeasibleG = G2;
          bestInfeasibleTs = ts;
        }

        ts.resize( ts.size() - 2 );
      }
    }
    return -1.0;
  }

  // The Lin-Kernighan-Helsgaun algorithm, which strings gainful sequences of (possibly non-gainful) k-opt moves together
  bool INLINE_ATTRIBUTE kOptImpl_( size_t k,
                                   vector< size_t >& tour,
                                   vector< size_t >& position,
                                   vector< bool >& dontLook,
                                   const vector< size_t >& betterTour,
                                   const VDistances& distances,
                                   const vector< vector< size_t > >& nearestNeighbors,
                                   bool linKernighan )
  {
    updatePosition( tour, position );
    auto tourChanged = false;
    vector< size_t > ts;
    vector< pair< size_t, size_t > > added;
    vector< pair< size_t, size_t > > removed;
    vector< size_t > bestTs;
    vector< size_t > bestInfeasibleTs;
    vector< size_t > tourCopy;
    vector< size_t > tsHistory;
    for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
      if ( dontLook[ t1 ] ) {
        continue;
      }
      auto found = false;
      for ( const auto t2choice : { true, false } ) {
        const auto t2 = t2choice ? previous_( t1, tour, position ) : next_( t1, tour, position );
        if ( inBetterTour( t1, t2, betterTour ) ) {
          continue;
        }

        ts = { t1, t2 };
        added.clear();
        removed = { { t1, t2 }, { t2, t1 } };
        bestTs.clear();
        bestInfeasibleTs.clear();
        tourCopy = tour;
        tsHistory.clear();
        auto G0 = distances( t1, t2 );
        auto bestG = tolerance;
        auto bestInfeasibleG = tolerance;
        auto testChange = false;
        size_t lkDepth = 0;
        do {
          ++lkDepth;
          bestG = tolerance;
          bestInfeasibleG = useInfeasibleMoves ? tolerance : numeric_limits< double >::max();
          const size_t currentDepth = 1;
          auto gain = getBestKOptMove_( currentDepth,
                                        k,
                                        G0,
                                        ts,
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
            // Improving move found, make it!
            performKOptMove_( ts, tour, position );
            assert( isTour_( tour, position ) );
            tsHistory.insert( tsHistory.end(), ts.begin(), ts.end() );
            for ( const auto& th : tsHistory ) {
              dontLook[ th ] = false;
            }
            found = true;
            tourChanged = true;
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
            // Feasible move with positive g found, test it!
            performKOptMove_( bestTs, tour, position );
            testChange = true;
            G0 = bestG;
            ts = { bestTs.front(), bestTs.back() };
            tsHistory.insert( tsHistory.end(), bestTs.begin(), bestTs.end() );
            removed = { { bestTs.front(), bestTs.back() }, { bestTs.back(), bestTs.front() } };
            added.clear();
            for ( size_t i = 1; i + 1 < tsHistory.size(); i += 2 ) {
              added.push_back( { tsHistory[ i ], tsHistory[ i + 1 ] } );
              added.push_back( { tsHistory[ i + 1 ], tsHistory[ i ] } );
            }
          }
          else if ( bestInfeasibleG > tolerance ) {
            // Infeasible move with positive g found, explore it further!
            G0 = bestInfeasibleG;
            ts = bestInfeasibleTs;
            removed.clear();
            added.clear();
            for ( size_t i = 0; i + 1 < ts.size(); i += 2 ) {
              removed.push_back( { ts[ i ], ts[ i + 1 ] } );
              removed.push_back( { ts[ i + 1 ], ts[ i ] } );
            }
            for ( size_t i = 1; i + 1 < ts.size(); i += 2 ) {
              added.push_back( { ts[ i ], ts[ i + 1 ] } );
              added.push_back( { ts[ i + 1 ], ts[ i ] } );
            }
          }
        } while ( ( bestG > tolerance && lkDepth < maxGainMoves ) || ( bestInfeasibleG > tolerance && ts.size() < maxNumOfInfeasibleMoves * 2 * k + 1 ) );
        if ( testChange ) {
          swap( tour, tourCopy );
          updatePosition( tour, position );
        }
      }
      if ( !found ) {
        dontLook[ t1 ] = true;
      }
    }
    return tourChanged;
  }

  // Perturb the current tour by swapping edges
  void kSwapKick_( size_t numOfPairs,
                   vector< size_t >& tour,
                   vector< size_t >& position,
                   vector< bool >& dontLook,
                   std::function<int()> random )
  {
    assert( numOfPairs > 0 );
    updatePosition( tour, position );

    vector< size_t > swapEdges;
    swapEdges.reserve( 2 * numOfPairs );
    while ( swapEdges.size() < 2 * numOfPairs && swapEdges.size() < 2 * ( tour.size() / 4 ) ) {
      size_t i = random() % tour.size();
      if ( find( swapEdges.begin(), swapEdges.end(), i ) == swapEdges.end() ) {
        swapEdges.push_back( i );
      }
    }
    sort( swapEdges.begin(), swapEdges.end(), [&] ( const auto& ta, const auto& tb ) { return position[ ta ] < position[ tb ]; } );
    vector< size_t > ts;
    ts.reserve( 2 * swapEdges.size() );
    for ( size_t i = 0; i < swapEdges.size(); ++i ) {
      ts.push_back( swapEdges[ i ] );
      ts.push_back( next_( swapEdges[ i ], tour, position ) );
    }
    const auto tourCopy = tour;
    size_t index = 0;
    for ( size_t i = 0; i + 1 < ts.size(); i += 2 ) {
      for ( size_t t = ts[ i ]; t != ( i == 0 ? ts.back() : ts[ i - 1 ] ); t = previous_( t, tourCopy, position ), ++index ) {
        tour[ index ] = t;
      }
      tour[ index ] = i == 0 ? ts.back() : ts[ i - 1 ];
      ++index;
    }
    assert( index == tour.size() );

    for ( const auto& i : ts ) {
      dontLook[ i ] = false;
    }
  }

  // Run run the k-opt algorithm multiple times with pertubations inbetween iterations
  bool INLINE_ATTRIBUTE improveTourIteratedKOpt_( const size_t k,
                                                  const bool linKernighan,
                                                  const size_t iterations,
                                                  const bool useGain23,
                                                  vector< size_t >& tour,
                                                  vector< size_t >& betterTour,
                                                  const VDistances& distances,
                                                  const vector< vector< size_t > >& nearestNeighbors,
                                                  std::function<int()> random )
  {
    assert( k >= 2 );
    const auto improveTour = [ k, linKernighan, &betterTour ] ( vector< size_t >& tour,
                                                                vector< size_t >& position,
                                                                vector< bool >& dontLook,
                                                                const VDistances& distances,
                                                                const vector< vector< size_t > >& nearestNeighbors ) {
      if ( linKernighan ) {
        return kOptImpl_( k, tour, position, dontLook, betterTour, distances, nearestNeighbors, linKernighan );
      }
      switch ( k ) {
        case 2: return twoOrThreeOptImpl_( true, tour, position, dontLook, distances, nearestNeighbors );
        case 3: return twoOrThreeOptImpl_( false, tour, position, dontLook, distances, nearestNeighbors );
        default: return kOptImpl_( k, tour, position, dontLook, betterTour, distances, nearestNeighbors, linKernighan );
      }
    };

    auto change = false;
    vector< bool > dontLook( tour.size(), false );
    vector< bool > dontLook4( tour.size(), false );
    vector< size_t > position( tour.size() );
    auto bestTour = tour;
    auto bestLength = getLength_( tour, distances );
    auto tour1 = tour;
    auto tour2 = tour;
    auto nn = nearestNeighbors;
    for ( size_t i = 0; i < iterations; ++i ) {
      auto change4 = true;
      while ( change4 ) {
        change4 = false;
        while ( improveTour( tour, position, dontLook, distances, nn ) ) {
          change = true;
        }

        if ( useGain23 ) {
          fill( dontLook4.begin(), dontLook4.end(), false );
          while ( twoThreeOptImpl_( tour, position, dontLook4, dontLook, distances, nn ) ) {
            change4 = true;
            change = true;
          }
        }
      }

      const auto length = getLength_( tour, distances );
      if ( length < bestLength ) {
        bestTour = tour;
        bestLength = length;
        betterTour = bestTour;
        updateNearest( tour, tour1, tour2, nn );
      }

      tour = bestTour;
      if ( i + 1 < iterations ) {
        kSwapKick_( min( tour.size() / 4, static_cast< size_t >( ( random() ) % 16 ) ) + 1, tour, position, dontLook, random );
      }
    }
    return change;
  }

  // Run the k-opt algorithm
  bool INLINE_ATTRIBUTE improveTourKOpt_( const size_t k,
                                          const bool linKernighan,
                                          vector< size_t >& tour,
                                          vector< size_t >& betterTour,
                                          const VDistances& distances,
                                          const vector< vector< size_t > >& nearestNeighbors,
                                          std::function<int()> random )
  {
    const size_t numberOfIterations = 1;
    const auto useGain23 = k > 3;
    return improveTourIteratedKOpt_( k,
                                     linKernighan,
                                     numberOfIterations,
                                     useGain23,
                                     tour,
                                     betterTour,
                                     distances,
                                     nearestNeighbors,
                                     random );
  }


  inline double distanceRound( double d ) { return double( long( d + 0.5 ) ); }
} // anonymous namespace

vector< size_t > INLINE_ATTRIBUTE TravelingSalespersonProblemSolver::computeTour( const VDistances& distances )
{
  {
    double start( clock() );
    cerr << "Starting sanity check ... ";
    assert( !distances.empty() );
    for ( size_t i = 0; i < distances.size(); ++i ) {
      assert( distances( i, i ) == 0.0 );
      for ( size_t j = i + 1; j < distances.size(); ++j ) {
        //std::cout << distances(i,j) << " ";
        assert( fabs( distances( i, j ) - distances( j, i ) ) < 1e-9 );
        assert( distances( i, j ) >= 0.0 );
        assert( distances( j, i ) >= 0.0 );
      }
      //std::cout << std::endl;
    }
    double timeSanity( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "done: " << timeSanity << endl;

    for ( size_t i = 0; i < distances.size(); ++i ) {
      for ( size_t j = 0; j < distances.size(); ++j ) {
       // std::cout << distances(i,j) << " ";
      }
      //std::cout << std::endl;
    }
  }

  std::mt19937 generator( 1338 );
  std::uniform_int_distribution<> uniformDist( 0, distances.size() - 1 );
  auto random = std::bind( uniformDist, std::ref( generator ) );

  double start( clock() );
  const vector< size_t > tourRand = getRandomTour_( distances, random );
  double timeRand( ( clock() - start ) / CLOCKS_PER_SEC );
  start = clock();
  const vector< size_t > tourNN = getNearestNeighborTour_( distances, random );
  double timeNN( ( clock() - start ) / CLOCKS_PER_SEC );
  cerr << setprecision( 8 );


  vector< vector< size_t > > nearestNeighbors5( distances.size() );
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
      nearestNeighbors5[ i ] = vector< size_t >( nearestNeighbors30[ i ].begin(), nearestNeighbors30[ i ].begin() + min( static_cast< size_t >( 5 ), nearestNeighbors30[ i ].size() ) );
      helsgaun5[ i ] = vector< size_t >( helsgaun10[ i ].begin(), helsgaun10[ i ].begin() + min( static_cast< size_t >( 5 ), helsgaun10[ i ].size() ) );
    }
    auto addVector = [&] ( vector< vector< size_t > >& nn, const vector< vector< size_t > >& hn ) {
      assert( nn.size() == hn.size() );
      for ( size_t i = 0; i < nn.size(); ++i ) {
        nn[ i ].insert( nn[ i ].begin(), hn[ i ].begin(), hn[ i ].end() );

        nn[ i ].erase( remove_if( nn[ i ].begin() + hn[ i ].size(), nn[ i ].end(),
                                  [&] ( size_t index ) { return find( hn[ i ].begin(), hn[ i ].end(), index ) != hn[ i ].end(); } ),
                       nn[ i ].end() );
      }
    };
    addVector( nearestNeighbors30, helsgaun10 );
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

  for ( const auto& k : { 2, 3, 5 } ) {
    tour = tourGreedy;
    vector< size_t > betterTour;
    bool linKernighan = false;
    double start( clock() );
    improveTourKOpt_( k, linKernighan, tour, betterTour, distances, nearestNeighbors30, random );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << k << "-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
    assert( isTour_( tour, distances ) );
  }

  for ( const auto& k : { 3 } ) {
    tour = tourGreedy;
    vector< size_t > betterTour;
    bool linKernighan = false;
    double start( clock() );
    improveTourIteratedKOpt_( k, linKernighan, 1000, false, tour, betterTour, distances, nearestNeighbors30, random );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "I" << k << "    tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  for ( const auto& k : { 5 } ) {
    tour = tourGreedy;
    vector< size_t > betterTour;
    vector< size_t > tour1( tour );
    vector< size_t > tour2( tour );
    vector< vector< size_t > > nn( nearestNeighbors5 );
    vector< size_t > bestTour( tour );
    vector< vector< double > > points;
    double start( clock() );
    for ( size_t i = 0; i < 10; ++i ) {
      tour = getHelsgaunInitialTour_( nn, helsgaun10, helsgaunDistances10, bestTour, random );
      bool useGain23 = true;
      while ( improveTourKOpt_( k, useGain23, tour, betterTour, distances, nn, random ) ) {
        updateNearest( tour, tour1, tour2, nn );
      }
      if ( getLength_( tour, distances ) < getLength_( bestTour, distances ) ) {
        bestTour = tour;
        betterTour = bestTour;
      }
    }
    tour = bestTour;

    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << k << "-LK  tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  for ( const auto& k : { 5 } ) {
    tour = tourGreedy;
    vector< size_t > betterTour;
    bool linKernighan = true;
    bool useGain23 = false;
    double start( clock() );
    improveTourIteratedKOpt_( k, linKernighan, 100, useGain23, tour, betterTour, distances, nearestNeighbors5, random );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << k << "-ILK tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
    /*
    for ( const auto& t : tour ) {
      cerr << t << " -> ";
    }
    cerr << std::endl;
    */
  }

  if ( false ) {
    printTour_( "a", tourGreedy );
  }

  return tour;
}
