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
  const auto useInfeasibleMoves = true;
  const auto maxGainMoves = static_cast< size_t >( 1000 );

  bool INLINE_ATTRIBUTE twoOptInnerLoop_( vector< size_t >& tour,
                                          vector< size_t >& position,
                                          vector< bool >& dontLook,
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
      for ( const auto t2choice : { true, false } ) {
        const auto t2 = t2choice ? previous_( t1, tour, position ) : next_( t1, tour, position );
        for ( const auto& t3 : nearestNeighbors[ t2 ] ) {
          if ( t3 == previous_( t2, tour, position ) || t3 == next_( t2, tour, position ) ) {
            continue;
          }
          if ( distances( t2, t3 ) >= distances( t1, t2 ) ) {
            continue;
          }
          const auto t4 = t2choice ? next_( t3, tour, position ) : previous_( t3, tour, position );
          const auto gain = distances( t1, t2 ) + distances( t3, t4 ) - distances( t1, t4 ) - distances( t2, t3 );
          if ( gain > maxGain ) {
            maxGain = gain;
            bestTs = { t1, t2, t3, t4 };
          }
        }
      }
      if ( maxGain > tolerance ) {
        if ( !dontLook.empty() ) {
          for ( const auto& i : bestTs ) {
            dontLook[ i ]= false;
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

  bool INLINE_ATTRIBUTE threeOptInnerLoop_( vector< size_t >& tour,
                                            vector< size_t >& position,
                                            vector< bool >& dontLook,
                                            const VDistances& distances,
                                            const vector< vector< size_t > >& nearestNeighbors )
  {
    auto changed = false;
    updatePosition( tour, position );
    vector< size_t > bestTs;
    for ( size_t t1 = 0; t1 < tour.size(); ++t1 ) {
      if ( !dontLook.empty() && dontLook[ t1 ] ) {
        continue;
      }
      auto found = false;
      auto G = 0.0;
      for ( const auto t2choice : { true, false } ) {
        const auto t2 = t2choice ? previous_( t1, tour, position ) : next_( t1, tour, position );
        for ( const auto& t3 : nearestNeighbors[ t2 ] ) {
          if ( t3 == previous_( t2, tour, position ) || t3 == next_( t2, tour, position ) ) {
            continue;
          }
          const auto g1 = distances( t1, t2 ) - distances( t2, t3 );
          if ( g1 <= tolerance ) {
            continue;
          }
          // First choice of t4
          auto t4 = t2choice ? next_( t3, tour, position ) : previous_( t3, tour, position );
          if ( t4 == previous_( t2, tour, position ) || t4 == next_( t2, tour, position ) ) {
            continue;
          }
          {
            // Test for improving 2-opt move
            const auto gain = g1 + distances( t3, t4 ) - distances( t4, t1 );
            if ( gain > G ) {
              G = gain;
              bestTs = { t1, t2, t3, t4 };
            }
          }
          for ( const auto& t5 : nearestNeighbors[ t4 ] ) {
            if ( t5 == previous_( t4, tour, position ) || t5 == next_( t4, tour, position ) ) {
              continue;
            }
            const auto g2 = distances( t3, t4 ) - distances( t4, t5 );
            if ( g1 + g2 <= tolerance ) {
              continue;
            }

            // Select t6 such that a valid tour is created
            const auto t6 = between_( t2, t4, t5, position ) ? next_( t5, tour, position ) : previous_( t5, tour, position );
            if ( t6 == t1 ) {
              continue;
            }
            const auto g3 = distances( t5, t6 ) - distances( t6, t1 );
            const auto gain = g1 + g2 + g3;
            if ( gain > G ) {
              G = gain;
              bestTs = { t1, t2, t3, t4, t5, t6 };
            }
          }

          // Second choice of t4
          t4 = t2choice ? previous_( t3, tour, position ) : next_( t3, tour, position );
          for ( const auto& t5 : nearestNeighbors[ t4 ] ) {
            if ( t5 == previous_( t4, tour, position ) || t5 == next_( t4, tour, position ) ) {
              continue;
            }
            if ( ( t2choice && !between_( t3, t2, t5, position ) ) ||
                 ( !t2choice && !between_( t2, t3, t5, position ) ) ) {
              continue;
            }
            const auto g2 = distances( t3, t4 ) - distances( t4, t5 );
            if ( g1 + g2 <= tolerance ) {
              continue;
            }
            // Only consider one choice of t6. The other choice is possible, but clutters the code and doesn't lead to a significant improvement.
            const auto t6choice = t2choice;
            const auto t6 = t6choice ? next_( t5, tour, position ) : previous_( t5, tour, position );
            if ( t6 == t3 || t6 == t2 || t6 == t1 ) {
              continue;
            }
            const auto g3 = distances( t5, t6 ) - distances( t6, t1 );
            const auto gain = g1 + g2 + g3;
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

  bool INLINE_ATTRIBUTE twoThreeOptInnerLoop_( vector< size_t >& tour,
                                               vector< size_t >& position,
                                               vector< bool >& dontLook,
                                               vector< bool >& otherDontLook,
                                               const VDistances& distances,
                                               const vector< vector< size_t > >& nearestNeighbors )
  {
    // Performs infeasible 2-opt moves followed by a 2- or 3-opt move.
    // Includes double bridge moves.
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
          const auto t4 = t2choice == 0 ? next_( t3, tour, position ) : previous_( t3, tour, position );
          if ( t3 == t1 || t4 == t1 ) {
            continue;
          }
          const auto gainFirstBridge = distances( t1, t2 ) + distances( t3, t4 ) - distances( t2, t3 ) - distances( t1, t4 );
          if ( gainFirstBridge <= tolerance ) {
            continue;
          }

          for ( size_t t5 = t2; t5 != t3; t5 = t2choice == 0 ? next_( t5, tour, position ) : previous_( t5, tour, position ) ) {
            for ( const auto t6choice : { true, false } ) {
              auto t6 = t6choice ? next_( t5, tour, position ) : previous_( t5, tour, position );
              if ( t5 == t2 && t6 == t1 ) {
                continue;
              }
              for ( const auto& t7 : nearestNeighbors[ t6 ] ) {
                for ( const auto t8choice : { true, false  } ) {
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
                  }

                  for ( const auto& t9 : nearestNeighbors[ t8 ] ) {
                    const auto t10 = t8choice || between_( t2, t3, t9, position ) ? previous_( t9, tour, position ) : next_( t9, tour, position );
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

                    const auto gain3 =
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

  double INLINE_ATTRIBUTE getBestKOptMove_( size_t depth,
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
    const auto tb = ts.back();
    vector< pair< size_t, size_t > > tcTdPairs;
    for ( const auto& tc : nearestNeighbors[ tb ] ) {
      const auto G1 = G0 - distances( tb, tc );
      if ( G1 <= tolerance ) {
        continue;
      }
      if ( tc == next_( tb, tour, position ) ||
           tc == previous_( tb, tour, position ) ) {
        // The added edge should not belong to T
        continue;
      }
      for ( const auto tdChoice : { true, false } ) {
        const auto td = tdChoice ? previous_( tc, tour, position ) : next_( tc, tour, position );
        if ( depth + 1 == k ) {
          const auto G2 = G1 + distances( tc, td );
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
      const auto tc1 = p1.first;
      const auto td1 = p1.second;
      const auto tc2 = p2.first;
      const auto td2 = p2.second;
      return distances( tc1, td1 ) - distances( tb, tc1 )
             < distances( tc2, td2 ) - distances( tb, tc2 );
    } );
    */
    reverse( tcTdPairs.begin(), tcTdPairs.end() );

    while ( !tcTdPairs.empty() ) {
      const auto tc = tcTdPairs.back().first;
      const auto td = tcTdPairs.back().second;
      tcTdPairs.pop_back();

      ts.push_back( tc );
      ts.push_back( td );
      const auto G1 = G0 - distances( tb, tc );
      const auto G2 = G1 + distances( tc, td );
      auto gain = G2 - distances( ts.front(), ts.back() );
      if ( gain > tolerance && makesTour_( ts, tour, position ) ) {
        return gain;
      }

      if ( depth + 1 < k ) {
        added.push_back( make_pair( tb, tc ) );
        added.push_back( make_pair( tc, tb ) );
        removed.push_back( make_pair( tc, td ) );
        removed.push_back( make_pair( td, tc ) );
        gain = getBestKOptMove_( depth + 1,
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
                                        vector< size_t >& position,
                                        vector< bool >& dontLook,
                                        const vector< size_t >& betterTour,
                                        const VDistances& distances,
                                        const vector< vector< size_t > >& nearestNeighbors,
                                        bool linKernighan )
  {
    updatePosition( tour, position );
    auto anyChange = false;
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
        removed = { make_pair( t1, t2 ), make_pair( t2, t1 ) };
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
          auto gain = getBestKOptMove_( 1,
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
            for ( const auto& th : tsHistory ) {
              dontLook[ th ] = false;
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
          swap( tour, tourCopy );
          updatePosition( tour, position );
        }
      }
      if ( !found ) {
        dontLook[ t1 ] = true;
      }
    }
    return anyChange;
  }

  void kSwapKick( size_t k,
                  vector< size_t >& tour,
                  vector< size_t >& position,
                  vector< bool >& dontLook )
  {
    updatePosition( tour, position );

    vector< size_t > swapEdges;
    swapEdges.reserve( k );
    while ( swapEdges.size() < k && swapEdges.size() < tour.size() / 2 ) {
      size_t i = rand() % tour.size();
      if ( find( swapEdges.begin(), swapEdges.end(), i ) == swapEdges.end() ) {
        swapEdges.push_back( i );
      }
    }
    sort( swapEdges.begin(), swapEdges.end(), [&] ( const auto& ta, const auto& tb ) { return position[ ta ] < position[ tb ]; } );
    vector< size_t > ts;
    ts.reserve( 2 * swapEdges.size() );
    for ( size_t i = 0; i < swapEdges.size(); ++i ) {
      ts.push_back( swapEdges[ i ] );
      ts.push_back( next_( ts.back(), tour, position ) );
    }
    auto tourCopy = tour;
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

  bool INLINE_ATTRIBUTE improveTourKOpt_( size_t k,
                                          bool linKernighan,
                                          vector< size_t >& tour,
                                          const vector< size_t >& betterTour,
                                          const VDistances& distances,
                                          const vector< vector< size_t > >& nearestNeighbors )
  {
    vector< size_t > position( tour.size() );
    vector< bool > dontLook( tour.size(), false );
    vector< bool > dontLook4( tour.size(), false );
    bool change = false;
    bool change4 = true;
    auto innerLoop = [ k, linKernighan, &tour, &betterTour, &distances, &nearestNeighbors, &position, &dontLook ] {
      if ( linKernighan ) {
        return kOptOuterLoop_( k, tour, position, dontLook, betterTour, distances, nearestNeighbors, linKernighan );
      }
      switch ( k ) {
        case 2: return twoOptInnerLoop_( tour, position, dontLook, distances, nearestNeighbors );
        case 3: return threeOptInnerLoop_( tour, position, dontLook, distances, nearestNeighbors );
        default: return kOptOuterLoop_( k, tour, position, dontLook, betterTour, distances, nearestNeighbors, linKernighan );
      }
    };
    while ( change4 ) {
      change4 = false;
      while ( innerLoop() ) {
        change = true;
      }

      if ( k > 3 ) {
        fill( dontLook4.begin(), dontLook4.end(), false );
        while ( twoThreeOptInnerLoop_( tour, position, dontLook4, dontLook, distances, nearestNeighbors ) ) {
          change4 = true;
          change = true;
        }
      }
    }
    return change;
  }

  bool INLINE_ATTRIBUTE improveTourIterated_( function< bool ( vector< size_t >&,
                                                               vector< size_t >&,
                                                               vector< bool >&,
                                                               const VDistances&,
                                                               const vector< vector< size_t > >& ) > improveTour,
                                              size_t iterations,
                                              bool useGain23,
                                              vector< size_t >& tour,
                                              vector< size_t >& betterTour,
                                              const VDistances& distances,
                                              const vector< vector< size_t > >& nearestNeighbors )
  {
    assert( tour.size() > 8 );
    auto change = false;
    vector< bool > dontLook( tour.size(), false );
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
          vector< bool > dontLook4( tour.size(), false );
          while ( twoThreeOptInnerLoop_( tour, position, dontLook4, dontLook, distances, nn ) ) {
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
      if ( i + 1 == iterations ) {
        break;
      }
      tour = bestTour;
      kSwapKick( 6, tour, position, dontLook );
    }
    tour = bestTour;
    return change;
  }

  bool INLINE_ATTRIBUTE improveTourIteratedKOpt_( size_t k,
                                                  bool linKernighan,
                                                  size_t iterations,
                                                  bool useGain23,
                                                  vector< size_t >& tour,
                                                  vector< size_t >& betterTour,
                                                  const VDistances& distances,
                                                  const vector< vector< size_t > >& nearestNeighbors )
  {
    auto improveTour = [ k, linKernighan, &betterTour ] ( vector< size_t >& tour,
                                                          vector< size_t >& position,
                                                          vector< bool >& dontLook,
                                                          const VDistances& distances,
                                                          const vector< vector< size_t > >& nearestNeighbors ) {
      if ( linKernighan ) {
        return kOptOuterLoop_( k, tour, position, dontLook, betterTour, distances, nearestNeighbors, linKernighan );
      }
      switch ( k ) {
        case 2: return twoOptInnerLoop_( tour, position, dontLook, distances, nearestNeighbors );
        case 3: return threeOptInnerLoop_( tour, position, dontLook, distances, nearestNeighbors );
        default: return kOptOuterLoop_( k, tour, position, dontLook, betterTour, distances, nearestNeighbors, linKernighan );
      }
    };
    return improveTourIterated_( improveTour, iterations, useGain23, tour, betterTour, distances, nearestNeighbors );
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
        assert( fabs( distances( i, j ) - distances( j, i ) ) < 1e-9 );
        assert( distances( i, j ) > 0.0 );
        assert( distances( j, i ) > 0.0 );
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

        nn[ i ].erase( remove_if( nn[ i ].begin() + hn[ i ].size(), nn[ i ].end(),
                                  [&] ( size_t index ) { return find( hn[ i ].begin(), hn[ i ].end(), index ) != hn[ i ].end(); } ),
                       nn[ i ].end() );
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

  for ( const auto& k : { 2, 3, 5 } ) {
    tour = tourGreedy;
    vector< size_t > betterTour;
    bool linKernighan = false;
    double start( clock() );
    improveTourKOpt_( k, linKernighan, tour, betterTour, distances, nearestNeighbors30 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << k << "-opt tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
    assert( isTour_( tour, distances ) );
  }

  for ( const auto& k : { 3 } ) {
    tour = tourGreedy;
    vector< size_t > betterTour;
    bool linKernighan = false;
    double start( clock() );
    improveTourIteratedKOpt_( k, linKernighan, 1000, false, tour, betterTour, distances, nearestNeighbors30 );
    double time( ( clock() - start ) / CLOCKS_PER_SEC );
    cerr << "I" << k << "    tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
  }

  for ( size_t k = 5; k < 6; ++k ) {
    if ( true ) {
      tour = tourGreedy;
      vector< size_t > betterTour;
      vector< size_t > tour1( tour );
      vector< size_t > tour2( tour );
      vector< vector< size_t > > nn( nearestNeighbors5 );
      vector< size_t > bestTour( tour );
      vector< vector< double > > points;
      double start( clock() );
      for ( size_t i = 0; i < 10; ++i ) {
        tour = getHelsgaunInitialTour_( nn, helsgaun10, helsgaunDistances10, bestTour );
        bool useGain23 = true;
        while ( improveTourKOpt_( k, useGain23, tour, betterTour, distances, nn ) ) {
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
  }

  for ( size_t k = 5; k < 6; ++k ) {
    if ( true ) {
      tour = tourGreedy;
      vector< size_t > betterTour;
      bool linKernighan = true;
      bool useGain23 = false;
      double start( clock() );
      improveTourIteratedKOpt_( k, linKernighan, 100, useGain23, tour, betterTour, distances, nearestNeighbors5 );
      double time( ( clock() - start ) / CLOCKS_PER_SEC );
      cerr << k << "-ILK tour distance: " << getLength_( tour, distances ) << ", time: " << time << endl;
    }
  }

  if ( false ) {
    printTour_( "a", tourGreedy );
  }

  return tour;
}
