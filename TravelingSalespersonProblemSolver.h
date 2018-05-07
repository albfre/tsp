#ifndef TRAVELING_SALESPERSON_PROBLEM_SOLVER_H
#define TRAVELING_SALESPERSON_PROBLEM_SOLVER_H

#include "Distances.h"

namespace TravelingSalespersonProblemSolver {
  std::vector< size_t > computeTour( const VDistances& distances );
}

#endif // TRAVLING_SALESPERSON_PROBLEM_SOLVER_H
