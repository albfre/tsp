#ifndef TRAVELING_SALESPERSON_H
#define TRAVELING_SALESPERSON_H

#include <cstddef>
#include <vector>

namespace TravelingSalespersonProblemSolver {
  std::vector< size_t > computePath( const std::vector< std::vector< double > >& distances );
}

#endif
