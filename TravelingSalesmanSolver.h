#ifndef TRAVELING_SALESMAN_H
#define TRAVELING_SALESMAN_H

#include <cstddef>
#include <vector>

namespace TravelingSalesmanSolver {
  std::vector< size_t > computePath( const std::vector< std::vector< double > >& distances );
}

#endif
