#ifndef TRAVELING_SALESPERSON_H
#define TRAVELING_SALESPERSON_H

#include <cstddef>
#include <vector>
#include <functional>

namespace TravelingSalespersonProblemSolver {
  class VDistances {
    public:
      VDistances( const std::vector< std::vector< double > >& points,
                  std::function< double ( double ) > rounding );
      virtual double operator()( size_t i, size_t j ) const = 0;
      virtual size_t size() const = 0;
      bool empty() const;
    protected:
      double computeDistance_( const std::vector< double >& point1, const std::vector< double >& point2 ) const;
      std::function< double ( double ) > rounding_;
  };

  class MatrixDistances : public VDistances {
    public:
      MatrixDistances( const std::vector< std::vector< double > >& points,
                       std::function< double ( double ) > rounding = [] ( double d ) { return d; } );
      virtual double operator()( size_t i, size_t j ) const;
      virtual size_t size() const;
      void setMatrix( std::vector< std::vector< double > >& distances );
    private:
      std::vector< std::vector< double > > distances_;
  };

  class MatrixRoundedDistances : public MatrixDistances {
    public:
      MatrixRoundedDistances( const std::vector< std::vector< double > >& points );
  };

  class OnTheFlyDistances : public VDistances {
    public:
      OnTheFlyDistances( const std::vector< std::vector< double > >& points,
                       std::function< double ( double ) > rounding = [] ( double d ) { return d; } );
      virtual double operator()( size_t i, size_t j ) const;
      virtual size_t size() const;

    private:
      const std::vector< std::vector< double > > points_;
  };

  class OnTheFlyRoundedDistances : public OnTheFlyDistances {
    public:
      OnTheFlyRoundedDistances( const std::vector< std::vector< double > >& points );
  };

  std::vector< size_t > computeTour( const VDistances& distances );
}

#endif
