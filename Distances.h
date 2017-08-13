#ifndef TRAVELING_SALESPERSON_PROBLEM_SOLVER_DISTANCES_H
#define TRAVELING_SALESPERSON_PROBLEM_SOLVER_DISTANCES_H

#include <cstddef>
#include <vector>
#include <functional>

namespace TravelingSalespersonProblemSolver {
  class VDistances {
    public:
      VDistances( const std::vector< std::vector< double > >& points,
                  std::function< double ( double ) > rounding );
      virtual ~VDistances() = 0;
      virtual double operator()( size_t i, size_t j ) const = 0;
      virtual size_t size() const = 0;
      bool empty() const;
    protected:
      double computeDistance_( const double* point1, const double* point2, size_t pointDimension ) const;
      double computeDistance_( const std::vector< double >& point1, const std::vector< double >& point2 ) const;
      std::function< double ( double ) > rounding_;
  };

  class MatrixDistances : public VDistances {
    public:
      MatrixDistances( const std::vector< std::vector< double > >& points,
                       std::function< double ( double ) > rounding = [] ( double d ) { return d; } );
      virtual ~MatrixDistances();
      virtual double operator()( size_t i, size_t j ) const override;
      virtual size_t size() const override;
      void setMatrix( std::vector< std::vector< double > >& distances );
    private:
      size_t size_;
      std::vector< double > distances_;
  };

  class MatrixRoundedDistances : public MatrixDistances {
    public:
      MatrixRoundedDistances( const std::vector< std::vector< double > >& points );
      virtual ~MatrixRoundedDistances();
  };

  class OnTheFlyDistances : public VDistances {
    public:
      OnTheFlyDistances( const std::vector< std::vector< double > >& points,
                         std::function< double ( double ) > rounding = [] ( double d ) { return d; } );
      virtual ~OnTheFlyDistances();
      virtual double operator()( size_t i, size_t j ) const override;
      virtual size_t size() const override;

    private:
      const size_t size_;
      const size_t pointDimension_;
      std::vector< double > points_;
  };

  class OnTheFlyRoundedDistances : public OnTheFlyDistances {
    public:
      OnTheFlyRoundedDistances( const std::vector< std::vector< double > >& points );
      virtual ~OnTheFlyRoundedDistances();
  };
}

#endif // TRAVELING_SALESPERSON_PROBLEM_SOLVER_DISTANCES_H
