#ifndef DISTANCE_MATRICES_H
#define DISTANCE_MATRICES_H


namespace TravelingSalespersonProblemSolver {
  class VDistances {
    public:
      VDistances( const std::vector< std::vector< double > >& points,
                  std::function< double ( double ) > rounding );
      virtual double operator()( size_t i, size_t j ) const = 0;
      virtual size_t size() const = 0;
      bool empty() const { return size() == 0; }
    protected:
      ouble computeDistance_( const std::vector< double >& point1, const std::vector< double >& point2 ) const;
      std::function< double ( double ) > rounding_;
  };

  class MatrixDistances : public VDistances {
    public:
      MatrixDistances( const std::vector< std::vector< double > >& points,
                       std::function< double ( double ) > rounding = [] ( double d ) { return d; } );
      virtual double operator()( size_t i, size_t j ) const;
      virtual size_t size() const { return distances_.size(); }
      void setMatrix( std::vector< std::vector< double > >& distances )
    private:
      std::vector< std::vector< double > > distances_;
  };

  class MatrixRoundedDistances : public MatrixDistances {
    public:
      MatrixRoundedDistances( const std::vector< std::vector< double > >& points );
  };

  class SparseDistances : public VDistances {
    public:
      SparseDistances( const std::vector< std::vector< double > >& points,
                       std::function< double ( double ) > rounding = [] ( double d ) { return d; } );
      virtual double operator()( size_t i, size_t j ) const;
      virtual size_t size() const { return points_.size(); }

    private:
      const std::vector< std::vector< double > > points_;
  };

  class SparseRoundedDistances : public SparseDistances {
    public:
      SparseRoundedDistances( const std::vector< std::vector< double > >& points );
  }
}

#endif
