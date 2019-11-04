#ifndef __IPOPTBOXCONSTRAINTSOLVER_H
#define __IPOPTBOXCONSTRAINTSOLVER_H

#include <ipoptIncludes.h>
#include <SolverInfo.h>


#ifdef USE_IPOPT


/* use IPOPT to solve a constraint minimization problem of type 
 *  minimize f(x) 
 * over all x s.t. 
 * x_l \leq x \leq x_u
 * g_l(x) \leq g(x) \leq g_u(x)
 */
template<typename DataTypeContainer>
class IpoptBoxConstraintFirstOrderSolverInterface : public Ipopt::TNLP {

protected:
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  
  const aol::NonlinearEnergyOp<DataTypeContainer> &_energyOp;

  const VectorType &_startingPoint;
  VectorType &_nextIter;
 
  unsigned int _numTotalDofs;
  unsigned int _numConstraints;
  
  const RealType _x_l, _x_u;
  RealType _nlpError;
  int _iter;

public:
  IpoptBoxConstraintFirstOrderSolverInterface ( const aol::NonlinearEnergyOp<DataTypeContainer> &energyOp,
                                                const VectorType &startingPoint, VectorType &nextIter,
                                                const RealType x_l = -2.e+19,   const RealType x_u = 2.e+19 )
  : _energyOp ( energyOp ),
    _startingPoint ( startingPoint ), _nextIter ( nextIter ),
    _numTotalDofs ( startingPoint.size () ), _numConstraints ( 0 ),
    _x_l ( x_l ), _x_u ( x_u ) { }
    
  //returns info about the nlp
  virtual bool get_nlp_info ( Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style ) {
    // The number of DOFs.
    n = _numTotalDofs;
    // The number of constraints.
    m = _numConstraints;
    // The number of nonzeros in the jacobian: Possibly every entry.
    nnz_jac_g = _numConstraints * _numTotalDofs;
    //Nonzeros in hessian: 0 because not implemented
    nnz_h_lag = 0;
    //Use C index style
    index_style = C_STYLE;

    return true;
  }

  virtual bool get_bounds_info ( Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, 
                                 Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u ) {
      
    for ( Ipopt::Index i = 0; i < n; ++i ) {
        x_l[i] = _x_l;
        x_u[i] = _x_u;
    }

    return true;
  }

  virtual bool get_starting_point ( Ipopt::Index n, bool init_x,  Ipopt::Number* x, 
                                    bool init_z, Ipopt::Number* /*z_L*/, Ipopt::Number* /*z_U*/, 
                                    Ipopt::Index /*m*/, bool init_lambda, Ipopt::Number* /*lambda*/ ) {
    init_x == true; init_z == false; init_lambda == false;
    for ( int j = 0; j < _startingPoint.size (); ++j ) x[j] = _startingPoint[j];
    return true;
  }

  virtual bool eval_f ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number& obj_value ) {

    // Convert
    VectorType v ( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    RealType s;
    _energyOp.evaluateEnergy( v, s );
    obj_value = s;
    return true;

  }

  virtual bool eval_grad_f ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f ) {

    // Convert
    VectorType v ( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

   VectorType res ( v.size() );
   _energyOp.evaluateJacobian ( v, res );
   for ( int i = 0; i < n; ++i )
     grad_f[i] = res[i];
   return true;

  }

  virtual bool eval_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number *g ) {
    return true;
  }

  virtual bool eval_jac_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                            Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values ) {
    return true;
  }

  virtual bool intermediate_callback ( Ipopt::AlgorithmMode /*mode*/, Ipopt::Index iter, 
                                       Ipopt::Number /*obj_value*/,  Ipopt::Number /*inf_pr*/,  Ipopt::Number /*inf_du*/, Ipopt::Number /*mu*/,  Ipopt::Number /*d_norm*/, Ipopt::Number /*regularization_size*/, Ipopt::Number /*d_du*/, Ipopt::Number /*d_pr*/,  Ipopt::Index /*ls_trials*/, 
                                       const Ipopt::IpoptData* ip_data,
                                       Ipopt::IpoptCalculatedQuantities* ip_cq ) {

    // Retrieve primal variables (code similar to code from the IPOPT documentation):
    Ipopt::TNLPAdapter* tnlpAdapter = NULL;

    if ( ip_cq != 0 ) {
      Ipopt::OrigIpoptNLP* origNLP;
      origNLP = dynamic_cast< Ipopt::OrigIpoptNLP*> ( Ipopt::GetRawPtr ( ip_cq->GetIpoptNLP () ) );
      // If in restoration mode, origNLP will be NULL. Quit method in this case
      if ( origNLP == NULL )  return true;
      tnlpAdapter = dynamic_cast < Ipopt::TNLPAdapter* > ( GetRawPtr ( origNLP->nlp () ) );
    }
    else throw std::invalid_argument ( aol::strprintf ( "ip_cq == NULL in FILE %s at line %d.", __FILE__, __LINE__ ).c_str() );
    double *x = new double[_startingPoint.size ()];   // n == _numDofs
    tnlpAdapter->ResortX ( *ip_data->curr ()->x (), x );
    delete[] x;

    // Plot solution.
    _iter = iter;
    
    return true;
  }

  virtual void finalize_solution ( Ipopt::SolverReturn /*status*/, Ipopt::Index /*n*/,  const Ipopt::Number* x,
                                   const Ipopt::Number* /*z_L*/,  const Ipopt::Number* /*z_U*/, 
                                   Ipopt::Index /*m*/,  const Ipopt::Number* /*g*/, const Ipopt::Number* /*lambda*/, 
                                   Ipopt::Number /*obj_value*/,  const Ipopt::IpoptData* /*ip_data*/,  Ipopt::IpoptCalculatedQuantities* ip_cq ) {
    for ( int i = 0; i < _nextIter.size (); ++i ) _nextIter[i] = x[i];
    _nlpError = ip_cq->curr_nlp_error();
  }
  
public : 
    RealType getNLPError( ) const { return _nlpError; }
    int getNumIterations() const { return _iter; }
  
};




template <typename DataTypeContainer>
class IpoptBoxConstraintFirstOrderSolver {
protected:
 
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
    
  const aol::NonlinearEnergyOp<DataTypeContainer> &_energyOp;
  RealType _ipoptTol;
  RealType _MaxIterations;
  const RealType _x_l, _x_u;
  const int _ipoptPrintLevel;

public:
  IpoptBoxConstraintFirstOrderSolver ( const aol::NonlinearEnergyOp<DataTypeContainer> &E,
                                       const int MaxIterations, const RealType Tolerance,
                                       const RealType x_l = - 2.e+19, const RealType x_u = 2.e+19,
                                       const int ipoptPrintLevel = 5
                 )
  : _energyOp ( E ),
    _ipoptTol ( Tolerance ), _MaxIterations ( MaxIterations ),
    _x_l( x_l), _x_u( x_u ),
    _ipoptPrintLevel ( ipoptPrintLevel ) { }

  ~IpoptBoxConstraintFirstOrderSolver() {}
  
  void solve ( const VectorType &StartPosition, VectorType &Solution, SolverInfo<DataTypeContainer> &solverInfo ) const {

    // Set up masking operators:
    Ipopt::SmartPtr<IpoptBoxConstraintFirstOrderSolverInterface<DataTypeContainer> > tsOpt 
    = new IpoptBoxConstraintFirstOrderSolverInterface<DataTypeContainer> ( _energyOp, StartPosition, Solution, _x_l, _x_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp = new Ipopt::IpoptApplication ();
    
    //Derivative Test
#ifdef IPOPTDERIVATIVETEST
    ipoptApp->Options ()->SetStringValue  ( "derivative_test", "first-order" );
//     ipoptApp->Options ()->SetStringValue  ( "derivative_test_print_all", "yes" );
    ipoptApp->Options ()->SetStringValue  ( "derivative_test_print_all", "no" );
    ipoptApp->Options ()->SetNumericValue ( "derivative_test_tol", 1.e-4 );
    ipoptApp->Options ()->SetNumericValue ( "derivative_test_perturbation", 1.e-6 );
#endif //IPOPTDERIVATIVETEST
    
//     ipoptApp->Options()->SetStringValue( "output_file", "ipoptOutputTest.txt" );
    
    // Tolerance
    ipoptApp->Options()->SetNumericValue ( "tol", _ipoptTol * 1e-2 );
    ipoptApp->Options()->SetNumericValue ( "acceptable_tol", _ipoptTol );

    // Enable quasi-Newton hessian approximation
    ipoptApp->Options ()->SetStringValue ( "hessian_approximation", "limited-memory" );
    ipoptApp->Options ()->SetIntegerValue ( "max_iter", _MaxIterations );
    ipoptApp->Options ()->SetStringValue ( "linear_solver", "MUMPS" );

    ipoptApp->Options()->SetIntegerValue( "print_level", _ipoptPrintLevel );
    
    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();

    Ipopt::SmartPtr<Ipopt::TNLP> nlp;
    
    //! Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
    
    outputIpoptStatus ( ipoptStatus, true );
    
//     Ipopt::SmartPtr< Ipopt::IpoptCalculatedQuantities >  ipCQ ( nlp, const Ipopt::SmartPtr< IpoptData > &ip_data)
//     double error = ipCQ->curr_nlp_error();
    solverInfo.setSolverStatus( getIpoptStatus( ipoptStatus ).c_str() );
    solverInfo.setError( tsOpt->getNLPError() );
    solverInfo.setNumIterations( tsOpt->getNumIterations() );
    
  }
  
  void solve ( const VectorType &StartPosition, VectorType &Solution ) const {
    SolverInfo<DataTypeContainer> solverInfo;
    this->solve( StartPosition, Solution, solverInfo );   
  }

};






/* use IPOPT to solve a constraint minimization problem of type 
 *  minimize f(x) 
 * over all x s.t. 
 * x_l \leq x \leq x_u
 * g_l(x) \leq g(x) \leq g_u(x)
 */
template<typename DataTypeContainer >
class IpoptBoxConstraintSecondOrderSolverInterface : public Ipopt::TNLP {

protected:
    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType VectorType;
    typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
    typedef typename DataTypeContainer::MaskType MaskType;
    typedef typename DataTypeContainer::TripletType TripletType;
    
  const aol::NonlinearEnergyOp<DataTypeContainer> & _energyOp;

  const VectorType &_startingPoint;
  VectorType &_nextIter;
  
  unsigned int _numTotalDofs;
  
  const RealType _x_l, _x_u;
public:
  IpoptBoxConstraintSecondOrderSolverInterface ( const aol::NonlinearEnergyOp<DataTypeContainer> &Energy,
                           const VectorType &startingPoint, VectorType &nextIter,
                           const RealType x_l = -2.e+19, const RealType x_u = 2.e+19 )
  : _energyOp ( Energy ),
    _startingPoint ( startingPoint ), 
    _nextIter ( nextIter ),
    _numTotalDofs ( startingPoint.size () ),
    _x_l ( x_l ), _x_u ( x_u ) { }
    

  //returns info about the nlp
  virtual bool get_nlp_info ( Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style ) {
    // The number of DOFs.
    n = _numTotalDofs;
    // The number of constraints.
    m = 0;
    // The number of nonzeros in the jacobian: Possibly every entry. 
    nnz_jac_g = 0;
    
    nnz_h_lag = 0;

   std::vector<TripletType> tripletListHessian;
   VectorType tmpVec ( n ); tmpVec.setZero();
   _energyOp.evaluateTripletListHessianSym( tmpVec, tripletListHessian );
   nnz_h_lag = tripletListHessian.size();
    
    //Use C index style
    index_style = C_STYLE;

    return true;
  }

  virtual bool get_bounds_info ( Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, 
                                 Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u ) {
   for ( Ipopt::Index i = 0; i < n; ++i ) {
      x_l[i] = _x_l;
      x_u[i] = _x_u;
    }
    return true;
  }

  virtual bool get_starting_point ( Ipopt::Index n, bool init_x,  Ipopt::Number* x, 
                                    bool init_z, Ipopt::Number* /*z_L*/, Ipopt::Number* /*z_U*/, 
                                    Ipopt::Index /*m*/, bool init_lambda, Ipopt::Number* /*lambda*/ ) {
    init_x == true; init_z == false; init_lambda == false;
    for ( int j = 0; j < _startingPoint.size (); ++j ) x[j] = _startingPoint[j];
    return true;
  }

  virtual bool eval_f ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number& obj_value ) {
    // Convert
    VectorType v ( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    //obj_value = _energyOp.evaluateEnergy( v );
    RealType energy = 0.;
    _energyOp.evaluateEnergy( v, energy );
    obj_value = energy;
    return true;
  }

  virtual bool eval_grad_f ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f ) {
    // Convert
    VectorType v ( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

   VectorType res ( v.size() );
   _energyOp.evaluateJacobian( v, res );
   for ( int i = 0; i < n; ++i ) grad_f[i] = res[i];
   return true;
  }

  virtual bool eval_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number *g ) {
    return true;
  }

  virtual bool eval_jac_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                            Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values ) {
    return true;
  }
  
  // Hessian of Lagrangian: L = sigma * D^2 E + lambda D^2 Constraint
  virtual bool eval_h( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor,
                       Ipopt::Index m, const Ipopt::Number* lambda, bool new_lambda, 
                       Ipopt::Index nele_hess, Ipopt::Index* iRow,Ipopt::Index* jCol, Ipopt::Number* values){
    if (values == NULL) {
        // return the structure.
        Ipopt::Index idx=0;

        std::vector<TripletType> tripletListHessian;
        VectorType tmpVec ( n ); tmpVec.setZero();
        _energyOp.evaluateTripletListHessianSym( tmpVec, tripletListHessian );
        for( int index = 0; index < tripletListHessian.size(); ++index ){
              iRow[idx] = tripletListHessian[index].row();
              jCol[idx] = tripletListHessian[index].col();
              idx++;
        }
    }
    else {
        // return the values
        // Convert
        VectorType v ( n ); for ( int i = 0; i < n; ++i ) v[i] = x[i];
      
        
        Ipopt::Index idx=0;
        std::vector<TripletType> tripletListHessian;
        _energyOp.evaluateTripletListHessianSym( v, tripletListHessian );
        for( int index = 0; index < tripletListHessian.size(); ++index ){
              values[idx] = obj_factor * tripletListHessian[index].value();
              idx++;
        }
    }


    return true;
 }

  virtual bool intermediate_callback ( Ipopt::AlgorithmMode /*mode*/, Ipopt::Index iter, 
                                       Ipopt::Number /*obj_value*/,  Ipopt::Number /*inf_pr*/,  Ipopt::Number /*inf_du*/, Ipopt::Number /*mu*/,  Ipopt::Number /*d_norm*/, Ipopt::Number /*regularization_size*/, Ipopt::Number /*d_du*/, Ipopt::Number /*d_pr*/,  Ipopt::Index /*ls_trials*/, 
                                       const Ipopt::IpoptData* ip_data,
                                       Ipopt::IpoptCalculatedQuantities* ip_cq ) {

    // Retrieve primal variables (code similar to code from the IPOPT documentation):
    Ipopt::TNLPAdapter* tnlpAdapter = NULL;

    if ( ip_cq != 0 ) {
      Ipopt::OrigIpoptNLP* origNLP;
      origNLP = dynamic_cast< Ipopt::OrigIpoptNLP*> ( Ipopt::GetRawPtr ( ip_cq->GetIpoptNLP () ) );
      // If in restoration mode, origNLP will be NULL. Quit method in this case
      if ( origNLP == NULL )  return true;
      tnlpAdapter = dynamic_cast < Ipopt::TNLPAdapter* > ( GetRawPtr ( origNLP->nlp () ) );
    }
    else throw std::invalid_argument ( aol::strprintf ( "ip_cq == NULL in FILE %s at line %d.", __FILE__, __LINE__ ).c_str() );
    double *x = new double[_startingPoint.size ()];   // n == _numDofs
    tnlpAdapter->ResortX ( *ip_data->curr ()->x (), x );
    delete[] x;

    // Plot solution.
    
    return true;
  }

  virtual void finalize_solution ( Ipopt::SolverReturn /*status*/, Ipopt::Index /*n*/,  const Ipopt::Number* x,
                                   const Ipopt::Number* /*z_L*/,  const Ipopt::Number* /*z_U*/, 
                                   Ipopt::Index /*m*/,  const Ipopt::Number* /*g*/, const Ipopt::Number* /*lambda*/, 
                                   Ipopt::Number /*obj_value*/,  const Ipopt::IpoptData* /*ip_data*/,  Ipopt::IpoptCalculatedQuantities* /*ip_cq*/ ) {
    for ( int i = 0; i < _nextIter.size (); ++i ) _nextIter[i] = x[i];
  }
  
};




template <typename DataTypeContainer >
class IpoptBoxConstraintSecondOrderSolver {
protected:
 
    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType VectorType;
    typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
    
  const aol::NonlinearEnergyOp<DataTypeContainer> &_energyOp;
  const RealType _ipoptTol;
  const RealType _MaxIterations;
  const RealType _x_l, _x_u;
  const int _ipoptPrintLevel;
public:
  IpoptBoxConstraintSecondOrderSolver ( const aol::NonlinearEnergyOp<DataTypeContainer> &linearEnergy,
                                        const int MaxIterations, const RealType Tolerance,
                                        const RealType x_l = - 2.e+19,const RealType x_u = 2.e+19,
                                        const int ipoptPrintLevel = 5  )
  : _energyOp( linearEnergy ),
    _ipoptTol ( Tolerance ),
    _MaxIterations ( MaxIterations ),
    _x_l( x_l), _x_u( x_u ),
    _ipoptPrintLevel ( ipoptPrintLevel ) { }

  virtual ~IpoptBoxConstraintSecondOrderSolver() {}

  void solve ( const VectorType &StartPosition, VectorType &Solution ) const {

    // Set up masking operators:
    Ipopt::SmartPtr<IpoptBoxConstraintSecondOrderSolverInterface<DataTypeContainer> > tsOpt 
    = new IpoptBoxConstraintSecondOrderSolverInterface<DataTypeContainer> ( _energyOp, StartPosition, Solution, _x_l, _x_u );
    
    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp = new Ipopt::IpoptApplication ();
    
#ifdef DERIVATIVETEST
//     ipoptApp->Options ()->SetStringValue  ( "derivative_test", "first-order" );
    ipoptApp->Options ()->SetStringValue  ( "derivative_test", "second-order" );
    ipoptApp->Options ()->SetStringValue  ( "derivative_test_print_all", "yes" );
    ipoptApp->Options ()->SetNumericValue ( "derivative_test_tol", 1.e-4 );
    ipoptApp->Options ()->SetNumericValue ( "derivative_test_perturbation", 1.e-6 );
#endif //DERIVATIVETEST
    
    // Tolerance
    ipoptApp->Options()->SetNumericValue ( "tol", _ipoptTol * 1.e-2 );
    ipoptApp->Options()->SetNumericValue ( "acceptable_tol", _ipoptTol );
    //ipoptApp->Options()->SetNumericValue( "acceptable_constr_viol_tol", 1.e-6 );

    // Enable quasi-Newton hessian approximation
    ipoptApp->Options()->SetStringValue ( "hessian_approximation", "exact" );
    
    ipoptApp->Options()->SetIntegerValue ( "max_iter", _MaxIterations );
    ipoptApp->Options()->SetStringValue ( "linear_solver", "MUMPS" );

    //set print_level (from 0 - 12 )
    ipoptApp->Options()->SetIntegerValue( "print_level", _ipoptPrintLevel );
    
    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();

    Ipopt::SmartPtr<Ipopt::TNLP> nlp;
    
    //! Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
    // true means, output only if solver failed
    outputIpoptStatus ( ipoptStatus, true );
    
  }

};


#endif //USE_IPOPT

#endif //__IPOPTBOXCONSTRAINTSOLVER_H
