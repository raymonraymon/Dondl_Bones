#ifndef __ENERGYDEFINES_H
#define __ENERGYDEFINES_H

namespace aol {


// Energy is quadratic: E(u) = Mu - F \cdot u, i.e. Hessian is constant
template <typename DataTypeContainer >
class LinearEnergyOp {
protected :
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  
  const SparseMatrixType &_Hessian;
  const VectorType &_rhs;
  
public:
  LinearEnergyOp( const SparseMatrixType & Hessian, const VectorType & rhs ) : _Hessian( Hessian ), _rhs( rhs ){ }

  RealType evaluateEnergy ( const VectorType &Arg ) const { return 0.5 * (_Hessian * Arg ).dot( Arg ) - _rhs.dot( Arg );}
  void evaluateEnergy( const VectorType &Arg, RealType & energy ) const { energy = this->evaluateEnergy(Arg);}
  void evaluateJacobian( const VectorType &Arg, VectorType & jacobian ) const { jacobian = _Hessian * Arg - _rhs;}
  void evaluateHessian( const VectorType &/*Arg*/, SparseMatrixType & hessian ) const { hessian = _Hessian;}
  const SparseMatrixType & getHessian( ) const {return _Hessian;}
};


template <typename DataTypeContainer >
class NonlinearEnergyOp {

public:
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  
public:
  
  NonlinearEnergyOp() { }

  // Destroy polymorphic Ops correctly, important!
  virtual ~NonlinearEnergyOp () {}
  
  virtual void evaluateEnergy( const VectorType &Arg, RealType & energy ) const = 0;
  RealType evaluateEnergy ( VectorType &Arg ) { RealType energy = 0; evaluateEnergy( Arg, energy ); return energy; }
  virtual void evaluateJacobian( const VectorType &Arg, VectorType & jacobian ) const = 0;
  virtual void evaluateHessian( const VectorType &Arg, SparseMatrixType & hessian ) const = 0;
  virtual void evaluateTripletListHessianSym( const VectorType &Arg, std::vector<typename DataTypeContainer::TripletType> & tripletListHessian ) const = 0;
  
};


template <typename DataTypeContainer >
class NonlinearConstraintOps {

protected:
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  
  const int _numConstraints;
  
public:
  
  NonlinearConstraintOps( const int numConstraints) : _numConstraints ( numConstraints ) { }

  // Destroy polymorphic Ops correctly, important!
  virtual ~NonlinearConstraintOps () {}
  
  virtual void evaluateEnergy( const int numConstraint, const VectorType &Arg, RealType & energy ) const = 0;
  RealType evaluateEnergy ( const int numConstraint, VectorType &Arg ) { RealType energy = 0; evaluateEnergy( numConstraint, Arg, energy ); return energy; }
  virtual void evaluateJacobian( const int numConstraint, const VectorType &Arg, VectorType & jacobian ) const = 0;
  virtual void evaluateHessian( const int numConstraint, const VectorType &Arg, SparseMatrixType & hessian ) const = 0;
  
  const int getNumConstraints( ) const { return _numConstraints;}
  
};


template <typename DataTypeContainer >
class NonlinearConstraintOps_FromVector : public NonlinearConstraintOps<DataTypeContainer> {

protected:
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  
  const std::vector< std::reference_wrapper<aol::NonlinearEnergyOp<DataTypeContainer>> > &_constraintOps;
  
public:
  
  NonlinearConstraintOps_FromVector( const std::vector< std::reference_wrapper<aol::NonlinearEnergyOp<DataTypeContainer>> > &constraintOps) : 
  NonlinearConstraintOps<DataTypeContainer> ( constraintOps.size() ),
  _constraintOps ( constraintOps ) { }

  // Destroy polymorphic Ops correctly, important!
  ~NonlinearConstraintOps_FromVector () {}
  
  void evaluateEnergy( const int numConstraint, const VectorType &Arg, RealType & energy ) const {
      _constraintOps[numConstraint].get().evaluateEnergy( Arg, energy );
  };
  void evaluateJacobian( const int numConstraint, const VectorType &Arg, VectorType & jacobian ) const {
      _constraintOps[numConstraint].get().evaluateJacobian( Arg, jacobian );
  };
  void evaluateHessian( const int numConstraint, const VectorType &Arg, SparseMatrixType & hessian ) const {
      _constraintOps[numConstraint].get().evaluateHessian( Arg, hessian );
  };
  
};





template <typename DataTypeContainer >
class SparseJacobianNonlinearConstraintOps {

protected:
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  typedef typename DataTypeContainer::TripletType TripletType;
  
  const int _numConstraints;
  
public:
  
  SparseJacobianNonlinearConstraintOps( const int numConstraints) : _numConstraints ( numConstraints ) { }

  // Destroy polymorphic Ops correctly, important!
  virtual ~SparseJacobianNonlinearConstraintOps () {}
  
  virtual void evaluate( const VectorType &Arg, VectorType & constraintVec ) const = 0;
  virtual void evaluateJacobian( const VectorType &Arg, std::vector<TripletType> & jacobian ) const = 0;
  virtual void evaluateHessian( const VectorType &Arg, std::vector<std::vector<TripletType>> &hessian ) const = 0;
  
  const int getNumConstraints( ) const { return _numConstraints;}
  virtual const int sizeJacobian( ) const = 0;
  virtual const int sizeHessian( ) const = 0;
};


}//end namespace


#endif
