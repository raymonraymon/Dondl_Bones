#ifndef LINESEARCHMETHODS_H
#define LINESEARCHMETHODS_H

#include <energyDefines.h>
#include "linearSolver.h"

//===================================================================================================================================
// Stepsize control
//===================================================================================================================================


//!==========================================================================================================
//! Class for Armijo stepsize control
template<typename StepsizeControlHelper>
class ArmijoStepsizeControl{
    
protected:
  typedef typename StepsizeControlHelper::DTContainer   DataTypeContainer;
  typedef typename DataTypeContainer::RealType          RealType;
  typedef typename DataTypeContainer::VectorType        VectorType;
  
  const StepsizeControlHelper _stepSizeControlHelper;
  
  // The tolerence parameters.
  RealType _sigma;

  mutable RealType _startTau;
  RealType _tauMin, _tauMax;
    
public:
  ArmijoStepsizeControl( const StepsizeControlHelper & stepSizeControlHelper, RealType sigma = 0.25, RealType startTau = 1.0, RealType tauMin = 1.e-7, RealType tauMax = 4. ) 
  : _stepSizeControlHelper ( stepSizeControlHelper ), _sigma(sigma), _startTau(startTau), _tauMin(tauMin), _tauMax(tauMax){}
  
  ~ArmijoStepsizeControl( ){}

  // Armijo timestepping, i.e. E[p + tau*d] < E[p] + _sigma * tau * <DE[p], d>
  RealType getTimestepWidth( VectorType &DescentDir, VectorType &CurrentPosition, const RealType OldTau ) const{
    // Makes sure that tauMin <= startTau <= tauMax. For OldTau == 0, one descent direction component in the last step didn't lead to an energy reduction. step is not done by setting tau to zero.
    RealType tau = std::min( std::max(OldTau, _tauMin), _tauMax );

    const RealType f  = _stepSizeControlHelper.evaluateLinesearchFnc( CurrentPosition );
    RealType fNew     = _stepSizeControlHelper.evaluateLinesearchFnc( CurrentPosition, DescentDir, tau );
    const RealType Df = _stepSizeControlHelper.evaluateLinesearchGrad( CurrentPosition, DescentDir );
    RealType G        = (fNew - f)/(tau * Df);

    if ( Df >= 0 ) return 0;

    //check
    if( G >= _sigma && f >= fNew){
      //time step too small
      do{
        tau *= 2.;
        fNew = _stepSizeControlHelper.evaluateLinesearchFnc( CurrentPosition, DescentDir, tau );
        G = (fNew - f)/(tau * Df);
      } while (G >= _sigma && f >= fNew && tau <= _tauMax);
      tau *= 0.5;
    }
    else{
      // time step too large
      do{
        if (tau > _tauMin)
          tau *= 0.5;
        fNew = _stepSizeControlHelper.evaluateLinesearchFnc( CurrentPosition, DescentDir, tau );
        G = (fNew - f)/(Df*tau);
      } while ( ( (G < _sigma || f < fNew) ) && (tau > _tauMin));
    }
    
      // update
      DescentDir *= tau;
      CurrentPosition += DescentDir;

    // If tau == tauMin it's very likely that the condition for the step size control is not satisfied. In such a case we don't want the step to be done.
    if( tau > _tauMin ) return tau;
    else return 0.;
  }
  
};
  
  
//!==========================================================================================================
//! Class for Powell-Wolfe stepsize control
template<typename StepsizeControlHelper>
class WolfeStepsizeControl{
    
protected:
  typedef typename StepsizeControlHelper::DTContainer   DataTypeContainer;
  typedef typename DataTypeContainer::RealType          RealType;
  typedef typename DataTypeContainer::VectorType        VectorType;
  
  const StepsizeControlHelper _stepSizeControlHelper;
  
  // The tolerence parameters.
  RealType _sigma;
  //! Additional parameter for getTimestepWidthWithPowellWolfeLineSearch.
  RealType _beta;

  mutable RealType _startTau;
  RealType _tauMin, _tauMax;
    
public:
  WolfeStepsizeControl( const StepsizeControlHelper & stepSizeControlHelper, RealType sigma, RealType beta, RealType startTau = 1.0, RealType tauMin = 1e-6, RealType tauMax = 8. ) 
  : _stepSizeControlHelper ( stepSizeControlHelper ), _sigma(sigma), _beta(beta), _startTau(startTau), _tauMin(tauMin), _tauMax(tauMax){}
  
  ~WolfeStepsizeControl( ){}

  // CurrentPosition contains the currently optimal point,
  // DescentDir the linesearch direction,
  // Algorithm taken from: Nocedal, Wright: Numerical Optimization, Alg.3.5-3.6
  RealType getTimestepWidth( VectorType &DescentDir, VectorType &CurrentPosition, const RealType /*oldTau*/ ) const {

    RealType Tau = 1.; // important for Newton-based methods for superlinear convergence, since the first steplength is chosen which satisfies the Wolfe conditions
    RealType TauOld = 0.;
    RealType dE = -1. * DescentDir.squaredNorm();

    RealType energy = _stepSizeControlHelper.evaluateLinesearchFnc( CurrentPosition );
    RealType energyBackup = energy;
  
    VectorType positionNew( CurrentPosition.size() );

    bool stepSizeFound = false;
    do {
      positionNew = CurrentPosition;
      positionNew += Tau * DescentDir;
      RealType energyNew = _stepSizeControlHelper.evaluateLinesearchFnc( positionNew );

      if ( energyNew > std::min( energy + _sigma * Tau * dE, energyBackup ) ) {
        // Armijo condition violated...
        Tau = zoomWolfe( DescentDir, CurrentPosition, energy, dE, energyBackup, TauOld, Tau );
        DescentDir *= Tau;
        stepSizeFound = true;

      } else {
        // Armijo condition fulfilled...
        RealType dENew = _stepSizeControlHelper.evaluateLinesearchGrad( positionNew, DescentDir );

        if ( std::abs( dENew ) <= - _beta * dE ) {
          // strong Wolfe condition (Armijo condition + curvature condition) fulfilled...
          DescentDir *= Tau;
          CurrentPosition = positionNew;
          energy = energyNew;
          stepSizeFound = true;

        } else {
          // curvature condition violated...
          if ( dENew >= 0. ) {
            // too long step
            Tau = zoomWolfe( DescentDir, CurrentPosition, energy, dE, energy, Tau, TauOld );
            DescentDir *= Tau;
            stepSizeFound = true;

          } else {
            // too short step
            TauOld = Tau;
            energyBackup = energyNew;
            Tau *= 2.;
          }
        }
      }
    } while ( !stepSizeFound );
   
   return Tau;
  }

  
protected:
  // auxiliary method for WOlfe condition
  RealType zoomWolfe( const VectorType &DescentDir,  VectorType &CurrentPosition,
                      RealType &energy, const RealType dE,
                      RealType ELo, RealType TauLo, RealType TauHi ) const {
    // either: TauLo satisfies the Armijo condition, but is too short for curvature condition, and TauHi > TauLo does not satisfy Armijo condition
    // or: both satisfy Armijo condition; TauLo is too long and TauHi too short for curvature condition
    RealType dENew, energyNew;
    VectorType positionNew( CurrentPosition.size() ), Derivative( CurrentPosition.size() );
    RealType Tau = 0.;
    
    bool stepSizeFound = false;
    int counter = 0; // for non-termination due to rounding errors, we should fix the maximum number of iterations
    do {
      Tau = ( TauLo + TauHi ) / 2.;
      positionNew = CurrentPosition;
      positionNew += Tau * DescentDir;
      energyNew = _stepSizeControlHelper.evaluateLinesearchFnc( positionNew );      

      if ( energyNew > std::min( energy + _sigma * Tau * dE, ELo ) )
        // Armijo condition violated...
        TauHi = Tau;
      else {
        // Armijo condition fulfilled...
        dENew = _stepSizeControlHelper.evaluateLinesearchGrad( positionNew, DescentDir );

        if ( std::abs( dENew ) <= - _beta * dE ){
          // strong Wolfe (Armijo + curvature) condition fulfilled...
          stepSizeFound = true;
        }
        else {
          // curvature condition violated...
          if ( dENew * ( TauHi - TauLo ) >= 0. )
            TauHi = TauLo;
          TauLo = Tau;
          ELo = energyNew;
        }
      }
    } while ( !stepSizeFound && (++counter < 30) );

    if ( counter >= 30 ) {
      Tau = 0.;
    }
    else {
      CurrentPosition = positionNew;
      energy = energyNew;
    }  
    
    return Tau;
  }
  
};




//===================================================================================================================================
// First Order Methods
//===================================================================================================================================


//!==========================================================================================================
template<typename DataTypeContainer >
class StepsizeControlHelperFirstOrder {
    
public:
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType    RealType;
  typedef typename DataTypeContainer::VectorType  VectorType;
  
  const aol::NonlinearEnergyOp<DataTypeContainer> &_energyOp;
    
public:
  StepsizeControlHelperFirstOrder( const aol::NonlinearEnergyOp<DataTypeContainer> &E ) : _energyOp(E) {}
  
    //Returns the scalar objective function evaluated at CurrentPosition.
  RealType evaluateLinesearchFnc( const VectorType &CurrentPosition ) const {
    RealType Energy;
    _energyOp.evaluateEnergy( CurrentPosition, Energy);
    return Energy;
  }

  //Returns the dot product of the energy derivative and the descent direction.
  RealType evaluateLinesearchGrad( const VectorType &Position, const VectorType &DescentDir ) const {
    VectorType tmp( DescentDir );
    _energyOp.evaluateJacobian( Position, tmp );
    return tmp.dot(DescentDir);
  } 
  
  //Returns the scalar objective function evaluated at CurrentPosition+DescentDir*timestepWidth.
  RealType evaluateLinesearchFnc( const VectorType &CurrentPosition, const VectorType &DescentDir, const RealType timestepWidth ) const {
      return evaluateLinesearchFnc( CurrentPosition + timestepWidth * DescentDir );
  }
  
};




//! \brief This class implements a quasi-Newton method with BFGS-Update of the matrix for minimizing a scalar function.
template <typename DataTypeContainer, template<class> class StepSizeControllerType >
class QuasiNewtonBFGS{

protected:
  typedef typename DataTypeContainer::RealType              RealType;
  typedef typename DataTypeContainer::VectorType            VectorType;
  typedef typename DataTypeContainer::SparseMatrixType      MatrixType;
     
  const aol::NonlinearEnergyOp<DataTypeContainer> &_energyOp;
  StepsizeControlHelperFirstOrder<DataTypeContainer> _stepSizeControlHelper;
  StepSizeControllerType< StepsizeControlHelperFirstOrder<DataTypeContainer> > _stepsizeControl;
  
  const int _maxIterations;
  const RealType _stopEpsilon;  
  const int _reset;
  bool _writeOutput;
  mutable std::vector<VectorType> _solutionTrajectory; 
  
  // BFGS stored information
  mutable int _counter;
  mutable std::vector<VectorType> _y;
  mutable std::vector<VectorType> _dx;
  mutable std::vector<RealType> _dx_y;
  mutable std::vector<VectorType> _auxMemory;  

public:
  QuasiNewtonBFGS ( const aol::NonlinearEnergyOp<DataTypeContainer> &EnergyOp,
                    const int MaxIterations = 50,
                    const RealType StopEpsilon = 1.e-6,
                    const int Reset = 50,                    
                    const bool writeOutput = true ) :  
    _energyOp(EnergyOp), _stepSizeControlHelper( EnergyOp ), _stepsizeControl( _stepSizeControlHelper ),  
    _maxIterations ( MaxIterations ),
    _stopEpsilon ( StopEpsilon ),
    _reset ( Reset ),
    _writeOutput ( writeOutput ),
    _counter(0),
    _dx_y( Reset ) { }


public:
    
  void solve( const VectorType &Arg, VectorType &Dest ) const {
      
    if( _writeOutput ){
      RealType energy;
      _energyOp.evaluateEnergy(Arg, energy );
      std::cout << "=======================================================================" << std::endl;
      std::cout << "Start BFGS with max " << _maxIterations << " iterations and eps = " << _stopEpsilon << "." << std::endl;
      std::cout << "Initial energy: " << energy << std::endl;
      std::cout << "=======================================================================" << std::endl; 
    }
    
    Dest = Arg;
    
    VectorType Dx ( Dest.size() );
    VectorType f ( Dest.size() );
    VectorType Df ( Dest.size() );
    VectorType tmp ( Dest.size() );
    VectorType descentDir ( Dest.size() );

    RealType FNorm = 1e+15;
    RealType tau = 1.0;
    int iterations = 0;
    bool forcedReset = false;

    while ( FNorm > _stopEpsilon && (iterations < _maxIterations) && tau > 0. ) {
      // Quasi-Newton-iteration given by x_{k+1} = x_k - tau_k B_k^{-1} f(x_k)
      // iteration number k (i.e. find x_{k+1} from current approximation x_k)
        
      // compute f and its norm
      _energyOp.evaluateJacobian( Dest, f );

      // save x_k and f(x_k) for the computation of Dx = x_{k+1}-x_k and Df = f(x_{k+1})-f(x_k)
      Dx = Dest;
      Dx *= -1.;
      Df = f;
      Df *= -1.;

      applyInverse( f, descentDir );
      descentDir *= -1.;

      // get tau and update position and descent direction
      tau = _stepsizeControl.getTimestepWidth( descentDir, Dest, tau );
                             
      if ( ( tau == 0 ) && ( forcedReset == false ) ) {
        reset();
        tau = 1.;
        forcedReset = true;
        continue;
      }
      forcedReset = false;

      _energyOp.evaluateJacobian(Dest, tmp);
      FNorm = tmp.norm();

      Dx += Dest;
      Df += tmp;

      // update of B, B_{k+1} = B_k - \frac{B_k Dx Dx^T B_k}{Dx\cdot (B_k Dx)}+\frac{Df Df^T}{Df \cdot Dx}
      update(Dx, Df);
      iterations++;
      
      if ( _writeOutput ){ //|| (10*iterations)%_maxIterations == 0
        RealType energy;
        _energyOp.evaluateEnergy(Dest, energy );
        std::cout << "step = " << iterations << " , \t stepsize = " << tau << ", \t \t energy = " << energy << ", \t error = " << FNorm << std::endl;
        _solutionTrajectory.push_back( Dest );
      }

    } // end while
    
    if( _writeOutput ){
        RealType energy;
        _energyOp.evaluateEnergy(Dest, energy );
        std::cout << "=======================================================================" << std::endl;
        std::cout << "Finished BFGS after " << iterations << " steps." << std::endl;
        std::cout << "Final stepsize = " << tau << ", energy = " << energy << ", error = " << FNorm << std::endl;
        std::cout << "=======================================================================" << std::endl;
    }
  }
  
  const std::vector<VectorType> & getSolutionTrajectory( ) const {return _solutionTrajectory;}
  
protected:    
  void update( const VectorType &DX, const VectorType &Y ) const {
    if ( _counter == _reset )
      reset();
    _counter++;
    
    if ( static_cast<int>(_y.size()) < _counter ) {
      _y.push_back ( Y );
      _dx.push_back( DX );
      _auxMemory.push_back( Y );
    } else {
      _y[_counter-1] = Y;
      _dx[_counter-1] = DX;
    }
    _dx_y[_counter-1] = DX.dot(Y);
  }

  void reset () const {
    _y.clear();
    _dx.clear();
    _counter = 0;
  }

  void apply ( const VectorType &Arg, VectorType &Dest ) const {
    Dest = Arg;
    for ( int i = 0; i < _counter; ++i ) {
      _auxMemory[i] = _dx[i]; // _auxMemory[k] = B_k \Delta dx_k
      for ( int j = 0; j < i; ++j ) {
        _auxMemory[i] += ( (_y[j]).dot(_dx[i]) / _dx_y[j] ) * _y[j];
        _auxMemory[i] -= (( (_auxMemory[j]).dot(_dx[i]) ) / ( _auxMemory[j] * _dx[j] ) ) *_auxMemory[j];
      }
      Dest += ( ( (_y[i]).dot(Arg) ) / _dx_y[i] ) * _y[i];
      Dest -= ( ( _auxMemory[i] * Arg ) / ( _auxMemory[i] * _dx[i] ) ) *_auxMemory[i];
    }
  }

  void applyInverse ( const VectorType &Arg, VectorType &Dest ) const {
    Dest = Arg;
    for ( int i = _counter - 1; i >= 0; --i ) {
      RealType dotProd = ( (_dx[i]).dot(Dest) ) / _dx_y[i];
      _auxMemory[i] = _dx[i];
      _auxMemory[i] *= dotProd; // _auxMemory[k] = ( dx_k dx_k^T / (y_k^T dx_k)) (I - y_{k+1} dx_{k+1}^T / (y_{k+1}^T dx_{k+1})) ... (I - y_K dx_K^T / (y_K^T dx_K)) Arg
      Dest -= dotProd * _y[i]; // Dest = (I - y_k dx_k^T / (y_k^T dx_k)) ... (I - y_K dx_K^T / (y_K^T dx_K)) Arg
    }
    for ( int i = 0; i < _counter; ++i ) {
      Dest -= (Dest.dot(_y[i])  / _dx_y[i]) * _dx[i]; // aux = (I - y_k dx_k^T / (y_k^T dx_k))^T H_{k-1} (I - y_k dx_k^T / (y_k^T dx_k)) ... (I - y_K dx_K^T / (y_K^T dx_K)) Arg
      Dest += _auxMemory[i];
    }
  }

};






//===================================================================================================================================
// Second Order Methods
//===================================================================================================================================

//!==========================================================================================================
//! Class for Newton stepsize control
template<typename StepsizeControlHelper>
class NewtonOptimalStepsizeControl{
    
protected:
  typedef typename StepsizeControlHelper::DTContainer   DataTypeContainer;
  typedef typename DataTypeContainer::RealType          RealType;
  typedef typename DataTypeContainer::VectorType        VectorType;
  
  const StepsizeControlHelper _stepSizeControlHelper;

  const RealType _startTau;
  const RealType _tauMin, _tauMax;
    
public:
  NewtonOptimalStepsizeControl( const StepsizeControlHelper & stepSizeControlHelper, const RealType startTau = 1.0, const RealType tauMin = 1.e-6, const RealType tauMax = 4. ) 
  : _stepSizeControlHelper ( stepSizeControlHelper ), _startTau(startTau), _tauMin(tauMin), _tauMax(tauMax){}
  
  ~NewtonOptimalStepsizeControl( ){}

  // Calculates the step size based on "Schaback, Werner - Numerische Mathematik, 4. Auflage, Seite 129". Is global convergent, as long as F fulfils certain regularity conditions.
  RealType getTimestepWidth( VectorType &DescentDir, VectorType &CurrentPosition, const RealType /*OldTau*/ ) const {
      
    const RealType DescentDirNormSqr = DescentDir.squaredNorm();
    
    if( DescentDirNormSqr > 0. ){
      // initial guess for L: L = ||F(y)-F(x)-F'(x)(y-x)|| / ||y-x||^2, where x = CurrentPosition, y = newPosition = x + DescentDir
      VectorType newPosition = CurrentPosition + DescentDir;

      VectorType pTmpFx (DescentDir.size()), pTmpFy (DescentDir.size());
      _stepSizeControlHelper.getEnergyOp( ).evaluateJacobian(CurrentPosition, pTmpFx);
      _stepSizeControlHelper.getEnergyOp( ).evaluateJacobian(newPosition, pTmpFy);
      const RealType fNorm = pTmpFx.norm();
      VectorType pTmpDFx (DescentDir.size());
      pTmpDFx =  _stepSizeControlHelper.getHessian() * DescentDir;
      VectorType pTmp (DescentDir.size());
      pTmp = pTmpFy - pTmpFx - pTmpDFx;
      //! \todo check sign
//       pTmp += pTmp2;
//       pTmp -= pTmp2;
      RealType L = pTmp.norm()/DescentDirNormSqr;

      // If F is locally linear and matches the gradient with numerical perfection take the full step
      if ( L < 1e-15 ){
          CurrentPosition += DescentDir;
          return 1.0;
      }

      RealType tau = std::min( fNorm/(2.*L*DescentDirNormSqr), 1.0 );
      RealType temp1 = fNorm - std::sqrt(2.* _stepSizeControlHelper.evaluateLinesearchFnc( CurrentPosition, DescentDir, tau ));
      RealType temp2 = tau * (fNorm - L*tau*DescentDirNormSqr);
      RealType tauCandidate;
      if( temp1 < temp2 ){
        do{
          L *= 2.;
          tau = std::min( fNorm/(2.*L*DescentDirNormSqr), 1.0 );
          temp1 = fNorm - sqrt(2.* _stepSizeControlHelper.evaluateLinesearchFnc( CurrentPosition, DescentDir, tau ));
          temp2 = tau*(fNorm - L*tau*DescentDirNormSqr);
          // Prevent the while loop from getting stuck if DescentDir is not a descent direction.
          if ( tau < this->_tauMin ){
            temp1 = temp2; tau = 0.;
          }
        } while( temp1 < temp2 );
      }
      else{
        do{
          // Compute tauCandidate, temp1 and temp2 with L/2 instead of L to find out if a bigger step size is still feasible.
          tauCandidate = std::min( fNorm/(L*DescentDirNormSqr), 1.0 );
          if(tauCandidate == 1.0)
            break;
          temp1 = fNorm - std::sqrt( 2.* _stepSizeControlHelper.evaluateLinesearchFnc( CurrentPosition, DescentDir, tauCandidate ) );
          temp2 = tauCandidate * (fNorm - 0.5*L*tauCandidate*DescentDirNormSqr);

          if ( temp1 >= temp2 ){
            L /= 2.;
            tau = tauCandidate;
          }
        } while( temp1 >= temp2 );
      }
      
      // update
      DescentDir *= tau;
      CurrentPosition += DescentDir;
      return tau;
    }
    else
      return 0.;
  }
  
};



template<typename DataTypeContainer>
class StepsizeControlHelperSecondOrder {   
public:
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType    RealType;
  typedef typename DataTypeContainer::VectorType  VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  
  const aol::NonlinearEnergyOp<DataTypeContainer> &_energyOp;
  const SparseMatrixType &_Hessian;
    
public:
  StepsizeControlHelperSecondOrder( const aol::NonlinearEnergyOp<DataTypeContainer> &E, const SparseMatrixType &Hessian ) : _energyOp(E), _Hessian ( Hessian ) {}
  
  //Returns the scalar objective function evaluated at CurrentPosition.
  RealType evaluateLinesearchFnc( const VectorType &Position ) const {
    VectorType temp ( Position.size() );
    _energyOp.evaluateJacobian( Position, temp );
    return 0.5 * temp.squaredNorm();
  }

  //Returns the dot product of the energy derivative and the descent direction (if _Hessian is filled with D2E(position) )
  RealType evaluateLinesearchGrad( const VectorType &Position, const VectorType &DescentDir ) const {
    VectorType outer = _Hessian * DescentDir;
    VectorType grad (DescentDir.size());
    _energyOp.evaluateJacobian(Position, grad);
    return outer.dot(grad);
  }
  
  //Returns the scalar objective function evaluated at CurrentPosition+DescentDir*timestepWidth.
  RealType evaluateLinesearchFnc( const VectorType &CurrentPosition, const VectorType &DescentDir, const RealType timestepWidth ) const {
      return evaluateLinesearchFnc( CurrentPosition + timestepWidth * DescentDir );
  }
  
  const aol::NonlinearEnergyOp<DataTypeContainer> & getEnergyOp( ) const { return _energyOp; }
  const SparseMatrixType & getHessian( ) const { return _Hessian; };
  
};





//!==========================================================================================================
//! Newton method to find a root of of F: \R^n -> \R^n, which is the Jacobian of an energy E :\R^n -> \R
template <typename DataTypeContainer, template<class> class StepSizeControllerType >
class NewtonMethod{
    
protected:
  typedef typename DataTypeContainer::RealType              RealType;
  typedef typename DataTypeContainer::VectorType            VectorType;
  typedef typename DataTypeContainer::SparseMatrixType      MatrixType;
  
  const aol::NonlinearEnergyOp<DataTypeContainer> &_energyOp; 
  const int _numDofs;
  mutable MatrixType _Hessian;
  StepsizeControlHelperSecondOrder<DataTypeContainer> _stepSizeControlHelper;
  StepSizeControllerType< StepsizeControlHelperSecondOrder<DataTypeContainer> > _stepsizeControl;
  
  LinearSolverMethod _solverType;
  const int _maxIterations;
  const RealType _tolerance, _acceptTolerance;
  const int _outputLevel;
  
  mutable std::vector<VectorType> _solutionTrajectory;

public:
  NewtonMethod( const aol::NonlinearEnergyOp<DataTypeContainer> &E, const int numDofs,
                const int MaxIterations = 50, 
                const RealType tolerance = 1.e-8,
                const int solverType = 1, const int outputLevel = 2 ) 
  : _energyOp(E), _numDofs( numDofs ), _Hessian( numDofs, numDofs ), _stepSizeControlHelper( E, _Hessian ), _stepsizeControl( _stepSizeControlHelper ), 
    _solverType( static_cast<LinearSolverMethod> ( solverType ) ), _maxIterations(MaxIterations), _tolerance(tolerance), _acceptTolerance( numDofs * tolerance ), _outputLevel ( outputLevel ) {}

  RealType computeErrorNorm ( const VectorType& delta_x_k ) const { return delta_x_k.norm();}
 
  // x^{k+1} = x^k + tau * d^k, where d^k solves D^2E[x^k] d^k = - DE[x^k]
  void solve( const VectorType &x_0, VectorType &x_k ) const {
   
    x_k = x_0;
    VectorType F_x_k ( _numDofs ), delta_x_k ( _numDofs );
    _energyOp.evaluateJacobian(x_k, F_x_k);
    
    RealType tau = 1.0;    
    RealType FNorm = computeErrorNorm ( F_x_k );
    
    if( _outputLevel > 2 ){
      std::cout << "=======================================================================" << std::endl;
      std::cout << "Start Newton method with " << _maxIterations << " iterations and tolerance = " << _tolerance << "." << std::endl;
      std::cout << "=======================================================================" << std::endl;
    }
    
    int iterations = 0;
    while ( FNorm > _tolerance && (iterations < _maxIterations) && tau > 0. ) {
      iterations++;

      // Newton iteration given by x^{k+1} = x^k - tau D2F(x^k)^{-1}(DF(x^k))
      _energyOp.evaluateHessian( x_k, _Hessian);
      F_x_k *= -1.;
      
      //LinearSolver<ConfiguratorType>(_solverType).solve( _Hessian, F_x_k, delta_x_k );
      solveLinearSystem<DataTypeContainer> ( _Hessian, delta_x_k, F_x_k, _solverType );
      
      // get tau and update position and descent direction
      tau = _stepsizeControl.getTimestepWidth( delta_x_k, x_k, tau );

      if( tau > 0. ){
          _energyOp.evaluateJacobian( x_k, F_x_k );
          FNorm = computeErrorNorm ( F_x_k );
      }

      if( _outputLevel > 2 ) {
          std::cout << "step = " << iterations << " , stepsize = " << tau << " , norm = " << FNorm << std::endl;
          _solutionTrajectory.push_back( x_k );
      }
    } // end while
    
    if( _outputLevel > 1 ){
      cout << endl 
           << "=======================================================================" << endl
           << "Finished Newton method after " << iterations << " steps." << std::endl
           << "Final stepsize = " << tau << ", error = " << FNorm << std::endl
           << "=======================================================================" << endl << endl;
    }
    if( _outputLevel == 1 ){
         if( FNorm > _acceptTolerance ) {
            cout << endl 
           << "=======================================================================" << endl
           << "Finished Newton method after " << iterations << " steps without satisfying tolerance." << std::endl
           << "Final stepsize = " << tau << ", error = " << FNorm << std::endl
           << "=======================================================================" << endl << endl;
         }
    }

  }
  
  const std::vector<VectorType> & getSolutionTrajectory( ) const {return _solutionTrajectory;}
  
};




//=============================================================================
#endif // LINESEARCHMETHODS_HH
//=============================================================================
