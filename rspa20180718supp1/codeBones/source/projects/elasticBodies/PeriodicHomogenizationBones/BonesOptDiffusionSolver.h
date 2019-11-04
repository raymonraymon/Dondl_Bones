#ifndef __PERIODICHOMOGENIZATIONBONESMULTIPLELOADOPTIMALDIFFUSIONSOLVERBONES_H
#define __PERIODICHOMOGENIZATIONBONESMULTIPLELOADOPTIMALDIFFUSIONSOLVERBONES_H

#include <general.h>
#include <loadAndSave.h>
#include <linearSolver.h>
#include <LinearSystemSolver.h>
#include <quocHandler.h>

#include "BonesLinDiffusionEnergies.h"


using namespace quocFE;

namespace shapeOptBonePolymerPeriodicHomogenization{

int counterOptimalDiffusionSolverInterfaceMultipleLoad = 0;
    
template <typename MatOptConfigurator>
class OptimalDiffusionSolverInterfaceMultipleLoad {
protected :
  typedef typename MatOptConfigurator::ConfiguratorType        ConfiguratorType;
  typedef typename ConfiguratorType::DTContainer               DataTypeContainer;
  typedef typename ConfiguratorType::RealType                  RealType;
  typedef typename ConfiguratorType::MaskType                  MaskType;
  typedef typename ConfiguratorType::InitType                  MeshType;
  typedef typename ConfiguratorType::DomVecType                DomVecType;
  typedef typename ConfiguratorType::VectorType                VectorType;
  typedef typename ConfiguratorType::SparseMatrixType          SparseMatrixType;
  typedef typename DataTypeContainer::ParameterParserType      ParameterParserType;
  
  const ParameterParserType &_parser;
  const MatOptConfigurator &_matOpConf;
  const ConfiguratorType & _conf;
  const QuocHandler<ConfiguratorType> & _quocHandler;
  const int _dimDomain;

  mutable VectorType _pfCollabsed, _pfExtended;
  const MaskType & _mask;
  const std::vector<int> & _periodicIndices;
  
  const std::vector<DomVecType> &_affineDiffusionBone, &_affineDiffusionPolymer;
  const int _numLoadsBone, _numLoadsPolymer;
  mutable std::vector<VectorType> _solDiffusionBonePeriodic, _solDiffusionBonePeriodicallyExtended, _solDiffusionBoneAndMultiplier, 
                                  _solDiffusionPolymerPeriodic, _solDiffusionPolymerPeriodicallyExtended, _solDiffusionPolymerAndMultiplier;
  
  mutable DirectLinearSystemSolver<DataTypeContainer> _directLinearSystemSolverBone, _directLinearSystemSolverPolymer;
  mutable IterativeLinearSystemSolver<DataTypeContainer> _iterativeLinearSystemSolverBone, _iterativeLinearSystemSolverPolymer;

  public:

   OptimalDiffusionSolverInterfaceMultipleLoad( const ParameterParserType &Parser,
                                        const MatOptConfigurator & matOpConf,
                                        const VectorType & Phasefield,
                                        const QuocHandler<ConfiguratorType> & quocHandler,
                                        const std::vector<DomVecType> &affineDiffusionBone,
                                        const std::vector<DomVecType> &affineDiffusionPolymer ) :
         _parser ( Parser ),
         _matOpConf ( matOpConf ),
         _conf ( matOpConf._conf ),
         _quocHandler ( quocHandler ),
         _dimDomain( _conf.dimDomain ),
         _pfCollabsed ( Phasefield ),
         _pfExtended ( Phasefield ),
         _mask ( quocHandler.getPeriodicMask() ),
         _periodicIndices ( quocHandler.getPeriodicIndices() ),
         _affineDiffusionBone ( affineDiffusionBone ), _affineDiffusionPolymer ( affineDiffusionPolymer ),
         _numLoadsBone( _affineDiffusionBone.size() ), _numLoadsPolymer( _affineDiffusionPolymer.size() ),
         _directLinearSystemSolverBone( "opt diffusion bone", _parser.template get<string> ("saving.saveDirectory" ).c_str() ), 
         _directLinearSystemSolverPolymer( "opt diffusion polymer", _parser.template get<string> ("saving.saveDirectory" ).c_str() ), 
         _iterativeLinearSystemSolverBone ( "opt diffusion bone", _parser.template get<string> ("saving.saveDirectory" ).c_str() ),
         _iterativeLinearSystemSolverPolymer ( "opt diffusion polymer", _parser.template get<string> ("saving.saveDirectory" ).c_str() )
    {
         _quocHandler.extendVectorPeriodically( _pfExtended );
             
         for( int i=0; i<_numLoadsBone; ++i ){
             _solDiffusionBonePeriodic.push_back ( VectorType( _conf.getNumGlobalDofs() ) );
             _solDiffusionBonePeriodicallyExtended.push_back ( VectorType( _conf.getNumGlobalDofs() ) );
             _solDiffusionBoneAndMultiplier.push_back ( VectorType( _conf.getNumGlobalDofs() + 1 ) );
             
             //! \note: this is necessery if we want to use solveWithGuess 
             _solDiffusionBoneAndMultiplier[i].setZero();
         }
         
        for( int i=0; i<_numLoadsPolymer; ++i ){;
             _solDiffusionPolymerPeriodic.push_back ( VectorType( _conf.getNumGlobalDofs() ) );
             _solDiffusionPolymerPeriodicallyExtended.push_back ( VectorType( _conf.getNumGlobalDofs() ) );
             _solDiffusionPolymerAndMultiplier.push_back ( VectorType( _conf.getNumGlobalDofs() + 1 ) );
             
             //! \note: this is necessery if we want to use solveWithGuess 
             _solDiffusionPolymerAndMultiplier[i].setZero(); 
         }
         
    }
         
    const ParameterParserType& getParser( ) const { return _parser; }
    const MatOptConfigurator & getMatOptConfigurator() const { return _matOpConf;}
    const QuocHandler<ConfiguratorType> & getQuocHandler ( ) const { return _quocHandler; }

    const VectorType& getPhaseFieldPeriodicallyExtended( ) const { return _pfExtended; }
    
    
    template<MaterialTypeBonePolymer MaterialType>
    DirectLinearSystemSolver<DataTypeContainer>& getDirectLinearSolver( ) const{ 
     switch( MaterialType ){
        case BONE :    return _directLinearSystemSolverBone;    break;
        case POLYMER : return _directLinearSystemSolverPolymer; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    IterativeLinearSystemSolver<DataTypeContainer>& getIterativeLinearSolver( ) const{ 
     switch( MaterialType ){
        case BONE :    return _iterativeLinearSystemSolverBone;    break;
        case POLYMER : return _iterativeLinearSystemSolverPolymer; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    const VectorType& getSolutionDiffusionPeriodic(int i) const { 
     switch( MaterialType ){
        case BONE :    return _solDiffusionBonePeriodic[i];    break;
        case POLYMER : return _solDiffusionPolymerPeriodic[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    VectorType& SolutionDiffusionPeriodic(int i) const { 
     switch( MaterialType ){
        case BONE :    return _solDiffusionBonePeriodic[i];    break;
        case POLYMER : return _solDiffusionPolymerPeriodic[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    void setSolutionDiffusionPeriodic ( int i, const VectorType &solution ) const{ 
     switch( MaterialType ){
        case BONE :    _solDiffusionBonePeriodic[i] = solution;    break;
        case POLYMER : _solDiffusionPolymerPeriodic[i] = solution; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    const VectorType& getSolutionDiffusionPeriodicallyExtended(int i) const { 
    switch( MaterialType ){
        case BONE :    return _solDiffusionBonePeriodicallyExtended[i];    break;
        case POLYMER : return _solDiffusionPolymerPeriodicallyExtended[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    VectorType& SolutionDiffusionPeriodicallyExtended(int i) const { 
    switch( MaterialType ){
        case BONE :    return _solDiffusionBonePeriodicallyExtended[i];    break;
        case POLYMER : return _solDiffusionPolymerPeriodicallyExtended[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    const VectorType& getSolutionDiffusionAndMultiplier(int i) const { 
    switch( MaterialType ){
        case BONE :    return _solDiffusionBoneAndMultiplier[i];    break;
        case POLYMER : return _solDiffusionPolymerAndMultiplier[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    
    template<MaterialTypeBonePolymer MaterialType>
    VectorType& SolutionDiffusionAndMultiplier(int i) const{ 
    switch( MaterialType ){
        case BONE :    return _solDiffusionBoneAndMultiplier[i];    break;
        case POLYMER : return _solDiffusionPolymerAndMultiplier[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    const DomVecType& getAffineDiffusion(int i) const { 
    switch( MaterialType ){
        case BONE :    return _affineDiffusionBone[i];    break;
        case POLYMER : return _affineDiffusionPolymer[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    const int getNumDiffusionDofs( ) const {return _conf.getNumGlobalDofs(); }
    const int getNumPeriodicDiffusionDofs( ) const {return _conf.getNumGlobalDofs(); }
    
    template<MaterialTypeBonePolymer MaterialType>
    const int getNumLoads() const { 
    switch( MaterialType ){
        case BONE :    return _numLoadsBone;    break;
        case POLYMER : return _numLoadsPolymer; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }

};



template <typename MatOptConfigurator >
class OptimalDiffusionSolverMultipleLoad : public OptimalDiffusionSolverInterfaceMultipleLoad<MatOptConfigurator> {

  typedef typename MatOptConfigurator::ConfiguratorType        ConfiguratorType;
  typedef typename ConfiguratorType::DTContainer               DataTypeContainer;
  typedef typename ConfiguratorType::RealType                  RealType;
  typedef typename ConfiguratorType::MaskType                  MaskType;
  typedef typename ConfiguratorType::InitType                  MeshType;
  typedef typename ConfiguratorType::DomVecType                DomVecType;
  typedef typename ConfiguratorType::VectorType                VectorType;
  typedef typename ConfiguratorType::SparseMatrixType          SparseMatrixType;
  typedef typename DataTypeContainer::ParameterParserType      ParameterParserType;
  typedef typename DataTypeContainer::TripletType              TripletType;
 
protected :
  mutable SparseMatrixType _HessianLinDiffusionBone, _HessianLinDiffusionPolymer;
  mutable SparseMatrixType _HessianLinDiffusionMixedBone, _HessianLinDiffusionMixedPolymer;
  mutable SparseMatrixType _HessianLinDiffusionAffineBone, _HessianLinDiffusionAffinePolymer;
  mutable std::vector<VectorType> _rhsBone, _rhsPolymer;
  
  public:
      OptimalDiffusionSolverMultipleLoad( const ParameterParserType &Parser,  
                                       const MatOptConfigurator & matOpConf,
                                       const VectorType & Phasefield,
                                       const QuocHandler<ConfiguratorType> & quocHandler,
                                       const std::vector<DomVecType> &affineDiffusionBone, 
                                       const std::vector<DomVecType> &affineDiffusionPolymer
                                     ) :
         OptimalDiffusionSolverInterfaceMultipleLoad<MatOptConfigurator> ( Parser, matOpConf, Phasefield, quocHandler, affineDiffusionBone, affineDiffusionPolymer ),
         _HessianLinDiffusionBone( ( matOpConf._conf.getNumGlobalDofs() + 1) , (matOpConf._conf.getNumGlobalDofs() + 1)  ),
         _HessianLinDiffusionPolymer( ( matOpConf._conf.getNumGlobalDofs() + 1) , (matOpConf._conf.getNumGlobalDofs() + 1)  ),
         _HessianLinDiffusionMixedBone( matOpConf._conf.getNumGlobalDofs(), this->_dimDomain ),
         _HessianLinDiffusionMixedPolymer( matOpConf._conf.getNumGlobalDofs(), this->_dimDomain ),
         _HessianLinDiffusionAffineBone( this->_dimDomain, this->_dimDomain ),
         _HessianLinDiffusionAffinePolymer( this->_dimDomain, this->_dimDomain )
          
        {
            for( int i=0; i<this->template getNumLoads<BONE>(); ++i ){
                _rhsBone.push_back ( VectorType( ( matOpConf._conf.getNumGlobalDofs() + 1 ) ) );
            }
            for( int i=0; i<this->template getNumLoads<POLYMER>(); ++i ){
                _rhsPolymer.push_back ( VectorType (( matOpConf._conf.getNumGlobalDofs() + 1 ) ) );
            }
            this->assembleLinDiffusionHessian();
            
            if( this->getParser().template get<bool>("ConstraintProblem.solveWithDirectSolver") ){  
                this->template getDirectLinearSolver<BONE>().analyzePattern( this->getHessianLinDiffusion<BONE> () );
                this->template getDirectLinearSolver<POLYMER>().analyzePattern( this->getHessianLinDiffusion<POLYMER> () );
            }
            
            solve();
        }

      
    void assembleLinDiffusionHessian( ) const {
        #ifdef USE_OPENMP
        omp_set_nested(1);    
        #pragma omp parallel sections
        #endif
        {
            #ifdef USE_OPENMP 
            #pragma omp section
            #endif
            { 
                LinDiffusionHessian<MatOptConfigurator,BONE> ( this->_matOpConf, this->_pfExtended ).assemblePeriodic( _HessianLinDiffusionBone, this->_mask, this->_periodicIndices );
                _HessianLinDiffusionBone.makeCompressed();
            }
            #ifdef USE_OPENMP
            #pragma omp section
            #endif
            { 
                LinDiffusionHessian<MatOptConfigurator,POLYMER> ( this->_matOpConf, this->_pfExtended ).assemblePeriodic( _HessianLinDiffusionPolymer, this->_mask, this->_periodicIndices );
                _HessianLinDiffusionPolymer.makeCompressed();
            }
            #ifdef USE_OPENMP
            #pragma omp section
            #endif
            { 
                LinDiffusionHessian<MatOptConfigurator,BONE> ( this->_matOpConf, this->_pfExtended ).assembleMixedPeriodic( _HessianLinDiffusionMixedBone,this->_mask, this->_periodicIndices );
                _HessianLinDiffusionMixedBone.makeCompressed();
            }
            #ifdef USE_OPENMP
            #pragma omp section
            #endif
            { 
                 LinDiffusionHessian<MatOptConfigurator,POLYMER> ( this->_matOpConf, this->_pfExtended ).assembleMixedPeriodic( _HessianLinDiffusionMixedPolymer, this->_mask, this->_periodicIndices );
                 _HessianLinDiffusionMixedPolymer.makeCompressed();
            }
            #ifdef USE_OPENMP
            #pragma omp section
            #endif
            { 
                LinDiffusionHessian<MatOptConfigurator,BONE> ( this->_matOpConf, this->_pfExtended ).assembleAffine( _HessianLinDiffusionAffineBone );
                _HessianLinDiffusionAffineBone.makeCompressed();
            }
            #ifdef USE_OPENMP
            #pragma omp section
            #endif
            { 
                LinDiffusionHessian<MatOptConfigurator,POLYMER> ( this->_matOpConf, this->_pfExtended ).assembleAffine( _HessianLinDiffusionAffinePolymer );
                _HessianLinDiffusionAffinePolymer.makeCompressed();
            }
        }
        for( int i=0; i<this->_numLoadsBone; ++i ){
            VectorType RHSBONE = _HessianLinDiffusionMixedBone * this->template getAffineDiffusion<BONE> ( i );
            for( int j=0; j < RHSBONE.size(); ++j ) (_rhsBone[i])[j] = - RHSBONE[j];
            (_rhsBone[i])[this->getNumDiffusionDofs()] = 0.0;
        }
        
        for( int i=0; i<this->_numLoadsPolymer; ++i ){
            VectorType RHSPOLYMER = _HessianLinDiffusionMixedPolymer * this->template getAffineDiffusion<POLYMER> ( i );
            for( int j=0; j < RHSPOLYMER.size(); ++j ) (_rhsPolymer[i])[j] = - RHSPOLYMER[j];
            (_rhsPolymer[i])[this->getNumDiffusionDofs()] = 0.0;
        }
      }
        
        
    template< MaterialTypeBonePolymer MaterialType>
    const SparseMatrixType& getHessianLinDiffusion() const { 
      switch( MaterialType ){
        case BONE :
            return _HessianLinDiffusionBone;
            break;
        case POLYMER :
            return _HessianLinDiffusionPolymer;
            break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      } 
    }
    
    
    template< MaterialTypeBonePolymer MaterialType>
    const SparseMatrixType& getHessianLinDiffusionMixedPart() const { 
      switch( MaterialType ){
        case BONE :
            return _HessianLinDiffusionMixedBone;
            break;
        case POLYMER :
            return _HessianLinDiffusionMixedPolymer;
            break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      } 
    }
    
    template< MaterialTypeBonePolymer MaterialType>
    const SparseMatrixType& getHessianLinDiffusionAffinePart() const { 
      switch( MaterialType ){
        case BONE :
            return _HessianLinDiffusionAffineBone;
            break;
        case POLYMER :
            return _HessianLinDiffusionAffinePolymer;
            break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      } 
    }
    
    template< MaterialTypeBonePolymer MaterialType>
    const VectorType& getRHS(int i) const { 
      switch( MaterialType ){
        case BONE : return _rhsBone[i]; break;
        case POLYMER : return _rhsPolymer[i]; break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      } 
    }
    
    
    template< MaterialTypeBonePolymer MaterialType>
    const SparseMatrixType& getSystemMatForAdjointProblem () const { return getHessianLinDiffusion<MaterialType> (); }

    
    bool updatePhasefield( const VectorType &pf ) const {
        bool NewPhasefieldDiffersFromOld = false;
      VectorType pfCollabsed( pf );
      this->_quocHandler.collabseVectorPeriodically( pfCollabsed );
      RealType diff =( this->_pfCollabsed - pfCollabsed ).squaredNorm();
      if ( diff > 1.e-15 ){
        NewPhasefieldDiffersFromOld = true;
        this->_pfCollabsed = pfCollabsed;
        this->_pfExtended = this->_pfCollabsed;
        this->_quocHandler.extendVectorPeriodically( this->_pfExtended );
        this->assembleLinDiffusionHessian();
        solve();
      }
      return NewPhasefieldDiffersFromOld;
    }

protected:
    
  template<MaterialTypeBonePolymer MaterialType>
  void solve( ) const {
      if( this->getParser().template get<bool>("ConstraintProblem.solveWithDirectSolver") ){  
        this->template getDirectLinearSolver<MaterialType>().prepareSolver( this->getHessianLinDiffusion<MaterialType> () );
      }
      
      for( int i=0; i< this->template getNumLoads<MaterialType> (); ++i ){
          
        if( this->getParser().template get<bool>("ConstraintProblem.solveWithDirectSolver") ){  
            this->template getDirectLinearSolver<MaterialType>().solve( this->template SolutionDiffusionAndMultiplier<MaterialType>(i), this->getRHS<MaterialType>(i) ); 
        }else{
            this->template getIterativeLinearSolver<MaterialType>().solve( this->getHessianLinDiffusion<MaterialType> (), this->template SolutionDiffusionAndMultiplier<MaterialType>(i), this->getRHS<MaterialType>(i), this->getParser().template get<RealType> ("ConstraintProblem.toleranceLinearSystem"), this->getParser().template get<RealType>("ConstraintProblem.maxItersFacLinearSystem") ); 
        }
        
        VectorType solDiffusionPeriodic ( this->template getSolutionDiffusionAndMultiplier<MaterialType>(i).segment( 0, this->_conf.getNumGlobalDofs() ) );
        this->template setSolutionDiffusionPeriodic<MaterialType> ( i, solDiffusionPeriodic );
        
        this->_quocHandler.extendVectorPeriodically( this->template getSolutionDiffusionPeriodic<MaterialType>(i), this->template SolutionDiffusionPeriodicallyExtended<MaterialType>(i) );
        
      }
  }
  
  
  void solve() const {
        this->solve<BONE>();
        this->solve<POLYMER>();
  }
    
  
};

}//end namespace

#endif
