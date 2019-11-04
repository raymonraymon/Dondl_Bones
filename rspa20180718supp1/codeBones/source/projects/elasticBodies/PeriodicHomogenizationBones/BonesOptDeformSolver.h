#ifndef __PERIODICHOMOGENIZATIONBONESMULTIPLELOADOPTIMALDEFORMSOLVERBONES_H
#define __PERIODICHOMOGENIZATIONBONESMULTIPLELOADOPTIMALDEFORMSOLVERBONES_H

#include <general.h>
#include <loadAndSave.h>
#include <linearSolver.h>
#include <LinearSystemSolver.h>
#include <quocHandler.h>

#include "BonesLinElastEnergies.h"


using namespace quocFE;

namespace shapeOptBonePolymerPeriodicHomogenization{

int counterOptimalDeformSolverInterfaceMultipleLoad = 0;
    
template <typename MatOptConfigurator>
class OptimalDeformSolverInterfaceMultipleLoad {
protected :
  typedef typename MatOptConfigurator::ConfiguratorType        ConfiguratorType;
  typedef typename ConfiguratorType::DTContainer               DataTypeContainer;
  typedef typename ConfiguratorType::RealType                  RealType;
  typedef typename ConfiguratorType::MaskType                  MaskType;
  typedef typename ConfiguratorType::InitType                  MeshType;
  typedef typename ConfiguratorType::VectorType                VectorType;
  typedef typename ConfiguratorType::SparseMatrixType          SparseMatrixType;
  typedef typename DataTypeContainer::ParameterParserType      ParameterParserType;
  
  const ParameterParserType &_parser;
  const MatOptConfigurator &_matOpConf;
  const ConfiguratorType & _conf;
  const QuocHandler<ConfiguratorType> & _quocHandler;
  const int _dimDomain; const int _numAffineSymGradDofs;

  mutable VectorType _pfCollabsed, _pfExtended;
  const MaskType & _mask;
  const std::vector<int> & _periodicIndices;
  
  const std::vector<VectorType> &_affineDispBone, &_affineDispPolymer;
  const int _numLoads;
  mutable std::vector<VectorType> _solDispBonePeriodic, _solDispBonePeriodicallyExtended, _solDispBoneAndMultiplier, _solDispPolymerPeriodic, _solDispPolymerPeriodicallyExtended, _solDispPolymerAndMultiplier;
  
  mutable DirectLinearSystemSolver<DataTypeContainer> _directLinearSystemSolverBone, _directLinearSystemSolverPolymer;
  mutable IterativeLinearSystemSolver<DataTypeContainer> _iterativeLinearSystemSolverBone, _iterativeLinearSystemSolverPolymer;

  public:

   OptimalDeformSolverInterfaceMultipleLoad( const ParameterParserType &Parser,
                                        const MatOptConfigurator & matOpConf,
                                        const VectorType & Phasefield,
                                        const QuocHandler<ConfiguratorType> & quocHandler,
                                        const std::vector<VectorType> &affineDispBone,
                                        const std::vector<VectorType> &affineDispPolymer ) :
         _parser ( Parser ),
         _matOpConf ( matOpConf ),
         _conf ( matOpConf._conf ),
         _quocHandler ( quocHandler ),
         _dimDomain( _conf.dimDomain ),
         _numAffineSymGradDofs ( _conf.numAffineSymGradDofs ),
         _pfCollabsed ( Phasefield ),
         _pfExtended ( Phasefield ),
         _mask ( quocHandler.getPeriodicMask() ),
         _periodicIndices ( quocHandler.getPeriodicIndices() ),
         _affineDispBone ( affineDispBone ), _affineDispPolymer ( affineDispPolymer ),
         _numLoads( _affineDispBone.size() ),
         _directLinearSystemSolverBone( "opt disp bone", _parser.template get<string> ("saving.saveDirectory" ).c_str() ), 
         _directLinearSystemSolverPolymer( "opt disp polymer", _parser.template get<string> ("saving.saveDirectory" ).c_str() ), 
         _iterativeLinearSystemSolverBone ( "opt disp bone", _parser.template get<string> ("saving.saveDirectory" ).c_str() ),
         _iterativeLinearSystemSolverPolymer ( "opt disp polymer", _parser.template get<string> ("saving.saveDirectory" ).c_str() )
    {
         _quocHandler.extendVectorPeriodically( _pfExtended );
             
         for( int i=0; i<_numLoads; ++i ){
             
             _solDispBonePeriodic.push_back ( VectorType( _dimDomain * _conf.getNumGlobalDofs() ) );
             _solDispBonePeriodicallyExtended.push_back ( VectorType( _dimDomain * _conf.getNumGlobalDofs() ) );
             _solDispBoneAndMultiplier.push_back ( VectorType( _dimDomain * (_conf.getNumGlobalDofs() + 1) ) );
             
             _solDispPolymerPeriodic.push_back ( VectorType( _dimDomain * _conf.getNumGlobalDofs() ) );
             _solDispPolymerPeriodicallyExtended.push_back ( VectorType( _dimDomain * _conf.getNumGlobalDofs() ) );
             _solDispPolymerAndMultiplier.push_back ( VectorType( _dimDomain * (_conf.getNumGlobalDofs() + 1) ) );
             
             //! \note: this is necessery if we want to use solveWithGuess 
             _solDispBoneAndMultiplier[i].setZero(); _solDispPolymerAndMultiplier[i].setZero(); 
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
    const VectorType& getSolutionDisplacementPeriodic(int i) const { 
     switch( MaterialType ){
        case BONE :    return _solDispBonePeriodic[i];    break;
        case POLYMER : return _solDispPolymerPeriodic[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    VectorType& SolutionDisplacementPeriodic(int i) const { 
     switch( MaterialType ){
        case BONE :    return _solDispBonePeriodic[i];    break;
        case POLYMER : return _solDispPolymerPeriodic[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    void setSolutionDisplacementPeriodic ( int i, const VectorType &solution ) const{ 
     switch( MaterialType ){
        case BONE :    _solDispBonePeriodic[i] = solution;    break;
        case POLYMER : _solDispPolymerPeriodic[i] = solution; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    const VectorType& getSolutionDisplacementPeriodicallyExtended(int i) const { 
    switch( MaterialType ){
        case BONE :    return _solDispBonePeriodicallyExtended[i];    break;
        case POLYMER : return _solDispPolymerPeriodicallyExtended[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    VectorType& SolutionDisplacementPeriodicallyExtended(int i) const { 
    switch( MaterialType ){
        case BONE :    return _solDispBonePeriodicallyExtended[i];    break;
        case POLYMER : return _solDispPolymerPeriodicallyExtended[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    const VectorType& getSolutionDisplacementAndMultiplier(int i) const { 
    switch( MaterialType ){
        case BONE :    return _solDispBoneAndMultiplier[i];    break;
        case POLYMER : return _solDispPolymerAndMultiplier[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    
    template<MaterialTypeBonePolymer MaterialType>
    VectorType& SolutionDisplacementAndMultiplier(int i) const{ 
    switch( MaterialType ){
        case BONE :    return _solDispBoneAndMultiplier[i];    break;
        case POLYMER : return _solDispPolymerAndMultiplier[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    const VectorType& getAffineDisplacement(int i) const { 
    switch( MaterialType ){
        case BONE :    return _affineDispBone[i];    break;
        case POLYMER : return _affineDispPolymer[i]; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
    const int getNumDisplacementDofs( ) const {return _dimDomain * _conf.getNumGlobalDofs(); }
    const int getNumPeriodicDisplacementDofs( ) const {return _dimDomain * _conf.getNumGlobalDofs(); }
    
    const int getNumLoads() const {return _numLoads;}

};



template <typename MatOptConfigurator >
class OptimalDeformSolverMultipleLoad : public OptimalDeformSolverInterfaceMultipleLoad<MatOptConfigurator> {

  typedef typename MatOptConfigurator::ConfiguratorType        ConfiguratorType;
  typedef typename ConfiguratorType::DTContainer               DataTypeContainer;
  typedef typename ConfiguratorType::RealType                  RealType;
  typedef typename ConfiguratorType::MaskType                  MaskType;
  typedef typename ConfiguratorType::InitType                  MeshType;
  typedef typename ConfiguratorType::VectorType                VectorType;
  typedef typename ConfiguratorType::SparseMatrixType          SparseMatrixType;
  typedef typename DataTypeContainer::ParameterParserType      ParameterParserType;
  typedef typename DataTypeContainer::TripletType              TripletType;
 
protected :
  mutable SparseMatrixType _HessianLinElastBone, _HessianLinElastPolymer;
  mutable SparseMatrixType _HessianLinElastMixedBone, _HessianLinElastMixedPolymer;
  mutable SparseMatrixType _HessianLinElastAffineBone, _HessianLinElastAffinePolymer;
  mutable std::vector<VectorType> _rhsBone, _rhsPolymer;
  
  public:
      OptimalDeformSolverMultipleLoad( const ParameterParserType &Parser,  
                                       const MatOptConfigurator & matOpConf,
                                       const VectorType & Phasefield,
                                       const QuocHandler<ConfiguratorType> & quocHandler,
                                       const std::vector<VectorType> &affineDispBone, 
                                       const std::vector<VectorType> &affineDispPolymer
                                     ) :
         OptimalDeformSolverInterfaceMultipleLoad<MatOptConfigurator> ( Parser, matOpConf, Phasefield, quocHandler, affineDispBone, affineDispPolymer ),
         _HessianLinElastBone( this->_dimDomain * ( matOpConf._conf.getNumGlobalDofs() + 1) , this->_dimDomain * (matOpConf._conf.getNumGlobalDofs() + 1)  ),
         _HessianLinElastPolymer( this->_dimDomain * ( matOpConf._conf.getNumGlobalDofs() + 1) , this->_dimDomain * (matOpConf._conf.getNumGlobalDofs() + 1)  ),
         _HessianLinElastMixedBone( this->_dimDomain * matOpConf._conf.getNumGlobalDofs(), this->_numAffineSymGradDofs ),
         _HessianLinElastMixedPolymer( this->_dimDomain * matOpConf._conf.getNumGlobalDofs(), this->_numAffineSymGradDofs ),
         _HessianLinElastAffineBone( this->_numAffineSymGradDofs, this->_numAffineSymGradDofs ),
         _HessianLinElastAffinePolymer( this->_numAffineSymGradDofs, this->_numAffineSymGradDofs )
          
        {
            for( int i=0; i<this->_numLoads; ++i ){
                _rhsBone.push_back ( VectorType( this->_dimDomain * ( matOpConf._conf.getNumGlobalDofs() + 1 ) ) );
                _rhsPolymer.push_back ( VectorType ( this->_dimDomain * ( matOpConf._conf.getNumGlobalDofs() + 1 ) ) );
            }
            this->assembleLinElastHessian();
            
            if( this->getParser().template get<bool>("ConstraintProblem.solveWithDirectSolver") ){  
                this->template getDirectLinearSolver<BONE>().analyzePattern( this->getHessianLinElast<BONE> () );
                this->template getDirectLinearSolver<POLYMER>().analyzePattern( this->getHessianLinElast<POLYMER> () );
            }
            
            solve();
        }

      
    void assembleLinElastHessian( ) const {
        #ifdef USE_OPENMP
        omp_set_nested(1);    
        #pragma omp parallel sections
        #endif
        {
            #ifdef USE_OPENMP 
            #pragma omp section
            #endif
            { 
                LinElastHessian<MatOptConfigurator,BONE> ( this->_matOpConf, this->_pfExtended ).assemblePeriodic( _HessianLinElastBone, this->_mask, this->_periodicIndices );
                _HessianLinElastBone.makeCompressed();
            }
            #ifdef USE_OPENMP
            #pragma omp section
            #endif
            { 
                LinElastHessian<MatOptConfigurator,POLYMER> ( this->_matOpConf, this->_pfExtended ).assemblePeriodic( _HessianLinElastPolymer, this->_mask, this->_periodicIndices );
                _HessianLinElastPolymer.makeCompressed();
            }
            #ifdef USE_OPENMP
            #pragma omp section
            #endif
            { 
                LinElastHessian<MatOptConfigurator,BONE> ( this->_matOpConf, this->_pfExtended ).assembleMixedPeriodic( _HessianLinElastMixedBone,this->_mask, this->_periodicIndices );
                _HessianLinElastMixedBone.makeCompressed();
            }
            #ifdef USE_OPENMP
            #pragma omp section
            #endif
            { 
                 LinElastHessian<MatOptConfigurator,POLYMER> ( this->_matOpConf, this->_pfExtended ).assembleMixedPeriodic( _HessianLinElastMixedPolymer, this->_mask, this->_periodicIndices );
                 _HessianLinElastMixedPolymer.makeCompressed();
            }
            #ifdef USE_OPENMP
            #pragma omp section
            #endif
            { 
                LinElastHessian<MatOptConfigurator,BONE> ( this->_matOpConf, this->_pfExtended ).assembleAffine( _HessianLinElastAffineBone );
                _HessianLinElastAffineBone.makeCompressed();
            }
            #ifdef USE_OPENMP
            #pragma omp section
            #endif
            { 
                LinElastHessian<MatOptConfigurator,POLYMER> ( this->_matOpConf, this->_pfExtended ).assembleAffine( _HessianLinElastAffinePolymer );
                _HessianLinElastAffinePolymer.makeCompressed();
            }
        }
       
        
        for( int i=0; i<this->_numLoads; ++i ){
            VectorType RHSBONE = _HessianLinElastMixedBone * this->template getAffineDisplacement<BONE> ( i );
            VectorType RHSPOLYMER = _HessianLinElastMixedPolymer * this->template getAffineDisplacement<POLYMER> ( i );
        
            for( int j=0; j < RHSBONE.size(); ++j ) (_rhsBone[i])[j] = - RHSBONE[j];
            for( int j=0; j<this->_dimDomain; ++j ) (_rhsBone[i])[j + this->getNumDisplacementDofs()] = 0.0;
            for( int j=0; j < RHSPOLYMER.size(); ++j ) (_rhsPolymer[i])[j] = - RHSPOLYMER[j];
            for( int j=0; j<this->_dimDomain; ++j ) (_rhsPolymer[i])[j + this->getNumDisplacementDofs()] = 0.0;

        }
      }
        
        
    template< MaterialTypeBonePolymer MaterialType>
    const SparseMatrixType& getHessianLinElast() const { 
      switch( MaterialType ){
        case BONE :
            return _HessianLinElastBone;
            break;
        case POLYMER :
            return _HessianLinElastPolymer;
            break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      } 
    }
    
    
    template< MaterialTypeBonePolymer MaterialType>
    const SparseMatrixType& getHessianLinElastMixedPart() const { 
      switch( MaterialType ){
        case BONE :
            return _HessianLinElastMixedBone;
            break;
        case POLYMER :
            return _HessianLinElastMixedPolymer;
            break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      } 
    }
    
    template< MaterialTypeBonePolymer MaterialType>
    const SparseMatrixType& getHessianLinElastAffinePart() const { 
      switch( MaterialType ){
        case BONE :
            return _HessianLinElastAffineBone;
            break;
        case POLYMER :
            return _HessianLinElastAffinePolymer;
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
    const SparseMatrixType& getSystemMatForAdjointProblem () const { return getHessianLinElast<MaterialType> (); }
    
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
        this->assembleLinElastHessian();
        solve();
      }
      return NewPhasefieldDiffersFromOld;
    }

protected:
    
  template<MaterialTypeBonePolymer MaterialType>
  void solve( ) const {
      
      if( this->getParser().template get<bool>("ConstraintProblem.solveWithDirectSolver") ){  
        this->template getDirectLinearSolver<MaterialType>().prepareSolver( this->getHessianLinElast<MaterialType> () );
      }
      
      for( int i=0; i<this->_numLoads; ++i ){
          
        if( this->getParser().template get<bool>("ConstraintProblem.solveWithDirectSolver") ){  
            this->template getDirectLinearSolver<MaterialType>().solve( this->template SolutionDisplacementAndMultiplier<MaterialType>(i), this->getRHS<MaterialType>(i) ); 
        }else{
            this->template getIterativeLinearSolver<MaterialType>().solve( this->getHessianLinElast<MaterialType> (), this->template SolutionDisplacementAndMultiplier<MaterialType>(i), this->getRHS<MaterialType>(i), this->getParser().template get<RealType> ("ConstraintProblem.toleranceLinearSystem"), this->getParser().template get<RealType>("ConstraintProblem.maxItersFacLinearSystem") ); 
        }
        
        VectorType solDispPeriodic ( this->template getSolutionDisplacementAndMultiplier<MaterialType>(i).segment( 0, this->_dimDomain * this->_conf.getNumGlobalDofs() ) );
        this->template setSolutionDisplacementPeriodic<MaterialType> ( i, solDispPeriodic );
        
        this->_quocHandler.extendMultiVectorPeriodically( this->template getSolutionDisplacementPeriodic<MaterialType>(i), this->template SolutionDisplacementPeriodicallyExtended<MaterialType>(i) );
        
      }
  }
  
  
  void solve() const {
        this->solve<BONE>();
        this->solve<POLYMER>();
  }
    
  
};

}//end namespace

#endif
