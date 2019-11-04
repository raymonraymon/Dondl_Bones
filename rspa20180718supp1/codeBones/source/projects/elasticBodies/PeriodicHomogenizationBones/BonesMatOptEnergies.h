#ifndef __PERIODICHOMOGENIZATIONBONESMULTIPLELOADMATERIALOPTIMIZATIONBONESENERGIES_H
#define __PERIODICHOMOGENIZATIONBONESMULTIPLELOADMATERIALOPTIMIZATIONBONESENERGIES_H

#include "../materialOptLinElast/quocEnergiesWithMaterial.h"
#include "BonesMatOptCompliance.h"


#define _REGULARIZEDMAXFUNCTION

using namespace quocFE;

namespace shapeOptBonePolymerPeriodicHomogenization{

    
template <typename MatOptConfigurator>
class MaterialOptimizationMultipleLoadEnergyInfo {
 
    typedef typename MatOptConfigurator::RealType RealType;
    typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
    typedef typename ConfiguratorType::PointType PointType;
    
public :
    
    int _numLoads;
    RealType _complianceWeightBone, _complianceWeightPolymer, _complianceEnergy, _interfaceEnergy;
    RealType _volumeBone, _volumePolymer;
    PointType _barycenterBone;
    std::vector<RealType> _complianceBone, _compliancePolymer;
    MaterialOptimizationMultipleLoadEnergyInfo ( ) {}
    
    void setNumLoads( const int numLoads ) { _numLoads = numLoads; }
    void setComplianceEnergy( const RealType complianceEnergy ) { _complianceEnergy = complianceEnergy; };
    void setComplianceWeightBone( const RealType compliance ) { _complianceWeightBone = compliance; };
    void setComplianceWeightPolymer( const RealType compliance ) { _complianceWeightPolymer = compliance; };
    void setInterfaceEnergy( const RealType interfaceEnergy ) { _interfaceEnergy = interfaceEnergy; };
    void setComplianceBone( const std::vector<RealType> & compliance ) { _complianceBone = compliance; };
    void setCompliancePolymer( const std::vector<RealType> & compliance ) { _compliancePolymer = compliance; };
    void setVolumeBone( const RealType vol ) { _volumeBone = vol; };
    void setVolumePolymer( const RealType vol ) { _volumePolymer = vol; };
    void setBarycenterBone( const PointType bar ) { _barycenterBone = bar; }
};
    
    
    
template <typename MatOptConfigurator>
class MaterialOptimizationMultipleLoadEnergyOp : public aol::NonlinearEnergyOp< typename MatOptConfigurator::ConfiguratorType::DTContainer > 
{
  protected :
    typedef typename MatOptConfigurator::ConfiguratorType::DTContainer DataTypeContainer;
    typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
    typedef typename MatOptConfigurator::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;

    const OptimalDeformSolverMultipleLoad<MatOptConfigurator> & _OptDeformFinder;
    const RealType _c_compl, _c_interface, _epsInterfaceLength;
    
    mutable RealType _lastAreaEnergy, _lastJphysEnergy, _lastTotalEnergy;
    mutable RealType _complianceWeightBone, _complianceWeightPolymer;
    mutable std::vector<RealType> _complianceBONE, _compliancePOLYMER;
    mutable RealType _lastAreaResidual, _lastJphysResidual, _lastTotalResidual;
    
  public:
    MaterialOptimizationMultipleLoadEnergyOp( const MatOptConfigurator &matOptConf,
                                       const OptimalDeformSolverMultipleLoad<MatOptConfigurator> & OptDeformFinder
                                     ) : 
         _OptDeformFinder( OptDeformFinder ), 
         _c_compl ( matOptConf._factorComplianceCost ),
         _c_interface( matOptConf._factorInterfaceCost ),
         _epsInterfaceLength( matOptConf._epsInterfaceLength ),
         _complianceBONE ( OptDeformFinder.getNumLoads() ),
         _compliancePOLYMER ( OptDeformFinder.getNumLoads() )
         { } 
   
   void evaluateEnergy( const VectorType & v, RealType& energy ) const {
         energy = 0;
         _OptDeformFinder.updatePhasefield(v);
         
         // compute area term
         PfOp_ModicaMortola<MatOptConfigurator>( _OptDeformFinder.getMatOptConfigurator(), _epsInterfaceLength ).apply( _OptDeformFinder.getPhaseFieldPeriodicallyExtended(), _lastAreaEnergy);
         energy += _lastAreaEnergy * _c_interface;
         
         // compute compliance
#ifdef _REGULARIZEDMAXFUNCTION
        ComplianceEnergyMultipleLoad_RegularizedMaxFunction<MatOptConfigurator> complianceOp ( _OptDeformFinder );
#endif
            
        complianceOp.apply( v, _lastJphysEnergy );
        energy += _lastJphysEnergy * _c_compl;
            
        _complianceWeightBone= complianceOp.getLastComplianceWeightBone();
        _complianceWeightPolymer = complianceOp.getLastComplianceWeightPolymer();
        for( int i=0; i<_OptDeformFinder.getNumLoads(); ++i ){
            _complianceBONE[i] = complianceOp.getLastComplianceBone(i);
            _compliancePOLYMER[i] = complianceOp.getLastCompliancePolymer(i);
        }
         
         _lastTotalEnergy = energy;
    }
    
    RealType getLastAreaEnergy()   const { return _lastAreaEnergy;}
    RealType getLastJphysEnergy()  const { return _lastJphysEnergy;}
    RealType getLastComplianceWeightBone()  const { return _complianceWeightBone;}
    RealType getLastComplianceWeightPolymer()  const { return _complianceWeightPolymer;}
    RealType getLastComplianceBone(const int i)  const { return _complianceBONE[i];}
    RealType getLastCompliancePolymer(const int i) const { return _compliancePOLYMER[i];}
    
    void getLastEnergies( RealType &Compliance, RealType &Interface ) const{ Compliance = _lastJphysEnergy; Interface = _lastAreaEnergy;  }
    void printLastEnergies() const {
     cout << "Energy:" << endl
          << std::fixed << std::setprecision( 12 )
          << "Interface  = " << getLastAreaEnergy()   << " , \t c_i Interface  = " << _c_interface * getLastAreaEnergy() << endl
          << "Compliance = " << getLastJphysEnergy()  << " , \t c_c Compliance = " << _c_compl * getLastJphysEnergy() << endl
          << "Total (weighted ) = " << _lastTotalEnergy << endl;
    }
   
   void evaluateJacobian( const VectorType &v, VectorType& Deriv ) const {  
         Deriv.setZero(); 
         _OptDeformFinder.updatePhasefield(v);
         VectorType grad ( v.size() );
         
         // compute derivative of area term:
         PfOp_ModicaMortolaDerivative<MatOptConfigurator>( _OptDeformFinder.getMatOptConfigurator(), _epsInterfaceLength ).apply( _OptDeformFinder.getPhaseFieldPeriodicallyExtended(), grad );
         _OptDeformFinder.getQuocHandler().collabseVectorPeriodicallyAdditive( grad );
         Deriv += _c_interface * grad;
         _lastAreaResidual = grad.norm();
         
         // compute gradient of Jphys w.r.t. phasefield Jphys,v 
#ifdef _REGULARIZEDMAXFUNCTION
//         ComplianceEnergyMultipleLoad_RegularizedMaxFunction_Derivative<MatOptConfigurator>( _OptDeformFinder, _adjointProblemBone, _adjointProblemPolymer ).apply( v, grad );
          ComplianceEnergyMultipleLoad_RegularizedMaxFunction_Derivative<MatOptConfigurator>( _OptDeformFinder ).apply( v, grad );
#endif
        Deriv += _c_compl * grad;
        _lastJphysResidual = grad.norm();
         
        _lastTotalResidual = Deriv.norm();
   }
   
    RealType getLastAreaResidual()   const { return _lastAreaResidual;}
    RealType getLastJphysResidual()  const { return _lastJphysResidual;}
    
    void printLastResiduals() const {
     cout << "Residuals :" << endl
          << "Interface  = " << getLastAreaResidual() << endl
          << "Compliance = " << getLastJphysResidual() << endl
          << "total      = " << _lastTotalResidual << endl; 
    }
   
   void evaluateHessian( const VectorType &/*v*/, SparseMatrixType& /*Hessian*/ ) const {  
       throw std::invalid_argument( aol::strprintf ( "Hessian is so far not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
   }
   
    void evaluateTripletListHessianSym( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const {
        throw std::invalid_argument( aol::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }

   void getEnergyInfo( MaterialOptimizationMultipleLoadEnergyInfo<MatOptConfigurator> & energyInfo ) const {
        energyInfo.setNumLoads( _complianceBONE.size() );
        energyInfo.setComplianceEnergy( _lastJphysEnergy  );
        energyInfo.setComplianceWeightBone( _complianceWeightBone );
        energyInfo.setComplianceWeightPolymer( _complianceWeightPolymer  );
        energyInfo.setInterfaceEnergy( _lastAreaEnergy ); 
        energyInfo.setComplianceBone( _complianceBONE );
        energyInfo.setCompliancePolymer( _compliancePOLYMER);
   }
   
};

}//end namespace

#endif //__PERIODICHOMOGENIZATIONBONESMULTIPLELOADMATERIALOPTIMIZATIONBONESENERGIES_H
