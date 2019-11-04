#ifndef __PERIODICHOMOGENIZATIONBONESMULTIPLELOADMATERIALOPTIMIZATIONBONESSOLVER_H
#define __PERIODICHOMOGENIZATIONBONESMULTIPLELOADMATERIALOPTIMIZATIONBONESSOLVER_H

#include <derivativeTest.h>
#include <ipoptNonlinearConstraintSolver.h>
#include <LineSearchMethods.h>
#include <quocHandler.h>

#include "BonesEnergiesWithMaterial.h"
#include "BonesOptDeformSolver.h"
#include "BonesMatOptCompliance.h"
#include "BonesMatOptEnergies.h"
#include "BonesPlotter.h"
#include "BonesResultsPlotter.h"

#include "BonesLinDiffusionConstraints.h"

using namespace quocFE;

namespace shapeOptBonePolymerPeriodicHomogenization{

template<typename MatOptConfigurator>
class MaterialOptimizationMultipleLoadSolver {
  
  typedef typename MatOptConfigurator::ConfiguratorType         ConfiguratorType;
  typedef typename ConfiguratorType::RealType                   RealType;
  typedef typename ConfiguratorType::MaskType                   MaskType;
  typedef typename ConfiguratorType::InitType                   MeshType;
  typedef typename ConfiguratorType::DTContainer                DataTypeContainer;
  typedef typename DataTypeContainer::PointType                 PointType;
  typedef typename ConfiguratorType::DomVecType                 DomVecType;
  typedef typename ConfiguratorType::VectorType                 VectorType;
  typedef typename ConfiguratorType::SparseMatrixType           SparseMatrixType;
  typedef typename DataTypeContainer::ParameterParserType       ParameterParserType;

  mutable ParameterParserType _parser;
  const int _maxIterations;
  const RealType _tolerance;

public:
  MaterialOptimizationMultipleLoadSolver ( const ParameterParserType &Parser ) :
  _parser( Parser ),
  _maxIterations ( Parser.template get<int> ("MaterialOptimization.maxIterations") ),
  _tolerance( Parser.template get<RealType> ("MaterialOptimization.tolerance") ) { } 
  
  void computeOptimalDesign( const MeshType & mesh, VectorType & initMaterial, VectorType &solutionMaterial, 
                             const int refinementStep = 0, 
                             const bool onlyComputeInitialEnergy = false,
                             const bool solveWithVolumeConstraint = false,
                             const bool solveWithDiffusionContraint = false
                           ) const {
 
    ParameterParserType parser; parser = _parser;
    const string saveDirectoryRefinementStep = parser.createSubDirectory( "refinementStep" + std::to_string(refinementStep) );
    parser.set ( "saving.saveDirectory", saveDirectoryRefinementStep );
    const string saveDirectoryDeformation = parser.createSubDirectory( "Deformation" );
    
    auto startTime = std::chrono::high_resolution_clock::now(); 
    parser.dumpToFile( aol::strprintf ( "ParameterParser.ini" ) );
    if( onlyComputeInitialEnergy ) aol::consoleOutput( aol::strprintf( "compute initial energy for refinementStep %d", refinementStep ).c_str() );
    else                           aol::consoleOutput( aol::strprintf( "compute optimal design for refinementStep %d", refinementStep ).c_str() );
    
    ConfiguratorType conf ( mesh );
    const int dimDomain = conf.dimDomain;
    const int numAffineSymGradDofs = conf.numAffineSymGradDofs;
    QuocHandler<ConfiguratorType> quocHandler ( parser, conf );
    
    //! make sure that initial material is periodic
    quocHandler.collabseVectorPeriodically( initMaterial );
    
    //!initialize eps_area
    const RealType eps_factor = parser.template get<RealType>( "MaterialOptimization.eps_factor" );
    const RealType eps_area = mesh.getInterfaceWith() * eps_factor;
    
    //read loads
    const int numLoads =parser.template get<int> ( "AffineDisp.numLoads" );
    std::vector<VectorType> affineDispBone, affineDispPolymer;
    for( int i=1; i <= numLoads; ++i ){
        VectorType affineDisp ( numAffineSymGradDofs );
        parser.template getFixSizeVector<RealType,VectorType> ( aol::strprintf("AffineDisp.Load%d", i ).c_str(), affineDisp );
        affineDispBone.push_back ( affineDisp ); affineDispPolymer.push_back ( affineDisp );
    }
    
    //!determine weightFunction for Multiple Load 
    std::vector<RealType> weightFct_Bone_weightVec ( numLoads ), weightFct_Polymer_weightVec ( numLoads );
    parser.template getFixSizeVector<RealType,std::vector<RealType>> ( "MaterialOptimization.weightFct_Bone_weightVec", weightFct_Bone_weightVec );
    parser.template getFixSizeVector<RealType,std::vector<RealType>> ( "MaterialOptimization.weightFct_Polymer_weightVec", weightFct_Polymer_weightVec );

    //initialize MaterialOptimizationConfigurator
    MatOptConfigurator matOptConf ( parser, conf, eps_area, eps_factor, weightFct_Bone_weightVec, weightFct_Polymer_weightVec );
    
    //! initialize state equation
    OptimalDeformSolverMultipleLoad<MatOptConfigurator> OptDeformFinder ( parser, matOptConf, initMaterial, quocHandler, affineDispBone, affineDispPolymer );

    //! initialize volume operators
    VolumeConstraintPeriodicBC<MatOptConfigurator,BONE> volumeOpBone ( matOptConf, quocHandler );
    VolumeConstraintPeriodicBC<MatOptConfigurator,POLYMER> volumeOpPolymer ( matOptConf, quocHandler );
    PointType barycenterPoint; for(int i=0; i<barycenterPoint.size(); ++i ) barycenterPoint[i] = 0.5;
    BaryCenterConstraintBonesPeriodicBC<MatOptConfigurator,BONE> barycenterConstraintOp ( OptDeformFinder.getMatOptConfigurator(), quocHandler, barycenterPoint );     
    
    //! construct energy operators
    MaterialOptimizationMultipleLoadEnergyOp<MatOptConfigurator> totalEnergyOp ( matOptConf, OptDeformFinder );
    MaterialOptimizationMultipleLoadEnergyInfo<MatOptConfigurator> energyInfoInit;
    this->saveResults( "InitMaterial", OptDeformFinder, totalEnergyOp, volumeOpBone, volumeOpPolymer, barycenterConstraintOp, energyInfoInit, initMaterial, saveDirectoryRefinementStep, saveDirectoryDeformation  );
    
    //=========================================================================================================== 
    if( onlyComputeInitialEnergy ){
        TikzPlotterMaterialOptimizationBonesMultipleLoad<DataTypeContainer> tikzPlotterInit ( parser );
        tikzPlotterInit.template plotInit<MatOptConfigurator>( matOptConf, energyInfoInit );
        std::ofstream outParameter ( aol::strprintf ( "%s/Parameter.tex", saveDirectoryRefinementStep.c_str() ).c_str ()  );
        tikzPlotterInit.template plotParameters<MatOptConfigurator>( outParameter, matOptConf );
        return;
    }
  
    
    //=========================================================================================================== 
    SolverInfo<DataTypeContainer> solverInfo;
   
    // Solver - IPOPT with barycenter constraint
    if( solveWithVolumeConstraint ){
        
        aol::consoleOutput( "solve with volume constraint" );
        
        VolumeAndBarycenterConstraintPeriodicBC<MatOptConfigurator,BONE> volumeAndBarycenterConstraintOp ( OptDeformFinder.getMatOptConfigurator(), quocHandler, barycenterPoint ); 
        std::vector<RealType> constraint_l ( 1 + barycenterPoint.size() ), constraint_u ( 1 + barycenterPoint.size() );
        constraint_l[0] = parser.template get<RealType>("MaterialOptimization.VolumeConstraint_l");
        constraint_u[0] = parser.template get<RealType>("MaterialOptimization.VolumeConstraint_u");
        for( int i=0; i < barycenterPoint.size(); ++i ){ constraint_l[i+1] = 0.0; constraint_u[i+1] = 0.0;}

        IpoptInterfaceMultipleNLC<DataTypeContainer> ipoptSolver ( totalEnergyOp, volumeAndBarycenterConstraintOp, _maxIterations, _tolerance, parser.template get<RealType>("MaterialOptimization.BoxConstraint_l"), parser.template get<RealType>("MaterialOptimization.BoxConstraint_u" ),
        constraint_l, constraint_u, parser.template get<int> ( "MaterialOptimization.linearSolverTypeIpopt") );
        ipoptSolver.solve( initMaterial, solutionMaterial, solverInfo );
    }else{
        
        if( solveWithDiffusionContraint ){
            
            aol::consoleOutput( "solve with diffusion constraint" );
            
            //read loads
            std::vector<DomVecType> affineDiffusionBone, affineDiffusionPolymer;
            const int numDiffusionLoadsBone = parser.template get<int> ( "AffineDiffusion.numLoadsBone" );
            const int numDiffusionLoadsPolymer = parser.template get<int> ( "AffineDiffusion.numLoadsPolymer" );
            for( int i=1; i <= numDiffusionLoadsBone; ++i ){
                DomVecType affineDiffusion;
                parser.template getFixSizeVector<RealType,DomVecType> ( aol::strprintf("AffineDiffusion.LoadBone%d", i ).c_str(), affineDiffusion );
                affineDiffusionBone.push_back ( affineDiffusion );
            }
            for( int i=1; i <= numDiffusionLoadsPolymer; ++i ){
                DomVecType affineDiffusion;
                parser.template getFixSizeVector<RealType,DomVecType> ( aol::strprintf("AffineDiffusion.LoadPolymer%d", i ).c_str(), affineDiffusion );
                affineDiffusionPolymer.push_back ( affineDiffusion );
            }
            
            OptimalDiffusionSolverMultipleLoad<MatOptConfigurator> OptDiffusionFinder ( parser, matOptConf, initMaterial, quocHandler, affineDiffusionBone, affineDiffusionPolymer );
            
            DiffusionAndBarycenterConstraintPeriodicBC<MatOptConfigurator> diffusionAndBarycenterConstraintOp ( OptDeformFinder.getMatOptConfigurator(), quocHandler, OptDiffusionFinder, barycenterPoint ); 
           
            std::vector<RealType> constraint_l ( affineDiffusionBone.size() + affineDiffusionPolymer.size() + barycenterPoint.size() ), 
                                  constraint_u ( affineDiffusionBone.size() + affineDiffusionPolymer.size() + barycenterPoint.size() );
            for( int i=0; i < affineDiffusionBone.size(); ++i ){ 
                constraint_l[i] = parser.template get<RealType> ( aol::strprintf("AffineDiffusion.lowerBoundBone%d", i+1 ).c_str() ); 
                constraint_u[i] = 2.e+19; 
            }
            for( int i=0; i < affineDiffusionPolymer.size(); ++i ){ 
                constraint_l[i+affineDiffusionBone.size()] = parser.template get<RealType> ( aol::strprintf("AffineDiffusion.lowerBoundPolymer%d", i+1 ).c_str() ); ; 
                constraint_u[i+affineDiffusionBone.size()] = 2.e+19;
            }
            for( int i=0; i < barycenterPoint.size(); ++i ){ 
                constraint_l[i+affineDiffusionBone.size() + affineDiffusionPolymer.size()] = 0.0; 
                constraint_u[i+affineDiffusionBone.size() + affineDiffusionPolymer.size()] = 0.0; 
            }
            
            IpoptInterfaceMultipleNLC<DataTypeContainer> ipoptSolver ( totalEnergyOp, diffusionAndBarycenterConstraintOp, _maxIterations, _tolerance, parser.template get<RealType>("MaterialOptimization.BoxConstraint_l"), parser.template get<RealType>("MaterialOptimization.BoxConstraint_u" ), constraint_l, constraint_u, parser.template get<int> ( "MaterialOptimization.linearSolverTypeIpopt") );
            
            //=========================================================================================================== 
            //First Derivative Test
            if( parser.template get<int>( "DerivativeTest.order" ) ==  1 ){
//                 ipoptSolver.testDerivative( initMaterial );
                aol::consoleOutput( "test first Derivative of diffusion constraint" );
                VectorType testMaterial = VectorType::Random( initMaterial.size() );
                quocHandler.collabseVectorPeriodically( testMaterial );
                OptDeformFinder.updatePhasefield( testMaterial );
                OptDiffusionFinder.updatePhasefield( testMaterial );
                DiffusionOpInMaterialPeriodicBC<MatOptConfigurator,BONE> diffusionOp ( OptDeformFinder.getMatOptConfigurator(), quocHandler, OptDiffusionFinder, 0 );
                DerivativeTester<DataTypeContainer> ( diffusionOp, parser.template get<RealType>( "DerivativeTest.stepSize" ), testMaterial, 1, parser.template get<RealType> ( "DerivativeTest.tolerance" ) ).testFirstDerivative_AllDirections( );
                return;
            }
            //=========================================================================================================== 
            
            
            ipoptSolver.solve( initMaterial, solutionMaterial, solverInfo );
            
            
            aol::consoleOutput( "diffusionAndBarycenterConstraintOp for solution" );
            for( int i=0; i < affineDiffusionBone.size(); ++i ){ 
                int numConstraint = i;
                RealType val;
                diffusionAndBarycenterConstraintOp.evaluateEnergy( numConstraint, solutionMaterial, val );
                cout << "bone diffusion for load " << i << " : " << val << endl;
            }
            for( int i=0; i < affineDiffusionPolymer.size(); ++i ){ 
                int numConstraint = i + affineDiffusionBone.size();
                RealType val;
                diffusionAndBarycenterConstraintOp.evaluateEnergy( numConstraint, solutionMaterial, val );
                cout << "polymer diffusion for load " << i << " : " << val << endl;
            }
            
            
        }else{
            std::vector<RealType> constraint_l ( barycenterPoint.size() ), constraint_u ( barycenterPoint.size() );
            for( int i=0; i < barycenterPoint.size(); ++i ){ constraint_l[i] = 0.0; constraint_u[i] = 0.0;}
            
            IpoptInterfaceMultipleNLC<DataTypeContainer> ipoptSolver ( totalEnergyOp, barycenterConstraintOp, _maxIterations, _tolerance, parser.template get<RealType>("MaterialOptimization.BoxConstraint_l"), parser.template get<RealType>("MaterialOptimization.BoxConstraint_u" ), constraint_l, constraint_u, parser.template get<int> ( "MaterialOptimization.linearSolverTypeIpopt") );
            ipoptSolver.solve( initMaterial, solutionMaterial, solverInfo );
        }
    }
    

    //===========================================================================================================    
    //! save solution results 
    MaterialOptimizationMultipleLoadEnergyInfo<MatOptConfigurator> energyInfoSolution;
    this->saveResults( "SolMaterial", OptDeformFinder, totalEnergyOp, volumeOpBone, volumeOpPolymer, barycenterConstraintOp, energyInfoSolution, solutionMaterial, saveDirectoryRefinementStep, saveDirectoryDeformation );
    
    TikzPlotterMaterialOptimizationBonesMultipleLoad<DataTypeContainer> tikzPlotter ( parser );
    tikzPlotter.template plotAll<MatOptConfigurator>( matOptConf, energyInfoInit, energyInfoSolution, solverInfo );
    
    //! print elapsed time
    std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << endl << "duration of refinementStep = " << diff.count() << " sec" << endl;
    aol::consoleOutput( aol::strprintf( "finished computation optimal design for refinementStep %d", refinementStep ).c_str() );
    
    //!plot results
    if( parser.template get<int> ( "saving.plotResults" ) == 1 ) this->plotResults( parser );
    
 }
  
 void computeOptimalDesign ( const MeshType &mesh ) const {
    aol::consoleOutput( "compute optimal design" );
    ConfiguratorType conf ( mesh );
    QuocHandler<ConfiguratorType> quocHandler ( _parser, conf );
    //! initialize material
    VectorType initMaterial ( conf.getNumGlobalDofs() );
    quocHandler.switchMaterialType ( initMaterial );
    VectorType solMaterial ( initMaterial.size() );
    computeOptimalDesign ( mesh, initMaterial, solMaterial ); 
  }
  
  void computeOptimalDesign_MultiLevel( const MeshType &meshStartLevel, const int numRefinementSteps, const int StartLevel = 0 ) const {
    aol::consoleOutput( "compute optimal design with multilevel scheme" );
    cout << "we use " << numRefinementSteps << " adaptive refinement steps" << endl;
    ConfiguratorType confStartLevel ( meshStartLevel );
    QuocHandler<ConfiguratorType> quocHandlerStartLevel ( _parser, confStartLevel );
    //! initialize material
    VectorType initMaterialStartLevel ( confStartLevel.getNumGlobalDofs() );
    quocHandlerStartLevel.switchMaterialType ( initMaterialStartLevel );
    VectorType solMaterialStartLevel ( initMaterialStartLevel.size() );
    VectorType initMaterialCurrentLevel, solMaterialCurrentLevel;
    typename MeshType::IntVecType numDofVec; meshStartLevel.getNumDofVec( numDofVec );
    typename MeshType::IntVecType numDofVecStartLevel; meshStartLevel.getNumDofVec( numDofVecStartLevel );
    PointType lenghtVec; meshStartLevel.getLenghtVec( lenghtVec );
    
    const RealType factorMultilevel = _parser.template get<RealType> ( "MaterialOptimization.factorMultiLevel" );
    
    bool useVolumeConstraint = false; if( _parser.template get<int> ( "MaterialOptimization.useVolumeConstraint" )  ) useVolumeConstraint = true;
    bool useDiffusionConstraint = false; if( _parser.template get<int> ( "MaterialOptimization.useDiffusionConstraint" )  ) useDiffusionConstraint = true;
    
    for( int level = StartLevel; level <= numRefinementSteps; ++level ){
        if( level > StartLevel ){   
            MeshType meshOldLevel( numDofVec, lenghtVec );
            ConfiguratorType confOldLevel ( meshOldLevel );
            QuocHandler<ConfiguratorType> quocHandlerOldLevel( _parser, confOldLevel );
            quocHandlerOldLevel.extendVectorPeriodically( solMaterialCurrentLevel );
            QuocDiscreteFunctionDefault<ConfiguratorType> discreteFctOld ( confOldLevel, solMaterialCurrentLevel ); //note: here solMaterialCurrentLevel is actually OldLevel
            // size_k = round( fac^k (size_{k-1} - 1 ) + 1 )
            for( int i=0; i<numDofVec.size(); ++i ) numDofVec[i] =  std::round( std::pow(factorMultilevel, level - StartLevel ) * ( numDofVecStartLevel[i] - 1.) + 1. );
            _parser.template setFixSizeVector<typename MeshType::IntVecType>("InputMesh.NumDofVec", numDofVec, true );
            _parser.template setFixSizeVector<PointType>("InputMesh.LengthVec", lenghtVec, false );
            MeshType meshCurrentLevel ( numDofVec, lenghtVec );
            initMaterialCurrentLevel.resize ( meshCurrentLevel.getNumVertices() );
            for( int nodeIdxCurrentLevel=0; nodeIdxCurrentLevel < meshCurrentLevel.getNumVertices(); nodeIdxCurrentLevel++ ){
                const PointType& GlobalCoords = meshCurrentLevel.getVertex ( nodeIdxCurrentLevel );
                int elementNumberOld; PointType LocalCoordOld;
                confOldLevel.getLocalCoords ( GlobalCoords, elementNumberOld, LocalCoordOld );
                initMaterialCurrentLevel[nodeIdxCurrentLevel] = discreteFctOld.evaluate( meshOldLevel.getElement(elementNumberOld), LocalCoordOld );
            }
            solMaterialCurrentLevel.resize( initMaterialCurrentLevel.size() );
            computeOptimalDesign ( meshCurrentLevel, initMaterialCurrentLevel, solMaterialCurrentLevel, level, false, useVolumeConstraint, useDiffusionConstraint ); 
        }else{
            MeshType meshCurrentLevel ( meshStartLevel );  
            initMaterialCurrentLevel.resize( initMaterialStartLevel.size() ); initMaterialCurrentLevel = initMaterialStartLevel;
            solMaterialCurrentLevel.resize( initMaterialStartLevel.size() );
            computeOptimalDesign ( meshCurrentLevel, initMaterialCurrentLevel, solMaterialCurrentLevel, level, false, useVolumeConstraint, useDiffusionConstraint ); 
        }
    }
    
    
    if( _parser.template get<int> ( "MaterialOptimization.comparePerimeterMinimizer" ) == 1 ){
        _parser.set( "MaterialOptimization.factorComplianceCost", 0.0 );
        _parser.set( "MaterialOptimization.factorInterfaceCost", 1.e-0 );
        MeshType meshCompare ( numDofVec, lenghtVec );
        VectorType initMaterialCompare ( solMaterialCurrentLevel.size() ), solutionMaterialCompare ( solMaterialCurrentLevel.size() );
        initMaterialCompare = solMaterialCurrentLevel;
        computeOptimalDesign( meshCompare, initMaterialCompare, solutionMaterialCompare, numRefinementSteps + 1, false, true, false );
    }
  }
  
private :
    
  void saveResults( const string name, const OptimalDeformSolverMultipleLoad<MatOptConfigurator> & OptDeformFinder, 
                    const MaterialOptimizationMultipleLoadEnergyOp<MatOptConfigurator> &totalEnergyOp,
                    const VolumeConstraintPeriodicBC<MatOptConfigurator,BONE> &volumeOpBone, const VolumeConstraintPeriodicBC<MatOptConfigurator,POLYMER> &volumeOpPolymer,
                    const BaryCenterConstraintBonesPeriodicBC <MatOptConfigurator,BONE> &barycenterConstraintOp,
                    MaterialOptimizationMultipleLoadEnergyInfo<MatOptConfigurator> &energyInfo,
                    const VectorType &material, 
                    const string saveDirectoryRefinementStep,
                    const string saveDirectoryDeformation
//                     , const bool saveDeformations
                  ) const{
      
    aol::consoleOutput( "Energy of " + name + " : " );
    RealType energy; VectorType residual ( material.size() );
    totalEnergyOp.evaluateEnergy( material, energy ); totalEnergyOp.evaluateJacobian( material, residual );
    totalEnergyOp.printLastEnergies(); totalEnergyOp.printLastResiduals();
    totalEnergyOp.getEnergyInfo( energyInfo );
    RealType volumeBone; volumeOpBone.evaluateEnergy( material, volumeBone );
    RealType volumePolymer; volumeOpPolymer.evaluateEnergy( material, volumePolymer );
    energyInfo.setVolumeBone( volumeBone ); energyInfo.setVolumePolymer( volumePolymer );
    PointType barycenterDiff;
    for( int i=0; i<barycenterConstraintOp.getNumConstraints(); ++i ) barycenterConstraintOp.evaluateEnergy( i, material, barycenterDiff[i] );
    PointType barycenter ( barycenterConstraintOp.getBarycenterPoint() + barycenterDiff );
    energyInfo.setBarycenterBone( barycenter ); 
     
    aol::printVectorToFile<VectorType> ( material, aol::strprintf( "%s/%s.txt", saveDirectoryRefinementStep.c_str(), name.c_str() ).c_str(), 10 ); 
    
    VectorType materialScaledTo01( material.size() );
    for( int i=0; i<material.size(); i++ ) materialScaledTo01[i] = 0.5 * ( material[i] + 1.0 );
    OptDeformFinder.getQuocHandler().extendVectorPeriodically( materialScaledTo01 );
    
        for( int loadIdx=0; loadIdx< OptDeformFinder.getNumLoads(); ++loadIdx ){
            aol::printVectorToFile<VectorType> ( OptDeformFinder.template getSolutionDisplacementPeriodic<BONE>(loadIdx), 
                                                aol::strprintf( "%s/%s_DisplacementBonePeriodic_Dir%d.txt", saveDirectoryDeformation.c_str(), name.c_str(), loadIdx ).c_str(), 10 ); 
            aol::printVectorToFile<VectorType> ( OptDeformFinder.template getAffineDisplacement<BONE>(loadIdx), 
                                                aol::strprintf( "%s/%s_DisplacementBoneAffine_Dir%d.txt", saveDirectoryDeformation.c_str(), name.c_str(), loadIdx ).c_str(), 10 ); 
            aol::printVectorToFile<VectorType> ( OptDeformFinder.template getSolutionDisplacementPeriodic<POLYMER>(loadIdx), 
                                                aol::strprintf( "%s/%s_DisplacementPolymerPeriodic_Dir%d.txt", saveDirectoryDeformation.c_str(), name.c_str(), loadIdx ).c_str(), 10 ); 
            aol::printVectorToFile<VectorType> ( OptDeformFinder.template getAffineDisplacement<POLYMER>(loadIdx), 
                                                aol::strprintf( "%s/%s_DisplacementPolymerAffine_Dir%d.txt", saveDirectoryDeformation.c_str(), name.c_str(), loadIdx ).c_str(), 10 ); 

            OptDeformFinder.getQuocHandler().plotPeriodicPlusAffineSymGradDispWithMaterial ( OptDeformFinder.template getSolutionDisplacementPeriodic<BONE>(loadIdx), OptDeformFinder.template getAffineDisplacement<BONE>(loadIdx),  materialScaledTo01, aol::strprintf( "%s/%s_DeformationBone_Dir%d", saveDirectoryDeformation.c_str(), name.c_str(), loadIdx ).c_str() );
            OptDeformFinder.getQuocHandler().plotPeriodicPlusAffineSymGradDispWithMaterial ( OptDeformFinder.template getSolutionDisplacementPeriodic<POLYMER>(loadIdx), OptDeformFinder.template getAffineDisplacement<POLYMER>(loadIdx), materialScaledTo01, aol::strprintf( "%s/%s_DeformationPolymer_Dir%d", saveDirectoryDeformation.c_str(), name.c_str(), loadIdx ).c_str() );
        }
    
    OptDeformFinder.getQuocHandler().plotUndeformedWithMaterial ( materialScaledTo01, aol::strprintf( "%s_Undeformed", name.c_str() ).c_str() );
  }
  
  void plotResults( ParameterParserType &parserRefinementStep ) const{
    
    BonesResultsPlotter<MatOptConfigurator> plotter( parserRefinementStep );
 
    plotter.plotSingleCell( );
  
    const int numBlocksToPlot = _parser.template get<int> ("saving.numBlocks");
    plotter.plotBlockOfCells( numBlocksToPlot );
    
    const int plotInterface = _parser.template get<int> ("saving.plotInterface");
    if( plotInterface == 1 ) plotter.plotInterface();
    
    const int plotVonMises = _parser.template get<int>("saving.plotVonMises");
    if( plotVonMises ) plotter.plotVonMisesStresses();
    
    plotter.computeHomogenizedTensor();
    
    const int InitOrSolution = 3; //init 1, sol 2, both 3
    plotter.plotVTKToPNG( InitOrSolution, plotVonMises );
    
    plotter.summarizeAllResults( );
  }
  
      
};

}//end namespace

#endif //__QUOCMATERIALOPTIMIZATIONSOLVER_H
