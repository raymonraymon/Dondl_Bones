#include <quocConfigurators.h>
#include <quocDefines.h>

#include "BonesOptDeformSolver.h"

using namespace quocFE;
using namespace shapeOptBonePolymerPeriodicHomogenization;

// #define _USEQUOC2D
#define _USEQUOC3D

#ifdef _USEQUOC2D
typedef Quoc2DDataTypeContainer                                                                                         DataTypeContainer;
typedef QuocConfigurator2D<DataTypeContainer>                                                                           ConfiguratorType;
#endif

#ifdef _USEQUOC3D
typedef Quoc3DDataTypeContainer                                                                                         DataTypeContainer;
typedef QuocConfigurator3D<DataTypeContainer>                                                                           ConfiguratorType;
#endif

typedef typename DataTypeContainer::RealType                                                                            RealType;
typedef typename DataTypeContainer::VectorType                                                                          VectorType;
typedef typename DataTypeContainer::IntVecType                                                                          IntVecType;
typedef QuocMaterialOptimizationConfiguratorBones< ConfiguratorType >                                                   MatOptConfigurator;
typedef typename DataTypeContainer::ParameterParserType                                                                 ParameterParserType;
typedef typename DataTypeContainer::PointType                                                                           PointType;
typedef typename ConfiguratorType::SparseMatrixType                                                                     SparseMatrixType;
typedef typename ConfiguratorType::MaskType                                                                             MaskType;
typedef typename ConfiguratorType::InitType                                                                             MeshType;

int main () {
    aol::consoleOutputStartProgramm( "Find Optimal Deformation (linear elasticity) for bone and polymer into all directions" );
    
    //! start watch
    auto startTime = std::chrono::high_resolution_clock::now();
    
    //! Parser
    ParameterParserType parser( "../../../../ParameterParser/shapeDesignBones/BonesAffinePeriodic.ini", "../../../../ParameterParser/counter.txt", "/OptDeformMultipleLoad" );
   
    //! Initialize mesh and configurator 
    IntVecType numDofVec; parser.template getFixSizeVector<int,IntVecType> ("InputMesh.NumDofVec", numDofVec );
    PointType lengthVec; parser.template getFixSizeVector<RealType, PointType> ("InputMesh.LengthVec", lengthVec );
    MeshType mesh ( numDofVec, lengthVec );
    
    ConfiguratorType conf ( mesh );
    MatOptConfigurator matOpConf ( parser, conf );
    QuocHandler<ConfiguratorType> quocHandler ( parser, conf );

    const int dimDomain = ConfiguratorType::dimDomain;
    const int numAffineSymGradDofs = conf.numAffineSymGradDofs;
    
    cout << endl << "We have " << endl 
         << mesh.getNumElements()  << " Elements" << endl
         << mesh.getNumVertices() << " NumVertices" << endl
         << "dim domain = " << dimDomain << endl 
         << "numAffineSymGradDofs = " << numAffineSymGradDofs << endl << endl;;
   
    //initialize material
    VectorType material ( conf.getNumGlobalDofs() );
    quocHandler.switchMaterialType( material );
    quocHandler.collabseVectorPeriodically( material );
    
    //read loads
    const int numLoads = parser.template get<int> ( "AffineDisp.numLoads" );
    std::vector<VectorType> affineDispBone, affineDispPolymer;
    for( int i=1; i <= numLoads; ++i ){
        VectorType affineDisp ( numAffineSymGradDofs );
        parser.template getFixSizeVector<RealType,VectorType> ( aol::strprintf("AffineDisp.Load%d", i ).c_str(), affineDisp );
        affineDispBone.push_back ( affineDisp ); affineDispPolymer.push_back ( affineDisp );
    }

    //! solve
    OptimalDeformSolverMultipleLoad<MatOptConfigurator> OptDeformFinder( parser, matOpConf, material, quocHandler, affineDispBone, affineDispPolymer );

    //PLOT
    VectorType materialScaledTo01( material.size() );
    for( int i=0; i<material.size(); i++ ) materialScaledTo01[i] = 0.5 * ( material[i] + 1.0 );
    quocHandler.extendVectorPeriodically( materialScaledTo01 );
   
    for( int i=0; i<numLoads; ++i ){
        quocHandler.plotPeriodicPlusAffineSymGradDispWithMaterial ( OptDeformFinder.template getSolutionDisplacementPeriodic<BONE>(i), OptDeformFinder.template getAffineDisplacement<BONE>(i),  materialScaledTo01, aol::strprintf( "DeformationBone_Direction%d", i).c_str() );
        quocHandler.plotPeriodicPlusAffineSymGradDispWithMaterial ( OptDeformFinder.template getSolutionDisplacementPeriodic<POLYMER>(i), OptDeformFinder.template getAffineDisplacement<POLYMER>(i), materialScaledTo01, aol::strprintf( "DeformationPolymer_Direction%d", i).c_str() );
    }
    
    //! print elapsed time
    std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << endl << "duration = " << diff.count() << " sec" << endl;

  return ( EXIT_SUCCESS );
}
