#include <quocConfigurators.h>
#include <quocDefines.h>

#include "BonesOptDiffusionSolver.h"

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
    aol::consoleOutputStartProgramm( "Find Optimal Diffusion (linear) for bone and polymer into all directions" );
    
    //! start watch
    auto startTime = std::chrono::high_resolution_clock::now();
    
    //! Parser
    ParameterParserType parser( "../../../../ParameterParser/shapeDesignBones/BonesAffinePeriodic3D.ini", "../../../../ParameterParser/counter.txt", "/OptDiffusionMultipleLoad" );
   
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
         << "dim domain = " << dimDomain << endl  << endl;
   
    //initialize material
    VectorType material ( conf.getNumGlobalDofs() );
    quocHandler.switchMaterialType( material );
    quocHandler.collabseVectorPeriodically( material );
    
    //read loads
    std::vector<PointType> affineDiffusionBone, affineDiffusionPolymer;
    const int numDiffusionLoadsBone = parser.template get<int> ( "AffineDiffusion.numLoadsBone" );
    const int numDiffusionLoadsPolymer = parser.template get<int> ( "AffineDiffusion.numLoadsPolymer" );
    for( int i=1; i <= numDiffusionLoadsBone; ++i ){
                PointType affineDiffusion;
                parser.template getFixSizeVector<RealType,PointType> ( aol::strprintf("AffineDiffusion.LoadBone%d", i ).c_str(), affineDiffusion );
                affineDiffusionBone.push_back ( affineDiffusion );
    }
    for( int i=1; i <= numDiffusionLoadsPolymer; ++i ){
                PointType affineDiffusion;
                parser.template getFixSizeVector<RealType,PointType> ( aol::strprintf("AffineDiffusion.LoadPolymer%d", i ).c_str(), affineDiffusion );
                affineDiffusionPolymer.push_back ( affineDiffusion );
    }
      
    //! solve     
    OptimalDiffusionSolverMultipleLoad<MatOptConfigurator> OptDiffusionFinder ( parser, matOpConf, material, quocHandler, affineDiffusionBone, affineDiffusionPolymer );
    
    //! print elapsed time
    std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << endl << "duration = " << diff.count() << " sec" << endl;

  return ( EXIT_SUCCESS );
}
