#include <quocConfigurators.h>
#include <quocDefines.h>

#include "BonesMatOptSolver.h"

using namespace quocFE;
using namespace shapeOptBonePolymerPeriodicHomogenization;

#define _USEQUOC2D
// #define _USEQUOC3D

#ifdef _USEQUOC2D
typedef Quoc2DDataTypeContainer                                                                                         DataTypeContainer;
typedef QuocConfigurator2D<DataTypeContainer>                                                                           ConfiguratorType;
#endif

#ifdef _USEQUOC3D
typedef Quoc3DDataTypeContainer                                                                                         DataTypeContainer;
typedef QuocConfigurator3D<DataTypeContainer>                                                                           ConfiguratorType;
#endif

typedef typename DataTypeContainer::RealType                                                                            RealType;
typedef typename DataTypeContainer::IntVecType                                                                          IntVecType;
typedef typename DataTypeContainer::VectorType                                                                          VectorType;
typedef QuocMaterialOptimizationConfiguratorBones< ConfiguratorType >                                                   MatOptConfigurator;
typedef typename DataTypeContainer::ParameterParserType                                                                 ParameterParserType;
typedef typename DataTypeContainer::PointType                                                                           PointType;
typedef typename ConfiguratorType::SparseMatrixType                                                                     SparseMatrixType;
typedef typename ConfiguratorType::MaskType                                                                             MaskType;
typedef typename ConfiguratorType::InitType                                                                             MeshType;


//--------------------------------------------------------
int main(int argc, char ** argv) {
    
    aol::consoleOutputStartProgramm( "Phase field optimization for Bones" );
   
    ParameterParserType parser;
    cout << "Usage: programm [ParameterFile] [saveDirectory][onlyComputeInitialEnergy (optional, default false)] [int designType] [fileName design (if designType = 0)]" << endl;
    if (argc == 1) {
      parser = ParameterParserType ( "../../../../ParameterParser/shapeDesignBones/BonesAffinePeriodic2D.ini", "../../../../ParameterParser/counter.txt", "/MaterialOptBones2D" );
    }
    if( argc > 1 ) parser = ParameterParserType ( argv[1] );
    if( argc > 2 ) parser.set ( "saving.saveDirectory", argv[2] );
    bool onlyComputeInitialEnergy = false; 
    if( argc > 3 ) { 
        if( stoi(argv[3]) != 0 ) onlyComputeInitialEnergy = true; 
    }
    bool setDesignType = false; int designType; 
    if( argc > 4 ){ 
        setDesignType = true; 
        designType = stoi(argv[4]); 
    }
    bool readMaterialFromFile = false;
    string fileNameDesign;
    if( argc > 5 ){
        readMaterialFromFile = true;
        fileNameDesign = argv[5];
        parser.set<int> ( "MaterialOptimization.initMaterialType", 0 );
        parser.set<string> ("MaterialOptimization.materialFile", fileNameDesign );
    }
    
    //! start watch
    auto startTime = std::chrono::high_resolution_clock::now(); 
    
    //! Initialize mesh and configurator 
    IntVecType numDofVec; parser.template getFixSizeVector<int,IntVecType> ("InputMesh.NumDofVec", numDofVec );
    PointType lengthVec; parser.template getFixSizeVector<RealType, PointType> ("InputMesh.LengthVec", lengthVec );
    MeshType mesh ( numDofVec, lengthVec );
    
    ConfiguratorType conf ( mesh );
    QuocHandler<ConfiguratorType> quocHandler ( parser, conf );
    
    
    if( setDesignType ){
        cout << "set design from type" << endl;
        VectorType initMaterial( mesh.getNumVertices() );
        if( readMaterialFromFile ){
            cout << "read from file" << endl;
            quocHandler.switchMaterialType( initMaterial );
        }else{
            cout << "set from designNumber = " << designType << endl;
            const RealType volBoneMaterial = parser.template get<RealType> ("MaterialOptimization.VolInitialMaterial");
            string designTypeName;
            quocHandler.switchMaterialTypeForFixedVolume( designType, volBoneMaterial, initMaterial, designTypeName );
        }
        VectorType solutionMaterial( initMaterial );
        
        MaterialOptimizationMultipleLoadSolver<MatOptConfigurator> matOptProblem ( parser );
        const int refinementStep = 0;
        matOptProblem.computeOptimalDesign( mesh, initMaterial, solutionMaterial, refinementStep, onlyComputeInitialEnergy );
        
    }else{
        MaterialOptimizationMultipleLoadSolver<MatOptConfigurator> matOptProblem ( parser );
        matOptProblem.computeOptimalDesign_MultiLevel( mesh, parser.template get<int>("MaterialOptimization.numAdaptiveRefinementSteps") );
    }
    
    //! print elapsed time
    std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << endl << "duration = " << diff.count() << " sec" << endl;

    aol::consoleOutputFinishProgramm( "Phase field optimization for Bones" );
    
  return 0;
}

