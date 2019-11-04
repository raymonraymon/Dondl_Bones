#ifndef __BONESRESULTSPLOTTER_H
#define __BONESRESULTSPLOTTER_H

#include <quocConfigurators.h>
#include <quocHandler.h>
#include <quocDiscreteFunction.h>
#include <quocDefines.h>
#include <meshWithData.h>
#include <VTKInterfaceExtractor.h>
#include <VTKMeshWithData.h>
#include "BonesPlotter.h"
#include "BonesInterface.h"

using namespace quocFE;


namespace shapeOptBonePolymerPeriodicHomogenization{

template< typename MatOptConfigurator >
class BonesResultsPlotter {
    
private:
  
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfigurator::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
  typedef typename ConfiguratorType::DTContainer DataTypeContainer;
  typedef typename DataTypeContainer::IntVecType                                                                          IntVecType;
  typedef typename DataTypeContainer::ParameterParserType                                                                 ParameterParserType;
  typedef typename DataTypeContainer::PointType                                                                           PointType;
  typedef typename ConfiguratorType::SparseMatrixType                                                                     SparseMatrixType;
  typedef typename ConfiguratorType::MaskType                                                                             MaskType;
  typedef typename ConfiguratorType::InitType                                                                             MeshType;
  
  ParameterParserType &_parser;
  const string _saveDirectory;
  
  IntVecType _numDofVec; 
  PointType _lengthVec;  
  const int _numDofsOffsetOutside = 3;
  
public :
    
    BonesResultsPlotter ( ParameterParserType &parser ) :
    _parser ( parser ), _saveDirectory( parser.template get<string>( "saving.saveDirectory" ).c_str() ),
    _numDofsOffsetOutside( 3 ) {
        
        _parser.template getFixSizeVector<int,IntVecType> ("InputMesh.NumDofVec", _numDofVec );
        _parser.template getFixSizeVector<RealType, PointType> ("InputMesh.LengthVec", _lengthVec );
    }
    
        //=====================================================================================================================
        //Plot vtk files to png images
        //=====================================================================================================================
        void plotVTKToPNG( const int InitialOrSolution, const int plotVonMises ) const {
            aol::consoleOutput(  "plot vtk files to images" );
            const int numLoads = _parser.template get<int> ("AffineDisp.numLoads");
            VTKPlotterMaterialOptimizationBonesMultipleLoad<DataTypeContainer> vtkPlotter( _parser, numLoads );
            switch( InitialOrSolution ){
                case 1:{
                    vtkPlotter.plotAll("InitMaterial");  
                }break;
                case 2:{
                    vtkPlotter.plotAll("SolMaterial"); 
                    vtkPlotter.plotInterface( "SolMaterial" );
                }break;
                case 3:{
                    vtkPlotter.plotAll("InitMaterial");
                    vtkPlotter.plotAll("SolMaterial");
                    vtkPlotter.plotInterface( "SolMaterial" );
                    if( plotVonMises ) vtkPlotter.plotVonMisesStressesOnInterface( "SolMaterial", numLoads );
                }break;
                default :
                    throw std::invalid_argument( aol::strprintf ( "Wrong type. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
                    break;
            }
        }
        
        
        //! ========================================================================================
        //! plot material (B or P) around single micro cell
        //! ========================================================================================
        void plotSingleCell( ) const{
            
            MeshType mesh ( _numDofVec, _lengthVec );
            ConfiguratorType conf ( mesh );
            QuocHandler<ConfiguratorType> quocHandler( _parser, conf );
            
            VectorType material ( mesh.getNumVertices() );
            aol::loadVectorFromFile<VectorType>( material, aol::strprintf ( "%s/SolMaterial.txt", _saveDirectory.c_str()  ) ); 
            quocHandler.extendVectorPeriodically( material );
            QuocDiscreteFunctionDefault<ConfiguratorType> discreteFctSingleCell ( conf, material ); //note: here solMaterialBlock is actually OldLevel
        
            PointType offsetOutside; for( int i=0; i<offsetOutside.size(); ++i ) offsetOutside[i] = _numDofsOffsetOutside * mesh.getMeshSize( i );
            
            aol::consoleOutput(  "plot single cell with material outside" );
            IntVecType numDofVecOutside;
            PointType lengthVecOutside;
            for( int i=0; i<_numDofVec.size(); ++i ){
                numDofVecOutside[i] = _numDofVec[i] + 2 * _numDofsOffsetOutside;
                lengthVecOutside[i] = _lengthVec[i] + 2 * _numDofsOffsetOutside * mesh.getMeshSize( i ); 
            }
            MeshType meshOutside ( numDofVecOutside, lengthVecOutside );
            VectorType materialOutsideBone ( meshOutside.getNumVertices() );
            for( int i=0; i<materialOutsideBone.size(); ++i ) materialOutsideBone[i] = 1.0;
            VectorType materialOutsidePolymer ( meshOutside.getNumVertices() );
            for( int i=0; i<materialOutsidePolymer.size(); ++i ) materialOutsidePolymer[i] = -1.0;
            for( int nodeIdxOffset=0; nodeIdxOffset < meshOutside.getNumVertices(); nodeIdxOffset++ ){
                const PointType& GlobalCoordsOffset = meshOutside.getVertex ( nodeIdxOffset );
               
                bool inside = true;
                for( int i=0; i<GlobalCoordsOffset.size(); ++i ){
                  if( GlobalCoordsOffset[i] < offsetOutside[i] ) inside = false;
                  if( GlobalCoordsOffset[i] > _lengthVec[i] + offsetOutside[i] ) inside = false;
                }
                
                if( inside ) {
                    PointType GlobalCoordsSingleCell;
                    for( int i=0; i<GlobalCoordsSingleCell.size(); ++i){
                        GlobalCoordsSingleCell[i] = GlobalCoordsOffset[i] - _numDofsOffsetOutside * mesh.getMeshSize(i);
                    }
                    int elementNumberSingleCell; PointType LocalCoordSingleCell;
                    conf.getLocalCoords ( GlobalCoordsSingleCell, elementNumberSingleCell, LocalCoordSingleCell );
                    materialOutsideBone[nodeIdxOffset] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
                    materialOutsidePolymer[nodeIdxOffset] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
                }
            }
            
//             VTKMeshWithData<MeshType> meshSaverOutsideBone ( meshOutside );
            MeshWithData<MeshType> meshSaverOutsideBone ( meshOutside );
            meshSaverOutsideBone.addScalarData ( materialOutsideBone, "material", VERTEX_DATA ); //"SolMaterialOutsideBone"
            meshSaverOutsideBone.saveAsVTKSTRUCTUREDPOINTS ( aol::strprintf( "%s/SolMaterial_OutsideBone.vtk", _saveDirectory.c_str() ), offsetOutside );
//             meshSaverOutsideBone.save ( aol::strprintf( "%s/SolMaterial_OutsideBone.vtk", _saveDirectory.c_str() ), offsetOutside, VTKSTRUCTUREDPOINTS );
            
            MeshWithData<MeshType> meshSaverOutsidePolymer ( meshOutside );
//             VTKMeshWithData<MeshType> meshSaverOutsidePolymer ( meshOutside );
            meshSaverOutsidePolymer.addScalarData ( materialOutsidePolymer, "material", VERTEX_DATA ); //"SolMaterialOutsidePolymer"
            meshSaverOutsidePolymer.saveAsVTKSTRUCTUREDPOINTS ( aol::strprintf( "%s/SolMaterial_OutsidePolymer.vtk", _saveDirectory.c_str() ), offsetOutside );
//             meshSaverOutsidePolymer.save ( aol::strprintf( "%s/SolMaterial_OutsidePolymer.vtk", _saveDirectory.c_str() ), offsetOutside, VTKSTRUCTUREDPOINTS );
        
        }
        
        
        //=====================================================================================================================
        //Plot block of mesh with material 
        //=====================================================================================================================
        void plotBlockOfCells( const int numBlocksPerDirection ) const{
            aol::consoleOutput( "plot block of cells" );
            
            MeshType mesh ( _numDofVec, _lengthVec );
            ConfiguratorType conf ( mesh );
            QuocHandler<ConfiguratorType> quocHandler( _parser, conf );
            
            VectorType material ( mesh.getNumVertices() );
            aol::loadVectorFromFile<VectorType>( material, aol::strprintf ( "%s/SolMaterial.txt", _saveDirectory.c_str()  ) ); 
            quocHandler.extendVectorPeriodically( material );
            QuocDiscreteFunctionDefault<ConfiguratorType> discreteFctSingleCell ( conf, material ); //note: here solMaterialBlock is actually OldLevel
        
            PointType offsetOutside; for( int i=0; i<offsetOutside.size(); ++i ) offsetOutside[i] = _numDofsOffsetOutside * mesh.getMeshSize( i );
            
            //! ========================================================================================
            //! plot block
            //! ========================================================================================
            IntVecType _numDofVecBlock;
            PointType _lengthVecBlock;
            for( int i=0; i<_numDofVec.size(); ++i ){
                int oldSize = _numDofVec[i];
                _numDofVecBlock[i] = static_cast<int> ( numBlocksPerDirection * (oldSize - 1) + 1 );
                _lengthVecBlock[i] = static_cast<RealType> ( numBlocksPerDirection ) * _lengthVec[i]; 
            }
            MeshType meshBlock ( _numDofVecBlock, _lengthVecBlock );
            VectorType materialBlock ( meshBlock.getNumVertices() );
            for( int nodeIdxBlock=0; nodeIdxBlock < meshBlock.getNumVertices(); nodeIdxBlock++ ){
                const PointType& GlobalCoordsBlock = meshBlock.getVertex ( nodeIdxBlock );
                PointType GlobalCoordsSingleCell;
                for( int i=0; i<GlobalCoordsSingleCell.size(); ++i){
                    int multDirection = static_cast<int> ( GlobalCoordsBlock[i] / _lengthVec[i] );
                    GlobalCoordsSingleCell[i] = GlobalCoordsBlock[i] - multDirection * _lengthVec[i];
                }
                int elementNumberSingleCell; PointType LocalCoordSingleCell;
                conf.getLocalCoords ( GlobalCoordsSingleCell, elementNumberSingleCell, LocalCoordSingleCell );
                materialBlock[nodeIdxBlock] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
            }
            MeshWithData<MeshType> meshSaver ( meshBlock );
//             VTKMeshWithData<MeshType> meshSaver ( meshBlock );
            meshSaver.addScalarData ( materialBlock, "material", VERTEX_DATA ); //"SolMaterialBlock"
            meshSaver.saveAsVTKSTRUCTUREDPOINTS ( aol::strprintf( "%s/SolMaterial_Block.vtk", _saveDirectory.c_str() ) );
//             meshSaver.save ( aol::strprintf( "%s/SolMaterial_Block.vtk", _saveDirectory.c_str() ), VTKSTRUCTUREDPOINTS );
            
            //! ========================================================================================
            //! plot material (B or P) around block of micro cells
            //! ========================================================================================
            IntVecType _numDofVecBlockOutside;
            PointType _lengthVecBlockOutside;
            for( int i=0; i<_numDofVec.size(); ++i ){
                _numDofVecBlockOutside[i] = _numDofVecBlock[i] + 2 * _numDofsOffsetOutside;
                _lengthVecBlockOutside[i] = _lengthVecBlock[i] + 2 * _numDofsOffsetOutside * mesh.getMeshSize( i ); 
            }
            MeshType meshBlockOutside ( _numDofVecBlockOutside, _lengthVecBlockOutside );
            VectorType materialBlockOutsideBone ( meshBlockOutside.getNumVertices() );
            for( int i=0; i<materialBlockOutsideBone.size(); ++i ) materialBlockOutsideBone[i] = 1.0;
            VectorType materialBlockOutsidePolymer ( meshBlockOutside.getNumVertices() );
            for( int i=0; i<materialBlockOutsidePolymer.size(); ++i ) materialBlockOutsidePolymer[i] = -1.0;
            for( int nodeIdxBlock=0; nodeIdxBlock < meshBlockOutside.getNumVertices(); nodeIdxBlock++ ){
                const PointType& GlobalCoordsBlockOffset = meshBlockOutside.getVertex ( nodeIdxBlock );
                
                bool inside = true;
                for( int i=0; i<GlobalCoordsBlockOffset.size(); ++i ){
                  if( GlobalCoordsBlockOffset[i] < offsetOutside[i] ) inside = false;
                  if( GlobalCoordsBlockOffset[i] > _lengthVecBlock[i] + offsetOutside[i] ) inside = false;
                }
                
                if( inside ) {
                    
                    PointType GlobalCoordsBlock;
                    for( int i=0; i<GlobalCoordsBlock.size(); ++i){
                        GlobalCoordsBlock[i] = GlobalCoordsBlockOffset[i] - _numDofsOffsetOutside * mesh.getMeshSize(i);
                    }
                    PointType GlobalCoordsSingleCell;
                    for( int i=0; i<GlobalCoordsSingleCell.size(); ++i){
                      int multDirection = static_cast<int> ( GlobalCoordsBlock[i] / _lengthVec[i] );
                      GlobalCoordsSingleCell[i] = GlobalCoordsBlock[i] - multDirection * _lengthVec[i];
                    }
                    int elementNumberSingleCell; PointType LocalCoordSingleCell;
                    conf.getLocalCoords ( GlobalCoordsSingleCell, elementNumberSingleCell, LocalCoordSingleCell );
                    materialBlockOutsideBone[nodeIdxBlock] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
                    materialBlockOutsidePolymer[nodeIdxBlock] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
                }
            }
            
            MeshWithData<MeshType> meshSaverBlockOutsideBone ( meshBlockOutside );
            meshSaverBlockOutsideBone.addScalarData ( materialBlockOutsideBone, "material", VERTEX_DATA ); //"SolMaterialBlockOutsideBone"
            meshSaverBlockOutsideBone.saveAsVTKSTRUCTUREDPOINTS ( aol::strprintf( "%s/SolMaterial_BlockOutsideBone.vtk", _saveDirectory.c_str() ), offsetOutside );

            MeshWithData<MeshType> meshSaverBlockOutsidePolymer ( meshBlockOutside );
            meshSaverBlockOutsidePolymer.addScalarData ( materialBlockOutsidePolymer, "material", VERTEX_DATA ); //"SolMaterialBlockOutsidePolymer"
            meshSaverBlockOutsidePolymer.saveAsVTKSTRUCTUREDPOINTS ( aol::strprintf( "%s/SolMaterial_BlockOutsidePolymer.vtk", _saveDirectory.c_str() ), offsetOutside );
        }
        
        //=====================================================================================================================
        //plot interface
        //=====================================================================================================================
        void plotInterface( const bool plotVonMises = true ) const {
            aol::consoleOutput( "plot interface" );
            _parser.createSubDirectory("Interface");
            
            MeshType mesh ( _numDofVec, _lengthVec );
            ConfiguratorType conf ( mesh );
            QuocHandler<ConfiguratorType> quocHandler( _parser, conf );
            
            VectorType material ( mesh.getNumVertices() );
            aol::loadVectorFromFile<VectorType>( material, aol::strprintf ( "%s/SolMaterial.txt", _saveDirectory.c_str()  ) ); 
            quocHandler.extendVectorPeriodically( material );
            QuocDiscreteFunctionDefault<ConfiguratorType> discreteFctSingleCell ( conf, material ); //note: here solMaterialBlock is actually OldLevel
        
            PointType offsetOutside; for( int i=0; i<offsetOutside.size(); ++i ) offsetOutside[i] = _numDofsOffsetOutside * mesh.getMeshSize( i );
            
            //extract surface
            vtkInterfaceExtractor surfaceExtractor;
            const RealType threshold = _parser.template get<RealType> ( "saving.thresholdInterface" );

            surfaceExtractor.getInterfaceBonePolymer( aol::strprintf( "%s/SolMaterial_OutsideBone.vtk", _saveDirectory.c_str() ),
                                                      aol::strprintf( "%s/Interface/SolMaterial_InterfacePolymer.vtk", _saveDirectory.c_str() ), "material", threshold );
            surfaceExtractor.getInterfaceBonePolymer( aol::strprintf( "%s/SolMaterial_OutsidePolymer.vtk", _saveDirectory.c_str() ),
                                         aol::strprintf( "%s/Interface/SolMaterial_InterfaceBone.vtk", _saveDirectory.c_str() ), "material", threshold );
            surfaceExtractor.getInterfaceBonePolymer( aol::strprintf( "%s/SolMaterial_BlockOutsideBone.vtk", _saveDirectory.c_str() ),
                                         aol::strprintf( "%s/Interface/SolMaterial_BlockInterfacePolymer.vtk", _saveDirectory.c_str() ), "material", threshold );
            surfaceExtractor.getInterfaceBonePolymer( aol::strprintf( "%s/SolMaterial_BlockOutsidePolymer.vtk", _saveDirectory.c_str() ),
                                         aol::strprintf( "%s/Interface/SolMaterial_BlockInterfaceBone.vtk", _saveDirectory.c_str() ), "material", threshold );
            
            //use loop subdivision for interfaces
            const int numSubdivLevels = 2;
            surfaceExtractor.loopSubdivision( aol::strprintf( "%s/Interface/SolMaterial_InterfacePolymer.vtk", _saveDirectory.c_str() ),
                                              aol::strprintf( "%s/Interface/SolMaterial_InterfacePolymer_Subdiv%d.vtk", _saveDirectory.c_str(), numSubdivLevels ),
                                              numSubdivLevels );
            surfaceExtractor.loopSubdivision( aol::strprintf( "%s/Interface/SolMaterial_InterfaceBone.vtk", _saveDirectory.c_str() ),
                                              aol::strprintf( "%s/Interface/SolMaterial_InterfaceBone_Subdiv%d.vtk", _saveDirectory.c_str(), numSubdivLevels ),
                                              numSubdivLevels );
            surfaceExtractor.loopSubdivision( aol::strprintf( "%s/Interface/SolMaterial_BlockInterfacePolymer.vtk", _saveDirectory.c_str() ),
                                              aol::strprintf( "%s/Interface/SolMaterial_BlockInterfacePolymer_Subdiv%d.vtk", _saveDirectory.c_str(), numSubdivLevels ),
                                              numSubdivLevels );
            surfaceExtractor.loopSubdivision( aol::strprintf( "%s/Interface/SolMaterial_BlockInterfaceBone.vtk", _saveDirectory.c_str() ),
                                              aol::strprintf( "%s/Interface/SolMaterial_BlockInterfaceBone_Subdiv%d.vtk", _saveDirectory.c_str(), numSubdivLevels ),
                                              numSubdivLevels );
                                              
            
            //! ========================================================================================
            //! plot von mises stresses on interface
            //! ========================================================================================
            if( plotVonMises ){
            aol::consoleOutput( "plot von mises stresses on interface" );
            _parser.createSubDirectory("StressOnInterface");
            
            MatOptConfigurator matOpConf ( _parser, conf );
            const int numLoads = _parser.template get<int> ( "AffineDisp.numLoads" );
            
            // for bone
            for( int loadIdx=0; loadIdx<numLoads; ++loadIdx ){
                
                BonesInterface<MatOptConfigurator,BONE> interface ( matOpConf );
                
                VectorType displacementPeriodic ( conf.dimDomain * conf.getNumGlobalDofs() ), displacementPeriodicExtended ( conf.dimDomain * conf.getNumGlobalDofs() );
                aol::loadVectorFromFile<VectorType>( displacementPeriodic, aol::strprintf( "%s/Deformation/SolMaterial_DisplacementBonePeriodic_Dir%d.txt", _saveDirectory.c_str(),  loadIdx ) ); 
                quocHandler.extendMultiVectorPeriodically( displacementPeriodic, displacementPeriodicExtended );
                
                VectorType displacementAffine ( conf.numAffineSymGradDofs );
                aol::loadVectorFromFile<VectorType>( displacementAffine, aol::strprintf( "%s/Deformation/SolMaterial_DisplacementBoneAffine_Dir%d.txt", _saveDirectory.c_str(),  loadIdx ) ); 
                
                std::vector<RealType> stressVec, stressVecBlock;
                interface.evaluateStress( aol::strprintf( "%s/Interface/SolMaterial_InterfaceBone.vtk", _saveDirectory.c_str() ).c_str(), material, displacementPeriodicExtended, displacementAffine, stressVec, false  );
                interface.evaluateStress( aol::strprintf( "%s/Interface/SolMaterial_BlockInterfaceBone.vtk", _saveDirectory.c_str() ).c_str(), material, displacementPeriodicExtended, displacementAffine, stressVecBlock, true  );
                
                interface.plotStress( aol::strprintf( "%s/Interface/SolMaterial_InterfaceBone.vtk", _saveDirectory.c_str() ).c_str(), 
                                      aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_InterfaceBone_Direction%d.vtp", _saveDirectory.c_str(), loadIdx ).c_str(),
                                      "VonMisesStresses", stressVec );
                interface.plotStress( aol::strprintf( "%s/Interface/SolMaterial_BlockInterfaceBone.vtk", _saveDirectory.c_str() ).c_str(), 
                                      aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_BlockInterfaceBone_Direction%d.vtp", _saveDirectory.c_str(), loadIdx ).c_str(),
                                      "VonMisesStresses", stressVecBlock );
            }
            // for polymer
            for( int loadIdx=0; loadIdx<numLoads; ++loadIdx ){
                
                BonesInterface<MatOptConfigurator,POLYMER> interface ( matOpConf );
                
                VectorType displacementPeriodic ( conf.dimDomain * conf.getNumGlobalDofs() ), displacementPeriodicExtended ( conf.dimDomain * conf.getNumGlobalDofs() );
                aol::loadVectorFromFile<VectorType>( displacementPeriodic, aol::strprintf( "%s/Deformation/SolMaterial_DisplacementPolymerPeriodic_Dir%d.txt", _saveDirectory.c_str(),  loadIdx ) ); 
                quocHandler.extendMultiVectorPeriodically( displacementPeriodic, displacementPeriodicExtended );
                
                VectorType displacementAffine ( conf.numAffineSymGradDofs );
                aol::loadVectorFromFile<VectorType>( displacementAffine, aol::strprintf( "%s/Deformation/SolMaterial_DisplacementPolymerAffine_Dir%d.txt", _saveDirectory.c_str(),  loadIdx ) ); 
                
                std::vector<RealType> stressVec, stressVecBlock;
                interface.evaluateStress( aol::strprintf( "%s/Interface/SolMaterial_InterfacePolymer.vtk", _saveDirectory.c_str() ).c_str(), material, displacementPeriodicExtended, displacementAffine, stressVec, false  );
                interface.evaluateStress( aol::strprintf( "%s/Interface/SolMaterial_BlockInterfacePolymer.vtk", _saveDirectory.c_str() ).c_str(), material, displacementPeriodicExtended, displacementAffine, stressVecBlock, true  );
                
                interface.plotStress( aol::strprintf( "%s/Interface/SolMaterial_InterfacePolymer.vtk", _saveDirectory.c_str() ).c_str(), 
                                      aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_InterfacePolymer_Direction%d.vtp", _saveDirectory.c_str(), loadIdx ).c_str(),
                                      "VonMisesStresses", stressVec );
                interface.plotStress( aol::strprintf( "%s/Interface/SolMaterial_BlockInterfacePolymer.vtk", _saveDirectory.c_str() ).c_str(), 
                                      aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_BlockInterfacePolymer_Direction%d.vtp", _saveDirectory.c_str(), loadIdx ).c_str(),
                                      "VonMisesStresses", stressVecBlock );
            }
            }
            
        }
        
        //=====================================================================================================================
        //plot von mises stresses
        //=====================================================================================================================
        
        void plotVonMisesStresses( ) const{
            aol::consoleOutput( "plot von mises stresses" );
            _parser.createSubDirectory("Stress");
            
            MeshType mesh ( _numDofVec, _lengthVec );
            ConfiguratorType conf ( mesh );
            QuocHandler<ConfiguratorType> quocHandler( _parser, conf );
            
            VectorType material ( mesh.getNumVertices() );
            aol::loadVectorFromFile<VectorType>( material, aol::strprintf ( "%s/SolMaterial.txt", _saveDirectory.c_str()  ) ); 
            quocHandler.extendVectorPeriodically( material );
            QuocDiscreteFunctionDefault<ConfiguratorType> discreteFctSingleCell ( conf, material ); //note: here solMaterialBlock is actually OldLevel
            MatOptConfigurator matOpConf ( _parser, conf );
            const int numLoads = _parser.template get<int> ( "AffineDisp.numLoads" );
            
            //for bone
            for( int loadIdx=0; loadIdx<numLoads; ++loadIdx ){
                
                BonesInterface<MatOptConfigurator,BONE> interface ( matOpConf );
                
                VectorType displacementPeriodic ( conf.dimDomain * conf.getNumGlobalDofs() ), displacementPeriodicExtended ( conf.dimDomain * conf.getNumGlobalDofs() );
                aol::loadVectorFromFile<VectorType>( displacementPeriodic, aol::strprintf( "%s/Deformation/SolMaterial_DisplacementBonePeriodic_Dir%d.txt", _saveDirectory.c_str(),  loadIdx ) ); 
                quocHandler.extendMultiVectorPeriodically( displacementPeriodic, displacementPeriodicExtended );
                
                VectorType displacementAffine ( conf.numAffineSymGradDofs );
                aol::loadVectorFromFile<VectorType>( displacementAffine, aol::strprintf( "%s/Deformation/SolMaterial_DisplacementBoneAffine_Dir%d.txt", _saveDirectory.c_str(),  loadIdx ) ); 
                
                std::vector<RealType> stressVec3d, stressVecBlock3d;
                interface.evaluateStress( aol::strprintf( "%s/SolMaterial_Undeformed.vtk", _saveDirectory.c_str() ).c_str(), material, displacementPeriodicExtended, displacementAffine, stressVec3d, false  );
                interface.evaluateStress( aol::strprintf( "%s/SolMaterial_Block.vtk", _saveDirectory.c_str() ).c_str(), material, displacementPeriodicExtended, displacementAffine, stressVecBlock3d, true );
                interface.plotStress( aol::strprintf( "%s/SolMaterial_Undeformed.vtk", _saveDirectory.c_str() ).c_str(), 
                                      aol::strprintf( "%s/Stress/SolMaterial_VonMises_Bone_Direction%d.vtk", _saveDirectory.c_str(), loadIdx ).c_str(),
                                      "VonMisesStresses", stressVec3d );
                interface.plotStress( aol::strprintf( "%s/SolMaterial_Block.vtk", _saveDirectory.c_str() ).c_str(), 
                                      aol::strprintf( "%s/Stress/SolMaterial_VonMises_BlockBone_Direction%d.vtk", _saveDirectory.c_str(), loadIdx ).c_str(),
                                      "VonMisesStresses", stressVecBlock3d );
            }
            
            //for polymer
            for( int loadIdx=0; loadIdx<numLoads; ++loadIdx ){
                
                BonesInterface<MatOptConfigurator,POLYMER> interface ( matOpConf );
                
                VectorType displacementPeriodic ( conf.dimDomain * conf.getNumGlobalDofs() ),displacementPeriodicExtended ( conf.dimDomain * conf.getNumGlobalDofs() ) ;
                aol::loadVectorFromFile<VectorType>( displacementPeriodic, aol::strprintf( "%s/Deformation/SolMaterial_DisplacementPolymerPeriodic_Dir%d.txt", _saveDirectory.c_str(),  loadIdx ) ); 
                quocHandler.extendMultiVectorPeriodically( displacementPeriodic, displacementPeriodicExtended );
                
                VectorType displacementAffine ( conf.numAffineSymGradDofs );
                aol::loadVectorFromFile<VectorType>( displacementAffine, aol::strprintf( "%s/Deformation/SolMaterial_DisplacementPolymerAffine_Dir%d.txt", _saveDirectory.c_str(),  loadIdx ) ); 
                
                std::vector<RealType> stressVec3d, stressVecBlock3d;
                interface.evaluateStress( aol::strprintf( "%s/SolMaterial_Undeformed.vtk", _saveDirectory.c_str() ).c_str(), material, displacementPeriodicExtended, displacementAffine, stressVec3d, false  );
                interface.evaluateStress( aol::strprintf( "%s/SolMaterial_Block.vtk", _saveDirectory.c_str() ).c_str(), material, displacementPeriodicExtended, displacementAffine, stressVecBlock3d, true );
                interface.plotStress( aol::strprintf( "%s/SolMaterial_Undeformed.vtk", _saveDirectory.c_str() ).c_str(), 
                                      aol::strprintf( "%s/Stress/SolMaterial_VonMises_Polymer_Direction%d.vtk", _saveDirectory.c_str(), loadIdx ).c_str(),
                                      "VonMisesStresses", stressVec3d );
                interface.plotStress( aol::strprintf( "%s/SolMaterial_Block.vtk", _saveDirectory.c_str() ).c_str(), 
                                      aol::strprintf( "%s/Stress/SolMaterial_VonMises_BlockPolymer_Direction%d.vtk", _saveDirectory.c_str(), loadIdx ).c_str(),
                                      "VonMisesStresses", stressVecBlock3d );
            }
            
        }
        
        
        //=====================================================================================================================
        //compute homogenized tensor
        //=====================================================================================================================
        void computeHomogenizedTensor( ) const {
            
            aol::consoleOutput( "compute homogenized tensor" );
            
            _parser.createSubDirectory("HomogenizedTensor");
            
            
            MeshType mesh ( _numDofVec, _lengthVec );
            ConfiguratorType conf ( mesh );
            QuocHandler<ConfiguratorType> quocHandler( _parser, conf );
            
            VectorType material ( mesh.getNumVertices() );
            aol::loadVectorFromFile<VectorType>( material, aol::strprintf ( "%s/SolMaterial.txt", _saveDirectory.c_str()  ) ); 
            quocHandler.extendVectorPeriodically( material );
            QuocDiscreteFunctionDefault<ConfiguratorType> discreteFctSingleCell ( conf, material ); //note: here solMaterialBlock is actually OldLevel
            
            const int dimDomain = ConfiguratorType::dimDomain;
            const int numAffineSymGradDofs = conf.numAffineSymGradDofs;
                       
            //define loads
            int numLoadsAll;
            if( dimDomain == 2) numLoadsAll = 2;
            if( dimDomain == 3 ) numLoadsAll = 6;
            std::vector<VectorType> affineDispBone, affineDispPolymer;
            for( int i=0; i < numLoadsAll; ++i ){
                VectorType affineDisp ( numAffineSymGradDofs ); affineDisp.setZero();
                affineDisp[i] = 0.1;                
                affineDispBone.push_back ( affineDisp ); affineDispPolymer.push_back ( affineDisp );
            }
            
            //!determine weightFunction for Multiple Load
            std::vector<RealType> weightFct_Bone_weightVec ( numLoadsAll ), weightFct_Polymer_weightVec ( numLoadsAll );
            for( int i=0; i<numLoadsAll; ++i){
              weightFct_Bone_weightVec[i] = 1.; weightFct_Polymer_weightVec[i] = 1.;   
            }
            
            RealType eps_area = 1.0, eps_factor = 1.0;
            MatOptConfigurator matOpConf ( _parser, conf, eps_area, eps_factor, weightFct_Bone_weightVec, weightFct_Polymer_weightVec );
            
            //! solve
            OptimalDeformSolverMultipleLoad<MatOptConfigurator> OptDeformFinderAll( _parser, matOpConf, material, quocHandler, affineDispBone, affineDispPolymer );

            ComplianceEnergyMultipleLoad_RegularizedMaxFunction<MatOptConfigurator> complianceOp ( OptDeformFinderAll );
            RealType complianceEnergy;
            complianceOp.apply( material, complianceEnergy );
            
            std::vector<RealType> homogenizedTensor_Bone ( numLoadsAll ), homogenizedTensor_Polymer ( numLoadsAll );
            for( int i=0; i<numLoadsAll; ++i ){
                homogenizedTensor_Bone[i] = 100. * complianceOp.getLastComplianceBone(i);
                homogenizedTensor_Polymer[i] = 100. * complianceOp.getLastCompliancePolymer(i);
            }
            
            
            //compute volume 
            shapeOptBonePolymerPeriodicHomogenization::VolumeConstraintPeriodicBC<MatOptConfigurator,BONE> volumeOpBone ( matOpConf, quocHandler );
            shapeOptBonePolymerPeriodicHomogenization::VolumeConstraintPeriodicBC<MatOptConfigurator,POLYMER> volumeOpPolymer ( matOpConf, quocHandler );
            RealType volumeBone; volumeOpBone.evaluateEnergy( material, volumeBone );
            RealType volumePolymer; volumeOpPolymer.evaluateEnergy( material, volumePolymer );
            
            MaterialOptimizationMultipleLoadEnergyInfo<MatOptConfigurator> energyInfoAll;
            energyInfoAll.setNumLoads( numLoadsAll );
            energyInfoAll.setComplianceBone( homogenizedTensor_Bone );
            energyInfoAll.setCompliancePolymer( homogenizedTensor_Polymer );
            energyInfoAll.setVolumeBone( volumeBone );
            energyInfoAll.setVolumePolymer( volumePolymer );
            
            TikzPlotterMaterialOptimizationBonesMultipleLoad<DataTypeContainer> tikzPlotterAll ( _parser );
            tikzPlotterAll.template generateTexCodeForPaper<MatOptConfigurator>( "HomogenizedTensor/codeForPaper", energyInfoAll );

        }
        
        
        //=====================================================================================================================
        //Generate PDF with all results
        //=====================================================================================================================
        void summarizeAllResults( ) const {
            std::ofstream BashFilePDF ( aol::strprintf ( "%s/%s.sh", _saveDirectory.c_str (), "generatePDF"  ) );
            BashFilePDF << "cd " << _saveDirectory.c_str() << endl;
            BashFilePDF << "pdflatex Results" << ".tex" << endl;
            string systemCommandPdfLatex;
            systemCommandPdfLatex += "nohup "; 
            systemCommandPdfLatex += "bash ";
            systemCommandPdfLatex += _saveDirectory + "/generatePDF.sh";
            systemCommandPdfLatex += " >" + _saveDirectory + "/ConsoleOutputPdfLatex.out 2>&1 ";
            cout << "systemCommandPdfLatex = " << systemCommandPdfLatex << endl;
            bool failed; failed = ( system ( systemCommandPdfLatex.c_str() ) != EXIT_SUCCESS );
            if ( failed ) cerr << "programm returned an error." << endl;
        }
        
    
    
};


}//end namespace

#endif
