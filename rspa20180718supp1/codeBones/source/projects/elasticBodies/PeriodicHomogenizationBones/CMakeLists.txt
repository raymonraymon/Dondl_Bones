#####################################################
# executables
#####################################################

# executables
add_executable(BonesOptDeform BonesOptDeform.cpp)
add_executable(BonesOptDiffusion BonesOptDiffusion.cpp)
add_executable(BonesMatOpt2D BonesMatOpt2D.cpp)
add_executable(BonesMatOpt3D BonesMatOpt3D.cpp)
add_executable(BashBonesCompareDesigns BashBonesCompareDesigns.cpp)
add_executable(BashBonesMatOpt BashBonesMatOpt.cpp)

# Link the utility library
set(ELASTBODIESPERIODICHOMOGENIZATIONBONESLIBS ${SUITESPARSE_LIBRARIES} ${Boost_LIBRARIES})

target_link_libraries(BonesOptDeform LINK_PUBLIC ${ELASTBODIESPERIODICHOMOGENIZATIONBONESLIBS})
target_link_libraries(BonesOptDiffusion LINK_PUBLIC ${ELASTBODIESPERIODICHOMOGENIZATIONBONESLIBS})
target_link_libraries(BonesMatOpt2D LINK_PUBLIC ${ELASTBODIESPERIODICHOMOGENIZATIONBONESLIBS})
target_link_libraries(BonesMatOpt3D LINK_PUBLIC ${ELASTBODIESPERIODICHOMOGENIZATIONBONESLIBS})
target_link_libraries(BashBonesMatOpt LINK_PUBLIC ${ELASTBODIESPERIODICHOMOGENIZATIONBONESLIBS})
target_link_libraries(BashBonesCompareDesigns LINK_PUBLIC ${ELASTBODIESPERIODICHOMOGENIZATIONBONESLIBS})


if( USE_IPOPT )
target_link_libraries(BonesMatOpt2D LINK_PUBLIC ${IPOPT_LIBRARIES})
target_link_libraries(BonesMatOpt3D LINK_PUBLIC ${IPOPT_LIBRARIES})
endif( USE_IPOPT )

if( USE_OPENMP )
 target_compile_definitions(BonesOptDeform PUBLIC USE_OPENMP)
 set_target_properties(BonesOptDeform PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
 target_compile_definitions(BonesOptDiffusion PUBLIC USE_OPENMP)
 set_target_properties(BonesOptDiffusion PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
 target_compile_definitions(BonesMatOpt2D PUBLIC USE_OPENMP)
 set_target_properties(BonesMatOpt2D PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
  target_compile_definitions(BonesMatOpt3D PUBLIC USE_OPENMP)
 set_target_properties(BonesMatOpt3D PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
endif( USE_OPENMP )

if( USE_VTK )
 target_link_libraries(BonesOptDeform LINK_PUBLIC ${VTK_LIBRARIES} )
 target_link_libraries(BonesOptDiffusion LINK_PUBLIC ${VTK_LIBRARIES} )
 target_link_libraries(BonesMatOpt2D LINK_PUBLIC ${VTK_LIBRARIES} )
  target_link_libraries(BonesMatOpt3D LINK_PUBLIC ${VTK_LIBRARIES} )
 target_link_libraries(BashBonesCompareDesigns LINK_PUBLIC ${VTK_LIBRARIES} )
endif( USE_VTK )
