
#-----------------------------------------------------------------------------
set(BASELINE ${CMAKE_CURRENT_SOURCE_DIR}/../../Data/Baseline)
set(INPUT ${CMAKE_CURRENT_SOURCE_DIR}/../../Data/Input)
set(TEMP "${CMAKE_BINARY_DIR}/Testing/Temporary")

set(CLP ${MODULE_NAME})

#-----------------------------------------------------------------------------
add_executable(${CLP}Test ${CLP}Test.cxx)
target_link_libraries(${CLP}Test ${CLP}Lib ${SlicerExecutionModel_EXTRA_EXECUTABLE_TARGET_LIBRARIES})
set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})


#-----------------------------------------------------------------------------
set(testname ${CLP}Test1)
add_test(NAME ${testname} COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPoint
  --inputMesh ${BASELINE}/brainMesh.vtk
  --inputLabeledImage ${BASELINE}/brain.nii.gz
  --numIterations 8
  --flexibility 1
  --nonConnectivity 3
  #--outputMesh ${BASELINE}/brainMeshDeformed.vtk
   )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})

set(testname ${CLP}Test2)
add_test(NAME ${testname} COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPoint
  --inputMesh ${BASELINE}/nidusMesh.vtk
  --inputLabeledImage ${BASELINE}/nidus.nii.gz
  --numIterations 6
  --flexibility 1
  --nonConnectivity 2
  #--outputMesh ${BASELINE}/nidus/nidusMeshDeformed.vtk
   )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})

set(testname ${CLP}Test3)
add_test(NAME ${testname} COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPoint
  --inputMesh ${BASELINE}/ventriclesMesh.vtk
  --inputLabeledImage ${BASELINE}/ventricles.nii.gz
  --numIterations 7
  --flexibility 1
  --nonConnectivity 3
  #--outputMesh ${BASELINE}/ventriclesMeshDeformed.vtk
   )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
#ExternalData_add_target(${CLP}Data)

