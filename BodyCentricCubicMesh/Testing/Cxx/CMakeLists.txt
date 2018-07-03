
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
  --inputLabeledImage ${BASELINE}/brain.nii.gz
  --size 8  
  --fidelity 0.9
  --resampleResolution 0
  #--outputMesh ${BASELINE}/brainMesh.vtk
   )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})

set(testname ${CLP}Test2)
add_test(NAME ${testname} COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPoint
  --inputLabeledImage ${BASELINE}/nidus.nii.gz
  --size 9  
  --fidelity 0.9
  --resampleResolution 0
  #--outputMesh ${BASELINE}/nidusMesh.vtk
   )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})

set(testname ${CLP}Test3)
add_test(NAME ${testname} COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPoint
  --inputLabeledImage ${BASELINE}/ventricles.nii.gz
  --size 7  
  --fidelity 0.9
  --resampleResolution 0
  #--outputMesh ${BASELINE}/ventriclesMesh.vtk
   )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})


#-----------------------------------------------------------------------------
#ExternalData_add_target(${CLP}Data)

