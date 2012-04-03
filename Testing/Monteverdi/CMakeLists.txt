IF(WIN32)
  ADD_DEFINITIONS(-DWIN32)
ENDIF(WIN32)


SET(BASELINE ${InSAR_DATA}/Baseline)
SET(INPUTDATA ${InSAR_DATA}/Input)

IF(INSAR_DATA_USE_LARGEINPUT)
  SET(LARGEINPUT ${INSAR_DATA_LARGEINPUT_ROOT} )
ENDIF(INSAR_DATA_USE_LARGEINPUT)

SET(InSARModules_TESTS1 ${CXX_TEST_PATH}/otbInSARModulesTests1)

#----------------------------------------------------------------------------------------------
#------------------------------------ otbModulesTests1 ----------------------------------------
#----------------------------------------------------------------------------------------------

#-------------- Dynamic Fringe Viewer Module Tests ---------------------------------
ADD_TEST(moTvDynamicFringeViewerModuleTest ${InSARModules_TESTS1} 
otbDynamicFringeViewerModuleTest
${INPUTDATA}/fringe.png
0
)


SET(InSARModules_SRCS1
otbInSARModulesTests1.cxx
otbDynamicFringeViewerModuleTest.cxx
)


ADD_EXECUTABLE(otbInSARModulesTests1 otbInSARModulesTests1.cxx ${InSARModules_SRCS1})
TARGET_LINK_LIBRARIES(otbInSARModulesTests1 ${LIST_MODULE_LIBRARIES} OTBIO OTBCommon OTBTesting otbMonteverdiCommon otbDynamicFringeViewerModule)

