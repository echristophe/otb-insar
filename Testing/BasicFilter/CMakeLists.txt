IF(WIN32)
  ADD_DEFINITIONS(-DWIN32)
ENDIF(WIN32)


SET(BASELINE ${InSAR_DATA}/Baseline)
SET(INPUTDATA ${InSAR_DATA}/Input)

IF(INSAR_DATA_USE_LARGEINPUT)
  SET(LARGEINPUT ${INSAR_DATA_LARGEINPUT_ROOT} )
ENDIF(INSAR_DATA_USE_LARGEINPUT)

SET(EPSILON 0.01)

SET(BASICFILTER_TESTS ${CXX_TEST_PATH}/otbBasicFilterTests)

# -------------  otb::ModuloImageFilter ----------------------------
ADD_TEST(bfTvModuloImageFilterTest ${BASICFILTER_TESTS}
otbModuloImageFilterTest
)


SET(BASICFILTER_SRCS
otbBasicFilterTests.cxx
otbModuloImageFilterTest.cxx
)


ADD_EXECUTABLE(otbBasicFilterTests ${BASICFILTER_SRCS})
TARGET_LINK_LIBRARIES(otbBasicFilterTests OTBTesting OTBIO OTBCommon  ITKCommon BasicFilter)
