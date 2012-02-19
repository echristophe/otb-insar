IF(WIN32)
  ADD_DEFINITIONS(-DWIN32)
ENDIF(WIN32)


SET(BASELINE ${InSAR_DATA}/Baseline)
SET(INPUTDATA ${InSAR_DATA}/Input)

SET(EPSILON 0.01)

SET(BASELINE_TESTS ${CXX_TEST_PATH}/otbBaselineTests)

# -------------  otb::BaselineFunctorBase ----------------------------
ADD_TEST(blTuBaselineFunctorBaseNew ${BASELINE_TESTS}
otbBaselineFunctorBaseNew
)

							       
SET(BASELINE_SRCS
otbBaselineTests.cxx
otbBaselineFunctorBaseNew.cxx
)


ADD_EXECUTABLE(otbBaselineTests ${BASELINE_SRCS})
TARGET_LINK_LIBRARIES(otbBaselineTests OTBTesting OTBIO OTBCommon  ITKCommon Baseline)
