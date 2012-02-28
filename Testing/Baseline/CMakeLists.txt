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

# -------------  otb::HorizontalVerticalBaselineFunctor---------------
ADD_TEST(blTuHorizontalVerticalBaselineFunctorNew ${BASELINE_TESTS}
otbHorizontalVerticalBaselineFunctorNew
)

ADD_TEST(blTuHorizontalVerticalBaselineFunctorTest ${BASELINE_TESTS}
otbHorizontalVerticalBaselineFunctorTest
)

# -------------  otb::LengthOrientationBaselineFunctor---------------
ADD_TEST(blTuLengthOrientationBaselineFunctorNew ${BASELINE_TESTS}
otbLengthOrientationBaselineFunctorNew
)

ADD_TEST(blTuLengthOrientationBaselineFunctorTest ${BASELINE_TESTS}
otbLengthOrientationBaselineFunctorTest
)

# -------------  otb::ParallelPerpendicularBaselineFunctor---------------
ADD_TEST(blTuParallelPerpendicularBaselineFunctorNew ${BASELINE_TESTS}
otbParallelPerpendicularBaselineFunctorNew
)

ADD_TEST(blTuParallelPerpendicularBaselineFunctorTest ${BASELINE_TESTS}
otbParallelPerpendicularBaselineFunctorTest
)

IF(OTB_DATA_USE_LARGEINPUT)

ENDIF(OTB_DATA_USE_LARGEINPUT)
							       
SET(BASELINE_SRCS
otbBaselineTests.cxx
otbBaselineFunctorBaseNew.cxx
otbHorizontalVerticalBaselineFunctorNew.cxx
otbHorizontalVerticalBaselineFunctorTest.cxx
otbLengthOrientationBaselineFunctorNew.cxx
otbLengthOrientationBaselineFunctorTest.cxx
otbParallelPerpendicularBaselineFunctorNew.cxx
otbParallelPerpendicularBaselineFunctorTest.cxx
)


ADD_EXECUTABLE(otbBaselineTests ${BASELINE_SRCS})
TARGET_LINK_LIBRARIES(otbBaselineTests OTBTesting OTBIO OTBCommon  ITKCommon Baseline)
