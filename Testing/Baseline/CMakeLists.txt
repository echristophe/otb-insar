IF(WIN32)
  ADD_DEFINITIONS(-DWIN32)
ENDIF(WIN32)


SET(BASELINE ${InSAR_DATA}/Baseline)
SET(INPUTDATA ${InSAR_DATA}/Input)

IF(INSAR_DATA_USE_LARGEINPUT)
  SET(LARGEINPUT ${INSAR_DATA_LARGEINPUT_ROOT} )
ENDIF(INSAR_DATA_USE_LARGEINPUT)

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

# -------------  otb::PlatformPositionToBaselineCalculator---------------
ADD_TEST(blTuPlatformPositionToBaselineCalculatorNew ${BASELINE_TESTS}
otbPlatformPositionToBaselineCalculatorNew
)

IF(INSAR_DATA_USE_LARGEINPUT)
ADD_TEST(blTvPlatformPositionToBaselineCalculatorTest ${BASELINE_TESTS}
otbPlatformPositionToBaselineCalculatorTest
	${LARGEINPUT}/TERRASARX/2008-03-10_GrandCanyon_SSC/TSX1_SAR__SSC______SM_S_SRA_20080310T133220_20080310T133228/TSX1_SAR__SSC______SM_S_SRA_20080310T133220_20080310T133228.xml
	${LARGEINPUT}/TERRASARX/2008-03-21_GrandCanyon_SSC/TSX1_SAR__SSC______SM_S_SRA_20080321T133220_20080321T133228/TSX1_SAR__SSC______SM_S_SRA_20080321T133220_20080321T133228.xml
)
ENDIF(INSAR_DATA_USE_LARGEINPUT)
							       
# -------------  otb::BaselineCalculator---------------
ADD_TEST(blTuBaselineCalculatorNew ${BASELINE_TESTS}
otbPlatformPositionToBaselineCalculatorNew
)

IF(INSAR_DATA_USE_LARGEINPUT)
ADD_TEST(blTvBaselineCalculatorTest ${BASELINE_TESTS}
otbBaselineCalculatorTest
	${LARGEINPUT}/TERRASARX/2008-03-10_GrandCanyon_SSC/TSX1_SAR__SSC______SM_S_SRA_20080310T133220_20080310T133228/TSX1_SAR__SSC______SM_S_SRA_20080310T133220_20080310T133228.xml
	${LARGEINPUT}/TERRASARX/2008-03-21_GrandCanyon_SSC/TSX1_SAR__SSC______SM_S_SRA_20080321T133220_20080321T133228/TSX1_SAR__SSC______SM_S_SRA_20080321T133220_20080321T133228.xml
)
ENDIF(INSAR_DATA_USE_LARGEINPUT)

#----- otb::MultivariateRationalTransform ------------------------
ADD_TEST(blTuMultivariateRationalTransformNew ${BASELINE_TESTS}
otbMultivariateRationalTransformNew)

ADD_TEST(prTvMultivariateRationalTransform ${BASELINE_TESTS}
#  --compare-ascii ${NOTOL}
#  ${BASELINE_FILES}/otbMultivariateRationalTransformOutput.txt
#  ${TEMP}/otbMultivariateRationalTransformOutput.txt
otbMultivariateRationalTransformTest
${TEMP}/otbMultivariateRationalTransformOutput.txt
0 0
1 1
10 10
-10 -10
)

SET(BASELINE_SRCS
otbBaselineTests.cxx
otbBaselineFunctorBaseNew.cxx
otbHorizontalVerticalBaselineFunctorNew.cxx
otbHorizontalVerticalBaselineFunctorTest.cxx
otbLengthOrientationBaselineFunctorNew.cxx
otbLengthOrientationBaselineFunctorTest.cxx
otbParallelPerpendicularBaselineFunctorNew.cxx
otbParallelPerpendicularBaselineFunctorTest.cxx
otbPlatformPositionToBaselineCalculatorNew.cxx
otbPlatformPositionToBaselineCalculatorTest.cxx
otbBaselineCalculatorNew.cxx
otbBaselineCalculatorTest.cxx
otbMultivariateRationalTransform.cxx
)


ADD_EXECUTABLE(otbBaselineTests ${BASELINE_SRCS})
TARGET_LINK_LIBRARIES(otbBaselineTests OTBTesting OTBIO OTBCommon  ITKCommon Baseline)
