IF(WIN32)
  ADD_DEFINITIONS(-DWIN32)
ENDIF(WIN32)


SET(BASELINE ${InSAR_DATA}/Baseline)
SET(INPUTDATA ${InSAR_DATA}/Input)

SET(EPSILON 0.01)

SET(INTERPOLATION_TESTS ${CXX_TEST_PATH}/otbInterpolationTests)

# -------------  otb::ComplexInterpolateImageFunction ----------------------------
ADD_TEST(coTuComplexInterpolateImageFunctionNew ${INTERPOLATION_TESTS}
otbComplexInterpolateImageFunctionNew
)

ADD_TEST(coTvComplexBlackmanWindowInterpolateImageFunction ${INTERPOLATION_TESTS}
#     --compare-ascii ${NOTOL}
#                          ${BASELINE_FILES}/coRealBlackmanWindowInterpolateImageFunctionOutput.txt
#                     ${TEMP}/coRealBlackmanWindowInterpolateImageFunctionOutput.txt
        otbComplexInterpolateImageFunction
                ${INPUTDATA}/prison_toulouse.tif
                ${TEMP}/coRealBlackmanWindowInterpolateImageFunctionOutput.txt
				0 # interpolation ID "Blackman"
				3 # radius
                0.01 # normalizeZeroFrequency
                0.5 0.5
                127.33 44.9
                259.67 21.43
                12.13 61.79
                89.5 11
                128 128
                127.255 128.73
                -1 -1
                )

ADD_TEST(coTvComplexHammingWindowInterpolateImageFunction ${INTERPOLATION_TESTS}
#     --compare-ascii ${NOTOL}
#                          ${BASELINE_FILES}/coComplexHammingWindowInterpolateImageFunctionOutput.txt
#                     ${TEMP}/coComplexHammingWindowInterpolateImageFunctionOutput.txt
        otbComplexInterpolateImageFunction
                ${INPUTDATA}/prison_toulouse.tif
                ${TEMP}/coComplexHammingWindowInterpolateImageFunctionOutput.txt
				1 # interpolation ID "Hamming"
				3 # radius
                0.01 # normalizeZeroFrequency
                0.5 0.5 
				127.33 44.9
                259.67 21.43
                12.13 61.79
                89.5 11
                128 128
                127.255 128.73
                -1 -1
                )

ADD_TEST(coTvComplexHammingWindowResampleImageFilter ${INTERPOLATION_TESTS}
#     --compare-ascii ${NOTOL}
#                          ${BASELINE_FILES}/coComplexHammingWindowResampleImageFilterOutput.tif
#                     ${TEMP}/coComplexHammingWindowResampleImageFilterOutput.tif
        otbComplexResampleImageFilter
                ${INPUTDATA}/prison_toulouse.tif
                ${TEMP}/coComplexHammingWindowResampleImageFilterOutput.tif
				0 # interpolation ID "Hamming"
				3 # radius
                0.00 # normalizeZeroFrequency
                )

ADD_TEST(coTvRealHammingWindowResampleImageFilter ${INTERPOLATION_TESTS}
#     --compare-ascii ${NOTOL}
#                          ${BASELINE_FILES}/coComplexHammingWindowResampleImageFilterOutput.tif
#                     ${TEMP}/coComplexHammingWindowResampleImageFilterOutput.tif
        otbRealResampleImageFilter
                ${INPUTDATA}/prison_toulouse.tif
                ${TEMP}/coRealHammingWindowResampleImageFilterOutput.tif
				0 # interpolation ID "Hamming"
				3 # radius
                )	

ADD_TEST(TvGridIntersectionPointSetSource ${INTERPOLATION_TESTS}
  otbGridIntersectionPointSetSourceTest)

							       
SET(INTERPOLATION_SRCS
otbInterpolationTests.cxx
otbComplexInterpolateImageFunctionNew.cxx
otbComplexInterpolateImageFunction.cxx
otbComplexResampleImageFilter.cxx
otbRealResampleImageFilter.cxx
otbGridIntersectionPointSetSourceTest.cxx
)


ADD_EXECUTABLE(otbInterpolationTests ${INTERPOLATION_SRCS})
TARGET_LINK_LIBRARIES(otbInterpolationTests OTBTesting OTBIO OTBCommon  ITKCommon Interpolation)
