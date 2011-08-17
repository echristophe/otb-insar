IF(WIN32)
  ADD_DEFINITIONS(-DWIN32)
ENDIF(WIN32)


SET(BASELINE ${ComplexInterpolation_DATA}/Baseline)
SET(INPUTDATA ${ComplexInterpolation_DATA}/Input)

SET(EPSILON 0.01)

SET(INTERPOLATION_TESTS ${CXX_TEST_PATH}/otbInterpolationTests)

# -------------  otb::ComplexInterpolateImageFunction ----------------------------
ADD_TEST(coTuComplexInterpolateImageFunctionNew ${INTERPOLATION_TESTS}
otbComplexInterpolateImageFunctionNew
)

ADD_TEST(coTvComplexInterpolateImageFunction ${INTERPOLATION_TESTS}
#     --compare-ascii ${NOTOL}
#                          ${BASELINE_FILES}/coComplexInterpolateImageFunctionOutput.txt
#                     ${TEMP}/coComplexInterpolateImageFunctionOutput.txt
        otbComplexInterpolateImageFunction
                ${INPUTDATA}/prison_toulouse.tif
                ${TEMP}/coComplexInterpolateImageFunctionOutput.txt
                3 # radius
                0.10 # normalizeZeroFrequency
                0.5 0.5
                127.33 44.9
                259.67 21.43
                12.13 61.79
                89.5 11
                128 128
                127.255 128.73
                -1 -1
                )


       
SET(INTERPOLATION_SRCS
otbInterpolationTests.cxx
otbComplexInterpolateImageFunctionNew.cxx
otbComplexInterpolateImageFunction.cxx
)


ADD_EXECUTABLE(otbInterpolationTests ${INTERPOLATION_SRCS})
TARGET_LINK_LIBRARIES(otbInterpolationTests OTBTesting OTBIO OTBCommon  ITKCommon Interpolation)
