IF(WIN32)
  ADD_DEFINITIONS(-DWIN32)
ENDIF(WIN32)


SET(BASELINE ${InSAR_DATA}/Baseline)
SET(INPUTDATA ${InSAR_DATA}/Input)

IF(INSAR_DATA_USE_LARGEINPUT)
  SET(LARGEINPUT ${INSAR_DATA_LARGEINPUT_ROOT} )
ENDIF(INSAR_DATA_USE_LARGEINPUT)

SET(EPSILON 0.01)

SET(PROJECTION_TESTS ${CXX_TEST_PATH}/otbProjectionTests)


#----- otb::MultivariateRationalTransform ------------------------
ADD_TEST(prTuMultivariateRationalTransformNew ${PROJECTION_TESTS}
otbMultivariateRationalTransformNew)

ADD_TEST(prTvMultivariateRationalTransform ${PROJECTION_TESTS}
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

# -------      otb::MultivariateRationalTransformXMLFileReader -----------------------------------
 
ADD_TEST(prTuMultivariateRationalTransformXMLFileReaderNew ${PROJECTION_TESTS} 
	otbMultivariateRationalTransformXMLFileReaderNew)

ADD_TEST(prTvMultivariateRationalTransformXMLFileReader ${PROJECTION_TESTS} 
	otbMultivariateRationalTransformXMLFileReader
        ${BASELINE_FILES}/otbMultivariateRationalTransform_test.xml
        )

# -------      otb::MultivariateRationalTransformXMLFileWriter -----------------------------------
 
ADD_TEST(prTuHpMultivariateRationalTransformXMLFileWriterNew ${PROJECTION_TESTS} 
	otbMultivariateRationalTransformXMLFileWriterNew)

ADD_TEST(prTvHpMultivariateRationalTransformXMLFileWriter ${PROJECTION_TESTS} 
	otbMultivariateRationalTransformXMLFileWriter
        ${TEMP}/otbTvHpMultivariateRationalTransformXMLFileWriter_output.xml
        )


SET(PROJECTION_SRCS
otbProjectionTests.cxx
otbMultivariateRationalTransform.cxx
otbMultivariateRationalTransformXMLFileReaderNew.cxx
otbMultivariateRationalTransformXMLFileReader.cxx
otbMultivariateRationalTransformXMLFileWriterNew.cxx
otbMultivariateRationalTransformXMLFileWriter.cxx
)


ADD_EXECUTABLE(otbProjectionTests ${BASELINE_SRCS})
TARGET_LINK_LIBRARIES(otbProjectionTests OTBTesting OTBIO OTBCommon  ITKCommon Projection)
