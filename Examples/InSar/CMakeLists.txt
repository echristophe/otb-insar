PROJECT(Sar)

cmake_minimum_required(VERSION 2.4)

FIND_PACKAGE(OTB)
IF(OTB_FOUND)
  INCLUDE(${OTB_USE_FILE})
ELSE(OTB_FOUND)
  MESSAGE(FATAL_ERROR
          "OTB not found. Please set OTB_DIR")
ENDIF(OTB_FOUND)



ADD_EXECUTABLE(BasicInterferogramComputation BasicInterferogramComputation.cxx)
TARGET_LINK_LIBRARIES(BasicInterferogramComputation OTBCommon OTBIO)

ADD_EXECUTABLE(SARCoregistrationComputation SARCoregistrationComputation.cxx)
TARGET_LINK_LIBRARIES(SARCoregistrationComputation OTBCommon OTBIO)

ADD_EXECUTABLE(SensorModelRegistration SensorModelRegistration.cxx)
TARGET_LINK_LIBRARIES(SensorModelRegistration OTBCommon OTBIO Interpolation)

ADD_EXECUTABLE(ComplexCoarseRegistrationComputation ComplexCoarseRegistrationComputation.cxx)
TARGET_LINK_LIBRARIES(ComplexCoarseRegistrationComputation OTBCommon OTBIO Interpolation)
