# Compile and build the OTB library.
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
#        cmake_policy(SET CMP0015 NEW)
endif(COMMAND cmake_policy)
     
     
PROJECT(ComplexInterpolation)

#-----------------------------------------------------------------------------
# OTB option to access large file (> 4Go).
# ADD_DEFINITIONS(
# 	-D_LARGEFILE_SOURCE
# 	-D_FILE_OFFSET_BITS=64
# 	-D_LARGEFILE64_SOURCE
# 	-D_LARGE_FILES
# 	)

#-------------------------------
# Librairie OTB
#-------------------------------
FIND_PACKAGE(OTB)
IF(OTB_FOUND)
  INCLUDE(${OTB_USE_FILE})
ELSE(OTB_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build ComplexInterpolation without OTB.  Please set OTB_DIR.")
ENDIF(OTB_FOUND)

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_BINARY_DIR} 
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/Code
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/Interpolation
  ${CMAKE_CURRENT_BINARY_DIR}/Code
)


# Find the data directory
# SET(OTB_SR_DIR /home2/arnaud/PROJETS/SimuCapteur/data/BD-CNES)
# SET(OTB_SATRSR_DIR /home2/arnaud/PROJETS/SimuCapteur/hg/SimuPerformances/Data/Rsr)
# FIND_PATH(OTB_SPECTRAL_RESPONSE otbSR.readme ${OTB_SR_DIR} $ENV{OTB_SPECTRAL_RESPONSE})
# FIND_PATH(OTB_SAT_RSR otbSATRSR.readme ${OTB_SATRSR_DIR} $ENV{OTB_SAT_RSR})

SET(ComplexInterpolation_DATA ${ComplexInterpolation_SOURCE_DIR}/Data)
SET(BASELINE_FILES ${ComplexInterpolation_DATA}/Baseline)
SET(INPUTDATA ${ComplexInterpolation_DATA}/Input)


# Output directories.
SET (LIBRARY_OUTPUT_PATH ${ComplexInterpolation_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${ComplexInterpolation_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)
SET(ComplexInterpolation_LIBRARY_PATH "${LIBRARY_OUTPUT_PATH}")
SET(ComplexInterpolation_EXECUTABLE_PATH "${EXECUTABLE_OUTPUT_PATH}")
SET (ComplexInterpolationTesting_BINARY_DIR ${ComplexInterpolation_BINARY_DIR}/Testing)
SET(TEMP ${ComplexInterpolationTesting_BINARY_DIR}/Temporary)

INCLUDE_DIRECTORIES( ${ComplexInterpolation_LIBRARY_PATH} )

# Put "." in the rpath of all ComplexInterpolation exe.
LINK_DIRECTORIES(.)

# IF(WIN32)
#    IF(BUILD_SHARED_LIBS)
#      SET(ComplexInterpolation_COMMON_BUILD_TYPE "SHARED")
#    ELSE(BUILD_SHARED_LIBS)
#      SET(ComplexInterpolation_COMMON_BUILD_TYPE "STATIC")
#    ENDIF(BUILD_SHARED_LIBS)
#    SET(BUILD_SHARED_LIBS OFF)
# ENDIF(WIN32)

SET(ComplexInterpolation_INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib/ComplexInterpolation")
SET(ComplexInterpolation_INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin")
SET(ComplexInterpolation_INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include/ComplexInterpolation")

ENABLE_TESTING()

SUBDIRS(Code Testing)




