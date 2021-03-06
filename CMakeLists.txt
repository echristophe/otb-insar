# Compile and build the InSAR project.
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
#  CMAKE_POLICY(SET CMP0015 NEW)
ENDIF(COMMAND CMAKE_POLICY)
     
PROJECT(InSAR)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules")


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
          "Cannot build InSAR without OTB.  Please set OTB_DIR.")
ENDIF(OTB_FOUND)

OPTION(INSAR_DATA_USE_LARGEINPUT "Use Large inputs images test." OFF)
MARK_AS_ADVANCED(INSAR_DATA_USE_LARGEINPUT)
IF(INSAR_DATA_USE_LARGEINPUT)
  FIND_PATH(INSAR_DATA_LARGEINPUT_ROOT OTBData.readme $ENV{INSAR_DATA_LARGEINPUT_ROOT})
  MARK_AS_ADVANCED(INSAR_DATA_LARGEINPUT_ROOT)
ENDIF(INSAR_DATA_USE_LARGEINPUT)

#-------------------------------
# External libraries
#-------------------------------
FIND_PACKAGE(VTK QUIET)
IF(VTK_FOUND)
	MESSAGE(STATUS "VTK Found")
	INCLUDE(${VTK_USE_FILE})
ENDIF(VTK_FOUND)

#-------------------------------
# Librairie FFTW
#-------------------------------
find_package(fftw3 REQUIRED)


INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/Code
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/Interpolation
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/InterferogramFormation
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/Coregistration
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/Baseline
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/BasicFilter
  ${CMAKE_CURRENT_BINARY_DIR}/Code
  ${FFTW3_INCLUDE_DIRS}
)

# Generate Monteverdi modules
OPTION(GENERATE_MONTEVERDI_APPLICATION "Generate the Monteverdi application." OFF)

IF(GENERATE_MONTEVERDI_APPLICATION)
  #Find Monteverdi and load its settings.
  FIND_PACKAGE(Monteverdi)
  IF(Monteverdi_FOUND)
    MESSAGE("Monteverdi use file: ${Monteverdi_USE_FILE}")
    INCLUDE(${Monteverdi_USE_FILE})
    INCLUDE_DIRECTORIES(${InSAR_SOURCE_DIR}/Monteverdi/Modules/DynamicFringeViewer/
	                    ${InSAR_SOURCE_DIR}/Monteverdi/Common/    
                        ${InSAR_BINARY_DIR}/Monteverdi/Modules/DynamicFringeViewer/
                        )
    SUBDIRS(Monteverdi)

  ELSE(Monteverdi_FOUND)
    MESSAGE("Monteverdi not found.  Please set Monteverdi_DIR")
  ENDIF(Monteverdi_FOUND)
ENDIF(GENERATE_MONTEVERDI_APPLICATION)


SET(InSAR_DATA ${InSAR_SOURCE_DIR}/Data)
SET(BASELINE_FILES ${InSAR_DATA}/Baseline)
SET(INPUTDATA ${InSAR_DATA}/Input)


# Output directories.
SET (LIBRARY_OUTPUT_PATH ${InSAR_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${InSAR_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)
SET(InSAR_LIBRARY_PATH "${LIBRARY_OUTPUT_PATH}")
SET(InSAR_EXECUTABLE_PATH "${EXECUTABLE_OUTPUT_PATH}")
SET (InSARTesting_BINARY_DIR ${InSAR_BINARY_DIR}/Testing)
SET(TEMP ${InSARTesting_BINARY_DIR}/Temporary)

INCLUDE_DIRECTORIES( ${InSAR_LIBRARY_PATH} )

# Put "." in the rpath of all InSAR exe.
LINK_DIRECTORIES(.)

# IF(WIN32)
#    IF(BUILD_SHARED_LIBS)
#      SET(InSAR_COMMON_BUILD_TYPE "SHARED")
#    ELSE(BUILD_SHARED_LIBS)
#      SET(InSAR_COMMON_BUILD_TYPE "STATIC")
#    ENDIF(BUILD_SHARED_LIBS)
#    SET(BUILD_SHARED_LIBS OFF)
# ENDIF(WIN32)

SET(InSAR_INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib/InSAR")
SET(InSAR_INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin")
SET(InSAR_INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include/InSAR")


ADD_SUBDIRECTORY(Code)

OPTION(BUILD_TESTING  "Build testing." ON)
IF (BUILD_TESTING)
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(Testing)
ENDIF (BUILD_TESTING)

OPTION(BUILD_EXAMPLES "Build the Examples directory." ON)
IF (BUILD_EXAMPLES)
  ADD_SUBDIRECTORY(Examples)
ENDIF (BUILD_EXAMPLES)

OPTION(GENERATE_DOCUMENTATION "Generate the LateX documentation." ON)
IF(GENERATE_DOCUMENTATION)
  ADD_SUBDIRECTORY(Documents)
ENDIF(GENERATE_DOCUMENTATION)

