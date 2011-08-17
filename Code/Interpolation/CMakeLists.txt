
# Sources of non-templated classes.
FILE(GLOB Interpolation_SRCS "*.cxx" )

ADD_LIBRARY(Interpolation ${Interpolation_SRCS})
TARGET_LINK_LIBRARIES (Interpolation OTBCommon)

INSTALL(TARGETS Interpolation
RUNTIME DESTINATION ${ComplexInterpolation_INSTALL_BIN_DIR} COMPONENT RuntimeLibraries
LIBRARY DESTINATION ${ComplexInterpolation_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
ARCHIVE DESTINATION ${ComplexInterpolation_INSTALL_LIB_DIR} COMPONENT Development)

FILE(GLOB __files1 "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB __files2 "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")

INSTALL(FILES ${__files1} ${__files2}
    DESTINATION ${ComplexInterpolation_INSTALL_INCLUDE_DIR}
    COMPONENT Development)