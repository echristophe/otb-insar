
# Sources of non-templated classes.
FILE(GLOB Baseline_SRCS "*.cxx" )

ADD_LIBRARY(Baseline ${Baseline_SRCS})
TARGET_LINK_LIBRARIES (Baseline OTBCommon OTBOssimAdapters)

INSTALL(TARGETS Baseline
RUNTIME DESTINATION ${InSAR_INSTALL_BIN_DIR} COMPONENT RuntimeLibraries
LIBRARY DESTINATION ${InSAR_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
ARCHIVE DESTINATION ${InSAR_INSTALL_LIB_DIR} COMPONENT Development)

FILE(GLOB __files1 "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB __files2 "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")

INSTALL(FILES ${__files1} ${__files2}
    DESTINATION ${InSAR_INSTALL_INCLUDE_DIR}
    COMPONENT Development)