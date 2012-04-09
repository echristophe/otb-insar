# Only if OTB is generate with Visu GUI option
IF(OTB_USE_VISU_GUI)

IF(WIN32)
  ADD_DEFINITIONS(-DWIN32)
ENDIF(WIN32)

SET(MODULE_SRCS
otbDynamicFringeViewerModule.cxx
)


FLTK_WRAP_UI(otbDynamicFringeViewerModule otbDynamicFringeViewerModuleGUI.fl)

# Suppress warning Fl_RGB_Image deferencing pointeur
IF(CMAKE_COMPILER_IS_GNUCXX)
      SET_SOURCE_FILES_PROPERTIES(
            ${MODULE_SRCS} 
            ${otbDynamicFringeViewerModule_FLTK_UI_SRCS}
            PROPERTIES COMPILE_FLAGS -fno-strict-aliasing)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

INCLUDE_DIRECTORIES(${InSAR_BINARY_DIR}/Monteverdi/Modules/DynamicFringeViewer)
INCLUDE_DIRECTORIES(${InSAR_SOURCE_DIR}/Monteverdi/Modules/DynamicFringeViewer)


ADD_LIBRARY(otbDynamicFringeViewerModule otbDynamicFringeViewerModule.cxx 
                                         ${otbDynamicFringeViewerModule_FLTK_UI_SRCS})
TARGET_LINK_LIBRARIES(otbDynamicFringeViewerModule OTBIO OTBVisualization OTBLearning OTBMonteverdiCommon OTBGui OTBVisualization flu OTBInSARMonteverdiCommon)

ENDIF(OTB_USE_VISU_GUI)
