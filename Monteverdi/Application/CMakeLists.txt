#Create Monteverdi application
ADD_EXECUTABLE(InSARMonteverdi InSARMonteverdi.cxx)
TARGET_LINK_LIBRARIES(InSARMonteverdi otbMonteverdi 
                      otbDynamicFringeViewerModule 
                      ${Monteverdi_MODULES_LIBRARIES})
