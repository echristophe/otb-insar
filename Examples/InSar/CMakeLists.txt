
ADD_EXECUTABLE(BasicInterferogramComputation BasicInterferogramComputation.cxx)
TARGET_LINK_LIBRARIES(BasicInterferogramComputation OTBCommon OTBIO)

ADD_EXECUTABLE(BaselineComputation BaselineComputation.cxx)
TARGET_LINK_LIBRARIES(BaselineComputation OTBCommon OTBIO Baseline)

IF(VTK_FOUND)
ADD_EXECUTABLE(BaselinePlot BaselinePlot.cxx)
TARGET_LINK_LIBRARIES(BaselinePlot OTBCommon OTBIO Baseline vtkHybrid vtkCharts)
ENDIF(VTK_FOUND)

ADD_EXECUTABLE(JustTheInterferogramComputation JustTheInterferogramComputation.cxx)
TARGET_LINK_LIBRARIES(JustTheInterferogramComputation OTBCommon OTBIO InterferogramFormation)

ADD_EXECUTABLE(JustFlatEarthRemovalInterferogramComputation JustFlatEarthRemovalInterferogramComputation.cxx)
TARGET_LINK_LIBRARIES(JustFlatEarthRemovalInterferogramComputation OTBCommon OTBIO InterferogramFormation)

#ADD_EXECUTABLE(SARCoregistrationComputation SARCoregistrationComputation.cxx)
#TARGET_LINK_LIBRARIES(SARCoregistrationComputation OTBCommon OTBIO)

ADD_EXECUTABLE(SensorModelRegistration SensorModelRegistration.cxx)
TARGET_LINK_LIBRARIES(SensorModelRegistration OTBCommon OTBIO Interpolation)

ADD_EXECUTABLE(ComplexCoarseRegistrationComputation ComplexCoarseRegistrationComputation.cxx)
TARGET_LINK_LIBRARIES(ComplexCoarseRegistrationComputation OTBCommon OTBIO Interpolation)

ADD_EXECUTABLE(FFTConvolutionTest FFTConvolutionTest.cxx)
TARGET_LINK_LIBRARIES(FFTConvolutionTest OTBCommon OTBIO Interpolation)

ADD_EXECUTABLE(CoherenceComputation CoherenceComputation.cxx)
TARGET_LINK_LIBRARIES(CoherenceComputation OTBCommon OTBIO)

#ADD_EXECUTABLE(ComplexCoarseFineStepTest ComplexCoarseFineStepTest.cxx)
#TARGET_LINK_LIBRARIES(ComplexCoarseFineStepTest OTBCommon OTBIO Interpolation)

# Note: couldn't compile this, had to deactivate (emmanuel)
#ADD_EXECUTABLE(InterferometricChain InterferometricChain.cxx)
#TARGET_LINK_LIBRARIES(InterferometricChain OTBCommon OTBIO Coregistration Interpolation InterferogramFormation)

ADD_EXECUTABLE(FlatEarthRemovalFFT FlatEarthRemovalFFT.cxx)
TARGET_LINK_LIBRARIES(FlatEarthRemovalFFT OTBCommon OTBIO InterferogramFormation)

ADD_EXECUTABLE(InterferogramHSVColorGeneration InterferogramHSVColorGeneration.cxx)
TARGET_LINK_LIBRARIES(InterferogramHSVColorGeneration OTBCommon OTBIO)

