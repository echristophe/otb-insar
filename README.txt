Original situation for the ESA-SOCIS project

* a basic standalone example for Insar is in
Examples/InSar/BasicInterferogramComputation.cxx
This does only integer pixel registration (no interpolation). See Emmanuel Christophe for more details.
Remember to use a version of OTB with Concept Checking in ITK turned off. To do this, you have to set (under cmake Advanced Options) ITK_USE_CONCEPT_CHECKING to OFF when compiling OTB.

* an implementation for complex interpolation. See Patrick Imbo for details.
