# Compile and build the InSAR project.
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
        cmake_policy(SET CMP0015 NEW)
endif(COMMAND cmake_policy)
     
PROJECT(InSARDocuments)

# Include latex cmake mode
INCLUDE(CMake/UseLATEX.cmake)

SUBDIRS(CR)
SUBDIRS(SOCISReport)

