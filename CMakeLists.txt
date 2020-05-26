cmake_minimum_required(VERSION 3.1)

project(PCBPLACER)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(THREADS_PREFER_PTHREAD_FLAG ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pthread -Ofast -funsafe-math-optimizations -fpermissive -fPIC")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb")

find_package(Boost REQUIRED)
find_library(SWIG_EXECUTABLE swig)
if(NOT SWIG_EXECUTABLE)
set(SWIG_EXECUTABLE "swig")
endif()
#find_package( PythonInterp 2.7 REQUIRED )
find_package( PythonInterp 3 REQUIRED )
#find_package( PythonLibs 2.7 REQUIRED )
find_package( PythonLibs 3 REQUIRED )

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARYDIR})

set(PCBPLACER_HOME ${PROJECT_SOURCE_DIR} )
set(KICADPARSER_HOME
  ${PROJECT_SOURCE_DIR}/module/KicadParser
)

add_subdirectory(${KICADPARSER_HOME})


set (PCBPLACER_SRC  
  src/GridBasedPlacer.cpp
  src/globalParam.cpp
  src/readFiles.cpp
  src/PcbPlacer_wrapper.cpp
  src/linreg.cpp
  src/Module.cpp
  src/Hierarchy.cpp
  src/Node.cpp
  #src/HPlacerUtils.cpp
  )

set (PCBPLACER_HEADER
  src/GridBasedPlacer.hpp
  src/globalParam.h
  src/util.h
  src/Node.hpp
  src/Hierarchy.hpp
  src/Module.hpp
  src/readFiles.hpp
  src/Pin.hpp
  src/linreg.h
  src/Logger.hpp
  #src/HPlacerUtils.hpp
  )

set (PLACER_SWIG_FILES
  src/PcbPlacer.i
)

############################################################
# include directories
############################################################
include_directories(
  ${KICADPARSER_HOME}/src
  ${PYTHON_INCLUDE_DIRS}
  )

############################################################
# Create a library
############################################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PCBPLACER_HOME}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PCBPLACER_HOME}/lib)

add_library( PCBPLACERlib SHARED
  ${PCBPLACER_SRC}
)

set_target_properties( PCBPLACERlib
  PROPERTIES OUTPUT_NAME "PCBPLACER"
)

############################################################
# Create an executable
############################################################
add_custom_command( OUTPUT ${PCBPLACER_HOME}/src/PcbPlacer_wrapper.cpp ${PCBPLACER_LIB_FILES}
  COMMAND ${SWIG_EXECUTABLE} -python -py3 -c++ -o ${PCBPLACER_HOME}/src/PcbPlacer_wrapper.cpp ${PCBPLACER_HOME}/src/PcbPlacer.i   
  WORKING_DIRECTORY ${PCBPLACER_HOME}
  DEPENDS ${PLACER_SWIG_FILES}
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PCBPLACER_HOME}/bin)

add_executable( PCBPLACER src/main.cpp )
add_dependencies( PCBPLACER kicadpcbparserlib )

target_link_libraries( PCBPLACER
  PCBPLACERlib
  kicadpcbparserlib
  ${Boost_LIBRARIES}
  ${PYTHON_LIBRARIES}
)

############################################################
# Install
############################################################

install(FILES ${PCBPLACER_HEADER} DESTINATION include)
install(TARGETS PCBPLACER DESTINATION bin)
install(TARGETS PCBPLACERlib DESTINATION lib)
