cmake_minimum_required(VERSION 3.1)

project(PCBPLACER)

#set(CXX g++)
#set(CC gcc)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(THREADS_PREFER_PTHREAD_FLAG ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pthread -O3 -fpermissive -Wl,--no-as-needed -ldl -lutil")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb")

#set(BOOST_ROOT /usr/local/boost-1.62.0/)
#set(Boost_INCLUDE_DIR /usr/local/boost-1.62.0/include/)
#set(BOOST_LIBRARYDIR /usr/local/boost-1.62.0/lib)
find_package(Boost REQUIRED)
#find_package(Swig REQUIRED)
set(SWIG_EXECUTABLE "/usr/bin/swig")
#find_package(Python REQUIRED)
#find_package(PythonLibs REQUIRED)
find_package( PythonInterp 2.7 REQUIRED )
find_package( PythonLibs 2.7 REQUIRED )

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARYDIR})

set(PCBPLACER_HOME ${PROJECT_SOURCE_DIR} )
set(KICADPARSER_HOME
  ${PROJECT_SOURCE_DIR}/module/KicadParser
)

add_subdirectory(${KICADPARSER_HOME})


set (PCBPLACER_SRC  
  src/BoardGrid.cpp
  src/GridBasedPlacer.cpp
  src/globalParam.cpp
  src/readFiles.cpp
  src/placer_wrapper.cpp
  )

set (PCBPLACER_HEADER
  src/BoardGrid.h
  src/GridBasedPlacer.h
  src/globalParam.h
  src/util.h
  src/Node.hpp
  src/readFiles.hpp
  src/Pin.hpp
  )

set (PLACER_SWIG_FILES
  src/placer.i
)

#set (PCBPLACER_LIB_FILES
#  bin/CMakeFiles/PCBPLACERlib.dir/src/*.o
#)

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

#add_library( PCBPLACERlib STATIC
#  ${PCBPLACER_SRC}
#)

add_library( PCBPLACERlib SHARED
  ${PCBPLACER_SRC}
)

set_target_properties( PCBPLACERlib
  PROPERTIES OUTPUT_NAME "PCBPLACER"
)

############################################################
# Create an executable
############################################################
add_custom_command( OUTPUT ${PCBPLACER_HOME}/src/placer_wrapper.cpp ${PCBPLACER_LIB_FILES}
  COMMAND ${SWIG_EXECUTABLE} -python -c++ -o ${PCBPLACER_HOME}/src/placer_wrapper.cpp ${PCBPLACER_HOME}/src/placer.i   
  WORKING_DIRECTORY ${PCBPLACER_HOME}
  DEPENDS ${PLACER_SWIG_FILES}
)

#add_custom_command( OUTPUT ${PCBPLACER_HOME}/src/_placer.so
#  COMMAND gcc -shared ${PCBPLACER_HOME}/bin/CMakeFiles/PCBPLACERlib.dir/src/GridBasedPlacer.cpp.o ${PCBPLACER_HOME}/bin/CMakeFiles/PCBPLACERlib.dir/src/BoardGrid.cpp.o  ${PCBPLACER_HOME}/bin/CMakeFiles/PCBPLACERlib.dir/src/globalParam.cpp.o ${PCBPLACER_HOME}/bin/CMakeFiles/PCBPLACERlib.dir/src/readFiles.cpp.o ${PCBPLACER_HOME}/bin/CMakeFiles/PCBPLACERlib.dir/src/placer_wrapper.cpp.o -o ${PCBPLACER_HOME}/src/_placer.so -lstdc++
#  WORKING_DIRECTORY ${PCBPLACER_HOME}
#  DEPENDS ${PCBPLACER_HOME}/src/placer_wrapper.cpp ${PCBPLACER_LIB_FILES}
#)

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
