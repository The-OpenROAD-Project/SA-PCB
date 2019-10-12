cmake_minimum_required(VERSION 3.1)

project(PCBPLACER)

set(CXX g++)
set(CC gcc)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pthread -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb")

set(BOOST_ROOT /usr/local/boost-1.62.0/)
set(Boost_INCLUDE_DIR /usr/local/boost-1.62.0/include/)
set(BOOST_LIBRARYDIR /usr/local/boost-1.62.0/lib)
find_package(Boost)

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

############################################################
# include directories
############################################################
include_directories(
  ${KICADPARSER_HOME}/src
  )

############################################################
# Create a library
############################################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PCBPLACER_HOME}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PCBPLACER_HOME}/lib)

add_library( PCBPLACERlib STATIC
  ${PCBPLACER_SRC}
)

set_target_properties( PCBPLACERlib
  PROPERTIES OUTPUT_NAME "PCBPLACER"
)

############################################################
# Create an executable
############################################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PCBPLACER_HOME}/bin)

add_executable( PCBPLACER src/main.cpp )
add_dependencies( PCBPLACER kicadpcbparserlib )

target_link_libraries( PCBPLACER
  PCBPLACERlib
  kicadpcbparserlib
  ${Boost_LIBRARIES}
)

############################################################
# Install
############################################################

install(FILES ${PCBPLACER_HEADER} DESTINATION include)
install(TARGETS PCBPLACER DESTINATION bin)
install(TARGETS PCBPLACERlib DESTINATION lib)
