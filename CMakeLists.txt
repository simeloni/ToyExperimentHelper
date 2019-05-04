cmake_minimum_required(VERSION 3.0.0)
project(Main VERSION 0.1.0 LANGUAGES CXX)

#-------------------------------------------------------------------------------------------
# Set the paths where to put the libraries
#-------------------------------------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(HEADER_OUTPUT_PATH ${CMAKE_BINARY_DIR}/include)

include(CTest)
enable_testing()

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)

set(main_SRC main.cpp)

#-------Executables-----------------
add_executable(Main ${main_SRC})
target_link_libraries(Main ${ROOT_LIBRARIES})


#------Libraries--------------------
## To actually compile the library
add_library(ToyExperiment SHARED ${PROJECT_SOURCE_DIR}/src/ToyExperiment.cpp)
target_link_libraries(Main ToyExperiment)

# TODO: Understand how you could compile this as a library and then import it 
#       in another program


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
