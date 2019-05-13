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

find_package(ROOT REQUIRED COMPONENTS MathMore RooFit)
include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/Examples/include)

set(main_SRC main.cpp)

#-------Executables-----------------
add_executable(Main ${main_SRC})
#----Subdirectories
add_subdirectory(Examples)

#------Libraries--------------------
## To actually compile the library
add_library(ToyModule SHARED ${PROJECT_SOURCE_DIR}/src/ToyModule.cpp)
add_library(ToyExperiment SHARED ${PROJECT_SOURCE_DIR}/src/ToyExperiment.cpp)
add_library(Parameter SHARED ${PROJECT_SOURCE_DIR}/src/Parameter.cpp)
add_library(FitParameters SHARED ${PROJECT_SOURCE_DIR}/src/FitParameters.cpp)
add_library(funcs SHARED ${PROJECT_SOURCE_DIR}/src/funcs.cpp)
target_link_libraries(funcs ${ROOT_LIBRARIES})
target_link_libraries(ToyModule ${ROOT_LIBRARIES})
target_link_libraries(Parameter funcs)
target_link_libraries(FitParameters Parameter)
target_link_libraries(Main ToyExperiment ToyModule ExampleParameters ExampleExperiment ExampleModule)

# TODO: Understand how you could compile this as a library and then import it 
#       in another program

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
