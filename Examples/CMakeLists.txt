cmake_minimum_required(VERSION 3.0.0)
project(Examples LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Examples/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Examples/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Examples/lib)
set(HEADER_OUTPUT_PATH ${CMAKE_BINARY_DIR}/include/Examples)

include_directories(${PROJECT_SOURCE_DIR}/include)

#-----Libraries
add_library(ExampleParameters SHARED ${PROJECT_SOURCE_DIR}/src/ExampleParameters.cpp)
add_library(ExampleExperiment SHARED ${PROJECT_SOURCE_DIR}/src/ExampleExperiment.cpp)
add_library(ExampleModule     SHARED ${PROJECT_SOURCE_DIR}/src/ExampleModule.cpp)
target_link_libraries(ExampleParameters FitParameters Parameter)
target_link_libraries(ExampleExperiment ${ROOT_LIBRARIES} FitParameters Parameter)
target_link_libraries(ExampleModule ExampleParameters ToyModule ExampleExperiment ${ROOT_LIBRARIES})
