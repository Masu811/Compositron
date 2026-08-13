cmake_minimum_required(VERSION 4.0)

project(compositron)

set(CMAKE_CXX_STANDARD           23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS        OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 REQUIRED)
find_package(Ceres REQUIRED)

add_library(${PROJECT_NAME} SHARED
    src/core/utils.cpp
    src/core/fitting.cpp
    src/dbs/fitting.cpp
    src/dbs/dbspectrum.cpp
    src/cdbs/cdbspectrum.cpp
    src/core/measurement.cpp
    src/importers/png_importer.cpp
    src/importers/SLOPE_n42_importer.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${PROJECT_NAME} PUBLIC
    Eigen3::Eigen
    Ceres::ceres
)

add_executable(exampleMeasurement examples/measurement.cpp)

target_link_libraries(exampleMeasurement PRIVATE ${PROJECT_NAME})
