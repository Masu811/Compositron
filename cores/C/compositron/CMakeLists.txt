cmake_minimum_required(VERSION 4.0)

project(compositron C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_options(
    -fsanitize=address
    -fno-omit-frame-pointer
    -g
)

add_link_options(
    -fsanitize=address
)

add_library(${PROJECT_NAME} SHARED
    src/cdbs/cdbspectrum.c
    src/core/measurement.c
    src/core/utils.c
    src/dbs/dbspectrum.c
    src/importers/png_importer.c
    src/importers/SLOPE_n42_importer.c
)

target_include_directories(${PROJECT_NAME} PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

add_executable(exampleMeasurement examples/measurement.c)

find_package(LibXml2 REQUIRED)

target_link_libraries(exampleMeasurement PRIVATE ${PROJECT_NAME} LibXml2::LibXml2)
