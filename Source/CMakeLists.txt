cmake_minimum_required(VERSION 3.9)
project(PolyhedralSplinesIGA)

# Required by OpenMesh
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Force make to stop after first error
if(UNIX)
    add_definitions(
        -Wfatal-errors # stop after first error
        -w
    )
endif()

add_compile_definitions(_USE_MATH_DEFINES)

################################################################################
# Clone and build install dependencies
################################################################################
include(ExternalProject)

# OpenMesh
set(OPENMESH_ROOT ${CMAKE_SOURCE_DIR}/External/OpenMesh)

ExternalProject_Add(
    OpenMesh
    PREFIX ${OPENMESH_ROOT}
    GIT_REPOSITORY https://www.graphics.rwth-aachen.de:9000/OpenMesh/OpenMesh.git
    GIT_TAG db8c0608830ae7324b2b2e1d91e4db565b2d3eae
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
               -DBUILD_APPS=FALSE
               -DCMAKE_BUILD_TYPE=RELEASE
)

# Eigen
set(Eigen_ROOT ${CMAKE_SOURCE_DIR}/External/Eigen)

ExternalProject_Add(
    Eigen
    PREFIX ${Eigen_ROOT}
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
               -DBUILD_APPS=FALSE
               -DCMAKE_BUILD_TYPE=RELEASE
)

################################################################################

link_directories(${OPENMESH_ROOT}/lib)

# Find source files
file(GLOB_RECURSE SOURCES src/*.cpp)
add_executable(PolyhedralSplinesIGA src/main.cpp ${SOURCES})

### Link OpenMesh ###
include_directories(SYSTEM ${OPENMESH_ROOT}/include)
target_link_libraries(PolyhedralSplinesIGA OpenMeshCore OpenMeshTools)

# OpenMesh need to be built before main program
add_dependencies(PolyhedralSplinesIGA OpenMesh)

### Link Eigen ###
add_dependencies(PolyhedralSplinesIGA Eigen) # Make sure Eigen is downloaded before building PolyhedralSplinesIGA
include_directories(SYSTEM ${Eigen_ROOT}/src/Eigen) # Adjust path if Eigen's headers are in a different subdirectory


# Copy tables to build folder
set(TABLE_SRC_PATH ${CMAKE_SOURCE_DIR}/src/Patch/Table)
set(SUPPORTS_SRC_PATH ${CMAKE_SOURCE_DIR}/src/Patch/SupportInfo)
set(TESTFILE_SRC_PATH ${CMAKE_SOURCE_DIR}/testfile)

if(UNIX)
    set(TABLE_TRG_PATH ${CMAKE_SOURCE_DIR}/../build/Table)
    set(SUPPORTS_TRG_PATH ${CMAKE_SOURCE_DIR}/../build/SupportInfo)
    set(TESTFILE_TRG_PATH ${CMAKE_SOURCE_DIR}/../build/testfile)
elseif(WIN32)
    set(TABLE_TRG_PATH ${CMAKE_SOURCE_DIR}/../build/Release/Table)
    set(SUPPORTS_TRG_PATH ${CMAKE_SOURCE_DIR}/../build/Release/SupportInfo)
    set(TESTFILE_TRG_PATH ${CMAKE_SOURCE_DIR}/../build/Release/testfile)
endif()

add_custom_command(
    TARGET PolyhedralSplinesIGA
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${TABLE_SRC_PATH} ${TABLE_TRG_PATH}
)

add_custom_command(
        TARGET PolyhedralSplinesIGA
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${SUPPORTS_SRC_PATH} ${SUPPORTS_TRG_PATH}
)

add_custom_command(
        TARGET PolyhedralSplinesIGA
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${TESTFILE_SRC_PATH} ${TESTFILE_TRG_PATH}
)
