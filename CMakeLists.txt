cmake_minimum_required(VERSION 3.22)

project(MovingCones)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(${CMAKE_SOURCE_DIR}/include)
set(SRC "${CMAKE_SOURCE_DIR}/src/Cholmod_Solver.cpp"
              "${CMAKE_SOURCE_DIR}/src/ConeOptimization.cpp"
              "${CMAKE_SOURCE_DIR}/src/Eigen_Solver.cpp"
              "${CMAKE_SOURCE_DIR}/src/MeshDefinition.cpp"
              "${CMAKE_SOURCE_DIR}/src/MeshLaplace.cpp"
              "${CMAKE_SOURCE_DIR}/src/MeshTopo.cpp")
add_library(${PROJECT_NAME}_lib ${SRC})
add_executable(${PROJECT_NAME} "${CMAKE_SOURCE_DIR}/src/main.cpp")

# Find Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(OpenMesh)
include_directories(${OPENMESH_INCLUDE_DIRS})

# Find CHOLMOD
find_package(CHOLMOD REQUIRED)
if (CHOLMOD_FOUND)
    message(STATUS "CHOLMOD found")
    include_directories(${CHOLMOD_INCLUDE_DIR})
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen OpenMeshCore OpenMeshTools ${CHOLMOD_LIBRARIES} ${PROJECT_NAME}_lib)