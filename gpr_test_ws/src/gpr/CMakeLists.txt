cmake_minimum_required(VERSION 2.8.3)
project(gpr)

set(CMAKE_VERBOSE_MAKEFILE "false")

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  geometry_msgs
  sensor_msgs
  mavros
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(Eigen3 REQUIRED)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN_INCLUDE_DIRS}
  include/${PROJECT_NAME}
  include
)

catkin_package(
    )

add_executable(gpr 
  src/gpr_node.cpp
  src/gpr.cpp
  src/input.cpp)
target_link_libraries(gpr
  ${catkin_LIBRARIES}
)
