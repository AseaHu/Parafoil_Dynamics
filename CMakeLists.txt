cmake_minimum_required(VERSION 3.1)
project(Parafoil_Dynamics)

find_package(OpenCV REQUIRED)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-std=c++14 -O3")
include_directories("/usr/local/include/eigen3")

add_executable(parafoil_system src/parafoil_system.cc src/Parafoil_Aerodynamics.cc)
target_link_libraries(parafoil_system ${OpenCV_LIBS})
