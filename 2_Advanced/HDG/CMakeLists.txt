
cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)

add_library(HDG SHARED HDGLaplace.cpp HDGElasticity)
target_link_libraries(HDG ngsolve)
