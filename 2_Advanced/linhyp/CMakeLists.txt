

cmake_minimum_required(VERSION 3.1)

include_directories(/usr/local/Cellar/open-mpi/4.0.0/include/)
find_package(NGSolve CONFIG REQUIRED
        HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Users/max/Dev/ngsolve-install/Contents/Resources/CMake C:/netgen
)

add_ngsolve_python_module(liblinhyp linhypDG.cpp)
# add_ngsolve_python_module(libequilibrate equilibrate_simple.cpp)
